#!/usr/bin/env Rscript
#
# get_API_data.R
#
# (Adapted from data_prep_optimized.R with fixes to ensure behavior matches the
#  original dplyr-based script while keeping data.table performance improvements)
#
# See original header in your code for usage notes.
#

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(readr)      # for safe read_csv fallback if needed
  library(readxl)     # for NZTCS
  library(stringr)
  library(stringdist)
  library(sf)
  library(insect)
  library(wilderlab)
  library(dplyr)      # used minimally for recode; could be removed
})

# ---------- Helpers ----------
`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nzchar(as.character(a))) a else b

msg <- function(...) cat(sprintf(...), "\n")

# Read keys from CSV or environment
read_keys <- function(keys_file = Sys.getenv("WILDER_KEYS_FILE", unset = "Data/wilder_keys.csv")) {
  if (file.exists(keys_file)) {
    kv <- tryCatch({
      dt <- fread(keys_file, header = FALSE, sep = ",", showProgress = FALSE)
      # support both comma and simple two-column files
      if (ncol(dt) < 2) stop("Keys CSV must have at least two columns")
      names_vec <- tolower(trimws(as.character(dt[[1]])))
      vals_vec  <- trimws(as.character(dt[[2]]))
      setNames(vals_vec, names_vec)
    }, error = function(e) {
      warning("Failed to read keys file: ", e$message)
      return(NULL)
    })
    if (!is.null(kv)) {
      key    <- kv[["key"]]     %||% kv[["api key"]] %||% NA_character_
      secret <- kv[["secret"]]  %||% kv[["api secret"]] %||% NA_character_
      xapikey <- kv[["xapikey"]] %||% kv[["x-api-key"]] %||% kv[["x_api_key"]] %||% NA_character_
      return(list(key = key, secret = secret, xapikey = xapikey))
    }
  }
  # fallback env
  list(
    key    = Sys.getenv("WILDER_KEY",    unset = NA_character_),
    secret = Sys.getenv("WILDER_SECRET", unset = NA_character_),
    xapikey= Sys.getenv("WILDER_XAPIKEY", unset = Sys.getenv("WILDER_X_API_KEY", unset = NA_character_))
  )
}

# Robust fetch with retries/backoff for get_wilderdata
fetch_wilder_safe <- function(type, ..., max_attempts = 4, base_sleep = 1) {
  attempt <- 1L
  repeat {
    res <- try(get_wilderdata(type, ...), silent = TRUE)
    if (!inherits(res, "try-error")) return(res)
    msg("Attempt %d for %s failed: %s", attempt, type, conditionMessage(attr(res, "condition")))
    if (attempt >= max_attempts) {
      warning("Giving up on ", type, " after ", attempt, " attempts")
      return(NULL)
    }
    Sys.sleep(base_sleep * attempt)
    attempt <- attempt + 1L
  }
}

# Safe per-taxid lineage fetch to avoid insect::get_lineage assertions
fetch_lineage_safe <- function(taxid, tdb) {
  if (is.na(taxid) || taxid == "") return(list())
  taxid_num <- suppressWarnings(as.numeric(taxid))
  if (is.na(taxid_num)) return(list())
  res <- tryCatch(insect::get_lineage(taxid_num, tdb), error = function(e) list())
  if (is.null(res)) res <- list()
  res
}

# species cleaning helper (vectorized)
clean_species_names <- function(x) {
  x <- tolower(x)
  x <- str_replace_all(x, "\\(.*?\\)", "")   # remove parentheses content
  x <- str_replace_all(x, "[[:punct:]]", " ") # remove punctuation to spaces
  x <- str_replace_all(x, "\\bsp\\.?\\b", " ")
  x <- str_replace_all(x, "\\bspecies\\b", " ")
  x <- str_replace_all(x, "\\s+", " ")
  trimws(x)
}

# sanitize layer names for GPKG (kept here for convenience if needed)
sanitize_name <- function(x, maxlen = 50) {
  nm <- iconv(as.character(x), to = "ASCII//TRANSLIT")
  nm <- str_replace_all(nm, "[^A-Za-z0-9_]", "_")
  nm <- str_replace_all(nm, "_+", "_")
  nm <- str_trim(nm, side = "both")
  nm <- str_replace_all(nm, "^_+|_+$", "")
  nm <- substr(nm, 1, maxlen)
  if (nzchar(nm)) nm else "unknown"
}

# ---------- Start ----------
keys <- read_keys()
key <- keys$key; secret <- keys$secret; xapikey <- keys$xapikey
if (is.na(key) || is.na(secret) || is.na(xapikey)) warning("API keys appear missing. Set Data/wilder_keys.csv or environment vars.")

today <- Sys.Date()
msg("Starting data prep: ", format(today))

# ---------- 1) Read NZTCS lookup ----------
msg("Reading NZTCS spreadsheet...")
nztcs_raw <- read_excel("Data/NZTCS.xlsx", sheet = "Exported Data")
nztcs_dt <- as.data.table(nztcs_raw)
# Normalize NZTCS columns we will use
norm_chr <- function(x) ifelse(is.na(x), NA_character_, str_squish(as.character(x)))
cap_case <- function(x) ifelse(is.na(x), NA_character_, str_to_title(x))
nztcs_dt[, species_nztcs := norm_chr(`Current Species Name`)]
nztcs_dt[, Genus := cap_case(norm_chr(Genus))]
nztcs_dt[, Family := cap_case(norm_chr(Family))]
nztcs_dt[, Order := cap_case(norm_chr(Order))]
nztcs_dt[, Class := cap_case(norm_chr(Class))]
nztcs_dt[, Phylum := cap_case(norm_chr(Phylum))]
nztcs_dt[, Status := norm_chr(Status)]
nztcs_dt[, Category := norm_chr(Category)]
nztcs_dt[, BioStatus := norm_chr(`Bio Status`)]
nztcs_dt[, YearAssessed := suppressWarnings(as.integer(`Year Assessed`))]
nztcs_dt[, ThreatReport := norm_chr(`Report Name`)]
# reduce to one record per species (most recent YearAssessed)
setorder(nztcs_dt, -YearAssessed)
nztcs_sp <- unique(nztcs_dt, by = "species_nztcs")
nztcs_sp <- nztcs_sp[, .(species_nztcs, Status, Category, BioStatus, ThreatReport, Genus, Family, Order, Class, Phylum)]
nztcs_sp[, species_clean := clean_species_names(species_nztcs)]

# ---------- 2) Fetch Wilderlab jobs/samples/taxa (with retries) ----------
msg("Fetching Wilderlab jobs/samples/taxa...")
jobs <- fetch_wilder_safe("jobs", key = key, secret = secret, xapikey = xapikey)
samples <- fetch_wilder_safe("samples", key = key, secret = secret, xapikey = xapikey)
taxa <- fetch_wilder_safe("taxa", key = key, secret = secret, xapikey = xapikey)

if (is.null(jobs) || is.null(samples) || is.null(taxa)) {
  stop("Failed to fetch essential Wilderlab endpoints (jobs/samples/taxa). Aborting.")
}

# Convert samples/taxa to data.table
setDT(samples); setDT(taxa)

# ---------- 3) Fetch per-job records robustly (list of data.tables) ----------
msg("Fetching per-job records (%d jobs)...", nrow(jobs))
fetch_records_for_job <- function(jid) {
  df <- fetch_wilder_safe("records", JobID = jid, key = key, secret = secret, xapikey = xapikey, max_attempts = 4)
  if (is.null(df)) return(NULL)
  # coerce to data.table and normalize factor/logical to character
  setDT(df)
  facs <- names(df)[vapply(df, is.factor, logical(1))]
  if (length(facs)) df[, (facs) := lapply(.SD, as.character), .SDcols = facs]
  logs <- names(df)[vapply(df, is.logical, logical(1))]
  if (length(logs)) df[, (logs) := lapply(.SD, as.character), .SDcols = logs]
  # ensure Rank exists
  if (!"Rank" %in% names(df)) df[, Rank := NA_character_]
  return(df)
}

# iterate jobs serially to avoid too many parallel requests (and hitting rate limits)
records_list <- vector("list", length = nrow(jobs))
for (i in seq_len(nrow(jobs))) {
  jid <- jobs$JobID[i]
  records_list[[i]] <- tryCatch(fetch_records_for_job(jid), error = function(e) {
    warning("JobID ", jid, " fetch failed: ", e$message); NULL
  })
  if (i %% 50 == 0) msg("Fetched %d/%d jobs...", i, nrow(jobs))
}
# compact
records_list <- Filter(Negate(is.null), records_list)
msg("Fetched records for %d jobs (non-null).", length(records_list))

# ---------- 4) Read public S3 CSVs (fast) ----------
msg("Reading public S3 CSVs (samples, records) with fread...")
# use data.table::fread as it's fastest; if network causes issues, fallback to read_csv
public_samples <- tryCatch(fread("http://s3.ap-southeast-2.amazonaws.com/wilderlab.publicdata/samples.csv", showProgress = FALSE),
                           error = function(e) { warning("fread public_samples failed: ", e$message); as.data.table(read_csv("http://s3.ap-southeast-2.amazonaws.com/wilderlab.publicdata/samples.csv")) })
public_records <- tryCatch(fread("http://s3.ap-southeast-2.amazonaws.com/wilderlab.publicdata/records.csv", showProgress = FALSE),
                           error = function(e) { warning("fread public_records failed: ", e$message); as.data.table(read_csv("http://s3.ap-southeast-2.amazonaws.com/wilderlab.publicdata/records.csv")) })

# ---------- 5) Combine per-job records fast with rbindlist and harmonise types ----------
msg("Combining per-job records (rbindlist)...")
if (length(records_list) == 0) {
  records_dt <- data.table()   # empty
} else {
  records_dt <- rbindlist(records_list, fill = TRUE, use.names = TRUE)
}

# Convert important columns to consistent types quickly.
num_cols <- intersect(c("Count", "TaxID", "Latitude", "Longitude", "TICIQuantile", "TICINoSeqs"), names(records_dt))
if (length(num_cols)) {
  for (c in num_cols) records_dt[, (c) := suppressWarnings(as.numeric(get(c)))]
}
char_cols <- intersect(c("UID", "HID", "ClientSampleID", "Report", "Rank", "Name", "CommonName", "Group"), names(records_dt))
if (length(char_cols)) records_dt[, (char_cols) := lapply(.SD, as.character), .SDcols = char_cols]

msg("Per-job records combined: %d rows, %d cols", nrow(records_dt), ncol(records_dt))

# ---------- 6) Harmonise public_records (convert to same schema) ----------
msg("Harmonising public_records into combined records table...")
setDT(public_records)
# ensure same columns exist: add missing columns in each table as NA
cols_union <- union(names(records_dt), names(public_records))
for (c in setdiff(cols_union, names(records_dt))) records_dt[, (c) := NA_character_]
for (c in setdiff(cols_union, names(public_records))) public_records[, (c) := NA_character_]

# coerce public_records columns to types matching records_dt: numeric or character
for (c in cols_union) {
  if (c %in% names(records_dt) && is.numeric(records_dt[[c]])) {
    public_records[, (c) := suppressWarnings(as.numeric(get(c)))]
  } else {
    public_records[, (c) := as.character(get(c))]
    records_dt[, (c) := as.character(get(c))]
  }
}

# Now rbindlist them
all_records_dt <- rbindlist(list(records_dt, public_records), use.names = TRUE, fill = TRUE)
msg("All records combined: %d rows", nrow(all_records_dt))

# Free memory of intermediates
rm(records_list, records_dt, public_records); gc()

# ---------- 7) Prepare samples and combine DOC + public samples ----------
msg("Combine samples...")
setDT(samples)
samples[, MakeDataPublic := ifelse(MakeDataPublic == 1, "Private", "Public")]
samples[, DOC_Data := "Yes"]
samples[, CollectionDate := as.IDate(CollectionDate)]
samples[, Latitude := suppressWarnings(as.numeric(Latitude))]
samples[, Longitude := suppressWarnings(as.numeric(Longitude))]

setDT(public_samples)
public_samples[, MakeDataPublic := "Public"]
public_samples[, DOC_Data := "No"]
public_samples[, CollectionDate := as.IDate(CollectionDate)]
public_samples[, Latitude := suppressWarnings(as.numeric(Latitude))]
public_samples[, Longitude := suppressWarnings(as.numeric(Longitude))]

# drop public samples already present in DOC samples by UID
if ("UID" %in% names(public_samples) && "UID" %in% names(samples)) {
  public_samples <- public_samples[!UID %in% samples$UID]
}

all_samples <- rbindlist(list(samples, public_samples), use.names = TRUE, fill = TRUE)
setDT(all_samples)
msg("All samples combined: %d rows", nrow(all_samples))

# ---------- 8) Join sample metadata into all_records_dt using keys (fast data.table join) ----------
# Ensure UID is character in both tables
if ("UID" %in% names(all_samples))  all_samples[, UID := as.character(UID)]
if ("UID" %in% names(all_records_dt)) all_records_dt[, UID := as.character(UID)]

# Use data.table update-join syntax to copy sample fields into all_records_dt
# This is the correct, efficient equivalent of merged_records$Latitude <- samples$Latitude[match(...)]
if ("UID" %in% names(all_samples) && "UID" %in% names(all_records_dt)) {
  setkey(all_samples, UID)
  setkey(all_records_dt, UID)
  # assign-from-i style update-join: i.<col> refers to all_samples' columns
  all_records_dt[all_samples, Latitude := i.Latitude, on = "UID"]
  all_records_dt[all_samples, Longitude := i.Longitude, on = "UID"]
  all_records_dt[all_samples, ClientSampleID := i.ClientSampleID, on = "UID"]
  all_records_dt[all_samples, Report := i.Report, on = "UID"]
} else {
  # ensure columns exist so downstream code doesn't break
  for (c in c("Latitude", "Longitude", "ClientSampleID", "Report")) {
    if (!c %in% names(all_records_dt)) all_records_dt[, (c) := NA_character_]
  }
}

msg("Joined sample metadata (coordinates, ClientSampleID) into records.")

# ---------- 9) Add taxonomy lineages (use tdb from taxa and safe per-taxid calls) ----------
msg("Adding taxonomic lineages (insect::get_lineage) ...")
setDT(taxa)
tdb <- taxa[, 1:4, with = FALSE]
setnames(tdb, old = names(tdb), new = c("taxID", "parent_taxID", "rank", "name"))
tdb[, taxID := suppressWarnings(as.numeric(taxID))]

# fetch lineage per unique TaxID to avoid duplicate expensive calls
if (!"TaxID" %in% names(all_records_dt)) {
  all_records_dt[, TaxID := NA_character_]
}
unique_taxids <- unique(na.omit(all_records_dt$TaxID))
msg("Unique TaxIDs to lookup: %d", length(unique_taxids))

msg("Fetching lineages for unique TaxIDs and building lookup table...")

lineage_rows <- lapply(unique_taxids, function(tid) {
  lin <- fetch_lineage_safe(tid, tdb)  # returns list or empty list on failure
  safe_get <- function(l, rank) {
    if (is.null(l) || length(l) == 0) return(NA_character_)
    nms <- names(l)
    if (!is.null(nms) && rank %in% nms) {
      val <- l[[rank]]
      return(if (is.null(val)) NA_character_ else as.character(val))
    }
    return(NA_character_)
  }
  list(
    TaxID = as.character(tid),
    phylum = safe_get(lin, "phylum"),
    class  = safe_get(lin, "class"),
    order  = safe_get(lin, "order"),
    family = safe_get(lin, "family"),
    genus  = safe_get(lin, "genus"),
    species = safe_get(lin, "species")
  )
})

lineage_dt <- rbindlist(lineage_rows, fill = TRUE)
lineage_dt <- unique(lineage_dt, by = "TaxID")
lineage_dt[, TaxID := as.character(TaxID)]

# Merge lineage_dt into all_records_dt using a left join (keep all records)
if ("TaxID" %in% names(all_records_dt)) all_records_dt[, TaxID := as.character(TaxID)]

setkey(lineage_dt, TaxID)
setkey(all_records_dt, TaxID)
all_records_dt <- merge(all_records_dt, lineage_dt, by = "TaxID", all.x = TRUE, sort = FALSE)

msg("Lineage lookup merged. Ensuring lineage columns exist...")

# Ensure lineage columns exist (if merge didn't create them for some reason)
expected_lineage_cols <- c("phylum", "class", "order", "family", "genus", "species")
for (col in expected_lineage_cols) {
  if (!col %in% names(all_records_dt)) {
    all_records_dt[, (col) := NA_character_]
  } else {
    all_records_dt[, (col) := ifelse(is.na(get(col)), NA_character_, str_squish(as.character(get(col))))]
  }
}

# ---------- 10) Clean taxa text (normalize capitals etc.) ----------
msg("Cleaning taxonomy text...")
all_records_dt[, phylum := ifelse(is.na(phylum), NA_character_, str_to_title(str_squish(as.character(phylum))))]
all_records_dt[, class  := ifelse(is.na(class),  NA_character_, str_to_title(str_squish(as.character(class))))]
all_records_dt[, order  := ifelse(is.na(order),  NA_character_, str_to_title(str_squish(as.character(order))))]
all_records_dt[, family := ifelse(is.na(family), NA_character_, str_to_title(str_squish(as.character(family))))]
all_records_dt[, genus  := ifelse(is.na(genus),  NA_character_, str_to_title(str_squish(as.character(genus))))]

class_synonyms <- list("Actinopteri" = "Actinopterygii")
if ("class" %in% names(all_records_dt)) {
  for (k in names(class_synonyms)) {
    all_records_dt[class == k, class := class_synonyms[[k]]]
  }
}
if ("species" %in% names(all_records_dt)) {
  all_records_dt[, species := as.character(species)]
  all_records_dt[species == "Galaxias sp. D (Allibone et al., 1996)", species := 'Galaxias "species D"']
}

# ---------- 11) Fuzzy match species names to NZTCS (efficient unique-based matching) ----------
msg("Fuzzy matching species names to NZTCS...")
all_records_dt[, species_clean := clean_species_names(species)]
unique_species_clean <- unique(na.omit(all_records_dt$species_clean))
msg("Unique cleaned species names in records: %d", length(unique_species_clean))
choices <- nztcs_sp$species_clean
max_dist <- 0.12

matched_choices <- vapply(unique_species_clean, function(x) {
  idx <- stringdist::amatch(x, choices, method = "jw", maxDist = max_dist)
  if (is.na(idx)) NA_character_ else choices[idx]
}, FUN.VALUE = character(1), USE.NAMES = FALSE)

lookup_dt <- data.table(species_clean = unique_species_clean, matched_species = matched_choices)

# Merge lookup back into all_records_dt (left join) so we keep all records
all_records_dt <- merge(all_records_dt, lookup_dt, by = "species_clean", all.x = TRUE, sort = FALSE)

# join NZTCS attributes for matched species (left join)
all_records_dt <- merge(all_records_dt, nztcs_sp[, .(species_clean, species_nztcs, Status, Category, BioStatus, ThreatReport, Genus, Family, Order, Class, Phylum)], 
                        by = "species_clean", all.x = TRUE, sort = FALSE)

# ---------- 12) Spatial joins for Nga Awa & Regional Council using sf (convert smaller sample table) ----------
msg("Spatial joins: Nga Awa catchment and Regional Council (samples-level)...")
NA_shp_path <- "Data/Nga Awa shapefiles/DOC_NgāAwa_RiverSites_20250122_n14.shp"
RC_shp_path <- "Data/Regional Council shapefiles/regional-council-2022-generalised.shp"
NA_polygons <- tryCatch(st_read(NA_shp_path, quiet = TRUE), error = function(e) { stop("Failed to read Nga Awa shapefile: ", e$message) })
RC_polygons <- tryCatch(st_read(RC_shp_path, quiet = TRUE), error = function(e) { stop("Failed to read Regional Council shapefile: ", e$message) })

samples_valid <- copy(all_samples)[!is.na(Longitude) & !is.na(Latitude)]
if (nrow(samples_valid) > 0) {
  samples_sf <- st_as_sf(samples_valid, coords = c("Longitude", "Latitude"), crs = 4167, remove = FALSE)
  # transform into polygon CRS and spatial join
  samples_sf <- st_transform(samples_sf, st_crs(NA_polygons))
  NA_joined <- st_join(samples_sf, NA_polygons, left = TRUE)
  samples_valid[, Nga_Awa_Catchment := ifelse(is.na(NA_joined$Waterway_N), "not Nga Awa", NA_joined$Waterway_N)]
  
  samples_sf2 <- st_as_sf(samples_valid, coords = c("Longitude", "Latitude"), crs = 4167, remove = FALSE)
  samples_sf2 <- st_transform(samples_sf2, st_crs(RC_polygons))
  RC_joined <- st_join(samples_sf2, RC_polygons, left = TRUE)
  samples_valid[, Regional_Council := ifelse(is.na(RC_joined$REGC2022_1), "None", RC_joined$REGC2022_1)]
} else {
  samples_valid <- copy(all_samples)
  samples_valid[, Nga_Awa_Catchment := "not Nga Awa"]
  samples_valid[, Regional_Council := "None"]
}

# Attach back to all_records_dt by SID (prefer) or UID (fallback)
if ("SID" %in% names(all_records_dt) && "SID" %in% names(samples_valid)) {
  setkey(samples_valid, SID)
  setkey(all_records_dt, SID)
  all_records_dt[samples_valid, Nga_Awa_Catchment := i.Nga_Awa_Catchment, on = "SID"]
  all_records_dt[samples_valid, Regional_Council := i.Regional_Council, on = "SID"]
} else if ("UID" %in% names(all_records_dt) && "UID" %in% names(samples_valid)) {
  setkey(samples_valid, UID)
  setkey(all_records_dt, UID)
  all_records_dt[samples_valid, Nga_Awa_Catchment := i.Nga_Awa_Catchment, on = "UID"]
  all_records_dt[samples_valid, Regional_Council := i.Regional_Council, on = "UID"]
} else {
  if (!"Nga_Awa_Catchment" %in% names(all_records_dt)) all_records_dt[, Nga_Awa_Catchment := NA_character_]
  if (!"Regional_Council" %in% names(all_records_dt)) all_records_dt[, Regional_Council := NA_character_]
}

all_records_dt[, `:=`(
  Latitude = suppressWarnings(as.numeric(Latitude)),
  Longitude = suppressWarnings(as.numeric(Longitude)),
  Count = suppressWarnings(as.numeric(Count))
)]

# ---------- 13) Create export summary (efficient data.table aggregation) ----------
msg("Summarising by Report and TaxID (data.table aggregation)...")
required_cols <- c("Report", "TaxID", "UID", "Count", "ClientSampleID", "Rank", "Name", "CommonName",
                   "Group", "Latitude", "Longitude", "CollectionDate", "Status", "Category", "ThreatReport",
                   "CollectedBy", "DOC_Data", "MakeDataPublic", "Nga_Awa_Catchment", "Regional_Council",
                   "phylum", "class", "order", "family", "genus", "species", "Wilderlab_Sp_name", "species_nztcs")
for (c in required_cols) if (!c %in% names(all_records_dt)) all_records_dt[, (c) := NA_character_]

DT <- all_records_dt  # working name

# uid_counts per Report
uid_counts_dt <- DT[, .(total_UID = uniqueN(UID)), by = Report]

# summary by Report + TaxID
summary_dt <- DT[, .(
  unique_UID_count = uniqueN(UID),
  UID_list = paste(unique(UID), collapse = "-"),
  sum_count = sum(if ("Count" %in% names(.SD)) as.numeric(Count) else 0, na.rm = TRUE),
  ClientSampleID = as.character(first(ClientSampleID)),
  Rank = as.character(first(Rank)),
  Name = as.character(first(Name)),
  CommonName = as.character(first(CommonName)),
  Group = as.character(first(Group)),
  Latitude = as.numeric(first(Latitude)),
  Longitude = as.numeric(first(Longitude)),
  CollectionDate = first(CollectionDate),
  Status = as.character(first(Status)),
  Category = as.character(first(Category)),
  ThreatReport = as.character(first(ThreatReport)),
  CollectedBy = as.character(first(CollectedBy)),
  DOC_Data = as.character(first(DOC_Data)),
  MakeDataPublic = as.character(first(MakeDataPublic)),
  Nga_Awa_Catchment = as.character(first(Nga_Awa_Catchment)),
  Regional_Council = as.character(first(Regional_Council)),
  phylum = as.character(first(phylum)),
  class = as.character(first(class)),
  order = as.character(first(order)),
  family = as.character(first(family)),
  genus = as.character(first(genus)),
  species = as.character(first(species)),
  Wilderlab_Sp_name = as.character(first(Wilderlab_Sp_name)),
  NZTC_Sp_name = as.character(first(species_nztcs))
), by = .(Report, TaxID)]

# join total_UID into summary_dt and compute mean_count
summary_dt <- merge(summary_dt, uid_counts_dt, by = "Report", all.x = TRUE, sort = FALSE)
summary_dt[, mean_count := sum_count / total_UID]

out_cols <- c("Name", "CommonName", "ClientSampleID", "sum_count", "mean_count",
              "unique_UID_count", "total_UID", "Rank", "Group", "Status", "Category",
              "Latitude", "Longitude", "CollectionDate", "ThreatReport", "CollectedBy",
              "DOC_Data", "MakeDataPublic", "Nga_Awa_Catchment", "Regional_Council",
              "phylum", "class", "order", "family", "genus", "species",
              "Wilderlab_Sp_name", "NZTC_Sp_name", "TaxID", "Report")
out_cols <- intersect(out_cols, names(summary_dt))
setcolorder(summary_dt, out_cols)

summary_df <- as.data.frame(summary_dt)

# Ensure exported column names match what the Shiny app expects.

# Map from current summary_df column names to the Shiny-friendly names used by the app
rename_map <- c(
  "Name" = "Latin Name",
  "CommonName" = "Common Name",
  "ClientSampleID" = "Sample Name",
  "sum_count" = "Total Reads",
  "mean_count" = "Average reads",
  "unique_UID_count" = "Hit number",
  "total_UID" = "Total hits",
  "Rank" = "Taxonomic Rank",
  "Group" = "Taxon Group",
  "Status" = "Threat Status",
  "Category" = "Threat Category",
  "CollectionDate" = "Date",
  "ThreatReport" = "Threat Document",
  "CollectedBy" = "Collector",
  "DOC_Data" = "DOC Data",
  "MakeDataPublic" = "Public/Private",
  "Nga_Awa_Catchment" = "Nga Awa",
  "Regional_Council" = "Regional Council",
  "Wilderlab_Sp_name" = "Wilderlab Sp Name",
  "NZTC_Sp_name" = "NZTC Sp Name",
  "Report" = "Wilderlab Report"
)

# Only rename columns that actually exist in summary_df
existing_old_names <- intersect(names(rename_map), names(summary_df))
if (length(existing_old_names) > 0) {
  new_names <- rename_map[existing_old_names]
  # perform rename in place
  names(summary_df)[match(existing_old_names, names(summary_df))] <- unname(new_names)
}

# Ensure the Shiny app's expected columns exist (create missing ones as NA so filters won't error)
expected_cols_for_shiny <- c(
  "Latin Name","Common Name","Sample Name","Total Reads","Average reads","Hit number","Total hits",
  "Taxonomic Rank","Taxon Group","Threat Status","Threat Category","Date","Threat Document","Collector",
  "DOC Data","Public/Private","Nga Awa","Regional Council","Wilderlab Sp Name","NZTC Sp Name",
  "TaxID","Wilderlab Report","Latitude","Longitude"
)
missing_cols <- setdiff(expected_cols_for_shiny, names(summary_df))
if (length(missing_cols) > 0) {
  for (c in missing_cols) summary_df[[c]] <- NA
}

# Optionally re-order columns to put the most-used columns first (not required, but keeps parity with original)
preferred_order <- intersect(expected_cols_for_shiny, names(summary_df))
remaining <- setdiff(names(summary_df), preferred_order)
summary_df <- summary_df[, c(preferred_order, remaining)]

# ---------- 14) Save outputs (dated and unversioned copy for the Shiny app) ----------
date_suffix <- format(today, "%d%m%y")
out_file <- file.path("Data", paste0("records_", date_suffix, ".RDS"))
msg("Saving dated RDS -> ", out_file)
saveRDS(summary_df, out_file)

unversioned <- file.path("Data", "records.rds")
msg("Saving unversioned RDS -> ", unversioned)
saveRDS(summary_df, unversioned)

msg("Completed. Summary rows: %d", nrow(summary_df))