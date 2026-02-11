#!/usr/bin/env Rscript
#
# get_API_data.R
#
# Produces Data/records.rds compatible with ShinyApp.R
#
# Ensures:
# - All records from DOC API are marked DOC Data = "Yes"
# - All records from public S3 are marked DOC Data = "No"
# - Public S3 rows are retained even when UIDs overlap DOC API rows
# - Final saved records.rds contains the exact column names the Shiny app expects
#
# Dependencies: data.table, readr, readxl, stringr, stringdist, sf, insect, wilderlab, dplyr
#

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(readxl)
  library(stringr)
  library(stringdist)
  library(sf)
  library(insect)
  library(wilderlab)
  library(dplyr)
  library(wilderlab)
})

# ---------- Helpers ----------
`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nzchar(as.character(a))) a else b
msg <- function(...) cat(sprintf(...), "\n")

# Read keys from CSV or environment
read_keys <- function(keys_file = Sys.getenv("WILDER_KEYS_FILE", unset = "Data/wilder_keys.csv")) {
  if (file.exists(keys_file)) {
    kv <- tryCatch({
      dt <- fread(keys_file, header = FALSE, sep = ",", showProgress = FALSE)
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
  x <- tolower(as.character(x))
  x <- str_replace_all(x, "\\(.*?\\)", "")   # remove parentheses content
  x <- str_replace_all(x, "[[:punct:]]", " ") # replace punctuation with spaces
  x <- str_replace_all(x, "\\bsp\\.?\\b", " ")
  x <- str_replace_all(x, "\\bspecies\\b", " ")
  x <- str_replace_all(x, "\\s+", " ")
  trimws(x)
}

# small normalisers
norm_chr <- function(x) ifelse(is.na(x), NA_character_, str_squish(as.character(x)))
cap_case <- function(x) ifelse(is.na(x), NA_character_, str_to_title(x))

# ---------- Start ----------
keys <- read_keys()
key <- keys$key; secret <- keys$secret; xapikey <- keys$xapikey
if (is.na(key) || is.na(secret) || is.na(xapikey)) msg("API keys appear missing — public S3 will still be used.")

today <- Sys.Date()
msg("Starting data prep: ", format(today))

# ---------- 1) Read NZTCS spreadsheet ----------
msg("Reading NZTCS spreadsheet...")
nztcs_raw <- read_excel("Data/NZTCS.xlsx", sheet = "Exported Data")
nztcs_dt <- as.data.table(nztcs_raw)
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
setorder(nztcs_dt, -YearAssessed)
nztcs_sp <- unique(nztcs_dt, by = "species_nztcs")
nztcs_sp <- nztcs_sp[, .(species_nztcs, Status, Category, BioStatus, ThreatReport, Genus, Family, Order, Class, Phylum)]
nztcs_sp[, species_clean := clean_species_names(species_nztcs)]

# ---------- 2) Fetch Wilderlab jobs/samples/taxa (with retries) ----------
msg("Fetching Wilderlab jobs/samples/taxa...")
jobs <- tryCatch(fetch_wilder_safe("jobs", key = key, secret = secret, xapikey = xapikey), error = function(e) NULL)
samples <- tryCatch(fetch_wilder_safe("samples", key = key, secret = secret, xapikey = xapikey), error = function(e) NULL)
taxa <- tryCatch(fetch_wilder_safe("taxa", key = key, secret = secret, xapikey = xapikey), error = function(e) NULL)

# If API calls failed, continue but rely on public S3
if (is.null(samples)) samples <- data.table()
if (is.null(taxa)) taxa <- data.table()
if (is.null(jobs)) jobs <- data.table()

setDT(samples); setDT(taxa); setDT(jobs)

# ensure samples from API have provenance and DOC flag
if (nrow(samples) > 0) {
  samples[, Source := "DOC_API"]
  samples[, DOC_Data := "Yes"]
  samples[, MakeDataPublic := ifelse(as.character(MakeDataPublic) == "1" | as.numeric(MakeDataPublic) == 1, "Public", "Private")]
  samples[, CollectionDate := as.IDate(CollectionDate)]
  samples[, Latitude := suppressWarnings(as.numeric(Latitude))]
  samples[, Longitude := suppressWarnings(as.numeric(Longitude))]
}

# ---------- 3) Fetch per-job records robustly ----------
msg("Fetching per-job records...")
records_list <- list()
if (nrow(jobs) > 0) {
  for (i in seq_len(nrow(jobs))) {
    jid <- jobs$JobID[i]
    df <- tryCatch(fetch_wilder_safe("records", JobID = jid, key = key, secret = secret, xapikey = xapikey),
                   error = function(e) { warning("Job ", jid, " fetch failed: ", e$message); NULL })
    if (!is.null(df)) {
      setDT(df)
      df[, Source := "DOC_API"]
      df[, DOC_Data := "Yes"]
      # normalize types
      if ("Count" %in% names(df)) df[, Count := suppressWarnings(as.numeric(Count))]
      records_list[[length(records_list) + 1]] <- df
    }
    if (i %% 50 == 0 && i > 0) msg("Fetched %d/%d jobs...", i, nrow(jobs))
  }
}

records_dt <- if (length(records_list) > 0) rbindlist(records_list, fill = TRUE, use.names = TRUE) else data.table()

# ---------- 4) Read public S3 CSVs ----------
msg("Reading public S3 CSVs (samples, records)...")
public_samples <- tryCatch(fread("http://s3.ap-southeast-2.amazonaws.com/wilderlab.publicdata/samples.csv", showProgress = FALSE),
                           error = function(e) as.data.table(read_csv("http://s3.ap-southeast-2.amazonaws.com/wilderlab.publicdata/samples.csv")))
public_records <- tryCatch(fread("http://s3.ap-southeast-2.amazonaws.com/wilderlab.publicdata/records.csv", showProgress = FALSE),
                           error = function(e) as.data.table(read_csv("http://s3.ap-southeast-2.amazonaws.com/wilderlab.publicdata/records.csv")))

setDT(public_samples); setDT(public_records)
if (nrow(public_samples) > 0) {
  public_samples[, Source := "Public_S3"]
  public_samples[, DOC_Data := "No"]
  public_samples[, MakeDataPublic := "Public"]
  public_samples[, CollectionDate := as.IDate(CollectionDate)]
  public_samples[, Latitude := suppressWarnings(as.numeric(Latitude))]
  public_samples[, Longitude := suppressWarnings(as.numeric(Longitude))]
}
if (nrow(public_records) > 0) {
  public_records[, Source := "Public_S3"]
  public_records[, DOC_Data := "No"]
  # ensure Count numeric where present
  if ("Count" %in% names(public_records)) public_records[, Count := suppressWarnings(as.numeric(Count))]
}

# ---------- 5) Combine records: DO NOT drop public rows that share UID with DOC API ----------
msg("Combining records from API and public S3...")
# Harmonise column sets prior to rbind
cols_union <- union(names(records_dt), names(public_records))
for (c in setdiff(cols_union, names(records_dt))) records_dt[, (c) := NA_character_]
for (c in setdiff(cols_union, names(public_records))) public_records[, (c) := NA_character_]

# coerce types: numeric where appropriate in the union
for (c in cols_union) {
  if (c %in% names(records_dt) && is.numeric(records_dt[[c]])) {
    public_records[, (c) := suppressWarnings(as.numeric(get(c)))]
  } else {
    public_records[, (c) := as.character(get(c))]
    records_dt[, (c) := as.character(get(c))]
  }
}

all_records_dt <- rbindlist(list(records_dt, public_records), use.names = TRUE, fill = TRUE)
msg("All records combined: %d rows", nrow(all_records_dt))

# Remove intermediates
rm(records_dt, public_records, records_list); gc()

# ---------- 6) Combine samples similarly ----------
msg("Combining samples from API and public S3...")
samples_union_cols <- union(names(samples), names(public_samples))
for (c in setdiff(samples_union_cols, names(samples))) samples[, (c) := NA_character_]
for (c in setdiff(samples_union_cols, names(public_samples))) public_samples[, (c) := NA_character_]
# coerce numeric/character as needed
for (c in samples_union_cols) {
  if (c %in% names(samples) && is.numeric(samples[[c]])) public_samples[, (c) := suppressWarnings(as.numeric(get(c)))] else {
    public_samples[, (c) := as.character(get(c))]
    samples[, (c) := as.character(get(c))]
  }
}
all_samples <- rbindlist(list(samples, public_samples), use.names = TRUE, fill = TRUE)
msg("All samples combined: %d rows", nrow(all_samples))
rm(samples, public_samples); gc()

# ---------- 7) Attach sample metadata into records ----------
msg("Joining sample metadata into records (UID + Source preferred, fallback to UID)...")
if ("UID" %in% names(all_samples)) all_samples[, UID := as.character(UID)]
if ("UID" %in% names(all_records_dt)) all_records_dt[, UID := as.character(UID)]
if (!"Source" %in% names(all_records_dt)) all_records_dt[, Source := NA_character_]
if (!"Source" %in% names(all_samples)) all_samples[, Source := NA_character_]

# Ensure CollectionDate column exists in all_samples
if(!"CollectionDate" %in% names(all_samples)) {
  all_samples[, CollectionDate := as.IDate(NA)]
}

# Join by UID+Source first
setkeyv(all_samples, c("UID", "Source"))
setkeyv(all_records_dt, c("UID", "Source"))
all_records_dt[all_samples, Latitude := i.Latitude, on = .(UID, Source)]
all_records_dt[all_samples, Longitude := i.Longitude, on = .(UID, Source)]
all_records_dt[all_samples, ClientSampleID := i.ClientSampleID, on = .(UID, Source)]
all_records_dt[all_samples, Report := i.Report, on = .(UID, Source)]
all_records_dt[all_samples, MakeDataPublic := i.MakeDataPublic, on = .(UID, Source)]
all_records_dt[all_samples, DOC_Data := i.DOC_Data, on = .(UID, Source)]
all_records_dt[all_samples, CollectionDate := i.CollectionDate, on = .(UID, Source)]  # ADD THIS LINE

# fallback join by UID only for missing metadata
missing_idx <- which(is.na(all_records_dt$ClientSampleID) & !is.na(all_records_dt$UID))
if (length(missing_idx) > 0 && "UID" %in% names(all_samples)) {
  tmp_samples <- copy(all_samples)
  setkey(tmp_samples, UID)
  setkey(all_records_dt, UID)
  all_records_dt[tmp_samples, Latitude := fifelse(is.na(Latitude), i.Latitude, Latitude), on = "UID"]
  all_records_dt[tmp_samples, Longitude := fifelse(is.na(Longitude), i.Longitude, Longitude), on = "UID"]
  all_records_dt[tmp_samples, ClientSampleID := fifelse(is.na(ClientSampleID), as.character(i.ClientSampleID), ClientSampleID), on = "UID"]
  all_records_dt[tmp_samples, Report := fifelse(is.na(Report), as.character(i.Report), Report), on = "UID"]
  all_records_dt[tmp_samples, MakeDataPublic := fifelse(is.na(MakeDataPublic), as.character(i.MakeDataPublic), MakeDataPublic), on = "UID"]
  all_records_dt[tmp_samples, DOC_Data := fifelse(is.na(DOC_Data), as.character(i.DOC_Data), DOC_Data), on = "UID"]
  all_records_dt[tmp_samples, CollectionDate := fifelse(is.na(CollectionDate), i.CollectionDate, CollectionDate), on = "UID"]  # ADD THIS LINE
}
msg("Sample metadata joined.")

# ---------- 8) Taxonomy lineage merge using get_lineages() ----------
msg("Adding taxonomic lineages using wilderlab::get_lineages()...")

# Extract unique TaxIDs and ensure they're numeric for get_lineages()
unique_taxids <- unique(na.omit(as.numeric(all_records_dt$TaxID)))

if (length(unique_taxids) > 0) {
  tryCatch({
    # Use the new get_lineages() function from wilderlab
    lineages <- get_lineages(unique_taxids)
    lineage_dt <- as.data.table(lineages)
    
    # Ensure TaxID column exists and is character in both tables for proper joining
    if (!"TaxID" %in% names(lineage_dt)) {
      lineage_dt[, TaxID := as.character(rownames(lineages))]
    } else {
      lineage_dt[, TaxID := as.character(TaxID)]
    }
    
    # Convert TaxID in records to character for consistent merging
    all_records_dt[, TaxID := as.character(TaxID)]
    
    # Merge lineages into records
    setkey(lineage_dt, TaxID)
    setkey(all_records_dt, TaxID)
    all_records_dt <- merge(all_records_dt, lineage_dt, by = "TaxID", all.x = TRUE, sort = FALSE)
    
    msg("Lineages merged successfully.")
  }, error = function(e) {
    msg("Warning: Failed to retrieve lineages: %s", e$message)
    msg("Proceeding without lineage data. Contact Wilderlab if taxonomy IDs are not recognized.")
  })
} else {
  msg("No TaxIDs found to retrieve lineages for.")
}

# Normalize taxa capitalization
if ("domain" %in% names(all_records_dt)) {
  all_records_dt[, domain := ifelse(is.na(domain), NA_character_, str_to_title(str_squish(as.character(domain))))]
}
if ("phylum" %in% names(all_records_dt)) {
  all_records_dt[, phylum := ifelse(is.na(phylum), NA_character_, str_to_title(str_squish(as.character(phylum))))]
}
if ("class" %in% names(all_records_dt)) {
  all_records_dt[, class := ifelse(is.na(class), NA_character_, str_to_title(str_squish(as.character(class))))]
}
if ("order" %in% names(all_records_dt)) {
  all_records_dt[, order := ifelse(is.na(order), NA_character_, str_to_title(str_squish(as.character(order))))]
}
if ("family" %in% names(all_records_dt)) {
  all_records_dt[, family := ifelse(is.na(family), NA_character_, str_to_title(str_squish(as.character(family))))]
}
if ("genus" %in% names(all_records_dt)) {
  all_records_dt[, genus := ifelse(is.na(genus), NA_character_, str_to_title(str_squish(as.character(genus))))]
}

# After the merge, consolidate the columns
all_records_dt[, c("Domain.x", "Phylum.x", "Class.x", "Order.x", "Family.x", "Genus.x") := NULL]

# Rename the .y columns to clean names
setnames(all_records_dt, 
         old = c("Phylum.y", "Class.y", "Order.y", "Family.y", "Genus.y"),
         new = c("Phylum", "Class", "Order", "Family", "Genus"))

# Create Species column: use Name when Rank == "species", otherwise NA
all_records_dt[, Species := fifelse(Rank == "species", Name, NA_character_)]

all_records_dt$Subfamily = NULL
all_records_dt$Superfamily = NULL
all_records_dt$species = NULL

# Now your original code will work
class_synonyms <- list("Actinopteri" = "Actinopterygii")
for (k in names(class_synonyms)) {
  all_records_dt[Class == k, Class := class_synonyms[[k]]]
}

all_records_dt[, Species := as.character(Species)]
all_records_dt[Species == "Galaxias sp. D (Allibone et al., 1996)", species := 'Galaxias "species D"']

class_synonyms <- list("Actinopteri" = "Actinopterygii")
for (k in names(class_synonyms)) {
  all_records_dt[as.character(Class) == k, Class := class_synonyms[[k]]]
}

# ---------- 9) Fuzzy-match to NZTCS ----------
msg("Fuzzy matching Species names to NZTCS...")
all_records_dt[, species_clean := clean_species_names(Species)]
unique_species_clean <- unique(na.omit(all_records_dt$species_clean))
choices <- nztcs_sp$species_clean
max_dist <- 0.12
matched_choices <- vapply(unique_species_clean, function(x) {
  idx <- stringdist::amatch(x, choices, method = "jw", maxDist = max_dist)
  if (is.na(idx)) NA_character_ else choices[idx]
}, FUN.VALUE = character(1), USE.NAMES = FALSE)
lookup_dt <- data.table(species_clean = unique_species_clean, matched_Species = matched_choices)
all_records_dt <- merge(all_records_dt, lookup_dt, by = "species_clean", all.x = TRUE, sort = FALSE)
all_records_dt <- merge(all_records_dt, nztcs_sp[, .(species_clean, species_nztcs, Status, Category, BioStatus, ThreatReport, Genus, Family, Order, Class, Phylum)], by = "species_clean", all.x = TRUE, sort = FALSE)

# ---------- 10) Spatial joins for Nga Awa & Regional Council ----------
msg("Spatial joins for Nga Awa & Regional Council...")
NA_shp_path <- "Data/Nga Awa shapefiles/DOC_NgāAwa_RiverSites_20250122_n14.shp"
RC_shp_path <- "Data/Regional Council shapefiles/regional-council-2022-generalised.shp"
NA_polygons <- tryCatch(st_read(NA_shp_path, quiet = TRUE), error = function(e) NULL)
RC_polygons <- tryCatch(st_read(RC_shp_path, quiet = TRUE), error = function(e) NULL)

# Use the combined samples table for spatial joins (faster)
samples_sf_source <- copy(all_samples)[!is.na(Longitude) & !is.na(Latitude)]
if (nrow(samples_sf_source) > 0 && !is.null(NA_polygons)) {
  samples_sf <- st_as_sf(samples_sf_source, coords = c("Longitude", "Latitude"), crs = 4167, remove = FALSE)
  samples_sf <- st_transform(samples_sf, st_crs(NA_polygons))
  NA_joined <- st_join(samples_sf, NA_polygons, left = TRUE)
  samples_sf_source[, Nga_Awa_Catchment := ifelse(is.na(NA_joined$Waterway_N), "not Nga Awa", NA_joined$Waterway_N)]
} else {
  samples_sf_source[, Nga_Awa_Catchment := "not Nga Awa"]
}

if (nrow(samples_sf_source) > 0 && !is.null(RC_polygons)) {
  samples_sf2 <- st_as_sf(samples_sf_source, coords = c("Longitude", "Latitude"), crs = 4167, remove = FALSE)
  samples_sf2 <- st_transform(samples_sf2, st_crs(RC_polygons))
  RC_joined <- st_join(samples_sf2, RC_polygons, left = TRUE)
  samples_sf_source[, Regional_Council := ifelse(is.na(RC_joined$REGC2022_1), "None", RC_joined$REGC2022_1)]
} else {
  samples_sf_source[, Regional_Council := "None"]
}

# Attach Nga_Awa & Regional_Council back to records by UID+Source first, then UID fallback
if ("SID" %in% names(all_records_dt) && "SID" %in% names(samples_sf_source)) {
  setkeyv(samples_sf_source, c("SID", "Source"))
  setkeyv(all_records_dt, c("SID", "Source"))
  all_records_dt[samples_sf_source, Nga_Awa_Catchment := i.Nga_Awa_Catchment, on = .(SID, Source)]
  all_records_dt[samples_sf_source, Regional_Council := i.Regional_Council, on = .(SID, Source)]
} else if ("UID" %in% names(all_records_dt) && "UID" %in% names(samples_sf_source)) {
  setkeyv(samples_sf_source, c("UID", "Source"))
  setkeyv(all_records_dt, c("UID", "Source"))
  all_records_dt[samples_sf_source, Nga_Awa_Catchment := i.Nga_Awa_Catchment, on = .(UID, Source)]
  all_records_dt[samples_sf_source, Regional_Council := i.Regional_Council, on = .(UID, Source)]
} else {
  if (!"Nga_Awa_Catchment" %in% names(all_records_dt)) all_records_dt[, Nga_Awa_Catchment := NA_character_]
  if (!"Regional_Council" %in% names(all_records_dt)) all_records_dt[, Regional_Council := NA_character_]
}

# ---------- 11) Prepare summary by Report + TaxID (matching original script's output) ----------
msg("Aggregating summary by Report and TaxID...")
required_cols <- c("Report", "TaxID", "UID", "Count", "ClientSampleID", "Rank", "Name", "CommonName",
                   "Group", "Latitude", "Longitude", "CollectionDate", "Status", "Category", "ThreatReport",
                   "CollectedBy", "DOC_Data", "MakeDataPublic", "Nga_Awa_Catchment", "Regional_Council",
                   "Phylum", "Class", "Order", "Family", "Genus", "Species", "Wilderlab_Sp_name", "species_nztcs")
for (c in required_cols) if (!c %in% names(all_records_dt)) all_records_dt[, (c) := NA_character_]

DT <- all_records_dt

uid_counts_dt <- DT[, .(total_UID = uniqueN(UID)), by = Report]

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
  Phylum = as.character(first(Phylum)),
  Class = as.character(first(Class)),
  Order = as.character(first(Order)),
  Family = as.character(first(Family)),
  Genus = as.character(first(Genus)),
  Species = as.character(first(Species)),
  Wilderlab_Sp_name = as.character(first(Wilderlab_Sp_name)),
  NZTC_Sp_name = as.character(first(species_nztcs))
), by = .(Report, TaxID)]

summary_dt <- merge(summary_dt, uid_counts_dt, by = "Report", all.x = TRUE, sort = FALSE)
summary_dt[, mean_count := sum_count / total_UID]

# Choose and Order columns that will feed into the Shiny app
out_cols <- c("Name", "CommonName", "ClientSampleID", "sum_count", "mean_count",
              "unique_UID_count", "total_UID", "Rank", "Group", "Status", "Category",
              "Latitude", "Longitude", "CollectionDate", "ThreatReport", "CollectedBy",
              "DOC_Data", "MakeDataPublic", "Nga_Awa_Catchment", "Regional_Council",
              "Phylum", "Class", "Order", "Family", "Genus", "Species",
              "Wilderlab_Sp_name", "NZTC_Sp_name", "TaxID", "Report")
out_cols <- intersect(out_cols, names(summary_dt))
setcolOrder(summary_dt, out_cols)

summary_df <- as.data.frame(summary_dt)

# ---------- 12) Ensure final schema & column names expected by ShinyApp.R ----------
msg("Normalising column names to Shiny schema...")

# mapping from current names to Shiny expected names
rename_map <- c(
  "Name" = "Latin Name",
  "CommonName" = "Common Name",
  "ClientSampleID" = "Sample Name",
  "sum_count" = "Total Reads",
  "mean_count" = "Average reads",
  "unique_UID_count" = "Hit number",
  "total_UID" = "Total replicates",
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

# rename only existing columns
existing_old_names <- intersect(names(rename_map), names(summary_df))
if (length(existing_old_names) > 0) {
  names(summary_df)[match(existing_old_names, names(summary_df))] <- unname(rename_map[existing_old_names])
}

# Ensure presence of all expected Shiny columns
expected_cols_for_shiny <- c(
  "Latin Name","Common Name","Sample Name","Total Reads","Average reads","Hit number","Total replicates",
  "Taxonomic Rank","Taxon Group","Threat Status","Threat Category","Date","Threat Document","Collector",
  "DOC Data","Public/Private","Nga Awa","Regional Council","Wilderlab Sp Name","NZTC Sp Name",
  "TaxID","Wilderlab Report","Latitude","Longitude",
  # also keep the raw taxonomic cols used by filters
  "Phylum","Class","Order","Family","Genus","Species"
)
missing_cols <- setdiff(expected_cols_for_shiny, names(summary_df))
if (length(missing_cols) > 0) {
  for (c in missing_cols) summary_df[[c]] <- NA
}

# Coerce types
if ("Date" %in% names(summary_df)) {
  summary_df$Date <- as.IDate(summary_df$Date)
}
if ("Latitude" %in% names(summary_df)) summary_df$Latitude <- suppressWarnings(as.numeric(summary_df$Latitude))
if ("Longitude" %in% names(summary_df)) summary_df$Longitude <- suppressWarnings(as.numeric(summary_df$Longitude))
if ("Total Reads" %in% names(summary_df)) summary_df$`Total Reads` <- suppressWarnings(as.numeric(summary_df$`Total Reads`))
if ("Average reads" %in% names(summary_df)) summary_df$`Average reads` <- suppressWarnings(as.numeric(summary_df$`Average reads`))

# Set final column Order: expected first (for Shiny) then the rest
preferred_Order <- c(intersect(expected_cols_for_shiny, names(summary_df)),
                     setdiff(names(summary_df), expected_cols_for_shiny))
summary_df <- summary_df[, preferred_Order, drop = FALSE]

# ---------- 13) Save outputs (dated and unversioned copy for the Shiny app) ----------
date_suffix <- format(today, "%d%m%y")
out_file <- file.path("Data", paste0("records_", date_suffix, ".RDS"))
msg("Saving dated RDS -> ", out_file)
saveRDS(summary_df, out_file)

unversioned <- file.path("Data", "records.rds")
msg("Saving unversioned RDS -> ", unversioned)
saveRDS(summary_df, unversioned)

msg("Completed. Summary rows: %d", nrow(summary_df))
