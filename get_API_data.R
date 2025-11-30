# Optimised & annotated data preparation script for Wilderlab + NZTCS integration
# - Consolidates repeated code
# - Uses purrr::map_dfr for API loops
# - Uses vectorised joins and safer lookups
# - Improves fuzzy matching with stringdist::amatch
# - Adds light error handling and recommendations to use env vars for secrets
#
# Run: source("data_prep_optimized.R") (ensure required packages installed and data paths / env vars set)

# -------------------------------
# Libraries
# -------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(readxl)
  library(stringr)
  library(stringdist)
  library(sf)
  library(insect)
  library(wilderlab)
})

# -------------------------------
# Configuration / secrets
# -------------------------------
# It's safer to store API keys in environment variables:
# Sys.setenv(WILDER_KEY="...", WILDER_SECRET="...", WILDER_XAPIKEY="...")
key    <- Sys.getenv("WILDER_KEY",    unset = NA_character_)
secret <- Sys.getenv("WILDER_SECRET", unset = NA_character_)
xapikey <- Sys.getenv("WILDER_XAPIKEY", unset = NA_character_)

if (is.na(key) || is.na(secret) || is.na(xapikey)) {
  warning("Wilderlab API keys not set in environment variables (WILDER_KEY / WILDER_SECRET / WILDER_XAPIKEY).",
          " Set them prior to running to avoid hardcoding secrets.")
}

today <- Sys.Date()

# -------------------------------
# Read local NZTCS spreadsheet
# -------------------------------
nztcs <- read_excel("Data/NZTCS.xlsx", sheet = "Exported Data")

# -------------------------------
# Pull Wilderlab data (DOC API + public data)
# - Use map_dfr for records per job (more concise + faster)
# - Use safe wrapper so one failing JobID won't stop the whole run
# -------------------------------
jobs <- get_wilderdata("jobs", key = key, secret = secret, xapikey = xapikey)
samples <- get_wilderdata("samples", key = key, secret = secret, xapikey = xapikey)
taxa <- get_wilderdata("taxa", key = key, secret = secret, xapikey = xapikey)

# records: we request records per job and bind rows
# Using purrr::possibly to continue on errors and return NULL for failed jobs
records <- purrr::map_dfr(jobs$JobID, function(jid) {
  tryCatch({
    get_wilderdata("records", JobID = jid, key = key, secret = secret, xapikey = xapikey)
  }, error = function(e) {
    warning("Failed to fetch records for JobID: ", jid, " - ", e$message)
    return(NULL)
  })
})

# -------------------------------
# Wilderlab public datasets (S3)
# -------------------------------
# Use readr::read_csv which is faster and more consistent than base read.csv
public_samples <- read_csv("http://s3.ap-southeast-2.amazonaws.com/wilderlab.publicdata/samples.csv",
                           show_col_types = FALSE)
public_records <- read_csv("http://s3.ap-southeast-2.amazonaws.com/wilderlab.publicdata/records.csv",
                           show_col_types = FALSE)

# -------------------------------
# Clean and merge samples (DOC + public)
# - Normalize columns in a compact pipeline
# -------------------------------
# Normalize flags and add DOC_Data indicator
DOC_samples <- samples %>%
  mutate(
    MakeDataPublic = ifelse(MakeDataPublic == 1, "Private", "Public"),
    DOC_Data = "Yes",
    CollectionDate = as.Date(CollectionDate),
    Latitude = as.numeric(Latitude),
    Longitude = as.numeric(Longitude)
  )

public_samples <- public_samples %>%
  mutate(
    MakeDataPublic = "Public",
    DOC_Data = "No",
    CollectionDate = as.Date(CollectionDate),
    Latitude = as.numeric(Latitude),
    Longitude = as.numeric(Longitude)
  ) %>%
  # only include public samples that the DOC (private) dataset doesn't already have
  filter(!UID %in% DOC_samples$UID)

all_samples <- bind_rows(DOC_samples, public_samples)

# Keep the useful sample metadata for joins later
merged_samples <- all_samples %>%
  select(SID, UID, JobID, CollectionDate, ClientSampleID,
         Latitude, Longitude, EnvironmentType, TICI, TICIVersion,
         TICINoSeqs, TICIQuantile, MakeDataPublic, DOC_Data, Report)

# -------------------------------
# Merge records (DOC + public)
# -------------------------------
public_records <- public_records %>%
  filter(!UID %in% records$UID)

merged_records <- bind_rows(records, public_records)

# -------------------------------
# Add taxonomy lineages using insect::get_lineage
# - Build tdb from taxa and query once
# -------------------------------
tdb <- taxa %>% select(1:4)
colnames(tdb) <- c("taxID", "parent_taxID", "rank", "name")

# get_lineage can return a list of named character vectors; wrap in try for safety
lineages <- insect::get_lineage(merged_records$TaxID, tdb)

# Extract ranks safely (return NA when not present)
get_rank <- function(lin, rank) {
  if (is.null(lin) || is.na(lin)) return(NA_character_)
  val <- lin[[rank]]
  if (is.null(val)) NA_character_ else val
}

merged_records <- merged_records %>%
  mutate(
    phylum = vapply(lineages, get_rank, CHARACTER(1), rank = "phylum"),
    class  = vapply(lineages, get_rank, CHARACTER(1), rank = "class"),
    order  = vapply(lineages, get_rank, CHARACTER(1), rank = "order"),
    family = vapply(lineages, get_rank, CHARACTER(1), rank = "family"),
    genus  = vapply(lineages, get_rank, CHARACTER(1), rank = "genus"),
    species = vapply(lineages, get_rank, CHARACTER(1), rank = "species")
  )

# Add coordinates and sample metadata to records by looking up in merged_samples
# Use match against merged_samples$UID (not samples which may be incomplete)
merged_records$Latitude <- merged_samples$Latitude[match(merged_records$UID, merged_samples$UID)]
merged_records$Longitude <- merged_samples$Longitude[match(merged_records$UID, merged_samples$UID)]
merged_records$ClientSampleID <- merged_samples$ClientSampleID[match(merged_records$UID, merged_samples$UID)]
merged_records$Report <- merged_samples$Report[match(merged_records$UID, merged_samples$UID)]

# -------------------------------
# Taxonomy cleaning utilities
# -------------------------------
norm_chr <- function(x) ifelse(is.na(x), NA_character_, str_squish(as.character(x)))
cap_case <- function(x) ifelse(is.na(x), NA_character_, str_to_title(x))

class_synonyms <- c("Actinopteri" = "Actinopterygii")
species_synonyms <- c("Galaxias sp. D (Allibone et al., 1996)" = 'Galaxias "species D"')

clean_taxa <- function(df) {
  df %>%
    mutate(
      across(c(phylum, class, order, family, genus, species), ~ norm_chr(.x)),
      class = recode(class, !!!class_synonyms),
      species = recode(species, !!!species_synonyms),
      phylum = cap_case(phylum),
      class  = cap_case(class),
      order  = cap_case(order),
      family = cap_case(family),
      genus  = cap_case(genus)
    )
}

records <- clean_taxa(merged_records)

# -------------------------------
# Prepare NZTCS table
# -------------------------------
nztcs_taxa <- nztcs %>%
  mutate(
    species_nztcs = norm_chr(`Current Species Name`),
    Genus  = cap_case(norm_chr(Genus)),
    Family = cap_case(norm_chr(Family)),
    Order  = cap_case(norm_chr(Order)),
    Class  = cap_case(norm_chr(Class)),
    Phylum = cap_case(norm_chr(Phylum)),
    Status = norm_chr(Status),
    Category = norm_chr(Category),
    BioStatus = norm_chr(`Bio Status`),
    YearAssessed = suppressWarnings(as.integer(`Year Assessed`)),
    ThreatReport = norm_chr(`Report Name`)
  ) %>%
  select(species_nztcs, Genus, Family, Order, Class, Phylum, Status, Category, BioStatus, ThreatReport, YearAssessed)

nztcs_sp <- nztcs_taxa %>%
  arrange(desc(YearAssessed)) %>%
  distinct(species_nztcs, .keep_all = TRUE) %>%
  select(species_nztcs, Status, Category, BioStatus, ThreatReport, Genus, Family, Order, Class, Phylum)

# -------------------------------
# Fuzzy match species names (records -> NZTCS)
# - Normalize strings, remove parentheses/punctuation, use amatch for speed
# -------------------------------
clean_species <- function(x) {
  x <- tolower(x %||% "")
  x <- str_replace_all(x, "\\(.*\\)", "")         # remove parentheses
  x <- str_replace_all(x, "[[:punct:]]", "")      # remove punctuation
  x <- str_replace_all(x, "\\bsp\\.?\\b", "")     # remove 'sp' or 'sp.'
  x <- str_replace_all(x, "\\bspecies\\b", "")    # remove 'species' word
  x <- str_replace_all(x, "\\s+", " ")            # normalize whitespace
  trimws(x)
}

# create cleaned name columns
records <- records %>% mutate(species_clean = clean_species(species))
nztcs_sp <- nztcs_sp %>% mutate(species_clean = clean_species(species_nztcs))

# fuzzy match using stringdist::amatch (fast, returns NA if no match within maxDist)
fuzzy_match_amatch <- function(name, choices, max_dist = 0.12) {
  if (is.na(name) || name == "") return(NA_character_)
  idx <- stringdist::amatch(name, choices, method = "jw", maxDist = max_dist)
  if (is.na(idx)) return(NA_character_) else return(choices[idx])
}

# Map cleaned choices back to NZTCS species_clean
choices <- nztcs_sp$species_clean
records <- records %>%
  mutate(matched_species = vapply(species_clean, fuzzy_match_amatch, CHARACTER(1), choices = choices, max_dist = 0.12))

# -------------------------------
# Join matched NZTCS attributes back to records
# -------------------------------
# Include species_nztcs in the join (so we can name NZTC_Sp_name explicitly)
merged <- records %>%
  left_join(
    nztcs_sp %>% select(species_clean, species_nztcs, Status, Category, BioStatus, ThreatReport, Genus, Family, Order, Class, Phylum),
    by = c("matched_species" = "species_clean")
  ) %>%
  rename(
    Wilderlab_Sp_name = species,
    NZTC_Sp_name = species_nztcs
  )

if (!"Status" %in% names(merged)) merged$Status <- NA_character_

# -------------------------------
# Spatial joins: Nga Awa catchment and Regional Council
# - Read shapefiles once, transform sample geometry once per shapefile
# -------------------------------
# helper to safely read a shapefile
safe_read_sf <- function(path) {
  tryCatch(
    st_read(path, quiet = TRUE),
    error = function(e) { stop("Failed to read shapefile: ", path, " - ", e$message) }
  )
}

NA_polygons <- safe_read_sf("Data/Nga Awa shapefiles/DOC_NgāAwa_RiverSites_20250122_n14.shp")
RC_polygons <- safe_read_sf("Data/Regional Council shapefiles/regional-council-2022-generalised.shp")

# Work only on valid coordinate rows and convert to sf once per polygon
all_samples_valid <- all_samples %>% filter(!is.na(Longitude) & !is.na(Latitude))

# function to join a polygon attribute to samples
join_polygon_attr <- function(samples_df, polygons, attr_col, new_name) {
  sample_sf <- st_as_sf(samples_df, coords = c("Longitude", "Latitude"), crs = 4167, remove = FALSE)
  sample_sf <- st_transform(sample_sf, st_crs(polygons))
  joined <- st_join(sample_sf, polygons, left = TRUE)
  # extract attribute (if missing, set to fallback)
  vals <- joined[[attr_col]]
  vals[is.na(vals)] <- NA_character_
  samples_df[[new_name]] <- vals
  samples_df
}

# Nga Awa: attribute column is "Waterway_N" in the original script
NA_joined <- join_polygon_attr(all_samples_valid, NA_polygons, "Waterway_N", "Nga_Awa_Catchment")
NA_joined$Nga_Awa_Catchment <- ifelse(is.na(NA_joined$Nga_Awa_Catchment), "not Nga Awa", NA_joined$Nga_Awa_Catchment)

# Regional Council: attribute column is "REGC2022_1" in the original script
RC_joined <- join_polygon_attr(all_samples_valid, RC_polygons, "REGC2022_1", "Regional_Council")
RC_joined$Regional_Council <- ifelse(is.na(RC_joined$Regional_Council), "None", RC_joined$Regional_Council)

# Attach these attributes back to all_samples (left join by SID)
all_samples <- all_samples %>%
  left_join(NA_joined %>% select(SID, Nga_Awa_Catchment), by = "SID") %>%
  left_join(RC_joined %>% select(SID, Regional_Council), by = "SID") %>%
  mutate(
    Nga_Awa_Catchment = ifelse(is.na(Nga_Awa_Catchment), "not Nga Awa", Nga_Awa_Catchment),
    Regional_Council = ifelse(is.na(Regional_Council), "None", Regional_Council)
  )

# -------------------------------
# Merge records with sample metadata and select columns
# -------------------------------
all_merged <- merged %>%
  left_join(all_samples, by = "UID") %>%
  mutate(
    Latitude = coalesce(as.numeric(Latitude.x), as.numeric(Latitude.y)),
    Longitude = coalesce(as.numeric(Longitude.x), as.numeric(Longitude.y)),
    ClientSampleID = coalesce(ClientSampleID.x, ClientSampleID.y)
  ) %>%
  select(-starts_with("Latitude."), -starts_with("Longitude."), -starts_with("ClientSampleID."))

# Select a consistent set of columns
selected_df <- all_merged %>%
  select(UID, TaxID, Rank, Name, CommonName, Group, Count, Latitude, Longitude, ClientSampleID, CollectionDate,
         Status, Category, ThreatReport, CollectedBy, DOC_Data, MakeDataPublic, Nga_Awa_Catchment, Regional_Council,
         phylum, class, order, family, genus, Wilderlab_Sp_name, NZTC_Sp_name, Report)

# -------------------------------
# Summarise by Wilderlab Report and TaxID
# -------------------------------
uid_counts <- selected_df %>%
  group_by(Report) %>%
  summarise(total_UID = n_distinct(UID), .groups = "drop")

summary_df <- selected_df %>%
  group_by(Report, TaxID) %>%
  summarise(
    unique_UID_count = n_distinct(UID),
    UID_list = paste(unique(UID), collapse = "-"),
    sum_count = sum(Count, na.rm = TRUE),
    ClientSampleID = first(ClientSampleID),
    across(c(Rank, Name, CommonName, Group, Latitude, Longitude, CollectionDate,
             Status, Category, ThreatReport, CollectedBy, DOC_Data, MakeDataPublic, Nga_Awa_Catchment, Regional_Council,
             phylum, class, order, family, genus, Wilderlab_Sp_name, NZTC_Sp_name), first),
    .groups = "drop"
  ) %>%
  left_join(uid_counts, by = "Report") %>%
  mutate(mean_count = sum_count / total_UID)

summary_df <- summary_df %>% 
  select(Name, CommonName, ClientSampleID, sum_count, mean_count, unique_UID_count, total_UID, Rank, Group, Status,
         Category, Latitude, Longitude, CollectionDate, ThreatReport, CollectedBy, DOC_Data, MakeDataPublic,
         Nga_Awa_Catchment, Regional_Council, phylum, class, order, family, genus, Wilderlab_Sp_name, NZTC_Sp_name,
         TaxID, Report)

# Friendly column names
summary_df <- summary_df %>%
  rename(
    `Latin Name` = Name,
    `Common Name` = CommonName,
    `Sample Name` = ClientSampleID,
    `Total Reads` = sum_count,
    `Average reads` = mean_count,
    `Hit number` = unique_UID_count,
    `Total hits` = total_UID,
    `Taxonomic Rank` = Rank,
    `Taxon Group` = Group,
    `Threat Status` = Status,
    `Threat Category` = Category,
    `Date` = CollectionDate,
    `Threat Document` = ThreatReport,
    `Collector` = CollectedBy,
    `DOC Data` = DOC_Data,
    `Public/Private` = MakeDataPublic,
    `Nga Awa` = Nga_Awa_Catchment,
    `Regional Council` = Regional_Council,
    `Wilderlab Sp Name` = Wilderlab_Sp_name,
    `NZTC Sp Name` = NZTC_Sp_name,
    `TaxID` = TaxID,
    `Wilderlab Report` = Report
  )

# -------------------------------
# Save for shiny app or downstream use
# -------------------------------
saveRDS(summary_df, "Data/records.rds")
message("Saved Data/records.rds (rows: ", nrow(summary_df), ")")