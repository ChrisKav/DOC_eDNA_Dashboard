#!/usr/bin/env Rscript
#
# run_nga_awa_by_catchment.R
#
# Driver script that:
#  - loads Data/records.rds (the summary exported by get_API_data.R / previous pipeline)
#  - normalises column names to the ones expected by analyze_community()
#  - iterates over each Nga Awa catchment and runs analyses for several taxonomic communities
#  - writes outputs to Output/<sanitized_NgaAwa>/<CommunityName>/...
#
# Usage:
#  Rscript run_nga_awa_by_catchment.R
#

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(sf)
  library(dplyr)
})

# -------------------------
# Robust script directory detection
# -------------------------
cmdargs <- commandArgs(trailingOnly = FALSE)
file_arg_prefix <- "--file="
script_path <- NA_character_

i <- grep(file_arg_prefix, cmdargs)
if (length(i) > 0) {
  script_path <- sub(file_arg_prefix, "", cmdargs[i][1])
  script_path <- tryCatch(normalizePath(script_path), error = function(e) NA_character_)
}

if (is.na(script_path) || is.null(script_path)) {
  sf <- sys.frames()
  if (!is.null(sf) && length(sf) > 0) {
    for (j in seq_along(sf)) {
      of <- attr(sf[[j]], "ofile")
      if (!is.null(of)) {
        script_path <- tryCatch(normalizePath(of), error = function(e) NA_character_)
        break
      }
    }
  }
}

if (is.na(script_path) || is.null(script_path)) {
  warning("Unable to auto-detect script path; using working directory as script_dir. If files are not found, run this script from its directory or set working directory accordingly.")
  script_dir <- getwd()
} else {
  script_dir <- dirname(script_path)
}

# Source the analysis function (assumes community_analysis_functions.R is in same dir)
func_path <- file.path(script_dir, "community_analysis_functions.R")
if (!file.exists(func_path)) {
  stop("Required file community_analysis_functions.R not found in ", script_dir, ". Place it next to this driver and re-run.")
}
source(func_path)

# Load records.rds
rds_path <- file.path("Data", "records.rds")
if (!file.exists(rds_path)) stop("Data/records.rds not found. Run get_API_data.R first to create it.")
merged <- readRDS(rds_path)
merged <- as_tibble(merged)

# Normalize column names so analyze_community finds expected names
col_renames <- c(
  "Sample Name" = "Sample.Name",
  "Latin Name" = "Latin.Name",
  "DOC Data" = "DOC.Data",
  "Threat Status" = "Threat.Status",
  "Threat Category" = "Threat.Category",
  "Taxon Group" = "Taxon.Group",
  "Taxonomic Rank" = "Taxonomic.Rank",
  "Threat Document" = "Threat.Document",
  "Nga Awa" = "Nga_Awa",
  "Nga_Awa_Catchment" = "Nga_Awa",
  "Regional Council" = "Regional_Council",
  "Public/Private" = "Public.Private",
  "Wilderlab Sp Name" = "Wilderlab_Sp_name",
  "NZTC Sp Name" = "NZTC_Sp_name",
  "UID" = "UID",
  "TaxID" = "TaxID",
  "Latitude" = "Latitude",
  "Longitude" = "Longitude",
  "species" = "species",
  "family" = "family",
  "phylum" = "phylum",
  "ClientSampleID" = "ClientSampleID"
)

for (old in names(col_renames)) {
  if (old %in% names(merged)) {
    names(merged)[names(merged) == old] <- col_renames[[old]]
  }
}

# Ensure Sample.Name exists (fallback to UID or ClientSampleID)
if (!"Sample.Name" %in% names(merged)) {
  if ("UID" %in% names(merged)) merged <- merged %>% mutate(Sample.Name = as.character(UID))
  else if ("ClientSampleID" %in% names(merged)) merged <- merged %>% mutate(Sample.Name = as.character(ClientSampleID))
}

# Ensure numeric coords
if ("Latitude" %in% names(merged)) merged$Latitude <- suppressWarnings(as.numeric(merged$Latitude))
if ("Longitude" %in% names(merged)) merged$Longitude <- suppressWarnings(as.numeric(merged$Longitude))

# -------------------------
# Find and load spatial helper layers from Data/ (robust)
# -------------------------
find_shapefile <- function(root = "Data", patterns = c()) {
  if (!dir.exists(root)) return(character(0))
  shp_files <- list.files(root, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(shp_files) == 0) return(character(0))
  if (length(patterns) == 0) return(shp_files)
  # score candidates by matching any of the provided patterns (case-insensitive)
  scores <- sapply(shp_files, function(p) {
    fn <- basename(p)
    sum(sapply(patterns, function(pt) as.integer(grepl(pt, fn, ignore.case = TRUE))))
  })
  # prefer files with highest score; if tie, return the first
  if (all(scores == 0)) return(character(0))
  shp_files[which.max(scores)]
}

# Try to locate REC2 rivers shapefile
rec2_candidates <- find_shapefile("Data", patterns = c("river_lines", "rec2", "river_line", "river_lines.shp", "river"))
rec2_path <- NA_character_
if (length(rec2_candidates) > 0) {
  rec2_path <- rec2_candidates[1]
} else {
  # fallback to the original expected path (relative)
  rec2_path_fallback <- file.path("REC2_Layers", "River_Lines.shp")
  rec2_path <- if (file.exists(rec2_path_fallback)) rec2_path_fallback else NA_character_
}

# Try to locate Nga Awa shapefile (names can contain unicode 'ā' or ascii variants)
nga_candidates <- find_shapefile("Data", patterns = c("ngaa", "ngaawa", "nga_awa", "ngāawa", "DOC_Ng"))
nga_awa_path <- NA_character_
if (length(nga_candidates) > 0) {
  nga_awa_path <- nga_candidates[1]
} else {
  nga_awa_path_fallback <- file.path("Nga Awa shapefiles", "DOC_NgāAwa_RiverSites_20250122_n14.shp")
  nga_awa_path <- if (file.exists(nga_awa_path_fallback)) nga_awa_path_fallback else NA_character_
}

rec2_rivers <- NULL
nga_awa <- NULL

if (!is.na(rec2_path) && nzchar(rec2_path) && file.exists(rec2_path)) {
  message("Reading REC2 rivers shapefile: ", rec2_path)
  rec2_rivers_try <- tryCatch(st_read(rec2_path, quiet = TRUE), error = function(e) { warning("Failed to read REC2 shapefile: ", e$message); NULL })
  if (!is.null(rec2_rivers_try)) {
    rec2_rivers <- tryCatch(st_transform(rec2_rivers_try, 2193), error = function(e) {
      warning("Failed to transform REC2 to EPSG:2193; proceeding with original CRS: ", e$message)
      rec2_rivers_try
    })
    message("REC2 rivers loaded (rows): ", ifelse(is.null(rec2_rivers), 0, nrow(rec2_rivers)))
  }
} else {
  warning("REC2 rivers shapefile not found. Spatial maps will be omitted for REC2 rivers.")
}

if (!is.na(nga_awa_path) && nzchar(nga_awa_path) && file.exists(nga_awa_path)) {
  message("Reading Nga Awa shapefile: ", nga_awa_path)
  nga_awa_try <- tryCatch(st_read(nga_awa_path, quiet = TRUE), error = function(e) { warning("Failed to read Nga Awa shapefile: ", e$message); NULL })
  if (!is.null(nga_awa_try)) {
    nga_awa <- tryCatch(st_transform(nga_awa_try, 2193), error = function(e) {
      warning("Failed to transform Nga Awa to EPSG:2193; proceeding with original CRS: ", e$message)
      nga_awa_try
    })
    message("Nga Awa shapefile loaded (rows): ", ifelse(is.null(nga_awa), 0, nrow(nga_awa)))
  }
} else {
  warning("Nga Awa shapefile not found. Spatial maps will be omitted or limited.")
}

# Define taxonomic group lists (same as your Script.R)
fish_families <- c("Retropinnidae", "Galaxiidae", "Anguillidae", "Cheimarrichthyidae",
                   "Eleotridae", "Mugilidae", "Tripterygiidae", "Pleuronectidae", 
                   "Microdesmidae", "Gobiidae", "Ictaluridae", "Cyprinidae", 
                   "Poeciliidae", "Salmonidae", "Percidae")

macroinvert_families <- c("Aeshnidae", "Ameletopsidae", "Atyidae", "Austroperlidae", "Calocidae", 
                          "Chathamiidae", "Chordodidae", "Chiltoniidae", "Coenagrionidae", "Coloburiscidae", 
                          "Conoesucidae", "Corduliidae", "Corixidae", "Corophiidae", "Corydalidae", 
                          "Crambidae", "Daphniidae", "Dalyellidae", "Dugesiidae", "Dytiscidae", "Ecnomidae", 
                          "Elmidae", "Ephemeridae", "Eustheniidae", "Gordiidae", "Gripopterygidae", 
                          "Helicophidae", "Helicopsychidae", "Hydraenidae", "Hydrobiidae", "Hydrobiosidae", 
                          "Hydrophilidae", "Hydropsychidae", "Hydroptilidae", "Hyalellidae", "Hymenosomatidae", 
                          "Hyriidae", "Idoteidae", "Kokiriidae", "Latiidae", "Leptoceridae", "Leptophlebiidae", 
                          "Libellulidae", "Limnadiidae", "Lestidae", "Lymnaeidae", "Melanopsidae", "Mesoveliidae", 
                          "Mononchidae", "Mysidae", "Nannochoristidae", "Nereididae", "Nesameletidae", 
                          "Notonectidae", "Notonemouridae", "Oeconesidae", "Oniscigastridae", "Paracalliopiidae", 
                          "Parastacidae", "Petaluridae", "Philopotamidae", "Philorheithridae", "Phoxocephalidae", 
                          "Phreatoicidae", "Planorbidae", "Polycentropodidae", "Pontogeneiidae", "Prorhynchidae", 
                          "Psychomyiidae", "Rallidentidae", "Saldidae", "Siphlaenigmatidae", "Sphaeriidae", 
                          "Sphaeromatidae", "Spionidae", "Talitridae", "Tateidae", "Temnocephalidae", 
                          "Triopsidae", "Typhloplanidae", "Veliidae")

macrophyte_families <- c("Characeae", "Hydatellaceae", "Haloragaceae", "Potamogetonaceae", "Hydrocharitaceae",
                         "Nymphaeaceae", "Cyperaceae", "Juncaceae", "Typhaceae", "Polygonaceae")

diatom_phyla <- c("Bacillariophyta")
ciliate_phyla <- c("Ciliophora")
rotifer_phyla <- c("Rotifera")

# Determine unique Nga Awa catchments in the data
if ("Nga_Awa" %in% names(merged)) {
  catch_col <- "Nga_Awa"
} else if ("Nga Awa" %in% names(merged)) {
  catch_col <- "Nga Awa"
} else if ("Nga_Awa_Catchment" %in% names(merged)) {
  catch_col <- "Nga_Awa_Catchment"
} else if ("Regional_Council" %in% names(merged)) {
  catch_col <- "Regional_Council"
} else {
  stop("No Nga Awa / Regional Council column found in records.rds. Ensure the exported RDS contains Nga Awa catchment info.")
}

catchments <- sort(unique(na.omit(merged[[catch_col]])))
if (length(catchments) == 0) stop("No Nga Awa catchments found in the records file.")

# ----- IGNORE 'not_Nga_Awa' catchment (and variants like 'not Nga Awa') -----
bad_pattern <- "^not[_ ]?nga[_ ]?awa$"
filtered_catchments <- catchments[!grepl(bad_pattern, catchments, ignore.case = TRUE) & nzchar(as.character(catchments))]
removed <- setdiff(catchments, filtered_catchments)
if (length(removed) > 0) {
  message("Ignoring catchment entries matching 'not Nga Awa' pattern: ", paste(removed, collapse = ", "))
}
catchments <- filtered_catchments

if (length(catchments) == 0) {
  stop("After filtering out 'not Nga Awa' catchments, no catchments remain to analyze.")
}

cat("Found", length(catchments), "Nga Awa catchments to process. They are:\n")
print(catchments)

# Helper to sanitize folder names
sanitize_folder <- function(x) {
  nm <- iconv(x, to = "ASCII//TRANSLIT")
  nm <- str_replace_all(nm, "[^A-Za-z0-9_ -]", "_")
  nm <- str_trim(nm)
  nm <- str_replace_all(nm, " ", "_")
  nm <- substr(nm, 1, 120)
  if (nchar(nm) == 0) nm <- "unknown"
  nm
}

# Loop through catchments
for (c in catchments) {
  safe_name <- sanitize_folder(as.character(c))
  base_out <- file.path("Output", safe_name)
  cat("\n========================================\nStarting analyses for catchment:", c, "->", base_out, "\n========================================\n")
  
  # Subset merged to this catchment
  subset_df <- merged %>% filter((!!as.name(catch_col)) == c)
  
  # If no rows, skip
  if (nrow(subset_df) == 0) {
    cat("No records for catchment:", c, " — skipping.\n")
    next
  }
  
  # For each community, create community-specific subset and call analyze_community
  communities <- list(
    Fish = subset_df %>% filter(family %in% fish_families),
    Macroinvertebrate = subset_df %>% filter(family %in% macroinvert_families),
    Macrophytes = subset_df %>% filter(family %in% macrophyte_families),
    Diatoms = subset_df %>% filter(phylum %in% diatom_phyla),
    Ciliates = subset_df %>% filter(phylum %in% ciliate_phyla),
    Rotifers = subset_df %>% filter(phylum %in% rotifer_phyla)
  )
  
  for (comm_name in names(communities)) {
    comm_df <- communities[[comm_name]]
    out_dir <- file.path(base_out, comm_name)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    # call analyze_community wrapped in tryCatch to continue on errors
    tryCatch({
      analyze_community(df = subset_df,
                        community_data = comm_df,
                        community_name = comm_name,
                        output_folder = out_dir,
                        rec2_rivers = rec2_rivers,
                        nga_awa = nga_awa,
                        grouping_distance = 200,
                        buffer_distance = 5000)
    }, error = function(e) {
      warning("Analysis failed for catchment=", c, " community=", comm_name, ": ", e$message)
      writeLines(paste0("ERROR: ", e$message), con = file.path(out_dir, "error.txt"))
    })
  }
  
  cat("Completed catchment:", c, "\n")
} # end catchment loop

cat("\nAll catchment runs complete. Outputs under Output/<NgaAwa>/\n")