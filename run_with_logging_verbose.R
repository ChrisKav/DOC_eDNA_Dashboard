# Safe batch runner for community analyses with detailed error logging.
# Place in repo root (next to run_community_analysis.R and community_analysis_functions.R)
# Run:
#   Rscript run_community_analysis_safe.R
# or from interactive R:
#   source("run_community_analysis_safe.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ggplot2)
  library(sf)
  library(ggrepel)
  library(ggnewscale)
  library(cluster)
  library(factoextra)
  library(ggspatial)
  library(patchwork)
  library(flextable)
  library(officer)
  library(indicspecies)
  library(pheatmap)
})

# create logging dirs
logdir <- file.path(".", "logs")
problems_dir <- file.path(logdir, "problems")
dir.create(problems_dir, recursive = TRUE, showWarnings = FALSE)

make_ts <- function() format(Sys.time(), "%Y%m%dT%H%M%S")
safe_name <- function(x) gsub("[^A-Za-z0-9_\\-\\.]", "_", x)

# Logging helper - saves rds and txt and appends summary csv
log_problem_full <- function(catchment = NA_character_, community = NA_character_, df_catchment = NULL, community_data = NULL, err = NULL, context = NULL) {
  ts <- make_ts()
  c_nm <- safe_name(ifelse(is.na(catchment) || catchment == "", "unknown_catchment", as.character(catchment)))
  comm_nm <- safe_name(ifelse(is.na(community) || community == "", "unknown_community", as.character(community)))
  base <- file.path(problems_dir, paste0("problem_", c_nm, "_", comm_nm, "_", ts))
  # capture traceback (best effort)
  tb <- tryCatch({
    # when called inside error handler, traceback() prints the last stack
    paste(capture.output(traceback(max.lines = 200)), collapse = "\n")
  }, error = function(e) paste0("traceback capture failed: ", e$message))
  saveRDS(list(
    time = Sys.time(),
    catchment = catchment,
    community = community,
    df_catchment = df_catchment,
    community_data = community_data,
    error_condition = err,
    error_message = if (!is.null(err)) conditionMessage(err) else NA_character_,
    error_print = if (!is.null(err)) paste(capture.output(print(err)), collapse = "\n") else NA_character_,
    traceback = tb,
    context = context,
    session = utils::sessionInfo()
  ), file = paste0(base, ".rds"))
  # human readable text
  txt_lines <- c(
    paste0("time: ", Sys.time()),
    paste0("catchment: ", c_nm),
    paste0("community: ", comm_nm),
    "",
    "error message:",
    if (!is.null(err)) conditionMessage(err) else "NULL",
    "",
    "error (printed):",
    if (!is.null(err)) paste(capture.output(print(err)), collapse = "\n") else "NULL",
    "",
    "traceback:",
    tb,
    "",
    "context:",
    if (!is.null(context)) paste(capture.output(str(context)), collapse = "\n") else "NULL",
    "",
    "sessionInfo:",
    paste(capture.output(utils::sessionInfo()), collapse = "\n")
  )
  writeLines(txt_lines, con = paste0(base, ".txt"))
  # append CSV summary row
  summary_row <- data.frame(
    time = Sys.time(),
    catchment = c_nm,
    community = comm_nm,
    message = if (!is.null(err)) conditionMessage(err) else NA_character_,
    rds = paste0(basename(base), ".rds"),
    txt = paste0(basename(base), ".txt"),
    stringsAsFactors = FALSE
  )
  csv_file <- file.path(logdir, "error_summary.csv")
  write.table(summary_row, file = csv_file, sep = ",", row.names = FALSE,
              col.names = !file.exists(csv_file), append = TRUE)
  invisible(list(rds = paste0(base, ".rds"), txt = paste0(base, ".txt"), summary = summary_row))
}

# Safe wrapper around analyze_community
safe_analyze <- function(df, community_data, community_name, output_folder, rec2_rivers, nga_awa, grouping_distance = 200, buffer_distance = 5000) {
  tryCatch({
    res <- analyze_community(df = df,
                             community_data = community_data,
                             community_name = community_name,
                             output_folder = output_folder,
                             rec2_rivers = rec2_rivers,
                             nga_awa = nga_awa,
                             grouping_distance = grouping_distance,
                             buffer_distance = buffer_distance)
    # return result if successful
    res
  }, error = function(e) {
    # Save full debug info and rethrow (so the batch keeps going)
    log_problem_full(
      catchment = if ("Nga.Awa" %in% names(df)) unique(na.omit(as.character(df$Nga.Awa))) else NA_character_,
      community = community_name,
      df_catchment = df,
      community_data = community_data,
      err = e,
      context = list(call = match.call())
    )
    # print useful info to console too
    cat("ERROR in analyze_community - community:", community_name, "\n")
    cat("Error message:", conditionMessage(e), "\n\n")
    cat("Full printed condition:\n")
    print(e)
    cat("\nTraceback (captured below):\n")
    try(traceback(), silent = TRUE)
    # return NULL so caller knows it failed
    NULL
  })
}

# --- Now replicate the parts of run_community_analysis.R that detect catchments and run them,
#     but using safe_analyze so we get detailed logs.

# Source existing analyze function
source("community_analysis_functions.R")  # must export analyze_community

# Load spatial data (same as original)
cat("\n########## LOADING SPATIAL DATA ##########\n")
rec2_rivers <- st_read("Data/REC2_Layers/River_Lines.shp") %>% st_transform(2193)
nga_awa <- st_read("Data/Nga Awa shapefiles/DOC_NgāAwa_RiverSites_20250122_n14.shp") %>% st_transform(2193)
cat("Spatial layers loaded successfully\n\n")

# Load records
cat("########## LOADING RECORDS DATA ##########\n")
df_all <- readRDS("Data/records_Waipoua.rds")

# Normalize column names to match
if("Sample Name" %in% colnames(df_all)) {
  colnames(df_all) <- gsub(" ", ".", colnames(df_all))
}

required_cols <- c("Sample.Name", "Nga.Awa", "species", "family", "phylum", 
                   "Latin.Name", "Threat.Status", "Date", "DOC.Data",
                   "Latitude", "Longitude")
missing_cols <- setdiff(required_cols, colnames(df_all))
if(length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

# identify catchments (same filter)
catchments <- df_all %>%
  filter(!is.na(Nga.Awa) & Nga.Awa != "not Nga Awa" & Nga.Awa != "") %>%
  pull(Nga.Awa) %>% unique() %>% sort()

cat(paste0("Found ", length(catchments), " Nga Awa catchments:\n"))
print(catchments)
cat("\n")

# Community filters (can reuse those defined in your functions file if present there, otherwise re-declare)
# For safety, try to reuse variables from the sourced file; fallback values below if not present.
if (!exists("fish_families")) fish_families <- c("Retropinnidae", "Galaxiidae", "Anguillidae", "Cheimarrichthyidae",
                                                 "Eleotridae", "Mugilidae", "Tripterygiidae", "Pleuronectidae", 
                                                 "Microdesmidae", "Gobiidae", "Ictaluridae", "Cyprinidae", 
                                                 "Poeciliidae", "Salmonidae", "Percidae")
if (!exists("macroinvert_families")) macroinvert_families <- c("Aeshnidae","Ameletopsidae","Atyidae")  # short fallback
if (!exists("macrophyte_families")) macrophyte_families <- c("Characeae","Hydatellaceae")
if (!exists("diatom_phyla")) diatom_phyla <- c("Bacillariophyta")
if (!exists("ciliate_phyla")) ciliate_phyla <- c("Ciliophora")
if (!exists("rotifer_phyla")) rotifer_phyla <- c("Rotifera")

# Ensure Output dir
if(!dir.exists("Output")) dir.create("Output", recursive = TRUE)

all_results <- list()

for (catchment in catchments) {
  cat("\n========================================\n")
  cat("Processing catchment:", catchment, "\n")
  cat("========================================\n\n")
  
  df_catchment <- df_all %>% filter(Nga.Awa == catchment)
  if (nrow(df_catchment) == 0) {
    cat("No records, skipping\n")
    next
  }
  
  catchment_folder <- file.path("Output", safe_name(catchment))
  if(!dir.exists(catchment_folder)) dir.create(catchment_folder, recursive = TRUE)
  
  results <- list()
  
  # -- Fish
  fish_data <- df_catchment %>% filter(family %in% fish_families)
  cat("Fish records:", nrow(fish_data), "\n")
  if (nrow(fish_data) > 0) {
    results$fish <- safe_analyze(df = df_catchment, community_data = fish_data, community_name = "Fish",
                                 output_folder = file.path(catchment_folder, "Fish"),
                                 rec2_rivers = rec2_rivers, nga_awa = nga_awa)
    if (is.null(results$fish)) cat("Fish analysis failed and was logged\n") else cat("Fish analysis completed\n")
  } else {
    cat("Fish: no data\n")
  }
  
  # -- Macroinvertebrates
  macro_data <- df_catchment %>% filter(family %in% macroinvert_families)
  cat("Macroinvertebrate records:", nrow(macro_data), "\n")
  if (nrow(macro_data) > 0) {
    results$macroinvert <- safe_analyze(df = df_catchment, community_data = macro_data, community_name = "Macroinvertebrate",
                                        output_folder = file.path(catchment_folder, "Macroinvertebrates"),
                                        rec2_rivers = rec2_rivers, nga_awa = nga_awa)
    if (is.null(results$macroinvert)) cat("Macroinvertebrate analysis failed and was logged\n") else cat("Macroinvertebrate analysis completed\n")
  } else {
    cat("Macroinvertebrates: no data\n")
  }
  
  # -- Macrophytes
  mac_data <- df_catchment %>% filter(family %in% macrophyte_families)
  cat("Macrophyte records:", nrow(mac_data), "\n")
  if (nrow(mac_data) > 0) {
    results$macrophyte <- safe_analyze(df = df_catchment, community_data = mac_data, community_name = "Macrophytes",
                                       output_folder = file.path(catchment_folder, "Macrophytes"),
                                       rec2_rivers = rec2_rivers, nga_awa = nga_awa)
    if (is.null(results$macrophyte)) cat("Macrophyte analysis failed and was logged\n") else cat("Macrophyte analysis completed\n")
  } else {
    cat("Macrophytes: no data\n")
  }
  
  # -- Diatoms
  diatom_data <- df_catchment %>% filter(phylum %in% diatom_phyla)
  cat("Diatom records:", nrow(diatom_data), "\n")
  if (nrow(diatom_data) > 0) {
    results$diatom <- safe_analyze(df = df_catchment, community_data = diatom_data, community_name = "Diatoms",
                                   output_folder = file.path(catchment_folder, "Diatoms"),
                                   rec2_rivers = rec2_rivers, nga_awa = nga_awa)
    if (is.null(results$diatom)) cat("Diatom analysis failed and was logged\n") else cat("Diatom analysis completed\n")
  } else {
    cat("Diatoms: no data\n")
  }
  
  # -- Ciliates
  ciliate_data <- df_catchment %>% filter(phylum %in% ciliate_phyla)
  cat("Ciliate records:", nrow(ciliate_data), "\n")
  if (nrow(ciliate_data) > 0) {
    results$ciliate <- safe_analyze(df = df_catchment, community_data = ciliate_data, community_name = "Ciliates",
                                    output_folder = file.path(catchment_folder, "Ciliates"),
                                    rec2_rivers = rec2_rivers, nga_awa = nga_awa)
    if (is.null(results$ciliate)) cat("Ciliate analysis failed and was logged\n") else cat("Ciliate analysis completed\n")
  } else {
    cat("Ciliates: no data\n")
  }
  
  # -- Rotifers
  rotifer_data <- df_catchment %>% filter(phylum %in% rotifer_phyla)
  cat("Rotifer records:", nrow(rotifer_data), "\n")
  if (nrow(rotifer_data) > 0) {
    results$rotifer <- safe_analyze(df = df_catchment, community_data = rotifer_data, community_name = "Rotifers",
                                    output_folder = file.path(catchment_folder, "Rotifers"),
                                    rec2_rivers = rec2_rivers, nga_awa = nga_awa)
    if (is.null(results$rotifer)) cat("Rotifer analysis failed and was logged\n") else cat("Rotifer analysis completed\n")
  } else {
    cat("Rotifers: no data\n")
  }
  
  all_results[[catchment]] <- results
  # Save intermediate results
  saveRDS(all_results, file = file.path("Output", "all_catchment_results_partial.RDS"))
}

# Final save
saveRDS(all_results, file = "Output/all_catchment_results.RDS")
cat("Safe batch run complete. Check logs/ and logs/problems for errors.\n")