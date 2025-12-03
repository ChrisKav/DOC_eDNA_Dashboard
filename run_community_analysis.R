# ============================================================
# BATCH ANALYSIS FOR ALL NGA AWA CATCHMENTS
# Runs community analysis for each catchment separately
# ============================================================

# -----------------------------
# 0. Load Required Libraries
# -----------------------------
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

# -----------------------------
# 1. SOURCE THE ANALYSIS FUNCTION
# -----------------------------
# Make sure community_analysis_functions.R contains ONLY the analyze_community function
# Remove any execution code from the bottom of community_analysis_functions.R before sourcing
source("community_analysis_functions.R")

# -----------------------------
# 2. LOAD SPATIAL DATA
# -----------------------------
cat("\n########## LOADING SPATIAL DATA ##########\n")
rec2_rivers <- st_read("Data/REC2_Layers/River_Lines.shp") %>% st_transform(2193)
nga_awa <- st_read("Data/Nga Awa shapefiles/DOC_NgÄAwa_RiverSites_20250122_n14.shp") %>% st_transform(2193)
cat("Spatial layers loaded successfully\n\n")

# -----------------------------
# 3. LOAD AND PREPARE DATA
# -----------------------------
cat("########## LOADING RECORDS DATA ##########\n")
df_all <- readRDS("Data/records.rds")

cat("\nColumn names in records.rds:\n")
print(colnames(df_all))
cat("\n")

# Convert spaces to dots if needed (to match community_analysis_functions.R format)
if("Sample Name" %in% colnames(df_all)) {
  colnames(df_all) <- gsub(" ", ".", colnames(df_all))
  cat("Converted spaces to dots in column names\n\n")
}

# Check required columns
required_cols <- c("Sample.Name", "Nga.Awa", "species", "family", "phylum", 
                   "Latin.Name", "Threat.Status", "Date", "DOC.Data",
                   "Latitude", "Longitude")

missing_cols <- setdiff(required_cols, colnames(df_all))
if(length(missing_cols) > 0) {
  cat("\nERROR - Missing required columns:\n")
  cat("  Missing:", paste(missing_cols, collapse = ", "), "\n")
  cat("  Available:", paste(colnames(df_all), collapse = ", "), "\n\n")
  stop("Cannot proceed")
}

cat("All required columns present\n\n")

# -----------------------------
# 4. DEFINE COMMUNITY FILTERS
# -----------------------------
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

# -----------------------------
# 5. IDENTIFY CATCHMENTS
# -----------------------------
catchments <- df_all %>%
  filter(!is.na(Nga.Awa) & Nga.Awa != "not Nga Awa" & Nga.Awa != "") %>%
  pull(Nga.Awa) %>%
  unique() %>%
  sort()

cat(paste0("Found ", length(catchments), " Nga Awa catchments:\n"))
for(i in seq_along(catchments)) {
  cat(sprintf("  %2d. %s\n", i, catchments[i]))
}
cat("\n")

# Diagnostic check
cat("Diagnostic check - data availability per catchment:\n")
cat(sprintf("%-30s %8s %8s %8s\n", "Catchment", "Total", "Fish", "Macroinv"))
cat(rep("-", 60), "\n", sep = "")

for(catchment in catchments) {
  n_total <- sum(df_all$Nga.Awa == catchment, na.rm = TRUE)
  n_fish <- sum(df_all$Nga.Awa == catchment & df_all$family %in% fish_families, na.rm = TRUE)
  n_macro <- sum(df_all$Nga.Awa == catchment & df_all$family %in% macroinvert_families, na.rm = TRUE)
  
  cat(sprintf("%-30s %8d %8d %8d\n", catchment, n_total, n_fish, n_macro))
}
cat("\n")

# -----------------------------
# 6. ANALYSIS FUNCTION
# -----------------------------
run_catchment_analysis <- function(catchment_name, df_all, 
                                   rec2_rivers, nga_awa,
                                   fish_families, macroinvert_families, 
                                   macrophyte_families, diatom_phyla,
                                   ciliate_phyla, rotifer_phyla) {
  
  cat("\n")
  cat(rep("#", 80), "\n", sep = "")
  cat(sprintf("### %s ###\n", toupper(catchment_name)))
  cat(rep("#", 80), "\n\n", sep = "")
  
  # Filter data
  df_catchment <- df_all %>% filter(Nga.Awa == catchment_name)
  
  if(nrow(df_catchment) == 0) {
    cat("WARNING: No data - skipping\n")
    return(NULL)
  }
  
  cat(sprintf("Processing %d records\n\n", nrow(df_catchment)))
  
  # Create output folder
  catchment_folder <- file.path("Output", gsub("[^A-Za-z0-9_]", "_", catchment_name))
  if(!dir.exists(catchment_folder)) dir.create(catchment_folder, recursive = TRUE)
  
  results <- list()
  
  # FISH
  tryCatch({
    fish_data <- df_catchment %>% filter(family %in% fish_families)
    cat(sprintf("Fish: %d records\n", nrow(fish_data)))
    if(nrow(fish_data) > 0) {
      results$fish <- analyze_community(
        df = df_catchment, community_data = fish_data, community_name = "Fish",
        output_folder = file.path(catchment_folder, "Fish"),
        rec2_rivers = rec2_rivers, nga_awa = nga_awa,
        grouping_distance = 200, buffer_distance = 5000
      )
      cat("  -> Complete\n")
    } else {
      cat("  -> Skipped (no data)\n")
    }
  }, error = function(e) {
    cat("  -> ERROR:", e$message, "\n")
  })
  
  # MACROINVERTEBRATES
  tryCatch({
    macro_data <- df_catchment %>% filter(family %in% macroinvert_families)
    cat(sprintf("Macroinvertebrates: %d records\n", nrow(macro_data)))
    if(nrow(macro_data) > 0) {
      results$macroinvert <- analyze_community(
        df = df_catchment, community_data = macro_data, community_name = "Macroinvertebrate",
        output_folder = file.path(catchment_folder, "Macroinvertebrates"),
        rec2_rivers = rec2_rivers, nga_awa = nga_awa,
        grouping_distance = 200, buffer_distance = 5000
      )
      cat("  -> Complete\n")
    } else {
      cat("  -> Skipped (no data)\n")
    }
  }, error = function(e) {
    cat("  -> ERROR:", e$message, "\n")
  })
  
  # MACROPHYTES
  tryCatch({
    macrophyte_data <- df_catchment %>% filter(family %in% macrophyte_families)
    cat(sprintf("Macrophytes: %d records\n", nrow(macrophyte_data)))
    if(nrow(macrophyte_data) > 0) {
      results$macrophyte <- analyze_community(
        df = df_catchment, community_data = macrophyte_data, community_name = "Macrophytes",
        output_folder = file.path(catchment_folder, "Macrophytes"),
        rec2_rivers = rec2_rivers, nga_awa = nga_awa,
        grouping_distance = 200, buffer_distance = 5000
      )
      cat("  -> Complete\n")
    } else {
      cat("  -> Skipped (no data)\n")
    }
  }, error = function(e) {
    cat("  -> ERROR:", e$message, "\n")
  })
  
  # DIATOMS
  tryCatch({
    diatom_data <- df_catchment %>% filter(phylum %in% diatom_phyla)
    cat(sprintf("Diatoms: %d records\n", nrow(diatom_data)))
    if(nrow(diatom_data) > 0) {
      results$diatom <- analyze_community(
        df = df_catchment, community_data = diatom_data, community_name = "Diatoms",
        output_folder = file.path(catchment_folder, "Diatoms"),
        rec2_rivers = rec2_rivers, nga_awa = nga_awa,
        grouping_distance = 200, buffer_distance = 5000
      )
      cat("  -> Complete\n")
    } else {
      cat("  -> Skipped (no data)\n")
    }
  }, error = function(e) {
    cat("  -> ERROR:", e$message, "\n")
  })
  
  # CILIATES
  tryCatch({
    ciliate_data <- df_catchment %>% filter(phylum %in% ciliate_phyla)
    cat(sprintf("Ciliates: %d records\n", nrow(ciliate_data)))
    if(nrow(ciliate_data) > 0) {
      results$ciliate <- analyze_community(
        df = df_catchment, community_data = ciliate_data, community_name = "Ciliates",
        output_folder = file.path(catchment_folder, "Ciliates"),
        rec2_rivers = rec2_rivers, nga_awa = nga_awa,
        grouping_distance = 200, buffer_distance = 5000
      )
      cat("  -> Complete\n")
    } else {
      cat("  -> Skipped (no data)\n")
    }
  }, error = function(e) {
    cat("  -> ERROR:", e$message, "\n")
  })
  
  # ROTIFERS
  tryCatch({
    rotifer_data <- df_catchment %>% filter(phylum %in% rotifer_phyla)
    cat(sprintf("Rotifers: %d records\n", nrow(rotifer_data)))
    if(nrow(rotifer_data) > 0) {
      results$rotifer <- analyze_community(
        df = df_catchment, community_data = rotifer_data, community_name = "Rotifers",
        output_folder = file.path(catchment_folder, "Rotifers"),
        rec2_rivers = rec2_rivers, nga_awa = nga_awa,
        grouping_distance = 200, buffer_distance = 5000
      )
      cat("  -> Complete\n")
    } else {
      cat("  -> Skipped (no data)\n")
    }
  }, error = function(e) {
    cat("  -> ERROR:", e$message, "\n")
  })
  
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")
  cat(sprintf("COMPLETED: %s\n", toupper(catchment_name)))
  cat(rep("=", 80), "\n\n", sep = "")
  
  return(results)
}

# -----------------------------
# 7. RUN ALL ANALYSES
# -----------------------------
cat("\n########## STARTING BATCH ANALYSIS ##########\n")
cat(sprintf("Processing %d catchments\n\n", length(catchments)))

all_results <- list()

for(i in seq_along(catchments)) {
  catchment <- catchments[i]
  cat(sprintf("\n[%d/%d] %s\n", i, length(catchments), catchment))
  
  all_results[[catchment]] <- run_catchment_analysis(
    catchment_name = catchment,
    df_all = df_all,
    rec2_rivers = rec2_rivers,
    nga_awa = nga_awa,
    fish_families = fish_families,
    macroinvert_families = macroinvert_families,
    macrophyte_families = macrophyte_families,
    diatom_phyla = diatom_phyla,
    ciliate_phyla = ciliate_phyla,
    rotifer_phyla = rotifer_phyla
  )
}

# -----------------------------
# 8. SUMMARY REPORT
# -----------------------------
cat("\n")
cat(rep("#", 80), "\n", sep = "")
cat("### BATCH ANALYSIS COMPLETE ###\n")
cat(rep("#", 80), "\n\n", sep = "")

cat("Summary of analyses:\n\n")
for(catchment in names(all_results)) {
  cat(sprintf("  %s:\n", catchment))
  if(!is.null(all_results[[catchment]])) {
    communities <- names(all_results[[catchment]])
    if(length(communities) > 0) {
      cat(sprintf("    - %s\n", paste(communities, collapse = ", ")))
    } else {
      cat("    - No communities analyzed\n")
    }
  } else {
    cat("    - Failed or no data\n")
  }
}

cat("\nOutputs saved to Output/ folder\n")
cat("Results object saved to Output/all_catchment_results.RDS\n")

saveRDS(all_results, file = "Output/all_catchment_results.RDS")

cat("\n########## DONE! ##########\n")