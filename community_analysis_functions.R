#!/usr/bin/env Rscript
#
# community_analysis_functions.R
#
# Contains analyze_community() adapted from your Script.R but written so it can
# be sourced and reused by a driver that runs analyses per Nga Awa catchment.
#
# The function expects a data.frame / tibble with columns:
#   Sample.Name, Latitude, Longitude, species, family, phylum, Latin.Name,
#   Threat.Status, DOC.Data, and any other metadata used in outputs.
#
# It writes a set of output files into the output_folder you provide.
#
# NOTE: This file only defines the function and supporting helpers. It does NOT
# run anything by itself. Source it from a driver script.
#

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

# Small helper to safely create folders
safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# Main analysis function (adapted, but behaviour kept similar to your Script.R)
analyze_community <- function(df,
                              community_data,
                              community_name,
                              output_folder,
                              rec2_rivers = NULL,
                              nga_awa = NULL,
                              grouping_distance = 200,
                              buffer_distance = 5000) {
  # df: full dataset (used for some summaries)
  # community_data: filtered dataset for the community (rows to analyse)
  # output_folder: path where outputs will be saved (function will create it)
  
  safe_dir_create(output_folder)
  cat("Running", community_name, "analysis ->", output_folder, "\n")
  
  # Ensure key columns exist
  required_cols <- c("Sample.Name", "Latitude", "Longitude", "species", "family", "phylum", "Latin.Name", "Threat.Status", "DOC.Data")
  for (c in required_cols) {
    if (!c %in% names(community_data)) community_data[[c]] <- NA
  }
  
  # STEP 1: Prepare community presence-absence matrix per Sample.Name
  site_coords <- community_data %>%
    select(Sample.Name, Latitude, Longitude) %>%
    distinct()
  
  # Build presence (1/0) matrix
  community_matrix <- community_data %>%
    filter(!is.na(Sample.Name)) %>%
    group_by(Sample.Name, species) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(presence = 1) %>%
    pivot_wider(names_from = species, values_from = presence, values_fill = 0)
  
  # Join coordinates (if available)
  community_matrix <- community_matrix %>%
    left_join(site_coords, by = "Sample.Name")
  
  # Remove any column named "NA" introduced by pivoting
  if ("NA" %in% names(community_matrix)) community_matrix <- community_matrix %>% select(-`NA`)
  
  # Replace NA numeric species cells with 0
  species_cols <- setdiff(names(community_matrix), c("Sample.Name", "Latitude", "Longitude"))
  community_matrix <- community_matrix %>%
    mutate(across(all_of(species_cols), ~replace_na(.x, 0)))
  
  # If there are no sites with any species, write an empty indicator and exit gracefully
  if (nrow(community_matrix) == 0 || length(species_cols) == 0) {
    cat("No data available for", community_name, "in this subset. Writing empty placeholder and returning.\n")
    saveRDS(list(message = paste("No data for", community_name)), file.path(output_folder, "no_data.rds"))
    return(invisible(NULL))
  }
  
  # Create spatial object if coords present
  coords_present <- all(c("Longitude", "Latitude") %in% names(community_matrix)) && any(!is.na(community_matrix$Longitude) & !is.na(community_matrix$Latitude))
  if (coords_present) {
    community_sf <- community_matrix %>%
      filter(!is.na(Longitude) & !is.na(Latitude)) %>%
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
      st_transform(2193)
  } else {
    # create an empty sf to avoid downstream failures
    community_sf <- st_sf(geometry = st_sfc(), crs = 2193)
  }
  
  # Group sites by distance (if we have at least one point)
  if (nrow(community_sf) > 0) {
    groups <- st_is_within_distance(community_sf, community_sf, dist = grouping_distance)
    community_sf$group_id <- sapply(groups, function(x) min(x))
  } else {
    community_sf$group_id <- integer(0)
  }
  
  # Attach group_id back to community_matrix by Sample.Name (keep first if multiple)
  group_df <- community_sf %>%
    st_drop_geometry() %>%
    select(Sample.Name, group_id)
  
  community_matrix <- community_matrix %>%
    left_join(group_df, by = "Sample.Name")
  
  # assign a group id to any samples without group_id
  if (!"group_id" %in% names(community_matrix)) community_matrix$group_id <- seq_len(nrow(community_matrix))
  
  # Aggregate by group
  aggregated <- community_matrix %>%
    group_by(group_id) %>%
    summarise(
      Sample.Name = first(Sample.Name),
      Latitude = first(Latitude),
      Longitude = first(Longitude),
      across(where(is.numeric) & !c(Latitude, Longitude), max),
      .groups = "drop"
    ) %>%
    mutate(
      Representative_Site = paste0("Site_", group_id),
      Richness = rowSums(select(., -group_id, -Sample.Name, -Latitude, -Longitude))
    )
  
  filtered_matrix <- aggregated %>% filter(Richness > 0)
  
  # species matrix for multivariate analyses (rows = sites, cols = species)
  species_matrix <- filtered_matrix %>%
    select(-group_id, -Sample.Name, -Representative_Site, -Latitude, -Longitude, -Richness) %>%
    as.data.frame()
  if (ncol(species_matrix) == 0) {
    cat("No species columns remain after aggregation for", community_name, "\n")
    saveRDS(list(message = paste("No species for", community_name)), file.path(output_folder, "no_species.rds"))
    return(invisible(NULL))
  }
  rownames(species_matrix) <- filtered_matrix$Representative_Site
  
  # Save a quick CSV of aggregated sites for inspection
  write_csv(filtered_matrix %>% st_drop_geometry(), file.path(output_folder, "aggregated_sites.csv"))
  
  # STEP 2: Sampling effort summary (by year and DOC vs Non-DOC)
  sampling_summary <- df %>%
    mutate(Year = suppressWarnings(as.integer(format(as.Date(Date), "%Y")))) %>%
    distinct(Sample.Name, Year, DOC.Data) %>%
    mutate(DOC.Data = ifelse(is.na(DOC.Data) | DOC.Data == "", "Unknown", DOC.Data)) %>%
    group_by(Year, DOC.Data) %>%
    summarise(n_samples = n(), .groups = "drop") %>%
    mutate(DOC.Data = factor(DOC.Data))
  
  sampling_totals <- sampling_summary %>%
    group_by(Year) %>%
    summarise(total = sum(n_samples), .groups = "drop")
  
  sampling_table <- sampling_summary %>%
    pivot_wider(names_from = DOC.Data, values_from = n_samples, values_fill = 0) %>%
    left_join(sampling_totals, by = "Year") %>%
    select(Year, everything())
  
  # save sampling table and figure
  write_csv(sampling_table, file.path(output_folder, "01_sampling_summary_by_year.csv"))
  
  # Plot (attempt; if plotting fails continue)
  tryCatch({
    sampling_plot <- ggplot(sampling_summary, aes(x = factor(Year), y = n_samples, fill = DOC.Data)) +
      geom_col(color = "black", width = 0.7, linewidth = 0.5) +
      geom_text(data = sampling_totals, aes(x = factor(Year), y = total, label = total, fill = NULL),
                vjust = -0.5, size = 4, fontface = "bold") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      theme_bw(base_size = 12) +
      labs(title = paste(community_name, "Sampling Effort Over Time"), x = "Year", y = "Number of Samples")
    ggsave(file.path(output_folder, "01_sampling_effort_by_year.png"), sampling_plot, width = 8, height = 6, dpi = 300, bg = "white")
  }, error = function(e) {
    warning("Sampling plot failed: ", e$message)
  })
  
  # STEP 3: Threatened species summary
  species_threat <- community_data %>%
    filter(!is.na(species) & species != "") %>%
    distinct(Latin.Name, species, Threat.Status) %>%
    mutate(Threat.Status = ifelse(is.na(Threat.Status) | Threat.Status == "", "Not Assessed", Threat.Status))
  
  threat_summary <- species_threat %>%
    group_by(Threat.Status) %>%
    summarise(n_species = n(), .groups = "drop") %>%
    arrange(desc(n_species))
  
  write_csv(threat_summary, file.path(output_folder, "02_threat_summary.csv"))
  
  # STEP 4: NMDS ordination (attempt; bail out if it fails)
  results <- list()
  if (nrow(species_matrix) >= 2 && ncol(species_matrix) >= 2) {
    set.seed(123)
    nmds2 <- tryCatch(metaMDS(species_matrix, distance = "jaccard", k = 2, trymax = 100), error = function(e) NULL)
    if (!is.null(nmds2)) {
      stress_val <- round(nmds2$stress, 3)
      scores_df <- as.data.frame(scores(nmds2)$sites)
      scores_df$Representative_Site <- rownames(scores_df)
      scores_df <- left_join(scores_df, filtered_matrix %>% select(Representative_Site, Latitude, Longitude, Richness), by = "Representative_Site")
      results$nmds2 <- nmds2
      # save plot
      tryCatch({
        nmds_plot <- ggplot(scores_df, aes(x = NMDS1, y = NMDS2)) +
          geom_point(aes(size = Richness), color = "#0072B2", alpha = 0.7) +
          geom_text_repel(aes(label = Representative_Site), size = 3, max.overlaps = 20) +
          theme_minimal() +
          labs(title = paste("NMDS of", community_name), subtitle = paste("Stress:", stress_val))
        ggsave(file.path(output_folder, "03_NMDS_2D.png"), nmds_plot, width = 10, height = 8, dpi = 300, bg = "white")
      }, error = function(e) warning("Failed saving NMDS plot: ", e$message))
    }
  }
  
  # STEP 5: Clustering (attempt if possible)
  if (nrow(species_matrix) >= 2 && ncol(species_matrix) >= 2) {
    dist_matrix <- vegdist(species_matrix, method = "jaccard")
    sil_width <- numeric(9)
    for (k in 2:10) {
      pam_fit <- tryCatch(pam(dist_matrix, k = k), error = function(e) NULL)
      if (!is.null(pam_fit)) sil_width[k - 1] <- pam_fit$silinfo$avg.width else sil_width[k - 1] <- NA
    }
    optimal_k <- which.max(sil_width, na.rm = TRUE) + 1
    if (is.finite(optimal_k) && optimal_k >= 2) {
      pam_result <- pam(dist_matrix, k = optimal_k)
      clusters <- pam_result$clustering
      results$clusters <- clusters
      # save cluster summary
      cluster_summary <- data.frame(Cluster = names(table(clusters)), Number_of_Sites = as.numeric(table(clusters)))
      write_csv(cluster_summary, file.path(output_folder, "05_cluster_summary.csv"))
    }
  }
  
  # STEP 6: Species frequency by cluster (if clusters exist)
  # (This step mirrors parts of your original pipeline but is defensive)
  # Save species matrix and results
  saveRDS(list(
    species_matrix = species_matrix,
    aggregated_sites = filtered_matrix,
    results = results,
    threat_summary = threat_summary
  ), file.path(output_folder, "analysis_objects.rds"))
  
  cat("Finished analysis for", community_name, "->", output_folder, "\n")
  return(invisible(list(
    species_matrix = species_matrix,
    aggregated_sites = filtered_matrix,
    results = results
  )))
}