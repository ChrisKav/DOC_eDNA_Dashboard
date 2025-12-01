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

# Main analysis function - matched to your Script.R outputs and filenames
analyze_community <- function(df, 
                              community_data,
                              community_name,
                              output_folder,
                              rec2_rivers = "Data/REC2_Layers/River_Lines.shp",
                              nga_awa = "Data/Nga Awa shapefiles/DOC_NgāAwa_RiverSites_20250122_n14.shp",
                              grouping_distance = 200,
                              buffer_distance = 5000) {
  
  safe_dir_create(output_folder)
  cat("\n", rep("=", 70), "\n", sep = "")
  cat(paste0("ANALYZING ", toupper(community_name), " COMMUNITIES\n"))
  cat(rep("=", 70), "\n\n", sep = "")
  
  # Ensure columns exist to avoid errors
  if (!"Sample.Name" %in% names(community_data)) community_data$Sample.Name <- NA
  if (!"Latitude" %in% names(community_data)) community_data$Latitude <- NA
  if (!"Longitude" %in% names(community_data)) community_data$Longitude <- NA
  if (!"species" %in% names(community_data)) community_data$species <- NA
  if (!"Latin.Name" %in% names(community_data)) community_data$Latin.Name <- NA
  if (!"Threat.Status" %in% names(community_data)) community_data$Threat.Status <- NA
  if (!"DOC.Data" %in% names(df)) df$DOC.Data <- NA
  if (!"Date" %in% names(df)) df$Date <- NA
  
  # STEP 1: Data Preparation
  cat("Step 1: Preparing data...\n")
  site_coords <- community_data %>%
    select(Sample.Name, Latitude, Longitude) %>%
    distinct()
  
  community_matrix <- community_data %>%
    filter(!is.na(Sample.Name)) %>%
    group_by(Sample.Name, species) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(presence = 1) %>%
    pivot_wider(names_from = species, values_from = presence, values_fill = 0) %>%
    select(-any_of("NA")) %>%
    right_join(site_coords, by = "Sample.Name") %>%
    mutate(across(where(is.numeric) & !c(Latitude, Longitude), ~replace_na(.x, 0)))
  
  # Create spatial object if coords exist
  coords_available <- all(c("Longitude", "Latitude") %in% names(community_matrix)) && any(!is.na(community_matrix$Longitude) & !is.na(community_matrix$Latitude))
  if (coords_available) {
    community_sf <- st_as_sf(community_matrix %>% filter(!is.na(Longitude) & !is.na(Latitude)),
                             coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)
    # transform to NZTM if possible - match Script.R behaviour
    community_sf <- tryCatch(st_transform(community_sf, 2193), error = function(e) { warning("Failed to transform community_sf to 2193: ", e$message); community_sf })
  } else {
    # create an empty sf with no geometry but keep data for aggregation
    community_sf <- st_as_sf(community_matrix[0, , drop = FALSE], coords = c("Longitude","Latitude"), crs = 4326, remove = FALSE)
    community_sf$geometry <- st_sfc(crs = 2193)
  }
  
  # Group sites by distance (within community_sf coordinate system; fallback uses grouping_distance in metres)
  if (nrow(community_sf) > 0) {
    groups <- st_is_within_distance(community_sf, community_sf, dist = grouping_distance)
    community_sf$group_id <- sapply(groups, function(x) min(x))
  } else {
    community_sf$group_id <- integer(0)
  }
  
  # Aggregate by group_id
  merged_matrix <- community_matrix %>%
    mutate(group_id = community_sf$group_id) %>%
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
  
  filtered_matrix <- merged_matrix %>% filter(Richness > 0)
  
  # Prepare species matrix for analyses
  species_matrix <- filtered_matrix %>%
    select(-group_id, -Sample.Name, -Representative_Site, -Latitude, -Longitude, -Richness) %>%
    as.data.frame()
  rownames(species_matrix) <- filtered_matrix$Representative_Site
  
  cat(paste0("  - ", nrow(filtered_matrix), " sites with ", community_name, " present\n"))
  cat(paste0("  - ", ncol(species_matrix), " species detected\n\n"))
  
  # STEP 2: Sampling Effort Analysis
  cat("Step 2: Analyzing sampling effort...\n")
  sampling_summary <- df %>%
    mutate(Year = suppressWarnings(as.numeric(format(as.Date(Date), "%Y")))) %>%
    distinct(Sample.Name, Year, DOC.Data) %>%
    mutate(DOC.Data = ifelse(DOC.Data == "Yes", "DOC", "Non-DOC")) %>%
    group_by(Year, DOC.Data) %>%
    summarise(n_samples = n(), .groups = "drop") %>%
    mutate(DOC.Data = factor(DOC.Data, levels = c("DOC", "Non-DOC")))
  
  sampling_totals <- sampling_summary %>%
    group_by(Year) %>%
    summarise(total = sum(n_samples), .groups = "drop")
  
  sampling_table <- sampling_summary %>%
    pivot_wider(names_from = DOC.Data, values_from = n_samples, values_fill = 0) %>%
    left_join(sampling_totals, by = "Year") %>%
    select(Year, DOC, `Non-DOC`, Total = total)
  
  # Save docx table like Script.R
  ft_sampling <- sampling_table %>%
    flextable() %>%
    theme_vanilla() %>%
    set_header_labels(Year = "Year", DOC = "DOC Samples", `Non-DOC` = "Non-DOC Samples", Total = "Total") %>%
    bold(part = "header") %>%
    autofit()
  tryCatch(save_as_docx(ft_sampling, path = file.path(output_folder, "01_Sampling_Summary_by_Year.docx")), error = function(e) warning("Failed saving sampling docx: ", e$message))
  
  # Save sampling figure PNG with same filename
  sampling_plot <- ggplot(sampling_summary, aes(x = factor(Year), y = n_samples, fill = fct_rev(DOC.Data))) +
    geom_col(color = "black", width = 0.7, linewidth = 0.5) +
    geom_text(data = sampling_totals, aes(x = factor(Year), y = total, label = total, fill = NULL),
              vjust = -0.5, size = 4, fontface = "bold") +
    scale_fill_manual(values = c("DOC" = "#2E7D32", "Non-DOC" = "#1976D2"), name = "Data Source") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5)) +
    labs(title = paste(community_name, "Sampling Effort Over Time"),
         subtitle = paste("Total unique sites sampled:", sum(sampling_summary$n_samples, na.rm = TRUE)),
         x = "Year", y = "Number of Samples")
  tryCatch(ggsave(file.path(output_folder, "01_Sampling_Effort_by_Year.png"), sampling_plot, width = 8, height = 6, dpi = 300, bg = "white"), error = function(e) warning("Failed saving sampling_plot: ", e$message))
  
  # STEP 3: Threatened Species Analysis
  cat("Step 3: Analyzing threatened species...\n")
  species_threat <- community_data %>%
    filter(!is.na(species) & species != "") %>%
    distinct(Latin.Name, species, Threat.Status) %>%
    mutate(Threat.Status = case_when(
      is.na(Threat.Status) | Threat.Status == "" ~ "Not Assessed",
      TRUE ~ Threat.Status
    ))
  
  threat_summary <- species_threat %>%
    group_by(Threat.Status) %>%
    summarise(n_species = n(), .groups = "drop") %>%
    arrange(desc(n_species))
  
  threat_details <- species_threat %>%
    select(`Threat Status` = Threat.Status, `Species Name` = Latin.Name) %>%
    arrange(`Threat Status`, `Species Name`)
  
  ft_threat <- threat_details %>%
    as_grouped_data(groups = "Threat Status") %>%
    as_flextable() %>%
    theme_vanilla() %>%
    bold(j = 1, i = ~ !is.na(`Threat Status`)) %>%
    bg(j = 1, i = ~ !is.na(`Threat Status`), bg = "#EFEFEF") %>%
    autofit()
  tryCatch(save_as_docx(ft_threat, path = file.path(output_folder, "02_Threatened_Species_List.docx")), error = function(e) warning("Failed saving threat docx: ", e$message))
  
  # Threat figure
  threat_plot <- ggplot(threat_summary, aes(x = reorder(Threat.Status, n_species), y = n_species)) +
    geom_col(aes(fill = Threat.Status), color = "black", width = 0.7) +
    geom_text(aes(label = n_species), hjust = -0.3, size = 4, fontface = "bold") +
    scale_fill_manual(values = c("Not Threatened" = "#4CAF50", "Declining" = "#FFC107",
                                 "At Risk" = "#FF9800", "Threatened" = "#F44336",
                                 "Nationally Critical" = "#B71C1C", "Nationally Endangered" = "#D32F2F",
                                 "Nationally Vulnerable" = "#E64A19", "Not Assessed" = "#9E9E9E")) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.y = element_blank(), legend.position = "none",
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
    labs(title = paste(community_name, "Species by Conservation Status"),
         subtitle = paste("Total species identified:", sum(threat_summary$n_species, na.rm = TRUE)),
         x = "Threat Status", y = "Number of Species")
  tryCatch(ggsave(file.path(output_folder, "02_Threatened_Species_Summary.png"), threat_plot, width = 8, height = 6, dpi = 300, bg = "white"), error = function(e) warning("Failed saving threat plot: ", e$message))
  
  # STEP 4: NMDS Ordination
  cat("Step 4: Running NMDS ordination...\n")
  nmds2 <- NULL; nmds3 <- NULL; scores_df <- NULL; scores_3d <- NULL; species_scores <- NULL; stress_val <- NA
  if (nrow(species_matrix) >= 2 && ncol(species_matrix) >= 1) {
    set.seed(123)
    nmds2 <- tryCatch(metaMDS(species_matrix, distance = "jaccard", k = 2, trymax = 100), error = function(e) { warning("metaMDS 2D failed: ", e$message); NULL })
    nmds3 <- tryCatch(metaMDS(species_matrix, distance = "jaccard", k = 3, trymax = 100), error = function(e) { warning("metaMDS 3D failed: ", e$message); NULL })
    if (!is.null(nmds2)) {
      scores_df <- as.data.frame(scores(nmds2)$sites)
      scores_df$Representative_Site <- rownames(scores_df)
      scores_df <- left_join(scores_df, filtered_matrix %>% select(Representative_Site, Latitude, Longitude, Richness), by = "Representative_Site")
      stress_val <- round(nmds2$stress, 3)
      cat(paste0("  - 2D NMDS stress: ", stress_val, "\n"))
      # species envfit
      species_fit <- tryCatch(envfit(nmds2, species_matrix, permutations = 999), error = function(e) NULL)
      if (!is.null(species_fit)) {
        species_scores <- as.data.frame(scores(species_fit, "vectors"))
        species_scores$species <- rownames(species_scores)
      }
      # plot NMDS
      nmds_plot <- ggplot(scores_df, aes(x = NMDS1, y = NMDS2)) +
        geom_point(aes(size = Richness), color = "#0072B2", alpha = 0.7) +
        geom_text_repel(aes(label = Representative_Site), size = 3, max.overlaps = 20)
      if (!is.null(species_scores)) {
        nmds_plot <- nmds_plot +
          geom_segment(data = species_scores, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                       arrow = arrow(length = unit(0.25, "cm")), color = "darkred") +
          geom_text_repel(data = species_scores, aes(x = NMDS1, y = NMDS2, label = species),
                          color = "darkred", size = 3)
      }
      nmds_plot <- nmds_plot +
        theme_minimal(base_size = 14) +
        labs(title = paste("NMDS Ordination of", community_name, "Communities"),
             subtitle = paste("Jaccard distance | Stress =", stress_val),
             x = "NMDS Axis 1", y = "NMDS Axis 2", size = "Richness") +
        theme(panel.grid = element_blank())
      tryCatch(ggsave(file.path(output_folder, "03_NMDS_Ordination_Plot.png"), nmds_plot, width = 10, height = 8, dpi = 300, bg = "white"), error = function(e) warning("Failed saving NMDS plot: ", e$message))
    }
    if (!is.null(nmds3)) {
      scores_3d <- as.data.frame(scores(nmds3, display = "sites"))
      scores_3d$Site <- rownames(scores_3d)
    }
  } else {
    warning("Not enough site/species data to perform NMDS.")
  }
  
  # STEP 5: Clustering Analysis
  cat("Step 5: Performing clustering analysis...\n")
  clusters <- NULL
  if (nrow(species_matrix) >= 2 && ncol(species_matrix) >= 1) {
    dist_matrix <- vegdist(species_matrix, method = "jaccard")
    # Silhouette analysis from k=2 to k=10
    sil_width <- numeric(9)
    for(k in 2:10) {
      pam_fit <- tryCatch(pam(dist_matrix, k = k), error = function(e) NULL)
      sil_width[k-1] <- if (!is.null(pam_fit)) pam_fit$silinfo$avg.width else NA
    }
    # Silhouette plot
    tryCatch({
      png(file.path(output_folder, "04_Silhouette_Analysis.png"), width = 2400, height = 1800, res = 300)
      plot(2:10, sil_width, type = "b", pch = 19, 
           xlab = "Number of clusters", ylab = "Average silhouette width",
           main = paste("Optimal Clusters via Silhouette Method -", community_name))
      abline(v = which.max(sil_width) + 1, lty = 2, col = "red")
      dev.off()
    }, error = function(e) warning("Failed saving silhouette plot: ", e$message))
    optimal_k <- which.max(sil_width) + 1
    if (!is.na(optimal_k) && optimal_k >= 2) {
      cat(paste0("  - Optimal clusters: ", optimal_k, "\n"))
      set.seed(123)
      pam_result <- pam(dist_matrix, k = optimal_k)
      clusters <- pam_result$clustering
      # Dendrogram
      tryCatch({
        png(file.path(output_folder, "05_Dendrogram.png"), width = 3000, height = 2000, res = 300)
        hclust_result <- hclust(dist_matrix, method = "average")
        plot(hclust_result, main = paste("Dendrogram of", community_name, "Communities (Average Linkage)"),
             xlab = "Sites", sub = "", hang = -1)
        rect.hclust(hclust_result, k = optimal_k, border = "red")
        dev.off()
      }, error = function(e) warning("Failed saving dendrogram: ", e$message))
      # Cluster summary docx (06_)
      cluster_summary <- data.frame(
        Cluster = names(table(clusters)),
        Number_of_Sites = as.numeric(table(clusters))
      )
      ft_cluster_summary <- cluster_summary %>%
        flextable() %>%
        theme_vanilla() %>%
        set_header_labels(Cluster = "Cluster", Number_of_Sites = "Number of Sites") %>%
        bold(part = "header") %>%
        autofit()
      tryCatch(save_as_docx(ft_cluster_summary, path = file.path(output_folder, "06_Cluster_Summary.docx")), error = function(e) warning("Failed saving cluster summary docx: ", e$message))
    } else {
      warning("Could not determine optimal_k for clustering.")
    }
  } else {
    warning("Not enough data for clustering analysis.")
  }
  
  # STEP 6: Species-Cluster Associations
  cat("Step 6: Analyzing species-cluster associations...\n")
  if (!is.null(clusters) && !is.null(scores_3d)) {
    filtered_matrix_clustered <- filtered_matrix %>%
      left_join(scores_3d %>% select(Site, Cluster = as.character(ifelse(is.null(clusters), NA, clusters[match(as.character(scores_3d$Site), names(clusters))])) ), by = c("Representative_Site" = "Site"))
    # The above join is conservative; easier approach below:
    filtered_matrix_clustered <- filtered_matrix
    filtered_matrix_clustered$Cluster <- as.character(clusters[match(filtered_matrix_clustered$Representative_Site, names(clusters))])
    
    # species frequency per cluster
    species_freq <- filtered_matrix_clustered %>%
      select(-group_id, -Sample.Name, -Representative_Site, -Latitude, -Longitude, -Richness) %>%
      pivot_longer(cols = -Cluster, names_to = "Species", values_to = "Presence") %>%
      filter(Presence > 0) %>%
      count(Cluster, Species) %>%
      group_by(Cluster) %>%
      mutate(total_sites = n_distinct(filtered_matrix_clustered$Representative_Site[filtered_matrix_clustered$Cluster == Cluster]),
             frequency = n / total_sites * 100) %>%
      arrange(Cluster, desc(frequency))
    
    clean_data <- species_freq %>%
      filter(!is.na(Species)) %>%
      select(Cluster, Species, n, frequency) %>%
      arrange(Cluster, desc(frequency)) %>%
      mutate(frequency = round(frequency, 2))
    
    ft_species_freq <- clean_data %>%
      as_grouped_data(groups = "Cluster") %>%
      as_flextable() %>%
      theme_vanilla() %>%
      set_header_labels(Species = "Species Name", n = "Count", frequency = "Frequency (%)") %>%
      autofit() %>%
      bg(j = 1, i = ~ !is.na(Cluster), bg = "#EFEFEF") %>% 
      bold(j = 1, i = ~ !is.na(Cluster), bold = TRUE)
    tryCatch(save_as_docx(ft_species_freq, path = file.path(output_folder, "07_Species_Frequency_by_Cluster.docx")), error = function(e) warning("Failed saving species frequency docx: ", e$message))
    
    # Indicator species analysis
    cluster_vector <- as.numeric(as.factor(filtered_matrix_clustered$Cluster))
    species_only <- filtered_matrix_clustered %>%
      select(-group_id, -Sample.Name, -Representative_Site, -Latitude, -Longitude, -Richness, -Cluster)
    if (ncol(species_only) > 0) {
      indval_result <- tryCatch(multipatt(species_only, cluster_vector, control = how(nperm = 999)), error = function(e) NULL)
      if (!is.null(indval_result)) {
        sign_table <- indval_result$sign
        if (nrow(sign_table) > 0) {
          significant_indicators <- sign_table %>%
            as.data.frame() %>%
            rownames_to_column(var = "Species") %>%
            mutate(Cluster = as.character(index),
                   Indicator_Value = round(stat, 3),
                   p_value = round(p.value, 4)) %>%
            select(Species, Cluster, Indicator_Value, p_value) %>%
            arrange(Cluster, desc(Indicator_Value))
          if (nrow(significant_indicators) > 0) {
            ft_indicator <- significant_indicators %>%
              flextable() %>%
              theme_vanilla() %>%
              set_header_labels(Species = "Species Name", Cluster = "Associated Cluster",
                                Indicator_Value = "IndVal Statistic", p_value = "p-value") %>%
              autofit() %>%
              bold(part = "header") %>%
              bg(i = ~ p_value < 0.01, bg = "#FFEB9C", part = "body") %>%
              add_footer_lines("IndVal ranges from 0 to 1. p < 0.05 indicates significant association.") %>%
              fontsize(size = 9, part = "footer")
            tryCatch(save_as_docx(ft_indicator, path = file.path(output_folder, "08_Indicator_Species_Analysis.docx")), error = function(e) warning("Failed saving indicator docx: ", e$message))
          }
        }
      }
    }
    
    # Heatmap
    freq_matrix <- species_freq %>%
      select(Cluster, Species, frequency) %>%
      pivot_wider(names_from = Cluster, values_from = frequency, values_fill = 0) %>%
      column_to_rownames("Species") %>%
      as.matrix()
    tryCatch({
      png(file.path(output_folder, "09_Species_Frequency_Heatmap.png"), width = 3000, height = 3000, res = 300)
      pheatmap(freq_matrix,
               color = colorRampPalette(c("white", "lightblue", "darkblue"))(50),
               cluster_rows = TRUE, cluster_cols = FALSE,
               main = paste("Species Frequency (%) by Cluster -", community_name),
               fontsize = 10, angle_col = 0, cellwidth = 30, cellheight = 12,
               display_numbers = round(freq_matrix, 0), number_color = "black", fontsize_number = 8)
      dev.off()
    }, error = function(e) warning("Failed saving heatmap: ", e$message))
    
    characteristic_species <- species_freq %>%
      filter(frequency >= 50) %>%
      group_by(Cluster) %>%
      summarise(n_characteristic = n(),
                characteristic_species = paste(Species, collapse = ", "),
                .groups = "drop")
    ft_characteristic <- characteristic_species %>%
      flextable() %>%
      theme_vanilla() %>%
      set_header_labels(Cluster = "Cluster", n_characteristic = "Number of Species",
                        characteristic_species = "Species (≥50% frequency)") %>%
      bold(part = "header") %>%
      autofit()
    tryCatch(save_as_docx(ft_characteristic, path = file.path(output_folder, "10_Characteristic_Species_Summary.docx")), error = function(e) warning("Failed saving characteristic docx: ", e$message))
  } else {
    warning("Cluster-based species analyses skipped because clustering did not produce results.")
  }
  
  # STEP 7: Spatial Maps (if spatial data provided)
  if (!is.null(rec2_rivers) && !is.null(nga_awa) && !is.null(scores_df)) {
    cat("Step 7: Creating spatial maps...\n")
    # ensure spatial layers are in 2193 like Script.R
    rec2_rivers_try <- tryCatch(st_transform(rec2_rivers, 2193), error = function(e) rec2_rivers)
    nga_awa_try <- tryCatch(st_transform(nga_awa, 2193), error = function(e) nga_awa)
    
    # Build community_sf in 2193 (if not already)
    community_sf_try <- tryCatch(st_transform(community_sf, 2193), error = function(e) community_sf)
    sampling_buffer <- st_buffer(community_sf_try, dist = buffer_distance)
    sampling_extent_nztm <- st_bbox(st_union(sampling_buffer))
    sampling_extent_sf <- st_as_sfc(sampling_extent_nztm) %>% st_transform(4326)
    sampling_extent <- st_bbox(sampling_extent_sf)
    
    # Crop spatial layers - use safe tryCatch around spatial ops
    nga_awa_crop <- tryCatch(st_crop(nga_awa_try, sampling_extent_nztm), error = function(e) nga_awa_try)
    rec2_trimmed <- tryCatch(st_intersection(rec2_rivers_try, nga_awa_crop), error = function(e) rec2_rivers_try)
    rec2_sampling <- tryCatch(st_crop(rec2_trimmed, sampling_extent_nztm), error = function(e) rec2_trimmed)
    
    # Prepare map data
    map_data <- merged_matrix %>%
      left_join(if (!is.null(scores_3d)) scores_3d %>% select(Site, Cluster) else tibble(Site = character(), Cluster = character()),
                by = c("Representative_Site" = "Site"))
    num_clusters <- length(unique(na.omit(map_data$Cluster)))
    palette_colors <- scales::hue_pal()(max(1, num_clusters))
    
    # Map: NMDS Axis 1
    map_nmds1 <- ggplot() +
      geom_sf(data = rec2_sampling, aes(size = StreamOrde), color = "steelblue", show.legend = FALSE) +
      scale_size(range = c(0.1, 1.5)) +
      geom_point(data = scores_df, aes(x = Longitude, y = Latitude, fill = NMDS1), 
                 size = 6, alpha = 0.8, shape = 21, color = "black", stroke = 0.5) +
      scale_fill_viridis_c(option = "plasma") +
      annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.9, line_width = 1) +
      annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering,
                             height = unit(1.2, "cm"), width = unit(1.2, "cm")) +
      coord_sf(crs = 4326, xlim = c(sampling_extent["xmin"], sampling_extent["xmax"]),
               ylim = c(sampling_extent["ymin"], sampling_extent["ymax"])) +
      theme_bw(base_size = 12) +
      theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
            panel.grid.minor = element_blank(), legend.position = "right",
            plot.title = element_text(size = 13, face = "bold")) +
      labs(title = "NMDS Axis 1", x = "Longitude", y = "Latitude", fill = "NMDS1")
    # Map: NMDS Axis 2
    map_nmds2 <- map_nmds1 %+% map_nmds1 + labs(title = "NMDS Axis 2")
    # Actually build map_nmds2 properly to use NMDS2
    map_nmds2 <- ggplot() +
      geom_sf(data = rec2_sampling, aes(size = StreamOrde), color = "steelblue", show.legend = FALSE) +
      scale_size(range = c(0.1, 1.5)) +
      geom_point(data = scores_df, aes(x = Longitude, y = Latitude, fill = NMDS2), 
                 size = 6, alpha = 0.8, shape = 21, color = "black", stroke = 0.5) +
      scale_fill_viridis_c(option = "plasma") +
      annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.9, line_width = 1) +
      annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering,
                             height = unit(1.2, "cm"), width = unit(1.2, "cm")) +
      coord_sf(crs = 4326, xlim = c(sampling_extent["xmin"], sampling_extent["xmax"]),
               ylim = c(sampling_extent["ymin"], sampling_extent["ymax"])) +
      theme_bw(base_size = 12) +
      theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
            panel.grid.minor = element_blank(), legend.position = "right",
            plot.title = element_text(size = 13, face = "bold")) +
      labs(title = "NMDS Axis 2", x = "Longitude", y = "Latitude", fill = "NMDS2")
    
    # Combine and save NMDS spatial gradient
    nmds_combined <- map_nmds1 + map_nmds2 +
      plot_annotation(title = paste("NMDS Ordination of", community_name, "Communities"),
                      subtitle = paste("Jaccard distance | Stress =", ifelse(is.na(stress_val), "", stress_val)),
                      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                                    plot.subtitle = element_text(size = 12, hjust = 0.5)))
    tryCatch(ggsave(file.path(output_folder, "11_NMDS_Spatial_Gradients.png"), nmds_combined, width = 14, height = 6, dpi = 300, bg = "white"), error = function(e) warning("Failed saving NMDS spatial gradients: ", e$message))
    
    # Sampling locations map
    sampling_map <- ggplot() +
      geom_sf(data = rec2_sampling, aes(size = StreamOrde), color = "steelblue", show.legend = FALSE) +
      scale_size(range = c(0.3, 2)) +
      geom_point(data = merged_matrix, aes(x = Longitude, y = Latitude),
                 fill = "darkred", color = "black", shape = 21, size = 5, alpha = 0.8, stroke = 0.8) +
      geom_point(data = merged_matrix %>% filter(Richness == 0), aes(x = Longitude, y = Latitude),
                 color = "black", shape = 1, size = 5, stroke = 1.2) +
      annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.9, line_width = 1) +
      annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering,
                             height = unit(1.2, "cm"), width = unit(1.2, "cm")) +
      coord_sf(crs = 4326, xlim = c(sampling_extent["xmin"], sampling_extent["xmax"]),
               ylim = c(sampling_extent["ymin"], sampling_extent["ymax"])) +
      theme_bw(base_size = 12) +
      theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 11)) +
      labs(title = "Sampling Locations",
           subtitle = paste("Filled circles =", tolower(community_name), "present, hollow circles = none detected"),
           x = "Longitude", y = "Latitude")
    tryCatch(ggsave(file.path(output_folder, "12_Sampling_Locations_Map.png"), sampling_map, width = 8, height = 7, dpi = 300, bg = "white"), error = function(e) warning("Failed saving sampling_map: ", e$message))
    
    # Cluster map
    cluster_map <- ggplot() +
      geom_sf(data = rec2_sampling, aes(size = StreamOrde), color = "steelblue", show.legend = FALSE) +
      scale_size(range = c(0.3, 2)) +
      geom_point(data = map_data %>% filter(Richness > 0),
                 aes(x = Longitude, y = Latitude, fill = Cluster),
                 size = 6, alpha = 0.9, shape = 21, color = "black", stroke = 0.8) +
      geom_point(data = map_data %>% filter(Richness == 0),
                 aes(x = Longitude, y = Latitude, color = Cluster),
                 shape = 1, size = 6, stroke = 1.2) +
      scale_fill_manual(values = palette_colors, name = "Cluster") +
      scale_color_manual(values = palette_colors, name = "Cluster") +
      annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.9, line_width = 1) +
      annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering,
                             height = unit(1.2, "cm"), width = unit(1.2, "cm")) +
      coord_sf(crs = 4326, xlim = c(sampling_extent["xmin"], sampling_extent["xmax"]),
               ylim = c(sampling_extent["ymin"], sampling_extent["ymax"])) +
      theme_bw(base_size = 12) +
      theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
            legend.position = "right", plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 11)) +
      labs(title = paste(community_name, "Community Clusters"),
           subtitle = paste("PAM clustering (k =", ifelse(exists("optimal_k"), optimal_k, NA), ") based on Jaccard distance"),
           x = "Longitude", y = "Latitude")
    tryCatch(ggsave(file.path(output_folder, "13_Cluster_Map.png"), cluster_map, width = 9, height = 7, dpi = 300, bg = "white"), error = function(e) warning("Failed saving cluster_map: ", e$message))
    
    cat("  - 3 spatial maps created\n")
  } else {
    cat("Spatial maps skipped (rec2_rivers, nga_awa or NMDS/site scores missing)\n")
  }
  
  # Completion message and return values (match Script.R return)
  cat("\n", rep("=", 70), "\n", sep = "")
  cat(paste0(toupper(community_name), " ANALYSIS COMPLETE!\n"))
  cat(paste0("All outputs saved to '", output_folder, "' folder\n"))
  cat(rep("=", 70), "\n\n", sep = "")
  
  return(list(
    species_matrix = species_matrix,
    nmds2 = nmds2,
    nmds3 = nmds3,
    clusters = clusters,
    filtered_matrix = filtered_matrix,
    merged_matrix = merged_matrix,
    community_sf = community_sf
  ))
}