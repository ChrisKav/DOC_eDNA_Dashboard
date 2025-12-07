# ============================================================
# MODULAR COMMUNITY ANALYSIS SCRIPT
# Analyze any biological community (Fish, Macroinvertebrates, etc.)
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
# 1. MAIN ANALYSIS FUNCTION
# -----------------------------
analyze_community <- function(df, 
                              community_data,
                              community_name,
                              output_folder,
                              rec2_rivers = NULL,
                              nga_awa = NULL,
                              grouping_distance = 200,
                              buffer_distance = 5000) {
  
  # Create output folder
  if(!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
    cat(paste0("Created '", output_folder, "' folder\n"))
  }
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat(paste0("ANALYZING ", toupper(community_name), " COMMUNITIES\n"))
  cat(rep("=", 70), "\n\n", sep = "")
  
  # -----------------------------
  # STEP 1: Data Preparation
  # -----------------------------
  cat("Step 1: Preparing data...\n")
  
  site_coords <- community_data %>%
    select(Sample.Name, Latitude, Longitude) %>%
    distinct()
  
  # Build presence-absence matrix
  community_matrix <- community_data %>%
    group_by(Sample.Name, species) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(presence = 1) %>%
    pivot_wider(names_from = species, values_from = presence, values_fill = 0) %>%
    select(-any_of("NA")) %>%
    right_join(site_coords, by = "Sample.Name") %>%
    mutate(across(where(is.numeric) & !c(Latitude, Longitude), ~replace_na(.x, 0)))
  
  # Create spatial object
  community_sf <- st_as_sf(community_matrix, 
                           coords = c("Longitude", "Latitude"), 
                           crs = 4326) %>%
    st_transform(2193)
  
  # Group sites
  groups <- st_is_within_distance(community_sf, community_sf, dist = grouping_distance)
  community_sf$group_id <- sapply(groups, function(x) min(x))
  
  # Aggregate by group
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
  
  # Prepare species matrix
  species_matrix <- filtered_matrix %>%
    select(-group_id, -Sample.Name, -Representative_Site, -Latitude, -Longitude, -Richness) %>%
    as.data.frame()
  rownames(species_matrix) <- filtered_matrix$Representative_Site
  
  cat(paste0("  - ", nrow(filtered_matrix), " sites with ", community_name, " present\n"))
  cat(paste0("  - ", ncol(species_matrix), " species detected\n\n"))
  
  # -----------------------------
  # STEP 2: Sampling Effort Analysis (FIXED)
  # -----------------------------
  cat("Step 2: Analyzing sampling effort...\n")
  
  # Convert Date to proper date format and extract year (try several common formats)
  sampling_summary <- df %>%
    mutate(
      Date = as.character(Date),
      Date = as.Date(Date, format = "%Y-%m-%d"),
      Date = if_else(is.na(Date), as.Date(Date, format = "%d/%m/%Y"), Date),
      Date = if_else(is.na(Date), as.Date(Date, format = "%m/%d/%Y"), Date),
      Year = as.numeric(format(Date, "%Y"))
    ) %>%
    filter(!is.na(Year)) %>%  # Remove any rows where year extraction failed
    distinct(Sample.Name, Year, DOC.Data) %>%
    mutate(DOC.Data = ifelse(DOC.Data == "Yes", "DOC", "Non-DOC")) %>%
    group_by(Year, DOC.Data) %>%
    summarise(n_samples = n(), .groups = "drop") %>%
    mutate(DOC.Data = factor(DOC.Data, levels = c("DOC", "Non-DOC")))
  
  # Totals by year
  sampling_totals <- sampling_summary %>%
    group_by(Year) %>%
    summarise(total = sum(n_samples), .groups = "drop")
  
  # Pivot to wide table and ensure both DOC and Non-DOC columns exist (fill missing with 0)
  sampling_table <- sampling_summary %>%
    pivot_wider(names_from = DOC.Data, values_from = n_samples, values_fill = 0) %>%
    { if (!"DOC" %in% names(.)) mutate(., DOC = 0) else . } %>%
    { if (!"Non-DOC" %in% names(.)) mutate(., `Non-DOC` = 0) else . } %>%
    left_join(sampling_totals, by = "Year") %>%
    select(Year, DOC, `Non-DOC`, Total = total)
  
  # Prepare flextable and plot as before (no change)
  ft_sampling <- sampling_table %>%
    flextable() %>%
    theme_vanilla() %>%
    set_header_labels(Year = "Year", DOC = "DOC Samples", `Non-DOC` = "Non-DOC Samples", Total = "Total") %>%
    bold(part = "header") %>%
    autofit()
  
  save_as_docx(ft_sampling, path = file.path(output_folder, "01_Sampling_Summary_by_Year.docx"))
  
  # Figure
  sampling_plot <- ggplot(
    sampling_summary %>% mutate(DOC.Data = factor(DOC.Data, levels = c("DOC", "Non-DOC"))),
    aes(x = factor(Year), y = n_samples, fill = fct_rev(DOC.Data))
  ) +
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
         subtitle = paste("Total unique sites sampled:", sum(sampling_summary$n_samples)),
         x = "Year", y = "Number of Samples")
  
  ggsave(file.path(output_folder, "01_Sampling_Effort_by_Year.png"), sampling_plot, 
         width = 8, height = 6, dpi = 300, bg = "white")
  
  # -----------------------------
  # STEP 3: Threatened Species Analysis
  # -----------------------------
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
  
  # Table
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
  
  save_as_docx(ft_threat, path = file.path(output_folder, "02_Threatened_Species_List.docx"))
  
  # Figure
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
         subtitle = paste("Total species identified:", sum(threat_summary$n_species)),
         x = "Threat Status", y = "Number of Species")
  
  ggsave(file.path(output_folder, "02_Threatened_Species_Summary.png"), threat_plot, 
         width = 8, height = 6, dpi = 300, bg = "white")
  
  # -----------------------------
  # STEP 4: NMDS Ordination
  # -----------------------------
  cat("Step 4: Running NMDS ordination...\n")
  
  set.seed(123)
  nmds2 <- metaMDS(species_matrix, distance = "jaccard", k = 2, trymax = 100)
  nmds3 <- metaMDS(species_matrix, distance = "jaccard", k = 3, trymax = 100)
  
  scores_df <- as.data.frame(scores(nmds2)$sites)
  scores_df$Representative_Site <- rownames(scores_df)
  scores_df <- left_join(scores_df, 
                         filtered_matrix %>% select(Representative_Site, Latitude, Longitude, Richness),
                         by = "Representative_Site")
  
  stress_val <- round(nmds2$stress, 3)
  cat(paste0("  - 2D NMDS stress: ", stress_val, "\n"))
  
  scores_3d <- as.data.frame(scores(nmds3, display = "sites"))
  scores_3d$Site <- rownames(scores_3d)
  
  # Fit species vectors
  species_fit <- envfit(nmds2, species_matrix, permutations = 999)
  species_scores <- as.data.frame(scores(species_fit, "vectors"))
  species_scores$species <- rownames(species_scores)
  
  # Figure
  nmds_plot <- ggplot(scores_df, aes(x = NMDS1, y = NMDS2)) +
    geom_point(aes(size = Richness), color = "#0072B2", alpha = 0.7) +
    geom_text_repel(aes(label = Representative_Site), size = 3, max.overlaps = 20) +
    geom_segment(data = species_scores, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                 arrow = arrow(length = unit(0.25, "cm")), color = "darkred") +
    geom_text_repel(data = species_scores, aes(x = NMDS1, y = NMDS2, label = species),
                    color = "darkred", size = 3) +
    theme_minimal(base_size = 14) +
    labs(title = paste("NMDS Ordination of", community_name, "Communities"),
         subtitle = paste("Jaccard distance | Stress =", stress_val),
         x = "NMDS Axis 1", y = "NMDS Axis 2", size = "Richness") +
    theme(panel.grid = element_blank())
  
  ggsave(file.path(output_folder, "03_NMDS_Ordination_Plot.png"), nmds_plot, 
         width = 10, height = 8, dpi = 300, bg = "white")
  
  # -----------------------------
  # STEP 5: Clustering Analysis
  # -----------------------------
  cat("Step 5: Performing clustering analysis...\n")
  
  dist_matrix <- vegdist(species_matrix, method = "jaccard")
  
  # Silhouette analysis
  sil_width <- numeric(9)
  for(k in 2:10) {
    pam_fit <- pam(dist_matrix, k = k)
    sil_width[k-1] <- pam_fit$silinfo$avg.width
  }
  
  # Figure: Silhouette
  png(file.path(output_folder, "04_Silhouette_Analysis.png"), width = 2400, height = 1800, res = 300)
  plot(2:10, sil_width, type = "b", pch = 19, 
       xlab = "Number of clusters", ylab = "Average silhouette width",
       main = paste("Optimal Clusters via Silhouette Method -", community_name))
  abline(v = which.max(sil_width) + 1, lty = 2, col = "red")
  dev.off()
  
  optimal_k <- which.max(sil_width) + 1
  cat(paste0("  - Optimal clusters: ", optimal_k, "\n"))
  
  # PAM clustering
  set.seed(123)
  pam_result <- pam(dist_matrix, k = optimal_k)
  clusters <- pam_result$clustering
  scores_3d$Cluster <- factor(clusters[match(scores_3d$Site, names(clusters))])
  
  # Figure: Dendrogram
  png(file.path(output_folder, "05_Dendrogram.png"), width = 3000, height = 2000, res = 300)
  hclust_result <- hclust(dist_matrix, method = "average")
  plot(hclust_result, main = paste("Dendrogram of", community_name, "Communities (Average Linkage)"),
       xlab = "Sites", sub = "", hang = -1)
  rect.hclust(hclust_result, k = optimal_k, border = "red")
  dev.off()
  
  # Table: Cluster Summary
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
  
  save_as_docx(ft_cluster_summary, path = file.path(output_folder, "06_Cluster_Summary.docx"))
  
  # -----------------------------
  # STEP 6: Species-Cluster Associations
  # -----------------------------
  cat("Step 6: Analyzing species-cluster associations...\n")
  
  filtered_matrix_clustered <- filtered_matrix %>%
    left_join(scores_3d %>% select(Site, Cluster), by = c("Representative_Site" = "Site"))
  
  # Calculate species frequency
  species_freq <- filtered_matrix_clustered %>%
    select(-group_id, -Sample.Name, -Representative_Site, -Latitude, -Longitude, -Richness) %>%
    pivot_longer(cols = -Cluster, names_to = "Species", values_to = "Presence") %>%
    filter(Presence > 0) %>%
    count(Cluster, Species) %>%
    group_by(Cluster) %>%
    mutate(total_sites = sum(table(filtered_matrix_clustered$Cluster)[as.character(Cluster)]),
           frequency = n / total_sites * 100) %>%
    arrange(Cluster, desc(frequency))
  
  # Table: Species Frequency
  clean_data <- species_freq %>%
    filter(Species != "count") %>%
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
  
  save_as_docx(ft_species_freq, path = file.path(output_folder, "07_Species_Frequency_by_Cluster.docx"))
  
  # Indicator species analysis
  cluster_vector <- as.numeric(filtered_matrix_clustered$Cluster)
  species_only <- filtered_matrix_clustered %>%
    select(-group_id, -Sample.Name, -Representative_Site, -Latitude, -Longitude, -Richness, -Cluster)
  
  indval_result <- multipatt(species_only, cluster_vector, control = how(nperm = 999))
  
  # Table: Indicator Species
  sign_table <- indval_result$sign
  
  if(nrow(sign_table) > 0) {
    significant_indicators <- sign_table %>%
      filter(p.value < 0.05) %>%
      mutate(Species = rownames(.), Cluster = as.character(index),
             Indicator_Value = round(stat, 3), p_value = round(p.value, 4)) %>%
      select(Species, Cluster, Indicator_Value, p_value) %>%
      arrange(Cluster, desc(Indicator_Value))
    
    if(nrow(significant_indicators) > 0) {
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
      
      save_as_docx(ft_indicator, path = file.path(output_folder, "08_Indicator_Species_Analysis.docx"))
    }
  }
  
  # Figure: Heatmap
  freq_matrix <- species_freq %>%
    select(Cluster, Species, frequency) %>%
    pivot_wider(names_from = Cluster, values_from = frequency, values_fill = 0) %>%
    column_to_rownames("Species") %>%
    as.matrix()
  
  png(file.path(output_folder, "09_Species_Frequency_Heatmap.png"), width = 3000, height = 3000, res = 300)
  pheatmap(freq_matrix,
           color = colorRampPalette(c("white", "lightblue", "darkblue"))(50),
           cluster_rows = TRUE, cluster_cols = FALSE,
           main = paste("Species Frequency (%) by Cluster -", community_name),
           fontsize = 10, angle_col = 0, cellwidth = 30, cellheight = 12,
           display_numbers = round(freq_matrix, 0), number_color = "black", fontsize_number = 8)
  dev.off()
  
  # Table: Characteristic Species
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
                      characteristic_species = "Species (â‰¥50% frequency)") %>%
    bold(part = "header") %>%
    autofit()
  
  save_as_docx(ft_characteristic, path = file.path(output_folder, "10_Characteristic_Species_Summary.docx"))
  
  # -----------------------------
  # STEP 7: Spatial Maps (if spatial data provided)
  # -----------------------------
  if(!is.null(rec2_rivers) && !is.null(nga_awa)) {
    cat("Step 7: Creating spatial maps...\n")
    
    # Setup map extent
    sampling_buffer <- st_buffer(community_sf, dist = buffer_distance)
    sampling_extent_nztm <- st_bbox(st_union(sampling_buffer))
    sampling_extent_sf <- st_as_sfc(sampling_extent_nztm) %>% st_transform(4326)
    sampling_extent <- st_bbox(sampling_extent_sf)
    
    # Crop spatial layers
    nga_awa_crop <- st_crop(nga_awa, sampling_extent_nztm)
    rec2_trimmed <- st_intersection(rec2_rivers, nga_awa_crop)
    rec2_sampling <- st_crop(rec2_trimmed, sampling_extent_nztm)
    
    # Prepare map data
    map_data <- merged_matrix %>%
      left_join(scores_3d %>% select(Site, Cluster), by = c("Representative_Site" = "Site"))
    
    num_clusters <- length(unique(na.omit(map_data$Cluster)))
    palette_colors <- scales::hue_pal()(num_clusters)
    
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
    
    # Combined NMDS map
    nmds_combined <- map_nmds1 + map_nmds2 +
      plot_annotation(title = paste("NMDS Ordination of", community_name, "Communities"),
                      subtitle = paste("Jaccard distance | Stress =", stress_val),
                      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                                    plot.subtitle = element_text(size = 12, hjust = 0.5)))
    
    ggsave(file.path(output_folder, "11_NMDS_Spatial_Gradients.png"), nmds_combined, 
           width = 14, height = 6, dpi = 300, bg = "white")
    
    # Map: Sampling Locations
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
    
    ggsave(file.path(output_folder, "12_Sampling_Locations_Map.png"), sampling_map, 
           width = 8, height = 7, dpi = 300, bg = "white")
    
    # Map: Cluster Map
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
           subtitle = paste("PAM clustering (k =", optimal_k, ") based on Jaccard distance"),
           x = "Longitude", y = "Latitude")
    
    ggsave(file.path(output_folder, "13_Cluster_Map.png"), cluster_map, 
           width = 9, height = 7, dpi = 300, bg = "white")
    
    cat("  - 3 spatial maps created\n")
  }
  
  # -----------------------------
  # Completion Message
  # -----------------------------
  cat("\n", rep("=", 70), "\n", sep = "")
  cat(paste0(toupper(community_name), " ANALYSIS COMPLETE!\n"))
  cat(paste0("All outputs saved to '", output_folder, "' folder\n"))
  cat(rep("=", 70), "\n\n", sep = "")
  
  # Return important objects for further use
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

# -----------------------------
# 2. DEFINE COMMUNITY FILTERS
# -----------------------------

# Fish families
fish_families <- c("Retropinnidae", "Galaxiidae", "Anguillidae", "Cheimarrichthyidae",
                   "Eleotridae", "Mugilidae", "Tripterygiidae", "Pleuronectidae", 
                   "Microdesmidae", "Gobiidae", "Ictaluridae", "Cyprinidae", 
                   "Poeciliidae", "Salmonidae", "Percidae")

# Macroinvertebrate families
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

#Macrophyte families
macrophyte_families <- c("Characeae", "Hydatellaceae", "Haloragaceae", "Potamogetonaceae", "Hydrocharitaceae",
                         "Nymphaeaceae", "Cyperaceae", "Juncaceae", "Typhaceae", "Polygonaceae")

#Diatom phylum
diatom_phyla <- c("Bacillariophyta")

#Ciliate phylum
ciliate_phyla <- c("Ciliophora")

#Rotifer phylum
rotifer_phyla <- c("Rotifera")

# -----------------------------
# 3. LOAD SPATIAL DATA (once for all analyses)
# -----------------------------
cat("\n########## LOADING SPATIAL DATA ##########\n")
rec2_rivers <- st_read("Data/REC2_Layers/River_Lines.shp") %>% st_transform(2193)
nga_awa <- st_read("Data/Nga Awa shapefiles/DOC_NgāAwa_RiverSites_20250122_n14.shp") %>% st_transform(2193)
cat("Spatial layers loaded successfully\n")

# -----------------------------
# 4. RUN ANALYSES
# -----------------------------

