library(shiny)
library(data.table)
library(shinycssloaders)
library(leaflet)
library(sf)            # for spatial write
library(dplyr)         # used in GPKG creation
library(stringr)       # used to sanitize layer names

# Load and convert to data.table
merged <- as.data.table(readRDS("Data/records.rds"))

# Helper function
unique_sorted <- function(x) sort(unique(na.omit(x)))

# UI
ui <- fluidPage(
  titlePanel("Biodiversity Records Filter"),
  sidebarLayout(
    sidebarPanel(
      selectInput("phylum", "Phylum", choices = c("All", unique_sorted(merged$phylum))),
      selectInput("class", "Class", choices = "All"),
      selectInput("order", "Order", choices = "All"),
      selectInput("family", "Family", choices = "All"),
      selectInput("genus", "Genus", choices = "All"),
      selectInput("species", "Species", choices = "All"),
      selectInput("Threat Status", "Threat Status", choices = c("All", unique_sorted(merged$`Threat Status`))),
      selectInput("Threat Category", "Threat Category", choices = c("All", unique_sorted(merged$`Threat Category`))),
      selectInput("Taxon Group", "Taxonomic Group", choices = c("All", unique_sorted(merged$`Taxon Group`))),
      selectInput("Taxonomic Rank", "Taxonomic Rank", choices = c("All", unique_sorted(merged$`Taxonomic Rank`))),
      selectInput("Threat Document", "Threat Classification Document", choices = c("All", unique_sorted(merged$`Threat Document`))),
      selectInput("DOC Data", "DOC's Data", choices = c("All", unique_sorted(merged$`DOC Data`))),
      selectInput("Public/Private", "Public/Private", choices = c("All", unique_sorted(merged$`Public/Private`))),
      selectInput("Nga Awa", "Nga Awa catchment", choices = c("All", unique_sorted(merged$`Nga Awa`))),
      selectInput("Regional Council", "Region", choices = c("All", unique_sorted(merged$`Regional Council`))),
      actionButton("resetFilters", "Reset Filters"),
      hr(),
      downloadButton("downloadData", "Download Filtered Data"),
      br(), br(),
      # New controls for Geodatabase export
      selectInput(
        "gdb_level",
        "Create GDB layers grouped by:",
        choices = c(
          "Species" = "species",
          "Genus" = "genus",
          "Family" = "family",
          "Latin Name" = "Latin Name",
          "Wilderlab Sp Name" = "Wilderlab Sp Name"
        ),
        selected = "species"
      ),
      downloadButton("downloadGDB", "Download Geopackage (per-taxon layers)")
    ),
    mainPanel(
      h4("Filtered data preview"),
      withSpinner(dataTableOutput("dataPreview"), type = 4, color = "#337ab7"),
      hr(),
      h4("Map Preview"),
      leafletOutput("mapPreview", height = 500),
      hr(),
      h4("Summary"),
      verbatimTextOutput("summaryPanel")
    )
  )
)

# Server
server <- function(input, output, session) {
  observeEvent(input$resetFilters, {
    updateSelectInput(session, "phylum", selected = "All")
    updateSelectInput(session, "class", choices = "All", selected = "All")
    updateSelectInput(session, "order", choices = "All", selected = "All")
    updateSelectInput(session, "family", choices = "All", selected = "All")
    updateSelectInput(session, "genus", choices = "All", selected = "All")
    updateSelectInput(session, "species", choices = "All", selected = "All")
    updateSelectInput(session, "Threat Status", selected = "All")
    updateSelectInput(session, "Threat Category", selected = "All")
    updateSelectInput(session, "Taxon Group", selected = "All")
    updateSelectInput(session, "Taxonomic Rank", selected = "All")
    updateSelectInput(session, "Threat Document", selected = "All")
    updateSelectInput(session, "DOC Data", selected = "All")
    updateSelectInput(session, "Public/Private", selected = "All")
    updateSelectInput(session, "Nga Awa", selected = "All")
    updateSelectInput(session, "Regional Council", selected = "All")
  })
  
  observeEvent(input$phylum, {
    classes <- unique_sorted(merged[phylum == input$phylum | input$phylum == "All", class])
    updateSelectInput(session, "class", choices = c("All", classes), selected = "All")
    updateSelectInput(session, "order", choices = "All")
    updateSelectInput(session, "family", choices = "All")
    updateSelectInput(session, "genus", choices = "All")
    updateSelectInput(session, "species", choices = "All")
  })
  
  observeEvent(input$class, {
    orders <- unique_sorted(merged[
      (phylum == input$phylum | input$phylum == "All") &
        (class == input$class | input$class == "All"), order])
    updateSelectInput(session, "order", choices = c("All", orders), selected = "All")
    updateSelectInput(session, "family", choices = "All")
    updateSelectInput(session, "genus", choices = "All")
    updateSelectInput(session, "species", choices = "All")
  })
  
  observeEvent(input$order, {
    families <- unique_sorted(merged[
      (phylum == input$phylum | input$phylum == "All") &
        (class == input$class | input$class == "All") &
        (order == input$order | input$order == "All"), family])
    updateSelectInput(session, "family", choices = c("All", families), selected = "All")
    updateSelectInput(session, "genus", choices = "All")
    updateSelectInput(session, "species", choices = "All")
  })
  
  observeEvent(input$family, {
    genera <- unique_sorted(merged[
      (phylum == input$phylum | input$phylum == "All") &
        (class == input$class | input$class == "All") &
        (order == input$order | input$order == "All") &
        (family == input$family | input$family == "All"), genus])
    updateSelectInput(session, "genus", choices = c("All", genera), selected = "All")
    updateSelectInput(session, "species", choices = "All")
  })
  
  observeEvent(input$genus, {
    spp <- unique_sorted(merged[
      (phylum == input$phylum | input$phylum == "All") &
        (class == input$class | input$class == "All") &
        (order == input$order | input$order == "All") &
        (family == input$family | input$family == "All") &
        (genus == input$genus | input$genus == "All"), species])
    updateSelectInput(session, "species", choices = c("All", spp), selected = "All")
  })
  
  filtered_data <- reactiveVal()
  
  get_filtered_data <- reactive({
    dt <- copy(merged)
    filters <- list(
      phylum = input$phylum,
      class = input$class,
      order = input$order,
      family = input$family,
      genus = input$genus,
      species = input$species,
      `Threat Status` = input$`Threat Status`,
      `Threat Category` = input$`Threat Category`,
      `Threat Document` = input$`Threat Document`,
      `Taxon Group` = input$`Taxon Group`,
      `Taxonomic Rank` = input$`Taxonomic Rank`,
      `DOC Data` = input$`DOC Data`,
      `Public/Private` = input$`Public/Private`,
      `Nga Awa` = input$`Nga Awa`,
      `Regional Council` = input$`Regional Council`
    )
    for (col in names(filters)) {
      val <- filters[[col]]
      if (!is.null(val) && val != "All") {
        dt <- dt[get(col) == val]
      }
    }
    dt
  })
  
  observe({
    filtered_data(get_filtered_data())
  })
  
  output$dataPreview <- renderDataTable({
    req(filtered_data())
  })
  
  output$summaryPanel <- renderText({
    dt <- filtered_data()
    record_count <- nrow(dt)
    taxid_count <- if ("TaxID" %in% names(dt)) uniqueN(dt$TaxID) else NA
    sample_count <- if ("Sample Name" %in% names(dt)) uniqueN(dt$`Sample Name`) else NA
    paste(
      "Number of records:", record_count,
      "
Unique TaxIDs:", taxid_count,
      "
Unique Sample Names:", sample_count
    )
  })
  
  output$mapPreview <- renderLeaflet({
    dt <- filtered_data()
    req(dt)
    
    if (!("Latitude" %in% names(dt)) || !("Longitude" %in% names(dt))) return(NULL)
    
    dt <- dt[!is.na(Latitude) & !is.na(Longitude)]
    
    leaflet(dt) %>%
      addTiles() %>%
      addMarkers(
        ~Longitude, ~Latitude,
        popup = ~paste("Sample:", `Sample Name`, "<br>Species:", `Latin Name`),
        clusterOptions = markerClusterOptions()
      )
  })
  
  output$downloadData <- downloadHandler(
    filename = function() paste0("filtered_data_", Sys.Date(), ".csv"),
    content = function(file) {
      fwrite(filtered_data(), file)
    }
  )
  
  # -------------------------------
  # New: Download Geopackage (GPKG) with one layer per taxa (by chosen taxonomic level)
  # -------------------------------
  sanitize_layer_name <- function(x) {
    # Replace non-alnum with underscore, collapse multiple underscores, trim length
    name <- iconv(as.character(x), to = "ASCII//TRANSLIT")
    name <- str_replace_all(name, "[^A-Za-z0-9_]", "_")
    name <- str_replace_all(name, "_+", "_")
    name <- str_trim(name, side = "both")
    name <- str_replace_all(name, "^_+|_+$", "")
    # enforce maximum reasonable length for layer names
    name <- substr(name, 1, 50)
    if (nchar(name) == 0) name <- "unknown"
    name
  }
  
  output$downloadGDB <- downloadHandler(
    filename = function() {
      paste0("filtered_gdb_", Sys.Date(), ".gpkg")
    },
    content = function(file) {
      dt <- copy(filtered_data())
      req(nrow(dt) >= 0) # allow zero rows but handle gracefully
      
      # Determine grouping column (use selected, but fall back if missing)
      group_col <- input$gdb_level
      available_cols <- names(dt)
      # if chosen column isn't present, try alternatives
      if (!(group_col %in% available_cols)) {
        # common alternatives tried in order
        alternatives <- c("Latin Name", "species", "Wilderlab Sp Name", "Wilderlab_Sp_name", "Genus", "genus", "Family", "family")
        found <- intersect(alternatives, available_cols)
        if (length(found) > 0) {
          group_col <- found[[1]]
        } else {
          # as last resort, create a single layer called "all_records"
          group_col <- NULL
        }
      }
      
      # ensure Longitude/Latitude exist
      has_coords <- all(c("Longitude", "Latitude") %in% names(dt))
      if (has_coords) {
        # check for numeric coords
        dt[, Longitude := as.numeric(Longitude)]
        dt[, Latitude := as.numeric(Latitude)]
        coords_valid <- any(!is.na(dt$Longitude) & !is.na(dt$Latitude))
      } else {
        coords_valid <- FALSE
      }
      
      # If no grouping column, write entire filtered table as single layer
      if (is.null(group_col)) {
        # write whole dataset as one layer
        if (coords_valid) {
          sf_all <- tryCatch(
            st_as_sf(dt[!is.na(Longitude) & !is.na(Latitude)], coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE),
            error = function(e) NULL
          )
          # write spatial subset (if any spatial) and non-spatial remainder as separate layers
          first <- TRUE
          if (!is.null(sf_all) && nrow(sf_all) > 0) {
            if (first) {
              st_write(sf_all, file, layer = "all_records_spatial", driver = "GPKG", delete_dsn = TRUE, quiet = TRUE)
              first <- FALSE
            } else {
              st_write(sf_all, file, layer = "all_records_spatial", driver = "GPKG", append = TRUE, quiet = TRUE)
            }
          }
          # non-spatial remainder
          non_spatial_dt <- dt[is.na(Longitude) | is.na(Latitude)]
          if (nrow(non_spatial_dt) > 0) {
            # write non-spatial table
            if (first) {
              st_write(non_spatial_dt, file, layer = "all_records_table", driver = "GPKG", delete_dsn = TRUE, quiet = TRUE)
              first <- FALSE
            } else {
              st_write(non_spatial_dt, file, layer = "all_records_table", driver = "GPKG", append = TRUE, quiet = TRUE)
            }
          }
        } else {
          # no coords at all, write single non-spatial layer
          st_write(dt, file, layer = "all_records_table", driver = "GPKG", delete_dsn = TRUE, quiet = TRUE)
        }
        return(NULL)
      }
      
      # Otherwise, group by chosen taxa level
      groups <- unique(na.omit(dt[[group_col]]))
      # If no non-NA groups exist, create a single "unknown" group
      if (length(groups) == 0) groups <- "unknown"
      
      first_write <- TRUE
      # Progress indicator
      withProgress(message = "Writing Geopackage...", value = 0, {
        n_groups <- length(groups)
        cnt <- 0
        for (g in groups) {
          cnt <- cnt + 1
          incProgress(1 / n_groups, detail = paste("Layer:", as.character(g)))
          subset_dt <- dt[get(group_col) == g]
          # create safe layer name (GPKG layer names should be short and safe)
          layer_name <- sanitize_layer_name(g)
          # avoid duplicate layer names by appending an index if necessary
          # track existing layers by a simple counter (if same sanitized name appears)
          # We'll append a suffix if multiple groups sanitize to same name
          # We can count previously used layer names in the file by keeping a vector
          if (!exists("used_layer_names", envir = session$userData)) {
            session$userData$used_layer_names <- character(0)
          }
          if (layer_name %in% session$userData$used_layer_names) {
            # append small suffix
            i <- 2
            new_name <- paste0(layer_name, "_", i)
            while (new_name %in% session$userData$used_layer_names) {
              i <- i + 1
              new_name <- paste0(layer_name, "_", i)
            }
            layer_name <- new_name
          }
          session$userData$used_layer_names <- c(session$userData$used_layer_names, layer_name)
          
          # Write layer: spatial if coords present & valid, otherwise non-spatial table
          if (coords_valid && any(!is.na(subset_dt$Longitude) & !is.na(subset_dt$Latitude))) {
            sf_obj <- tryCatch(
              st_as_sf(subset_dt[!is.na(Longitude) & !is.na(Latitude)], coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE),
              error = function(e) NULL
            )
            if (!is.null(sf_obj) && nrow(sf_obj) > 0) {
              if (first_write) {
                st_write(sf_obj, file, layer = layer_name, driver = "GPKG", delete_dsn = TRUE, quiet = TRUE)
                first_write <- FALSE
              } else {
                st_write(sf_obj, file, layer = layer_name, driver = "GPKG", append = TRUE, quiet = TRUE)
              }
            } else {
              # If conversion failed, fallback to non-spatial write
              if (first_write) {
                st_write(subset_dt, file, layer = layer_name, driver = "GPKG", delete_dsn = TRUE, quiet = TRUE)
                first_write <- FALSE
              } else {
                st_write(subset_dt, file, layer = layer_name, driver = "GPKG", append = TRUE, quiet = TRUE)
              }
            }
          } else {
            # no coordinates: write non-spatial table
            if (first_write) {
              st_write(subset_dt, file, layer = layer_name, driver = "GPKG", delete_dsn = TRUE, quiet = TRUE)
              first_write <- FALSE
            } else {
              st_write(subset_dt, file, layer = layer_name, driver = "GPKG", append = TRUE, quiet = TRUE)
            }
          }
        } # end for groups
      }) # end withProgress
    } # end content
  ) # end downloadHandler
}

shinyApp(ui = ui, server = server)