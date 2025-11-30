
library(shiny)
library(data.table)
library(shinycssloaders)
library(leaflet)

# Load and convert to data.table
merged <- as.data.table(readRDS("records.rds"))

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
      downloadButton("downloadData", "Download Filtered Data")
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
}

shinyApp(ui = ui, server = server)
