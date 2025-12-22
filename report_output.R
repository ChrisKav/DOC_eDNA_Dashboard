library(officer)
library(magrittr)
library(fs)
library(docxtractr)
library(flextable)

# Configuration ----
PARENT_DIR <- "Output"
IMG_WIDTH <- 5
IMG_HEIGHT <- 4

# Helper Functions ----

#' Get unique taxa across all catchments
#' @param catchments Character vector of catchment paths
#' @return Character vector of unique taxa names
get_taxa <- function(catchments) {
  taxa_dirs <- lapply(catchments, function(ct) {
    dir_ls(ct, type = "directory", recurse = FALSE)
  })
  unique(basename(unlist(taxa_dirs)))
}

#' Add figures to document
#' @param doc officer document object
#' @param fig_paths Character vector of figure file paths
#' @return Modified officer document
add_figures <- function(doc, fig_paths) {
  if (length(fig_paths) == 0) return(doc)
  
  doc <- doc %>% body_add_par("Figures:", style = "heading 3")
  
  for (fig in fig_paths) {
    doc <- doc %>%
      body_add_par(basename(fig), style = "Normal") %>%
      body_add_img(src = fig, width = IMG_WIDTH, height = IMG_HEIGHT)
  }
  doc
}

#' Add tables from docx files to document
#' @param doc officer document object
#' @param tbl_paths Character vector of table file paths
#' @return Modified officer document
add_tables <- function(doc, tbl_paths) {
  if (length(tbl_paths) == 0) return(doc)
  
  doc <- doc %>% body_add_par("Tables:", style = "heading 3")
  
  for (tbl_path in tbl_paths) {
    tryCatch({
      docx_in <- docxtractr::read_docx(tbl_path)
      tbl_list <- docxtractr::docx_extract_all_tbls(docx_in)
      
      for (tab in tbl_list) {
        doc <- doc %>% body_add_flextable(flextable(tab))
      }
    }, error = function(e) {
      warning(sprintf("Failed to extract table from %s: %s", 
                      basename(tbl_path), e$message))
    })
  }
  doc
}

#' Add catchment content to document
#' @param doc officer document object
#' @param catchment Path to catchment directory
#' @param taxon Taxon name
#' @return Modified officer document or NULL if no content
add_catchment_content <- function(doc, catchment, taxon) {
  river_name <- basename(catchment)
  taxon_dir <- path(catchment, taxon)
  
  if (!dir_exists(taxon_dir)) return(NULL)
  
  # Find figures and tables
  figs <- dir_ls(taxon_dir, regexp = "\\.png$", type = "file")
  tbls <- dir_ls(taxon_dir, regexp = "\\.docx$", type = "file")
  
  # Skip if no content
  if (length(figs) + length(tbls) == 0) return(NULL)
  
  # Add catchment heading and content
  doc <- doc %>%
    body_add_par(river_name, style = "heading 2") %>%
    add_figures(figs) %>%
    add_tables(tbls) %>%
    body_add_par("", style = "Normal")
  
  doc
}

# Main Script ----

# Validate parent directory
if (!dir_exists(PARENT_DIR)) {
  stop(sprintf("Parent directory '%s' does not exist", PARENT_DIR))
}

# Get all catchments
catchments <- dir_ls(PARENT_DIR, type = "directory")

if (length(catchments) == 0) {
  stop(sprintf("No catchment directories found in '%s'", PARENT_DIR))
}

message(sprintf("Found %d catchments", length(catchments)))

# Get all taxa
taxa <- get_taxa(catchments)

if (length(taxa) == 0) {
  stop("No taxa directories found in catchments")
}

message(sprintf("Found %d taxa: %s", length(taxa), paste(taxa, collapse = ", ")))

# Generate summary document for each taxon
for (taxon in taxa) {
  message(sprintf("Processing taxon: %s", taxon))
  
  doc <- officer::read_docx() %>%
    body_add_par(sprintf("Figures and Tables for %s", taxon), style = "heading 1")
  
  catchment_count <- 0
  
  for (catchment in catchments) {
    result <- add_catchment_content(doc, catchment, taxon)
    if (!is.null(result)) {
      doc <- result
      catchment_count <- catchment_count + 1
    }
  }
  
  # Save document
  output_file <- path(PARENT_DIR, sprintf("%s_summary.docx", taxon))
  print(doc, target = output_file)
  
  message(sprintf("  Created %s with content from %d catchments", 
                  basename(output_file), catchment_count))
}

message("Report generation complete!")