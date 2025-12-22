library(officer)
library(magrittr)
library(fs)
library(docxtractr)
library(flextable)

# Set parent directory that contains all the river folders
parent_dir <- "Output"  # Change appropriately

# List all river catchments
catchments <- dir_ls(parent_dir, type = "directory")

# Find all taxa (assuming all catchments have the same taxa structure)
taxa <- unique(unlist(lapply(catchments, function(ct) dir_ls(ct, type="directory", recurse=FALSE))))
taxa <- basename(taxa)

# Loop for each taxon to create a summary docx
for (taxon in taxa) {
  doc <- officer::read_docx()
  doc <- doc %>% body_add_par(paste("Figures and Tables for", taxon), style = "heading 1")
  
  for (catchment in catchments) {
    river_name <- basename(catchment)
    taxon_dir <- file.path(catchment, taxon)
    if (dir_exists(taxon_dir)) {
      figs <- dir_ls(taxon_dir, regexp = "\\.png$", type = "file")
      tbls <- dir_ls(taxon_dir, regexp = "\\.docx$", type = "file")
      
      if (length(figs) + length(tbls) == 0) next  # Skip empty catchments
      
      doc <- doc %>% body_add_par(river_name, style = "heading 2")
      
      # Add Figures
      if (length(figs) > 0) {
        doc <- doc %>% body_add_par("Figures:", style = "heading 3")
        for (fig in figs) {
          doc <- doc %>%
            body_add_par(basename(fig), style = "Normal") %>%
            body_add_img(src = fig, width = 5, height = 4)
        }
      }
      # Add Tables
      if (length(tbls) > 0) {
        doc <- doc %>% body_add_par("Tables:", style = "heading 3")
        for (tbl in tbls) {
          docx_in <- docxtractr::read_docx(tbl)
          tbl_list <- docxtractr::docx_extract_all_tbls(docx_in)
          for (tab in tbl_list) {
            doc <- doc %>% body_add_flextable(flextable::flextable(tab))
          }
        }
      }
      doc <- doc %>% body_add_par("", style = "Normal")
    }
  }
  # Save the summary file for the taxon
  print(doc, target = file.path(parent_dir, paste0(taxon, "_summary.docx")))
}