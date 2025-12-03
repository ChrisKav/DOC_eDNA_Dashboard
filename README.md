# Wilderlab eDNA Processing & Analysis

A comprehensive R-based pipeline for processing Wilderlab eDNA data, conducting community ecological analyses, and exploring biodiversity records interactively.

## Overview

This repository provides tools for:

1. **Data Preparation** - Fetch and harmonize eDNA records from Wilderlab API and public sources, enrich with taxonomy and spatial attributes
2. **Community Analysis** - Run standardized ecological analyses (NMDS ordination, clustering, indicator species) across multiple taxonomic groups and catchments
3. **Interactive Exploration** - Filter, visualize, and export biodiversity records through a Shiny web interface

---

## Features

### Data Processing (`get_API_data.R`)
- Fetches records from Wilderlab API (private DOC data) and public S3
- Enriches taxonomy using NCBI lineages via `insect` package
- Fuzzy matches species to NZ Threat Classification System (NZTCS)
- Performs spatial joins to assign Ngā Awa catchments and Regional Council boundaries
- Produces analysis-ready dataset: `Data/records.rds`

### Community Analysis (`community_analysis_functions.R`, `run_community_analysis.R`)
- Modular analysis framework supporting multiple community types:
  - Fish (15 families)
  - Macroinvertebrates (90+ families)
  - Macrophytes (10 families)
  - Diatoms, Ciliates, Rotifers (by phylum)
- Automated workflows for each Ngā Awa catchment:
  - NMDS ordination with species vectors
  - PAM clustering with silhouette optimization
  - Indicator species analysis (IndVal)
  - Species frequency heatmaps
  - Spatial distribution maps
- Generates publication-ready outputs (PNG figures, Word tables)

### Interactive Viewer (`ShinyApp.R`)
- Hierarchical taxonomic filters (phylum → species)
- Conservation status and data provenance filters
- Interactive map preview with clustering
- CSV export of filtered records
- GeoPackage export with taxa-specific layers

---

## Installation

### Prerequisites

**R packages:**
```r
install.packages(c(
  "data.table", "readr", "readxl", "stringr", "stringdist", 
  "sf", "insect", "wilderlab", "dplyr", "tidyverse", 
  "vegan", "ggplot2", "ggrepel", "ggnewscale", "cluster", 
  "factoextra", "ggspatial", "patchwork", "flextable", 
  "officer", "indicspecies", "pheatmap", "shiny", 
  "shinycssloaders", "leaflet"
))
```

**System dependencies for `sf` package:**
- **Ubuntu/Debian:** `sudo apt install libgdal-dev libgeos-dev libproj-dev`
- **macOS:** `brew install gdal geos proj`
- **Windows:** Binary packages usually include dependencies

### Setup

1. Clone repository:
   ```bash
   git clone https://github.com/yourusername/wilderlab-edna.git
   cd wilderlab-edna
   ```

2. Create required directories:
   ```bash
   mkdir -p Data Output logs
   ```

3. Add required data files (see [Data Requirements](#data-requirements))

---

## Data Requirements

Place these files in the `Data/` directory:

### Essential Files
- **NZTCS.xlsx** - NZ Threat Classification System spreadsheet (sheet: "Exported Data")
- **Shapefiles:**
  - `Nga Awa shapefiles/DOC_NgāAwa_RiverSites_20250122_n14.shp` (+ .shx, .dbf, .prj)
  - `REC2_Layers/River_Lines.shp` (River Environment Classification v2)
  - `Regional Council shapefiles/regional-council-2022-generalised.shp`

### API Credentials (Optional)

For accessing DOC private data, create `Data/wilder_keys.csv`:
```csv
Key,your_api_key_here
Secret,your_api_secret_here
xapikey,your_x_api_key_here
```

**Alternative:** Set environment variables:
```bash
export WILDER_KEY="your_api_key"
export WILDER_SECRET="your_api_secret"
export WILDER_XAPIKEY="your_x_api_key"
```

*Note: Public S3 data will be used if credentials are unavailable*

---

## Usage

### 1. Data Preparation

Fetch latest records and create `Data/records.rds`:

```bash
Rscript get_API_data.R
```

**Expected columns in output:**
- Sample metadata: `Sample Name`, `Date`, `Latitude`, `Longitude`
- Taxonomy: `species`, `genus`, `family`, `order`, `class`, `phylum`
- Conservation: `Threat Status`, `Threat Category`
- Spatial: `Nga Awa`, `Regional Council`
- Provenance: `DOC Data` (Yes/No)

### 2. Community Analysis

Run analyses for all Ngā Awa catchments:

```bash
Rscript run_community_analysis.R
```

**For detailed logging:**
```bash
Rscript run_per_catchment_logged.R
```

**Debug single catchment:**
```bash
Rscript run_per_catchment_logged.R "Arahura River"
```

Outputs saved to: `Output/<Catchment>/<CommunityType>/`

### 3. Shiny Application

Launch interactive viewer:

```bash
Rscript -e "shiny::runApp('ShinyApp.R', port = 8100, launch.browser = TRUE)"
```

Or in RStudio: Open `ShinyApp.R` → Click "Run App"

---

## Output Structure

```
Output/
├── Arahura_River/
│   ├── Fish/
│   │   ├── 01_Sampling_Effort_by_Year.png
│   │   ├── 02_Threatened_Species_Summary.png
│   │   ├── 03_NMDS_Ordination_Plot.png
│   │   ├── 04_Silhouette_Analysis.png
│   │   ├── 05_Dendrogram.png
│   │   ├── 06_Cluster_Summary.docx
│   │   ├── 07_Species_Frequency_by_Cluster.docx
│   │   ├── 08_Indicator_Species_Analysis.docx
│   │   ├── 09_Species_Frequency_Heatmap.png
│   │   ├── 10_Characteristic_Species_Summary.docx
│   │   ├── 11_NMDS_Spatial_Gradients.png
│   │   ├── 12_Sampling_Locations_Map.png
│   │   └── 13_Cluster_Map.png
│   ├── Macroinvertebrates/
│   ├── Macrophytes/
│   └── ...
├── Clutha_River/
└── ...

logs/
└── run_20250203_143022.log
```

---

## Workflow Examples

### Complete Analysis Pipeline

```r
# 1. Fetch and prepare data
source("get_API_data.R")

# 2. Run all catchment analyses
source("run_community_analysis.R")

# 3. Launch Shiny app to explore results
shiny::runApp("ShinyApp.R")
```

### Custom Analysis for Single Community

```r
source("community_analysis_functions.R")

# Load data
df <- readRDS("Data/records.rds")
df_arahura <- df[df$`Nga Awa` == "Arahura River", ]

# Filter to fish
fish_data <- df_arahura[df_arahura$family %in% fish_families, ]

# Run analysis
results <- analyze_community(
  df = df_arahura,
  community_data = fish_data,
  community_name = "Fish",
  output_folder = "Output/Custom/Arahura_Fish",
  rec2_rivers = NULL,  # Set to NULL to skip maps
  nga_awa = NULL
)
```

---

## Troubleshooting

### Common Issues

**Missing `records.rds`:**
```bash
# Run data preparation first
Rscript get_API_data.R
```

**Date parsing warnings:**
- Ensure `Date` column uses ISO format (YYYY-MM-DD)
- Script handles DD/MM/YYYY and MM/DD/YYYY as fallbacks

**GDAL/sf installation errors:**
- Ubuntu: `sudo apt install libgdal-dev libgeos-dev libproj-dev`
- macOS: `brew install gdal`
- Verify: `sf::sf_extSoftVersion()`

**Empty community analysis:**
- Check species column names match filter lists
- Verify catchment name in `Nga Awa` column exactly matches script input
- Inspect filtering: `sum(df$family %in% fish_families)`

**Memory issues (large catchments):**
- Process catchments individually using `run_per_catchment_logged.R`
- Increase R memory limit: `options(java.parameters = "-Xmx8g")`

### Debugging Tools

**Check data structure:**
```r
df <- readRDS("Data/records.rds")
str(df)
table(df$`Nga Awa`)
```

**Test single community:**
```r
source("community_analysis_functions.R")
df <- readRDS("Data/records.rds")
fish <- df[df$family %in% fish_families, ]
nrow(fish)  # Should be > 0
```

**Review logs:**
```bash
cat logs/run_*.log | grep ERROR
```

---

## Configuration

### Customizing Community Filters

Edit `community_analysis_functions.R` to modify taxonomic groups:

```r
# Add custom fish family
fish_families <- c(fish_families, "Myxinidae")

# Define new community
plankton_orders <- c("Calanoida", "Cyclopoida", "Harpacticoida")
```

### Adjusting Analysis Parameters

In `analyze_community()` function calls:

```r
analyze_community(
  ...,
  grouping_distance = 200,   # Spatial clustering threshold (m)
  buffer_distance = 5000     # Map buffer around sites (m)
)
```

---

## Citation

If using this pipeline in publications, please cite:

- Wilderlab eDNA platform: [https://wilderlab.co.nz](https://wilderlab.co.nz)
- R packages: `vegan`, `sf`, `indicspecies` (see `citation("package_name")`)

---

## Contributing

Contributions welcome! Please:
1. Fork repository
2. Create feature branch (`git checkout -b feature/new-analysis`)
3. Commit changes (`git commit -m 'Add new analysis type'`)
4. Push to branch (`git push origin feature/new-analysis`)
5. Open Pull Request

---

## License

[Specify license - e.g., MIT, GPL-3.0]

---

## Contact

For questions or issues:
- Open a GitHub issue
- Email: [your.email@domain.com]

---

## Acknowledgments

- Department of Conservation (DOC) for data access
- Wilderlab for eDNA processing and API
- REC2 and LINZ for spatial datasets
- NZ Threat Classification System contributors