# Wilderlab eDNA Data Processing & Shiny App

This repository contains two main R scripts:

- `get_API_data.R` — Pulls data from the Wilderlab API and public S3, cleans and enriches the data (taxonomy, spatial joins, fuzzy matching to NZTCS), and saves a preprocessed RDS file.
- `ShinyApp.R` — A Shiny application that loads the preprocessed records RDS file, provides interactive filters, a map preview, and download options (CSV and a GeoPackage with one layer per taxon).

This README explains what each script does, required dependencies, configuration (API keys), how to run them, expected outputs and troubleshooting tips.

---

## Table of contents

- [Overview](#overview)
- [Files](#files)
- [Requirements](#requirements)
- [Configuration — API keys](#configuration---api-keys)
- [Running the data preparation script](#running-the-data-preparation-script)
- [Running the Shiny app](#running-the-shiny-app)
- [GeoPackage export (from the app)](#geopackage-export-from-the-app)
- [Troubleshooting & diagnostics](#troubleshooting--diagnostics)
- [Development notes & suggestions](#development-notes--suggestions)
- [License](#license)

---

## Overview

- `get_API_data.R` is intended to be run periodically (manually or in a scheduled job). It:
  - Reads Wilderlab credentials (from a CSV or environment variables).
  - Fetches jobs, samples, taxa and records from the Wilderlab API.
  - Loads Wilderlab public S3 samples/records.
  - Harmonises per-job data frames to avoid type conflicts.
  - Adds taxonomic lineages using `insect::get_lineage`.
  - Performs fuzzy matching of species names to NZTCS taxonomy.
  - Adds spatial attributes (Nga Awa catchment, Regional Council) using shapefiles.
  - Writes an output RDS to `Data/records_DDMMYY.RDS` (dated with ddmmyy, e.g. `records_301125.RDS`).

- `ShinyApp.R` is a Shiny app that:
  - Loads a preprocessed RDS (by default it loads `records.rds` in the working directory).
  - Provides cascading filters (phylum → class → order → family → genus → species and additional metadata fields).
  - Shows a preview table and Leaflet map.
  - Allows downloading the filtered results as CSV.
  - Provides a button to download a GeoPackage (GPKG) where each selected taxon (at the chosen taxonomic level) is written as a separate layer.

---

## Files

- `get_API_data.R` — main preprocessing script.
- `ShinyApp.R` — Shiny UI + server script.
- `Data/` — expected folder with:
  - `wilder_keys.csv` (recommended; see below).
  - `NTZCS.xlsx` and necessary shapefiles referenced by the script.
  - Output records files such as `records_DDMMYY.RDS` (created by `get_API_data.R`).
- `README.md` — this file.

---

## Requirements

R (tested with R 4.x) and the following CRAN packages:

- dplyr
- purrr
- readr
- readxl
- stringr
- stringdist
- sf
- insect
- wilderlab
- data.table
- shiny
- shinycssloaders
- leaflet

You can install them with:

```r
install.packages(c(
  "dplyr", "purrr", "readr", "readxl", "stringr", "stringdist",
  "sf", "insect", "wilderlab", "data.table", "shiny",
  "shinycssloaders", "leaflet"
))
```

Note: `sf` may require system libraries (GDAL/GEOS/PROJ). Consult the `sf` installation notes for your OS.

---

## Configuration — API keys

The data prep script reads Wilderlab credentials from a simple two-column CSV by default at:

```
Data/wilder_keys.csv
```

Expected layout (no header required — the script reads the first two columns):

Column 1 | Column 2
--- | ---
Key | your-key-value
Secret | your-secret-value
xapikey | your-x-api-key-value

Examples (CSV rows):
```
Key,xxxxxxxxxxxxxxxxxxxxxxxxx
Secret,yyyyyyyyyyyyyyyyyyyyyyyy
xapikey,zzzzzzzzzzzzzzzzzzzzzzzz
```

Alternative: set environment variables (useful for CI or servers)

- WILDER_KEYS_FILE — path to an alternative CSV,
- WILDER_KEY, WILDER_SECRET, WILDER_XAPIKEY (or WILDER_X_API_KEY).

The script will prefer the CSV file if it exists; otherwise it will fall back to environment variables and warn if no keys are found.

---

## Running the data preparation script

From an R session or command line:

R interactive:
```r
# from project root
source("get_API_data.R")
```

Command line with Rscript:
```bash
Rscript get_API_data.R
```

Notes:
- The script expects `Data/NZTCS.xlsx` and shapefiles under `Data/` (see script header for exact filenames).
- By default the script writes a dated RDS to `Data/records_DDMMYY.RDS` (ddmmyy format). Example: `Data/records_301125.RDS`.
- If you want the Shiny app to use the newly created dated file, either:
  - Copy/rename the dated file to the name the app expects (`records.rds`), or
  - Edit `ShinyApp.R` to load the latest `records_*.RDS` file (a small change to the `readRDS()` call).

Example to create a `records.rds` symlink (on *nix/macOS):
```bash
# find newest and copy
cp Data/records_$(date +%d%m%y).RDS Data/records.rds
# or create a symlink
ln -sf Data/records_$(date +%d%m%y).RDS Data/records.rds
```

Windows PowerShell example to copy:
```powershell
Copy-Item -Path Data/records_$(Get-Date -Format ddMMyy).RDS -Destination Data/records.rds -Force
```

---

## Running the Shiny app

Start the app from command line:

```bash
Rscript -e "shiny::runApp('ShinyApp.R', port = 8100, launch.browser = TRUE)"
```

Or from within R/RStudio:
```r
source("ShinyApp.R")
# or
shiny::runApp("ShinyApp.R")
```

Important: `ShinyApp.R` expects an RDS file (default `records.rds`) in the working directory. Ensure `Data/records.rds` or `records.rds` is present or change the `readRDS()` path in `ShinyApp.R`.

---

## GeoPackage export (from the app)

The app provides a "Download Geopackage (per-taxon layers)" button and a selector `Create GDB layers grouped by:`. Choose the grouping level (species/genus/family/Latin name etc.) and click the button — the app will:

- Group the filtered records by the selected taxon column.
- Create one layer per taxon in a GeoPackage (*.gpkg).
- If coordinates are present for records, the layer will be spatial (geometry column using Latitude/Longitude, EPSG:4326). Non-spatial rows (missing coords) are written as attribute tables.
- Layer names are sanitized for safety (ASCII-only, underscores, trimmed).

Notes / limitations:
- GeoPackage layer name length and character rules differ by drivers — the script truncates/sanitizes layer names and will disambiguate duplicates by appending suffixes.
- Large numbers of taxon groups will result in long write times and a large file; use the filters to reduce the export if necessary.

---

## Troubleshooting & diagnostics

1. curl / connection issues when fetching API data:
   - If you see errors like `Recv failure: Connection was reset`, try:
     - Test with `curl -v https://connect.wilderlab.co.nz` from a terminal to isolate network/TLS issues.
     - Check proxy/VPN/firewall settings — set `http_proxy` / `https_proxy` in the environment if necessary.
     - Increase timeouts or retry (you can run the script later if the server is transiently unavailable).

2. Type mismatch when combining per-job data:
   - The script harmonises types (logical/factor → character) before binding rows. If you add new fields upstream, you may need to extend the harmonisation.

3. sf / GDAL errors when writing GeoPackage:
   - Ensure `sf` has been installed with appropriate system dependencies (GDAL/PROJ/GEOS).
   - If `st_write()` errors with driver issues, check your GDAL/PROJ versions and driver availability.

4. App cannot find `records.rds`:
   - Ensure you either copy the dated `records_DDMMYY.RDS` to `records.rds`, change the `readRDS()` path in `ShinyApp.R`, or update the working directory.

---

## Development notes & suggestions

- Consider adding a small step in `get_API_data.R` to always write an unversioned `Data/records.rds` (in addition to the dated file) — this will simplify running `ShinyApp.R` without changes.
- Add logging (e.g., `logger` package) and caching for API calls to avoid repeated API traffic during development.
- Add unit tests for core transformations (fuzzy matching, spatial joins).
- For very large datasets, consider switching heavy grouping operations to `data.table` or using a database backend (PostGIS) for spatial joins and exports.
- If you need to share credentials securely in CI, prefer environment variables or secret stores rather than checked-in CSV files.

---

## Contributing

Contributions welcome. Please:
- Open an issue to describe the change or bug.
- Send a PR with clear commit messages and tests (if applicable).

---

## License

Choose an appropriate license for your project (e.g. MIT). If you want, add a `LICENSE` file to this repo and note it here.

---

If you want, I can:
- Add a short example GitHub Actions workflow for running `get_API_data.R` nightly and storing the dated RDS as an artifact,
- Modify `get_API_data.R` to also write `Data/records.rds` (unversioned) automatically,
- Provide a sample `wilder_keys.csv` template file to include in the repo (but DO NOT commit real keys).

Which follow-up would you prefer?  