# Ensure exported column names match what the Shiny app expects.

# Map from current summary_df column names to the Shiny-friendly names used by the app
rename_map <- c(
  "Name" = "Latin Name",
  "CommonName" = "Common Name",
  "ClientSampleID" = "Sample Name",
  "sum_count" = "Total Reads",
  "mean_count" = "Average reads",
  "unique_UID_count" = "Hit number",
  "total_UID" = "Total hits",
  "Rank" = "Taxonomic Rank",
  "Group" = "Taxon Group",
  "Status" = "Threat Status",
  "Category" = "Threat Category",
  "CollectionDate" = "Date",
  "ThreatReport" = "Threat Document",
  "CollectedBy" = "Collector",
  "DOC_Data" = "DOC Data",
  "MakeDataPublic" = "Public/Private",
  "Nga_Awa_Catchment" = "Nga Awa",
  "Regional_Council" = "Regional Council",
  "Wilderlab_Sp_name" = "Wilderlab Sp Name",
  "NZTC_Sp_name" = "NZTC Sp Name",
  "Report" = "Wilderlab Report"
)

# Only rename columns that actually exist in summary_df
existing_old_names <- intersect(names(rename_map), names(summary_df))
if (length(existing_old_names) > 0) {
  new_names <- rename_map[existing_old_names]
  # perform rename in place
  names(summary_df)[match(existing_old_names, names(summary_df))] <- unname(new_names)
}

# Ensure the Shiny app's expected columns exist (create missing ones as NA so filters won't error)
expected_cols_for_shiny <- c(
  "Latin Name","Common Name","Sample Name","Total Reads","Average reads","Hit number","Total hits",
  "Taxonomic Rank","Taxon Group","Threat Status","Threat Category","Date","Threat Document","Collector",
  "DOC Data","Public/Private","Nga Awa","Regional Council","Wilderlab Sp Name","NZTC Sp Name",
  "TaxID","Wilderlab Report","Latitude","Longitude"
)
missing_cols <- setdiff(expected_cols_for_shiny, names(summary_df))
if (length(missing_cols) > 0) {
  for (c in missing_cols) summary_df[[c]] <- NA
}

# Optionally re-order columns to put the most-used columns first (not required, but keeps parity with original)
preferred_order <- intersect(expected_cols_for_shiny, names(summary_df))
remaining <- setdiff(names(summary_df), preferred_order)
summary_df <- summary_df[, c(preferred_order, remaining)]
# --- End insertion block ---