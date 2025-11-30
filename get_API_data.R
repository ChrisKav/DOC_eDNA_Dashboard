library(dplyr)
library(insect)
library(readr)
library(readxl)
library(stringr)
library(stringdist)
library(sf)
library(wilderlab)

today <- Sys.Date()# Get today's date

nztcs   <- read_excel("Data/NZTCS.xlsx", sheet = "Exported Data")

# Pull DOC data from API
key <- "AKIATVYXGCYLUO3ORORC"
secret <- "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"
xapikey <- "aiZhm77XKe7hk0KLOzEut9fNV3xRjOoc9IlWqEuz"

jobs <- get_wilderdata("jobs", key = key, secret = secret, xapikey = xapikey)
jobs

samples <- get_wilderdata("samples", key = key, secret = secret, xapikey = xapikey)
samples

taxa <- get_wilderdata("taxa", key = key, secret = secret, xapikey = xapikey)
head(taxa, 10)

records <- vector(mode = "list", length = nrow(jobs))
for(i in seq_along(records)){
  records[[i]] <- get_wilderdata("records", JobID = jobs$JobID[i],
                                 key = key, secret = secret, xapikey = xapikey)
}
records <- do.call("rbind", records)

#  Pull Wilderlab public dataset
public_samples <- read.csv("http://s3.ap-southeast-2.amazonaws.com/wilderlab.publicdata/samples.csv")
public_records <- read.csv("http://s3.ap-southeast-2.amazonaws.com/wilderlab.publicdata/records.csv")

#Merge DOC and public samples datasets

DOC_samples <- samples %>%
  mutate(MakeDataPublic = ifelse(MakeDataPublic == 1, "Private", "Public"))

DOC_samples <- DOC_samples %>%
  mutate(DOC_Data = "Yes")

public_samples <- public_samples %>%
  mutate(DOC_Data = "No")

public_samples <- public_samples %>%
  mutate(MakeDataPublic = "Public")

public_samples <- public_samples %>%
  filter(!UID %in% DOC_samples$UID)

DOC_samples <- DOC_samples %>%
  mutate(CollectionDate = as.Date(CollectionDate))

public_samples <- public_samples %>%
  mutate(CollectionDate = as.Date(CollectionDate))

DOC_samples$Latitude <- as.numeric(DOC_samples$Latitude)
DOC_samples$Longitude <- as.numeric(DOC_samples$Longitude)

public_samples$Latitude <- as.numeric(public_samples$Latitude)
public_samples$Longitude <- as.numeric(public_samples$Longitude)

all_samples <- bind_rows(DOC_samples, public_samples)

merged_samples <- all_samples %>%
  select(SID, UID, JobID, CollectionDate, ClientSampleID,
         Latitude, Longitude, EnvironmentType, TICI, TICIVersion,
         TICINoSeqs, TICIQuantile, MakeDataPublic, DOC_Data, Report)

#Merge DOC and public records datasets
public_records <- public_records %>%
  filter(!UID %in% records$UID)

merged_records <- bind_rows(records, public_records)

# Add taxonomy to records
tdb <- taxa[, 1:4]
colnames(tdb) <- c("taxID", "parent_taxID", "rank", "name")
lineages = insect::get_lineage(merged_records$TaxID, tdb)
merged_records$phylum <- sapply(lineages, "[", "phylum")
merged_records$class <- sapply(lineages, "[", "class")
merged_records$order <- sapply(lineages, "[", "order")
merged_records$family <- sapply(lineages, "[", "family")
merged_records$genus <- sapply(lineages, "[", "genus")
merged_records$species <- sapply(lineages, "[", "species")
merged_records$Latitude <- samples$Latitude[match(merged_records$UID, merged_samples$UID)]
merged_records$Longitude <- samples$Longitude[match(merged_records$UID, merged_samples$UID)]
merged_records$ClientSampleID <- samples$ClientSampleID[match(merged_records$UID, merged_samples$UID)]

# Run this once in a separate R script
clean_taxa <- function(df) {
  df %>%
    mutate(
      across(c(phylum, class, order, family, genus, species), norm_chr),
      # Normalise known synonyms to match NZTCS spelling
      class = dplyr::recode(class, !!!class_synonyms),
      species = dplyr::recode(species, !!!species_synonyms),
      # (Optional) Pretty capitalization for ranks (does not change species binomials)
      phylum = cap_case(phylum),
      class  = cap_case(class),
      order  = cap_case(order),
      family = cap_case(family),
      genus  = cap_case(genus)
    )
}

# -------------------------------
# Utilities
# -------------------------------
norm_chr <- function(x) ifelse(is.na(x), NA_character_, str_squish(as.character(x)))
cap_case <- function(x) ifelse(is.na(x), NA_character_, str_to_title(x))

class_synonyms <- c("Actinopteri" = "Actinopterygii")
species_synonyms <- c("Galaxias sp. D (Allibone et al., 1996)" = "Galaxias \"species D\"")

# -------------------------------
# Load and clean as before
# -------------------------------

# Clean eDNA records
clean_taxa <- function(df) {
  df %>%
    mutate(
      across(c(phylum, class, order, family, genus, species), norm_chr),
      class = recode(class, !!!class_synonyms),
      species = recode(species, !!!species_synonyms),
      phylum = cap_case(phylum),
      class  = cap_case(class),
      order  = cap_case(order),
      family = cap_case(family),
      genus  = cap_case(genus)
    )
}
records <- clean_taxa(merged_records)

# Prepare NZTCS
nztcs_taxa <- nztcs %>%
  mutate(
    species_nztcs = norm_chr(`Current Species Name`),
    Genus  = cap_case(norm_chr(Genus)),
    Family = cap_case(norm_chr(Family)),
    Order  = cap_case(norm_chr(Order)),
    Class  = cap_case(norm_chr(Class)),
    Phylum = cap_case(norm_chr(Phylum)),
    Status = norm_chr(Status),
    Category = norm_chr(Category),
    BioStatus = norm_chr(`Bio Status`),
    YearAssessed = suppressWarnings(as.integer(`Year Assessed`)),
    ThreatReport = norm_chr(`Report Name`)
  ) %>%
  select(species_nztcs, Genus, Family, Order, Class, Phylum, Status, Category, BioStatus, ThreatReport, YearAssessed)

nztcs_sp <- nztcs_taxa %>%
  arrange(desc(YearAssessed)) %>%
  distinct(species_nztcs, .keep_all = TRUE) %>%
  select(species_nztcs, Status, Category, BioStatus, ThreatReport, Genus, Family, Order, Class, Phylum)

# -------------------------------
# Fuzzy match species names
# -------------------------------

clean_species <- function(x) {
  x <- tolower(x)
  x <- str_replace_all(x, '\\(.*\\)', '')         # Remove anything in parentheses
  x <- str_replace_all(x, '[[:punct:]]', '')      # Remove all punctuation
  x <- str_replace_all(x, 'sp\\.?', '')           # Remove 'sp.'
  x <- str_replace_all(x, 'species', '')          # Remove species
  x <- str_replace_all(x, '\\s+', ' ')            # Normalize whitespace
  x <- trimws(x)                                  # Trim leading/trailing spaces
  return(x)
}

records <- records %>% mutate(species_clean = clean_species(species))
nztcs_sp <- nztcs_sp %>% mutate(species_clean = clean_species(species_nztcs))

fuzzy_match <- function(name, choices, max_dist = 0.12) {
  distances <- stringdist(name, choices, method = "jw")
  best <- which.min(distances)
  if (length(best) == 0 || distances[best] > max_dist) return(NA)
  return(choices[best])
}

records$matched_species <- sapply(records$species_clean, function(x) fuzzy_match(x, nztcs_sp$species_clean))

# -------------------------------
# Join matched records
# -------------------------------

merged <- records %>%
  left_join(nztcs_sp %>% select(species_clean, Status, Category, BioStatus, ThreatReport),
            by = c("matched_species" = "species_clean"))

# Safety net
if (!"Status" %in% names(merged)) merged$Status <- NA_character_

# -------------------------------
# Add Nga Awa catchment data
# -------------------------------

NA_polygons <- st_read("Data/Nga Awa shapefiles/DOC_NgāAwa_RiverSites_20250122_n14.shp")

all_samples_valid <- all_samples %>%
  filter(!is.na(Longitude) & !is.na(Latitude))

NA_sf <- st_as_sf(all_samples_valid, coords = c("Longitude", "Latitude"), crs = 4167)
NA_sf <- st_transform(NA_sf, st_crs(NA_polygons))

NA_joined <- st_join(NA_sf, NA_polygons, left = TRUE)

all_samples_valid$Nga_Awa_Catchment <- ifelse(is.na(NA_joined$Waterway_N), "not Nga Awa", NA_joined$Waterway_N)

all_samples <- all_samples %>%
  left_join(all_samples_valid %>% select(SID, Nga_Awa_Catchment), by = "SID") %>%
  mutate(Nga_Awa_Catchment = ifelse(is.na(Nga_Awa_Catchment), "not Nga Awa", Nga_Awa_Catchment))

# -------------------------------
# Add Regional Council
# -------------------------------

RC_polygons <- st_read("Data/Regional Council shapefiles/regional-council-2022-generalised.shp")

all_samples_valid <- all_samples %>%
  filter(!is.na(Longitude) & !is.na(Latitude))

RC_sf <- st_as_sf(all_samples_valid, coords = c("Longitude", "Latitude"), crs = 4167)
RC_sf <- st_transform(RC_sf, st_crs(RC_polygons))

RC_joined <- st_join(RC_sf, RC_polygons, left = TRUE)

all_samples_valid$Regional_Council <- ifelse(is.na(RC_joined$REGC2022_1), "None", RC_joined$REGC2022_1)

all_samples <- all_samples %>%
  left_join(all_samples_valid %>% select(SID, Regional_Council), by = "SID") %>%
  mutate(Regional_Council = ifelse(is.na(Regional_Council), "None", Regional_Council))
                               
# -------------------------------
# Save preprocessed data
# -------------------------------

all_merged <- merged %>%
  left_join(all_samples, by = 'UID') %>%
  mutate(
    Latitude = coalesce(as.numeric(Latitude.x), as.numeric(Latitude.y)),
    Longitude = coalesce(as.numeric(Longitude.x), as.numeric(Longitude.y)),
    ClientSampleID = coalesce(ClientSampleID.x, ClientSampleID.y)
  ) %>%
  select(-Latitude.x, -Latitude.y, -Longitude.x, -Longitude.y, -ClientSampleID.x, -ClientSampleID.y)

colnames(all_merged)[c(20,21)] <- c("Wilderlab_Sp_name", "NZTC_Sp_name")

selected_df <- all_merged %>%
  select(UID, TaxID, Rank, Name, CommonName, Group, Count, Latitude, Longitude, ClientSampleID, CollectionDate, 
         Status, Category, ThreatReport, CollectedBy, DOC_Data, MakeDataPublic, Nga_Awa_Catchment, Regional_Council, phylum,
         class, order, family, genus, species, Wilderlab_Sp_name, NZTC_Sp_name, Report)

#Merge selected_df by Wilderlab Report.

# Step 1: Get number of unique UID per Report
uid_counts <- selected_df %>%
  group_by(Report) %>%
  summarise(total_UID = n_distinct(UID), .groups = "drop")

# Step 2: Summarise by Report and TaxID
summary_df <- selected_df %>%
  group_by(Report, TaxID) %>%
  summarise(
    unique_UID_count = n_distinct(UID),
    UID_list = paste(unique(UID), collapse = "-"),
    sum_count = sum(Count, na.rm = TRUE),
    ClientSampleID = first(ClientSampleID),  # representative
    across(c(Rank, Name, CommonName, Group,Latitude, Longitude, CollectionDate, 
             Status, Category, ThreatReport, CollectedBy, DOC_Data, MakeDataPublic, Nga_Awa_Catchment, Regional_Council, phylum,
             class, order, family, genus, species, Wilderlab_Sp_name, NZTC_Sp_name), first), # keep other columns
    .groups = "drop"
  ) %>%
  left_join(uid_counts, by = "Report") %>%
  mutate(mean_count = sum_count / total_UID)

summary_df <- summary_df %>% 
  select(Name, CommonName, ClientSampleID, sum_count, mean_count, unique_UID_count, total_UID, Rank, Group, Status,
         Category, Latitude, Longitude, CollectionDate,ThreatReport, CollectedBy,DOC_Data, MakeDataPublic, 
         Nga_Awa_Catchment, Regional_Council, phylum, class, order, family, genus, species, Wilderlab_Sp_name, NZTC_Sp_name,
                                    TaxID, Report)

summary_df <- summary_df %>%
  rename(
    `Latin Name` = Name,
    `Common Name` = CommonName,
    `Sample Name` = ClientSampleID,
    `Total Reads` = sum_count,
    `Average reads` = mean_count,
    `Hit number` = unique_UID_count,
    `Total hits` = total_UID,
    `Taxonomic Rank` = Rank,
    `Taxon Group` = Group,
    `Threat Status` = Status,
    `Threat Category` = Category,
    `Date` = CollectionDate,
    `Threat Document` = ThreatReport,
    `Collector` = CollectedBy,
    `DOC Data` = DOC_Data,
    `Public/Private` = MakeDataPublic, 
    `Nga Awa` = Nga_Awa_Catchment,
    `Regional Council` = Regional_Council,
    `Wilderlab Sp Name` = Wilderlab_Sp_name,
    `NZTC Sp Name` = NZTC_Sp_name,
    `TaxID` = TaxID,
    `Wilderlab Report` = Report
  )

# Save data for export into shinyapp

saveRDS(summary_df, "Data/records.rds")
