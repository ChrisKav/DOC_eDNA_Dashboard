#!/usr/bin/env Rscript
#
# build_taxon_tree.R
#
# For a user-defined taxonomic group (typically Family), this script:
#   1. Loads Data/records.rds (produced by get_API_data.R)
#   2. Finds all Wilderlab jobs (Reports) that detected taxa in that group
#   3. Looks up each job's ResultsLink in the Wilderlab `jobs` table
#      (re-fetched via the wilderlab API; cached to Data/jobs.rds)
#   4. Downloads the per-job .xlsx result files into Data/xlsx_cache/
#   5. Extracts sequences from the `full` sheet, filtered to the target group
#      plus its sisters (everything in the same higher rank for tree context)
#   6. Splits sequences by genetic marker (the `Target` column: CI, WV, TP, ...)
#   7. Aligns each marker independently with DECIPHER
#   8. Builds neighbour-joining trees per marker (ape) and an optional
#      concatenated supermatrix tree across markers
#   9. Annotates tips with Species + Genus labels, and flags tips whose
#      reference rank is *coarser* than species (genus / family / no rank)
#      as "candidate undescribed taxa"
#  10. Writes per-marker FASTA, alignment, Newick tree and an annotated PDF
#      figure to Output/Phylogenetics/<group_name>/
#
# Usage (interactive):
#   source("build_taxon_tree.R")
#   build_taxon_tree(rank = "Family", taxon = "Galaxiidae")
#
# Usage (command line):
#   Rscript build_taxon_tree.R --rank Family --taxon Galaxiidae
#   Rscript build_taxon_tree.R --rank Genus  --taxon Galaxias --markers CI,WV
#
# Conventions follow get_API_data.R: data.table for tabular work, msg() for
# logging, fault-tolerant API calls, and Data/wilder_keys.csv (or WILDER_*
# env vars) for credentials.
#

options(stringsAsFactors = FALSE)
# Stop ggtree from erroring out on small NJ trees whose midpoint-rooted
# branches go slightly negative (a normal artefact of NJ + midpoint).
options(ignore.negative.edge = TRUE)
# R's default 60s timeout is too short for the >300 MB public records.csv;
# bump it so the UID -> JobID bridge has a chance to complete.
options(timeout = max(600, getOption("timeout", 60)))
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(stringr)
  library(dplyr)
  library(httr)
  library(Biostrings)
  library(DECIPHER)
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(wilderlab)
})

# ============================================================================
# 0. Helpers (mirrored from get_API_data.R for consistency)
# ============================================================================

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a) && nzchar(as.character(a))) a else b
msg    <- function(...) cat(sprintf(...), "\n")

# Read Wilderlab API credentials from CSV or env, identical pattern to
# get_API_data.R::read_keys() so a single keys file works for both scripts.
read_keys <- function(keys_file = Sys.getenv("WILDER_KEYS_FILE",
                                             unset = "Data/wilder_keys.csv")) {
  if (file.exists(keys_file)) {
    kv <- tryCatch({
      dt <- fread(keys_file, header = FALSE, sep = ",", showProgress = FALSE)
      if (ncol(dt) < 2) stop("Keys CSV must have at least two columns")
      setNames(trimws(as.character(dt[[2]])),
               tolower(trimws(as.character(dt[[1]]))))
    }, error = function(e) { warning("Failed to read keys: ", e$message); NULL })
    if (!is.null(kv)) {
      return(list(
        key     = kv[["key"]]     %||% kv[["api key"]]    %||% NA_character_,
        secret  = kv[["secret"]]  %||% kv[["api secret"]] %||% NA_character_,
        xapikey = kv[["xapikey"]] %||% kv[["x-api-key"]]  %||%
                  kv[["x_api_key"]] %||% NA_character_
      ))
    }
  }
  list(
    key     = Sys.getenv("WILDER_KEY",     unset = NA_character_),
    secret  = Sys.getenv("WILDER_SECRET",  unset = NA_character_),
    xapikey = Sys.getenv("WILDER_XAPIKEY",
                         unset = Sys.getenv("WILDER_X_API_KEY",
                                            unset = NA_character_))
  )
}

# Retrying fetch wrapper around get_wilderdata(), same pattern as
# get_API_data.R::fetch_wilder_safe(). Returns NULL on persistent failure.
fetch_wilder_safe <- function(type, ..., max_attempts = 4, base_sleep = 1) {
  attempt <- 1L
  repeat {
    res <- try(get_wilderdata(type, ...), silent = TRUE)
    if (!inherits(res, "try-error")) return(res)
    msg("Attempt %d for %s failed: %s", attempt, type,
        conditionMessage(attr(res, "condition")))
    if (attempt >= max_attempts) {
      warning("Giving up on ", type, " after ", attempt, " attempts")
      return(NULL)
    }
    Sys.sleep(base_sleep * attempt)
    attempt <- attempt + 1L
  }
}

# ============================================================================
# 1. Load records.rds and filter to the target taxonomic group
# ============================================================================

#' Load Data/records.rds and return the rows matching the chosen taxon.
#'
#' @param rank  Taxonomic level column to match against. One of
#'              "Phylum", "Class", "Order", "Family", "Genus", "Species".
#' @param taxon Taxon name (case-insensitive); e.g. "Galaxiidae".
#' @param records_path Path to the records RDS file.
#' @return data.table of records with one row per UID.
filter_records_to_taxon <- function(rank, taxon,
                                    records_path = "Data/records.rds") {
  if (!file.exists(records_path)) {
    stop("records.rds not found at ", records_path,
         "\n  Run get_API_data.R first to build it.")
  }
  records <- as.data.table(readRDS(records_path))

  if (!rank %in% names(records)) {
    stop("Rank column '", rank, "' not found in records.rds.\n",
         "  Available taxonomic columns: ",
         paste(intersect(c("Phylum","Class","Order","Family","Genus","Species"),
                         names(records)), collapse = ", "))
  }

  hits <- records[tolower(get(rank)) == tolower(taxon)]
  msg("Records matching %s = '%s': %d aggregated rows",
      rank, taxon, nrow(hits))
  if (nrow(hits) == 0) {
    stop("No records matched. Try a different rank or check spelling.")
  }

  # get_API_data.R aggregates DT by (Report, TaxID) and stores UIDs as a
  # "-"-joined string in `UID_list`. Variant pipelines may rename this or
  # keep a per-row `UID` column instead. Try all the likely names.
  uid_src_col <- NULL
  for (candidate in c("UID_list", "UID list", "UID_List", "uid_list", "UID")) {
    if (candidate %in% names(hits)) { uid_src_col <- candidate; break }
  }
  if (is.null(uid_src_col)) {
    stop("records.rds has no UID-like column. Looked for: ",
         "UID_list, 'UID list', UID_List, uid_list, UID.\n",
         "  Available columns: ",
         paste(head(names(hits), 30), collapse = ", "),
         if (length(names(hits)) > 30) ", ..." else "")
  }
  msg("Using '%s' as UID source.", uid_src_col)

  # If there's a stale "UID" column hanging around, drop it before exploding
  # so the new UID column we produce is unambiguous.
  if ("UID" %in% names(hits) && uid_src_col != "UID") hits[, UID := NULL]

  # Explode the "-"-joined UID list into one row per UID. We use a per-row
  # id and merge approach — avoids the data.table scoping quirk that broke
  # the previous `by = setdiff(names(hits), 'UID')` pattern when UID wasn't
  # actually a column in scope at the time j was evaluated.
  hits[, .uid_raw := as.character(get(uid_src_col))]
  hits <- hits[!is.na(.uid_raw) & nzchar(.uid_raw)]
  hits[, .row_id := .I]

  exploded <- hits[, .(UID = unlist(strsplit(.uid_raw, "-", fixed = TRUE))),
                   by = .row_id]

  hits[, .uid_raw := NULL]
  hits <- merge(hits, exploded, by = ".row_id", allow.cartesian = TRUE)
  hits[, .row_id := NULL]

  hits[, UID := as.character(UID)]
  hits <- hits[!is.na(UID) & nzchar(UID)]

  msg("After expanding UIDs: %d rows, %d unique UIDs",
      nrow(hits), length(unique(hits$UID)))
  hits
}

# ============================================================================
# 2. Build a UID -> JobID map from the samples tables
# ============================================================================

#' Return a data.table with columns (UID, JobID, Report_URL, Source) merging
#' samples from BOTH the authenticated Wilderlab API AND the public S3
#' dataset — the same combine-both pattern get_API_data.R uses for
#' `all_samples`. Rows are tagged with Source = "DOC_API" or "Public_S3"
#' so you can tell where each (UID, JobID) pairing came from.
#'
#' NOTE: the public S3 samples.csv does NOT carry a JobID column — only
#' `Report` (the openwaters HTML viewer URL). We capture that as
#' Report_URL so downstream code can still group/track those UIDs even
#' though they can't be resolved to a ResultsLink without manual help.
get_samples_table <- function() {
  parts <- list()

  # 1. Authenticated API samples (your account; has JobID directly)
  keys <- read_keys()
  if (!any(is.na(c(keys$key, keys$secret, keys$xapikey)))) {
    msg("Fetching samples from Wilderlab API...")
    s_api <- fetch_wilder_safe("samples",
                               key = keys$key, secret = keys$secret,
                               xapikey = keys$xapikey)
    if (!is.null(s_api) && nrow(s_api) > 0 && "UID" %in% names(s_api)) {
      sa <- as.data.table(s_api)
      parts$api <- data.table(
        UID        = as.character(sa$UID),
        JobID      = if ("JobID"  %in% names(sa)) as.character(sa$JobID)  else NA_character_,
        Report_URL = if ("Report" %in% names(sa)) as.character(sa$Report) else NA_character_,
        Source     = "DOC_API"
      )
      msg("  API samples: %d rows", nrow(parts$api))
    }
  } else {
    msg("API credentials not set — skipping authenticated samples fetch.")
  }

  # 2. Public S3 samples — these have UID and Report (openwaters HTML URL)
  # but typically NO JobID column. We still pull what's there.
  pub_url <- "http://s3.ap-southeast-2.amazonaws.com/wilderlab.publicdata/samples.csv"
  msg("Fetching public samples from %s", pub_url)
  s_pub <- tryCatch(fread(pub_url, showProgress = FALSE),
                    error = function(e) {
                      warning("Public samples.csv fetch failed: ", e$message)
                      NULL
                    })
  if (!is.null(s_pub) && nrow(s_pub) > 0 && "UID" %in% names(s_pub)) {
    sp <- as.data.table(s_pub)
    has_jobid  <- "JobID"  %in% names(sp)
    has_report <- "Report" %in% names(sp)
    parts$public <- data.table(
      UID        = as.character(sp$UID),
      JobID      = if (has_jobid)  as.character(sp$JobID)  else NA_character_,
      Report_URL = if (has_report) as.character(sp$Report) else NA_character_,
      Source     = "Public_S3"
    )
    msg("  Public samples: %d rows (JobID column present: %s, Report column present: %s)",
        nrow(parts$public), has_jobid, has_report)
  }

  # 3. Public S3 records.csv — UID -> JobID bridge for public data.
  # Wilderlab's live API records table has no JobID column (the API
  # requires JobID as an input parameter, so it's implicit). But when
  # those records are exported to the flat public CSV, the JobID column
  # is typically stamped onto each row. Use it as the bridge that the
  # public samples.csv can no longer provide.
  pub_records_url <- "http://s3.ap-southeast-2.amazonaws.com/wilderlab.publicdata/records.csv"
  msg("Fetching public records.csv for UID -> JobID bridge from %s",
      pub_records_url)
  r_pub <- tryCatch(fread(pub_records_url, showProgress = FALSE),
                    error = function(e) {
                      warning("Public records.csv fetch failed: ", e$message)
                      NULL
                    })
  if (!is.null(r_pub) && nrow(r_pub) > 0) {
    rp <- as.data.table(r_pub)
    job_col <- intersect(c("JobID", "Job_ID", "Job", "jobid", "job_id"),
                         names(rp))[1]
    if (!is.na(job_col) && "UID" %in% names(rp)) {
      parts$public_records <- unique(rp[, .(
        UID        = as.character(get("UID")),
        JobID      = as.character(get(job_col)),
        Report_URL = NA_character_,
        Source     = "Public_S3_records"
      )])
      msg("  Public records.csv: %d unique (UID, JobID) pairings via '%s'",
          nrow(parts$public_records), job_col)
    } else {
      msg("  Public records.csv has no JobID-like column. Cols: %s",
          paste(head(names(rp), 15), collapse = ", "))
    }
  }

  if (length(parts) == 0) {
    stop("Could not obtain a samples table from any source.")
  }

  # 4. Merge — rbindlist + dedup per UID, but prefer rows that actually
  # have a JobID (otherwise dedup could keep a JobID-less public_samples
  # row over the JobID-bearing public_records row for the same UID).
  combined <- rbindlist(parts, fill = TRUE, use.names = TRUE)
  combined[, has_job := !is.na(JobID) & nzchar(JobID)]
  # Sort: rows with JobID first, then API > public_records > public_samples
  source_rank <- c("DOC_API" = 1L, "Public_S3_records" = 2L, "Public_S3" = 3L)
  combined[, src_rank := source_rank[Source]]
  setorder(combined, UID, -has_job, src_rank)
  out <- combined[, .SD[1], by = UID]
  out[, c("has_job", "src_rank") := NULL]

  msg("UID lookup table: %d unique UIDs (with JobID: %d, with Report_URL: %d)",
      nrow(out),
      sum(!is.na(out$JobID) & nzchar(out$JobID)),
      sum(!is.na(out$Report_URL) & nzchar(out$Report_URL)))
  src_tbl <- out[, .N, by = Source]
  msg("  Provenance: %s",
      paste(sprintf("%s=%d", src_tbl$Source, src_tbl$N), collapse = ", "))
  out
}

# ============================================================================
# 3. Obtain the jobs table (with ResultsLink) — try multiple sources
# ============================================================================

#' Return a data.table with at minimum columns (JobID, ResultsLink).
#' Sources, tried in priority order:
#'   1. Data/manual_jobs.csv      — user override (best for stubborn cases).
#'                                  Expected columns: JobID, ResultsLink.
#'   2. Wilderlab API "jobs"      — only your account's jobs.
#'   3. Public S3 jobs.csv        — only present if Wilderlab publishes one.
#'   4. Data/jobs.rds             — cached result of a previous run.
#' The combined table is cached back to Data/jobs.rds.
get_jobs_table <- function(jobs_cache_path = "Data/jobs.rds",
                           manual_csv      = "Data/manual_jobs.csv",
                           force_refresh   = FALSE) {

  parts <- list()

  # 3a. manual_jobs.csv: simplest way for the user to supply ResultsLinks
  # that no API can find. Two-column CSV is enough.
  if (file.exists(manual_csv)) {
    mj <- tryCatch(fread(manual_csv, showProgress = FALSE),
                   error = function(e) {
                     warning("Could not read ", manual_csv, ": ", e$message)
                     NULL
                   })
    if (!is.null(mj) && all(c("JobID", "ResultsLink") %in% names(mj))) {
      mj[, JobID := as.character(JobID)]
      parts$manual <- mj[, .(JobID, ResultsLink)]
      msg("Loaded %d manual JobID -> ResultsLink rows from %s",
          nrow(mj), manual_csv)
    }
  }

  # 3b. authenticated API jobs (covers your own account's jobs)
  keys <- read_keys()
  if (!any(is.na(c(keys$key, keys$secret, keys$xapikey)))) {
    msg("Fetching jobs from Wilderlab API...")
    j_api <- fetch_wilder_safe("jobs",
                               key = keys$key, secret = keys$secret,
                               xapikey = keys$xapikey)
    if (!is.null(j_api) && nrow(j_api) > 0 &&
        all(c("JobID", "ResultsLink") %in% names(j_api))) {
      j_api_dt <- as.data.table(j_api)
      j_api_dt[, JobID := as.character(JobID)]
      parts$api <- j_api_dt[, .(JobID, ResultsLink)]
    }
  } else {
    msg("API credentials not set — skipping authenticated jobs fetch.")
  }

  # 3c. public S3 jobs.csv (may or may not exist; try anyway)
  pub_url <- "http://s3.ap-southeast-2.amazonaws.com/wilderlab.publicdata/jobs.csv"
  j_pub <- tryCatch(suppressWarnings(fread(pub_url, showProgress = FALSE)),
                    error = function(e) NULL)
  if (!is.null(j_pub) && nrow(j_pub) > 0 &&
      all(c("JobID", "ResultsLink") %in% names(j_pub))) {
    msg("Loaded public jobs.csv (%d rows)", nrow(j_pub))
    jp <- as.data.table(j_pub)
    jp[, JobID := as.character(JobID)]
    parts$public <- jp[, .(JobID, ResultsLink)]
  } else {
    msg("No public jobs.csv available (404 or missing columns) — skipping.")
  }

  # 3d. previous cache (only used if we got nothing fresh)
  if (!force_refresh && length(parts) == 0 && file.exists(jobs_cache_path)) {
    msg("Falling back to cached jobs at %s", jobs_cache_path)
    cached <- as.data.table(readRDS(jobs_cache_path))
    if (all(c("JobID", "ResultsLink") %in% names(cached))) {
      cached[, JobID := as.character(JobID)]
      parts$cached <- cached[, .(JobID, ResultsLink)]
    }
  }

  if (length(parts) == 0) {
    stop("Could not obtain ANY jobs/ResultsLink data.\n",
         "  Provide Wilderlab API credentials, or supply Data/manual_jobs.csv\n",
         "  with two columns: JobID, ResultsLink.")
  }

  # Combine: priority order is manual > api > public > cached. Use first
  # non-NA ResultsLink per JobID.
  priority <- c("manual" = 1L, "api" = 2L, "public" = 3L, "cached" = 4L)
  combined <- rbindlist(lapply(names(parts), function(nm) {
    d <- copy(parts[[nm]]); d[, src := nm]; d
  }), use.names = TRUE, fill = TRUE)
  combined[, prio := priority[src]]
  setorder(combined, JobID, prio)
  jobs <- combined[, .SD[1], by = JobID][, .(JobID, ResultsLink, source = src)]

  dir.create(dirname(jobs_cache_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(jobs, jobs_cache_path)
  msg("Combined jobs table: %d rows (manual=%d, api=%d, public=%d, cached=%d) -> %s",
      nrow(jobs),
      nrow(jobs[source == "manual"]),
      nrow(jobs[source == "api"]),
      nrow(jobs[source == "public"]),
      nrow(jobs[source == "cached"]),
      jobs_cache_path)
  jobs
}

# ============================================================================
# 3. Resolve a job's ResultsLink to the latest .xlsx URL and download it
# ============================================================================

#' Pick the most recent xlsx URL from a comma-separated ResultsLink string.
#' Wilderlab embeds an ISO-ish timestamp ("YYMMDDHHMM") into the filename,
#' so the last sortable string is the freshest report.
pick_latest_results_url <- function(results_link) {
  if (is.na(results_link) || !nzchar(results_link)) return(NA_character_)
  urls <- str_split(results_link, ",\\s*")[[1]]
  urls <- urls[grepl("\\.xlsx$", urls, ignore.case = TRUE)]
  if (length(urls) == 0) return(NA_character_)
  # Sort by filename — the trailing timestamp orders chronologically
  urls[order(basename(urls), decreasing = TRUE)][1]
}

#' Download a single xlsx file, caching to Data/xlsx_cache/. Returns the
#' local path or NA on failure.
download_xlsx <- function(url, cache_dir = "Data/xlsx_cache", quiet = FALSE) {
  if (is.na(url) || !nzchar(url)) return(NA_character_)
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  local_path <- file.path(cache_dir, basename(url))
  if (file.exists(local_path) && file.size(local_path) > 1024) {
    if (!quiet) msg("  cached: %s", basename(local_path))
    return(local_path)
  }
  res <- tryCatch({
    resp <- GET(url, write_disk(local_path, overwrite = TRUE), timeout(120))
    if (status_code(resp) != 200) {
      warning("HTTP ", status_code(resp), " for ", url)
      file.remove(local_path)
      return(NA_character_)
    }
    if (!quiet) msg("  downloaded: %s (%d kB)", basename(local_path),
                    round(file.size(local_path) / 1024))
    local_path
  }, error = function(e) {
    warning("Download failed for ", url, ": ", e$message)
    if (file.exists(local_path)) file.remove(local_path)
    NA_character_
  })
  res
}

#' For every JobID in `job_ids`, find its ResultsLink in the jobs table and
#' download the latest xlsx. Returns a data.table: JobID, url, local_path,
#' status, source.
download_xlsx_for_jobs <- function(job_ids, jobs,
                                   cache_dir = "Data/xlsx_cache") {
  digits_only <- function(x) gsub("[^0-9]", "", as.character(x))
  jobs <- copy(jobs)
  jobs[, JobID_key := digits_only(as.character(JobID))]

  reps <- unique(as.character(job_ids))
  reps <- reps[!is.na(reps) & nzchar(reps)]
  msg("Resolving ResultsLink for %d unique JobIDs...", length(reps))

  out <- data.table(JobID = reps, url = NA_character_,
                    local_path = NA_character_,
                    status = "missing", source = NA_character_)
  for (i in seq_along(reps)) {
    rid <- reps[i]
    rid_key <- digits_only(rid)
    j <- jobs[JobID_key == rid_key]
    if (nrow(j) == 0) {
      out[i, status := "no_job_record"]
      next
    }
    url <- pick_latest_results_url(j$ResultsLink[1])
    if (is.na(url)) {
      out[i, status := "no_results_link"]; out[i, source := j$source[1]]
      next
    }
    out[i, `:=`(url = url, source = j$source[1])]
    lp <- download_xlsx(url, cache_dir = cache_dir)
    if (is.na(lp)) {
      out[i, status := "download_failed"]
    } else {
      out[i, `:=`(local_path = lp, status = "ok")]
    }
  }
  unreachable <- out[status != "ok"]
  msg("Download summary: ok=%d, no_job=%d, no_link=%d, failed=%d",
      sum(out$status == "ok"),
      sum(out$status == "no_job_record"),
      sum(out$status == "no_results_link"),
      sum(out$status == "download_failed"))
  if (nrow(unreachable) > 0) {
    miss_csv <- file.path(cache_dir, "_unreachable_jobs.csv")
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    fwrite(unreachable, miss_csv)
    msg("  -> %d unreachable JobIDs logged to %s",
        nrow(unreachable), miss_csv)
    msg("  To recover them: put a CSV at Data/manual_jobs.csv with two")
    msg("  columns (JobID, ResultsLink). One row per missing JobID.")
  }
  out
}

# ============================================================================
# 4. Extract sequences from each xlsx and filter to the target taxon
# ============================================================================

#' Read the `full` sheet from one xlsx and return rows whose taxonomic
#' lineage intersects the target taxon. Because the `full` sheet only
#' carries ScientificName + Group (not full lineage), we look up the
#' lineage for each TaxID against the per-job records subset already
#' filtered to the taxon, plus name-based fallbacks.
#'
#' @param xlsx_path Path to one Wilderlab results xlsx.
#' @param records_subset Rows of records.rds that are already filtered to
#'        the target taxon, restricted to this Report. We use this to know
#'        which (Report, UID, TaxID) tuples are "in scope".
#' @param rank One of "Phylum","Class","Order","Family","Genus","Species".
#' @param taxon The user-specified taxon name.
#' @return data.table of sequences in scope.
extract_sequences_one_job <- function(xlsx_path, records_subset, rank, taxon) {
  # Read the `full` sheet — this is the one with one row per sequence variant
  # per marker (Target).
  full <- tryCatch(read_excel(xlsx_path, sheet = "full"),
                   error = function(e) {
                     warning("Could not read 'full' sheet of ",
                             basename(xlsx_path), ": ", e$message)
                     NULL
                   })
  if (is.null(full) || nrow(full) == 0) return(NULL)
  full <- as.data.table(full)

  # Read metadata to extract the JobID (handles cases where filename and
  # internal JobID disagree).
  meta <- tryCatch(suppressMessages(read_excel(xlsx_path, sheet = "metadata",
                                               col_names = FALSE)),
                   error = function(e) NULL)
  job_id <- if (!is.null(meta)) {
    m <- as.data.table(meta)
    val <- m[grepl("^JobID", as.character(m[[1]]), ignore.case = TRUE)][[2]][1]
    as.character(val)
  } else NA_character_

  # The UID columns are whatever extra integer-named columns appear after
  # the standard taxonomy columns.
  fixed_cols <- c("Sequence", "Target", "ScientificName", "Rank",
                  "TaxID", "CommonName", "Group")
  uid_cols   <- setdiff(names(full), fixed_cols)

  # In-scope TaxIDs = those present in this report's slice of records.rds
  # that matched the target taxon. records.rds stores TaxID as character;
  # the xlsx full sheet stores it as numeric (e.g. 109671.0). Normalise both
  # to integer-string so the match works.
  norm_taxid <- function(x) {
    n <- suppressWarnings(as.integer(round(as.numeric(x))))
    ifelse(is.na(n), NA_character_, as.character(n))
  }
  in_scope_taxids <- unique(norm_taxid(records_subset$TaxID))
  in_scope_taxids <- in_scope_taxids[!is.na(in_scope_taxids)]

  # Also keep rows whose ScientificName is or starts with the target taxon
  # (handles "Galaxiidae", "Galaxias sp.", etc. that may not carry a
  # resolved TaxID matching records.rds).
  full[, TaxID_chr := norm_taxid(TaxID)]
  full[, match_taxid := !is.na(TaxID_chr) & TaxID_chr %in% in_scope_taxids]
  full[, match_name  := grepl(paste0("(^|\\s)", taxon), ScientificName,
                              ignore.case = TRUE)]

  out <- full[match_taxid | match_name]
  if (nrow(out) == 0) return(NULL)

  # Pivot UID columns to a list-column of sample read-counts per row, so
  # one Sequence-per-Target row carries every sample it was detected in.
  detect <- out[, ..uid_cols]
  detect_long <- lapply(seq_len(nrow(detect)), function(i) {
    vals <- as.numeric(unlist(detect[i, ]))
    nm   <- uid_cols
    keep <- !is.na(vals) & vals > 0
    if (!any(keep)) return(data.table(UID = character(), reads = numeric()))
    data.table(UID = nm[keep], reads = vals[keep])
  })

  result <- data.table(
    Report         = if (!is.na(job_id)) job_id
                     else tools::file_path_sans_ext(basename(xlsx_path)),
    Target         = out$Target,
    ScientificName = out$ScientificName,
    Rank           = out$Rank,
    TaxID          = out$TaxID_chr,
    CommonName     = out$CommonName,
    Group          = out$Group,
    Sequence       = out$Sequence,
    detections     = detect_long
  )
  result
}

# ============================================================================
# 5. Merge sequences across jobs and build per-marker datasets
# ============================================================================

#' Combine sequences from all jobs and dedupe so each unique
#' (Target, Sequence) pair appears once, with detections aggregated.
merge_sequences <- function(sequence_tables) {
  sequence_tables <- sequence_tables[!sapply(sequence_tables, is.null)]
  if (length(sequence_tables) == 0) {
    stop("No sequences extracted from any job.")
  }
  combined <- rbindlist(sequence_tables, use.names = TRUE, fill = TRUE)
  msg("Total sequence rows across jobs: %d (unique sequences: %d)",
      nrow(combined),
      length(unique(paste(combined$Target, combined$Sequence))))

  # Dedupe by (Target, Sequence). Keep the most resolved taxonomy for the
  # representative row (species > genus > family > ... > no rank).
  rank_order <- c("species" = 1L, "subspecies" = 1L, "subgenus" = 2L,
                  "genus" = 3L, "subfamily" = 4L, "family" = 5L,
                  "tribe" = 6L, "subtribe" = 6L, "suborder" = 6L,
                  "order" = 7L, "subclass" = 8L, "class" = 9L,
                  "subphylum" = 10L, "phylum" = 11L, "domain" = 12L,
                  "kingdom" = 13L, "clade" = 14L, "no rank" = 99L)
  rp <- rank_order[as.character(combined$Rank)]
  rp[is.na(rp)] <- 99L
  combined[, rank_priority := rp]
  setorder(combined, Target, Sequence, rank_priority)

  # Roll up: concatenate the unique detections across duplicate rows.
  dedup <- combined[, .(
    ScientificName = ScientificName[1],
    Rank           = Rank[1],
    TaxID          = TaxID[1],
    CommonName     = CommonName[1],
    Group          = Group[1],
    Reports        = paste(unique(Report), collapse = ";"),
    detections     = list(unique(rbindlist(detections, fill = TRUE)))
  ), by = .(Target, Sequence)]

  msg("After dedup by (Target, Sequence): %d rows", nrow(dedup))
  dedup
}

# ============================================================================
# 6. Align each marker (Target) with DECIPHER
# ============================================================================

#' Build a DNAStringSet for one marker and align with DECIPHER::AlignSeqs().
align_one_marker <- function(marker_dt, marker) {
  if (nrow(marker_dt) < 3) {
    msg("Marker %s has only %d sequences; skipping alignment.",
        marker, nrow(marker_dt))
    return(NULL)
  }
  # Sanitize sequences: uppercase, ATCGN only — DECIPHER tolerates ambiguity
  # codes but it's worth catching obvious garbage early.
  seqs <- toupper(gsub("[^ACGTUNRYSWKMBDHV-]", "N", marker_dt$Sequence))
  dnass <- DNAStringSet(seqs)

  # Build informative names: TipID = paste(Target, idx) — short, unique,
  # safe for Newick. We carry the long label separately for plotting.
  names(dnass) <- sprintf("%s_%03d", marker, seq_len(nrow(marker_dt)))

  msg("Aligning %d sequences for marker %s...", length(dnass), marker)
  aln <- tryCatch(AlignSeqs(dnass, verbose = FALSE, processors = NULL),
                  error = function(e) {
                    warning("Alignment failed for ", marker, ": ", e$message)
                    NULL
                  })
  aln
}

# ============================================================================
# 7. Build a neighbour-joining tree and decorate tips
# ============================================================================

#' Build a NJ tree from an alignment, returning an `ape::phylo` object with
#' tip metadata as an attribute "tip_data".
build_tree_one_marker <- function(aln, marker_dt, marker) {
  if (is.null(aln) || length(aln) < 3) return(NULL)
  # Convert DECIPHER alignment to ape DNAbin for dist.dna
  aln_dnabin <- as.DNAbin(aln)
  # K80 distance is a sensible default for short barcode markers
  d <- tryCatch(dist.dna(aln_dnabin, model = "K80",
                         pairwise.deletion = TRUE),
                error = function(e) {
                  warning("dist.dna failed for ", marker, ": ", e$message)
                  NULL
                })
  if (is.null(d) || any(is.nan(d)) || all(is.na(d))) {
    # Fall back to raw distance which tolerates large gaps
    d <- dist.dna(aln_dnabin, model = "raw", pairwise.deletion = TRUE)
  }
  tree <- tryCatch(nj(d), error = function(e) {
    warning("NJ failed for ", marker, ": ", e$message); NULL
  })
  if (is.null(tree)) return(NULL)

  # Midpoint-root for readability if phangorn is available, else leave
  # unrooted.
  tree <- tryCatch(phangorn::midpoint(tree),
                   error = function(e) tree)

  # NJ + midpoint on small/strangely-distanced trees can produce negative
  # edge lengths that crash ggtree. Clamp to 0 (preserves topology, just
  # tidies up the visual). Also guard against NA or non-finite edges.
  if (!is.null(tree$edge.length)) {
    tree$edge.length[!is.finite(tree$edge.length)] <- 0
    tree$edge.length <- pmax(tree$edge.length, 0)
  }

  # Build tip annotations aligned to tree$tip.label
  tip_ids <- tree$tip.label
  idx     <- match(tip_ids, names(aln))
  # detection-list -> sample counts (safe for NULL / empty entries)
  n_samp <- vapply(marker_dt$detections[idx], function(x) {
    if (is.null(x) || length(x) == 0) 0L else as.integer(nrow(x))
  }, integer(1))
  # Per-tip UID list (sorted, comma-joined). Truncated for the on-tree
  # label below, but the full list is kept in `uid_list` for the CSV.
  uid_list <- vapply(marker_dt$detections[idx], function(x) {
    if (is.null(x) || length(x) == 0) ""
    else paste(sort(unique(as.character(x$UID))), collapse = ",")
  }, character(1))
  uid_short <- vapply(strsplit(uid_list, ",", fixed = TRUE), function(u) {
    if (length(u) == 0 || all(!nzchar(u))) return("")
    head_u <- head(u, 3)
    if (length(u) > 3) paste0(paste(head_u, collapse = ","),
                              ",+", length(u) - 3, " more")
    else paste(head_u, collapse = ",")
  }, character(1))
  reads_total <- vapply(marker_dt$detections[idx], function(x) {
    if (is.null(x) || length(x) == 0 || !"reads" %in% names(x)) 0L
    else as.integer(sum(as.numeric(x$reads), na.rm = TRUE))
  }, integer(1))
  meta <- data.table(
    label          = tip_ids,
    ScientificName = marker_dt$ScientificName[idx],
    Rank           = marker_dt$Rank[idx],
    Group          = marker_dt$Group[idx],
    TaxID          = marker_dt$TaxID[idx],
    Reports        = marker_dt$Reports[idx],
    n_samples      = n_samp,
    uid_list       = uid_list,
    total_reads    = reads_total
  )
  # Flag tips that are *not* identified to species — those are the eDNA
  # variants of greatest interest for "undescribed taxa" follow-up.
  meta[, is_candidate := !tolower(Rank) %in% c("species", "subspecies")]
  # On-tree label: species + rank + UIDs (truncated). The full UID list
  # is exported separately to `tip_lookup_<marker>.csv`.
  meta[, tip_label := paste0(
    ifelse(is.na(ScientificName), "Unknown", ScientificName),
    " [", Rank, "] ",
    "n=", n_samples,
    ifelse(nzchar(uid_short), paste0(" UIDs: ", uid_short), "")
  )]
  attr(tree, "tip_data") <- meta
  tree
}

# ============================================================================
# 8. Identify candidate undescribed taxa
# ============================================================================

#' Summarise tips that look like candidate undescribed taxa: rank coarser
#' than species AND detected in at least one sample. Returns a data.table
#' ready to save as CSV alongside the tree.
summarise_candidates <- function(tree_per_marker) {
  rows <- list()
  for (marker in names(tree_per_marker)) {
    tree <- tree_per_marker[[marker]]
    if (is.null(tree)) next
    md <- attr(tree, "tip_data")
    cand <- md[is_candidate == TRUE & n_samples > 0]
    if (nrow(cand) == 0) next
    cand[, Target := marker]
    rows[[marker]] <- cand
  }
  if (length(rows) == 0) return(data.table())
  rbindlist(rows, fill = TRUE)
}

# ============================================================================
# 9. Plot the trees with ggtree, colour by Rank and highlight candidates
# ============================================================================

plot_tree_one_marker <- function(tree, marker, group_name, out_dir) {
  if (is.null(tree)) return(invisible(NULL))
  md <- attr(tree, "tip_data")

  # Safe x-axis upper bound (handles all-zero / NA edge lengths).
  max_depth <- tryCatch(max(node.depth.edgelength(tree), na.rm = TRUE),
                        error = function(e) 1)
  if (!is.finite(max_depth) || max_depth <= 0) max_depth <- 1

  pdf_path <- file.path(out_dir,
                        sprintf("tree_%s_%s.pdf",
                                gsub("[^A-Za-z0-9]+", "_", group_name),
                                marker))

  # --- Try ggtree first (richer plot) -----------------------------------
  # We attach metadata via inner_join on label rather than `%<+%`, since
  # %<+% can trip a `group_info` recycle error on small trees whose tip
  # data has fewer rows than the total node count.
  plotted <- tryCatch({
    tree_data <- ggtree::fortify(tree)
    tree_data <- merge(tree_data, md, by = "label", all.x = TRUE,
                       sort = FALSE)
    p <- ggtree(tree, layout = "rectangular", linewidth = 0.3) +
      geom_tiplab(data = tree_data,
                  aes(label = tip_label, colour = is_candidate),
                  size = 2.2, hjust = -0.05, na.rm = TRUE) +
      geom_tippoint(data = tree_data,
                    aes(colour = is_candidate),
                    size = 1.8, na.rm = TRUE) +
      scale_colour_manual(values = c(`TRUE` = "#D7263D",
                                     `FALSE` = "#1B998B"),
                          name = "Candidate undescribed",
                          na.value = "grey50",
                          na.translate = FALSE) +
      ggtitle(sprintf("%s — marker %s (NJ tree, K80)",
                      group_name, marker)) +
      theme_tree2() +
      theme(legend.position = "right",
            plot.title = element_text(size = 11, face = "bold")) +
      xlim(0, max_depth * 1.7)

    ggsave(pdf_path, p,
           width = 11, height = max(4, nrow(md) * 0.20),
           limitsize = FALSE)
    TRUE
  }, error = function(e) {
    warning("ggtree plot failed for marker ", marker,
            " (falling back to ape::plot.phylo): ", e$message)
    FALSE
  })

  # --- ape fallback: simpler but unconditionally robust ----------------
  if (!plotted) {
    # Re-label tips with the informative label for the ape plot
    t2 <- tree
    new_labels <- md$tip_label[match(t2$tip.label, md$label)]
    t2$tip.label <- ifelse(is.na(new_labels), t2$tip.label, new_labels)
    tip_col <- ifelse(md$is_candidate[match(tree$tip.label, md$label)],
                      "#D7263D", "#1B998B")
    tip_col[is.na(tip_col)] <- "grey50"

    pdf(pdf_path,
        width = 11,
        height = max(4, length(tree$tip.label) * 0.20))
    ape::plot.phylo(t2, type = "phylogram", cex = 0.7,
                    label.offset = max_depth * 0.01,
                    tip.color = tip_col,
                    main = sprintf("%s — marker %s (NJ tree, K80) [fallback]",
                                   group_name, marker))
    ape::add.scale.bar(cex = 0.6)
    legend("bottomleft", legend = c("Candidate undescribed", "Identified to species"),
           text.col = c("#D7263D", "#1B998B"), bty = "n", cex = 0.7)
    dev.off()
  }

  msg("Wrote %s", pdf_path)
  pdf_path
}

# ============================================================================
# 10. Top-level orchestration
# ============================================================================

#' Build per-marker phylogenetic trees for a taxonomic group.
#'
#' @param rank   Taxonomic rank to filter on (default "Family").
#' @param taxon  Taxon name to match in that rank (e.g. "Galaxiidae").
#' @param markers Optional character vector of Target codes to restrict to
#'        (e.g. c("CI","WV")). If NULL (default) every marker present is used.
#' @param min_seqs Minimum sequences a marker must have to attempt
#'        alignment + tree (default 3).
#' @param force_refresh_jobs If TRUE, ignore Data/jobs.rds and re-fetch.
#' @param merge_by_species If TRUE (default), clades that end up sharing
#'        the same species name (paraphyletic species on the NJ tree)
#'        are merged into one. Pre-merge IDs are kept as
#'        `genetic_subclade`; `n_subclades > 1` in the summary CSV is
#'        your cryptic-species shortlist.
#' @return Invisibly returns a list with all intermediate objects so the
#'         caller can inspect or post-process.
#' @export
build_taxon_tree <- function(rank = "Family",
                             taxon,
                             markers = NULL,
                             min_seqs = 3,
                             force_refresh_jobs = FALSE,
                             out_root = "Output/Phylogenetics",
                             merge_by_species = TRUE) {

  stopifnot(is.character(rank), length(rank) == 1,
            is.character(taxon), length(taxon) == 1)
  rank <- tools::toTitleCase(tolower(rank))

  # --- 1. records.rds ----------------------------------------------------
  msg("\n==== build_taxon_tree: %s = %s ====", rank, taxon)
  records_hit <- filter_records_to_taxon(rank, taxon)

  # --- 2. UID -> JobID map ----------------------------------------------
  samples_tbl <- get_samples_table()
  records_hit[, UID := as.character(UID)]
  records_hit[samples_tbl, `:=`(JobID      = i.JobID,
                                Report_URL = i.Report_URL),
              on = "UID"]
  matched_uids <- records_hit[!is.na(JobID) & nzchar(JobID)]
  unresolved   <- unique(records_hit[is.na(JobID) | !nzchar(JobID),
                                     .(UID, Report_URL)])
  msg("UIDs resolved to a JobID: %d / %d",
      length(unique(matched_uids$UID)),
      length(unique(records_hit$UID)))
  if (nrow(unresolved) > 0) {
    msg("  Unresolved UIDs: %d (Wilderlab's public samples.csv no longer ",
        nrow(unresolved))
    msg("    carries JobID, and the public records.csv didn't bridge them).")
    miss_csv <- "Data/unresolved_uids.csv"
    dir.create("Data", showWarnings = FALSE, recursive = TRUE)
    fwrite(unresolved, miss_csv)
    msg("    Wrote unresolved UID list to %s", miss_csv)
  }
  if (nrow(matched_uids) == 0) {
    stop("Could not map ANY taxon-matching UID to a JobID.\n",
         "  Wilderlab dropped the JobID column from public samples.csv, and\n",
         "  public records.csv didn't carry it either. You can recover by\n",
         "  populating Data/manual_jobs.csv with rows of (JobID, ResultsLink)\n",
         "  obtained from https://publicdata.wilderlab.co/ (use Region/UID\n",
         "  filters to find the jobs covering your taxon).")
  }

  # --- 3. jobs table with ResultsLink ------------------------------------
  jobs <- get_jobs_table(force_refresh = force_refresh_jobs)

  # --- 4. download xlsx --------------------------------------------------
  job_ids <- unique(as.character(matched_uids$JobID))
  msg("Distinct JobIDs to fetch: %d", length(job_ids))
  dls <- download_xlsx_for_jobs(job_ids, jobs)
  ok  <- dls[status == "ok"]
  if (nrow(ok) == 0) {
    stop("No xlsx files could be downloaded for any matching job.\n",
         "  See Data/xlsx_cache/_unreachable_jobs.csv for the list, then\n",
         "  populate Data/manual_jobs.csv with rows of (JobID, ResultsLink).")
  }

  # --- 5. extract sequences per job --------------------------------------
  msg("\nExtracting sequences from %d xlsx files...", nrow(ok))
  seq_tables <- vector("list", nrow(ok))
  for (i in seq_len(nrow(ok))) {
    jid <- ok$JobID[i]
    rec_sub <- matched_uids[as.character(JobID) == jid]
    # Use single-bracket + list() so NULL returns DON'T silently delete
    # the slot from seq_tables. `seq_tables[[i]] <- NULL` would shrink the
    # list and crash the loop a few dozen iterations later.
    seq_tables[i] <- list(tryCatch(
      extract_sequences_one_job(ok$local_path[i], rec_sub, rank, taxon),
      error = function(e) {
        warning("Job ", jid, " sequence extraction failed: ", e$message)
        NULL
      }))
    msg("  job %s -> %d sequences", jid,
        if (is.null(seq_tables[[i]])) 0 else nrow(seq_tables[[i]]))
  }

  # --- 5. merge ----------------------------------------------------------
  merged <- merge_sequences(seq_tables)
  if (!is.null(markers)) {
    merged <- merged[Target %in% markers]
    msg("Restricted to markers %s: %d rows",
        paste(markers, collapse = ","), nrow(merged))
  }

  # --- 6. per-marker alignment + tree -----------------------------------
  out_dir <- file.path(out_root, gsub("[^A-Za-z0-9]+", "_", taxon))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  marker_results <- list()
  trees <- list()
  for (mk in unique(merged$Target)) {
    sub <- merged[Target == mk]
    if (nrow(sub) < min_seqs) {
      msg("Marker %s: only %d sequences (< %d) — skipping.",
          mk, nrow(sub), min_seqs)
      next
    }

    aln <- align_one_marker(sub, mk)
    if (is.null(aln)) next

    # Write FASTA and alignment for downstream use
    fa_path <- file.path(out_dir, sprintf("seqs_%s.fasta", mk))
    aln_path <- file.path(out_dir, sprintf("aln_%s.fasta", mk))
    writeXStringSet(DNAStringSet(setNames(sub$Sequence,
                                          sprintf("%s_%03d", mk,
                                                  seq_len(nrow(sub))))),
                    fa_path)
    writeXStringSet(aln, aln_path)

    tree <- tryCatch(build_tree_one_marker(aln, sub, mk),
                     error = function(e) {
                       warning("Tree build failed for marker ", mk, ": ",
                               e$message)
                       NULL
                     })
    if (is.null(tree)) next

    nwk_path <- file.path(out_dir, sprintf("tree_%s.nwk", mk))
    tryCatch(write.tree(tree, nwk_path),
             error = function(e) warning("Newick write failed for ", mk,
                                         ": ", e$message))

    # Per-marker tip lookup: every tree tip mapped to its full UID list,
    # ScientificName, Rank, JobIDs, total reads. This is the table to
    # join back to records.rds / samples to find collection date, site,
    # lat/lon, etc. for each tip on the tree.
    lookup_path <- file.path(out_dir, sprintf("tip_lookup_%s.csv", mk))
    tryCatch({
      md_full <- attr(tree, "tip_data")
      if (!is.null(md_full) && nrow(md_full) > 0) {
        fwrite(md_full[, .(tip_id = label, ScientificName, Rank, Group,
                           TaxID, is_candidate, n_samples, total_reads,
                           JobIDs = Reports, UIDs = uid_list)],
               lookup_path)
        msg("  Wrote %s (%d tips)", lookup_path, nrow(md_full))
      }
    }, error = function(e)
       warning("Tip lookup CSV failed for ", mk, ": ", e$message))

    # Plot failures (commonly from ggtree quirks on small trees) must not
    # abort the rest of the pipeline — the FASTA, alignment, Newick and
    # tip_lookup CSV are still useful even when the PDF can't be drawn.
    tryCatch(plot_tree_one_marker(tree, mk, taxon, out_dir),
             error = function(e) {
               warning("Tree plot failed for marker ", mk,
                       " (FASTA / alignment / Newick / tip_lookup still saved): ",
                       e$message)
             })

    trees[[mk]] <- tree
    marker_results[[mk]] <- list(sequences = sub, alignment = aln,
                                 tree = tree, fasta = fa_path,
                                 alignment_path = aln_path,
                                 newick = nwk_path)
  }

  # --- 7. candidate-undescribed summary ----------------------------------
  candidates <- summarise_candidates(trees)
  if (nrow(candidates) > 0) {
    cand_path <- file.path(out_dir, "candidate_undescribed_taxa.csv")
    # Include the full UID list and total_reads so each candidate can be
    # traced directly to the source samples in records.rds / samples.
    out_cols <- intersect(c("Target", "tip_id" = "label", "ScientificName",
                            "Rank", "Group", "TaxID", "n_samples",
                            "total_reads", "JobIDs", "UIDs"),
                          names(candidates))
    fwrite(candidates[, .(Target, tip_id = label, ScientificName, Rank,
                          Group, TaxID, n_samples, total_reads,
                          JobIDs = Reports, UIDs = uid_list)],
           cand_path)
    msg("\nWrote %d candidate undescribed-taxon tips to %s",
        nrow(candidates), cand_path)
  } else {
    msg("\nNo candidate undescribed taxa flagged (all tips identified to species).")
  }

  # --- 8. geographic mapping --------------------------------------------
  geo <- tryCatch(map_taxon_tree(
                    result = list(records = records_hit,
                                  marker_results = marker_results,
                                  output_dir = out_dir),
                    taxon = taxon),
                  error = function(e) {
                    warning("map_taxon_tree failed: ", e$message)
                    NULL
                  })

  # --- 9. clade-level geographic mapping --------------------------------
  clades <- tryCatch(map_clades_geographically(
                       result = list(records = records_hit,
                                     marker_results = marker_results,
                                     output_dir = out_dir),
                       taxon = taxon,
                       merge_by_species = merge_by_species),
                     error = function(e) {
                       warning("map_clades_geographically failed: ",
                               e$message)
                       NULL
                     })

  msg("\nDone. Outputs in %s", out_dir)
  invisible(list(records = records_hit,
                 jobs = jobs,
                 downloads = dls,
                 merged_sequences = merged,
                 marker_results = marker_results,
                 trees = trees,
                 candidates = candidates,
                 geo = geo,
                 clades = clades,
                 output_dir = out_dir))
}

# ============================================================================
# 10b. Geographic mapping of tree tips
# ============================================================================

#' Map every tip on every marker tree to the geographic sites where it was
#' detected. Produces, per marker:
#'   - geo_detections_<MK>.csv : one row per (tip, UID) detection with lat/lon
#'   - geo_summary_<MK>.csv    : per-tip lat/lon range, catchments, regions
#'   - map_<taxon>_<MK>.html   : interactive Leaflet map (if leaflet installed)
#'   - map_static_<taxon>_<MK>.pdf : static ggplot scatter
#' Plus a consolidated geo_summary_all_markers.csv across all markers.
#'
#' Designed to be called directly on a `build_taxon_tree()` result, or
#' automatically as the final step of that function.
#'
#' @param result A list with $records (filtered records.rds), $marker_results
#'        (each with $tree carrying tip_data), $output_dir.
#' @param taxon Display name to use in titles / filenames. Defaults to
#'        basename(output_dir).
#' @param markers Optional subset of marker codes to map. NULL = all.
#' @export
map_taxon_tree <- function(result, taxon = NULL, markers = NULL) {

  if (is.null(result$records) || nrow(result$records) == 0) {
    msg("map_taxon_tree: no records — skipping.")
    return(invisible(NULL))
  }
  if (length(result$marker_results) == 0) {
    msg("map_taxon_tree: no marker results — skipping.")
    return(invisible(NULL))
  }

  out_dir <- result$output_dir
  if (is.null(taxon)) taxon <- basename(out_dir)
  taxon_clean <- gsub("[^A-Za-z0-9]+", "_", taxon)

  records <- as.data.table(result$records)
  if (!all(c("Latitude", "Longitude", "UID") %in% names(records))) {
    msg("map_taxon_tree: records lacks Latitude/Longitude/UID — skipping.")
    return(invisible(NULL))
  }

  # Build UID -> location lookup. Preserve other useful site metadata
  # (sample id, date, catchment, region) where present.
  loc_cols <- c("UID",
                intersect(c("Latitude","Longitude","ClientSampleID",
                            "CollectionDate","Nga_Awa_Catchment",
                            "Regional_Council"),
                          names(records)))
  uid_loc <- unique(records[, loc_cols, with = FALSE], by = "UID")
  uid_loc[, Latitude  := suppressWarnings(as.numeric(Latitude))]
  uid_loc[, Longitude := suppressWarnings(as.numeric(Longitude))]
  uid_loc <- uid_loc[!is.na(Latitude) & !is.na(Longitude)]
  msg("Mapping: %d UIDs have valid coordinates", nrow(uid_loc))
  if (nrow(uid_loc) == 0) return(invisible(NULL))

  marker_codes <- names(result$marker_results)
  if (!is.null(markers)) marker_codes <- intersect(marker_codes, markers)

  geo_summary_list <- list()
  map_paths <- list()

  for (mk in marker_codes) {
    tree <- result$marker_results[[mk]]$tree
    if (is.null(tree)) next
    md <- attr(tree, "tip_data")
    if (is.null(md) || nrow(md) == 0 || !"uid_list" %in% names(md)) next

    # Long form: one row per (tip, UID detection)
    tip_det <- md[, .(UID = unlist(strsplit(as.character(uid_list),
                                            ",", fixed = TRUE))),
                  by = .(tip_id = label, ScientificName, Rank,
                         is_candidate, tip_label)]
    tip_det <- tip_det[nzchar(UID)]
    tip_det <- merge(tip_det, uid_loc, by = "UID", all.x = TRUE)
    tip_det <- tip_det[!is.na(Latitude) & !is.na(Longitude)]
    if (nrow(tip_det) == 0) {
      msg("  Marker %s: 0 geocoded detections — skipping.", mk)
      next
    }
    msg("  Marker %s: %d geocoded detections across %d tips",
        mk, nrow(tip_det), length(unique(tip_det$tip_id)))

    # Per-tip geographic summary
    tip_geo <- tip_det[, .(
      n_detections = .N,
      unique_sites = if ("ClientSampleID" %in% names(.SD))
                       uniqueN(ClientSampleID) else NA_integer_,
      lat_min = min(Latitude, na.rm = TRUE),
      lat_max = max(Latitude, na.rm = TRUE),
      lon_min = min(Longitude, na.rm = TRUE),
      lon_max = max(Longitude, na.rm = TRUE),
      catchments = if ("Nga_Awa_Catchment" %in% names(.SD))
                     paste(unique(na.omit(Nga_Awa_Catchment)),
                           collapse = ";") else "",
      regions    = if ("Regional_Council" %in% names(.SD))
                     paste(unique(na.omit(Regional_Council)),
                           collapse = ";") else ""
    ), by = .(tip_id, ScientificName, Rank, is_candidate)]
    tip_geo[, Target := mk]
    setcolorder(tip_geo, c("Target", setdiff(names(tip_geo), "Target")))
    geo_summary_list[[mk]] <- tip_geo

    # CSVs
    sum_csv <- file.path(out_dir, sprintf("geo_summary_%s.csv", mk))
    det_csv <- file.path(out_dir, sprintf("geo_detections_%s.csv", mk))
    fwrite(tip_geo, sum_csv)
    fwrite(tip_det, det_csv)
    msg("    Wrote %s and %s", basename(sum_csv), basename(det_csv))

    # Interactive Leaflet map (clusters at high density, popups per point)
    if (requireNamespace("leaflet", quietly = TRUE) &&
        requireNamespace("htmlwidgets", quietly = TRUE)) {
      tryCatch(
        map_paths[[paste0(mk, "_html")]] <-
          .build_leaflet_map(tip_det, mk, taxon, taxon_clean, out_dir),
        error = function(e)
          warning("Leaflet map failed for ", mk, ": ", e$message))
    } else {
      msg("    (leaflet/htmlwidgets not installed — skipping HTML map)")
    }

    # Static PDF map (always)
    tryCatch(
      map_paths[[paste0(mk, "_pdf")]] <-
        .build_static_map(tip_det, mk, taxon, taxon_clean, out_dir),
      error = function(e)
        warning("Static map failed for ", mk, ": ", e$message))
  }

  if (length(geo_summary_list) > 0) {
    all_path <- file.path(out_dir, "geo_summary_all_markers.csv")
    fwrite(rbindlist(geo_summary_list, fill = TRUE), all_path)
    msg("Wrote consolidated %s", basename(all_path))
  }

  invisible(list(summaries = geo_summary_list, maps = map_paths))
}

# ---- internal: Leaflet HTML map ---------------------------------------------
.build_leaflet_map <- function(tip_det, mk, taxon, taxon_clean, out_dir) {
  td <- copy(tip_det)
  td[, color := fifelse(is_candidate == TRUE, "#D7263D", "#1B998B")]
  td[, popup := paste0(
    "<b>", tip_id, "</b><br>",
    "Species: ", ScientificName,
    " <i>(", Rank, ")</i><br>",
    "Candidate undescribed: <b>", is_candidate, "</b><br>",
    "UID: ", UID, "<br>",
    if ("ClientSampleID" %in% names(td)) paste0("Sample: ", ClientSampleID, "<br>") else "",
    if ("CollectionDate" %in% names(td)) paste0("Date: ", CollectionDate, "<br>") else "",
    if ("Nga_Awa_Catchment" %in% names(td)) paste0("Catchment: ", Nga_Awa_Catchment, "<br>") else "",
    if ("Regional_Council"  %in% names(td)) paste0("Region: ", Regional_Council) else ""
  )]

  m <- leaflet::leaflet(td) |>
    leaflet::addTiles() |>
    leaflet::addCircleMarkers(
      lng = ~Longitude, lat = ~Latitude,
      color = ~color, fillColor = ~color,
      radius = 5, stroke = TRUE, weight = 1, opacity = 0.85,
      fillOpacity = 0.65,
      popup = ~popup,
      label = ~paste0(tip_id, ": ", ScientificName),
      clusterOptions = leaflet::markerClusterOptions(
        spiderfyOnMaxZoom = TRUE,
        showCoverageOnHover = FALSE)) |>
    leaflet::addLegend("bottomright",
      colors = c("#D7263D", "#1B998B"),
      labels = c("Candidate undescribed", "Identified to species"),
      title = sprintf("%s — marker %s", taxon, mk),
      opacity = 0.9)

  html_path <- file.path(out_dir,
                         sprintf("map_%s_%s.html", taxon_clean, mk))
  htmlwidgets::saveWidget(m, html_path, selfcontained = TRUE)
  msg("    Wrote %s", basename(html_path))
  html_path
}

# ---- internal: static ggplot PDF map ----------------------------------------
.build_static_map <- function(tip_det, mk, taxon, taxon_clean, out_dir) {
  # Optional basemap (rnaturalearth) if installed — falls back to plain
  # scatter on white otherwise.
  basemap <- NULL
  if (requireNamespace("rnaturalearth", quietly = TRUE) &&
      requireNamespace("sf",            quietly = TRUE)) {
    basemap <- tryCatch(
      rnaturalearth::ne_countries(scale = "medium", returnclass = "sf"),
      error = function(e) NULL)
  }

  # Sensible map limits with a small buffer
  lon_lim <- range(tip_det$Longitude, na.rm = TRUE) + c(-0.5, 0.5)
  lat_lim <- range(tip_det$Latitude,  na.rm = TRUE) + c(-0.5, 0.5)

  p <- ggplot()
  if (!is.null(basemap)) {
    p <- p + ggplot2::geom_sf(data = basemap,
                              fill = "grey95", colour = "grey70",
                              linewidth = 0.2)
  }
  p <- p +
    ggplot2::geom_point(data = tip_det[is_candidate == FALSE],
                        ggplot2::aes(x = Longitude, y = Latitude),
                        colour = "#1B998B", size = 1.4, alpha = 0.55) +
    # Candidates drawn on top, bigger, red
    ggplot2::geom_point(data = tip_det[is_candidate == TRUE],
                        ggplot2::aes(x = Longitude, y = Latitude),
                        colour = "#D7263D", size = 2.2, alpha = 0.85) +
    ggplot2::coord_sf(xlim = lon_lim, ylim = lat_lim,
                      default_crs = sf::st_crs(4326),
                      expand = FALSE) +
    ggplot2::labs(
      title = sprintf("%s — geographic distribution (marker %s)",
                      taxon, mk),
      subtitle = sprintf("%d geocoded detections across %d tips (red = candidate undescribed)",
                         nrow(tip_det),
                         length(unique(tip_det$tip_id))),
      x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "grey92"),
                   plot.title = ggplot2::element_text(face = "bold"))

  # Fallback if rnaturalearth/sf missing: use coord_quickmap
  if (is.null(basemap)) {
    p <- p + ggplot2::coord_quickmap(xlim = lon_lim, ylim = lat_lim)
  }

  pdf_path <- file.path(out_dir,
                        sprintf("map_static_%s_%s.pdf", taxon_clean, mk))
  ggplot2::ggsave(pdf_path, p, width = 8, height = 10)
  msg("    Wrote %s", basename(pdf_path))
  pdf_path
}

# ============================================================================
# 10c. Clade-level geographic mapping
# ============================================================================

#' Define clades by cutting each marker tree at a distance threshold, then
#' map each clade geographically with one consistent colour per clade.
#'
#' Rationale: a clade of similar sequences detected in one region is the
#' signal of a real cryptic taxon; the same clade scattered nationwide is
#' more likely a widespread species with multiple haplotypes that the
#' reference database happens to lack. Looking at clades on a map
#' (rather than individual tips) makes this distinction visible.
#'
#' Method: hierarchical clustering (complete linkage by default) on the
#' tree's cophenetic distances, cut at `distance_threshold`. This is
#' equivalent to "all tips within the threshold of each other form a
#' clade". For NJ + K80 in eDNA barcode markers, ~0.03 is a reasonable
#' default; vary it to explore.
#'
#' @param result Output of build_taxon_tree().
#' @param taxon  Display name; defaults to basename(result$output_dir).
#' @param distance_threshold Cut height for clade definition (default 0.03).
#' @param markers Optional subset of marker codes.
#' @param linkage `hclust` linkage; one of "complete" (default), "single",
#'        "average".
#' @param min_clade_size Drop clades with fewer than this many tips from
#'        plots (still listed in CSVs).
#' @param min_clade_size Drop clades with fewer than this many tips from
#'        plots (still listed in CSVs).
#' @param max_clades_to_map Deprecated under the new default — every clade
#'        is drawn with a distinct colour from the cycling palette. Kept
#'        for backward compatibility; default `Inf` (no cap).
#' @param focus_clades Deprecated under the new default — every clade is
#'        drawn. Kept for backward compatibility; passing a value just
#'        prints an informational message.
#' Map every clade on every marker tree to the geographic sites where it
#' was detected, with clades defined by tree topology.
#'
#' Algorithm (the only one — no thresholds, no fixed K):
#'   1. Post-order traverse the tree; precompute the set of unique
#'      species names found in the subtree below every internal node.
#'   2. Top-down walk from root: if a node's subtree contains 0 or 1
#'      unique species, the whole subtree is one clade. Otherwise the
#'      node is a clade-split boundary and we recurse into its children.
#'   3. Each clade is named by its species (or "Unnamed clade N" when no
#'      species-level tip is inside).
#'   4. Optional: clades that end up sharing the same species name
#'      (paraphyletic species on the NJ tree) are merged via union-find,
#'      with the pre-merge IDs preserved on every tip as `genetic_subclade`.
#'      `n_subclades > 1` in the summary flags cryptic-species candidates.
#'
#' The number of clades is purely data-driven — there's no cap.
#'
#' @param result Output of build_taxon_tree().
#' @param taxon Display name; defaults to basename(result$output_dir).
#' @param markers Optional subset of marker codes.
#' @param min_clade_size Drop clades smaller than this from plots
#'        (still listed in CSVs). Default 1 (keep everything).
#' @param max_clades_to_map Deprecated; kept for API compatibility.
#'        Every clade is drawn distinctly.
#' @param focus_clades Deprecated; kept for API compatibility.
#' @param merge_by_species If TRUE (default), clades that end up sharing
#'        the same species name are merged back together. The pre-merge
#'        clade IDs are retained on every tip as `genetic_subclade`, so
#'        clades with `n_subclades > 1` in the summary are exactly the
#'        candidate cryptic species. Set FALSE to see the raw paraphyletic
#'        splits.
#' @export
map_clades_geographically <- function(result, taxon = NULL,
                                       markers = NULL,
                                       min_clade_size = 1L,
                                       max_clades_to_map = Inf,
                                       focus_clades = NULL,
                                       merge_by_species = TRUE) {

  if (length(result$marker_results) == 0) {
    msg("map_clades_geographically: no marker results — skipping."); return(invisible(NULL))
  }
  if (is.null(result$records) || nrow(result$records) == 0) {
    msg("map_clades_geographically: no records — skipping."); return(invisible(NULL))
  }

  out_dir <- result$output_dir
  if (is.null(taxon)) taxon <- basename(out_dir)
  taxon_clean <- gsub("[^A-Za-z0-9]+", "_", taxon)

  # UID -> location lookup (same approach as map_taxon_tree)
  records <- as.data.table(result$records)
  loc_cols <- c("UID",
                intersect(c("Latitude","Longitude","ClientSampleID",
                            "CollectionDate","Nga_Awa_Catchment",
                            "Regional_Council"), names(records)))
  uid_loc <- unique(records[, loc_cols, with = FALSE], by = "UID")
  uid_loc[, Latitude  := suppressWarnings(as.numeric(Latitude))]
  uid_loc[, Longitude := suppressWarnings(as.numeric(Longitude))]
  uid_loc <- uid_loc[!is.na(Latitude) & !is.na(Longitude)]

  marker_codes <- names(result$marker_results)
  if (!is.null(markers)) marker_codes <- intersect(marker_codes, markers)

  msg("\nClade-level mapping: topology recursion (one clade per subtree with ≤1 species)")
  out <- list()

  for (mk in marker_codes) {
    tree <- result$marker_results[[mk]]$tree
    if (is.null(tree)) next
    md <- attr(tree, "tip_data")
    if (is.null(md) || nrow(md) == 0 || !"uid_list" %in% names(md)) next
    if (length(tree$tip.label) < 2) {
      msg("  Marker %s: <2 tips, skipping.", mk); next
    }

    # 1. Define clades by tree topology only — no thresholds, no fixed K.
    #    Each clade is a connected subtree containing ≤1 unique species.
    #    The number of clades is purely data-driven.
    coph <- cophenetic.phylo(tree)
    hc   <- hclust(as.dist(coph), method = "complete")
    # Compute fine-grained sub-clusters at d=0.03 — surfaces cryptic
    # species candidates inside each named clade via `n_subclades`.
    raw_clade_id <- cutree(hc, h = 0.03)

    asg <- .assign_clades_by_topology(tree, md)
    clade_id  <- asg$clade_id[names(raw_clade_id)]
    tip_to_key <- asg$tip_to_name[names(raw_clade_id)]
    n_named   <- sum(!startsWith(asg$clade_name_per_id, "Unnamed clade"))
    n_unnamed <- sum( startsWith(asg$clade_name_per_id, "Unnamed clade"))
    msg("  Marker %s: topology -> %d clades (%d named, %d unnamed)",
        mk, length(unique(clade_id)), n_named, n_unnamed)

    # Optional: merge clades that ended up sharing a species name
    # (paraphyletic species on the NJ tree). Pre-merge IDs become
    # `genetic_subclade`.
    n_split_species  <- 0L
    n_raw_clades_pre <- length(unique(clade_id))
    if (isTRUE(merge_by_species)) {
      # Get the name of each clade
      clade_to_name <- unique(data.table(c = clade_id, n = tip_to_key))
      by_name <- split(clade_to_name$c, clade_to_name$n)
      by_name <- by_name[!startsWith(names(by_name), "Unnamed clade")]
      all_c  <- as.character(unique(clade_id))
      parent <- setNames(all_c, all_c)
      find_root <- function(x) {
        x <- as.character(x)
        while (parent[[x]] != x) {
          parent[[x]] <<- parent[[parent[[x]]]]
          x <- parent[[x]]
        }
        x
      }
      union_c <- function(a, b) {
        ra <- find_root(a); rb <- find_root(b)
        if (ra != rb) parent[[ra]] <<- rb
      }
      for (cs in by_name) {
        cs <- unique(cs)
        if (length(cs) > 1) {
          n_split_species <- n_split_species + 1L
          anchor <- cs[1]
          for (other in cs[-1]) union_c(anchor, other)
        }
      }
      if (n_split_species > 0) {
        merged <- vapply(clade_id,
                         function(c) as.integer(find_root(c)),
                         integer(1))
        names(merged) <- names(clade_id)
        clade_id <- merged
        msg("    merge_by_species: %d species spanned multiple subtrees; %d -> %d clades after merge",
            n_split_species, n_raw_clades_pre, length(unique(clade_id)))
      }
    }

    # 2. Relabel clades by total n_samples so clade 1 = biggest
    tmp_md <- copy(md)[, .(label, n_samples = as.integer(n_samples))]
    tmp_md[, clade_raw := clade_id[label]]
    clade_order <- tmp_md[, .(total_n = sum(n_samples)), by = clade_raw][
                          order(-total_n)]
    clade_order[, new_id := seq_len(.N)]
    relabel <- setNames(clade_order$new_id, clade_order$clade_raw)
    clade_assign <- data.table(
      label            = names(clade_id),
      clade            = as.integer(relabel[as.character(clade_id)]),
      genetic_subclade = as.integer(raw_clade_id[names(clade_id)])
    )
    # Attach the topology-derived clade name to every tip — this is the
    # canonical name from the recursion (species name or "Unnamed clade N")
    # and what gets used in legends/maps.
    clade_assign[, clade_key := tip_to_key[label]]

    sizes <- as.integer(table(clade_assign$clade))
    msg("  Marker %s: %d clades (sizes: %d–%d tips; %d singletons)",
        mk, length(unique(clade_assign$clade)),
        min(sizes), max(sizes), sum(sizes == 1))

    # 3. Build per-tip table with clade label + per-clade summary
    md_clade <- merge(md, clade_assign, by = "label", sort = FALSE)
    # The clade name IS the topology-derived label. Renumber "Unnamed
    # clade N" by the new clade IDs so the numbering stays consistent
    # after renumbering by total samples.
    clade_names <- unique(md_clade[, .(clade, clade_name = clade_key)])
    clade_names[, has_species_label := !startsWith(clade_name, "Unnamed clade")]
    # Renumber unnamed clades to match their final clade IDs
    clade_names[!has_species_label,
                clade_name := paste0("Unnamed clade ", clade)]
    clade_summary <- md_clade[, .(
      n_tips           = .N,
      n_samples_total  = sum(as.integer(n_samples)),
      reads_total      = sum(as.integer(total_reads)),
      taxa             = paste(sort(unique(ScientificName)), collapse = ";"),
      ranks            = paste(sort(unique(Rank)), collapse = ";"),
      n_candidate_tips = sum(is_candidate == TRUE),
      # n_subclades > 1 = cryptic-species candidate inside this named clade
      n_subclades      = uniqueN(genetic_subclade),
      tip_ids          = paste(label, collapse = ";")
    ), by = clade][order(clade)]
    clade_summary <- merge(clade_summary, clade_names, by = "clade",
                            sort = FALSE)
    clade_summary[, Target := mk]
    setcolorder(clade_summary, c("Target", "clade", "clade_name",
                                  "has_species_label", "n_subclades"))

    # 3. Long form: one row per (clade, tip, UID) with lat/lon
    det <- md_clade[, .(UID = unlist(strsplit(as.character(uid_list),
                                              ",", fixed = TRUE))),
                    by = .(clade, tip_id = label, ScientificName, Rank,
                           is_candidate)]
    det <- det[nzchar(UID)]
    det <- merge(det, uid_loc, by = "UID", all.x = TRUE)
    det <- det[!is.na(Latitude) & !is.na(Longitude)]

    # 4. Per-clade geographic stats
    clade_geo <- det[, .(
      n_geocoded_dets = .N,
      lat_min = min(Latitude), lat_max = max(Latitude),
      lon_min = min(Longitude), lon_max = max(Longitude),
      lat_range = max(Latitude) - min(Latitude),
      lon_range = max(Longitude) - min(Longitude),
      catchments = if ("Nga_Awa_Catchment" %in% names(.SD))
                     paste(sort(unique(na.omit(Nga_Awa_Catchment))),
                           collapse = ";") else "",
      regions    = if ("Regional_Council" %in% names(.SD))
                     paste(sort(unique(na.omit(Regional_Council))),
                           collapse = ";") else "",
      n_unique_sites = if ("ClientSampleID" %in% names(.SD))
                         uniqueN(ClientSampleID) else NA_integer_
    ), by = clade]
    clade_summary_full <- merge(clade_summary, clade_geo,
                                 by = "clade", all.x = TRUE)

    # 5. CSVs
    md_clade_named <- merge(md_clade, clade_names, by = "clade",
                             sort = FALSE)
    fwrite(md_clade_named[, .(Target = mk, clade, clade_name,
                              genetic_subclade,
                              tip_id = label,
                              ScientificName, Rank, is_candidate, n_samples,
                              total_reads, uid_list)],
           file.path(out_dir, sprintf("clades_%s.csv", mk)))
    fwrite(clade_summary_full,
           file.path(out_dir, sprintf("clades_summary_%s.csv", mk)))
    msg("    Wrote clades_%s.csv and clades_summary_%s.csv", mk, mk)

    # 6. Per-clade colour palette. Now every clade gets a distinct colour
    #    from the cycling palette — no "Other clades" bucket. If the user
    #    set focus_clades, those are highlighted by ordering but all
    #    clades remain coloured.
    all_clades <- sort(unique(clade_summary$clade))
    palette <- .clade_palette(length(all_clades))
    clade_colors <- setNames(palette, as.character(all_clades))

    # focus_clades retained as a parameter for backward compatibility but
    # now used only for messaging — every clade is still drawn.
    if (!is.null(focus_clades)) {
      msg("    focus_clades requested (%s) — every clade still drawn; focus is informational only.",
          paste(as.integer(focus_clades), collapse = ", "))
    }
    top_clades <- all_clades  # pass all to helpers

    # 7b. clade_id -> clade_name lookup, used to make legend / tip labels
    #    biologically readable ("C5: Aoteapsyche colonica" instead of "C5").
    name_lookup <- setNames(clade_names$clade_name,
                             as.character(clade_names$clade))

    # 8. Tree PDF with tips coloured by clade
    tryCatch({
      .plot_tree_by_clade(tree, md_clade, clade_colors, mk, taxon,
                          taxon_clean, out_dir, top_clades = top_clades,
                          name_lookup = name_lookup)
    }, error = function(e)
       warning("Clade tree plot failed for ", mk, ": ", e$message))

    # 9. Interactive Leaflet map with clades as toggleable layers
    if (requireNamespace("leaflet",    quietly = TRUE) &&
        requireNamespace("htmlwidgets", quietly = TRUE) && nrow(det) > 0) {
      tryCatch({
        .map_clades_leaflet(det, md_clade, clade_colors, mk, taxon,
                            taxon_clean, out_dir, min_clade_size,
                            top_clades = top_clades,
                            name_lookup = name_lookup)
      }, error = function(e)
         warning("Clade map (leaflet) failed for ", mk, ": ", e$message))
    }

    # 10. Static PDF map by clade
    if (nrow(det) > 0) {
      tryCatch({
        .map_clades_static(det, clade_colors, mk, taxon, taxon_clean,
                           out_dir, min_clade_size,
                           top_clades = top_clades,
                           name_lookup = name_lookup)
      }, error = function(e)
         warning("Clade map (static) failed for ", mk, ": ", e$message))
    }

    out[[mk]] <- list(clade_assign = md_clade_named,
                      clade_summary = clade_summary_full,
                      detections = det,
                      colors = clade_colors,
                      top_clades = top_clades,
                      name_lookup = name_lookup)
  }

  # Consolidated across all markers
  if (length(out) > 0) {
    all_summaries <- rbindlist(lapply(out, `[[`, "clade_summary"),
                                fill = TRUE)
    fwrite(all_summaries,
           file.path(out_dir, "clades_summary_all_markers.csv"))
  }
  invisible(out)
}

# ---- internal: distinguishable colour palette of any size -------------------
# For small N, picks from a curated qualitative palette; for larger N,
# falls back to evenly-spaced HSV hues which stay distinguishable even at
# 50+ clades (still imperfect for the eye past ~30, but readable in a
# legend and unique per clade).
.clade_palette <- function(n) {
  base <- c("#D7263D","#1B998B","#3D5A80","#F46036","#9C27B0","#FF9F1C",
            "#34A853","#E91E63","#4285F4","#FDDB3A","#5C415D","#EE6C4D",
            "#293241","#7B287D","#FBBC04","#0EAD69","#EA4335","#98C1D9",
            "#A33B20","#0B7A75","#603A40","#F18F01","#7768AE","#3B6064",
            "#C73E1D","#F46197","#37123C","#71A2B6","#235789","#F4D35E")
  if (n <= length(base)) return(base[seq_len(n)])
  # For more than length(base) clades, generate evenly-spaced HSV hues.
  # Vary saturation/value across blocks so adjacent clade numbers don't
  # look identical.
  sat <- rep(c(0.80, 0.55, 0.95), length.out = n)
  val <- rep(c(0.85, 0.95, 0.70), length.out = n)
  hue <- (seq_len(n) - 1L) / n
  grDevices::hsv(h = hue, s = sat, v = val)
}

# ---- internal: find the natural "barcode gap" threshold ---------------------
# Looks at the distribution of pairwise tree distances, fits a kernel
# density, and locates the deepest valley between the within-species peak
# (low distances) and the between-species peak (higher distances). That
# valley is the canonical "barcode gap". If no clear valley exists, falls
# back to a sensible default for COI-style data.
.find_barcode_gap <- function(d_matrix, default = 0.03,
                              search_quantile = 0.6, min_dist = 0.005) {
  d_vec <- d_matrix[upper.tri(d_matrix)]
  d_vec <- d_vec[is.finite(d_vec) & d_vec > 0]
  if (length(d_vec) < 20L) return(default)

  dens <- tryCatch(
    stats::density(d_vec, n = 1024L, from = 0,
                   to = stats::quantile(d_vec, 0.99, names = FALSE),
                   bw = "nrd0"),
    error = function(e) NULL)
  if (is.null(dens)) return(default)

  upper_x <- stats::quantile(d_vec, search_quantile, names = FALSE)
  mask <- dens$x >= min_dist & dens$x <= upper_x
  if (sum(mask) < 5L) return(default)
  x <- dens$x[mask]; y <- dens$y[mask]

  if (length(y) < 3L) return(default)
  # Local minima (strict): y[i] < both neighbours
  is_min <- c(FALSE,
              y[-c(1L, length(y))] < y[-c(length(y) - 1L, length(y))] &
              y[-c(1L, length(y))] < y[-c(1L, 2L)],
              FALSE)
  valley_idx <- which(is_min)
  if (length(valley_idx) == 0L) return(default)

  deepest <- valley_idx[which.min(y[valley_idx])]
  cut_x <- x[deepest]
  if (!is.finite(cut_x) || cut_x <= 0) return(default)
  cut_x
}

# ---- internal: assign clades by tree topology -------------------------------
# Algorithm (the only one — no thresholds, no fixed K):
#   1. Post-order traverse the rooted tree; precompute, for every node,
#      the set of unique species-level ScientificNames found in the
#      subtree rooted at that node.
#   2. Top-down walk from the root: a subtree with ≤1 unique species
#      becomes one clade. If a subtree has ≥2 species, this is a
#      clade-split boundary — recurse into the children.
#   3. Each resulting clade is named by its single species (if any) or
#      "Unnamed clade N" if no species-level tip is inside.
#
# Returns a list with:
#   $clade_id          int vector, named by tip label
#   $tip_to_name       char vector, named by tip label
#   $clade_name_per_id char vector, indexed by clade ID
.assign_clades_by_topology <- function(tree, md) {
  n_tips     <- length(tree$tip.label)
  n_internal <- tree$Nnode
  n_total    <- n_tips + n_internal

  if (n_tips < 2) {
    return(list(
      clade_id          = setNames(1L, tree$tip.label),
      tip_to_name       = setNames("Singleton", tree$tip.label),
      clade_name_per_id = setNames("Singleton", "1")
    ))
  }

  # Build tip -> species lookup (NA when the tip isn't species-level)
  tip_species <- setNames(rep(NA_character_, n_tips), tree$tip.label)
  if (!is.null(md) && nrow(md) > 0) {
    for (i in seq_len(nrow(md))) {
      lbl <- md$label[i]
      if (!(lbl %in% tree$tip.label)) next
      rk <- tolower(as.character(md$Rank[i]))
      nm <- as.character(md$ScientificName[i])
      if (rk %in% c("species", "subspecies") &&
          !is.na(nm) && nzchar(nm)) {
        tip_species[lbl] <- nm
      }
    }
  }

  # Build child relationships
  children <- vector("list", n_total)
  for (i in seq_len(nrow(tree$edge))) {
    p <- tree$edge[i, 1]
    c <- tree$edge[i, 2]
    children[[p]] <- c(children[[p]], c)
  }

  # Find root (the only node that's never a child)
  all_children <- unique(tree$edge[, 2])
  root_candidates <- setdiff(seq_len(n_total), all_children)
  root <- if (length(root_candidates) == 1L) root_candidates else (n_tips + 1L)

  # Post-order traversal (iterative; safe for deep trees)
  node_species   <- vector("list", n_total)
  node_tips_list <- vector("list", n_total)

  visited <- rep(FALSE, n_total)
  pending <- c(root)
  order   <- integer(0)
  while (length(pending) > 0) {
    top <- pending[length(pending)]
    if (top <= n_tips || visited[top]) {
      order   <- c(order, top)
      pending <- pending[-length(pending)]
    } else {
      visited[top] <- TRUE
      kids <- children[[top]]
      pending <- c(pending, kids)
    }
  }

  for (node in order) {
    if (node <= n_tips) {
      sp <- tip_species[tree$tip.label[node]]
      node_tips_list[[node]] <- node
      node_species[[node]]   <- if (is.na(sp)) character(0) else sp
    } else {
      kids <- children[[node]]
      node_tips_list[[node]] <- unlist(lapply(kids,
                                              function(k) node_tips_list[[k]]))
      node_species[[node]]   <- unique(unlist(lapply(kids,
                                              function(k) node_species[[k]])))
    }
  }

  # Top-down clade assignment using an explicit queue
  clade_id          <- integer(n_tips)
  clade_name_per_id <- character(0)
  next_id           <- 1L
  unnamed_counter   <- 1L

  to_process <- list(root)
  while (length(to_process) > 0) {
    node <- to_process[[1]]
    to_process <- to_process[-1]

    sp   <- node_species[[node]]
    tips <- node_tips_list[[node]]

    if (length(sp) <= 1) {
      cid <- next_id
      next_id <- next_id + 1L
      clade_id[tips] <- cid
      if (length(sp) == 1) {
        clade_name_per_id[cid] <- sp[1]
      } else {
        clade_name_per_id[cid] <- paste0("Unnamed clade ", unnamed_counter)
        unnamed_counter <- unnamed_counter + 1L
      }
    } else {
      if (node > n_tips) {
        for (k in children[[node]]) to_process[[length(to_process) + 1L]] <- k
      } else {
        # A tip with >1 species shouldn't happen, but cope gracefully
        cid <- next_id
        next_id <- next_id + 1L
        clade_id[node] <- cid
        clade_name_per_id[cid] <- sp[1]
      }
    }
  }

  names(clade_id) <- tree$tip.label
  tip_to_name <- setNames(clade_name_per_id[clade_id], tree$tip.label)

  list(
    clade_id          = clade_id,
    tip_to_name       = tip_to_name,
    clade_name_per_id = clade_name_per_id
  )
}

# ---- internal: clade name = dominant species-level name, NA if none ---------
# Used to label clades after tree-based cutting. Returns NA when the clade
# has no tip identified to species/subspecies — caller then names it
# "Unnamed clade N" and flags it as a candidate undescribed lineage.
.compute_clade_name <- function(dt_clade) {
  d <- copy(dt_clade)
  d[, rk := tolower(as.character(Rank))]
  sp <- d[rk %in% c("species", "subspecies") &
          !is.na(ScientificName) & nzchar(ScientificName)]
  if (nrow(sp) == 0L) return(NA_character_)
  sp[, n := suppressWarnings(as.integer(n_samples))]
  sp[is.na(n), n := 0L]
  agg <- sp[, .(total_n = sum(n)), keyby = ScientificName]
  setorder(agg, -total_n, ScientificName)
  agg$ScientificName[1L]
}

# ---- internal: tree PDF with tips coloured by clade -------------------------
# Every clade gets a colour from the cycling palette; no greying-out of
# "non-top" clades any more.
.plot_tree_by_clade <- function(tree, md_clade, clade_colors, mk, taxon,
                                taxon_clean, out_dir, top_clades = NULL,
                                name_lookup = NULL) {
  tip_clade <- md_clade$clade[match(tree$tip.label, md_clade$label)]
  tip_col   <- clade_colors[as.character(tip_clade)]
  tip_col[is.na(tip_col)] <- "grey60"

  tip_lab <- md_clade$tip_label[match(tree$tip.label, md_clade$label)]
  tip_lab <- ifelse(is.na(tip_lab), tree$tip.label,
                    sprintf("[C%d] %s", tip_clade, tip_lab))

  t2 <- tree
  t2$tip.label <- tip_lab

  pdf_path <- file.path(out_dir,
                        sprintf("tree_clades_%s_%s.pdf", taxon_clean, mk))
  pdf(pdf_path, width = 11,
      height = max(4, length(tree$tip.label) * 0.20))
  ape::plot.phylo(t2, type = "phylogram", cex = 0.7,
                  label.offset = max(node.depth.edgelength(tree),
                                     na.rm = TRUE) * 0.01,
                  tip.color = tip_col, no.margin = FALSE,
                  main = sprintf("%s — marker %s (%d clades)",
                                 taxon, mk,
                                 length(unique(md_clade$clade))))
  ape::add.scale.bar(cex = 0.6)

  # Inline legend — cap at ~30 entries for readability, with "..." marker
  # if more. (Full list is in clades_summary_<MK>.csv.)
  clade_n <- table(md_clade$clade)
  legend_clades <- as.integer(names(sort(clade_n, decreasing = TRUE)))
  show_n <- min(length(legend_clades), 30L)
  show <- legend_clades[seq_len(show_n)]
  .lab <- function(cid) {
    nm <- if (!is.null(name_lookup)) name_lookup[as.character(cid)] else NA
    nm <- if (is.null(nm) || is.na(nm) || !nzchar(nm)) "" else paste0(": ", nm)
    sprintf("C%d%s (%d tips)", cid, nm,
            as.integer(clade_n[as.character(cid)]))
  }
  legend_text <- vapply(show, .lab, character(1))
  if (length(legend_clades) > show_n) {
    legend_text <- c(legend_text,
                     sprintf("… +%d more (see CSV)",
                             length(legend_clades) - show_n))
  }
  legend("bottomleft",
         legend = legend_text,
         text.col = c(clade_colors[as.character(show)],
                      if (length(legend_clades) > show_n) "grey40"),
         bty = "n", cex = 0.6)
  dev.off()
  msg("    Wrote %s", basename(pdf_path))
  pdf_path
}

# ---- internal: Leaflet map with clades as toggleable groups -----------------
# Every clade gets its own toggleable layer and a distinct colour from the
# cycling palette. No "Other clades" bucket — when the data has many clades,
# colours simply recycle in a way that adjacent clades stay distinguishable.
.map_clades_leaflet <- function(det, md_clade, clade_colors, mk, taxon,
                                taxon_clean, out_dir, min_clade_size,
                                top_clades = NULL, name_lookup = NULL) {
  td <- copy(det)

  .clade_label <- function(cid) {
    nm <- if (!is.null(name_lookup)) name_lookup[as.character(cid)] else NA
    if (is.null(nm) || is.na(nm) || !nzchar(nm)) sprintf("C%d", cid)
    else sprintf("C%d: %s", cid, nm)
  }

  td[, color := clade_colors[as.character(clade)]]
  td[, group := vapply(clade, .clade_label, character(1))]

  # Filter small clades from the map (still in CSVs). Default min=1 keeps
  # everything.
  size_per_clade <- md_clade[, .N, by = clade]
  keep <- size_per_clade[N >= min_clade_size]$clade
  td <- td[clade %in% keep]
  if (nrow(td) == 0) return(invisible(NULL))

  td[, popup := paste0(
    "<b>", group, "</b><br>",
    tip_id, " — ", ScientificName, " (", Rank, ")<br>",
    "UID: ", UID, "<br>",
    if ("ClientSampleID" %in% names(td))
       paste0("Sample: ", ClientSampleID, "<br>") else "",
    if ("CollectionDate" %in% names(td))
       paste0("Date: ", CollectionDate, "<br>") else "",
    if ("Nga_Awa_Catchment" %in% names(td))
       paste0("Catchment: ", Nga_Awa_Catchment) else "")]

  # Order: clades by numeric ID so the panel reads top-down sensibly
  present <- sort(unique(td$clade))
  group_order <- vapply(present, .clade_label, character(1))

  m <- leaflet::leaflet() |> leaflet::addTiles()
  for (g in group_order) {
    sub <- td[group == g]
    if (nrow(sub) == 0) next
    m <- m |> leaflet::addCircleMarkers(
      data = sub,
      lng = ~Longitude, lat = ~Latitude,
      color = ~color, fillColor = ~color,
      radius = 5, weight = 1,
      opacity = 0.85, fillOpacity = 0.65,
      popup = ~popup,
      label = ~paste0(group, " — ", ScientificName),
      group = g)
  }

  m <- m |>
    leaflet::addLayersControl(overlayGroups = group_order,
      options = leaflet::layersControlOptions(collapsed = FALSE)) |>
    leaflet::addLegend("bottomright",
      colors = clade_colors[as.character(present)],
      labels = group_order,
      title = sprintf("%s — marker %s (%d clades)",
                      taxon, mk, length(present)),
      opacity = 0.9)

  html_path <- file.path(out_dir,
                         sprintf("map_clades_%s_%s.html",
                                 taxon_clean, mk))
  htmlwidgets::saveWidget(m, html_path, selfcontained = TRUE)
  msg("    Wrote %s", basename(html_path))
  html_path
}

# ---- internal: static ggplot map coloured by clade --------------------------
# All clades coloured distinctly via the cycling palette; no "Other" bucket.
.map_clades_static <- function(det, clade_colors, mk, taxon,
                                taxon_clean, out_dir, min_clade_size,
                                top_clades = NULL, name_lookup = NULL) {
  td <- copy(det)

  .clade_label <- function(cid) {
    nm <- if (!is.null(name_lookup)) name_lookup[as.character(cid)] else NA
    if (is.null(nm) || is.na(nm) || !nzchar(nm)) sprintf("C%d", cid)
    else sprintf("C%d: %s", cid, nm)
  }

  present <- sort(unique(td$clade))
  lvls    <- vapply(present, .clade_label, character(1))
  td[, clade_lbl := factor(vapply(clade, .clade_label, character(1)),
                            levels = lvls)]

  sizes <- table(td$clade)
  keep <- as.integer(names(sizes[sizes >= min_clade_size]))
  td <- td[clade %in% keep]
  if (nrow(td) == 0) return(invisible(NULL))

  pal <- clade_colors[as.character(present)]
  names(pal) <- lvls

  basemap <- NULL
  if (requireNamespace("rnaturalearth", quietly = TRUE) &&
      requireNamespace("sf",            quietly = TRUE)) {
    basemap <- tryCatch(rnaturalearth::ne_countries(scale = "medium",
                                                    returnclass = "sf"),
                        error = function(e) NULL)
  }
  lon_lim <- range(td$Longitude) + c(-0.5, 0.5)
  lat_lim <- range(td$Latitude)  + c(-0.5, 0.5)

  p <- ggplot()
  if (!is.null(basemap)) {
    p <- p + ggplot2::geom_sf(data = basemap, fill = "grey95",
                              colour = "grey70", linewidth = 0.2)
  }
  p <- p +
    ggplot2::geom_point(data = td,
                  ggplot2::aes(x = Longitude, y = Latitude,
                               colour = clade_lbl),
                  size = 1.7, alpha = 0.8) +
    ggplot2::scale_colour_manual(values = pal, name = "Clade",
                                  drop = FALSE) +
    ggplot2::labs(
      title = sprintf("%s — clades on the map (marker %s)", taxon, mk),
      subtitle = sprintf("%d geocoded detections across %d clades",
                         nrow(td), length(present))) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "grey92"),
                   plot.title = ggplot2::element_text(face = "bold"),
                   legend.position = "right",
                   legend.text = ggplot2::element_text(size = 6),
                   legend.key.height = ggplot2::unit(0.3, "cm"))
  if (!is.null(basemap)) {
    p <- p + ggplot2::coord_sf(xlim = lon_lim, ylim = lat_lim,
                                default_crs = sf::st_crs(4326),
                                expand = FALSE)
  } else {
    p <- p + ggplot2::coord_quickmap(xlim = lon_lim, ylim = lat_lim)
  }

  pdf_path <- file.path(out_dir,
                        sprintf("map_clades_static_%s_%s.pdf",
                                taxon_clean, mk))
  # Scale page width with clade count to accommodate the legend
  page_w <- min(20, 9 + max(0, length(present) - 20) * 0.05)
  ggplot2::ggsave(pdf_path, p, width = page_w, height = 10)
  msg("    Wrote %s", basename(pdf_path))
  pdf_path
}

# ============================================================================
# 11. Command-line interface
# ============================================================================

# Lightweight CLI so the script can be invoked via Rscript. Recognises
# --rank, --taxon, --markers and --refresh-jobs.
.parse_cli_args <- function(args) {
  out <- list(rank = "Family", taxon = NA_character_,
              markers = NULL, force_refresh_jobs = FALSE)
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (a == "--rank")     { out$rank  <- args[i + 1]; i <- i + 2; next }
    if (a == "--taxon")    { out$taxon <- args[i + 1]; i <- i + 2; next }
    if (a == "--markers")  {
      out$markers <- strsplit(args[i + 1], "[,\\s]+", perl = TRUE)[[1]]
      i <- i + 2; next
    }
    if (a == "--refresh-jobs") { out$force_refresh_jobs <- TRUE; i <- i + 1; next }
    if (a %in% c("-h", "--help")) {
      cat("Usage: Rscript build_taxon_tree.R --rank <Rank> --taxon <Name> ",
          "[--markers CI,WV] [--refresh-jobs]\n", sep = "")
      quit(status = 0)
    }
    warning("Unknown argument: ", a); i <- i + 1
  }
  out
}

# Only run the CLI if executed directly (not when source()'d)
if (!interactive() && length(commandArgs(trailingOnly = TRUE)) > 0) {
  args <- .parse_cli_args(commandArgs(trailingOnly = TRUE))
  if (is.na(args$taxon)) stop("--taxon is required.")
  build_taxon_tree(rank = args$rank,
                   taxon = args$taxon,
                   markers = args$markers,
                   force_refresh_jobs = args$force_refresh_jobs)
}
