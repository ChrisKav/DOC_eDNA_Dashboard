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
#' @return Invisibly returns a list with all intermediate objects so the
#'         caller can inspect or post-process.
#' @export
build_taxon_tree <- function(rank = "Family",
                             taxon,
                             markers = NULL,
                             min_seqs = 3,
                             force_refresh_jobs = FALSE,
                             out_root = "Output/Phylogenetics") {

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
    seq_tables[[i]] <- tryCatch(
      extract_sequences_one_job(ok$local_path[i], rec_sub, rank, taxon),
      error = function(e) {
        warning("Job ", jid, " sequence extraction failed: ", e$message)
        NULL
      })
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

  msg("\nDone. Outputs in %s", out_dir)
  invisible(list(records = records_hit,
                 jobs = jobs,
                 downloads = dls,
                 merged_sequences = merged,
                 marker_results = marker_results,
                 trees = trees,
                 candidates = candidates,
                 output_dir = out_dir))
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
