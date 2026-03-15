#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  if (!requireNamespace("GISAIDR", quietly = TRUE)) {
    stop(
      paste(
        "Package 'GISAIDR' is not installed.",
        "Install it with devtools::install_github(\"Wytamma/GISAIDR\")."
      ),
      call. = FALSE
    )
  }
}))

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  defaults <- list(
    input = "H3N2-WIC/metadata/H3N2_WIC_gisaid_remaining_query.csv",
    out_dir = "H3N2-WIC/output/gisaidr",
    username = Sys.getenv("GISAIDR_USERNAME", unset = ""),
    password = Sys.getenv("GISAIDR_PASSWORD", unset = ""),
    database = Sys.getenv("GISAIDR_DATABASE", unset = ""),
    batch_size = 250L,
    download_sequences = FALSE,
    sequence_batch_size = 1000L
  )

  if (length(args) == 0L) {
    return(defaults)
  }

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument: %s", key), call. = FALSE)
    }
    key <- sub("^--", "", key)
    if (!(key %in% names(defaults))) {
      stop(sprintf("Unknown argument: --%s", key), call. = FALSE)
    }
    if (is.logical(defaults[[key]])) {
      defaults[[key]] <- TRUE
      i <- i + 1L
      next
    }
    if (i == length(args)) {
      stop(sprintf("Missing value for --%s", key), call. = FALSE)
    }
    defaults[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }

  defaults$batch_size <- as.integer(defaults$batch_size)
  defaults$sequence_batch_size <- as.integer(defaults$sequence_batch_size)
  defaults
}

opts <- parse_args(args)

if (!nzchar(opts$username) || !nzchar(opts$password)) {
  stop(
    paste(
      "GISAID credentials are required.",
      "Set GISAIDR_USERNAME and GISAIDR_PASSWORD or pass --username and --password."
    ),
    call. = FALSE
  )
}

read_query_names <- function(path) {
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  parts <- unlist(strsplit(text, "[,\n\r\t]+"))
  unique(trimws(parts[nzchar(trimws(parts))]))
}

normalize_strain <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("[[:space:]]+", "", x)
  x <- gsub("'", "", x, fixed = TRUE)
  x
}

first_present <- function(df, candidates) {
  hits <- candidates[candidates %in% colnames(df)]
  if (length(hits) == 0L) {
    return(NA_character_)
  }
  hits[[1L]]
}

write_fasta_chunks <- function(credentials, accession_ids, out_fasta, batch_size) {
  if (length(accession_ids) == 0L) {
    return(invisible(NULL))
  }

  if (file.exists(out_fasta)) {
    file.remove(out_fasta)
  }

  batches <- split(accession_ids, ceiling(seq_along(accession_ids) / batch_size))
  for (i in seq_along(batches)) {
    batch_ids <- batches[[i]]
    message(sprintf("Downloading sequences batch %d/%d (%d accessions)", i, length(batches), length(batch_ids)))
    res <- GISAIDR::download_files(
      credentials = credentials,
      list_of_accession_ids = batch_ids,
      sequences = TRUE
    )
    fasta_text <- res$sequences
    cat(fasta_text, file = out_fasta, append = TRUE)
    if (!endsWith(fasta_text, "\n")) {
      cat("\n", file = out_fasta, append = TRUE)
    }
  }
}

query_names <- read_query_names(opts$input)
if (length(query_names) == 0L) {
  stop(sprintf("No query strain names found in %s", opts$input), call. = FALSE)
}

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

login_args <- list(username = opts$username, password = opts$password)
if (nzchar(opts$database)) {
  login_args$database <- opts$database
}

message("Logging into GISAID via GISAIDR")
credentials <- do.call(GISAIDR::login, login_args)

batches <- split(query_names, ceiling(seq_along(query_names) / opts$batch_size))
all_hits <- list()
for (i in seq_along(batches)) {
  batch <- batches[[i]]
  message(sprintf("Querying batch %d/%d (%d strain names)", i, length(batches), length(batch)))
  hits <- tryCatch(
    GISAIDR::query(
      credentials = credentials,
      text = paste(batch, collapse = "\n"),
      load_all = TRUE
    ),
    error = function(e) {
      message(sprintf("Batch %d failed: %s", i, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(hits) || nrow(hits) == 0L) {
    next
  }
  hits$query_batch <- i
  all_hits[[length(all_hits) + 1L]] <- hits
}

if (length(all_hits) == 0L) {
  stop("No GISAID query hits were returned.", call. = FALSE)
}

hits_df <- do.call(rbind, all_hits)
strain_col <- first_present(hits_df, c("strain", "name"))
accession_col <- first_present(hits_df, c("accession_id", "gisaid_epi_isl"))
genbank_col <- first_present(hits_df, c("genbank_accession", "genbank"))
date_col <- first_present(hits_df, c("date", "collection_date"))
segment_col <- first_present(hits_df, c("segment", "gene"))

if (is.na(strain_col) || is.na(accession_col)) {
  stop("Could not find expected strain/accession columns in GISAIDR query output.", call. = FALSE)
}

query_map <- data.frame(
  query_strain = query_names,
  query_key = normalize_strain(query_names),
  stringsAsFactors = FALSE
)

hits_df$strain_key <- normalize_strain(hits_df[[strain_col]])
matched <- merge(
  query_map,
  hits_df,
  by.x = "query_key",
  by.y = "strain_key",
  all.x = FALSE,
  all.y = FALSE
)

accession_table <- data.frame(
  query_strain = matched$query_strain,
  matched_strain = matched[[strain_col]],
  accession_id = matched[[accession_col]],
  genbank_accession = if (!is.na(genbank_col)) matched[[genbank_col]] else NA_character_,
  collection_date = if (!is.na(date_col)) matched[[date_col]] else NA_character_,
  segment = if (!is.na(segment_col)) matched[[segment_col]] else NA_character_,
  query_batch = matched$query_batch,
  stringsAsFactors = FALSE
)

accession_table <- accession_table[!duplicated(accession_table[c("query_strain", "accession_id")]), ]
accession_table <- accession_table[order(accession_table$query_strain, accession_table$accession_id), ]

matched_queries <- unique(accession_table$query_strain)
unresolved_queries <- query_names[!(query_names %in% matched_queries)]

accession_path <- file.path(opts$out_dir, "wic_gisaidr_accessions.tsv")
unresolved_path <- file.path(opts$out_dir, "wic_gisaidr_unresolved.txt")
fasta_path <- file.path(opts$out_dir, "wic_gisaidr_sequences.fasta")

write.table(
  accession_table,
  file = accession_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
writeLines(unresolved_queries, unresolved_path)

message(sprintf("Queried strains: %d", length(query_names)))
message(sprintf("Matched query strains: %d", length(matched_queries)))
message(sprintf("Resolved accession rows: %d", nrow(accession_table)))
message(sprintf("Unresolved query strains: %d", length(unresolved_queries)))
message(sprintf("Accession table: %s", accession_path))
message(sprintf("Unresolved list: %s", unresolved_path))

if (isTRUE(opts$download_sequences)) {
  accession_ids <- unique(accession_table$accession_id[nzchar(accession_table$accession_id)])
  write_fasta_chunks(credentials, accession_ids, fasta_path, opts$sequence_batch_size)
  message(sprintf("Sequence FASTA: %s", fasta_path))
}
