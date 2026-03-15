#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required.", call. = FALSE)
  }
}))

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  defaults <- list(
    query_input = "H3N2-WIC/metadata/H3N2_WIC_gisaid_remaining_query.csv",
    metadata_glob = "H3N2-WIC/data/raw/gisaid/gisaid_meta_*.xls",
    out_dir = "H3N2-WIC/metadata"
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
    if (i == length(args)) {
      stop(sprintf("Missing value for --%s", key), call. = FALSE)
    }
    defaults[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }

  defaults
}

normalize_name <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("&", "AND", x, fixed = TRUE)
  x <- gsub("['`]", "", x)
  x <- gsub("[[:space:]]+", "", x)
  x <- gsub("[^A-Z0-9]", "", x)
  x
}

read_query_names <- function(path) {
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  parts <- unlist(strsplit(text, "[,\n\r\t]+"))
  unique(trimws(parts[nzchar(trimws(parts))]))
}

read_metadata_table <- function(path) {
  raw <- suppressMessages(
    readxl::read_excel(path, sheet = "Tabelle1", col_names = FALSE, .name_repair = "minimal")
  )
  if (nrow(raw) < 2L) {
    return(data.frame())
  }

  headers <- as.character(unlist(raw[1, ], use.names = FALSE))
  keep <- !is.na(headers) & nzchar(headers)
  headers <- headers[keep]
  dat <- raw[-1, keep, drop = FALSE]
  colnames(dat) <- headers

  needed <- c("Isolate_Id", "Isolate_Name", "Collection_Date", "Location", "Host")
  missing <- setdiff(c("Isolate_Id", "Isolate_Name"), colnames(dat))
  if (length(missing) > 0L) {
    stop(sprintf("%s is missing required columns: %s", basename(path), paste(missing, collapse = ", ")), call. = FALSE)
  }

  keep_cols <- intersect(needed, colnames(dat))
  dat <- as.data.frame(dat[, keep_cols, drop = FALSE], stringsAsFactors = FALSE)
  dat$metadata_file <- basename(path)
  dat
}

opts <- parse_args(args)
dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

query_names <- read_query_names(opts$query_input)
if (length(query_names) == 0L) {
  stop(sprintf("No query names found in %s", opts$query_input), call. = FALSE)
}

metadata_files <- sort(Sys.glob(opts$metadata_glob))
if (length(metadata_files) == 0L) {
  stop(sprintf("No metadata files matched %s", opts$metadata_glob), call. = FALSE)
}

metadata_parts <- lapply(metadata_files, read_metadata_table)
metadata <- do.call(rbind, metadata_parts)

metadata <- metadata[!is.na(metadata$Isolate_Id) & nzchar(metadata$Isolate_Id) &
                     !is.na(metadata$Isolate_Name) & nzchar(metadata$Isolate_Name), , drop = FALSE]

metadata$query_key <- normalize_name(metadata$Isolate_Name)
metadata <- metadata[!duplicated(metadata[c("Isolate_Id", "query_key")]), , drop = FALSE]

query_df <- data.frame(
  query_strain = query_names,
  query_key = normalize_name(query_names),
  stringsAsFactors = FALSE
)

matches <- merge(
  query_df,
  metadata,
  by = "query_key",
  all.x = TRUE,
  sort = FALSE
)

matched_rows <- matches[!is.na(matches$Isolate_Id) & nzchar(matches$Isolate_Id), , drop = FALSE]
matched_rows <- matched_rows[order(matched_rows$query_strain, matched_rows$Isolate_Id), , drop = FALSE]

mapping_path <- file.path(opts$out_dir, "H3N2_WIC_gisaid_isolate_id_mapping.tsv")
unresolved_path <- file.path(opts$out_dir, "H3N2_WIC_gisaid_unresolved_after_metadata.txt")
download_csv_path <- file.path(opts$out_dir, "H3N2_WIC_gisaid_isolate_ids.csv")

write.table(
  matched_rows[, c("query_strain", "Isolate_Name", "Isolate_Id", "Collection_Date", "Location", "Host", "metadata_file")],
  file = mapping_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

matched_queries <- unique(matched_rows$query_strain)
unresolved_queries <- query_df$query_strain[!(query_df$query_strain %in% matched_queries)]
writeLines(unresolved_queries, unresolved_path)

download_ids <- unique(matched_rows$Isolate_Id)
writeLines(paste(download_ids, collapse = ","), download_csv_path)

cat(sprintf("Query strains: %d\n", nrow(query_df)))
cat(sprintf("Matched query strains: %d\n", length(matched_queries)))
cat(sprintf("Matched Isolate_Id rows: %d\n", nrow(matched_rows)))
cat(sprintf("Unique Isolate_Id values: %d\n", length(download_ids)))
cat(sprintf("Unresolved query strains: %d\n", length(unresolved_queries)))
cat(sprintf("Mapping table: %s\n", mapping_path))
cat(sprintf("Unresolved list: %s\n", unresolved_path))
cat(sprintf("One-line Isolate_Id file: %s\n", download_csv_path))
