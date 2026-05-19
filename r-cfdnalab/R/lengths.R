#' Read cfDNAlab fragment-length counts.
#'
#' Loads a `<prefix>.length_counts.tsv.zst` file created with the `cfdna
#' lengths` CLI command.
#'
#' @param path Path to a cfDNAlab length-count `.tsv.zst` file.
#'
#' @return A mode-specific `cfdnalab_length_counts` object.
#' @export
#'
#' @examples
#' \dontrun{
#' lengths <- read_lengths("sample.length_counts.tsv.zst")
#' length_bins(lengths)
#' length_data_frame(lengths, value = "fraction")
#' }
read_lengths <- function(path) {
  path <- cf_validate_length_counts_path(path)
  table <- cf_read_length_counts_table(path)
  if (ncol(table) == 0L) {
    stop("Length-count TSV must contain at least one count column", call. = FALSE)
  }

  # The length-count TSV is self-describing enough for the public loader. Keep
  # settings JSON out of this path so command provenance stays separate from
  # the tabular data contract.
  count_columns <- cf_length_count_columns(names(table))
  bin_metadata <- cf_parse_length_count_columns(count_columns)
  mode <- cf_infer_length_count_mode(names(table), count_columns)
  row_metadata <- cf_length_row_metadata(table, mode, count_columns[[1L]])
  counts <- cf_length_counts_matrix_from_table(table, count_columns)

  if (identical(mode, "global") && nrow(counts) != 1L) {
    stop("Global length-count output must contain exactly one row", call. = FALSE)
  }
  if (!identical(nrow(row_metadata), nrow(counts))) {
    stop("Length-count row metadata does not match count row count", call. = FALSE)
  }

  object <- list(
    path = path,
    mode = mode,
    length_bin_idx0 = bin_metadata$length_bin_idx0,
    length_start_bp = bin_metadata$length_start_bp,
    length_end_bp = bin_metadata$length_end_bp,
    count_column = bin_metadata$count_column,
    counts = counts,
    row_metadata = row_metadata
  )

  class(object) <- c(
    switch(
      mode,
      global = "cfdnalab_global_length_counts",
      windowed = "cfdnalab_windowed_length_counts",
      grouped = "cfdnalab_grouped_length_counts"
    ),
    "cfdnalab_length_counts"
  )
  object
}

#' @export
#' @rdname length_bins
length_bins.cfdnalab_length_counts <- function(x, ...) {
  cf_reject_unused_arguments(...)
  data.frame(
    length_bin_idx = cf_index0_to_r_index(x$length_bin_idx0),
    length_start_bp = x$length_start_bp,
    length_end_bp = x$length_end_bp,
    length_midpoint_bp = (x$length_start_bp + x$length_end_bp) / 2,
    length_width_bp = x$length_end_bp - x$length_start_bp,
    count_column = x$count_column,
    stringsAsFactors = FALSE
  )
}

#' @export
#' @rdname length_bin_idx
#' @param length Fragment length in base pairs.
length_bin_idx.cfdnalab_length_counts <- function(x, length, ...) {
  cf_reject_unused_arguments(...)
  length <- cf_validate_fragment_length(length)
  matches <- which(x$length_start_bp <= length & length < x$length_end_bp)
  if (length(matches) == 0L) {
    stop("No length-count bin contains length ", length, call. = FALSE)
  }
  if (length(matches) > 1L) {
    stop("Multiple length-count bins contain length ", length, call. = FALSE)
  }
  cf_index0_to_r_index(x$length_bin_idx0[[matches]])
}

#' @export
#' @rdname length_counts_matrix
length_counts_matrix.cfdnalab_length_counts <- function(x, ...) {
  cf_reject_unused_arguments(...)
  x$counts
}

#' @export
#' @rdname length_counts_vector
length_counts_vector.cfdnalab_global_length_counts <- function(x, ...) {
  cf_reject_unused_arguments(...)
  stats::setNames(as.numeric(x$counts[1L, ]), x$count_column)
}

#' @export
#' @rdname window_metadata
window_metadata.cfdnalab_windowed_length_counts <- function(x, ...) {
  cf_reject_unused_arguments(...)
  x$row_metadata
}

#' @export
#' @rdname group_metadata
group_metadata.cfdnalab_grouped_length_counts <- function(x, ...) {
  cf_reject_unused_arguments(...)
  x$row_metadata
}

#' @export
#' @rdname group_idx
#' @param group_name Group name to look up.
group_idx.cfdnalab_grouped_length_counts <- function(x, group_name, ...) {
  cf_reject_unused_arguments(...)
  group_name <- cf_validate_scalar_string(group_name, "group_name")
  matched_index <- cf_find_unique_value_index(
    x$row_metadata$group_name,
    group_name,
    "Unknown length-count group name: ",
    "Length-count group name is not unique: "
  )
  x$row_metadata$group_idx[[matched_index]]
}

#' @export
#' @rdname length_data_frame
#' @param value Which value to return:
#'   - `"count"` returns raw counts.
#'   - `"fraction"` returns counts divided by the row total.
#'   - `"density"` returns fractions divided by `length_width_bp`, giving
#'     fraction per base pair so bins with different widths are comparable.
#' @param keep_wide If `TRUE`, return one row per output unit with one value
#'   column per length bin. If `FALSE`, return one row per output unit and
#'   length bin.
#' @param max_blacklisted_fraction Optional maximum `blacklisted_fraction` in
#'   0..1 to retain before reshaping. The default `1.0` keeps all rows.
length_data_frame.cfdnalab_global_length_counts <- function(
  x,
  value = c("count", "fraction", "density"),
  keep_wide = FALSE,
  max_blacklisted_fraction = 1.0,
  ...
) {
  cf_reject_unused_arguments(...)
  row_indices <- cf_length_row_indices_after_blacklist_filter(x, max_blacklisted_fraction)
  cf_length_data_frame_for_rows(x, row_indices, match.arg(value), keep_wide)
}

#' @export
#' @rdname length_data_frame
#' @param window_idx Optional one-based window index vector.
length_data_frame.cfdnalab_windowed_length_counts <- function(
  x,
  window_idx = NULL,
  value = c("count", "fraction", "density"),
  keep_wide = FALSE,
  max_blacklisted_fraction = 1.0,
  ...
) {
  cf_reject_unused_arguments(...)
  row_indices <- if (is.null(window_idx)) {
    seq_len(nrow(x$counts))
  } else {
    cf_validate_r_indices(window_idx, nrow(x$counts), "window_idx")
  }
  row_indices <- cf_apply_length_blacklist_filter(x, row_indices, max_blacklisted_fraction)
  cf_length_data_frame_for_rows(x, row_indices, match.arg(value), keep_wide)
}

#' @export
#' @rdname length_data_frame
#' @param group Optional group name or one-based group index vector. Use either
#'   `group` or `group_idx`, not both.
#' @param group_idx Optional one-based group index vector.
length_data_frame.cfdnalab_grouped_length_counts <- function(
  x,
  group = NULL,
  group_idx = NULL,
  value = c("count", "fraction", "density"),
  keep_wide = FALSE,
  max_blacklisted_fraction = 1.0,
  ...
) {
  cf_reject_unused_arguments(...)
  if (!is.null(group) && !is.null(group_idx)) {
    stop("Use either group or group_idx, not both", call. = FALSE)
  }
  row_indices <- if (!is.null(group)) {
    cf_resolve_length_group_indices(x, group)
  } else if (!is.null(group_idx)) {
    cf_validate_r_indices(group_idx, nrow(x$counts), "group_idx")
  } else {
    seq_len(nrow(x$counts))
  }
  row_indices <- cf_apply_length_blacklist_filter(x, row_indices, max_blacklisted_fraction)
  cf_length_data_frame_for_rows(x, row_indices, match.arg(value), keep_wide)
}

#' Print a length-count object.
#'
#' @param x A cfDNAlab length-count object.
#' @param ... Ignored.
#'
#' @return Invisibly returns `x`.
#' @export
#' @keywords internal
print.cfdnalab_length_counts <- function(x, ...) {
  cat("<cfDNAlab length counts>\n")
  cat("Path: ", x$path, "\n", sep = "")
  cat("Mode: ", x$mode, "\n", sep = "")
  cat("Rows: ", nrow(x$counts), "\n", sep = "")
  cat("Length bins: ", length(x$length_bin_idx0), "\n", sep = "")
  invisible(x)
}

#' Validate a length-count TSV path.
#'
#' @param path User-supplied path.
#'
#' @return A normalized path string.
#' @noRd
cf_validate_length_counts_path <- function(path) {
  if (length(path) != 1L || !is.character(path) || is.na(path)) {
    stop("Length-count path must be a single path string", call. = FALSE)
  }
  path <- normalizePath(path, mustWork = FALSE)
  if (!file.exists(path)) {
    stop("Length-count TSV does not exist: ", path, call. = FALSE)
  }
  if (dir.exists(path)) {
    stop("Length-count path exists but is a directory: ", path, call. = FALSE)
  }
  if (!grepl("\\.tsv(\\.zst)?$", path, ignore.case = TRUE)) {
    stop("Length-count path must end in '.tsv' or '.tsv.zst': ", path, call. = FALSE)
  }
  path
}

#' Read a length-count TSV table.
#'
#' @param path Validated TSV path.
#'
#' @return A base data frame.
#' @noRd
cf_read_length_counts_table <- function(path) {
  if (grepl("\\.zst$", path, ignore.case = TRUE)) {
    zstd <- Sys.which("zstd")
    if (!nzchar(zstd)) {
      stop("Reading .tsv.zst length-count files requires the zstd command-line tool", call. = FALSE)
    }
    # data.table::fread() is the fast TSV parser here; zstd streams compressed
    # outputs to it without materializing an intermediate decompressed file.
    command <- paste(shQuote(zstd), "-dc", shQuote(path))
    return(data.table::fread(cmd = command, data.table = FALSE, check.names = FALSE))
  }
  # Plain TSV support is mainly for small fixtures and local debugging.
  data.table::fread(path, data.table = FALSE, check.names = FALSE)
}

#' Return and validate length-count columns.
#'
#' @param column_names TSV column names.
#'
#' @return Character vector of count columns.
#' @noRd
cf_length_count_columns <- function(column_names) {
  count_columns <- grep("^count_[0-9]+(_[0-9]+)?$", column_names, value = TRUE)
  if (length(count_columns) == 0L) {
    stop("Length-count TSV must contain count columns named count_<length> or count_<start>_<end>", call. = FALSE)
  }
  if (anyDuplicated(column_names)) {
    stop("Length-count TSV column names must be unique", call. = FALSE)
  }
  first_count_column <- match(count_columns[[1L]], column_names)
  expected_count_columns <- column_names[first_count_column:length(column_names)]
  if (!identical(count_columns, expected_count_columns)) {
    stop("Length-count TSV count columns must be contiguous and follow metadata columns", call. = FALSE)
  }
  count_columns
}

#' Infer the length-count output mode.
#'
#' @param column_names TSV column names.
#' @param count_columns Count column names.
#'
#' @return One of `global`, `windowed`, or `grouped`.
#' @noRd
cf_infer_length_count_mode <- function(column_names, count_columns) {
  first_count_column <- match(count_columns[[1L]], column_names)
  metadata_columns <- column_names[seq_len(first_count_column - 1L)]
  if (identical(metadata_columns, character(0))) {
    return("global")
  }
  if (identical(metadata_columns, c("chrom", "start", "end")) ||
    identical(metadata_columns, c("chrom", "start", "end", "blacklisted_fraction"))) {
    return("windowed")
  }
  if (identical(metadata_columns, c("group_name", "eligible_windows")) ||
    identical(metadata_columns, c("group_name", "eligible_windows", "blacklisted_fraction"))) {
    return("grouped")
  }
  stop(
    "Could not infer length-count output mode from metadata columns: ",
    paste(metadata_columns, collapse = ", "),
    call. = FALSE
  )
}

#' Parse length-bin metadata from count column names.
#'
#' @param count_columns Count column names.
#'
#' @return A data frame with internal zero-based bin indices.
#' @noRd
cf_parse_length_count_columns <- function(count_columns) {
  starts <- integer(length(count_columns))
  ends <- integer(length(count_columns))
  for (column_index in seq_along(count_columns)) {
    suffix <- sub("^count_", "", count_columns[[column_index]])
    parts <- strsplit(suffix, "_", fixed = TRUE)[[1L]]
    if (length(parts) == 1L) {
      start <- cf_parse_nonnegative_integer_string(parts[[1L]], count_columns[[column_index]])
      end <- start + 1L
    } else if (length(parts) == 2L) {
      start <- cf_parse_nonnegative_integer_string(parts[[1L]], count_columns[[column_index]])
      end <- cf_parse_nonnegative_integer_string(parts[[2L]], count_columns[[column_index]])
    } else {
      stop("Invalid length-count column name: ", count_columns[[column_index]], call. = FALSE)
    }
    starts[[column_index]] <- start
    ends[[column_index]] <- end
  }
  cf_validate_half_open_intervals(starts, ends, "length bin start", "length bin end")
  intervals <- paste(starts, ends, sep = "-")
  if (anyDuplicated(intervals)) {
    stop("Length-count TSV contains duplicate length bins", call. = FALSE)
  }
  data.frame(
    length_bin_idx0 = seq_along(count_columns) - 1L,
    length_start_bp = starts,
    length_end_bp = ends,
    count_column = count_columns,
    stringsAsFactors = FALSE
  )
}

#' Parse a non-negative integer string.
#'
#' @param value String value.
#' @param name Value name for error messages.
#'
#' @return An integer.
#' @noRd
cf_parse_nonnegative_integer_string <- function(value, name) {
  if (!grepl("^[0-9]+$", value)) {
    stop("Invalid non-negative integer in ", name, call. = FALSE)
  }
  parsed <- as.numeric(value)
  if (!is.finite(parsed) || parsed > .Machine$integer.max) {
    stop("Length bin bound in ", name, " must fit in R integer range", call. = FALSE)
  }
  as.integer(parsed)
}

#' Build row metadata for a length-count object.
#'
#' @param table TSV table.
#' @param mode Inferred mode.
#' @param first_count_column First count column name.
#'
#' @return A data frame.
#' @noRd
cf_length_row_metadata <- function(table, mode, first_count_column) {
  metadata_columns <- names(table)[seq_len(match(first_count_column, names(table)) - 1L)]
  metadata <- table[metadata_columns]
  if (identical(mode, "global")) {
    return(data.frame(row.names = seq_len(nrow(table))))
  }
  if (identical(mode, "windowed")) {
    cf_validate_character_vector(metadata$chrom, "chrom")
    cf_validate_half_open_intervals(metadata$start, metadata$end, "start", "end")
    out <- data.frame(
      window_idx = seq_len(nrow(table)),
      chrom = as.character(metadata$chrom),
      start = as.integer(metadata$start),
      end = as.integer(metadata$end),
      stringsAsFactors = FALSE
    )
  } else {
    cf_validate_character_vector(metadata$group_name, "group_name")
    cf_validate_nonnegative_integer_vector(metadata$eligible_windows, "eligible_windows")
    out <- data.frame(
      group_idx = seq_len(nrow(table)),
      group_name = as.character(metadata$group_name),
      eligible_windows = as.integer(metadata$eligible_windows),
      stringsAsFactors = FALSE
    )
  }
  if ("blacklisted_fraction" %in% names(metadata)) {
    cf_validate_fraction_vector(metadata$blacklisted_fraction, "blacklisted_fraction")
    out$blacklisted_fraction <- metadata$blacklisted_fraction
  }
  out
}

#' Build a numeric count matrix from a length-count TSV.
#'
#' @param table TSV table.
#' @param count_columns Count column names.
#'
#' @return A numeric matrix.
#' @noRd
cf_length_counts_matrix_from_table <- function(table, count_columns) {
  for (count_column in count_columns) {
    cf_validate_nonnegative_numeric_vector(table[[count_column]], count_column)
  }
  counts <- as.matrix(table[count_columns])
  storage.mode(counts) <- "numeric"
  colnames(counts) <- count_columns
  rownames(counts) <- NULL
  counts
}

#' Validate a character vector.
#'
#' @param values Values to validate.
#' @param value_name Human-readable value name.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_character_vector <- function(values, value_name) {
  if (!is.character(values) || any(is.na(values)) || any(!nzchar(values))) {
    stop(value_name, " must contain non-empty character strings", call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate a fragment length.
#'
#' @param length Fragment length value.
#'
#' @return The validated integer length.
#' @noRd
cf_validate_fragment_length <- function(length) {
  if (
    length(length) != 1L ||
      !is.numeric(length) ||
      is.na(length) ||
      !is.finite(length) ||
      length > .Machine$integer.max ||
      length != as.integer(length) ||
      length < 0L
  ) {
    stop("Fragment length must be a single non-negative integer", call. = FALSE)
  }
  as.integer(length)
}

#' Validate public R indices.
#'
#' @param indices User-supplied indices.
#' @param size Axis size.
#' @param name Human-readable index name.
#'
#' @return Integer vector of one-based indices.
#' @noRd
cf_validate_r_indices <- function(indices, size, name) {
  if (
    !is.numeric(indices) ||
      any(is.na(indices)) ||
      any(!is.finite(indices)) ||
      any(indices > .Machine$integer.max) ||
      any(indices != as.integer(indices))
  ) {
    stop(name, " must contain integer values", call. = FALSE)
  }
  indices <- as.integer(indices)
  if (any(indices < 1L | indices > size)) {
    stop(name, " contains values outside 1..", size, call. = FALSE)
  }
  indices
}

#' Resolve grouped length-count selectors to row indices.
#'
#' @param x Grouped length-count object.
#' @param group Group names or one-based indices.
#'
#' @return One-based row indices.
#' @noRd
cf_resolve_length_group_indices <- function(x, group) {
  if (is.character(group)) {
    matches <- match(group, x$row_metadata$group_name)
    if (any(is.na(matches))) {
      stop("Unknown length-count group name: ", sQuote(group[is.na(matches)][[1L]]), call. = FALSE)
    }
    if (anyDuplicated(x$row_metadata$group_name)) {
      stop("Length-count group names are not unique", call. = FALSE)
    }
    return(as.integer(matches))
  }
  cf_validate_r_indices(group, nrow(x$counts), "group")
}

#' Apply a blacklist fraction filter to row indices.
#'
#' @param x Length-count object.
#' @param row_indices One-based row indices.
#' @param max_blacklisted_fraction Optional maximum blacklist fraction.
#'
#' @return Filtered one-based row indices.
#' @noRd
cf_apply_length_blacklist_filter <- function(x, row_indices, max_blacklisted_fraction) {
  cf_apply_row_blacklist_filter(x$row_metadata, row_indices, max_blacklisted_fraction)
}

#' Return row indices after a blacklist filter.
#'
#' @param x Length-count object.
#' @param max_blacklisted_fraction Optional maximum blacklist fraction.
#'
#' @return One-based row indices.
#' @noRd
cf_length_row_indices_after_blacklist_filter <- function(x, max_blacklisted_fraction) {
  cf_apply_length_blacklist_filter(x, seq_len(nrow(x$counts)), max_blacklisted_fraction)
}

#' Build a long or wide length data frame for selected rows.
#'
#' @param x Length-count object.
#' @param row_indices One-based row indices.
#' @param value Value type.
#' @param keep_wide Whether to keep wide shape.
#'
#' @return A data frame.
#' @noRd
cf_length_data_frame_for_rows <- function(x, row_indices, value, keep_wide) {
  cf_validate_scalar_logical(keep_wide, "keep_wide")
  values <- cf_length_value_matrix(x, row_indices, value)
  if (isTRUE(keep_wide)) {
    return(cf_length_wide_data_frame(x, row_indices, values, value))
  }
  cf_length_long_data_frame(x, row_indices, values, value)
}

#' Validate a scalar logical.
#'
#' @param value Logical value.
#' @param name Human-readable value name.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_scalar_logical <- function(value, name) {
  if (length(value) != 1L || !is.logical(value) || is.na(value)) {
    stop(name, " must be TRUE or FALSE", call. = FALSE)
  }
  invisible(TRUE)
}

#' Compute count-derived values for selected rows.
#'
#' @param x Length-count object.
#' @param row_indices One-based row indices.
#' @param value Value type.
#'
#' @return Numeric matrix.
#' @noRd
cf_length_value_matrix <- function(x, row_indices, value) {
  counts <- x$counts[row_indices, , drop = FALSE]
  if (identical(value, "count")) {
    return(counts)
  }

  row_totals <- rowSums(counts)
  fractions <- counts
  positive_rows <- row_totals > 0
  fractions[positive_rows, ] <- counts[positive_rows, , drop = FALSE] / row_totals[positive_rows]
  fractions[!positive_rows, ] <- NA_real_

  if (identical(value, "fraction")) {
    return(fractions)
  }

  widths <- x$length_end_bp - x$length_start_bp
  sweep(fractions, 2L, widths, "/")
}

#' Build a wide length data frame.
#'
#' @param x Length-count object.
#' @param row_indices One-based row indices.
#' @param values Value matrix.
#' @param value Value type.
#'
#' @return A data frame.
#' @noRd
cf_length_wide_data_frame <- function(x, row_indices, values, value) {
  value_columns <- cf_length_value_column_names(x$count_column, value)
  colnames(values) <- value_columns
  value_frame <- as.data.frame(values, check.names = FALSE)
  metadata <- x$row_metadata[row_indices, , drop = FALSE]
  data.frame(metadata, value_frame, row.names = NULL, stringsAsFactors = FALSE)
}

#' Build a long length data frame.
#'
#' @param x Length-count object.
#' @param row_indices One-based row indices.
#' @param values Value matrix.
#' @param value Value type.
#'
#' @return A data frame.
#' @noRd
cf_length_long_data_frame <- function(x, row_indices, values, value) {
  num_rows <- length(row_indices)
  num_bins <- length(x$length_bin_idx0)
  metadata <- x$row_metadata[row_indices, , drop = FALSE]
  metadata <- metadata[rep(seq_len(num_rows), each = num_bins), , drop = FALSE]
  bins <- length_bins(x)
  bins <- bins[rep(seq_len(num_bins), times = num_rows), , drop = FALSE]
  out <- data.frame(metadata, bins, stringsAsFactors = FALSE, row.names = NULL)
  out[[value]] <- as.vector(t(values))
  out
}

#' Build wide value column names.
#'
#' @param count_columns Source count columns.
#' @param value Value type.
#'
#' @return Character vector.
#' @noRd
cf_length_value_column_names <- function(count_columns, value) {
  if (identical(value, "count")) {
    return(count_columns)
  }
  sub("^count_", paste0(value, "_"), count_columns)
}
