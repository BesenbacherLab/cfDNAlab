#' Supported end-motif storage modes.
#'
#' These values mirror the cfDNAlab Zarr `storage_mode` root attribute.
#'
#' @noRd
END_MOTIF_VALID_STORAGE_MODES <- c("dense", "sparse_coo")

#' Supported end-motif row modes.
#'
#' These values mirror the cfDNAlab Zarr `row_mode` root attribute.
#'
#' @noRd
END_MOTIF_VALID_ROW_MODES <- c("global", "size", "bed", "grouped_bed")

#' Read cfDNAlab end-motif counts.
#'
#' Loads a `<prefix>.end_motifs.zarr` store created with the \code{cfdna ends}
#' CLI tool from the main \code{cfDNAlab} rust package.
#' It validates the cfDNAlab schema, row metadata, motif metadata,
#' and dense or sparse count layout.
#'
#' @param path Path to a cfDNAlab end-motif `.zarr` directory.
#'
#' @return One of `cfdnalab_global_end_motif_counts`,
#'   `cfdnalab_windowed_end_motif_counts`, or
#'   `cfdnalab_grouped_end_motif_counts`, depending on the row mode.
#' @export
#'
#' @examples
#' \dontrun{
#' ends <- read_end_motifs("sample.end_motifs.zarr")
#' motifs(ends)
#' sparse_counts_matrix(ends)
#' }
read_end_motifs <- function(path) {
  path <- cf_validate_zarr_store_path(path, "End-motif")
  root_attributes <- cf_root_attributes(path)
  cf_validate_schema(root_attributes, "end_motif_counts", "end-motif")

  storage_mode <- cf_validate_allowed_string(
    root_attributes$storage_mode,
    END_MOTIF_VALID_STORAGE_MODES,
    "end-motif storage mode"
  )
  row_mode <- cf_validate_allowed_string(
    root_attributes$row_mode,
    END_MOTIF_VALID_ROW_MODES,
    "end-motif row mode"
  )

  store <- cf_open_zarr(path, "end-motif")
  cf_required_arrays(store, cf_end_motif_required_arrays(storage_mode, row_mode), "End-motif")
  cf_validate_dimension_names(path, "motif_index", "motif")
  cf_validate_dimension_names(path, "motif_byte", "motif_byte")
  cf_validate_dimension_names(path, "motif_ascii", c("motif", "motif_byte"))
  cf_validate_dimension_names(path, "row", "row")

  motif_axis <- cf_read_vector(store, "motif_index", "End-motif")
  motif_byte <- cf_read_vector(store, "motif_byte", "End-motif")
  row <- cf_read_vector(store, "row", "End-motif")
  motif_ascii <- cf_read_array(store, "motif_ascii", "End-motif")
  motif <- cf_decode_motif_ascii(motif_ascii, length(motif_axis), length(motif_byte))

  cf_validate_axis(motif_axis, "motif_index")
  cf_validate_axis(motif_byte, "motif_byte")
  cf_validate_axis(row, "row")
  if (identical(row_mode, "global") && length(row) != 1L) {
    stop("global end-motif stores must contain exactly one row", call. = FALSE)
  }

  counts <- NULL
  sparse <- NULL
  if (identical(storage_mode, "dense")) {
    cf_validate_dimension_names(path, "counts", c("row", "motif"))
    counts <- cf_get_array(store, "counts", "End-motif")
    expected_shape <- c(length(row), length(motif_axis))
    if (!identical(cf_array_shape(counts), expected_shape)) {
      stop(
        "dense counts shape does not match row and motif axes: counts=",
        paste(cf_array_shape(counts), collapse = "x"),
        ", coordinates=",
        paste(expected_shape, collapse = "x"),
        call. = FALSE
      )
    }
  } else {
    sparse <- list(
      row_idx0 = cf_read_vector(store, "sparse/row", "End-motif"),
      motif_idx0 = cf_read_vector(store, "sparse/motif", "End-motif"),
      count = cf_read_vector(store, "sparse/count", "End-motif"),
      shape = cf_read_vector(store, "sparse/shape", "End-motif")
    )
    cf_validate_dimension_names(path, "sparse/row", "nnz")
    cf_validate_dimension_names(path, "sparse/motif", "nnz")
    cf_validate_dimension_names(path, "sparse/count", "nnz")
    cf_validate_dimension_names(path, "sparse/shape", "sparse_dimension")
    cf_validate_dimension_names(path, "sparse/sparse_dimension", "sparse_dimension")
    sparse_dimension_labels <- cf_read_labels(
      path,
      "sparse/sparse_dimension",
      "sparse_dimension_name",
      2L
    )
    if (!identical(sparse_dimension_labels, c("row", "motif"))) {
      stop("sparse_dimension labels must be row, motif", call. = FALSE)
    }
    cf_validate_nonnegative_integer_vector(sparse$shape, "sparse/shape")
    if (any(sparse$shape > .Machine$integer.max)) {
      stop("sparse/shape values must fit in R integer range", call. = FALSE)
    }
    if (!identical(as.integer(sparse$shape), c(length(row), length(motif_axis)))) {
      stop("sparse/shape does not match row and motif axes", call. = FALSE)
    }
    cf_validate_same_length(sparse$motif_idx0, sparse$row_idx0, "sparse/motif", "sparse/row")
    cf_validate_same_length(sparse$count, sparse$row_idx0, "sparse/count", "sparse/row")
    cf_validate_index_vector(sparse$row_idx0, length(row), "sparse/row")
    cf_validate_index_vector(sparse$motif_idx0, length(motif_axis), "sparse/motif")
    cf_validate_nonnegative_numeric_vector(sparse$count, "sparse/count")
  }

  row_metadata <- cf_read_end_motif_row_metadata(
    path = path,
    store = store,
    row = row,
    row_mode = row_mode
  )

  object <- list(
    path = path,
    store = store,
    root_attributes = root_attributes,
    storage_mode = storage_mode,
    row_mode = row_mode,
    motif_idx0 = as.integer(motif_axis),
    motif = motif,
    row_idx0 = as.integer(row),
    counts = counts,
    sparse = sparse,
    row_metadata = row_metadata
  )

  class(object) <- c(
    switch(
      row_mode,
      global = "cfdnalab_global_end_motif_counts",
      size = "cfdnalab_windowed_end_motif_counts",
      bed = "cfdnalab_windowed_end_motif_counts",
      grouped_bed = "cfdnalab_grouped_end_motif_counts"
    ),
    "cfdnalab_end_motif_counts",
    "cfdnalab_zarr_store"
  )
  object
}

#' Return arrays required by an end-motif store.
#'
#' @param storage_mode End-motif storage mode.
#' @param row_mode End-motif row mode.
#'
#' @return Character vector of required array paths.
#' @noRd
cf_end_motif_required_arrays <- function(storage_mode, row_mode) {
  required <- c("motif_index", "motif_byte", "motif_ascii", "row")
  if (identical(storage_mode, "dense")) {
    required <- c(required, "counts")
  } else {
    required <- c(
      required,
      "sparse/row",
      "sparse/motif",
      "sparse/count",
      "sparse/shape",
      "sparse/sparse_dimension"
    )
  }
  if (row_mode %in% c("size", "bed")) {
    required <- c(
      required,
      "chromosome",
      "row_chromosome",
      "row_start_bp",
      "row_end_bp",
      "blacklisted_fraction"
    )
  } else if (identical(row_mode, "grouped_bed")) {
    required <- c(required, "group", "eligible_windows", "blacklisted_fraction")
  }
  required
}

#' Read end-motif row metadata.
#'
#' @param path Path to a Zarr store.
#' @param store Open Zarr store.
#' @param row Row-axis values.
#' @param row_mode End-motif row mode.
#'
#' @return A data frame describing count rows.
#' @noRd
cf_read_end_motif_row_metadata <- function(path, store, row, row_mode) {
  if (identical(row_mode, "global")) {
    cf_validate_dimension_names(path, "row", "row")
    labels <- cf_read_labels(path, "row", "row_label", length(row))
    return(data.frame(row_label = labels, stringsAsFactors = FALSE))
  }

  if (row_mode %in% c("size", "bed")) {
    cf_validate_dimension_names(path, "chromosome", "chromosome")
    cf_validate_dimension_names(path, "row_chromosome", "row")
    cf_validate_dimension_names(path, "row_start_bp", "row")
    cf_validate_dimension_names(path, "row_end_bp", "row")
    cf_validate_dimension_names(path, "blacklisted_fraction", "row")
    chromosome <- cf_read_vector(store, "chromosome", "End-motif")
    chromosome_name <- cf_read_labels(path, "chromosome", "chromosome_name", length(chromosome))
    row_chromosome <- cf_read_vector(store, "row_chromosome", "End-motif")
    row_start_bp <- cf_read_vector(store, "row_start_bp", "End-motif")
    row_end_bp <- cf_read_vector(store, "row_end_bp", "End-motif")
    blacklisted_fraction <- cf_read_vector(store, "blacklisted_fraction", "End-motif")
    cf_validate_axis(chromosome, "chromosome")
    cf_validate_same_length(row_chromosome, row, "row_chromosome", "row")
    cf_validate_same_length(row_start_bp, row, "row_start_bp", "row")
    cf_validate_same_length(row_end_bp, row, "row_end_bp", "row")
    cf_validate_same_length(blacklisted_fraction, row, "blacklisted_fraction", "row")
    cf_validate_index_vector(row_chromosome, length(chromosome), "row_chromosome")
    cf_validate_half_open_intervals(row_start_bp, row_end_bp, "row_start_bp", "row_end_bp")
    cf_validate_fraction_vector(blacklisted_fraction, "blacklisted_fraction")
    return(data.frame(
      window_idx = cf_index0_to_r_index(row),
      chromosome_idx = cf_index0_to_r_index(row_chromosome),
      chromosome_name = chromosome_name[as.integer(row_chromosome) + 1L],
      window_start_bp = row_start_bp,
      window_end_bp = row_end_bp,
      blacklisted_fraction = blacklisted_fraction,
      stringsAsFactors = FALSE
    ))
  }

  group <- cf_read_vector(store, "group", "End-motif")
  cf_validate_dimension_names(path, "group", "row")
  cf_validate_dimension_names(path, "eligible_windows", "row")
  cf_validate_dimension_names(path, "blacklisted_fraction", "row")
  group_name <- cf_read_labels(path, "group", "group_name", length(group))
  eligible_windows <- cf_read_vector(store, "eligible_windows", "End-motif")
  blacklisted_fraction <- cf_read_vector(store, "blacklisted_fraction", "End-motif")
  cf_validate_axis(group, "group")
  cf_validate_same_length(group, row, "group", "row")
  cf_validate_same_length(group_name, row, "group_name", "row")
  cf_validate_same_length(eligible_windows, row, "eligible_windows", "row")
  cf_validate_same_length(blacklisted_fraction, row, "blacklisted_fraction", "row")
  cf_validate_nonnegative_integer_vector(eligible_windows, "eligible_windows")
  cf_validate_fraction_vector(blacklisted_fraction, "blacklisted_fraction")
  data.frame(
    group_idx = cf_index0_to_r_index(group),
    group_name = group_name,
    eligible_windows = eligible_windows,
    blacklisted_fraction = blacklisted_fraction,
    stringsAsFactors = FALSE
  )
}

#' @export
#' @rdname storage_mode
storage_mode.cfdnalab_end_motif_counts <- function(x, ...) {
  cf_reject_unused_arguments(...)
  x$storage_mode
}

#' @export
#' @rdname row_mode
row_mode.cfdnalab_end_motif_counts <- function(x, ...) {
  cf_reject_unused_arguments(...)
  x$row_mode
}

#' @export
#' @rdname motifs
motifs.cfdnalab_end_motif_counts <- function(x, ...) {
  cf_reject_unused_arguments(...)
  data.frame(
    motif_idx = cf_index0_to_r_index(x$motif_idx0),
    motif = x$motif,
    stringsAsFactors = FALSE
  )
}

#' @export
#' @rdname motif_idx
#' @param motif Motif label to look up.
motif_idx.cfdnalab_end_motif_counts <- function(x, motif, ...) {
  cf_reject_unused_arguments(...)
  motif <- cf_validate_scalar_string(motif, "motif")
  matches <- which(x$motif == motif)
  if (length(matches) == 0L) {
    stop("Unknown end-motif label: ", sQuote(motif), call. = FALSE)
  }
  if (length(matches) > 1L) {
    stop("End-motif label is not unique: ", sQuote(motif), call. = FALSE)
  }
  cf_index0_to_r_index(x$motif_idx0[[matches]])
}

#' @export
#' @rdname has_motif
#' @param motif Motif label to test.
has_motif.cfdnalab_end_motif_counts <- function(x, motif, ...) {
  cf_reject_unused_arguments(...)
  motif <- cf_validate_scalar_string(motif, "motif")
  any(x$motif == motif)
}

#' @export
#' @rdname windows
windows.cfdnalab_windowed_end_motif_counts <- function(x, ...) {
  cf_reject_unused_arguments(...)
  x$row_metadata
}

#' @export
#' @rdname group_metadata
group_metadata.cfdnalab_grouped_end_motif_counts <- function(x, ...) {
  cf_reject_unused_arguments(...)
  x$row_metadata
}

#' @export
#' @rdname group_idx
#' @param group_name Group name to look up.
group_idx.cfdnalab_grouped_end_motif_counts <- function(x, group_name, ...) {
  cf_reject_unused_arguments(...)
  group_name <- cf_validate_scalar_string(group_name, "group_name")
  matches <- which(x$row_metadata$group_name == group_name)
  if (length(matches) == 0L) {
    stop("Unknown end-motif group name: ", sQuote(group_name), call. = FALSE)
  }
  if (length(matches) > 1L) {
    stop("End-motif group name is not unique: ", sQuote(group_name), call. = FALSE)
  }
  x$row_metadata$group_idx[[matches]]
}

#' @export
#' @rdname sparse_counts_matrix
sparse_counts_matrix.cfdnalab_end_motif_counts <- function(x, ...) {
  cf_reject_unused_arguments(...)
  if (identical(x$storage_mode, "dense")) {
    return(Matrix::Matrix(x$counts$read(), sparse = TRUE))
  }
  Matrix::sparseMatrix(
    i = cf_index0_to_r_index(x$sparse$row_idx0),
    j = cf_index0_to_r_index(x$sparse$motif_idx0),
    x = as.numeric(x$sparse$count),
    dims = as.integer(x$sparse$shape)
  )
}

#' @export
#' @rdname dense_counts_matrix
#' @param allow_densify If `TRUE`, allow sparse stores to be converted to a dense
#'   in-memory matrix. Sparse stores error by default.
dense_counts_matrix.cfdnalab_end_motif_counts <- function(x, allow_densify = FALSE, ...) {
  cf_reject_unused_arguments(...)
  if (identical(x$storage_mode, "dense")) {
    return(cf_read_array(x$store, "counts", "End-motif"))
  }
  if (!isTRUE(allow_densify)) {
    stop(
      "This end-motif store is sparse. Use sparse_counts_matrix() or set allow_densify = TRUE.",
      call. = FALSE
    )
  }
  as.matrix(sparse_counts_matrix(x))
}

#' @export
#' @rdname dense_counts_vector
#' @param allow_densify If `TRUE`, allow sparse stores to be converted to dense
#'   in memory before returning the vector.
dense_counts_vector.cfdnalab_global_end_motif_counts <- function(
  x,
  allow_densify = FALSE,
  ...
) {
  cf_reject_unused_arguments(...)
  counts <- as.vector(dense_counts_matrix(x, allow_densify = allow_densify)[1L, ])
  stats::setNames(counts, x$motif)
}

#' @export
#' @rdname dense_data_frame
#' @param allow_densify If `TRUE`, allow sparse stores to be converted to dense
#'   in memory before returning the data frame.
dense_data_frame.cfdnalab_global_end_motif_counts <- function(
  x,
  allow_densify = FALSE,
  ...
) {
  cf_reject_unused_arguments(...)
  data.frame(
    motifs(x),
    count = unname(dense_counts_vector(x, allow_densify = allow_densify)),
    stringsAsFactors = FALSE
  )
}

#' @export
#' @rdname dense_data_frame_for_window
#' @param window_idx One-based window index.
#' @param allow_densify If `TRUE`, allow sparse stores to be converted to dense
#'   in memory before returning the data frame.
dense_data_frame_for_window.cfdnalab_windowed_end_motif_counts <- function(
  x,
  window_idx,
  allow_densify = FALSE,
  ...
) {
  cf_reject_unused_arguments(...)
  window_idx <- cf_validate_r_index(window_idx, length(x$row_idx0), "window_idx")
  counts <- dense_counts_matrix(x, allow_densify = allow_densify)[window_idx, ]
  metadata <- x$row_metadata[window_idx, , drop = FALSE]
  data.frame(
    motifs(x),
    count = as.vector(counts),
    metadata[rep(1L, length(x$motif)), , drop = FALSE],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' @export
#' @rdname dense_data_frame_for_group
#' @param group Group name or one-based group index.
#' @param allow_densify If `TRUE`, allow sparse stores to be converted to dense
#'   in memory before returning the data frame.
dense_data_frame_for_group.cfdnalab_grouped_end_motif_counts <- function(
  x,
  group,
  allow_densify = FALSE,
  ...
) {
  cf_reject_unused_arguments(...)
  group_idx0 <- cf_resolve_end_motif_group_idx0(x, group)
  group_r_index <- cf_index0_to_r_index(group_idx0)
  counts <- dense_counts_matrix(x, allow_densify = allow_densify)[group_r_index, ]
  metadata <- x$row_metadata[group_r_index, , drop = FALSE]
  data.frame(
    motifs(x),
    count = as.vector(counts),
    metadata[rep(1L, length(x$motif)), , drop = FALSE],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' @export
#' @rdname dense_data_frame_for_motif
#' @param motif Motif label.
#' @param allow_densify If `TRUE`, allow sparse stores to be converted to dense
#'   in memory before returning the data frame.
dense_data_frame_for_motif.cfdnalab_end_motif_counts <- function(
  x,
  motif,
  allow_densify = FALSE,
  ...
) {
  cf_reject_unused_arguments(...)
  motif_idx0 <- cf_resolve_end_motif_idx0(x, motif)
  counts <- dense_counts_matrix(x, allow_densify = allow_densify)[, cf_index0_to_r_index(motif_idx0)]
  data.frame(
    x$row_metadata,
    motif_idx = cf_index0_to_r_index(motif_idx0),
    motif = motif,
    count = as.vector(counts),
    stringsAsFactors = FALSE
  )
}

#' @export
#' @rdname sparse_data_frame
sparse_data_frame.cfdnalab_end_motif_counts <- function(x, ...) {
  cf_reject_unused_arguments(...)
  if (!identical(x$storage_mode, "sparse_coo")) {
    stop("sparse_data_frame() is only available for sparse_coo output", call. = FALSE)
  }
  row_idx0 <- as.integer(x$sparse$row_idx0)
  motif_idx0 <- as.integer(x$sparse$motif_idx0)
  data.frame(
    row_idx = cf_index0_to_r_index(row_idx0),
    motif_idx = cf_index0_to_r_index(motif_idx0),
    motif = x$motif[cf_index0_to_r_index(motif_idx0)],
    count = as.numeric(x$sparse$count),
    stringsAsFactors = FALSE
  )
}

#' @export
#' @rdname sparse_data_frame_for_window
#' @param window_idx One-based window index.
sparse_data_frame_for_window.cfdnalab_windowed_end_motif_counts <- function(x, window_idx, ...) {
  cf_reject_unused_arguments(...)
  window_idx <- cf_validate_r_index(window_idx, length(x$row_idx0), "window_idx")
  cf_sparse_data_frame_for_row_idx0(x, cf_r_index_to_index0(window_idx))
}

#' @export
#' @rdname sparse_data_frame_for_group
#' @param group Group name or one-based group index.
sparse_data_frame_for_group.cfdnalab_grouped_end_motif_counts <- function(x, group, ...) {
  cf_reject_unused_arguments(...)
  group_idx0 <- cf_resolve_end_motif_group_idx0(x, group)
  cf_sparse_data_frame_for_row_idx0(x, group_idx0)
}

#' @export
#' @rdname sparse_data_frame_for_motif
#' @param motif Motif label.
sparse_data_frame_for_motif.cfdnalab_end_motif_counts <- function(x, motif, ...) {
  cf_reject_unused_arguments(...)
  motif_idx0 <- cf_resolve_end_motif_idx0(x, motif)
  if (!identical(x$storage_mode, "sparse_coo")) {
    stop("sparse_data_frame_for_motif() is only available for sparse_coo output", call. = FALSE)
  }
  matches <- as.integer(x$sparse$motif_idx0) == motif_idx0
  row_idx0 <- as.integer(x$sparse$row_idx0[matches])
  data.frame(
    x$row_metadata[cf_index0_to_r_index(row_idx0), , drop = FALSE],
    motif_idx = rep(cf_index0_to_r_index(motif_idx0), length(row_idx0)),
    motif = rep(motif, length(row_idx0)),
    count = as.numeric(x$sparse$count[matches]),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Return sparse non-zero count rows for one count row.
#'
#' @param x A cfDNAlab end-motif object.
#' @param row_idx0 Internal zero-based row index.
#'
#' @return A data frame with one row per stored non-zero count.
#' @noRd
cf_sparse_data_frame_for_row_idx0 <- function(x, row_idx0) {
  if (!identical(x$storage_mode, "sparse_coo")) {
    stop("Sparse row data frames are only available for sparse_coo output", call. = FALSE)
  }
  row_idx0 <- cf_validate_index0(row_idx0, length(x$row_idx0), "row_idx0")
  matches <- as.integer(x$sparse$row_idx0) == row_idx0
  motif_idx0 <- as.integer(x$sparse$motif_idx0[matches])
  data.frame(
    x$row_metadata[cf_index0_to_r_index(row_idx0), , drop = FALSE][rep(1L, sum(matches)), , drop = FALSE],
    motif_idx = cf_index0_to_r_index(motif_idx0),
    motif = x$motif[cf_index0_to_r_index(motif_idx0)],
    count = as.numeric(x$sparse$count[matches]),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Resolve a grouped end-motif group selector.
#'
#' @param x A `cfdnalab_grouped_end_motif_counts` object.
#' @param group Group name or one-based group index.
#'
#' @return A zero-based group index.
#' @noRd
cf_resolve_end_motif_group_idx0 <- function(x, group) {
  if (is.character(group)) {
    return(cf_r_index_to_index0(group_idx(x, group)))
  }
  cf_r_index_to_index0(cf_validate_r_index(group, length(x$row_idx0), "group"))
}

#' Resolve an end-motif selector.
#'
#' @param x A `cfdnalab_end_motif_counts` object.
#' @param motif Motif label.
#'
#' @return A zero-based motif index.
#' @noRd
cf_resolve_end_motif_idx0 <- function(x, motif) {
  cf_r_index_to_index0(motif_idx(x, motif))
}

#' Print an end-motif object.
#'
#' @param x A cfDNAlab end-motif object.
#' @param ... Ignored.
#'
#' @return Invisibly returns `x`.
#' @export
#' @keywords internal
print.cfdnalab_end_motif_counts <- function(x, ...) {
  cat("<cfDNAlab end-motif counts>\n")
  cat("Path: ", x$path, "\n", sep = "")
  cat("Storage mode: ", x$storage_mode, "\n", sep = "")
  cat("Row mode: ", x$row_mode, "\n", sep = "")
  cat("Rows: ", length(x$row_idx0), "\n", sep = "")
  cat("Motifs: ", length(x$motif), "\n", sep = "")
  invisible(x)
}
