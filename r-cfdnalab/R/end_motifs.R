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

#' Supported end-motif column-axis kinds.
#'
#' @noRd
END_MOTIF_VALID_AXIS_KINDS <- c("motif", "motif_group")

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
  motif_axis_kind <- root_attributes$motif_axis_kind
  if (is.null(motif_axis_kind)) {
    if (!identical(as.integer(root_attributes$cfdnalab_schema_version), 1L)) {
      stop("end-motif schema v2 stores must declare motif_axis_kind", call. = FALSE)
    }
    motif_axis_kind <- "motif"
  }
  motif_axis_kind <- cf_validate_allowed_string(
    motif_axis_kind,
    END_MOTIF_VALID_AXIS_KINDS,
    "end-motif motif axis kind"
  )

  cf_reject_empty_end_motif_counts(path, storage_mode)

  store <- cf_open_zarr(path, "end-motif")
  cf_required_arrays(
    store,
    cf_end_motif_required_arrays(storage_mode, row_mode, motif_axis_kind),
    "End-motif"
  )
  cf_validate_dimension_names(path, "motif_index", "motif")
  cf_validate_dimension_names(path, "row", "row")

  motif_axis <- cf_read_vector(store, "motif_index", "End-motif")
  row <- cf_read_vector(store, "row", "End-motif")
  motif <- NULL
  if (identical(motif_axis_kind, "motif")) {
    cf_validate_dimension_names(path, "motif_byte", "motif_byte")
    cf_validate_dimension_names(path, "motif_ascii", c("motif", "motif_byte"))
    motif_byte <- cf_read_vector(store, "motif_byte", "End-motif")
    motif_ascii <- cf_read_array(store, "motif_ascii", "End-motif")
    motif <- cf_decode_motif_ascii(motif_ascii, length(motif_axis), length(motif_byte))
    cf_validate_axis(motif_byte, "motif_byte")
  } else {
    motif <- cf_read_labels(path, "motif_index", "motif_group", length(motif_axis))
  }

  cf_validate_axis(motif_axis, "motif_index")
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
    motif_axis_kind = motif_axis_kind,
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
#' @param motif_axis_kind End-motif column-axis kind.
#'
#' @return Character vector of required array paths.
#' @noRd
cf_end_motif_required_arrays <- function(storage_mode, row_mode, motif_axis_kind) {
  required <- c("motif_index", "row")
  if (identical(motif_axis_kind, "motif")) {
    required <- c(required, "motif_byte", "motif_ascii")
  }
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

#' Reject sparse stores without count entries before opening the Zarr hierarchy.
#'
#' The CRAN `zarr` reader rejects zero-length array dimensions while opening a
#' store. Sparse end-motif stores with no observed counts therefore need a
#' schema-level error before `zarr::open_zarr()` constructs every array node.
#'
#' @param path Path to an end-motif Zarr store.
#' @param storage_mode End-motif storage mode.
#'
#' @return Invisibly returns `TRUE` when count arrays are non-empty.
#' @noRd
cf_reject_empty_end_motif_counts <- function(path, storage_mode) {
  motif_shape <- cf_zarr_array_metadata_shape(path, "motif_index")
  if (!is.null(motif_shape) && length(motif_shape) >= 1L && motif_shape[[1L]] == 0L) {
    cf_stop_no_end_motif_counts_available()
  }

  if (identical(storage_mode, "sparse_coo")) {
    sparse_count_shape <- cf_zarr_array_metadata_shape(path, "sparse/count")
    if (
      !is.null(sparse_count_shape) &&
        length(sparse_count_shape) >= 1L &&
        sparse_count_shape[[1L]] == 0L
    ) {
      cf_stop_no_end_motif_counts_available()
    }
  }

  invisible(TRUE)
}

#' Return the raw metadata shape for one Zarr array.
#'
#' Missing arrays are left to the normal required-array validation. This helper
#' only exists to catch valid-but-empty sparse stores before `zarr::open_zarr()`
#' raises a generic shape error.
#'
#' @param path Path to a Zarr store.
#' @param array_name Slash-separated array path within the store.
#'
#' @return Integer shape vector, or `NULL` when the array metadata is absent.
#' @noRd
cf_zarr_array_metadata_shape <- function(path, array_name) {
  metadata_path <- do.call(
    file.path,
    as.list(c(path, strsplit(array_name, "/", fixed = TRUE)[[1L]], "zarr.json"))
  )
  if (!file.exists(metadata_path)) {
    return(NULL)
  }

  metadata <- cf_read_json_file(metadata_path)
  shape <- unlist(metadata$shape, use.names = FALSE)
  if (
    !is.numeric(shape) ||
      any(is.na(shape)) ||
      any(!is.finite(shape)) ||
      any(shape != floor(shape)) ||
      any(shape > .Machine$integer.max) ||
      any(shape < 0)
  ) {
    stop(array_name, " metadata shape must be a non-negative integer vector", call. = FALSE)
  }
  as.integer(shape)
}

#' Raise the public no-counts error for end-motif stores.
#'
#' @return Never returns.
#' @noRd
cf_stop_no_end_motif_counts_available <- function() {
  stop(
    "No end-motif counts are available in this store. ",
    "If you expected motifs or groups from `--motifs-file` with zero counts to remain in the output, ",
    "rerun `cfdna ends` with `--all-motifs`.",
    call. = FALSE
  )
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
    if (any(row_start_bp > .Machine$integer.max) || any(row_end_bp > .Machine$integer.max)) {
      stop("End-motif window coordinates must fit in R integer range", call. = FALSE)
    }
    return(data.frame(
      window_idx = cf_index0_to_r_index(row),
      chrom = chromosome_name[as.integer(row_chromosome) + 1L],
      start = as.integer(row_start_bp),
      end = as.integer(row_end_bp),
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
  if (any(eligible_windows > .Machine$integer.max)) {
    stop("eligible_windows values must fit in R integer range", call. = FALSE)
  }
  data.frame(
    group_idx = cf_index0_to_r_index(group),
    group_name = group_name,
    eligible_windows = as.integer(eligible_windows),
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
  matched_index <- cf_find_unique_value_index(
    x$motif,
    motif,
    "Unknown end-motif label: ",
    "End-motif label is not unique: "
  )
  cf_index0_to_r_index(x$motif_idx0[[matched_index]])
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
#' @rdname window_metadata
window_metadata.cfdnalab_windowed_end_motif_counts <- function(x, ...) {
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
  matched_index <- cf_find_unique_value_index(
    x$row_metadata$group_name,
    group_name,
    "Unknown end-motif group name: ",
    "End-motif group name is not unique: "
  )
  x$row_metadata$group_idx[[matched_index]]
}

#' @export
#' @rdname sparse_counts_matrix
#' @param motifs Optional motif label vector. Use either `motifs` or
#'   `motif_idxs`, not both.
#' @param motif_idxs Optional one-based motif index vector.
sparse_counts_matrix.cfdnalab_global_end_motif_counts <- function(
  x,
  motifs = NULL,
  motif_idxs = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_sparse_end_motif_matrix_for_indices(
    x,
    row_indices = seq_len(length(x$row_idx0)),
    motif_indices = cf_resolve_end_motif_indices(x, motifs, motif_idxs)
  )
}

#' @export
#' @rdname sparse_counts_matrix
#' @param window_idxs Optional one-based window index vector for windowed output.
sparse_counts_matrix.cfdnalab_windowed_end_motif_counts <- function(
  x,
  window_idxs = NULL,
  motifs = NULL,
  motif_idxs = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_sparse_end_motif_matrix_for_indices(
    x,
    row_indices = cf_resolve_end_motif_window_indices(x, window_idxs),
    motif_indices = cf_resolve_end_motif_indices(x, motifs, motif_idxs)
  )
}

#' @export
#' @rdname sparse_counts_matrix
#' @param groups Optional group name vector for grouped output. Use either
#'   `groups` or `group_idxs`, not both.
#' @param group_idxs Optional one-based group index vector for grouped output.
sparse_counts_matrix.cfdnalab_grouped_end_motif_counts <- function(
  x,
  groups = NULL,
  group_idxs = NULL,
  motifs = NULL,
  motif_idxs = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_sparse_end_motif_matrix_for_indices(
    x,
    row_indices = cf_resolve_end_motif_group_indices(x, groups, group_idxs),
    motif_indices = cf_resolve_end_motif_indices(x, motifs, motif_idxs)
  )
}

#' Build a sparse end-motif matrix for selected rows and motifs.
#'
#' @param x End-motif object.
#' @param row_indices One-based row indices.
#' @param motif_indices One-based motif indices.
#'
#' @return A `Matrix` sparse matrix.
#' @noRd
cf_sparse_end_motif_matrix_for_indices <- function(x, row_indices, motif_indices) {
  if (identical(x$storage_mode, "dense")) {
    return(Matrix::Matrix(
      x$counts$read()[row_indices, motif_indices, drop = FALSE],
      sparse = TRUE
    ))
  }
  selected_row_idx0 <- cf_r_index_to_index0(row_indices)
  selected_motif_idx0 <- cf_r_index_to_index0(motif_indices)
  sparse_row_idx0 <- as.integer(x$sparse$row_idx0)
  sparse_motif_idx0 <- as.integer(x$sparse$motif_idx0)
  matches <- sparse_row_idx0 %in% selected_row_idx0 &
    sparse_motif_idx0 %in% selected_motif_idx0
  Matrix::sparseMatrix(
    i = match(sparse_row_idx0[matches], selected_row_idx0),
    j = match(sparse_motif_idx0[matches], selected_motif_idx0),
    x = as.numeric(x$sparse$count[matches]),
    dims = as.integer(c(length(row_indices), length(motif_indices)))
  )
}

#' @export
#' @rdname dense_counts_matrix
#' @param motifs Optional motif label vector. Use either `motifs` or
#'   `motif_idxs`, not both.
#' @param motif_idxs Optional one-based motif index vector.
#' @param allow_densify If `TRUE`, allow sparse stores to be converted to a dense
#'   in-memory matrix. Sparse stores error by default.
dense_counts_matrix.cfdnalab_global_end_motif_counts <- function(
  x,
  allow_densify = FALSE,
  motifs = NULL,
  motif_idxs = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_dense_end_motif_matrix_for_indices(
    x,
    row_indices = seq_len(length(x$row_idx0)),
    motif_indices = cf_resolve_end_motif_indices(x, motifs, motif_idxs),
    allow_densify = allow_densify
  )
}

#' @export
#' @rdname dense_counts_matrix
#' @param window_idxs Optional one-based window index vector for windowed output.
dense_counts_matrix.cfdnalab_windowed_end_motif_counts <- function(
  x,
  allow_densify = FALSE,
  window_idxs = NULL,
  motifs = NULL,
  motif_idxs = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_dense_end_motif_matrix_for_indices(
    x,
    row_indices = cf_resolve_end_motif_window_indices(x, window_idxs),
    motif_indices = cf_resolve_end_motif_indices(x, motifs, motif_idxs),
    allow_densify = allow_densify
  )
}

#' @export
#' @rdname dense_counts_matrix
#' @param groups Optional group name vector for grouped output. Use either
#'   `groups` or `group_idxs`, not both.
#' @param group_idxs Optional one-based group index vector for grouped output.
dense_counts_matrix.cfdnalab_grouped_end_motif_counts <- function(
  x,
  allow_densify = FALSE,
  groups = NULL,
  group_idxs = NULL,
  motifs = NULL,
  motif_idxs = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_dense_end_motif_matrix_for_indices(
    x,
    row_indices = cf_resolve_end_motif_group_indices(x, groups, group_idxs),
    motif_indices = cf_resolve_end_motif_indices(x, motifs, motif_idxs),
    allow_densify = allow_densify
  )
}

#' Build a dense end-motif matrix for selected rows and motifs.
#'
#' @param x End-motif object.
#' @param row_indices One-based row indices.
#' @param motif_indices One-based motif indices.
#' @param allow_densify Whether to allow sparse-store densification.
#'
#' @return A dense numeric matrix.
#' @noRd
cf_dense_end_motif_matrix_for_indices <- function(
  x,
  row_indices,
  motif_indices,
  allow_densify
) {
  if (identical(x$storage_mode, "dense")) {
    return(cf_read_array(x$store, "counts", "End-motif")[
      row_indices,
      motif_indices,
      drop = FALSE
    ])
  }
  if (!isTRUE(allow_densify)) {
    stop(
      "This end-motif store is sparse. Use sparse_counts_matrix() or set allow_densify = TRUE.",
      call. = FALSE
    )
  }
  as.matrix(cf_sparse_end_motif_matrix_for_indices(x, row_indices, motif_indices))
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
#' @rdname end_motif_data_frame
#' @param densify If `TRUE`, sparse outputs add explicit zero-count rows for
#'   selected observed motifs. Dense outputs ignore this option.
#' @param motifs Optional motif label vector. Use either `motifs` or
#'   `motif_idxs`, not both.
#' @param motif_idxs Optional one-based motif index vector.
end_motif_data_frame.cfdnalab_global_end_motif_counts <- function(
  x,
  densify = FALSE,
  motifs = NULL,
  motif_idxs = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_end_motif_data_frame(
    x,
    row_indices = seq_len(length(x$row_idx0)),
    motif_indices = cf_resolve_end_motif_indices(x, motifs, motif_idxs),
    densify = densify,
    max_blacklisted_fraction = 1.0
  )
}

#' @export
#' @rdname end_motif_data_frame
#' @param window_idxs Optional one-based window index vector for windowed output.
#' @param max_blacklisted_fraction Maximum row `blacklisted_fraction` in 0..1
#'   to retain before returning counts. The default `1.0` keeps all selected
#'   rows.
end_motif_data_frame.cfdnalab_windowed_end_motif_counts <- function(
  x,
  window_idxs = NULL,
  densify = FALSE,
  motifs = NULL,
  motif_idxs = NULL,
  max_blacklisted_fraction = 1.0,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_end_motif_data_frame(
    x,
    row_indices = cf_resolve_end_motif_window_indices(x, window_idxs),
    motif_indices = cf_resolve_end_motif_indices(x, motifs, motif_idxs),
    densify = densify,
    max_blacklisted_fraction = max_blacklisted_fraction
  )
}

#' @export
#' @rdname end_motif_data_frame
#' @param groups Optional group name vector for grouped output. Use either
#'   `groups` or `group_idxs`, not both.
#' @param group_idxs Optional one-based group index vector for grouped output.
end_motif_data_frame.cfdnalab_grouped_end_motif_counts <- function(
  x,
  groups = NULL,
  group_idxs = NULL,
  densify = FALSE,
  motifs = NULL,
  motif_idxs = NULL,
  max_blacklisted_fraction = 1.0,
  ...
) {
  cf_reject_unused_arguments(...)
  row_indices <- cf_resolve_end_motif_group_indices(x, groups, group_idxs)
  cf_end_motif_data_frame(
    x,
    row_indices = row_indices,
    motif_indices = cf_resolve_end_motif_indices(x, motifs, motif_idxs),
    densify = densify,
    max_blacklisted_fraction = max_blacklisted_fraction
  )
}

#' Shared implementation for mode-specific end-motif data frame methods.
#'
#' @param x End-motif object.
#' @param row_indices One-based row indices.
#' @param motif_indices One-based motif indices.
#' @param densify Whether to densify sparse stores.
#' @param max_blacklisted_fraction Maximum blacklist fraction.
#'
#' @return A data frame.
#' @noRd
cf_end_motif_data_frame <- function(
  x,
  row_indices,
  motif_indices,
  densify,
  max_blacklisted_fraction
) {
  cf_validate_scalar_logical(densify, "densify")
  row_indices <- cf_apply_end_motif_blacklist_filter(x, row_indices, max_blacklisted_fraction)
  if (identical(x$storage_mode, "sparse_coo") && !isTRUE(densify)) {
    return(cf_stored_end_motif_data_frame_for_indices(x, row_indices, motif_indices))
  }
  cf_complete_end_motif_data_frame_for_indices(x, row_indices, motif_indices, densify)
}

#' Apply a blacklist fraction filter to end-motif row indices.
#'
#' @param x End-motif object.
#' @param row_indices One-based row indices.
#' @param max_blacklisted_fraction Maximum blacklist fraction.
#'
#' @return Filtered one-based row indices.
#' @noRd
cf_apply_end_motif_blacklist_filter <- function(x, row_indices, max_blacklisted_fraction) {
  cf_apply_row_blacklist_filter(x$row_metadata, row_indices, max_blacklisted_fraction)
}

#' Resolve grouped end-motif selectors to one-based row indices.
#'
#' @param x Grouped end-motif object.
#' @param groups Optional group names.
#' @param group_idxs Optional one-based group indices.
#'
#' @return One-based row indices.
#' @noRd
cf_resolve_end_motif_group_indices <- function(x, groups, group_idxs) {
  if (!is.null(groups) && !is.null(group_idxs)) {
    stop("Use either groups or group_idxs, not both", call. = FALSE)
  }
  if (!is.null(groups)) {
    cf_validate_character_vector(groups, "groups")
    cf_validate_unique_values(groups, "groups")
    return(vapply(
      groups,
      function(group_name) {
        cf_find_unique_value_index(
          x$row_metadata$group_name,
          group_name,
          "Unknown end-motif group name: ",
          "End-motif group name is not unique: "
        )
      },
      integer(1L),
      USE.NAMES = FALSE
    ))
  }
  if (!is.null(group_idxs)) {
    group_indices <- cf_validate_r_indices(
      group_idxs,
      length(x$row_idx0),
      "group_idxs"
    )
    cf_validate_unique_values(group_indices, "group_idxs")
    return(group_indices)
  }
  seq_len(length(x$row_idx0))
}

#' Resolve windowed end-motif selectors to one-based row indices.
#'
#' @param x Windowed end-motif object.
#' @param window_idxs Optional one-based window indices.
#'
#' @return One-based row indices.
#' @noRd
cf_resolve_end_motif_window_indices <- function(x, window_idxs) {
  if (is.null(window_idxs)) {
    return(seq_len(length(x$row_idx0)))
  }
  window_indices <- cf_validate_r_indices(
    window_idxs,
    length(x$row_idx0),
    "window_idxs"
  )
  cf_validate_unique_values(window_indices, "window_idxs")
  window_indices
}

#' Resolve end-motif selectors to one-based motif indices.
#'
#' @param x End-motif object.
#' @param motifs Optional motif labels.
#' @param motif_idxs Optional one-based motif indices.
#'
#' @return One-based motif indices.
#' @noRd
cf_resolve_end_motif_indices <- function(x, motifs, motif_idxs) {
  if (!is.null(motifs) && !is.null(motif_idxs)) {
    stop("Use either motifs or motif_idxs, not both", call. = FALSE)
  }
  if (!is.null(motifs)) {
    cf_validate_character_vector(motifs, "motifs")
    cf_validate_unique_values(motifs, "motifs")
    return(vapply(
      motifs,
      function(motif) {
        motif_idx(x, motif)
      },
      integer(1L),
      USE.NAMES = FALSE
    ))
  }
  if (!is.null(motif_idxs)) {
    motif_indices <- cf_validate_r_indices(
      motif_idxs,
      length(x$motif_idx0),
      "motif_idxs"
    )
    cf_validate_unique_values(motif_indices, "motif_idxs")
    return(motif_indices)
  }
  seq_along(x$motif_idx0)
}

#' Build a complete end-motif data frame for selected rows and motifs.
#'
#' @param x End-motif object.
#' @param row_indices One-based row indices.
#' @param motif_indices One-based motif indices.
#' @param densify Whether to allow sparse-store densification.
#'
#' @return A data frame with one row per selected row and motif.
#' @noRd
cf_complete_end_motif_data_frame_for_indices <- function(x, row_indices, motif_indices, densify) {
  if (length(row_indices) == 0L || length(motif_indices) == 0L) {
    return(cf_empty_end_motif_data_frame(x))
  }
  counts <- dense_counts_matrix(x, allow_densify = densify)[
    row_indices,
    motif_indices,
    drop = FALSE
  ]
  num_rows <- length(row_indices)
  num_motifs <- length(motif_indices)
  metadata <- x$row_metadata[row_indices, , drop = FALSE]
  metadata <- metadata[rep(seq_len(num_rows), each = num_motifs), , drop = FALSE]
  motif_metadata <- motifs(x)[motif_indices, , drop = FALSE]
  motif_metadata <- motif_metadata[rep(seq_len(num_motifs), times = num_rows), , drop = FALSE]
  data.frame(
    metadata,
    motif_metadata,
    count = as.vector(t(counts)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Build an end-motif data frame from stored COO rows.
#'
#' @param x End-motif object.
#' @param row_indices One-based row indices.
#' @param motif_indices One-based motif indices.
#'
#' @return A data frame with one row per stored non-zero count.
#' @noRd
cf_stored_end_motif_data_frame_for_indices <- function(x, row_indices, motif_indices) {
  if (length(row_indices) == 0L || length(motif_indices) == 0L) {
    return(cf_empty_end_motif_data_frame(x))
  }
  selected_row_idx0 <- cf_r_index_to_index0(row_indices)
  selected_motif_idx0 <- cf_r_index_to_index0(motif_indices)
  sparse_row_idx0 <- as.integer(x$sparse$row_idx0)
  sparse_motif_idx0 <- as.integer(x$sparse$motif_idx0)
  matches <- sparse_row_idx0 %in% selected_row_idx0 &
    sparse_motif_idx0 %in% selected_motif_idx0
  if (!any(matches)) {
    return(cf_empty_end_motif_data_frame(x))
  }
  matched_row_idx0 <- sparse_row_idx0[matches]
  matched_motif_idx0 <- sparse_motif_idx0[matches]
  sort_order <- order(
    match(matched_row_idx0, selected_row_idx0),
    match(matched_motif_idx0, selected_motif_idx0)
  )
  matched_row_idx0 <- matched_row_idx0[sort_order]
  matched_motif_idx0 <- matched_motif_idx0[sort_order]
  matched_counts <- as.numeric(x$sparse$count[matches])[sort_order]
  matched_row_indices <- cf_index0_to_r_index(matched_row_idx0)
  matched_motif_indices <- cf_index0_to_r_index(matched_motif_idx0)
  data.frame(
    x$row_metadata[matched_row_indices, , drop = FALSE],
    motifs(x)[matched_motif_indices, , drop = FALSE],
    count = matched_counts,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Build an empty end-motif data frame with public columns.
#'
#' @param x End-motif object.
#'
#' @return A zero-row data frame.
#' @noRd
cf_empty_end_motif_data_frame <- function(x) {
  data.frame(
    x$row_metadata[integer(0), , drop = FALSE],
    motifs(x)[integer(0), , drop = FALSE],
    count = numeric(),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
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
