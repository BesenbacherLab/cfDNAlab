#' Supported cfDNAlab schema-version ranges.
#'
#' Version ranges are keyed by schema name so one output schema can evolve
#' without making unrelated loaders accept versions they do not understand.
#'
#' @noRd
CFDNALAB_SCHEMA_VERSION_RANGES <- list(
  midpoint_profiles = c(min = 1L, max = 1L),
  end_motif_counts = c(min = 1L, max = 1L)
)

#' Validate a cfDNAlab Zarr store path.
#'
#' @param path Path supplied by the user.
#' @param label Human-readable store label used in error messages.
#'
#' @return A normalized path string.
#' @noRd
cf_validate_zarr_store_path <- function(path, label) {
  if (length(path) != 1L || !is.character(path) || is.na(path)) {
    stop(label, " path must be a single path string", call. = FALSE)
  }
  path <- normalizePath(path, mustWork = FALSE)
  if (!file.exists(path)) {
    stop(label, " Zarr store does not exist: ", path, call. = FALSE)
  }
  if (!dir.exists(path)) {
    stop(label, " Zarr store path exists but is not a directory: ", path, call. = FALSE)
  }
  if (!grepl("\\.zarr$", path, ignore.case = TRUE)) {
    stop(label, " Zarr store path must end in '.zarr': ", path, call. = FALSE)
  }
  metadata_path <- file.path(path, "zarr.json")
  if (!file.exists(metadata_path)) {
    stop(label, " Zarr store is missing root zarr.json: ", metadata_path, call. = FALSE)
  }
  path
}

#' Read a JSON file as a list.
#'
#' @param path Path to a JSON file.
#'
#' @return A list parsed from JSON.
#' @noRd
cf_read_json_file <- function(path) {
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

#' Read root Zarr attributes.
#'
#' @param path Path to a Zarr store.
#'
#' @return A list of root attributes.
#' @noRd
cf_root_attributes <- function(path) {
  metadata <- cf_read_json_file(file.path(path, "zarr.json"))
  attrs <- metadata$attributes
  if (is.null(attrs)) {
    stop("Zarr root metadata is missing attributes: ", path, call. = FALSE)
  }
  attrs
}

#' Read array metadata from `zarr.json`.
#'
#' @param path Path to a Zarr store.
#' @param array_name Slash-separated array path within the store.
#'
#' @return A list parsed from the array metadata JSON.
#' @noRd
cf_array_metadata <- function(path, array_name) {
  metadata_path <- do.call(
    file.path,
    as.list(c(path, strsplit(array_name, "/", fixed = TRUE)[[1L]], "zarr.json"))
  )
  if (!file.exists(metadata_path)) {
    stop("Zarr array metadata is missing: ", metadata_path, call. = FALSE)
  }
  cf_read_json_file(metadata_path)
}

#' Read array attributes.
#'
#' @param path Path to a Zarr store.
#' @param array_name Slash-separated array path within the store.
#'
#' @return A list of array attributes.
#' @noRd
cf_array_attributes <- function(path, array_name) {
  metadata <- cf_array_metadata(path, array_name)
  attrs <- metadata$attributes
  if (is.null(attrs)) {
    attrs <- list()
  }
  attrs
}

#' Validate Zarr v3 dimension names for an array.
#'
#' @param path Path to a Zarr store.
#' @param array_name Slash-separated array path within the store.
#' @param expected_dimension_names Expected dimension names in order.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_dimension_names <- function(path, array_name, expected_dimension_names) {
  metadata <- cf_array_metadata(path, array_name)
  dimension_names <- unlist(metadata$dimension_names, use.names = FALSE)
  if (!is.character(dimension_names)) {
    stop(array_name, " dimensions must be character strings", call. = FALSE)
  }
  if (!identical(dimension_names, expected_dimension_names)) {
    stop(
      array_name,
      " dimensions must be ",
      paste(expected_dimension_names, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Validate cfDNAlab schema root attributes.
#'
#' @param attrs Root Zarr attributes.
#' @param expected_schema Expected cfDNAlab schema name.
#' @param label Human-readable schema label used in error messages.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_schema <- function(attrs, expected_schema, label) {
  schema <- attrs$cfdnalab_schema
  if (!identical(schema, expected_schema)) {
    stop(
      "Expected cfdnalab_schema='", expected_schema, "', found ",
      deparse(schema),
      call. = FALSE
    )
  }

  version_range <- CFDNALAB_SCHEMA_VERSION_RANGES[[expected_schema]]
  if (is.null(version_range)) {
    stop(
      "No supported schema-version range is registered for ",
      sQuote(expected_schema),
      call. = FALSE
    )
  }

  version <- attrs$cfdnalab_schema_version
  if (
    length(version) != 1L ||
      !is.numeric(version) ||
      is.na(version) ||
      !is.finite(version) ||
      version != as.integer(version) ||
      version < version_range[["min"]] ||
      version > version_range[["max"]]
  ) {
    stop(
      "Unsupported ", label, " schema version: ",
      deparse(version),
      ". Supported range: ",
      version_range[["min"]],
      "..",
      version_range[["max"]],
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Open a Zarr store with a package-specific error.
#'
#' @param path Path to a Zarr store.
#' @param label Human-readable store label used in error messages.
#'
#' @return A `zarr` store object.
#' The CRAN `zarr` reader loads `qs2` lazily while opening the V3 stores written
#' by cfDNAlab and uses `bit64` for `int64` coordinate arrays. Import one symbol
#' from each package so the runtime dependencies are installed with this package
#' and visible to R package checks.
#' @importFrom bit64 as.integer64
#' @importFrom qs2 qs_deserialize
#' @noRd
cf_open_zarr <- function(path, label) {
  tryCatch(
    zarr::open_zarr(path, read_only = TRUE),
    error = function(error) {
      stop(
        "Could not open ",
        label,
        " Zarr store at ",
        path,
        ": ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )
}

#' Convert an array name to a Zarr absolute path.
#'
#' @param array_name Slash-separated array path within the store.
#'
#' @return A slash-prefixed Zarr path.
#' @noRd
cf_zarr_path <- function(array_name) {
  paste0("/", array_name)
}

#' Find a unique exact-match index.
#'
#' @param values Vector to search.
#' @param value Scalar value to find.
#' @param unknown_message Error prefix for missing values.
#' @param duplicate_message Error prefix for duplicate values.
#'
#' @return A scalar one-based index into `values`.
#' @noRd
cf_find_unique_value_index <- function(values, value, unknown_message, duplicate_message) {
  matches <- which(values == value)
  if (length(matches) == 0L) {
    stop(unknown_message, sQuote(value), call. = FALSE)
  }
  if (length(matches) > 1L) {
    stop(duplicate_message, sQuote(value), call. = FALSE)
  }
  matches[[1L]]
}

#' Get a Zarr array node.
#'
#' @param store Open Zarr store.
#' @param array_name Slash-separated array path within the store.
#' @param label Human-readable store label used in error messages.
#'
#' @return A Zarr array node.
#' @noRd
cf_get_array <- function(store, array_name, label) {
  node <- tryCatch(
    cf_get_zarr_node(store, array_name),
    error = function(error) {
      stop(
        label,
        " Zarr store is missing array '",
        array_name,
        "': ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )
  if (is.null(node)) {
    stop(label, " Zarr store is missing array '", array_name, "'", call. = FALSE)
  }
  node
}

#' Get a Zarr node by path.
#'
#' The CRAN `zarr` package exposes nested V3 arrays through `get_node()` after
#' their parent group metadata has been discovered. Direct `[[` lookup works for
#' top-level arrays but returns `NULL` for nested paths such as `sparse/row`.
#'
#' @param store Open Zarr store.
#' @param array_name Slash-separated array path within the store.
#'
#' @return A Zarr node or `NULL`.
#' @noRd
cf_get_zarr_node <- function(store, array_name) {
  node_path <- cf_zarr_path(array_name)
  if (!is.null(store$get_node) && is.function(store$get_node)) {
    return(store$get_node(node_path))
  }
  store[[node_path]]
}

#' Read a full Zarr array.
#'
#' @param store Open Zarr store.
#' @param array_name Slash-separated array path within the store.
#' @param label Human-readable store label used in error messages.
#'
#' @return An R vector, matrix, or array.
#' @noRd
cf_read_array <- function(store, array_name, label) {
  array <- cf_get_array(store, array_name, label)
  if (!is.null(array$read) && is.function(array$read)) {
    return(array$read())
  }
  array[]
}

#' Read a full Zarr array as a vector.
#'
#' @param store Open Zarr store.
#' @param array_name Slash-separated array path within the store.
#' @param label Human-readable store label used in error messages.
#'
#' @return A vector.
#' @noRd
cf_read_vector <- function(store, array_name, label) {
  values <- cf_read_array(store, array_name, label)
  dimensions <- dim(values)
  if (!is.null(dimensions) && length(dimensions) > 1L) {
    stop(array_name, " must be a vector array", call. = FALSE)
  }
  if (!is.null(dimensions)) {
    dim(values) <- NULL
  }
  values
}

#' Read a selected Zarr array slice.
#'
#' @param array Zarr array node.
#' @param selection List of one-based R index selections.
#' @param label Human-readable store label used in error messages.
#'
#' @return An R vector, matrix, or array.
#' @noRd
cf_read_slice <- function(array, selection, label) {
  if (!is.null(array$read) && is.function(array$read)) {
    return(array$read(selection = selection))
  }
  stop(label, " Zarr backend does not expose array$read()", call. = FALSE)
}

#' Ensure arrays are present in a Zarr store.
#'
#' @param store Open Zarr store.
#' @param array_names Array paths that must exist.
#' @param label Human-readable store label used in error messages.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_required_arrays <- function(store, array_names, label) {
  missing <- character()
  for (array_name in array_names) {
    node <- tryCatch(cf_get_zarr_node(store, array_name), error = function(error) NULL)
    if (is.null(node)) {
      missing <- c(missing, array_name)
    }
  }
  if (length(missing) > 0L) {
    stop(label, " Zarr store is missing arrays: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

#' Reject unused public method arguments.
#'
#' @param ... Arguments captured by an S3 method.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_reject_unused_arguments <- function(...) {
  arguments <- list(...)
  if (length(arguments) > 0L) {
    argument_names <- names(arguments)
    argument_names[is.na(argument_names) | !nzchar(argument_names)] <- "<unnamed>"
    stop(
      "Unused argument(s): ",
      paste(argument_names, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Return a Zarr array shape as integers.
#'
#' @param array Zarr array node.
#'
#' @return Integer vector of dimensions.
#' @noRd
cf_array_shape <- function(array) {
  as.integer(array$shape)
}

#' Validate a contiguous zero-based axis.
#'
#' @param values Axis values.
#' @param axis_name Human-readable axis name.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_axis <- function(values, axis_name) {
  cf_validate_integer_vector(values, axis_name)
  if (any(values > .Machine$integer.max)) {
    stop(axis_name, " axis values must fit in R integer range", call. = FALSE)
  }
  expected <- if (length(values) == 0L) {
    integer()
  } else {
    seq.int(0L, length(values) - 1L)
  }
  if (!identical(as.integer(values), expected)) {
    stop(axis_name, " axis must contain contiguous zero-based indices", call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate integer-like numeric vector values.
#'
#' @param values Values to validate.
#' @param value_name Human-readable value name.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_integer_vector <- function(values, value_name) {
  if (
    !is.numeric(values) ||
      any(is.na(values)) ||
      any(!is.finite(values)) ||
      any(values != floor(values))
  ) {
    stop(value_name, " must contain integer values", call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate non-negative integer-like numeric vector values.
#'
#' @param values Values to validate.
#' @param value_name Human-readable value name.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_nonnegative_integer_vector <- function(values, value_name) {
  cf_validate_integer_vector(values, value_name)
  if (any(values < 0L)) {
    stop(value_name, " must contain non-negative integer values", call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate zero-based sparse index values.
#'
#' @param values Index values to validate.
#' @param size Axis size.
#' @param value_name Human-readable value name.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_index_vector <- function(values, size, value_name) {
  cf_validate_nonnegative_integer_vector(values, value_name)
  if (any(values > .Machine$integer.max)) {
    stop(value_name, " values must fit in R integer range", call. = FALSE)
  }
  if (any(values >= size)) {
    stop(value_name, " contains an index outside 0..", size - 1L, call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate finite non-negative numeric values.
#'
#' @param values Values to validate.
#' @param value_name Human-readable value name.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_nonnegative_numeric_vector <- function(values, value_name) {
  if (
    !is.numeric(values) ||
      any(is.na(values)) ||
      any(!is.finite(values)) ||
      any(values < 0)
  ) {
    stop(value_name, " must contain finite non-negative values", call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate paired half-open intervals.
#'
#' @param starts Start coordinates.
#' @param ends End coordinates.
#' @param start_name Human-readable start-coordinate name.
#' @param end_name Human-readable end-coordinate name.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_half_open_intervals <- function(starts, ends, start_name, end_name) {
  cf_validate_nonnegative_integer_vector(starts, start_name)
  cf_validate_nonnegative_integer_vector(ends, end_name)
  cf_validate_same_length(starts, ends, start_name, end_name)
  if (any(starts >= ends)) {
    stop(start_name, " must be smaller than ", end_name, call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate finite fractions.
#'
#' @param values Fraction values.
#' @param value_name Human-readable value name.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_fraction_vector <- function(values, value_name) {
  if (
    !is.numeric(values) ||
      any(is.na(values)) ||
      any(!is.finite(values)) ||
      any(values < 0 | values > 1)
  ) {
    stop(value_name, " must contain finite fractions in 0..1", call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate a scalar fraction.
#'
#' @param value Fraction value.
#' @param name Human-readable value name.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_scalar_fraction <- function(value, name) {
  if (
    length(value) != 1L ||
      !is.numeric(value) ||
      is.na(value) ||
      !is.finite(value) ||
      value < 0 ||
      value > 1
  ) {
    stop(name, " must be a single finite fraction in 0..1", call. = FALSE)
  }
  invisible(TRUE)
}

#' Apply a blacklist fraction filter to row indices.
#'
#' @param row_metadata Data frame with row metadata.
#' @param row_indices One-based row indices.
#' @param max_blacklisted_fraction Maximum blacklist fraction.
#'
#' @return Filtered one-based row indices.
#' @noRd
cf_apply_row_blacklist_filter <- function(row_metadata, row_indices, max_blacklisted_fraction) {
  cf_validate_scalar_fraction(max_blacklisted_fraction, "max_blacklisted_fraction")
  if (!"blacklisted_fraction" %in% names(row_metadata)) {
    if (max_blacklisted_fraction == 1) {
      return(row_indices)
    }
    stop("Cannot filter by max_blacklisted_fraction because this output has no blacklisted_fraction column", call. = FALSE)
  }
  row_indices[row_metadata$blacklisted_fraction[row_indices] <= max_blacklisted_fraction]
}

#' Validate two vectors have the same length.
#'
#' @param values Values being checked.
#' @param reference Reference vector.
#' @param value_name Human-readable value name.
#' @param reference_name Human-readable reference name.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_same_length <- function(values, reference, value_name, reference_name) {
  if (length(values) != length(reference)) {
    stop(
      value_name,
      " length (",
      length(values),
      ") does not match ",
      reference_name,
      " length (",
      length(reference),
      ")",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Validate a scalar R index.
#'
#' @param index Index supplied by the user.
#' @param size Axis size.
#' @param name Human-readable index name.
#'
#' @return The validated integer index.
#' @noRd
cf_validate_r_index <- function(index, size, name) {
  if (
    length(index) != 1L ||
      !is.numeric(index) ||
      is.na(index) ||
      !is.finite(index) ||
      index != as.integer(index)
  ) {
    stop(name, " must be a single integer", call. = FALSE)
  }
  index <- as.integer(index)
  if (index < 1L || index > size) {
    stop(name, " ", index, " is outside 1..", size, call. = FALSE)
  }
  index
}

#' Validate an internal zero-based index.
#'
#' @param index0 Internal index read from a cfDNAlab Zarr array.
#' @param size Axis size.
#' @param name Human-readable index name.
#'
#' @return The validated integer index.
#' @noRd
cf_validate_index0 <- function(index0, size, name) {
  if (
    length(index0) != 1L ||
      !is.numeric(index0) ||
      is.na(index0) ||
      !is.finite(index0) ||
      index0 != as.integer(index0)
  ) {
    stop(name, " must be a single integer", call. = FALSE)
  }
  index0 <- as.integer(index0)
  if (index0 < 0L || index0 >= size) {
    stop(name, " ", index0, " is outside 0..", size - 1L, call. = FALSE)
  }
  index0
}

#' Convert a public R index to an internal zero-based index.
#'
#' @param index Public one-based R index.
#'
#' @return A zero-based integer index.
#' @noRd
cf_r_index_to_index0 <- function(index) {
  as.integer(index) - 1L
}

#' Convert an internal zero-based index to a public R index.
#'
#' @param index0 Internal zero-based index.
#'
#' @return A one-based integer index.
#' @noRd
cf_index0_to_r_index <- function(index0) {
  as.integer(index0) + 1L
}

#' Validate a scalar string.
#'
#' @param value Value supplied by the user or read from metadata.
#' @param name Human-readable value name.
#'
#' @return The validated string.
#' @noRd
cf_validate_scalar_string <- function(value, name) {
  if (length(value) != 1L || is.na(value) || !is.character(value)) {
    stop(name, " must be a single character string", call. = FALSE)
  }
  value
}

#' Validate a scalar string against an allowed set.
#'
#' @param value Value supplied by the user or read from metadata.
#' @param allowed Allowed character values.
#' @param name Human-readable value name.
#'
#' @return The validated string.
#' @noRd
cf_validate_allowed_string <- function(value, allowed, name) {
  value <- cf_validate_scalar_string(value, name)
  if (!value %in% allowed) {
    stop("Unsupported ", name, ": ", deparse(value), call. = FALSE)
  }
  value
}

#' Read labels from Zarr array attributes.
#'
#' @param path Path to a Zarr store.
#' @param array_name Slash-separated array path within the store.
#' @param expected_field Expected `label_field` attribute.
#' @param expected_length Expected number of labels.
#'
#' @return Character vector of labels.
#' @noRd
cf_read_labels <- function(path, array_name, expected_field, expected_length) {
  attrs <- cf_array_attributes(path, array_name)
  if (!identical(attrs$label_field, expected_field)) {
    stop(
      array_name,
      " metadata must declare label_field = ",
      expected_field,
      call. = FALSE
    )
  }
  if (is.null(attrs$labels)) {
    stop(array_name, " metadata is missing labels", call. = FALSE)
  }
  labels <- unlist(attrs$labels, use.names = FALSE)
  if (!is.character(labels)) {
    stop(array_name, " labels must be character strings", call. = FALSE)
  }
  if (length(labels) != expected_length) {
    stop(
      array_name,
      " labels length (",
      length(labels),
      ") does not match axis length (",
      expected_length,
      ")",
      call. = FALSE
    )
  }
  as.character(labels)
}

#' Decode fixed-width motif labels from an ASCII byte matrix.
#'
#' @param bytes Numeric matrix of ASCII bytes.
#' @param n_motifs Expected number of motif rows.
#' @param motif_width Expected number of bytes per motif.
#'
#' @return Character vector of motif labels.
#' @noRd
cf_decode_motif_ascii <- function(bytes, n_motifs, motif_width) {
  if (is.null(dim(bytes))) {
    stop("motif_ascii must be a rank-2 array", call. = FALSE)
  }
  if (
    !is.numeric(bytes) ||
      any(is.na(bytes)) ||
      any(!is.finite(bytes)) ||
      any(bytes != as.integer(bytes)) ||
      any(bytes < 0L | bytes > 127L)
  ) {
    stop("motif_ascii must contain ASCII byte values in 0..127", call. = FALSE)
  }
  if (nrow(bytes) != n_motifs) {
    stop(
      "motif_ascii row count (",
      nrow(bytes),
      ") does not match motif_index length (",
      n_motifs,
      ")",
      call. = FALSE
    )
  }
  if (ncol(bytes) != motif_width) {
    stop(
      "motif_ascii column count (",
      ncol(bytes),
      ") does not match motif_byte length (",
      motif_width,
      ")",
      call. = FALSE
    )
  }
  unname(apply(bytes, 1L, function(row) {
    row <- as.integer(row)
    row <- row[row != 0L]
    rawToChar(as.raw(row))
  }))
}

#' @export
#' @rdname schema_version
schema_version.cfdnalab_zarr_store <- function(x, ...) {
  cf_reject_unused_arguments(...)
  as.integer(x$root_attributes$cfdnalab_schema_version)
}
