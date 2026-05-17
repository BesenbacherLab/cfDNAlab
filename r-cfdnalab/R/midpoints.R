#' Required arrays for midpoint-profile Zarr stores.
#'
#' The loader checks these before reading data so missing-array errors are
#' reported together instead of failing later during a specific metadata read.
#'
#' @noRd
MIDPOINT_REQUIRED_ARRAYS <- c(
  "counts",
  "group",
  "eligible_intervals",
  "length_bin",
  "length_start_bp",
  "length_end_bp",
  "position",
  "position_bin_start_bp",
  "position_bin_end_bp"
)

#' Read cfDNAlab midpoint profiles.
#'
#' Loads a \code{<prefix>.midpoint_profiles.zarr} store created with the
#' \code{cfdna midpoints} CLI tool from the main \code{cfDNAlab} rust package.
#' It validates the schema, coordinate axes, labels, and count-array shape.
#'
#' @param path Path to a cfDNAlab midpoint-profile `.zarr` directory.
#'
#' @return A `cfdnalab_midpoint_profiles` object.
#' @export
#'
#' @examples
#' \dontrun{
#' midpoints <- read_midpoints("sample.midpoint_profiles.zarr")
#' groups(midpoints)
#' profile_data_frame(midpoints, group = "LYL1", length_bin_idx = 1)
#' }
read_midpoints <- function(path) {
  path <- cf_validate_zarr_store_path(path, "Midpoint profile")
  root_attributes <- cf_root_attributes(path)
  cf_validate_schema(root_attributes, "midpoint_profiles", "midpoint profile")

  store <- cf_open_zarr(path, "midpoint profile")
  cf_required_arrays(store, MIDPOINT_REQUIRED_ARRAYS, "Midpoint profile")
  cf_validate_dimension_names(path, "counts", c("group", "length_bin", "position"))
  cf_validate_dimension_names(path, "group", "group")
  cf_validate_dimension_names(path, "eligible_intervals", "group")
  cf_validate_dimension_names(path, "length_bin", "length_bin")
  cf_validate_dimension_names(path, "length_start_bp", "length_bin")
  cf_validate_dimension_names(path, "length_end_bp", "length_bin")
  cf_validate_dimension_names(path, "position", "position")
  cf_validate_dimension_names(path, "position_bin_start_bp", "position")
  cf_validate_dimension_names(path, "position_bin_end_bp", "position")

  counts <- cf_get_array(store, "counts", "Midpoint profile")
  group <- cf_read_vector(store, "group", "Midpoint profile")
  length_bin <- cf_read_vector(store, "length_bin", "Midpoint profile")
  position <- cf_read_vector(store, "position", "Midpoint profile")
  group_name <- cf_read_labels(path, "group", "group_name", length(group))
  eligible_intervals <- cf_read_vector(store, "eligible_intervals", "Midpoint profile")
  length_start_bp <- cf_read_vector(store, "length_start_bp", "Midpoint profile")
  length_end_bp <- cf_read_vector(store, "length_end_bp", "Midpoint profile")
  position_bin_start_bp <- cf_read_vector(store, "position_bin_start_bp", "Midpoint profile")
  position_bin_end_bp <- cf_read_vector(store, "position_bin_end_bp", "Midpoint profile")

  cf_validate_axis(group, "group")
  cf_validate_axis(length_bin, "length_bin")
  cf_validate_axis(position, "position")
  cf_validate_same_length(eligible_intervals, group, "eligible_intervals", "group")
  cf_validate_same_length(length_start_bp, length_bin, "length_start_bp", "length_bin")
  cf_validate_same_length(length_end_bp, length_bin, "length_end_bp", "length_bin")
  cf_validate_same_length(position_bin_start_bp, position, "position_bin_start_bp", "position")
  cf_validate_same_length(position_bin_end_bp, position, "position_bin_end_bp", "position")
  cf_validate_nonnegative_integer_vector(eligible_intervals, "eligible_intervals")
  cf_validate_half_open_intervals(
    length_start_bp,
    length_end_bp,
    "length_start_bp",
    "length_end_bp"
  )
  cf_validate_half_open_intervals(
    position_bin_start_bp,
    position_bin_end_bp,
    "position_bin_start_bp",
    "position_bin_end_bp"
  )

  expected_shape <- c(length(group), length(length_bin), length(position))
  if (!identical(cf_array_shape(counts), expected_shape)) {
    stop(
      "counts shape does not match midpoint coordinate arrays: counts=",
      paste(cf_array_shape(counts), collapse = "x"),
      ", coordinates=",
      paste(expected_shape, collapse = "x"),
      call. = FALSE
    )
  }

  structure(
    list(
      path = path,
      store = store,
      root_attributes = root_attributes,
      counts = counts,
      group_idx0 = as.integer(group),
      group_name = group_name,
      eligible_intervals = eligible_intervals,
      length_bin_idx0 = as.integer(length_bin),
      length_start_bp = length_start_bp,
      length_end_bp = length_end_bp,
      position_idx0 = as.integer(position),
      position_bin_start_bp = position_bin_start_bp,
      position_bin_end_bp = position_bin_end_bp
    ),
    class = c("cfdnalab_midpoint_profiles", "cfdnalab_zarr_store")
  )
}

#' @export
#' @rdname groups
groups.cfdnalab_midpoint_profiles <- function(x, ...) {
  cf_reject_unused_arguments(...)
  data.frame(
    group_idx = cf_index0_to_r_index(x$group_idx0),
    group_name = x$group_name,
    eligible_intervals = x$eligible_intervals,
    stringsAsFactors = FALSE
  )
}

#' @export
#' @rdname length_bins
length_bins.cfdnalab_midpoint_profiles <- function(x, ...) {
  cf_reject_unused_arguments(...)
  data.frame(
    length_bin_idx = cf_index0_to_r_index(x$length_bin_idx0),
    length_start_bp = x$length_start_bp,
    length_end_bp = x$length_end_bp,
    stringsAsFactors = FALSE
  )
}

#' @export
#' @rdname positions
positions.cfdnalab_midpoint_profiles <- function(x, ...) {
  cf_reject_unused_arguments(...)
  data.frame(
    position_idx = cf_index0_to_r_index(x$position_idx0),
    position_bin_start_bp = x$position_bin_start_bp,
    position_bin_end_bp = x$position_bin_end_bp,
    stringsAsFactors = FALSE
  )
}

#' @export
#' @rdname group_idx
#' @param group_name Group name to look up.
group_idx.cfdnalab_midpoint_profiles <- function(x, group_name, ...) {
  cf_reject_unused_arguments(...)
  group_name <- cf_validate_scalar_string(group_name, "group_name")
  matches <- which(x$group_name == group_name)
  if (length(matches) == 0L) {
    stop("Unknown midpoint group name: ", sQuote(group_name), call. = FALSE)
  }
  if (length(matches) > 1L) {
    stop("Midpoint group name is not unique: ", sQuote(group_name), call. = FALSE)
  }
  cf_index0_to_r_index(x$group_idx0[[matches]])
}

#' @export
#' @rdname length_bin_idx
#' @param length Fragment length in base pairs.
length_bin_idx.cfdnalab_midpoint_profiles <- function(x, length, ...) {
  cf_reject_unused_arguments(...)
  if (
    length(length) != 1L ||
      !is.numeric(length) ||
      is.na(length) ||
      !is.finite(length) ||
      length != as.integer(length) ||
      length < 0L
  ) {
    stop("Fragment length must be a single non-negative integer", call. = FALSE)
  }
  matches <- which(x$length_start_bp <= length & length < x$length_end_bp)
  if (length(matches) == 0L) {
    stop("No midpoint length bin contains length ", length, call. = FALSE)
  }
  if (length(matches) > 1L) {
    stop("Multiple midpoint length bins contain length ", length, call. = FALSE)
  }
  cf_index0_to_r_index(x$length_bin_idx0[[matches]])
}

#' @export
#' @rdname profile_array
#' @param group_idx One-based group index. Use either `group_idx` or `group`.
#' @param length_bin_idx One-based length-bin index. Use either `length_bin_idx` or `length`.
#' @param group Group name. Use either `group_idx` or `group`.
#' @param length Fragment length in base pairs. Use either `length_bin_idx` or `length`.
profile_array.cfdnalab_midpoint_profiles <- function(
  x,
  group_idx = NULL,
  length_bin_idx = NULL,
  group = NULL,
  length = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  group_idx0 <- cf_resolve_midpoint_group_idx0(x, group_idx, group)
  length_bin_idx0 <- cf_resolve_midpoint_length_bin_idx0(x, length_bin_idx, length)
  as.vector(cf_read_slice(
    x$counts,
    list(
      cf_index0_to_r_index(group_idx0),
      cf_index0_to_r_index(length_bin_idx0),
      seq_along(x$position_idx0)
    ),
    "Midpoint profile"
  ))
}

#' @export
#' @rdname profile_data_frame
#' @param group_idx One-based group index. Use either `group_idx` or `group`.
#' @param length_bin_idx One-based length-bin index. Use either `length_bin_idx` or `length`.
#' @param group Group name. Use either `group_idx` or `group`.
#' @param length Fragment length in base pairs. Use either `length_bin_idx` or `length`.
profile_data_frame.cfdnalab_midpoint_profiles <- function(
  x,
  group_idx = NULL,
  length_bin_idx = NULL,
  group = NULL,
  length = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  group_idx0 <- cf_resolve_midpoint_group_idx0(x, group_idx, group)
  length_bin_idx0 <- cf_resolve_midpoint_length_bin_idx0(x, length_bin_idx, length)
  group_r_index <- cf_index0_to_r_index(group_idx0)
  length_bin_r_index <- cf_index0_to_r_index(length_bin_idx0)
  counts <- profile_array(
    x,
    group_idx = group_r_index,
    length_bin_idx = length_bin_r_index
  )
  data.frame(
    group_idx = cf_index0_to_r_index(x$group_idx0[[group_r_index]]),
    group_name = x$group_name[[group_r_index]],
    eligible_intervals = x$eligible_intervals[[group_r_index]],
    length_bin_idx = cf_index0_to_r_index(x$length_bin_idx0[[length_bin_r_index]]),
    length_start_bp = x$length_start_bp[[length_bin_r_index]],
    length_end_bp = x$length_end_bp[[length_bin_r_index]],
    position_idx = cf_index0_to_r_index(x$position_idx0),
    position_bin_start_bp = x$position_bin_start_bp,
    position_bin_end_bp = x$position_bin_end_bp,
    count = counts,
    stringsAsFactors = FALSE
  )
}

#' @export
#' @rdname midpoint_array
midpoint_array.cfdnalab_midpoint_profiles <- function(x, ...) {
  cf_reject_unused_arguments(...)
  cf_read_array(x$store, "counts", "Midpoint profile")
}

#' Resolve a midpoint group selector.
#'
#' @param x A `cfdnalab_midpoint_profiles` object.
#' @param group_idx_value Optional one-based group index.
#' @param group_value Optional group name.
#'
#' @return A zero-based group index.
#' @noRd
cf_resolve_midpoint_group_idx0 <- function(x, group_idx_value, group_value) {
  if (!is.null(group_idx_value) && !is.null(group_value)) {
    stop("Use either group_idx or group, not both", call. = FALSE)
  }
  if (!is.null(group_value)) {
    return(cf_r_index_to_index0(group_idx(x, group_value)))
  }
  if (is.null(group_idx_value)) {
    stop("group_idx or group is required", call. = FALSE)
  }
  cf_r_index_to_index0(
    cf_validate_r_index(group_idx_value, length(x$group_idx0), "group_idx")
  )
}

#' Resolve a midpoint length-bin selector.
#'
#' @param x A `cfdnalab_midpoint_profiles` object.
#' @param length_bin_idx_value Optional one-based length-bin index.
#' @param length_value Optional fragment length in base pairs.
#'
#' @return A zero-based length-bin index.
#' @noRd
cf_resolve_midpoint_length_bin_idx0 <- function(x, length_bin_idx_value, length_value) {
  if (!is.null(length_bin_idx_value) && !is.null(length_value)) {
    stop("Use either length_bin_idx or length, not both", call. = FALSE)
  }
  if (!is.null(length_value)) {
    return(cf_r_index_to_index0(length_bin_idx(x, length_value)))
  }
  if (is.null(length_bin_idx_value)) {
    stop("length_bin_idx or length is required", call. = FALSE)
  }
  cf_r_index_to_index0(
    cf_validate_r_index(
      length_bin_idx_value,
      length(x$length_bin_idx0),
      "length_bin_idx"
    )
  )
}

#' Print a midpoint-profile object.
#'
#' @param x A `cfdnalab_midpoint_profiles` object.
#' @param ... Ignored.
#'
#' @return Invisibly returns `x`.
#' @export
#' @keywords internal
print.cfdnalab_midpoint_profiles <- function(x, ...) {
  cat("<cfDNAlab midpoint profiles>\n")
  cat("Path: ", x$path, "\n", sep = "")
  cat("Groups: ", length(x$group_idx0), "\n", sep = "")
  cat("Length bins: ", length(x$length_bin_idx0), "\n", sep = "")
  cat("Positions: ", length(x$position_idx0), "\n", sep = "")
  invisible(x)
}
