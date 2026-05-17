#' Return the cfDNAlab schema version.
#'
#' @param x A cfDNAlab loader object.
#' @param ... Reserved for future methods.
#'
#' @return An integer schema version.
#' @export
schema_version <- function(x, ...) {
  UseMethod("schema_version")
}

#' Return the end-motif storage mode.
#'
#' @param x A cfDNAlab end-motif object.
#' @param ... Reserved for future methods.
#'
#' @return A scalar character value, currently `"dense"` or `"sparse_coo"`.
#' @export
storage_mode <- function(x, ...) {
  UseMethod("storage_mode")
}

#' Return the end-motif row mode.
#'
#' @param x A cfDNAlab end-motif object.
#' @param ... Reserved for future methods.
#'
#' @return A scalar character value describing the row axis.
#' @export
row_mode <- function(x, ...) {
  UseMethod("row_mode")
}

#' Return group metadata.
#'
#' @param x A cfDNAlab object with a group axis.
#' @param ... Reserved for future methods.
#'
#' @return A data frame with one row per group.
#' @export
groups <- function(x, ...) {
  UseMethod("groups")
}

#' Return midpoint length-bin metadata.
#'
#' @param x A cfDNAlab midpoint-profile object.
#' @param ... Reserved for future methods.
#'
#' @return A data frame with one row per length bin.
#' @export
length_bins <- function(x, ...) {
  UseMethod("length_bins")
}

#' Return midpoint position-bin metadata.
#'
#' @param x A cfDNAlab midpoint-profile object.
#' @param ... Reserved for future methods.
#'
#' @return A data frame with one row per position bin.
#' @export
positions <- function(x, ...) {
  UseMethod("positions")
}

#' Return window metadata.
#'
#' @param x A cfDNAlab object with window rows.
#' @param ... Reserved for future methods.
#'
#' @return A data frame with one row per window.
#' @export
windows <- function(x, ...) {
  UseMethod("windows")
}

#' Look up a group index.
#'
#' @param x A cfDNAlab object with group labels.
#' @param ... Method-specific lookup arguments.
#'
#' @return A scalar one-based integer group index.
#' @export
group_idx <- function(x, ...) {
  UseMethod("group_idx")
}

#' Look up the midpoint length-bin index containing a fragment length.
#'
#' @param x A cfDNAlab midpoint-profile object.
#' @param ... Method-specific lookup arguments.
#'
#' @return A scalar one-based integer length-bin index.
#' @export
length_bin_idx <- function(x, ...) {
  UseMethod("length_bin_idx")
}

#' Return end-motif metadata.
#'
#' @param x A cfDNAlab end-motif object.
#' @param ... Reserved for future methods.
#'
#' @return A data frame with one row per motif.
#' @export
motifs <- function(x, ...) {
  UseMethod("motifs")
}

#' Look up a motif index.
#'
#' @param x A cfDNAlab end-motif object.
#' @param ... Method-specific lookup arguments.
#'
#' @return A scalar one-based integer motif index.
#' @export
motif_idx <- function(x, ...) {
  UseMethod("motif_idx")
}

#' Test whether an end-motif label exists.
#'
#' @param x A cfDNAlab end-motif object.
#' @param ... Method-specific lookup arguments.
#'
#' @return A scalar logical.
#' @export
has_motif <- function(x, ...) {
  UseMethod("has_motif")
}

#' Return one midpoint profile as an array vector.
#'
#' @param x A cfDNAlab midpoint-profile object.
#' @param ... Method-specific profile selection arguments.
#'
#' @return A numeric vector with one value per position bin.
#' @export
profile_array <- function(x, ...) {
  UseMethod("profile_array")
}

#' Return one midpoint profile as a data frame.
#'
#' @param x A cfDNAlab midpoint-profile object.
#' @param ... Method-specific profile selection arguments.
#'
#' @return A data frame with one row per position bin.
#' @export
profile_data_frame <- function(x, ...) {
  UseMethod("profile_data_frame")
}

#' Return the full midpoint count array.
#'
#' @param x A cfDNAlab midpoint-profile object.
#' @param ... Reserved for future methods.
#'
#' @return A three-dimensional numeric array.
#' @export
midpoint_array <- function(x, ...) {
  UseMethod("midpoint_array")
}

#' Return end-motif counts as a dense matrix.
#'
#' Sparse stores are not densified unless the method explicitly receives
#' `allow_densify = TRUE`.
#'
#' @param x A cfDNAlab end-motif object.
#' @param ... Method-specific arguments.
#'
#' @return A dense numeric matrix.
#' @export
dense_counts_matrix <- function(x, ...) {
  UseMethod("dense_counts_matrix")
}

#' Return global end-motif counts as a named vector.
#'
#' Sparse stores are not densified unless the method explicitly receives
#' `allow_densify = TRUE`.
#'
#' @param x A cfDNAlab global end-motif object.
#' @param ... Method-specific arguments.
#'
#' @return A named numeric vector with one value per motif.
#' @export
dense_counts_vector <- function(x, ...) {
  UseMethod("dense_counts_vector")
}

#' Return dense end-motif counts as a data frame.
#'
#' Sparse stores are not densified unless the method explicitly receives
#' `allow_densify = TRUE`.
#'
#' @param x A cfDNAlab end-motif object.
#' @param ... Method-specific arguments.
#'
#' @return A data frame containing motif metadata and counts.
#' @export
dense_data_frame <- function(x, ...) {
  UseMethod("dense_data_frame")
}

#' Return dense end-motif counts for one window.
#'
#' Sparse stores are not densified unless the method explicitly receives
#' `allow_densify = TRUE`.
#'
#' @param x A cfDNAlab windowed end-motif object.
#' @param ... Method-specific window selection arguments.
#'
#' @return A data frame with one row per motif.
#' @export
dense_data_frame_for_window <- function(x, ...) {
  UseMethod("dense_data_frame_for_window")
}

#' Return dense end-motif counts for one group.
#'
#' Sparse stores are not densified unless the method explicitly receives
#' `allow_densify = TRUE`.
#'
#' @param x A cfDNAlab grouped end-motif object.
#' @param ... Method-specific group selection arguments.
#'
#' @return A data frame with one row per motif.
#' @export
dense_data_frame_for_group <- function(x, ...) {
  UseMethod("dense_data_frame_for_group")
}

#' Return dense end-motif counts for one motif.
#'
#' Sparse stores are not densified unless the method explicitly receives
#' `allow_densify = TRUE`.
#'
#' @param x A cfDNAlab end-motif object.
#' @param ... Method-specific motif selection arguments.
#'
#' @return A data frame with one row per count row.
#' @export
dense_data_frame_for_motif <- function(x, ...) {
  UseMethod("dense_data_frame_for_motif")
}

#' Return end-motif counts as a sparse matrix.
#'
#' Sparse stores are converted directly from their stored COO arrays. Dense
#' stores are read into memory before conversion.
#'
#' @param x A cfDNAlab end-motif object.
#' @param ... Reserved for future methods.
#'
#' @return A `Matrix` sparse matrix.
#' @export
sparse_counts_matrix <- function(x, ...) {
  UseMethod("sparse_counts_matrix")
}

#' Return stored sparse end-motif counts as a data frame.
#'
#' @param x A cfDNAlab sparse end-motif object.
#' @param ... Reserved for future methods.
#'
#' @return A data frame with one row per stored non-zero count.
#' @export
sparse_data_frame <- function(x, ...) {
  UseMethod("sparse_data_frame")
}

#' Return stored sparse end-motif counts for one window.
#'
#' @param x A cfDNAlab windowed sparse end-motif object.
#' @param ... Method-specific window selection arguments.
#'
#' @return A data frame with one row per stored non-zero count.
#' @export
sparse_data_frame_for_window <- function(x, ...) {
  UseMethod("sparse_data_frame_for_window")
}

#' Return stored sparse end-motif counts for one group.
#'
#' @param x A cfDNAlab grouped sparse end-motif object.
#' @param ... Method-specific group selection arguments.
#'
#' @return A data frame with one row per stored non-zero count.
#' @export
sparse_data_frame_for_group <- function(x, ...) {
  UseMethod("sparse_data_frame_for_group")
}

#' Return stored sparse end-motif counts for one motif.
#'
#' @param x A cfDNAlab sparse end-motif object.
#' @param ... Method-specific motif selection arguments.
#'
#' @return A data frame with one row per stored non-zero count.
#' @export
sparse_data_frame_for_motif <- function(x, ...) {
  UseMethod("sparse_data_frame_for_motif")
}
