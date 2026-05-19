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
group_metadata <- function(x, ...) {
  UseMethod("group_metadata")
}

#' Return length-bin metadata.
#'
#' @param x A cfDNAlab object with length bins.
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
#' @return A data frame with one row per window. Public genomic window
#'   metadata uses `window_idx`, `chrom`, `start`, and `end` columns.
#' @export
window_metadata <- function(x, ...) {
  UseMethod("window_metadata")
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

#' Look up the length-bin index containing a fragment length.
#'
#' @param x A cfDNAlab object with length bins.
#' @param ... Method-specific lookup arguments.
#'
#' @return A scalar one-based integer length-bin index.
#'
#' @details
#' Errors if no length bin contains the requested fragment length.
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

#' Return midpoint profiles as a data frame.
#'
#' @param x A cfDNAlab midpoint-profile object.
#' @param ... Method-specific profile selection arguments.
#'
#' @return A data frame with one row per selected group, length bin, and
#'   position bin.
#' @export
midpoint_data_frame <- function(x, ...) {
  UseMethod("midpoint_data_frame")
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

#' Return length-count values as a matrix.
#'
#' @param x A cfDNAlab length-count object.
#' @param ... Reserved for future methods.
#'
#' @return A numeric matrix with one row per output unit and one column per
#'   length bin.
#' @export
length_counts_matrix <- function(x, ...) {
  UseMethod("length_counts_matrix")
}

#' Return global length-count values as a vector.
#'
#' @param x A cfDNAlab global length-count object.
#' @param ... Reserved for future methods.
#'
#' @return A named numeric vector with one value per length bin.
#' @export
length_counts_vector <- function(x, ...) {
  UseMethod("length_counts_vector")
}

#' Return length-count values as a data frame.
#'
#' @param x A cfDNAlab length-count object.
#' @param ... Method-specific selection arguments.
#'
#' @return A data frame with length-bin metadata and count-derived values.
#' @export
length_data_frame <- function(x, ...) {
  UseMethod("length_data_frame")
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

#' Return end-motif counts as a data frame.
#'
#' Sparse outputs return stored non-zero rows unless the method explicitly
#' receives `densify = TRUE`. Densifying adds explicit zero-count rows for
#' selected observed motifs. Dense outputs always include zero counts.
#'
#' @param x A cfDNAlab end-motif object.
#' @param ... Method-specific selection arguments.
#'
#' @return A data frame containing row metadata, motif metadata, and counts.
#' @export
end_motif_data_frame <- function(x, ...) {
  UseMethod("end_motif_data_frame")
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
