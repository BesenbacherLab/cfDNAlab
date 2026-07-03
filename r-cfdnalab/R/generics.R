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

#' Return the output storage mode.
#'
#' @param x A cfDNAlab object with dense or sparse matrix output.
#' @param ... Reserved for future methods.
#'
#' @return A scalar character value, currently `"dense"` or `"sparse_coo"`.
#' @export
storage_mode <- function(x, ...) {
  UseMethod("storage_mode")
}

#' Return the output row mode.
#'
#' @param x A cfDNAlab object with matrix rows.
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

#' Return motif-axis metadata.
#'
#' For ordinary end-motif stores, the `motif` column contains concrete motif
#' labels. For grouped motifs-file output, the same column contains user-defined
#' group names from the motif axis. Reference k-mer stores use the same column
#' for concrete k-mers or k-mer group names. For observed-only reference k-mer
#' output, this is the combined set of motifs or motifs-file targets observed
#' anywhere in the output.
#'
#' @param x A cfDNAlab object with a motif axis.
#' @param ... Reserved for future methods.
#'
#' @return A data frame with one row per motif-axis label.
#' @export
motifs <- function(x, ...) {
  UseMethod("motifs")
}

#' Look up a motif index.
#'
#' @param x A cfDNAlab object with motif labels.
#' @param ... Method-specific lookup arguments.
#'
#' @return A scalar one-based integer motif index.
#' @export
motif_idx <- function(x, ...) {
  UseMethod("motif_idx")
}

#' Test whether a motif label exists.
#'
#' @param x A cfDNAlab object with motif labels.
#' @param ... Method-specific lookup arguments.
#'
#' @return A scalar logical.
#' @export
has_motif <- function(x, ...) {
  UseMethod("has_motif")
}

#' Return the reference k-mer motif-axis kind.
#'
#' @param x A cfDNAlab reference k-mer object.
#' @param ... Reserved for future methods.
#'
#' @return A scalar character value, either `"motif"` or `"motif_group"`.
#' @export
motif_axis_kind <- function(x, ...) {
  UseMethod("motif_axis_kind")
}

#' Return the reference k-mer size.
#'
#' @param x A cfDNAlab reference k-mer object.
#' @param ... Reserved for future methods.
#'
#' @return A scalar integer k-mer size.
#' @export
kmer_size <- function(x, ...) {
  UseMethod("kmer_size")
}

#' Return whether reference k-mers were canonicalized.
#'
#' @param x A cfDNAlab reference k-mer object.
#' @param ... Reserved for future methods.
#'
#' @return A scalar logical.
#' @export
canonical <- function(x, ...) {
  UseMethod("canonical")
}

#' Return whether all requested reference k-mer motifs were kept.
#'
#' For full k-mer output, this means every A/C/G/T k-mer for the requested k.
#' For motifs-file output, this means every target from the motifs file.
#'
#' @param x A cfDNAlab reference k-mer object.
#' @param ... Reserved for future methods.
#'
#' @return A scalar logical.
#' @export
all_motifs <- function(x, ...) {
  UseMethod("all_motifs")
}

#' Return the reference k-mer window assignment rule.
#'
#' @param x A cfDNAlab reference k-mer object.
#' @param ... Reserved for future methods.
#'
#' @return A scalar character value.
#' @export
assign_by <- function(x, ...) {
  UseMethod("assign_by")
}

#' Return the reference contig footprint.
#'
#' @param x A cfDNAlab reference k-mer object.
#' @param ... Reserved for future methods.
#'
#' @return JSON-decoded reference contig footprint metadata.
#' @export
reference_contig_footprint <- function(x, ...) {
  UseMethod("reference_contig_footprint")
}

#' Return reference k-mer row scaling factors.
#'
#' Reference k-mer outputs store frequencies. Multiplying a row's frequency by
#' its `row_scaling_factor` gives the reconstructed count for that row.
#'
#' @param x A cfDNAlab reference k-mer object.
#' @param ... Reserved for future methods.
#'
#' @return A data frame with row metadata and `row_scaling_factor`.
#' @export
row_scaling_factors <- function(x, ...) {
  UseMethod("row_scaling_factors")
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

#' Return counts as a dense matrix.
#'
#' Sparse output stores only non-zero values. These methods do not create a
#' zero-filled dense matrix from sparse output unless they explicitly receive
#' `allow_densify = TRUE`. For objects with a motif axis, densifying fills
#' zeroes only across the labels returned by `motifs(x)`.
#'
#' @param x A cfDNAlab object with count values or reconstructable counts.
#' @param ... Method-specific arguments.
#'
#' @return A dense numeric matrix.
#' @export
dense_counts_matrix <- function(x, ...) {
  UseMethod("dense_counts_matrix")
}

#' Return global counts as a named vector.
#'
#' Sparse output stores only non-zero values. These methods do not create a
#' zero-filled dense vector from sparse output unless they explicitly receive
#' `allow_densify = TRUE`. For objects with a motif axis, densifying fills
#' zeroes only across the labels returned by `motifs(x)`.
#'
#' @param x A cfDNAlab global object with count values or reconstructable counts.
#' @param ... Method-specific arguments.
#'
#' @return A named numeric vector with one value per motif.
#' @export
dense_counts_vector <- function(x, ...) {
  UseMethod("dense_counts_vector")
}

#' Return reference k-mer frequencies as a dense matrix.
#'
#' Sparse output stores only non-zero frequencies. This method does not create
#' a zero-filled dense matrix from sparse output unless it explicitly receives
#' `allow_densify = TRUE`. Densifying fills zeroes only across the motif axis
#' returned by `motifs(x)`.
#'
#' @param x A cfDNAlab reference k-mer object.
#' @param ... Method-specific arguments.
#'
#' @return A dense numeric matrix.
#' @export
dense_frequencies_matrix <- function(x, ...) {
  UseMethod("dense_frequencies_matrix")
}

#' Return global reference k-mer frequencies as a named vector.
#'
#' Sparse output stores only non-zero frequencies. This method does not create
#' a zero-filled dense vector from sparse output unless it explicitly receives
#' `allow_densify = TRUE`. Densifying fills zeroes only across the motif axis
#' returned by `motifs(x)`.
#'
#' @param x A cfDNAlab global reference k-mer object.
#' @param ... Method-specific arguments.
#'
#' @return A named numeric vector with one value per motif.
#' @export
dense_frequencies_vector <- function(x, ...) {
  UseMethod("dense_frequencies_vector")
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

#' Return reference k-mer frequencies and reconstructed counts as a data frame.
#'
#' Sparse output stores only non-zero frequencies. By default, sparse output
#' returns those stored rows only. With `densify = TRUE`, the data frame also
#' includes zero-frequency rows for the selected rows and the selected motifs
#' returned by `motifs(x)`. For observed-only output, those selected labels are
#' the combined set observed anywhere in the output. Densifying does not add
#' every possible k-mer unless `all_motifs(x)` is `TRUE`. Dense output always
#' includes zeroes.
#'
#' @param x A cfDNAlab reference k-mer object.
#' @param ... Method-specific selection arguments.
#'
#' @return A data frame containing row metadata, motif metadata, `frequency`,
#'   and reconstructed `count`.
#' @export
ref_kmer_data_frame <- function(x, ...) {
  UseMethod("ref_kmer_data_frame")
}

#' Return counts as a sparse matrix.
#'
#' Sparse output is returned without building a zero-filled dense matrix. Dense
#' output is read into memory before conversion to a sparse matrix.
#'
#' @param x A cfDNAlab object with count values or reconstructable counts.
#' @param ... Reserved for future methods.
#'
#' @return A `Matrix` sparse matrix.
#' @export
sparse_counts_matrix <- function(x, ...) {
  UseMethod("sparse_counts_matrix")
}

#' Return reference k-mer frequencies as a sparse matrix.
#'
#' Sparse output is returned without building a zero-filled dense matrix. Dense
#' output is read into memory before conversion to a sparse matrix.
#'
#' @param x A cfDNAlab reference k-mer object.
#' @param ... Method-specific arguments.
#'
#' @return A `Matrix` sparse matrix.
#' @export
sparse_frequencies_matrix <- function(x, ...) {
  UseMethod("sparse_frequencies_matrix")
}
