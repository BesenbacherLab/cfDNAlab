#' Correct end-motif counts for reference k-mer composition.
#'
#' This helper starts from `end_motif_data_frame()` and adds `corrected_count`
#' and `corrected_frequency`. Formula internals such as reference frequencies
#' and correction factors are not returned.
#'
#' Motif labels are matched to reference k-mers by removing the `_` separator,
#' for example `AT_CG -> ATCG`. Motif-group outputs are matched directly by
#' group label.
#'
#' Reference k-mer output is read without densifying. For sparse reference
#' output, omitted row/motif pairs are treated as zero frequency.
#'
#' Reference correction divides each observed end-motif count by a
#' reference-based correction factor for the matched row. This factor is
#' computed from the motif frequencies in the reference k-mer output and
#' normalized so a uniform reference composition leaves counts unchanged.
#' Motifs that are common in the reference row are scaled down. Motifs that are
#' rare in the reference row are scaled up. Only motifs with a positive
#' reference frequency contribute to the row's correction support.
#'
#' Two-sided correction modes:
#'
#' When motif labels contain both outside and inside bases, such as `"AC_GT"`,
#' `two_sided_correction` chooses both the motif labels in the result and the
#' correction factor used for each returned count.
#'
#' - `"joint"` keeps full labels such as `"AC_GT"` and corrects each count
#'   using the exact reference k-mer `"ACGT"`.
#'
#' - `"split"` keeps full labels such as `"AC_GT"`, but calculates the
#'   correction factor from the two sides separately. For `"AC_GT"`, separate
#'   correction factors are calculated for outside label `"AC"` and inside
#'   label `"GT"`. Those two correction factors are multiplied and applied to
#'   the observed `"AC_GT"` count. Use this when you want full two-sided motif
#'   labels in the result, but the exact full reference k-mers are too sparse or
#'   you want the reference correction to treat outside and inside sequence
#'   composition separately.
#'
#' - `"outside"` returns outside labels such as `"AC_"`. For each outside
#'   label, all full motif counts with that outside label are summed first. For
#'   example, `"AC_AA"` and `"AC_GT"` both contribute to the `"AC_"` count. That
#'   summed count is corrected using the outside label `"AC"`.
#'
#' - `"inside"` returns inside labels such as `"_GT"`. For each inside label,
#'   all full motif counts with that inside label are summed first. For example,
#'   `"AA_GT"` and `"AC_GT"` both contribute to the `"_GT"` count. That summed
#'   count is corrected using the inside label `"GT"`.
#'
#' One-sided outputs do not accept an explicit mode.
#'
#' For `"split"`, `"outside"`, and `"inside"`, side-specific reference frequencies
#' are calculated from the loaded full-length reference k-mers. For example, the
#' outside frequency for `"AC"` is the sum of frequencies for loaded k-mers with
#' prefix `"AC"`, such as `"ACTG"` and `"ACAA"`. The inside frequency for `"TG"` is
#' the corresponding sum over loaded k-mers with suffix `"TG"`. Separate shorter
#' reference k-mer runs are not required.
#'
#' A motifs file used for the reference output restricts these sums to the k-mers
#' in that file. Without a motifs file, all k-mers in the reference output can
#' contribute, including k-mers absent from the sample end-motif output.
#'
#' `corrected_frequency` is normalized from `corrected_count` over the full
#' correction-mode motif axis for each output row. Motif selection filters those
#' frequencies afterward and does not renormalize them. With
#' `unsupported_motifs = "keep_na"`, one undefined positive corrected count makes
#' all frequencies in that output row `NA`. Correction fails if division by a
#' positive reference factor would produce a non-finite corrected count.
#'
#' An observed sample motif with a positive count is unsupported when it has no
#' positive correction factor under the selected mode. By default this is an
#' error. Set
#' `unsupported_motifs = "drop"` to omit those rows, or
#' `unsupported_motifs = "keep_na"` to keep them with `NA` corrected counts.
#'
#' By default, end-motif and reference k-mer rows must match exactly. If
#' `ref_kmers` is global and `ends` is windowed or grouped, pass
#' `use_global_bias = TRUE` to apply the global reference composition to every
#' end-motif row.
#'
#' `cfdna ends` and `cfdna ref-kmers` both write forward-oriented motif labels.
#' Right-end motifs have already been reverse-complemented by `cfdna ends`.
#'
#' Window, group, and motif selectors follow the same rules as
#' `end_motif_data_frame()`. Motif selectors choose the returned end-motif
#' rows. They do not change the reference support used for correction, so
#' selecting a motif gives the same corrected value as filtering the full
#' corrected data frame afterward.
#'
#' @param ends A loaded cfDNAlab end-motif object.
#' @param ref_kmers A loaded cfDNAlab reference k-mer object built with matching
#'   reference settings.
#' @param window_idxs Optional one-based window index vector for windowed
#'   output.
#' @param groups Optional group name vector for grouped output. Use either
#'   `groups` or `group_idxs`, not both.
#' @param group_idxs Optional one-based group index vector for grouped output.
#' @param densify Whether sparse end-motif output should include explicit
#'   zero-count rows before correction.
#' @param motifs Optional end-motif label vector. Use either `motifs` or
#'   `motif_idxs`, not both.
#' @param motif_idxs Optional one-based end-motif index vector.
#' @param max_blacklisted_fraction Maximum row `blacklisted_fraction` in 0..1
#'   to keep before correction. For matched windowed or grouped references, this
#'   filters end-motif rows first. Reference rows are then matched to the
#'   remaining end-motif rows, not filtered independently by the reference
#'   file's blacklist fractions.
#' @param use_global_bias Whether a global reference k-mer output may be applied
#'   to every non-global end-motif row.
#' @param unsupported_motifs What to do when an observed sample motif has no
#'   positive correction factor under the selected mode. Use `"error"`,
#'   `"drop"`, or `"keep_na"`.
#' @param two_sided_correction Required when motif labels contain both outside
#'   and inside bases, such as `"AC_GT"`. Use `"joint"`, `"split"`,
#'   `"outside"`, or `"inside"`. Leave as `NULL` for one-sided motifs or motif
#'   groups.
#'
#' @return An end-motif data frame with `corrected_count` and
#'   `corrected_frequency`.
#' @noRd
cf_reference_corrected_end_motif_data_frame <- function(
  ends,
  ref_kmers,
  window_idxs = NULL,
  groups = NULL,
  group_idxs = NULL,
  densify = FALSE,
  motifs = NULL,
  motif_idxs = NULL,
  max_blacklisted_fraction = 1.0,
  use_global_bias = FALSE,
  unsupported_motifs = "error",
  two_sided_correction = NULL
) {
  unsupported_motifs <- cf_match_exact_argument(
    unsupported_motifs,
    c("error", "drop", "keep_na"),
    "unsupported_motifs"
  )
  cf_validate_scalar_logical(densify, "densify")
  context <- cf_prepare_reference_correction(
    ends,
    ref_kmers,
    motifs,
    motif_idxs,
    use_global_bias,
    two_sided_correction
  )
  cf_reference_corrected_end_motif_data_frame_from_context(
    ends = ends,
    ref_kmers = ref_kmers,
    context = context,
    window_idxs = window_idxs,
    groups = groups,
    group_idxs = group_idxs,
    densify = densify,
    max_blacklisted_fraction = max_blacklisted_fraction,
    unsupported_motifs = unsupported_motifs
  )
}

#' Prepare and validate shared reference-correction state.
#'
#' @return A list containing the correction mode, row columns, reference row
#'   columns, and selected output motif labels.
#' @noRd
cf_prepare_reference_correction <- function(
  ends,
  ref_kmers,
  motifs,
  motif_idxs,
  use_global_bias,
  two_sided_correction
) {
  # Validate the relationship between the stored outputs before deriving any
  # selectors, so failures do not depend on which result form the caller uses
  two_sided_correction <- cf_validate_two_sided_correction(two_sided_correction)
  if (!inherits(ends, "cfdnalab_end_motif_counts")) {
    stop("ends must be a cfDNAlab end-motif object", call. = FALSE)
  }
  if (!inherits(ref_kmers, "cfdnalab_ref_kmer_frequencies")) {
    stop("ref_kmers must be a cfDNAlab reference k-mer object", call. = FALSE)
  }
  cf_validate_scalar_logical(use_global_bias, "use_global_bias")
  if (isTRUE(use_global_bias) && !identical(ref_kmers$row_mode, "global")) {
    stop(
      "use_global_bias = TRUE requires a global reference k-mer output",
      call. = FALSE
    )
  }
  if (!identical(ends$row_mode, ref_kmers$row_mode)) {
    if (
      identical(ref_kmers$row_mode, "global") &&
        !identical(ends$row_mode, "global")
    ) {
      if (!isTRUE(use_global_bias)) {
        stop(
          "Reference k-mer output is global but end-motif output is ",
          ends$row_mode,
          ". Pass use_global_bias = TRUE to apply the global reference bias ",
          "to every end-motif row.",
          call. = FALSE
        )
      }
    } else {
      stop(
        "End-motif and reference k-mer row modes must match: ",
        ends$row_mode,
        " != ",
        ref_kmers$row_mode,
        call. = FALSE
      )
    }
  }
  cf_validate_reference_correction_motif_axes(ends, ref_kmers)

  # Resolve both the correction formula and the motif axis it produces. Side
  # modes derive a new axis, while exact and split modes retain the stored axis
  correction_mode <- cf_reference_correction_mode(
    ends,
    ref_kmers,
    two_sided_correction
  )
  row_columns <- cf_reference_correction_row_columns(ends$row_mode)
  reference_row_columns <- cf_reference_correction_reference_row_columns(
    ends,
    ref_kmers,
    use_global_bias
  )
  # Row validation uses the full stored metadata. Selection happens later, after
  # establishing that keyed correction cannot silently mismatch the two outputs
  cf_validate_reference_correction_rows(
    ends,
    ref_kmers,
    row_columns,
    reference_row_columns
  )

  list(
    correction_mode = correction_mode,
    row_columns = row_columns,
    reference_row_columns = reference_row_columns,
    selected_mode_labels = cf_selected_reference_correction_mode_labels(
      ends,
      correction_mode,
      motifs,
      motif_idxs
    )
  )
}

#' Build reference-corrected rows from validated shared state.
#'
#' @return A reference-corrected end-motif data frame.
#' @noRd
cf_reference_corrected_end_motif_data_frame_from_context <- function(
  ends,
  ref_kmers,
  context,
  window_idxs,
  groups,
  group_idxs,
  densify,
  max_blacklisted_fraction,
  unsupported_motifs
) {
  # Load every motif for the selected sample rows. Corrected frequencies need the
  # full result-axis total, and side modes need all full motifs before aggregation
  end_row_indices <- cf_reference_correction_end_row_indices(
    ends,
    window_idxs,
    groups,
    group_idxs
  )
  end_rows <- cf_end_motif_data_frame(
    ends,
    row_indices = end_row_indices,
    motif_indices = seq_along(ends$motif_idx0),
    densify = densify,
    max_blacklisted_fraction = max_blacklisted_fraction
  )
  if (nrow(end_rows) == 0L) {
    return(cf_add_empty_reference_correction_columns(end_rows))
  }

  # Preserve the public columns and first-seen sample-row order across base R
  # merges, which otherwise reorder their inputs by join keys
  output_columns <- names(end_rows)
  end_rows$.cfdnalab_row_order <- cf_reference_correction_row_order(
    end_rows,
    context$row_columns
  )

  ref_row_indices <- cf_reference_correction_ref_row_indices_from_end_rows(
    ref_kmers,
    end_rows,
    context$reference_row_columns
  )
  ref_rows <- cf_ref_kmer_data_frame(
    ref_kmers,
    row_indices = ref_row_indices,
    motif_indices = seq_along(ref_kmers$motif_idx0),
    densify = FALSE,
    max_blacklisted_fraction = 1.0
  )

  # Keep the complete reference motif set for the chosen reference rows
  # because support counts and side marginals must not depend on sample selection
  ref_rows <- cf_prepare_reference_correction_ref_rows(
    ref_rows,
    context$reference_row_columns
  )

  # Convert reference frequencies into correction factors according to whether
  # the result retains full motifs or collapses onto a selected motif side
  if (identical(context$correction_mode$mode, "exact")) {
    corrected <- cf_exact_reference_corrected_rows(
      ends,
      end_rows,
      ref_rows,
      context$reference_row_columns,
      unsupported_motifs
    )
  } else if (identical(context$correction_mode$mode, "split")) {
    corrected <- cf_split_reference_corrected_rows(
      end_rows,
      ref_rows,
      context$reference_row_columns,
      context$correction_mode,
      unsupported_motifs
    )
  } else {
    corrected <- cf_side_reference_corrected_rows(
      end_rows,
      ref_rows,
      context$reference_row_columns,
      context$correction_mode,
      output_columns,
      unsupported_motifs
    )
  }

  # Normalize before motif filtering so requesting a subset does not change the
  # frequencies that those motifs would have in the complete corrected result
  corrected <- cf_add_corrected_frequency(
    corrected,
    context$row_columns
  )
  corrected <- corrected[
    corrected$motif %in% context$selected_mode_labels,
    ,
    drop = FALSE
  ]

  # Restore requested motif order within original sample-row order after all
  # joins, aggregation, and filtering have finished
  if (length(context$selected_mode_labels) > 0L && nrow(corrected) > 0L) {
    corrected$.cfdnalab_motif_order <- match(
      corrected$motif,
      context$selected_mode_labels
    )
    corrected <- corrected[order(
      corrected$.cfdnalab_row_order,
      corrected$.cfdnalab_motif_order
    ), , drop = FALSE]
    corrected$.cfdnalab_motif_order <- NULL
  }
  row.names(corrected) <- NULL
  corrected[c(output_columns, "corrected_count", "corrected_frequency")]
}

#' Return reference-corrected counts as a dense matrix.
#'
#' @param ends End-motif object.
#' @param ref_kmers Reference k-mer object.
#' @param window_idxs Optional one-based window indices.
#' @param groups Optional group names.
#' @param group_idxs Optional one-based group indices.
#' @param motifs Optional motif labels.
#' @param motif_idxs Optional one-based motif indices.
#' @param allow_densify Whether sparse end-motif output may be densified.
#' @param max_blacklisted_fraction Maximum blacklist fraction.
#' @param use_global_bias Whether a global reference can be broadcast.
#' @param unsupported_motifs Unsupported motif policy.
#' @param two_sided_correction `NULL` for one-sided or motif-group axes.
#'   For two-sided axes, one of `"joint"`, `"split"`, `"outside"`, or
#'   `"inside"`, which determines both factor construction and the motif axis.
#'
#' @return A dense numeric matrix.
#' @noRd
cf_reference_corrected_counts_matrix <- function(
  ends,
  ref_kmers,
  window_idxs = NULL,
  groups = NULL,
  group_idxs = NULL,
  motifs = NULL,
  motif_idxs = NULL,
  allow_densify = FALSE,
  max_blacklisted_fraction = 1.0,
  use_global_bias = FALSE,
  unsupported_motifs = "error",
  two_sided_correction = NULL
) {
  cf_validate_fixed_shape_reference_correction_policy(
    unsupported_motifs,
    "dense_corrected_counts_matrix"
  )
  cf_validate_scalar_logical(allow_densify, "allow_densify")
  context <- cf_prepare_reference_correction(
    ends,
    ref_kmers,
    motifs,
    motif_idxs,
    use_global_bias,
    two_sided_correction
  )
  if (identical(ends$storage_mode, "sparse_coo") && !isTRUE(allow_densify)) {
    stop(
      "This end-motif store is sparse. Use sparse_corrected_counts_matrix() ",
      "or set allow_densify = TRUE.",
      call. = FALSE
    )
  }
  row_indices <- cf_reference_correction_end_row_indices(ends, window_idxs, groups, group_idxs)
  row_indices <- cf_apply_end_motif_blacklist_filter(ends, row_indices, max_blacklisted_fraction)
  corrected <- cf_reference_corrected_end_motif_data_frame_from_context(
    ends = ends,
    ref_kmers = ref_kmers,
    context = context,
    window_idxs = window_idxs,
    groups = groups,
    group_idxs = group_idxs,
    densify = TRUE,
    max_blacklisted_fraction = max_blacklisted_fraction,
    unsupported_motifs = unsupported_motifs
  )
  corrected_matrix <- matrix(
    corrected$corrected_count,
    nrow = length(row_indices),
    ncol = length(context$selected_mode_labels),
    byrow = TRUE
  )
  if (context$correction_mode$mode %in% c("outside", "inside")) {
    colnames(corrected_matrix) <- context$selected_mode_labels
  }
  corrected_matrix
}

#' Return reference-corrected counts as a sparse matrix.
#'
#' @param ends End-motif object.
#' @param ref_kmers Reference k-mer object.
#' @param window_idxs Optional one-based window indices.
#' @param groups Optional group names.
#' @param group_idxs Optional one-based group indices.
#' @param motifs Optional motif labels.
#' @param motif_idxs Optional one-based motif indices.
#' @param max_blacklisted_fraction Maximum blacklist fraction.
#' @param use_global_bias Whether a global reference can be broadcast.
#' @param unsupported_motifs Unsupported motif policy.
#' @param two_sided_correction `NULL` for one-sided or motif-group axes.
#'   For two-sided axes, one of `"joint"`, `"split"`, `"outside"`, or
#'   `"inside"`, which determines both factor construction and the motif axis.
#'
#' @return A `Matrix` sparse matrix.
#' @noRd
cf_sparse_reference_corrected_counts_matrix <- function(
  ends,
  ref_kmers,
  window_idxs = NULL,
  groups = NULL,
  group_idxs = NULL,
  motifs = NULL,
  motif_idxs = NULL,
  max_blacklisted_fraction = 1.0,
  use_global_bias = FALSE,
  unsupported_motifs = "error",
  two_sided_correction = NULL
) {
  cf_validate_fixed_shape_reference_correction_policy(
    unsupported_motifs,
    "sparse_corrected_counts_matrix"
  )
  context <- cf_prepare_reference_correction(
    ends,
    ref_kmers,
    motifs,
    motif_idxs,
    use_global_bias,
    two_sided_correction
  )
  row_indices <- cf_reference_correction_end_row_indices(ends, window_idxs, groups, group_idxs)
  row_indices <- cf_apply_end_motif_blacklist_filter(ends, row_indices, max_blacklisted_fraction)
  if (
    length(row_indices) == 0L ||
      length(context$selected_mode_labels) == 0L
  ) {
    return(cf_empty_sparse_reference_corrected_counts_matrix(
      length(row_indices),
      context$selected_mode_labels,
      context$correction_mode
    ))
  }
  corrected <- cf_reference_corrected_end_motif_data_frame_from_context(
    ends = ends,
    ref_kmers = ref_kmers,
    context = context,
    window_idxs = window_idxs,
    groups = groups,
    group_idxs = group_idxs,
    densify = FALSE,
    max_blacklisted_fraction = max_blacklisted_fraction,
    unsupported_motifs = unsupported_motifs
  )
  if (nrow(corrected) == 0L) {
    return(cf_empty_sparse_reference_corrected_counts_matrix(
      length(row_indices),
      context$selected_mode_labels,
      context$correction_mode
    ))
  }
  row_positions <- cf_reference_correction_row_positions(ends, row_indices)
  motif_positions <- stats::setNames(
    seq_along(context$selected_mode_labels),
    context$selected_mode_labels
  )
  corrected_values <- corrected$corrected_count
  stored <- corrected_values != 0 | is.na(corrected_values)
  if (!any(stored)) {
    return(cf_empty_sparse_reference_corrected_counts_matrix(
      length(row_indices),
      context$selected_mode_labels,
      context$correction_mode
    ))
  }
  Matrix::sparseMatrix(
    i = unname(row_positions[cf_reference_correction_row_keys(
      corrected,
      context$row_columns
    )])[stored],
    j = unname(motif_positions[corrected$motif])[stored],
    x = corrected_values[stored],
    dims = as.integer(c(length(row_indices), length(context$selected_mode_labels))),
    dimnames = if (context$correction_mode$mode %in% c("outside", "inside")) {
      list(NULL, context$selected_mode_labels)
    } else {
      NULL
    }
  )
}

#' Construct an empty fixed-shape sparse correction matrix.
#'
#' @return A sparse matrix with the requested dimensions.
#' @noRd
cf_empty_sparse_reference_corrected_counts_matrix <- function(
  row_count,
  mode_labels,
  correction_mode
) {
  Matrix::sparseMatrix(
    i = integer(),
    j = integer(),
    dims = as.integer(c(row_count, length(mode_labels))),
    dimnames = if (correction_mode$mode %in% c("outside", "inside")) {
      list(NULL, mode_labels)
    } else {
      NULL
    }
  )
}

#' Reject policies that would change a matrix axis.
#'
#' Dense and sparse corrected matrices have fixed row and motif dimensions, so
#' they cannot represent `"drop"` without removing cells or columns. They accept
#' `"error"` and `"keep_na"`. Data frame output remains the route for dropping
#' unsupported motif rows.
#'
#' @param unsupported_motifs Unsupported motif policy.
#' @param method_name Public method name.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_fixed_shape_reference_correction_policy <- function(unsupported_motifs, method_name) {
  unsupported_motifs <- cf_match_exact_argument(
    unsupported_motifs,
    c("error", "drop", "keep_na"),
    "unsupported_motifs"
  )
  if (identical(unsupported_motifs, "drop")) {
    stop(
      "unsupported_motifs = \"drop\" cannot be represented in a fixed-shape ",
      method_name,
      "() result. Use end_motif_data_frame(..., ref_kmers = ref_kmers, ",
      "unsupported_motifs = \"drop\") or unsupported_motifs = \"keep_na\".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Validate a nullable two-sided correction choice.
#'
#' `NULL` is retained for axes that do not need a two-sided choice. Non-null
#' values must match `"joint"`, `"split"`, `"outside"`, or `"inside"`
#' exactly.
#'
#' @param two_sided_correction User-supplied mode or `NULL`.
#'
#' @return `NULL` or the matched mode string.
#' @noRd
cf_validate_two_sided_correction <- function(two_sided_correction) {
  if (is.null(two_sided_correction)) {
    return(NULL)
  }
  cf_match_exact_argument(
    two_sided_correction,
    c("joint", "split", "outside", "inside"),
    "two_sided_correction"
  )
}

#' Resolve factor construction and the returned motif axis.
#'
#' Motif groups use exact matching and reject a two-sided choice. Empty motif
#' axes return exact mode without needing side widths. One-sided motifs also use
#' exact matching and reject a two-sided choice. A motif with bases on both
#' sides requires an explicit choice. `"joint"` uses exact full labels,
#' `"split"` records the inferred side widths while preserving full labels,
#' and side modes additionally build the derived outside or inside label axis.
#'
#' @param ends End-motif object.
#' @param ref_kmers Reference k-mer object.
#' @param two_sided_correction Nullable public mode value.
#'
#' @return List describing the internal correction mode.
#' @noRd
cf_reference_correction_mode <- function(ends, ref_kmers, two_sided_correction) {
  # Motif groups already define the complete correction axis and have no
  # outside/inside boundary that a two-sided mode could reinterpret
  if (identical(ends$motif_axis_kind, "motif_group")) {
    if (!is.null(two_sided_correction)) {
      stop("Motif-group end-motif outputs do not accept two_sided_correction", call. = FALSE)
    }
    return(list(mode = "exact"))
  }
  # An empty stored axis has no labels from which to infer side widths. Exact
  # mode still lets empty selections return their expected shape
  if (length(ends$motif) == 0L) {
    return(list(mode = "exact"))
  }

  # Resolve the boundary from sample labels and require every reference label to
  # use the same total width before choosing a correction formula
  widths <- cf_infer_end_motif_side_widths(ends$motif, ref_kmers$kmer_size)
  cf_validate_reference_labels_split_cleanly(ref_kmers$motif, widths$outside, widths$inside)
  if (widths$outside == 0L || widths$inside == 0L) {
    if (!is.null(two_sided_correction)) {
      stop("One-sided end-motif outputs do not accept two_sided_correction", call. = FALSE)
    }
    return(list(mode = "exact"))
  }
  # A true two-sided axis is ambiguous without an explicit choice because the
  # choice controls both the correction formula and the side-mode result shape
  if (is.null(two_sided_correction)) {
    stop("two-sided end-motif labels with both outside and inside bases require two_sided_correction", call. = FALSE)
  }
  if (identical(two_sided_correction, "joint")) {
    return(list(mode = "exact"))
  }
  # Split retains full motifs. Outside and inside modes additionally need a
  # stable derived axis for selection, matrix columns, and result ordering
  mode <- list(
    mode = two_sided_correction,
    outside_width = widths$outside,
    inside_width = widths$inside
  )
  if (two_sided_correction %in% c("outside", "inside")) {
    mode$side_labels <- cf_side_axis_labels(ends$motif, two_sided_correction)
  }
  mode
}

#' Infer and validate the outside and inside motif widths.
#'
#' Each label is split at its underscore. Its outside and inside widths must
#' sum to the reference k-mer size, and every loaded label must use the same
#' pair of widths so reference prefixes and suffixes are unambiguous.
#'
#' @param motif_labels End-motif labels.
#' @param reference_kmer_size Reference k-mer size.
#'
#' @return List with `outside` and `inside` widths.
#' @noRd
cf_infer_end_motif_side_widths <- function(motif_labels, reference_kmer_size) {
  outside_width <- NULL
  inside_width <- NULL
  for (motif_label in motif_labels) {
    parts <- cf_split_end_motif_label(motif_label)
    widths <- c(nchar(parts$outside), nchar(parts$inside))
    if (sum(widths) != reference_kmer_size) {
      stop(
        "End-motif width must match reference k-mer size (",
        reference_kmer_size,
        "): ",
        motif_label,
        call. = FALSE
      )
    }
    if (is.null(outside_width)) {
      outside_width <- widths[[1L]]
      inside_width <- widths[[2L]]
    } else if (!identical(c(outside_width, inside_width), widths)) {
      stop(
        "All end-motif labels must use the same outside and inside widths",
        call. = FALSE
      )
    }
  }
  list(outside = outside_width, inside = inside_width)
}

#' Split an end-motif label at its single underscore.
#'
#' Exactly one underscore is required. An empty outside or inside component is
#' valid, allowing one-sided labels such as `"_GT"` and `"AC_"`.
#'
#' @param motif_label End-motif label.
#'
#' @return List with `outside` and `inside`.
#' @noRd
cf_split_end_motif_label <- function(motif_label) {
  separator_positions <- gregexpr("_", motif_label, fixed = TRUE)[[1L]]
  if (
    length(separator_positions) != 1L ||
      separator_positions[[1L]] < 0L
  ) {
    stop(
      "End-motif label must contain exactly one '_' to separate outside and inside bases: ",
      motif_label,
      call. = FALSE
    )
  }
  separator_position <- separator_positions[[1L]]
  motif_width <- nchar(motif_label)
  outside <- if (separator_position == 1L) {
    ""
  } else {
    substr(motif_label, 1L, separator_position - 1L)
  }
  inside <- if (separator_position == motif_width) {
    ""
  } else {
    substr(motif_label, separator_position + 1L, motif_width)
  }
  list(outside = outside, inside = inside)
}

#' Require every reference label to cover the inferred two sides.
#'
#' A reference label has no underscore, so its length must equal
#' `outside_width + inside_width` before it can be split by prefix and suffix.
#'
#' @param reference_labels Reference motif labels.
#' @param outside_width Outside width.
#' @param inside_width Inside width.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_reference_labels_split_cleanly <- function(
  reference_labels,
  outside_width,
  inside_width
) {
  expected_width <- outside_width + inside_width
  bad_labels <- reference_labels[nchar(reference_labels) != expected_width]
  if (length(bad_labels) > 0L) {
    stop(
      "Reference motif label must split into outside width ",
      outside_width,
      " and inside width ",
      inside_width,
      ": ",
      bad_labels[[1L]],
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Build a derived side axis in loaded full-motif order.
#'
#' Outside mode converts each full label to forms such as `"AC_"`, while inside
#' mode produces forms such as `"_GT"`. Repeated side labels are removed while
#' preserving the first loaded full motif that introduced each label.
#'
#' @param motif_labels Concrete joint motif labels.
#' @param side_mode `"outside"` or `"inside"`.
#'
#' @return Side labels.
#' @noRd
cf_side_axis_labels <- function(motif_labels, side_mode) {
  motif_parts <- lapply(motif_labels, cf_split_end_motif_label)
  side_labels <- if (identical(side_mode, "outside")) {
    paste0(vapply(motif_parts, `[[`, character(1), "outside"), "_")
  } else {
    paste0("_", vapply(motif_parts, `[[`, character(1), "inside"))
  }
  unique(side_labels)
}

#' Resolve labels against the axis returned by the correction mode.
#'
#' Exact and split modes select from the stored full-motif axis. Outside and
#' inside modes select from a derived side-label axis, reject source motif
#' indices because they do not identify that axis, validate requested labels,
#' and preserve the caller's requested label order.
#'
#' @param ends End-motif object.
#' @param correction_mode Internal correction mode.
#' @param motifs Optional motif labels.
#' @param motif_idxs Optional motif indices.
#'
#' @return Selected mode labels.
#' @noRd
cf_selected_reference_correction_mode_labels <- function(
  ends,
  correction_mode,
  motifs,
  motif_idxs
) {
  # Exact and split results retain the stored motif axis, so stored labels and
  # indices keep their usual meaning
  if (correction_mode$mode %in% c("exact", "split")) {
    return(ends$motif[cf_resolve_end_motif_indices(ends, motifs, motif_idxs)])
  }
  # A side label represents several stored motifs, so a stored motif index cannot
  # identify a derived side-mode column
  if (!is.null(motif_idxs)) {
    stop(
      "motif index selectors are not supported for outside or inside reference correction",
      call. = FALSE
    )
  }
  if (is.null(motifs)) {
    return(correction_mode$side_labels)
  }
  # Validate requested labels against the derived axis without reordering them
  if (!is.character(motifs) || anyNA(motifs)) {
    stop("motifs must contain character strings", call. = FALSE)
  }
  if (any(duplicated(motifs))) {
    stop("motifs contains duplicate values", call. = FALSE)
  }
  unknown <- setdiff(motifs, correction_mode$side_labels)
  if (length(unknown) > 0L) {
    stop("Side-mode motif axis has no label ", sQuote(unknown[[1L]]), call. = FALSE)
  }
  motifs
}

#' Prepare only the reference rows and columns needed for correction.
#'
#' Rename reference `motif` and `frequency` columns to avoid collisions with
#' sample columns and retain the row-identifying join columns. The caller has
#' already loaded only the reference rows needed by the selected sample rows.
#'
#' @param ref_rows Reference rows.
#' @param reference_row_columns Reference row columns.
#'
#' @return Reference rows with correction column names.
#' @noRd
cf_prepare_reference_correction_ref_rows <- function(
  ref_rows,
  reference_row_columns
) {
  names(ref_rows)[names(ref_rows) == "motif"] <- "reference_motif"
  names(ref_rows)[names(ref_rows) == "frequency"] <- "reference_frequency"
  ref_rows[c(reference_row_columns, "reference_motif", "reference_frequency")]
}

#' Count positive frequencies within each reference row.
#'
#' Global reference composition receives the same scalar count on every row.
#' Keyed output is grouped by its row-identifying metadata. Empty input returns
#' an aligned empty integer vector.
#'
#' @return Integer support count for each input frequency.
#' @noRd
cf_positive_support_counts <- function(frequencies, data_frame, row_columns) {
  positive_frequency <- as.integer(frequencies > 0)
  if (length(positive_frequency) == 0L) {
    return(integer())
  }
  if (length(row_columns) == 0L) {
    return(rep.int(sum(positive_frequency), length(positive_frequency)))
  }
  row_keys <- cf_reference_correction_row_keys(data_frame, row_columns)
  as.integer(stats::ave(positive_frequency, row_keys, FUN = sum))
}

#' Correct each sample motif with its matching full reference label.
#'
#' Motif-group labels match directly. Sequence labels match after removing the
#' underscore. For each reference row, the correction factor is the matching
#' frequency times the number of positive reference motifs, which expresses the
#' frequency relative to a uniform reference row. Each sample count is divided
#' by that factor after applying the unsupported-reference policy.
#'
#' @return Corrected rows.
#' @noRd
cf_exact_reference_corrected_rows <- function(
  ends,
  end_rows,
  ref_rows,
  reference_row_columns,
  unsupported_motifs
) {
  # Translate concrete sample labels to reference labels. Motif-group labels are
  # already shared verbatim by both outputs
  if (identical(ends$motif_axis_kind, "motif_group")) {
    end_rows$reference_motif <- end_rows$motif
  } else {
    end_rows$reference_motif <- gsub("_", "", end_rows$motif, fixed = TRUE)
  }
  merge_columns <- c(reference_row_columns, "reference_motif")
  if (any(duplicated(ref_rows[merge_columns]))) {
    stop("Reference k-mer rows are not unique for row and motif labels", call. = FALSE)
  }

  # Attach support before joining so frequency and support arrive together
  ref_rows$number_of_supported_motifs <- cf_positive_support_counts(
    ref_rows$reference_frequency,
    ref_rows,
    reference_row_columns
  )

  # Attach each motif's frequency and support count while restoring sample order
  # after merge(), then represent absent sparse reference pairs as unsupported
  end_rows$.cfdnalab_order <- seq_len(nrow(end_rows))
  corrected <- merge(
    end_rows,
    ref_rows,
    by = merge_columns,
    all.x = TRUE,
    sort = FALSE
  )
  corrected <- corrected[order(corrected$.cfdnalab_order), , drop = FALSE]
  corrected$.cfdnalab_order <- NULL
  corrected$reference_frequency[is.na(corrected$reference_frequency)] <- 0
  missing_support <- is.na(corrected$number_of_supported_motifs)
  corrected$number_of_supported_motifs[missing_support] <- 0L
  # With N supported motifs, frequency 1/N produces a denominator of 1
  corrected$reference_denominator <- (
    corrected$reference_frequency * corrected$number_of_supported_motifs
  )
  cf_apply_reference_denominator_policy(
    corrected,
    unsupported_motifs
  )
}

#' Correct full motifs with independently calculated side factors.
#'
#' Parse each sample motif into its outside and inside labels, calculate factors
#' from aggregated reference prefix and suffix frequencies, and join both factors
#' to every sample row. The full motif axis is retained. Its correction factor is
#' the product of the matching outside and inside factors.
#'
#' @return Corrected rows.
#' @noRd
cf_split_reference_corrected_rows <- function(
  end_rows,
  ref_rows,
  reference_row_columns,
  correction_mode,
  unsupported_motifs
) {
  # Add the two labels used to join each full sample motif to side marginals
  end_rows <- cf_add_end_motif_sides(
    end_rows,
    correction_mode$outside_width,
    correction_mode$inside_width
  )
  # Derive outside and inside factors independently from the full reference axis
  outside_denominator <- cf_side_reference_denominator(
    ref_rows,
    reference_row_columns,
    "outside",
    correction_mode$outside_width,
    correction_mode$inside_width
  )
  inside_denominator <- cf_side_reference_denominator(
    ref_rows,
    reference_row_columns,
    "inside",
    correction_mode$outside_width,
    correction_mode$inside_width
  )
  # Retain the full sample axis and combine the two side effects multiplicatively
  corrected <- cf_merge_side_denominator(
    end_rows,
    outside_denominator,
    reference_row_columns,
    "outside"
  )
  corrected <- cf_merge_side_denominator(
    corrected,
    inside_denominator,
    reference_row_columns,
    "inside"
  )
  corrected$reference_denominator <- (
    corrected$outside_reference_denominator *
      corrected$inside_reference_denominator
  )
  cf_apply_reference_denominator_policy(
    corrected,
    unsupported_motifs
  )
}

#' Aggregate sample counts to a side axis before correcting them.
#'
#' Relabel each full motif as an outside label such as `"AC_"` or an inside
#' label such as `"_GT"`, then sum counts with the same selected-row metadata
#' and side label. The returned motif axis is this derived side axis. Each
#' aggregated count is divided by the matching side correction factor.
#'
#' @return Corrected rows.
#' @noRd
cf_side_reference_corrected_rows <- function(
  end_rows,
  ref_rows,
  reference_row_columns,
  correction_mode,
  output_columns,
  unsupported_motifs
) {
  # Replace each full motif with the outside or inside label retained by this mode
  end_rows <- cf_add_end_motif_sides(
    end_rows,
    correction_mode$outside_width,
    correction_mode$inside_width
  )
  side_column <- if (identical(correction_mode$mode, "outside")) "outside" else "inside"
  end_rows$motif <- if (identical(correction_mode$mode, "outside")) {
    paste0(end_rows$outside, "_")
  } else {
    paste0("_", end_rows$inside)
  }
  end_rows$motif_idx <- match(end_rows$motif, correction_mode$side_labels)

  # Save row metadata once because aggregation below keeps only row order, the
  # derived side label, and its summed sample count
  row_metadata_columns <- setdiff(output_columns, c("motif_idx", "motif", "count"))
  row_metadata <- unique(end_rows[c(".cfdnalab_row_order", row_metadata_columns)])
  aggregate_columns <- c(".cfdnalab_row_order", "motif_idx", "motif", side_column)
  # Sum all full motifs that share the selected side before applying correction
  aggregated <- stats::aggregate(
    end_rows["count"],
    by = end_rows[aggregate_columns],
    FUN = sum
  )
  aggregated <- merge(
    aggregated,
    row_metadata,
    by = ".cfdnalab_row_order",
    all.x = TRUE,
    sort = FALSE
  )
  aggregated <- aggregated[order(
    aggregated$.cfdnalab_row_order,
    aggregated$motif_idx
  ), , drop = FALSE]
  # Calculate marginals from the complete reference axis and attach only the side
  # factor represented by the derived result axis
  side_denominator <- cf_side_reference_denominator(
    ref_rows,
    reference_row_columns,
    side_column,
    correction_mode$outside_width,
    correction_mode$inside_width
  )
  corrected <- cf_merge_side_denominator(
    aggregated,
    side_denominator,
    reference_row_columns,
    side_column
  )
  denominator_column <- paste0(side_column, "_reference_denominator")
  corrected$reference_denominator <- corrected[[denominator_column]]
  cf_apply_reference_denominator_policy(
    corrected,
    unsupported_motifs
  )
}

#' Parse end-motif labels into validated outside and inside columns.
#'
#' The parsed components must match the widths inferred for the complete motif
#' axis. This catches inconsistent labels before side aggregation or joins can
#' silently assign them to the wrong reference prefix or suffix.
#'
#' @return End rows with side columns.
#' @noRd
cf_add_end_motif_sides <- function(end_rows, outside_width, inside_width) {
  split_labels <- lapply(end_rows$motif, cf_split_end_motif_label)
  end_rows$outside <- vapply(split_labels, `[[`, character(1), "outside")
  end_rows$inside <- vapply(split_labels, `[[`, character(1), "inside")
  invalid_width <- (
    nchar(end_rows$outside) != outside_width |
      nchar(end_rows$inside) != inside_width
  )
  if (any(invalid_width)) {
    stop("End-motif label does not match inferred side widths", call. = FALSE)
  }
  end_rows
}

#' Build the reference correction-factor table for one motif side.
#'
#' Reduce each full reference label to the requested prefix or suffix, aggregate
#' its frequencies within each reference row, and multiply each marginal by the
#' number of positive-frequency side labels in that row. Uniform side composition
#' therefore gives correction factor one.
#'
#' @return A side correction-factor data frame.
#' @noRd
cf_side_reference_denominator <- function(
  ref_rows,
  reference_row_columns,
  side_column,
  outside_width,
  inside_width
) {
  ref_rows[[side_column]] <- if (identical(side_column, "outside")) {
    substr(ref_rows$reference_motif, 1L, outside_width)
  } else {
    substr(
      ref_rows$reference_motif,
      outside_width + 1L,
      outside_width + inside_width
    )
  }

  group_columns <- c(reference_row_columns, side_column)
  denominator_column <- paste0(side_column, "_reference_denominator")
  # Preserve the expected join columns when sparse input stores no reference rows
  if (nrow(ref_rows) == 0L) {
    denominators <- ref_rows[integer(0), group_columns, drop = FALSE]
    denominators[[denominator_column]] <- numeric()
    return(denominators)
  }

  # Collapse full reference motifs into marginal frequencies
  side_frequencies <- stats::aggregate(
    ref_rows["reference_frequency"],
    by = ref_rows[group_columns],
    FUN = sum
  )
  frequency_column <- ".cfdnalab_side_reference_frequency"
  names(side_frequencies)[ncol(side_frequencies)] <- frequency_column

  side_support_count <- cf_positive_support_counts(
    side_frequencies[[frequency_column]],
    side_frequencies,
    reference_row_columns
  )
  side_frequencies[[denominator_column]] <- (
    side_frequencies[[frequency_column]] * side_support_count
  )
  side_frequencies[c(group_columns, denominator_column)]
}

#' Join one side factor without changing sample-row order.
#'
#' A left join retains every sample row. The saved input position is restored
#' after `merge()`, and a missing factor becomes zero so unsupported positive
#' counts are handled by the selected policy.
#'
#' @return End rows with the selected side correction-factor column.
#' @noRd
cf_merge_side_denominator <- function(
  end_rows,
  side_denominator,
  reference_row_columns,
  side_column
) {
  end_rows$.cfdnalab_order <- seq_len(nrow(end_rows))
  denominator_column <- paste0(side_column, "_reference_denominator")
  side_columns <- c(reference_row_columns, side_column, denominator_column)
  corrected <- merge(
    end_rows,
    side_denominator[side_columns],
    by = c(reference_row_columns, side_column),
    all.x = TRUE,
    sort = FALSE
  )
  corrected <- corrected[order(corrected$.cfdnalab_order), , drop = FALSE]
  corrected$.cfdnalab_order <- NULL
  corrected[[denominator_column]][is.na(corrected[[denominator_column]])] <- 0
  row.names(corrected) <- NULL
  corrected
}

#' Apply unsupported-reference policy and divide supported counts.
#'
#' A non-positive correction factor is unsupported. Under `"error"`, only an
#' unsupported positive sample count is an error. `"drop"` removes every row
#' lacking support, including zero-count rows. `"keep_na"` retains the axis and
#' marks unsupported positive counts as `NA`, while unsupported zero counts
#' remain zero. Supported counts are divided by their correction factor.
#'
#' @return Corrected rows with `corrected_count`.
#' @noRd
cf_apply_reference_denominator_policy <- function(
  corrected,
  unsupported_motifs
) {
  denominator_column <- "reference_denominator"
  # A missing denominator means the left join found no reference support. Treat
  # it exactly like an explicit zero denominator
  corrected[[denominator_column]][is.na(corrected[[denominator_column]])] <- 0

  # Keep separate masks because unsupported zero counts are harmless under error
  # and keep_na, while drop intentionally removes the full unsupported axis
  unsupported_reference <- corrected[[denominator_column]] <= 0
  positive_unsupported_reference <- unsupported_reference & corrected$count > 0

  # Fail before changing the rows so the message includes every affected motif
  if (any(positive_unsupported_reference) && identical(unsupported_motifs, "error")) {
    unsupported_labels <- sort(unique(corrected$motif[positive_unsupported_reference]))
    stop(
      "Positive-count end motifs have no positive reference-based correction factor: ",
      paste(unsupported_labels, collapse = ", "),
      ". Pass unsupported_motifs = \"drop\" to omit those rows, or ",
      "unsupported_motifs = \"keep_na\" to keep them with NA corrected counts.",
      call. = FALSE
    )
  }

  # Drop uses reference support rather than observed count, so zero-count rows
  # without support are also removed and the result may have a variable shape
  if (identical(unsupported_motifs, "drop")) {
    corrected <- corrected[!unsupported_reference, , drop = FALSE]
  }

  # Initialize to zero so unsupported zero observations keep their exact result
  # and division only runs for rows with a positive denominator
  corrected$corrected_count <- 0
  supported_reference <- corrected[[denominator_column]] > 0
  corrected_values <- (
    corrected$count[supported_reference] / corrected[[denominator_column]][supported_reference]
  )
  # A positive denominator may still be too small for a finite floating-point
  # quotient. Reject infinities here before frequency normalization
  if (any(!is.finite(corrected_values))) {
    non_finite_labels <- sort(unique(
      corrected$motif[supported_reference][!is.finite(corrected_values)]
    ))
    stop(
      "Reference correction produced non-finite corrected counts for motifs: ",
      paste(non_finite_labels, collapse = ", "),
      call. = FALSE
    )
  }
  corrected$corrected_count[supported_reference] <- corrected_values
  # Only positive unsupported observations are unknown. Unsupported zero counts
  # remain zero. Frequency normalization later expands any resulting NA to its row
  if (identical(unsupported_motifs, "keep_na")) {
    corrected$corrected_count[positive_unsupported_reference] <- NA_real_
  }
  corrected
}

#' Normalize corrected counts within each selected sample row.
#'
#' Counts are scaled by the largest corrected count in each row before summing,
#' which preserves the normalized result without overflowing a finite total. A
#' zero-total row remains all zero. Under `"keep_na"`, any unsupported positive
#' count makes the entire row's corrected frequencies `NA`, because its
#' normalized composition is not defined.
#'
#' @return Corrected rows with `corrected_frequency`.
#' @noRd
cf_add_corrected_frequency <- function(corrected, row_columns) {
  corrected$corrected_frequency <- 0
  if (nrow(corrected) == 0L) {
    return(corrected)
  }

  # Normalize a row as a unit because missing counts affect the complete row and
  # scaling must happen before its counts are summed
  normalize_row <- function(row_counts) {
    row_has_unknown_count <- anyNA(row_counts)
    row_counts[is.na(row_counts)] <- 0
    row_maximum <- max(row_counts)
    if (row_maximum <= 0) {
      frequencies <- numeric(length(row_counts))
    } else {
      scaled_counts <- row_counts / row_maximum
      frequencies <- scaled_counts / sum(scaled_counts)
    }
    if (row_has_unknown_count) {
      frequencies[] <- NA_real_
    }
    frequencies
  }

  # Row-identifying metadata determines which counts share a normalization total
  row_keys <- cf_reference_correction_row_keys(corrected, row_columns)
  row_groups <- match(row_keys, unique(row_keys))
  corrected$corrected_frequency <- stats::ave(
    corrected$corrected_count,
    row_groups,
    FUN = normalize_row
  )
  corrected
}

#' Assign row-order indices from first occurrence of each selected row key.
#'
#' Repeated motif rows receive the same index. Unique row keys retain their
#' first-occurrence order, allowing joins and aggregation to restore sample-row
#' order independently of motif order.
#'
#' @return Integer row order vector.
#' @noRd
cf_reference_correction_row_order <- function(data_frame, row_columns) {
  row_keys <- cf_reference_correction_row_keys(data_frame, row_columns)
  match(row_keys, unique(row_keys))
}

#' Build row-position lookup for corrected sparse matrices.
#'
#' @param ends End-motif object.
#' @param row_indices Selected one-based row indices.
#'
#' @return Named integer vector mapping row keys to selected row positions.
#' @noRd
cf_reference_correction_row_positions <- function(ends, row_indices) {
  row_metadata <- ends$row_metadata[row_indices, , drop = FALSE]
  row_keys <- cf_reference_correction_row_keys(
    row_metadata,
    cf_reference_correction_row_columns(ends$row_mode)
  )
  stats::setNames(seq_along(row_indices), row_keys)
}

#' Build row keys from row metadata.
#'
#' @param data_frame Data frame with row-key columns.
#' @param row_columns Row-key columns or key names for lookup.
#'
#' @return Character row-key vector.
#' @noRd
cf_reference_correction_row_keys <- function(data_frame, row_columns) {
  missing_columns <- setdiff(row_columns, names(data_frame))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing row-key columns for reference correction: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(row_columns) == 1L && row_columns[[1L]] %in% names(data_frame)) {
    return(as.character(data_frame[[row_columns[[1L]]]]))
  }
  do.call(
    paste,
    c(data_frame[row_columns], sep = "\r")
  )
}

#' Resolve end-motif row selectors for reference correction.
#'
#' @param ends End-motif object.
#' @param window_idxs Optional one-based window indices.
#' @param groups Optional group names.
#' @param group_idxs Optional one-based group indices.
#'
#' @return One-based row indices for `cf_end_motif_data_frame()`.
#' @noRd
cf_reference_correction_end_row_indices <- function(ends, window_idxs, groups, group_idxs) {
  if (identical(ends$row_mode, "bed") || identical(ends$row_mode, "size")) {
    if (!is.null(groups) || !is.null(group_idxs)) {
      stop("Grouped selectors can only be used with grouped output", call. = FALSE)
    }
    return(cf_resolve_end_motif_window_indices(ends, window_idxs))
  }
  if (identical(ends$row_mode, "grouped_bed")) {
    if (!is.null(window_idxs)) {
      stop("window_idxs can only be used with windowed output", call. = FALSE)
    }
    return(cf_resolve_end_motif_group_indices(ends, groups, group_idxs))
  }
  if (!is.null(window_idxs) || !is.null(groups) || !is.null(group_idxs)) {
    stop("Row selectors cannot be used with global output", call. = FALSE)
  }
  seq_len(length(ends$row_idx0))
}

#' Resolve reference row indices from selected end-motif row keys.
#'
#' @param ref_kmers Reference k-mer object.
#' @param end_rows Selected end-motif rows.
#' @param reference_row_columns Reference row-key columns.
#'
#' @return One-based row indices for `cf_ref_kmer_data_frame()`.
#' @noRd
cf_reference_correction_ref_row_indices_from_end_rows <- function(
  ref_kmers,
  end_rows,
  reference_row_columns
) {
  if (length(reference_row_columns) == 0L) {
    return(seq_len(length(ref_kmers$row_idx0)))
  }

  reference_keys <- cf_reference_correction_row_keys(
    ref_kmers$row_metadata,
    reference_row_columns
  )
  selected_keys <- unique(cf_reference_correction_row_keys(end_rows, reference_row_columns))
  matched_rows <- match(selected_keys, reference_keys)
  if (anyNA(matched_rows)) {
    stop("Selected end-motif row has no matching reference k-mer row", call. = FALSE)
  }
  as.integer(matched_rows)
}

#' Validate motif axes for reference correction.
#'
#' @param ends End-motif object.
#' @param ref_kmers Reference k-mer object.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_reference_correction_motif_axes <- function(ends, ref_kmers) {
  if (isTRUE(ref_kmers$canonical)) {
    stop("Reference correction requires non-canonical reference k-mer output", call. = FALSE)
  }
  if (identical(ends$motif_axis_kind, "motif_group")) {
    if (!identical(ref_kmers$motif_axis_kind, "motif_group")) {
      stop("Grouped end-motif output requires grouped reference k-mer output", call. = FALSE)
    }
    return(invisible(TRUE))
  }

  if (!identical(ref_kmers$motif_axis_kind, "motif")) {
    stop("End-motif output with motif labels requires reference k-mer output with motif labels", call. = FALSE)
  }

  reference_motif_widths <- nchar(gsub("_", "", ends$motif, fixed = TRUE))
  if (any(reference_motif_widths != ref_kmers$kmer_size)) {
    stop(
      "End-motif width must match reference k-mer size (",
      ref_kmers$kmer_size,
      ")",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Return row columns used to match end-motif and reference k-mer rows.
#'
#' @param row_mode Output row mode.
#'
#' @return Character vector of row-key columns.
#' @noRd
cf_reference_correction_row_columns <- function(row_mode) {
  if (identical(row_mode, "global")) {
    return("row_label")
  }
  if (row_mode %in% c("size", "bed")) {
    return(c("window_idx", "chrom", "start", "end"))
  }
  if (identical(row_mode, "grouped_bed")) {
    return("group_name")
  }
  stop("Unsupported end-motif row mode for correction: ", row_mode, call. = FALSE)
}

#' Return reference row columns used for correction.
#'
#' @param ends End-motif object.
#' @param ref_kmers Reference k-mer object.
#' @param use_global_bias Whether a global reference can be broadcast.
#'
#' @return Character vector of reference row-key columns.
#' @noRd
cf_reference_correction_reference_row_columns <- function(ends, ref_kmers, use_global_bias) {
  if (
    isTRUE(use_global_bias) &&
      identical(ref_kmers$row_mode, "global") &&
      !identical(ends$row_mode, "global")
  ) {
    return(character())
  }
  cf_reference_correction_row_columns(ref_kmers$row_mode)
}

#' Validate row compatibility before reference correction.
#'
#' @param ends End-motif object.
#' @param ref_kmers Reference k-mer object.
#' @param row_columns End-motif row-key columns.
#' @param reference_row_columns Reference row-key columns.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_reference_correction_rows <- function(
  ends,
  ref_kmers,
  row_columns,
  reference_row_columns
) {
  if (length(reference_row_columns) == 0L) {
    return(invisible(TRUE))
  }

  end_row_keys <- cf_sorted_unique_rows(ends$row_metadata[row_columns])
  ref_row_keys <- cf_sorted_unique_rows(ref_kmers$row_metadata[reference_row_columns])
  if (nrow(end_row_keys) != nrow(ends$row_metadata)) {
    stop("End-motif row labels are not unique enough for correction", call. = FALSE)
  }
  if (nrow(ref_row_keys) != nrow(ref_kmers$row_metadata)) {
    stop("Reference k-mer row labels are not unique enough for correction", call. = FALSE)
  }
  if (!identical(end_row_keys, ref_row_keys)) {
    stop(
      "End-motif and reference k-mer rows do not match. Run ref-kmers with ",
      "the same windowing or grouping.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Sort and deduplicate row-key data frames.
#'
#' @param data_frame Data frame.
#'
#' @return Sorted unique data frame.
#' @noRd
cf_sorted_unique_rows <- function(data_frame) {
  data_frame <- unique(data_frame)
  data_frame <- data_frame[do.call(order, data_frame), , drop = FALSE]
  row.names(data_frame) <- NULL
  data_frame
}

#' Add reference-correction columns to an empty data frame.
#'
#' @param data_frame Empty end-motif data frame.
#'
#' @return Empty data frame with correction columns.
#' @noRd
cf_add_empty_reference_correction_columns <- function(data_frame) {
  data_frame$corrected_count <- numeric()
  data_frame$corrected_frequency <- numeric()
  data_frame
}
