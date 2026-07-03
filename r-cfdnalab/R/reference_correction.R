#' Correct end-motif counts for reference k-mer composition.
#'
#' This helper starts from `end_motif_data_frame()`. Each count is divided by
#' `reference_frequency * correction_motif_count`.
#' `correction_motif_count` is computed separately for each reference row from
#' motifs with positive reference frequency. A uniform reference within that row
#' leaves counts unchanged.
#'
#' Concrete end-motif labels are matched to reference k-mers by removing the
#' `_` separator, for example `AT_CG -> ATCG`. Motif-group outputs are matched
#' directly by group label.
#'
#' Reference k-mer output is read without densifying. For sparse reference
#' output, omitted row/motif pairs are treated as zero frequency.
#'
#' Positive end-motif counts with zero or missing reference frequency cannot be
#' divided by a reference bias. By default this is an error. Set
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
#' rows. They do not change the reference support used for scaling, so
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
#' @param unsupported_motifs What to do with positive end-motif counts that have
#'   no positive reference frequency. Use `"error"`, `"drop"`, or `"keep_na"`.
#'
#' @return An end-motif data frame with `reference_motif`,
#'   `reference_frequency`, `correction_motif_count`, `reference_scale`, and
#'   `reference_corrected_count`.
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
  unsupported_motifs = "error"
) {
  unsupported_motifs <- cf_match_exact_argument(
    unsupported_motifs,
    c("error", "drop", "keep_na"),
    "unsupported_motifs"
  )
  if (!inherits(ends, "cfdnalab_end_motif_counts")) {
    stop("ends must be a cfDNAlab end-motif object", call. = FALSE)
  }
  if (!inherits(ref_kmers, "cfdnalab_ref_kmer_frequencies")) {
    stop("ref_kmers must be a cfDNAlab reference k-mer object", call. = FALSE)
  }
  cf_validate_scalar_logical(densify, "densify")
  cf_validate_scalar_logical(use_global_bias, "use_global_bias")
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

  row_columns <- cf_reference_correction_row_columns(ends$row_mode)
  reference_row_columns <- cf_reference_correction_reference_row_columns(
    ends,
    ref_kmers,
    use_global_bias
  )
  cf_validate_reference_correction_rows(
    ends,
    ref_kmers,
    row_columns,
    reference_row_columns
  )

  end_row_indices <- cf_reference_correction_end_row_indices(
    ends,
    window_idxs,
    groups,
    group_idxs
  )
  ref_row_indices <- cf_reference_correction_ref_row_indices(
    ref_kmers,
    window_idxs,
    groups,
    group_idxs,
    use_global_bias
  )
  end_motif_indices <- cf_resolve_end_motif_indices(ends, motifs, motif_idxs)

  end_rows <- cf_end_motif_data_frame(
    ends,
    row_indices = end_row_indices,
    motif_indices = end_motif_indices,
    densify = densify,
    max_blacklisted_fraction = max_blacklisted_fraction
  )
  if (nrow(end_rows) == 0L) {
    return(cf_add_empty_reference_correction_columns(end_rows))
  }

  ref_rows <- cf_ref_kmer_data_frame(
    ref_kmers,
    row_indices = ref_row_indices,
    motif_indices = seq_along(ref_kmers$motif_idx0),
    densify = FALSE,
    max_blacklisted_fraction = 1.0
  )

  if (identical(ends$motif_axis_kind, "motif_group")) {
    end_rows$reference_motif <- end_rows$motif
  } else {
    end_rows$reference_motif <- gsub("_", "", end_rows$motif, fixed = TRUE)
  }
  end_column_names <- names(end_rows)
  names(ref_rows)[names(ref_rows) == "motif"] <- "reference_motif"
  names(ref_rows)[names(ref_rows) == "frequency"] <- "reference_frequency"
  ref_rows <- ref_rows[c(reference_row_columns, "reference_motif", "reference_frequency")]
  ref_rows <- cf_reference_correction_filter_ref_rows(
    ref_rows,
    end_rows,
    reference_row_columns
  )

  merge_columns <- c(reference_row_columns, "reference_motif")
  if (any(duplicated(ref_rows[merge_columns]))) {
    stop("Reference k-mer rows are not unique for row and motif labels", call. = FALSE)
  }
  reference_support_counts <- cf_reference_correction_support_counts(
    ref_rows,
    reference_row_columns
  )

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
  row.names(corrected) <- NULL

  missing_reference <- is.na(corrected$reference_frequency)
  if (any(missing_reference)) {
    corrected$reference_frequency[missing_reference] <- 0
  }
  corrected <- cf_add_reference_correction_motif_count(
    corrected,
    reference_support_counts,
    reference_row_columns
  )

  unsupported_reference <- corrected$reference_frequency <= 0
  positive_unsupported_reference <- unsupported_reference & corrected$count > 0
  if (any(positive_unsupported_reference) && identical(unsupported_motifs, "error")) {
    unsupported_labels <- sort(unique(corrected$reference_motif[positive_unsupported_reference]))
    stop(
      "Positive-count end motifs have no positive reference frequency: ",
      paste(unsupported_labels, collapse = ", "),
      ". Pass unsupported_motifs = \"drop\" to omit those rows, or ",
      "unsupported_motifs = \"keep_na\" to keep them with NA corrected counts.",
      call. = FALSE
    )
  }
  if (identical(unsupported_motifs, "drop")) {
    corrected <- corrected[!unsupported_reference, , drop = FALSE]
    positive_unsupported_reference <- positive_unsupported_reference[!unsupported_reference]
  }

  corrected$reference_scale <- corrected$reference_frequency * corrected$correction_motif_count
  corrected$reference_corrected_count <- 0
  supported_reference <- corrected$reference_scale > 0
  corrected$reference_corrected_count[supported_reference] <- (
    corrected$count[supported_reference] / corrected$reference_scale[supported_reference]
  )
  if (identical(unsupported_motifs, "keep_na")) {
    corrected$reference_corrected_count[positive_unsupported_reference] <- NA_real_
  }
  corrected[
    c(
      end_column_names,
      "reference_frequency",
      "correction_motif_count",
      "reference_scale",
      "reference_corrected_count"
    )
  ]
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
  unsupported_motifs = "error"
) {
  cf_validate_fixed_shape_reference_correction_policy(
    unsupported_motifs,
    "dense_corrected_counts_matrix"
  )
  cf_validate_scalar_logical(allow_densify, "allow_densify")
  if (identical(ends$storage_mode, "sparse_coo") && !isTRUE(allow_densify)) {
    stop(
      "This end-motif store is sparse. Use sparse_corrected_counts_matrix() ",
      "or set allow_densify = TRUE.",
      call. = FALSE
    )
  }
  row_indices <- cf_reference_correction_end_row_indices(ends, window_idxs, groups, group_idxs)
  row_indices <- cf_apply_end_motif_blacklist_filter(ends, row_indices, max_blacklisted_fraction)
  motif_indices <- cf_resolve_end_motif_indices(ends, motifs, motif_idxs)
  corrected <- cf_reference_corrected_end_motif_data_frame(
    ends,
    ref_kmers,
    window_idxs = window_idxs,
    groups = groups,
    group_idxs = group_idxs,
    densify = TRUE,
    motifs = motifs,
    motif_idxs = motif_idxs,
    max_blacklisted_fraction = max_blacklisted_fraction,
    use_global_bias = use_global_bias,
    unsupported_motifs = unsupported_motifs
  )
  matrix(
    corrected$reference_corrected_count,
    nrow = length(row_indices),
    ncol = length(motif_indices),
    byrow = TRUE
  )
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
  unsupported_motifs = "error"
) {
  cf_validate_fixed_shape_reference_correction_policy(
    unsupported_motifs,
    "sparse_corrected_counts_matrix"
  )
  row_indices <- cf_reference_correction_end_row_indices(ends, window_idxs, groups, group_idxs)
  row_indices <- cf_apply_end_motif_blacklist_filter(ends, row_indices, max_blacklisted_fraction)
  motif_indices <- cf_resolve_end_motif_indices(ends, motifs, motif_idxs)
  if (length(row_indices) == 0L || length(motif_indices) == 0L) {
    return(Matrix::sparseMatrix(
      i = integer(),
      j = integer(),
      dims = as.integer(c(length(row_indices), length(motif_indices)))
    ))
  }
  corrected <- cf_reference_corrected_end_motif_data_frame(
    ends,
    ref_kmers,
    window_idxs = window_idxs,
    groups = groups,
    group_idxs = group_idxs,
    densify = FALSE,
    motifs = motifs,
    motif_idxs = motif_idxs,
    max_blacklisted_fraction = max_blacklisted_fraction,
    use_global_bias = use_global_bias,
    unsupported_motifs = unsupported_motifs
  )
  if (nrow(corrected) == 0L) {
    return(Matrix::sparseMatrix(
      i = integer(),
      j = integer(),
      dims = as.integer(c(length(row_indices), length(motif_indices)))
    ))
  }
  row_positions <- cf_reference_correction_row_positions(ends, row_indices)
  motif_positions <- stats::setNames(seq_along(motif_indices), ends$motif[motif_indices])
  corrected_values <- corrected$reference_corrected_count
  stored <- corrected_values != 0 | is.na(corrected_values)
  if (!any(stored)) {
    return(Matrix::sparseMatrix(
      i = integer(),
      j = integer(),
      dims = as.integer(c(length(row_indices), length(motif_indices)))
    ))
  }
  Matrix::sparseMatrix(
    i = unname(row_positions[cf_reference_correction_row_keys(
      corrected,
      cf_reference_correction_row_columns(ends$row_mode)
    )])[stored],
    j = unname(motif_positions[corrected$motif])[stored],
    x = corrected_values[stored],
    dims = as.integer(c(length(row_indices), length(motif_indices)))
  )
}

#' Validate unsupported motif policy for fixed-shape outputs.
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

#' Resolve reference row selectors for reference correction.
#'
#' @param ref_kmers Reference k-mer object.
#' @param window_idxs Optional one-based window indices.
#' @param groups Optional group names.
#' @param group_idxs Optional one-based group indices.
#' @param use_global_bias Whether a global reference may be applied to every
#'   end-motif row.
#'
#' @return One-based row indices for `cf_ref_kmer_data_frame()`.
#' @noRd
cf_reference_correction_ref_row_indices <- function(
  ref_kmers,
  window_idxs,
  groups,
  group_idxs,
  use_global_bias
) {
  if (isTRUE(use_global_bias) && identical(ref_kmers$row_mode, "global")) {
    return(seq_len(length(ref_kmers$row_idx0)))
  }
  if (identical(ref_kmers$row_mode, "bed") || identical(ref_kmers$row_mode, "size")) {
    if (!is.null(groups) || !is.null(group_idxs)) {
      stop("Grouped selectors can only be used with grouped output", call. = FALSE)
    }
    return(cf_resolve_ref_kmer_window_indices(ref_kmers, window_idxs))
  }
  if (identical(ref_kmers$row_mode, "grouped_bed")) {
    if (!is.null(window_idxs)) {
      stop("window_idxs can only be used with windowed output", call. = FALSE)
    }
    return(cf_resolve_ref_kmer_group_indices(ref_kmers, groups, group_idxs))
  }
  if (!is.null(window_idxs) || !is.null(groups) || !is.null(group_idxs)) {
    stop("Row selectors cannot be used with global output", call. = FALSE)
  }
  seq_len(length(ref_kmers$row_idx0))
}

#' Validate motif axes for reference correction.
#'
#' @param ends End-motif object.
#' @param ref_kmers Reference k-mer object.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_reference_correction_motif_axes <- function(ends, ref_kmers) {
  if (identical(ends$motif_axis_kind, "motif_group")) {
    if (!identical(ref_kmers$motif_axis_kind, "motif_group")) {
      stop("Grouped end-motif output requires grouped reference k-mer output", call. = FALSE)
    }
    return(invisible(TRUE))
  }

  if (!identical(ref_kmers$motif_axis_kind, "motif")) {
    stop("Concrete end-motif output requires concrete reference k-mer output", call. = FALSE)
  }
  if (isTRUE(ref_kmers$canonical)) {
    stop("Reference correction requires non-canonical reference k-mer output", call. = FALSE)
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

#' Count positive reference-frequency motifs per correction row.
#'
#' @param ref_rows Reference k-mer rows with correction column names.
#' @param reference_row_columns Reference row-key columns.
#'
#' @return A scalar for broadcast global references, otherwise a data frame.
#' @noRd
cf_reference_correction_support_counts <- function(ref_rows, reference_row_columns) {
  positive_ref_rows <- ref_rows[ref_rows$reference_frequency > 0, , drop = FALSE]
  if (length(reference_row_columns) == 0L) {
    return(nrow(positive_ref_rows))
  }
  if (nrow(positive_ref_rows) == 0L) {
    return(data.frame(
      positive_ref_rows[reference_row_columns],
      correction_motif_count = integer(),
      stringsAsFactors = FALSE
    ))
  }
  counts <- stats::aggregate(
    rep(1L, nrow(positive_ref_rows)),
    by = positive_ref_rows[reference_row_columns],
    FUN = sum
  )
  names(counts)[ncol(counts)] <- "correction_motif_count"
  counts
}

#' Keep reference rows whose row keys are present in selected end rows.
#'
#' @param ref_rows Reference k-mer rows with correction column names.
#' @param end_rows Selected end-motif rows.
#' @param reference_row_columns Reference row-key columns.
#'
#' @return Filtered `ref_rows`.
#' @noRd
cf_reference_correction_filter_ref_rows <- function(
  ref_rows,
  end_rows,
  reference_row_columns
) {
  if (length(reference_row_columns) == 0L) {
    return(ref_rows)
  }
  selected_row_keys <- unique(
    cf_reference_correction_row_keys(end_rows, reference_row_columns)
  )
  ref_row_keys <- cf_reference_correction_row_keys(ref_rows, reference_row_columns)
  ref_rows[ref_row_keys %in% selected_row_keys, , drop = FALSE]
}

#' Add per-row correction motif counts to corrected end-motif rows.
#'
#' @param corrected End-motif rows after joining reference frequencies.
#' @param reference_support_counts Counts from
#'   `cf_reference_correction_support_counts()`.
#' @param reference_row_columns Reference row-key columns.
#'
#' @return `corrected` with a `correction_motif_count` column.
#' @noRd
cf_add_reference_correction_motif_count <- function(
  corrected,
  reference_support_counts,
  reference_row_columns
) {
  if (length(reference_row_columns) == 0L) {
    corrected$correction_motif_count <- as.integer(reference_support_counts)
    return(corrected)
  }

  corrected$.cfdnalab_order <- seq_len(nrow(corrected))
  corrected <- merge(
    corrected,
    reference_support_counts,
    by = reference_row_columns,
    all.x = TRUE,
    sort = FALSE
  )
  corrected <- corrected[order(corrected$.cfdnalab_order), , drop = FALSE]
  corrected$.cfdnalab_order <- NULL
  row.names(corrected) <- NULL
  corrected$correction_motif_count[is.na(corrected$correction_motif_count)] <- 0L
  corrected$correction_motif_count <- as.integer(corrected$correction_motif_count)
  corrected
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
  data_frame$reference_motif <- character()
  data_frame$reference_frequency <- numeric()
  data_frame$correction_motif_count <- integer()
  data_frame$reference_scale <- numeric()
  data_frame$reference_corrected_count <- numeric()
  data_frame
}
