#' Supported reference k-mer storage modes.
#'
#' These values mirror the cfDNAlab Zarr `storage_mode` root attribute.
#'
#' @noRd
REF_KMER_VALID_STORAGE_MODES <- c("dense", "sparse_coo")

#' Supported reference k-mer row modes.
#'
#' These values mirror the cfDNAlab Zarr `row_mode` root attribute.
#'
#' @noRd
REF_KMER_VALID_ROW_MODES <- c("global", "size", "bed", "grouped_bed")

#' Supported reference k-mer motif-axis kinds.
#'
#' @noRd
REF_KMER_VALID_AXIS_KINDS <- c("motif", "motif_group")

#' Read cfDNAlab reference k-mer frequencies.
#'
#' Loads a `<prefix>.ref_kmer_counts.zarr` output directory created with the
#' \code{cfdna ref-kmers} CLI tool from the main \code{cfDNAlab} rust package.
#' The directory is a Zarr store on disk, but ordinary workflows can use the
#' data frame and matrix helpers without working with Zarr directly.
#'
#' Reference k-mer outputs store frequencies. Count helpers reconstruct counts
#' by multiplying each frequency row by its `row_scaling_factor`. A row can
#' describe the whole reference, a genomic window, a BED interval, or a grouped
#' BED entry depending on how the command was run.
#'
#' The loader validates the cfDNAlab schema, row metadata, motif metadata,
#' frequency layout, and the metadata needed to reconstruct counts.
#'
#' @param path Path to a cfDNAlab reference k-mer `.zarr` directory.
#'
#' @return One of `cfdnalab_global_ref_kmer_frequencies`,
#'   `cfdnalab_windowed_ref_kmer_frequencies`, or
#'   `cfdnalab_grouped_ref_kmer_frequencies`, depending on the row mode.
#' @export
#'
#' @examples
#' \dontrun{
#' ref_kmers <- read_ref_kmers("sample.ref_kmer_counts.zarr")
#' motifs(ref_kmers)
#' sparse_frequencies_matrix(ref_kmers)
#' }
read_ref_kmers <- function(path) {
  path <- cf_validate_zarr_store_path(path, "Reference k-mer")
  root_attributes <- cf_root_attributes(path)
  cf_validate_schema(root_attributes, "ref_kmer_frequencies", "reference k-mer")

  storage_mode <- cf_validate_allowed_string(
    root_attributes$storage_mode,
    REF_KMER_VALID_STORAGE_MODES,
    "reference k-mer storage mode"
  )
  row_mode <- cf_validate_allowed_string(
    root_attributes$row_mode,
    REF_KMER_VALID_ROW_MODES,
    "reference k-mer row mode"
  )
  motif_axis_kind <- cf_validate_allowed_string(
    root_attributes$motif_axis_kind,
    REF_KMER_VALID_AXIS_KINDS,
    "reference k-mer motif axis kind"
  )
  ref_kmer_metadata <- cf_validate_ref_kmer_root_attributes(root_attributes, storage_mode)

  store <- cf_open_zarr(path, "reference k-mer")
  cf_required_arrays(
    store,
    cf_ref_kmer_required_arrays(storage_mode, row_mode, motif_axis_kind),
    "Reference k-mer"
  )
  cf_validate_dimension_names(path, "motif_index", "motif")
  cf_validate_dimension_names(path, "row", "row")

  motif_axis <- cf_read_vector(store, "motif_index", "Reference k-mer")
  row <- cf_read_vector(store, "row", "Reference k-mer")
  motif <- NULL
  if (identical(motif_axis_kind, "motif")) {
    cf_validate_dimension_names(path, "motif_byte", "motif_byte")
    cf_validate_dimension_names(path, "motif_ascii", c("motif", "motif_byte"))
    motif_byte <- cf_read_vector(store, "motif_byte", "Reference k-mer")
    motif_ascii <- cf_read_array(store, "motif_ascii", "Reference k-mer")
    motif <- cf_decode_motif_ascii(motif_ascii, length(motif_axis), length(motif_byte))
    cf_validate_axis(motif_byte, "motif_byte")
    cf_validate_ref_kmer_motifs(motif, ref_kmer_metadata$kmer_size, ref_kmer_metadata$canonical)
  } else {
    motif <- cf_read_labels(path, "motif_index", "motif_group", length(motif_axis))
    cf_validate_unique_axis_labels(motif, "reference k-mer motif-group")
  }

  cf_validate_axis(motif_axis, "motif_index")
  cf_validate_axis(row, "row")
  if (identical(row_mode, "global") && length(row) != 1L) {
    stop("global reference k-mer output must contain exactly one row", call. = FALSE)
  }

  cf_validate_dimension_names(path, "row_scaling_factor", "row")
  row_scaling_factor <- cf_read_vector(store, "row_scaling_factor", "Reference k-mer")
  cf_validate_same_length(row_scaling_factor, row, "row_scaling_factor", "row")
  cf_validate_nonnegative_numeric_vector(row_scaling_factor, "row_scaling_factor")

  reference_contig_footprint <- cf_read_ref_kmer_reference_footprint(path, store)

  frequencies <- NULL
  sparse <- NULL
  if (identical(storage_mode, "dense")) {
    cf_validate_dimension_names(path, "frequencies", c("row", "motif"))
    frequencies <- cf_get_array(store, "frequencies", "Reference k-mer")
    expected_shape <- c(length(row), length(motif_axis))
    if (!identical(cf_array_shape(frequencies), expected_shape)) {
      stop(
        "dense frequencies shape does not match row and motif axes: frequencies=",
        paste(cf_array_shape(frequencies), collapse = "x"),
        ", coordinates=",
        paste(expected_shape, collapse = "x"),
        call. = FALSE
      )
    }
  } else {
    sparse <- list(
      row_idx0 = cf_read_vector(store, "sparse/row", "Reference k-mer"),
      motif_idx0 = cf_read_vector(store, "sparse/motif", "Reference k-mer"),
      frequency = cf_read_vector(store, "sparse/frequency", "Reference k-mer"),
      shape = cf_read_vector(store, "sparse/shape", "Reference k-mer")
    )
    cf_validate_dimension_names(path, "sparse/row", "nnz")
    cf_validate_dimension_names(path, "sparse/motif", "nnz")
    cf_validate_dimension_names(path, "sparse/frequency", "nnz")
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
    cf_validate_same_length(sparse$frequency, sparse$row_idx0, "sparse/frequency", "sparse/row")
    cf_validate_index_vector(sparse$row_idx0, length(row), "sparse/row")
    cf_validate_index_vector(sparse$motif_idx0, length(motif_axis), "sparse/motif")
    cf_validate_ref_kmer_sparse_coordinates(sparse$row_idx0, sparse$motif_idx0)
    cf_validate_ref_kmer_frequency_vector(sparse$frequency, "sparse/frequency")
  }

  row_metadata <- cf_read_ref_kmer_row_metadata(
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
    kmer_size = ref_kmer_metadata$kmer_size,
    canonical = ref_kmer_metadata$canonical,
    all_motifs = ref_kmer_metadata$all_motifs,
    assign_by = ref_kmer_metadata$assign_by,
    motif_idx0 = as.integer(motif_axis),
    motif = motif,
    row_idx0 = as.integer(row),
    row_scaling_factor = as.numeric(row_scaling_factor),
    reference_contig_footprint = reference_contig_footprint,
    frequencies = frequencies,
    sparse = sparse,
    row_metadata = row_metadata
  )

  class(object) <- c(
    switch(
      row_mode,
      global = "cfdnalab_global_ref_kmer_frequencies",
      size = "cfdnalab_windowed_ref_kmer_frequencies",
      bed = "cfdnalab_windowed_ref_kmer_frequencies",
      grouped_bed = "cfdnalab_grouped_ref_kmer_frequencies"
    ),
    "cfdnalab_ref_kmer_frequencies",
    "cfdnalab_zarr_store"
  )
  object
}

#' Validate reference k-mer root attributes that define interpretation.
#'
#' @param attrs Root Zarr attributes.
#' @param storage_mode Validated storage mode.
#'
#' @return A list with validated metadata values.
#' @noRd
cf_validate_ref_kmer_root_attributes <- function(attrs, storage_mode) {
  cf_require_root_attribute(attrs, "value_units", "reference_kmer_frequency")
  cf_require_root_attribute(attrs, "count_units", "reference_kmer_count")
  cf_require_root_attribute(attrs, "row_scaling_factor_array", "row_scaling_factor")
  cf_require_root_attribute(
    attrs,
    "count_reconstruction",
    "reference_kmer_count = frequency * row_scaling_factor[row]"
  )
  if (identical(storage_mode, "dense")) {
    cf_require_root_attribute(attrs, "primary_array", "frequencies")
  } else {
    cf_require_root_attribute(attrs, "primary_group", "sparse")
    cf_require_root_attribute(attrs, "sparse_format", "coo")
    sparse_indices_base <- attrs$sparse_indices_base
    if (
      length(sparse_indices_base) != 1L ||
        !is.numeric(sparse_indices_base) ||
        is.na(sparse_indices_base) ||
        !is.finite(sparse_indices_base) ||
        sparse_indices_base != 0
    ) {
      stop("sparse_indices_base must be 0 for reference k-mer sparse COO output", call. = FALSE)
    }
  }

  kmer_size <- attrs$kmer_size
  if (
    length(kmer_size) != 1L ||
      !is.numeric(kmer_size) ||
      is.na(kmer_size) ||
      !is.finite(kmer_size) ||
      kmer_size != as.integer(kmer_size) ||
      kmer_size < 1L
  ) {
    stop("kmer_size must be a positive integer", call. = FALSE)
  }
  if (length(attrs$canonical) != 1L || !is.logical(attrs$canonical) || is.na(attrs$canonical)) {
    stop("canonical must be TRUE or FALSE", call. = FALSE)
  }
  if (length(attrs$all_motifs) != 1L || !is.logical(attrs$all_motifs) || is.na(attrs$all_motifs)) {
    stop("all_motifs must be TRUE or FALSE", call. = FALSE)
  }
  assign_by <- cf_validate_scalar_string(attrs$assign_by, "assign_by")
  if (!nzchar(assign_by)) {
    stop("assign_by must be a non-empty string", call. = FALSE)
  }

  list(
    kmer_size = as.integer(kmer_size),
    canonical = isTRUE(attrs$canonical),
    all_motifs = isTRUE(attrs$all_motifs),
    assign_by = assign_by
  )
}

#' Require a root attribute to equal an expected scalar value.
#'
#' @param attrs Root Zarr attributes.
#' @param name Attribute name.
#' @param expected Expected value.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_require_root_attribute <- function(attrs, name, expected) {
  if (!identical(attrs[[name]], expected)) {
    stop(
      "Reference k-mer root attribute ",
      sQuote(name),
      " must be ",
      sQuote(expected),
      ", found ",
      deparse(attrs[[name]]),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Return arrays required by a reference k-mer store.
#'
#' @param storage_mode Reference k-mer storage mode.
#' @param row_mode Reference k-mer row mode.
#' @param motif_axis_kind Reference k-mer motif-axis kind.
#'
#' @return Character vector of required array paths.
#' @noRd
cf_ref_kmer_required_arrays <- function(storage_mode, row_mode, motif_axis_kind) {
  required <- c("motif_index", "row", "row_scaling_factor", "reference_contig_footprint_json")
  if (identical(motif_axis_kind, "motif")) {
    required <- c(required, "motif_byte", "motif_ascii")
  }
  if (identical(storage_mode, "dense")) {
    required <- c(required, "frequencies")
  } else {
    required <- c(
      required,
      "sparse/row",
      "sparse/motif",
      "sparse/frequency",
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

#' Read reference k-mer row metadata.
#'
#' @param path Path to a Zarr store.
#' @param store Open Zarr store.
#' @param row Row-axis values.
#' @param row_mode Reference k-mer row mode.
#'
#' @return A data frame describing frequency rows.
#' @noRd
cf_read_ref_kmer_row_metadata <- function(path, store, row, row_mode) {
  if (identical(row_mode, "global")) {
    labels <- cf_read_labels(path, "row", "row_label", length(row))
    if (!identical(labels, "global")) {
      stop("global reference k-mer output must contain exactly one row labeled 'global'", call. = FALSE)
    }
    return(data.frame(row_label = labels, stringsAsFactors = FALSE))
  }

  if (row_mode %in% c("size", "bed")) {
    cf_validate_dimension_names(path, "chromosome", "chromosome")
    cf_validate_dimension_names(path, "row_chromosome", "row")
    cf_validate_dimension_names(path, "row_start_bp", "row")
    cf_validate_dimension_names(path, "row_end_bp", "row")
    cf_validate_dimension_names(path, "blacklisted_fraction", "row")
    chromosome <- cf_read_vector(store, "chromosome", "Reference k-mer")
    chromosome_name <- cf_read_labels(path, "chromosome", "chromosome_name", length(chromosome))
    row_chromosome <- cf_read_vector(store, "row_chromosome", "Reference k-mer")
    row_start_bp <- cf_read_vector(store, "row_start_bp", "Reference k-mer")
    row_end_bp <- cf_read_vector(store, "row_end_bp", "Reference k-mer")
    blacklisted_fraction <- cf_read_vector(store, "blacklisted_fraction", "Reference k-mer")
    cf_validate_axis(chromosome, "chromosome")
    cf_validate_same_length(row_chromosome, row, "row_chromosome", "row")
    cf_validate_same_length(row_start_bp, row, "row_start_bp", "row")
    cf_validate_same_length(row_end_bp, row, "row_end_bp", "row")
    cf_validate_same_length(blacklisted_fraction, row, "blacklisted_fraction", "row")
    cf_validate_index_vector(row_chromosome, length(chromosome), "row_chromosome")
    cf_validate_half_open_intervals(row_start_bp, row_end_bp, "row_start_bp", "row_end_bp")
    cf_validate_fraction_vector(blacklisted_fraction, "blacklisted_fraction")
    if (any(row_start_bp > .Machine$integer.max) || any(row_end_bp > .Machine$integer.max)) {
      stop("Reference k-mer window coordinates must fit in R integer range", call. = FALSE)
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

  group <- cf_read_vector(store, "group", "Reference k-mer")
  cf_validate_dimension_names(path, "group", "row")
  cf_validate_dimension_names(path, "eligible_windows", "row")
  cf_validate_dimension_names(path, "blacklisted_fraction", "row")
  group_name <- cf_read_labels(path, "group", "group_name", length(group))
  eligible_windows <- cf_read_vector(store, "eligible_windows", "Reference k-mer")
  blacklisted_fraction <- cf_read_vector(store, "blacklisted_fraction", "Reference k-mer")
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

#' Read the JSON reference contig footprint.
#'
#' @param path Path to a Zarr store.
#' @param store Open Zarr store.
#'
#' @return Decoded JSON value.
#' @noRd
cf_read_ref_kmer_reference_footprint <- function(path, store) {
  cf_validate_dimension_names(path, "reference_contig_footprint_json", "json_byte")
  json_bytes <- cf_read_vector(store, "reference_contig_footprint_json", "Reference k-mer")
  if (
    !is.numeric(json_bytes) ||
      any(is.na(json_bytes)) ||
      any(!is.finite(json_bytes)) ||
      any(json_bytes != as.integer(json_bytes)) ||
      any(json_bytes < 0L | json_bytes > 255L)
  ) {
    stop("reference_contig_footprint_json must contain byte values in 0..255", call. = FALSE)
  }
  json_text <- rawToChar(as.raw(as.integer(json_bytes)))
  jsonlite::fromJSON(json_text, simplifyVector = FALSE)
}

#' Validate concrete reference k-mer labels.
#'
#' @param motif Motif labels.
#' @param kmer_size Expected k-mer size.
#' @param canonical Whether labels must be canonical representatives.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_ref_kmer_motifs <- function(motif, kmer_size, canonical) {
  cf_validate_unique_axis_labels(motif, "reference k-mer motif")
  for (motif_index in seq_along(motif)) {
    label <- motif[[motif_index]]
    if (!identical(nchar(label, type = "bytes"), kmer_size)) {
      stop(
        "reference k-mer motif label ",
        sQuote(label),
        " at motif_idx ",
        motif_index,
        " has length ",
        nchar(label, type = "bytes"),
        ", expected ",
        kmer_size,
        call. = FALSE
      )
    }
    invalid_bases <- setdiff(strsplit(label, "", fixed = TRUE)[[1L]], c("A", "C", "G", "T"))
    if (length(invalid_bases) > 0L) {
      stop(
        "reference k-mer motif label ",
        sQuote(label),
        " at motif_idx ",
        motif_index,
        " contains invalid base(s): ",
        paste(invalid_bases, collapse = ", "),
        call. = FALSE
      )
    }
    if (isTRUE(canonical)) {
      canonical_label <- cf_canonical_ref_kmer(label)
      if (!identical(label, canonical_label)) {
        stop(
          "canonical reference k-mer motif label ",
          sQuote(label),
          " at motif_idx ",
          motif_index,
          " should be ",
          sQuote(canonical_label),
          call. = FALSE
        )
      }
    }
  }
  invisible(TRUE)
}

#' Validate that public axis labels are unique.
#'
#' @param labels Axis labels.
#' @param label_name Human-readable label type.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_unique_axis_labels <- function(labels, label_name) {
  duplicate_index <- anyDuplicated(labels)
  if (duplicate_index > 0L) {
    stop("duplicate ", label_name, " label: ", sQuote(labels[[duplicate_index]]), call. = FALSE)
  }
  invisible(TRUE)
}

#' Return the canonical reference k-mer label.
#'
#' @param motif Concrete A/C/G/T motif.
#'
#' @return Canonical motif label.
#' @noRd
cf_canonical_ref_kmer <- function(motif) {
  reverse_complement <- cf_reverse_complement(motif)
  if (nchar(motif, type = "bytes") %% 2L == 0L) {
    return(min(motif, reverse_complement))
  }
  bases <- strsplit(motif, "", fixed = TRUE)[[1L]]
  middle_base <- bases[[ceiling(length(bases) / 2L)]]
  if (middle_base %in% c("A", "C")) {
    return(motif)
  }
  reverse_complement
}

#' Return the reverse complement of a concrete A/C/G/T motif.
#'
#' @param motif Concrete A/C/G/T motif.
#'
#' @return Reverse-complement motif.
#' @noRd
cf_reverse_complement <- function(motif) {
  bases <- strsplit(motif, "", fixed = TRUE)[[1L]]
  paste(rev(chartr("ACGT", "TGCA", bases)), collapse = "")
}

#' Validate reference k-mer frequency values.
#'
#' @param values Frequency values.
#' @param value_name Human-readable value name.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_ref_kmer_frequency_vector <- function(values, value_name) {
  cf_validate_fraction_vector(values, value_name)
}

#' Validate sparse reference k-mer COO coordinate order.
#'
#' @param row_idx0 Zero-based row coordinate vector.
#' @param motif_idx0 Zero-based motif coordinate vector.
#'
#' @return Invisibly returns `TRUE`.
#' @noRd
cf_validate_ref_kmer_sparse_coordinates <- function(row_idx0, motif_idx0) {
  if (length(row_idx0) < 2L) {
    return(invisible(TRUE))
  }
  row_idx0 <- as.integer(row_idx0)
  motif_idx0 <- as.integer(motif_idx0)
  previous_row_idx0 <- row_idx0[-length(row_idx0)]
  previous_motif_idx0 <- motif_idx0[-length(motif_idx0)]
  current_row_idx0 <- row_idx0[-1L]
  current_motif_idx0 <- motif_idx0[-1L]
  is_sorted_unique <- current_row_idx0 > previous_row_idx0 |
    (current_row_idx0 == previous_row_idx0 & current_motif_idx0 > previous_motif_idx0)
  if (!all(is_sorted_unique)) {
    stop("sparse COO entries must be sorted and unique by row, motif", call. = FALSE)
  }
  invisible(TRUE)
}

#' @export
#' @rdname storage_mode
storage_mode.cfdnalab_ref_kmer_frequencies <- function(x, ...) {
  cf_reject_unused_arguments(...)
  x$storage_mode
}

#' @export
#' @rdname row_mode
row_mode.cfdnalab_ref_kmer_frequencies <- function(x, ...) {
  cf_reject_unused_arguments(...)
  x$row_mode
}

#' @export
#' @rdname motif_axis_kind
motif_axis_kind.cfdnalab_ref_kmer_frequencies <- function(x, ...) {
  cf_reject_unused_arguments(...)
  x$motif_axis_kind
}

#' @export
#' @rdname kmer_size
kmer_size.cfdnalab_ref_kmer_frequencies <- function(x, ...) {
  cf_reject_unused_arguments(...)
  x$kmer_size
}

#' @export
#' @rdname canonical
canonical.cfdnalab_ref_kmer_frequencies <- function(x, ...) {
  cf_reject_unused_arguments(...)
  x$canonical
}

#' @export
#' @rdname all_motifs
all_motifs.cfdnalab_ref_kmer_frequencies <- function(x, ...) {
  cf_reject_unused_arguments(...)
  x$all_motifs
}

#' @export
#' @rdname assign_by
assign_by.cfdnalab_ref_kmer_frequencies <- function(x, ...) {
  cf_reject_unused_arguments(...)
  x$assign_by
}

#' @export
#' @rdname reference_contig_footprint
reference_contig_footprint.cfdnalab_ref_kmer_frequencies <- function(x, ...) {
  cf_reject_unused_arguments(...)
  x$reference_contig_footprint
}

#' @export
#' @rdname row_scaling_factors
row_scaling_factors.cfdnalab_ref_kmer_frequencies <- function(x, ...) {
  cf_reject_unused_arguments(...)
  data.frame(
    x$row_metadata,
    row_scaling_factor = x$row_scaling_factor,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' @export
#' @rdname motifs
motifs.cfdnalab_ref_kmer_frequencies <- function(x, ...) {
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
motif_idx.cfdnalab_ref_kmer_frequencies <- function(x, motif, ...) {
  cf_reject_unused_arguments(...)
  motif <- cf_validate_scalar_string(motif, "motif")
  matched_index <- cf_find_unique_value_index(
    x$motif,
    motif,
    "Unknown reference k-mer motif label: ",
    "Reference k-mer motif label is not unique: "
  )
  cf_index0_to_r_index(x$motif_idx0[[matched_index]])
}

#' @export
#' @rdname has_motif
#' @param motif Motif label to test.
has_motif.cfdnalab_ref_kmer_frequencies <- function(x, motif, ...) {
  cf_reject_unused_arguments(...)
  motif <- cf_validate_scalar_string(motif, "motif")
  any(x$motif == motif)
}

#' @export
#' @rdname window_metadata
window_metadata.cfdnalab_windowed_ref_kmer_frequencies <- function(x, ...) {
  cf_reject_unused_arguments(...)
  x$row_metadata
}

#' @export
#' @rdname group_metadata
group_metadata.cfdnalab_grouped_ref_kmer_frequencies <- function(x, ...) {
  cf_reject_unused_arguments(...)
  x$row_metadata
}

#' @export
#' @rdname group_idx
#' @param group_name Group name to look up.
group_idx.cfdnalab_grouped_ref_kmer_frequencies <- function(x, group_name, ...) {
  cf_reject_unused_arguments(...)
  group_name <- cf_validate_scalar_string(group_name, "group_name")
  matched_index <- cf_find_unique_value_index(
    x$row_metadata$group_name,
    group_name,
    "Unknown reference k-mer group name: ",
    "Reference k-mer group name is not unique: "
  )
  x$row_metadata$group_idx[[matched_index]]
}

#' @export
#' @rdname sparse_frequencies_matrix
#' @param motifs Optional motif label vector. Use either `motifs` or
#'   `motif_idxs`, not both.
#' @param motif_idxs Optional one-based motif index vector.
sparse_frequencies_matrix.cfdnalab_global_ref_kmer_frequencies <- function(
  x,
  motifs = NULL,
  motif_idxs = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_sparse_ref_kmer_frequency_matrix_for_indices(
    x,
    row_indices = seq_len(length(x$row_idx0)),
    motif_indices = cf_resolve_ref_kmer_motif_indices(x, motifs, motif_idxs)
  )
}

#' @export
#' @rdname sparse_frequencies_matrix
#' @param window_idxs Optional one-based window index vector for windowed output.
sparse_frequencies_matrix.cfdnalab_windowed_ref_kmer_frequencies <- function(
  x,
  window_idxs = NULL,
  motifs = NULL,
  motif_idxs = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_sparse_ref_kmer_frequency_matrix_for_indices(
    x,
    row_indices = cf_resolve_ref_kmer_window_indices(x, window_idxs),
    motif_indices = cf_resolve_ref_kmer_motif_indices(x, motifs, motif_idxs)
  )
}

#' @export
#' @rdname sparse_frequencies_matrix
#' @param groups Optional group name vector for grouped output. Use either
#'   `groups` or `group_idxs`, not both.
#' @param group_idxs Optional one-based group index vector for grouped output.
sparse_frequencies_matrix.cfdnalab_grouped_ref_kmer_frequencies <- function(
  x,
  groups = NULL,
  group_idxs = NULL,
  motifs = NULL,
  motif_idxs = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_sparse_ref_kmer_frequency_matrix_for_indices(
    x,
    row_indices = cf_resolve_ref_kmer_group_indices(x, groups, group_idxs),
    motif_indices = cf_resolve_ref_kmer_motif_indices(x, motifs, motif_idxs)
  )
}

#' Build a sparse reference k-mer frequency matrix for selected rows and motifs.
#'
#' @param x Reference k-mer object.
#' @param row_indices One-based row indices.
#' @param motif_indices One-based motif indices.
#'
#' @return A `Matrix` sparse matrix.
#' @noRd
cf_sparse_ref_kmer_frequency_matrix_for_indices <- function(x, row_indices, motif_indices) {
  if (identical(x$storage_mode, "dense")) {
    frequencies <- cf_ref_kmer_dense_frequency_matrix_for_indices(x, row_indices, motif_indices)
    return(Matrix::Matrix(frequencies, sparse = TRUE))
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
    x = as.numeric(x$sparse$frequency[matches]),
    dims = as.integer(c(length(row_indices), length(motif_indices)))
  )
}

#' @export
#' @rdname sparse_counts_matrix
sparse_counts_matrix.cfdnalab_global_ref_kmer_frequencies <- function(
  x,
  motifs = NULL,
  motif_idxs = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_sparse_ref_kmer_count_matrix_for_indices(
    x,
    row_indices = seq_len(length(x$row_idx0)),
    motif_indices = cf_resolve_ref_kmer_motif_indices(x, motifs, motif_idxs)
  )
}

#' @export
#' @rdname sparse_counts_matrix
sparse_counts_matrix.cfdnalab_windowed_ref_kmer_frequencies <- function(
  x,
  window_idxs = NULL,
  motifs = NULL,
  motif_idxs = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_sparse_ref_kmer_count_matrix_for_indices(
    x,
    row_indices = cf_resolve_ref_kmer_window_indices(x, window_idxs),
    motif_indices = cf_resolve_ref_kmer_motif_indices(x, motifs, motif_idxs)
  )
}

#' @export
#' @rdname sparse_counts_matrix
sparse_counts_matrix.cfdnalab_grouped_ref_kmer_frequencies <- function(
  x,
  groups = NULL,
  group_idxs = NULL,
  motifs = NULL,
  motif_idxs = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_sparse_ref_kmer_count_matrix_for_indices(
    x,
    row_indices = cf_resolve_ref_kmer_group_indices(x, groups, group_idxs),
    motif_indices = cf_resolve_ref_kmer_motif_indices(x, motifs, motif_idxs)
  )
}

#' Build a sparse reference k-mer count matrix for selected rows and motifs.
#'
#' @param x Reference k-mer object.
#' @param row_indices One-based row indices.
#' @param motif_indices One-based motif indices.
#'
#' @return A `Matrix` sparse matrix.
#' @noRd
cf_sparse_ref_kmer_count_matrix_for_indices <- function(x, row_indices, motif_indices) {
  if (identical(x$storage_mode, "dense")) {
    counts <- cf_dense_ref_kmer_count_matrix_for_indices(
      x,
      row_indices,
      motif_indices,
      allow_densify = TRUE
    )
    return(Matrix::Matrix(counts, sparse = TRUE))
  }
  selected_row_idx0 <- cf_r_index_to_index0(row_indices)
  selected_motif_idx0 <- cf_r_index_to_index0(motif_indices)
  sparse_row_idx0 <- as.integer(x$sparse$row_idx0)
  sparse_motif_idx0 <- as.integer(x$sparse$motif_idx0)
  matches <- sparse_row_idx0 %in% selected_row_idx0 &
    sparse_motif_idx0 %in% selected_motif_idx0
  matched_row_indices <- cf_index0_to_r_index(sparse_row_idx0[matches])
  Matrix::sparseMatrix(
    i = match(sparse_row_idx0[matches], selected_row_idx0),
    j = match(sparse_motif_idx0[matches], selected_motif_idx0),
    x = as.numeric(x$sparse$frequency[matches]) * x$row_scaling_factor[matched_row_indices],
    dims = as.integer(c(length(row_indices), length(motif_indices)))
  )
}

#' @export
#' @rdname dense_frequencies_matrix
#' @param allow_densify If `TRUE`, allow sparse output to be converted to a
#'   zero-filled dense matrix in memory. Zero filling uses the selected motifs
#'   from `motifs(x)`, not every possible k-mer unless `all_motifs(x)` is
#'   `TRUE`. Sparse output errors by default.
#' @param motifs Optional motif label vector. Use either `motifs` or
#'   `motif_idxs`, not both.
#' @param motif_idxs Optional one-based motif index vector.
dense_frequencies_matrix.cfdnalab_global_ref_kmer_frequencies <- function(
  x,
  allow_densify = FALSE,
  motifs = NULL,
  motif_idxs = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_dense_ref_kmer_frequency_matrix_for_indices(
    x,
    row_indices = seq_len(length(x$row_idx0)),
    motif_indices = cf_resolve_ref_kmer_motif_indices(x, motifs, motif_idxs),
    allow_densify = allow_densify
  )
}

#' @export
#' @rdname dense_frequencies_matrix
#' @param window_idxs Optional one-based window index vector for windowed output.
dense_frequencies_matrix.cfdnalab_windowed_ref_kmer_frequencies <- function(
  x,
  allow_densify = FALSE,
  window_idxs = NULL,
  motifs = NULL,
  motif_idxs = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_dense_ref_kmer_frequency_matrix_for_indices(
    x,
    row_indices = cf_resolve_ref_kmer_window_indices(x, window_idxs),
    motif_indices = cf_resolve_ref_kmer_motif_indices(x, motifs, motif_idxs),
    allow_densify = allow_densify
  )
}

#' @export
#' @rdname dense_frequencies_matrix
#' @param groups Optional group name vector for grouped output. Use either
#'   `groups` or `group_idxs`, not both.
#' @param group_idxs Optional one-based group index vector for grouped output.
dense_frequencies_matrix.cfdnalab_grouped_ref_kmer_frequencies <- function(
  x,
  allow_densify = FALSE,
  groups = NULL,
  group_idxs = NULL,
  motifs = NULL,
  motif_idxs = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_dense_ref_kmer_frequency_matrix_for_indices(
    x,
    row_indices = cf_resolve_ref_kmer_group_indices(x, groups, group_idxs),
    motif_indices = cf_resolve_ref_kmer_motif_indices(x, motifs, motif_idxs),
    allow_densify = allow_densify
  )
}

#' Build a dense reference k-mer frequency matrix.
#'
#' @param x Reference k-mer object.
#' @param row_indices One-based row indices.
#' @param motif_indices One-based motif indices.
#' @param allow_densify Whether to allow sparse output to become a dense matrix.
#'
#' @return A dense numeric matrix.
#' @noRd
cf_dense_ref_kmer_frequency_matrix_for_indices <- function(
  x,
  row_indices,
  motif_indices,
  allow_densify
) {
  if (identical(x$storage_mode, "dense")) {
    return(cf_ref_kmer_dense_frequency_matrix_for_indices(x, row_indices, motif_indices))
  }
  if (!isTRUE(allow_densify)) {
    stop(
      "This reference k-mer output is sparse. Use sparse_frequencies_matrix() or set allow_densify = TRUE.",
      call. = FALSE
    )
  }
  as.matrix(cf_sparse_ref_kmer_frequency_matrix_for_indices(x, row_indices, motif_indices))
}

#' Read and validate selected dense frequencies.
#'
#' @param x Reference k-mer object.
#' @param row_indices One-based row indices.
#' @param motif_indices One-based motif indices.
#'
#' @return A dense numeric matrix.
#' @noRd
cf_ref_kmer_dense_frequency_matrix_for_indices <- function(x, row_indices, motif_indices) {
  frequencies <- cf_read_array(x$store, "frequencies", "Reference k-mer")[
    row_indices,
    motif_indices,
    drop = FALSE
  ]
  cf_validate_ref_kmer_frequency_vector(as.numeric(frequencies), "frequencies")
  frequencies
}

#' @export
#' @rdname dense_counts_matrix
dense_counts_matrix.cfdnalab_global_ref_kmer_frequencies <- function(
  x,
  allow_densify = FALSE,
  motifs = NULL,
  motif_idxs = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_dense_ref_kmer_count_matrix_for_indices(
    x,
    row_indices = seq_len(length(x$row_idx0)),
    motif_indices = cf_resolve_ref_kmer_motif_indices(x, motifs, motif_idxs),
    allow_densify = allow_densify
  )
}

#' @export
#' @rdname dense_counts_matrix
dense_counts_matrix.cfdnalab_windowed_ref_kmer_frequencies <- function(
  x,
  allow_densify = FALSE,
  window_idxs = NULL,
  motifs = NULL,
  motif_idxs = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_dense_ref_kmer_count_matrix_for_indices(
    x,
    row_indices = cf_resolve_ref_kmer_window_indices(x, window_idxs),
    motif_indices = cf_resolve_ref_kmer_motif_indices(x, motifs, motif_idxs),
    allow_densify = allow_densify
  )
}

#' @export
#' @rdname dense_counts_matrix
dense_counts_matrix.cfdnalab_grouped_ref_kmer_frequencies <- function(
  x,
  allow_densify = FALSE,
  groups = NULL,
  group_idxs = NULL,
  motifs = NULL,
  motif_idxs = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_dense_ref_kmer_count_matrix_for_indices(
    x,
    row_indices = cf_resolve_ref_kmer_group_indices(x, groups, group_idxs),
    motif_indices = cf_resolve_ref_kmer_motif_indices(x, motifs, motif_idxs),
    allow_densify = allow_densify
  )
}

#' Build a dense reference k-mer count matrix.
#'
#' @param x Reference k-mer object.
#' @param row_indices One-based row indices.
#' @param motif_indices One-based motif indices.
#' @param allow_densify Whether to allow sparse output to become a dense matrix.
#'
#' @return A dense numeric matrix.
#' @noRd
cf_dense_ref_kmer_count_matrix_for_indices <- function(
  x,
  row_indices,
  motif_indices,
  allow_densify
) {
  frequencies <- cf_dense_ref_kmer_frequency_matrix_for_indices(
    x,
    row_indices = row_indices,
    motif_indices = motif_indices,
    allow_densify = allow_densify
  )
  scaling <- matrix(
    x$row_scaling_factor[row_indices],
    nrow = length(row_indices),
    ncol = length(motif_indices)
  )
  frequencies * scaling
}

#' @export
#' @rdname dense_frequencies_vector
#' @param allow_densify If `TRUE`, allow sparse output to be converted to a
#'   zero-filled dense vector in memory. Zero filling uses the motif axis from
#'   `motifs(x)`, not every possible k-mer unless `all_motifs(x)` is `TRUE`.
dense_frequencies_vector.cfdnalab_global_ref_kmer_frequencies <- function(
  x,
  allow_densify = FALSE,
  ...
) {
  cf_reject_unused_arguments(...)
  frequencies <- as.vector(dense_frequencies_matrix(x, allow_densify = allow_densify)[1L, ])
  stats::setNames(frequencies, x$motif)
}

#' @export
#' @rdname dense_counts_vector
dense_counts_vector.cfdnalab_global_ref_kmer_frequencies <- function(
  x,
  allow_densify = FALSE,
  ...
) {
  cf_reject_unused_arguments(...)
  counts <- as.vector(dense_counts_matrix(x, allow_densify = allow_densify)[1L, ])
  stats::setNames(counts, x$motif)
}

#' @export
#' @rdname ref_kmer_data_frame
#' @param densify If `TRUE`, sparse output adds explicit zero-frequency rows
#'   for selected motifs in `motifs(x)`. For observed-only output, this is the
#'   combined set observed anywhere in the output. Densifying does not add
#'   every possible k-mer unless `all_motifs(x)` is `TRUE`. Dense outputs
#'   ignore this option.
#' @param motifs Optional motif label vector. Use either `motifs` or
#'   `motif_idxs`, not both.
#' @param motif_idxs Optional one-based motif index vector.
ref_kmer_data_frame.cfdnalab_global_ref_kmer_frequencies <- function(
  x,
  densify = FALSE,
  motifs = NULL,
  motif_idxs = NULL,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_ref_kmer_data_frame(
    x,
    row_indices = seq_len(length(x$row_idx0)),
    motif_indices = cf_resolve_ref_kmer_motif_indices(x, motifs, motif_idxs),
    densify = densify,
    max_blacklisted_fraction = 1.0
  )
}

#' @export
#' @rdname ref_kmer_data_frame
#' @param window_idxs Optional one-based window index vector for windowed output.
#' @param max_blacklisted_fraction Maximum row `blacklisted_fraction` in 0..1
#'   to retain before returning values. The default `1.0` keeps all selected
#'   rows.
ref_kmer_data_frame.cfdnalab_windowed_ref_kmer_frequencies <- function(
  x,
  window_idxs = NULL,
  densify = FALSE,
  motifs = NULL,
  motif_idxs = NULL,
  max_blacklisted_fraction = 1.0,
  ...
) {
  cf_reject_unused_arguments(...)
  cf_ref_kmer_data_frame(
    x,
    row_indices = cf_resolve_ref_kmer_window_indices(x, window_idxs),
    motif_indices = cf_resolve_ref_kmer_motif_indices(x, motifs, motif_idxs),
    densify = densify,
    max_blacklisted_fraction = max_blacklisted_fraction
  )
}

#' @export
#' @rdname ref_kmer_data_frame
#' @param groups Optional group name vector for grouped output. Use either
#'   `groups` or `group_idxs`, not both.
#' @param group_idxs Optional one-based group index vector for grouped output.
ref_kmer_data_frame.cfdnalab_grouped_ref_kmer_frequencies <- function(
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
  cf_ref_kmer_data_frame(
    x,
    row_indices = cf_resolve_ref_kmer_group_indices(x, groups, group_idxs),
    motif_indices = cf_resolve_ref_kmer_motif_indices(x, motifs, motif_idxs),
    densify = densify,
    max_blacklisted_fraction = max_blacklisted_fraction
  )
}

#' Shared implementation for mode-specific reference k-mer data frame methods.
#'
#' @param x Reference k-mer object.
#' @param row_indices One-based row indices.
#' @param motif_indices One-based motif indices.
#' @param densify Whether to add explicit zero-frequency rows for sparse output.
#' @param max_blacklisted_fraction Maximum blacklist fraction.
#'
#' @return A data frame.
#' @noRd
cf_ref_kmer_data_frame <- function(
  x,
  row_indices,
  motif_indices,
  densify,
  max_blacklisted_fraction
) {
  cf_validate_scalar_logical(densify, "densify")
  row_indices <- cf_apply_ref_kmer_blacklist_filter(x, row_indices, max_blacklisted_fraction)
  if (identical(x$storage_mode, "sparse_coo") && !isTRUE(densify)) {
    return(cf_stored_ref_kmer_data_frame_for_indices(x, row_indices, motif_indices))
  }
  cf_complete_ref_kmer_data_frame_for_indices(x, row_indices, motif_indices, densify)
}

#' Apply a blacklist fraction filter to reference k-mer row indices.
#'
#' @param x Reference k-mer object.
#' @param row_indices One-based row indices.
#' @param max_blacklisted_fraction Maximum blacklist fraction.
#'
#' @return Filtered one-based row indices.
#' @noRd
cf_apply_ref_kmer_blacklist_filter <- function(x, row_indices, max_blacklisted_fraction) {
  cf_apply_row_blacklist_filter(x$row_metadata, row_indices, max_blacklisted_fraction)
}

#' Resolve grouped reference k-mer selectors to one-based row indices.
#'
#' @param x Grouped reference k-mer object.
#' @param groups Optional group names.
#' @param group_idxs Optional one-based group indices.
#'
#' @return One-based row indices.
#' @noRd
cf_resolve_ref_kmer_group_indices <- function(x, groups, group_idxs) {
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
          "Unknown reference k-mer group name: ",
          "Reference k-mer group name is not unique: "
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

#' Resolve windowed reference k-mer selectors to one-based row indices.
#'
#' @param x Windowed reference k-mer object.
#' @param window_idxs Optional one-based window indices.
#'
#' @return One-based row indices.
#' @noRd
cf_resolve_ref_kmer_window_indices <- function(x, window_idxs) {
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

#' Resolve reference k-mer motif selectors to one-based motif indices.
#'
#' @param x Reference k-mer object.
#' @param motifs Optional motif labels.
#' @param motif_idxs Optional one-based motif indices.
#'
#' @return One-based motif indices.
#' @noRd
cf_resolve_ref_kmer_motif_indices <- function(x, motifs, motif_idxs) {
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

#' Build a complete reference k-mer data frame for selected rows and motifs.
#'
#' @param x Reference k-mer object.
#' @param row_indices One-based row indices.
#' @param motif_indices One-based motif indices.
#' @param densify Whether to add explicit zero-frequency rows for sparse output.
#'
#' @return A data frame with one row per selected row and motif.
#' @noRd
cf_complete_ref_kmer_data_frame_for_indices <- function(x, row_indices, motif_indices, densify) {
  if (length(row_indices) == 0L || length(motif_indices) == 0L) {
    return(cf_empty_ref_kmer_data_frame(x))
  }
  frequencies <- cf_dense_ref_kmer_frequency_matrix_for_indices(
    x,
    row_indices = row_indices,
    motif_indices = motif_indices,
    allow_densify = densify
  )
  counts <- frequencies * matrix(
    x$row_scaling_factor[row_indices],
    nrow = length(row_indices),
    ncol = length(motif_indices)
  )
  num_rows <- length(row_indices)
  num_motifs <- length(motif_indices)
  metadata <- x$row_metadata[row_indices, , drop = FALSE]
  metadata <- metadata[rep(seq_len(num_rows), each = num_motifs), , drop = FALSE]
  motif_metadata <- motifs(x)[motif_indices, , drop = FALSE]
  motif_metadata <- motif_metadata[rep(seq_len(num_motifs), times = num_rows), , drop = FALSE]
  data.frame(
    metadata,
    motif_metadata,
    frequency = as.vector(t(frequencies)),
    count = as.vector(t(counts)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Build a reference k-mer data frame from stored COO rows.
#'
#' @param x Reference k-mer object.
#' @param row_indices One-based row indices.
#' @param motif_indices One-based motif indices.
#'
#' @return A data frame with one row per stored non-zero frequency.
#' @noRd
cf_stored_ref_kmer_data_frame_for_indices <- function(x, row_indices, motif_indices) {
  if (length(row_indices) == 0L || length(motif_indices) == 0L) {
    return(cf_empty_ref_kmer_data_frame(x))
  }
  selected_row_idx0 <- cf_r_index_to_index0(row_indices)
  selected_motif_idx0 <- cf_r_index_to_index0(motif_indices)
  sparse_row_idx0 <- as.integer(x$sparse$row_idx0)
  sparse_motif_idx0 <- as.integer(x$sparse$motif_idx0)
  matches <- sparse_row_idx0 %in% selected_row_idx0 &
    sparse_motif_idx0 %in% selected_motif_idx0
  if (!any(matches)) {
    return(cf_empty_ref_kmer_data_frame(x))
  }
  matched_row_idx0 <- sparse_row_idx0[matches]
  matched_motif_idx0 <- sparse_motif_idx0[matches]
  sort_order <- order(
    match(matched_row_idx0, selected_row_idx0),
    match(matched_motif_idx0, selected_motif_idx0)
  )
  matched_row_idx0 <- matched_row_idx0[sort_order]
  matched_motif_idx0 <- matched_motif_idx0[sort_order]
  matched_frequencies <- as.numeric(x$sparse$frequency[matches])[sort_order]
  matched_row_indices <- cf_index0_to_r_index(matched_row_idx0)
  matched_motif_indices <- cf_index0_to_r_index(matched_motif_idx0)
  data.frame(
    x$row_metadata[matched_row_indices, , drop = FALSE],
    motifs(x)[matched_motif_indices, , drop = FALSE],
    frequency = matched_frequencies,
    count = matched_frequencies * x$row_scaling_factor[matched_row_indices],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Build an empty reference k-mer data frame with public columns.
#'
#' @param x Reference k-mer object.
#'
#' @return A zero-row data frame.
#' @noRd
cf_empty_ref_kmer_data_frame <- function(x) {
  data.frame(
    x$row_metadata[integer(0), , drop = FALSE],
    motifs(x)[integer(0), , drop = FALSE],
    frequency = numeric(),
    count = numeric(),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Print a reference k-mer object.
#'
#' @param x A cfDNAlab reference k-mer object.
#' @param ... Ignored.
#'
#' @return Invisibly returns `x`.
#' @export
#' @keywords internal
print.cfdnalab_ref_kmer_frequencies <- function(x, ...) {
  cat("<cfDNAlab reference k-mer frequencies>\n")
  cat("Path: ", x$path, "\n", sep = "")
  cat("Storage mode: ", x$storage_mode, "\n", sep = "")
  cat("Row mode: ", x$row_mode, "\n", sep = "")
  cat("Motif axis kind: ", x$motif_axis_kind, "\n", sep = "")
  cat("Rows: ", length(x$row_idx0), "\n", sep = "")
  cat("Motifs: ", length(x$motif), "\n", sep = "")
  invisible(x)
}
