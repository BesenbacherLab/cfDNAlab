test_that("dense windowed reference k-mers expose metadata, frequencies, and counts", {
  ref_kmers <- read_ref_kmers(make_dense_windowed_ref_kmer_zarr_fixture())

  expect_s3_class(ref_kmers, "cfdnalab_windowed_ref_kmer_frequencies")
  expect_equal(schema_version(ref_kmers), 1L)
  expect_equal(storage_mode(ref_kmers), "dense")
  expect_equal(row_mode(ref_kmers), "bed")
  expect_equal(motif_axis_kind(ref_kmers), "motif")
  expect_equal(kmer_size(ref_kmers), 2L)
  expect_false(canonical(ref_kmers))
  expect_false(all_motifs(ref_kmers))
  expect_equal(assign_by(ref_kmers), "count-overlap")
  expect_output(print(ref_kmers), "<cfDNAlab reference k-mer frequencies>", fixed = TRUE)
  expect_equal(
    vapply(reference_contig_footprint(ref_kmers), `[[`, character(1L), "name"),
    c("chr2", "chr10")
  )
  expect_equal(
    motifs(ref_kmers),
    data.frame(
      motif_idx = 1:3,
      motif = c("AA", "AC", "GT"),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_true(has_motif(ref_kmers, "AC"))
  expect_false(has_motif(ref_kmers, "TT"))
  expect_equal(motif_idx(ref_kmers, "GT"), 3L)
  expect_equal(
    window_metadata(ref_kmers),
    data.frame(
      window_idx = c(1L, 2L),
      chrom = c("chr2", "chr10"),
      start = c(10L, 40L),
      end = c(20L, 60L),
      blacklisted_fraction = c(0.25, 0),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_equal(
    row_scaling_factors(ref_kmers)$row_scaling_factor,
    c(4, 2),
    tolerance = 1e-8
  )
  expect_equal(
    dense_frequencies_matrix(ref_kmers),
    matrix(c(0.25, 0, 0.75, 0.5, 0.5, 0), nrow = 2L, byrow = TRUE),
    tolerance = 1e-8
  )
  expect_equal(
    dense_counts_matrix(ref_kmers),
    matrix(c(1, 0, 3, 1, 1, 0), nrow = 2L, byrow = TRUE),
    tolerance = 1e-8
  )
  expect_equal(
    dense_counts_matrix(ref_kmers, window_idxs = 2L, motifs = c("AC", "AA")),
    matrix(c(1, 1), nrow = 1L),
    tolerance = 1e-8
  )
  expect_equal(
    as.matrix(sparse_frequencies_matrix(ref_kmers, motifs = "GT")),
    matrix(c(0.75, 0), nrow = 2L),
    tolerance = 1e-8
  )
  expect_equal(
    as.matrix(sparse_counts_matrix(ref_kmers, motifs = "GT")),
    matrix(c(3, 0), nrow = 2L),
    tolerance = 1e-8
  )
})

test_that("dense windowed reference k-mer data frames use one-based selectors", {
  ref_kmers <- read_ref_kmers(make_dense_windowed_ref_kmer_zarr_fixture())

  first_window <- ref_kmer_data_frame(ref_kmers, window_idxs = 1L)
  expect_equal(names(first_window), c(
    "window_idx",
    "chrom",
    "start",
    "end",
    "blacklisted_fraction",
    "motif_idx",
    "motif",
    "frequency",
    "count"
  ))
  expect_equal(first_window$window_idx, c(1L, 1L, 1L))
  expect_equal(first_window$motif, c("AA", "AC", "GT"))
  expect_equal(first_window$frequency, c(0.25, 0, 0.75), tolerance = 1e-8)
  expect_equal(first_window$count, c(1, 0, 3), tolerance = 1e-8)

  selected <- ref_kmer_data_frame(ref_kmers, motifs = "AC", max_blacklisted_fraction = 0.1)
  expect_equal(selected$window_idx, 2L)
  expect_equal(selected$frequency, 0.5, tolerance = 1e-8)
  expect_equal(selected$count, 1, tolerance = 1e-8)
  expect_error(
    ref_kmer_data_frame(ref_kmers, motifs = "AA", motif_idxs = 1L),
    "Use either motifs or motif_idxs",
    fixed = TRUE
  )
  expect_error(
    ref_kmer_data_frame(ref_kmers, window_idxs = c(1L, 1L)),
    "window_idxs contains duplicate values",
    fixed = TRUE
  )
  expect_error(
    ref_kmer_data_frame(ref_kmers, window_idxs = 0L),
    "window_idxs contains values outside 1..2",
    fixed = TRUE
  )
})

test_that("sparse grouped reference k-mers densify only when requested", {
  ref_kmers <- read_ref_kmers(make_sparse_grouped_ref_kmer_zarr_fixture())

  expect_s3_class(ref_kmers, "cfdnalab_grouped_ref_kmer_frequencies")
  expect_equal(storage_mode(ref_kmers), "sparse_coo")
  expect_equal(row_mode(ref_kmers), "grouped_bed")
  expect_equal(group_idx(ref_kmers, "long_group"), 2L)
  expect_equal(
    group_metadata(ref_kmers),
    data.frame(
      group_idx = 1:3,
      group_name = c("A", "long_group", "empty"),
      eligible_windows = c(1L, 2L, 0L),
      blacklisted_fraction = c(0, 0.125, 0),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_error(dense_frequencies_matrix(ref_kmers), "Use sparse_frequencies_matrix")
  expect_equal(
    dense_frequencies_matrix(ref_kmers, allow_densify = TRUE),
    matrix(c(0.25, 0, 0.75, 0, 1, 0, 0, 0, 0), nrow = 3L, byrow = TRUE),
    tolerance = 1e-8
  )
  expect_equal(
    dense_counts_matrix(ref_kmers, allow_densify = TRUE),
    matrix(c(1, 0, 3, 0, 2, 0, 0, 0, 0), nrow = 3L, byrow = TRUE),
    tolerance = 1e-8
  )
  expect_equal(
    as.matrix(sparse_counts_matrix(ref_kmers, motifs = "GT")),
    matrix(c(3, 0, 0), nrow = 3L),
    tolerance = 1e-8
  )
  expect_equal(
    dense_counts_matrix(
      ref_kmers,
      groups = c("long_group", "A"),
      motifs = c("GT", "AA"),
      allow_densify = TRUE
    ),
    matrix(c(0, 0, 3, 1), nrow = 2L, byrow = TRUE),
    tolerance = 1e-8
  )
})

test_that("sparse grouped reference k-mer frames preserve selector order", {
  ref_kmers <- read_ref_kmers(make_sparse_grouped_ref_kmer_zarr_fixture())

  stored <- ref_kmer_data_frame(ref_kmers)
  expect_equal(stored$group_idx, c(1L, 1L, 2L))
  expect_equal(stored$group_name, c("A", "A", "long_group"))
  expect_equal(stored$motif_idx, c(1L, 3L, 2L))
  expect_equal(stored$motif, c("AA", "GT", "AC"))
  expect_equal(stored$frequency, c(0.25, 0.75, 1), tolerance = 1e-8)
  expect_equal(stored$count, c(1, 3, 2), tolerance = 1e-8)
  expect_false(any(grepl("idx0|index0", names(stored))))

  dense <- ref_kmer_data_frame(
    ref_kmers,
    groups = c("long_group", "A"),
    motifs = c("GT", "AA"),
    densify = TRUE
  )
  expect_equal(dense$group_name, c("long_group", "long_group", "A", "A"))
  expect_equal(dense$motif, c("GT", "AA", "GT", "AA"))
  expect_equal(dense$count, c(0, 0, 3, 1), tolerance = 1e-8)

  empty <- ref_kmer_data_frame(ref_kmers, groups = "empty", densify = TRUE)
  expect_equal(empty$motif, c("AA", "AC", "GT"))
  expect_equal(empty$frequency, c(0, 0, 0), tolerance = 1e-8)
  expect_equal(empty$count, c(0, 0, 0), tolerance = 1e-8)
})

test_that("global reference k-mer motif-group output uses motif selectors", {
  ref_kmers <- read_ref_kmers(make_dense_global_ref_kmer_group_zarr_fixture())

  expect_s3_class(ref_kmers, "cfdnalab_global_ref_kmer_frequencies")
  expect_equal(row_mode(ref_kmers), "global")
  expect_equal(motif_axis_kind(ref_kmers), "motif_group")
  expect_equal(
    motifs(ref_kmers),
    data.frame(
      motif_idx = c(1L, 2L),
      motif = c("left", "right"),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_equal(motif_idx(ref_kmers, "right"), 2L)
  expect_equal(
    dense_frequencies_vector(ref_kmers),
    c(left = 0.25, right = 0.75),
    tolerance = 1e-8
  )
  expect_equal(
    dense_counts_vector(ref_kmers),
    c(left = 1, right = 3),
    tolerance = 1e-8
  )
  frame <- ref_kmer_data_frame(ref_kmers)
  expect_equal(frame$row_label, c("global", "global"))
  expect_equal(frame$motif, c("left", "right"))
  expect_equal(frame$frequency, c(0.25, 0.75), tolerance = 1e-8)
  expect_equal(frame$count, c(1, 3), tolerance = 1e-8)
})

test_that("reference k-mer loader rejects schema and shape problems", {
  wrong_schema <- make_dense_windowed_ref_kmer_zarr_fixture(
    root_attributes = ref_kmer_root_attributes(
      storage_mode = "dense",
      row_mode = "bed",
      schema = "other"
    )
  )
  wrong_version <- make_dense_windowed_ref_kmer_zarr_fixture(
    root_attributes = ref_kmer_root_attributes(
      storage_mode = "dense",
      row_mode = "bed",
      schema_version = 99L
    )
  )
  missing_scaling <- make_dense_windowed_ref_kmer_zarr_fixture()
  unlink(file.path(missing_scaling, "row_scaling_factor"), recursive = TRUE)
  wrong_dimensions <- make_dense_windowed_ref_kmer_zarr_fixture()
  patch_zarr_metadata(wrong_dimensions, "frequencies", dimension_names = c("motif", "row"))
  wrong_motif_axis_dimensions <- make_dense_windowed_ref_kmer_zarr_fixture()
  patch_zarr_metadata(wrong_motif_axis_dimensions, "motif_index", dimension_names = "row")
  wrong_row_metadata_dimensions <- make_dense_windowed_ref_kmer_zarr_fixture()
  patch_zarr_metadata(wrong_row_metadata_dimensions, "row_start_bp", dimension_names = "motif")
  shape_mismatch <- make_sparse_grouped_ref_kmer_zarr_fixture(sparse_shape = c(3L, 2L))
  duplicate_sparse_coordinate <- make_sparse_grouped_ref_kmer_zarr_fixture(
    sparse_row = c(0L, 0L, 0L),
    sparse_motif = c(0L, 0L, 2L)
  )
  unsorted_sparse_coordinate <- make_sparse_grouped_ref_kmer_zarr_fixture(
    sparse_row = c(0L, 0L, 0L),
    sparse_motif = c(2L, 0L, 1L)
  )
  wrong_sparse_dimensions <- make_sparse_grouped_ref_kmer_zarr_fixture(
    sparse_dimension_labels = list("motif", "row")
  )

  expect_error(read_ref_kmers(wrong_schema), "Expected cfdnalab_schema")
  expect_error(read_ref_kmers(wrong_version), "Unsupported reference k-mer schema version")
  expect_error(read_ref_kmers(missing_scaling), "missing arrays: row_scaling_factor")
  expect_error(read_ref_kmers(wrong_dimensions), "frequencies dimensions must be")
  expect_error(
    read_ref_kmers(wrong_motif_axis_dimensions),
    "motif_index dimensions must be"
  )
  expect_error(
    read_ref_kmers(wrong_row_metadata_dimensions),
    "row_start_bp dimensions must be"
  )
  expect_error(read_ref_kmers(shape_mismatch), "sparse/shape does not match")
  expect_error(read_ref_kmers(duplicate_sparse_coordinate), "sorted and unique")
  expect_error(read_ref_kmers(unsorted_sparse_coordinate), "sorted and unique")
  expect_error(read_ref_kmers(wrong_sparse_dimensions), "sparse_dimension labels")
})

test_that("reference k-mer loader rejects metadata that changes count meaning", {
  wrong_units <- make_dense_windowed_ref_kmer_zarr_fixture(
    root_attributes = ref_kmer_root_attributes(
      storage_mode = "dense",
      row_mode = "bed",
      value_units = "other"
    )
  )
  wrong_scaling_array <- make_dense_windowed_ref_kmer_zarr_fixture(
    root_attributes = ref_kmer_root_attributes(
      storage_mode = "dense",
      row_mode = "bed",
      row_scaling_factor_array = "other"
    )
  )
  wrong_reconstruction <- make_dense_windowed_ref_kmer_zarr_fixture(
    root_attributes = ref_kmer_root_attributes(
      storage_mode = "dense",
      row_mode = "bed",
      count_reconstruction = "count = frequency"
    )
  )

  expect_error(read_ref_kmers(wrong_units), "value_units")
  expect_error(read_ref_kmers(wrong_scaling_array), "row_scaling_factor_array")
  expect_error(read_ref_kmers(wrong_reconstruction), "count_reconstruction")
})

test_that("reference k-mer loader rejects invalid values and motif labels", {
  bad_sparse_frequency <- make_sparse_grouped_ref_kmer_zarr_fixture(
    sparse_frequency = c(0.25, 1.25, 1)
  )
  bad_scaling <- make_dense_windowed_ref_kmer_zarr_fixture(
    row_scaling_factor = c(4, NaN)
  )
  invalid_base <- make_dense_windowed_ref_kmer_zarr_fixture(
    motifs = c("AA", "AN"),
    frequencies = matrix(0, nrow = 2L, ncol = 2L)
  )
  noncanonical <- make_dense_windowed_ref_kmer_zarr_fixture(
    motifs = c("AA", "GT"),
    canonical = TRUE,
    frequencies = matrix(0, nrow = 2L, ncol = 2L)
  )
  bad_dense_frequency <- make_dense_windowed_ref_kmer_zarr_fixture(
    frequencies = matrix(c(0.25, 1.25, 0, 0.5, 0.5, 0), nrow = 2L, byrow = TRUE)
  )

  expect_error(read_ref_kmers(bad_sparse_frequency), "sparse/frequency")
  expect_error(read_ref_kmers(bad_scaling), "row_scaling_factor")
  expect_error(read_ref_kmers(invalid_base), "invalid base")
  expect_error(read_ref_kmers(noncanonical), "canonical reference k-mer motif label")

  ref_kmers <- read_ref_kmers(bad_dense_frequency)
  expect_error(dense_frequencies_matrix(ref_kmers), "frequencies")
})

test_that("reference k-mer loader rejects invalid row metadata", {
  bad_interval <- make_dense_windowed_ref_kmer_zarr_fixture(
    row_start_bp = c(10L, 60L),
    row_end_bp = c(20L, 40L)
  )
  bad_window_fraction <- make_dense_windowed_ref_kmer_zarr_fixture(
    blacklisted_fraction = c(0.25, 1.25)
  )
  bad_chromosome_index <- make_dense_windowed_ref_kmer_zarr_fixture(
    row_chromosome = c(0L, 2L)
  )
  bad_eligible_windows <- make_sparse_grouped_ref_kmer_zarr_fixture(
    eligible_windows = c(1L, -1L, 0L)
  )
  bad_group_fraction <- make_sparse_grouped_ref_kmer_zarr_fixture(
    blacklisted_fraction = c(0, NaN, 0)
  )

  expect_error(
    read_ref_kmers(bad_interval),
    "row_start_bp must be smaller than row_end_bp"
  )
  expect_error(read_ref_kmers(bad_window_fraction), "blacklisted_fraction")
  expect_error(
    read_ref_kmers(bad_chromosome_index),
    "row_chromosome contains an index outside"
  )
  expect_error(read_ref_kmers(bad_eligible_windows), "eligible_windows")
  expect_error(read_ref_kmers(bad_group_fraction), "blacklisted_fraction")
})

test_that("reference k-mer loader rejects invalid JSON labels", {
  numeric_group_labels <- make_dense_global_ref_kmer_group_zarr_fixture()
  patch_zarr_metadata(
    numeric_group_labels,
    "motif_index",
    attributes = list(label_field = "motif_group", labels = list(1L, 2L))
  )
  control_character_label <- make_sparse_grouped_ref_kmer_zarr_fixture(
    group_labels = list("A", "bad\nlabel", "empty")
  )

  expect_error(read_ref_kmers(numeric_group_labels), "labels must be character strings")
  expect_error(
    read_ref_kmers(control_character_label),
    "labels must not contain control characters"
  )
})
