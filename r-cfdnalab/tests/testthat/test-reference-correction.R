test_that("reference correction keeps end-motif counts on count scale", {
  end_path <- make_dense_windowed_end_motif_zarr_fixture()
  ref_path <- make_dense_windowed_ref_kmer_zarr_fixture(
    motifs = c("A", "G", "T"),
    frequencies = matrix(
      c(1 / 3, 1 / 6, 1 / 2, 1 / 2, 1 / 4, 1 / 4),
      nrow = 2L,
      byrow = TRUE
    ),
    chromosome_names = c("chr1", "chr2"),
    row_chromosome = c(0L, 1L),
    row_start_bp = c(10L, 30L),
    row_end_bp = c(12L, 36L),
    blacklisted_fraction = c(0, 0.125),
    row_scaling_factor = c(6, 4)
  )
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  corrected <- end_motif_data_frame(ends, ref_kmers = ref_kmers)

  expect_equal(corrected$reference_motif, c("A", "G", "T", "A", "G", "T"))
  expect_equal(corrected$correction_motif_count, c(3L, 3L, 3L, 3L, 3L, 3L))
  expect_equal(corrected$reference_scale, c(1, 0.5, 1.5, 1.5, 0.75, 0.75))
  expect_equal(
    corrected$reference_corrected_count,
    c(0, 4, 0, 1, 0, 16 / 3)
  )
})

test_that("global reference correction keeps end-motif counts on count scale", {
  end_path <- make_dense_global_end_motif_zarr_fixture()
  ref_path <- make_dense_global_ref_kmer_zarr_fixture(
    motifs = c("A", "C", "G", "T"),
    frequencies = matrix(c(1 / 4, 1 / 8, 1 / 2, 1 / 8), nrow = 1L),
    row_scaling_factor = 8
  )
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  corrected <- end_motif_data_frame(ends, ref_kmers = ref_kmers)

  expect_equal(corrected$reference_motif, c("A", "C", "G", "T"))
  expect_equal(corrected$correction_motif_count, c(4L, 4L, 4L, 4L))
  expect_equal(corrected$reference_scale, c(1, 0.5, 2, 0.5))
  expect_equal(corrected$reference_corrected_count, c(1, 0, 1.25, 0))
})

test_that("reference correction selectors match filtering full correction", {
  end_path <- make_dense_windowed_end_motif_zarr_fixture()
  ref_path <- make_dense_windowed_ref_kmer_zarr_fixture(
    motifs = c("A", "G", "T"),
    frequencies = matrix(
      c(1 / 3, 1 / 6, 1 / 2, 1 / 2, 1 / 4, 1 / 4),
      nrow = 2L,
      byrow = TRUE
    ),
    chromosome_names = c("chr1", "chr2"),
    row_chromosome = c(0L, 1L),
    row_start_bp = c(10L, 30L),
    row_end_bp = c(12L, 36L),
    blacklisted_fraction = c(0, 0.125),
    row_scaling_factor = c(6, 4)
  )
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)
  full_correction <- end_motif_data_frame(ends, ref_kmers = ref_kmers)

  selected <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    window_idxs = 2L,
    motifs = "_A"
  )

  expected <- full_correction[
    full_correction$window_idx == 2L & full_correction$motif == "_A",
    ,
    drop = FALSE
  ]
  row.names(expected) <- NULL
  expect_equal(selected, expected)
  expect_equal(selected$correction_motif_count, 3L)
  expect_equal(selected$reference_corrected_count, 1)
})

test_that("reference correction blacklist filtering uses selected end rows", {
  end_path <- make_dense_windowed_end_motif_zarr_fixture()
  ref_path <- make_dense_windowed_ref_kmer_zarr_fixture(
    motifs = c("A", "G", "T"),
    frequencies = matrix(
      c(1 / 3, 1 / 6, 1 / 2, 1 / 2, 1 / 4, 1 / 4),
      nrow = 2L,
      byrow = TRUE
    ),
    chromosome_names = c("chr1", "chr2"),
    row_chromosome = c(0L, 1L),
    row_start_bp = c(10L, 30L),
    row_end_bp = c(12L, 36L),
    blacklisted_fraction = c(0.25, 0),
    row_scaling_factor = c(6, 4)
  )
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  corrected <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    max_blacklisted_fraction = 0.1
  )

  expect_equal(corrected$window_idx, c(1L, 1L, 1L))
  expect_equal(corrected$reference_corrected_count, c(0, 4, 0))
})

test_that("reference-corrected matrix extractors keep selected shape", {
  end_path <- make_dense_windowed_end_motif_zarr_fixture()
  ref_path <- make_dense_windowed_ref_kmer_zarr_fixture(
    motifs = c("A", "G", "T"),
    frequencies = matrix(
      c(1 / 3, 1 / 6, 1 / 2, 1 / 2, 1 / 4, 1 / 4),
      nrow = 2L,
      byrow = TRUE
    ),
    chromosome_names = c("chr1", "chr2"),
    row_chromosome = c(0L, 1L),
    row_start_bp = c(10L, 30L),
    row_end_bp = c(12L, 36L),
    blacklisted_fraction = c(0, 0.125),
    row_scaling_factor = c(6, 4)
  )
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)
  expected <- matrix(c(1, 0, 0, 4), nrow = 2L, byrow = TRUE)

  corrected_dense <- dense_corrected_counts_matrix(
    ends,
    ref_kmers,
    window_idxs = c(2L, 1L),
    motifs = c("_A", "_G")
  )
  corrected_sparse <- sparse_corrected_counts_matrix(
    ends,
    ref_kmers,
    window_idxs = c(2L, 1L),
    motifs = c("_A", "_G")
  )

  expect_equal(corrected_dense, expected)
  expect_equal(as.matrix(corrected_sparse), expected)
})

test_that("sparse reference-corrected matrix uses sparse end-motif input", {
  end_path <- make_sparse_windowed_end_motif_zarr_fixture()
  ref_path <- make_dense_windowed_ref_kmer_zarr_fixture(
    motifs = c("A", "G", "T"),
    frequencies = matrix(
      c(1 / 3, 1 / 6, 1 / 2, 1 / 2, 1 / 4, 1 / 4, 1 / 3, 1 / 3, 1 / 3),
      nrow = 3L,
      byrow = TRUE
    ),
    chromosome_names = c("chr1", "chr2"),
    row_chromosome = c(0L, 0L, 1L),
    row_start_bp = c(10L, 20L, 30L),
    row_end_bp = c(12L, 25L, 36L),
    blacklisted_fraction = c(0, 0.25, 0),
    row_scaling_factor = c(6, 4, 3)
  )
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  corrected_sparse <- sparse_corrected_counts_matrix(ends, ref_kmers)

  expect_equal(
    as.matrix(corrected_sparse),
    matrix(c(0, 4, 0, 1, 0, 16 / 3, 0, 0, 3), nrow = 3L, byrow = TRUE)
  )
  expect_error(
    dense_corrected_counts_matrix(ends, ref_kmers),
    "sparse_corrected_counts_matrix",
    fixed = TRUE
  )
})

test_that("sparse end and sparse reference correction uses sparse support", {
  end_path <- make_sparse_grouped_end_motif_zarr_fixture()
  ref_path <- make_sparse_grouped_ref_kmer_zarr_fixture(
    motifs = c("A", "C", "G"),
    sparse_row = c(0L, 0L),
    sparse_motif = c(0L, 2L),
    sparse_frequency = c(1 / 3, 2 / 3),
    group_labels = list("alpha", "beta"),
    sparse_shape = c(2L, 3L),
    row_scaling_factor = c(6, 0),
    eligible_windows = c(2L, 0L),
    blacklisted_fraction = c(0.125, 0)
  )
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  corrected_sparse <- sparse_corrected_counts_matrix(ends, ref_kmers)

  expect_equal(
    as.matrix(corrected_sparse),
    matrix(c(1.5, 0, 3.75, 0, 0, 0), nrow = 2L, byrow = TRUE)
  )
})

test_that("reference correction rejects positive counts at zero reference frequency", {
  end_path <- make_dense_windowed_end_motif_zarr_fixture()
  ref_path <- make_dense_windowed_ref_kmer_zarr_fixture(
    motifs = c("A", "G", "T"),
    frequencies = matrix(
      c(1 / 3, 1 / 6, 1 / 2, 1 / 2, 1 / 2, 0),
      nrow = 2L,
      byrow = TRUE
    ),
    chromosome_names = c("chr1", "chr2"),
    row_chromosome = c(0L, 1L),
    row_start_bp = c(10L, 30L),
    row_end_bp = c(12L, 36L),
    blacklisted_fraction = c(0, 0.125),
    row_scaling_factor = c(6, 4)
  )
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  expect_error(
    end_motif_data_frame(ends, ref_kmers = ref_kmers),
    "Positive-count end motifs have no positive reference frequency",
    fixed = TRUE
  )
})

test_that("reference correction uses row sparse reference support", {
  end_path <- make_dense_grouped_end_motif_zarr_fixture()
  ref_path <- make_sparse_grouped_ref_kmer_zarr_fixture(
    motifs = c("A", "C", "G"),
    sparse_row = c(0L, 0L),
    sparse_motif = c(0L, 2L),
    sparse_frequency = c(1 / 3, 2 / 3),
    group_labels = list("alpha", "beta"),
    sparse_shape = c(2L, 3L),
    row_scaling_factor = c(6, 0),
    eligible_windows = c(2L, 0L),
    blacklisted_fraction = c(0.125, 0)
  )
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  corrected <- end_motif_data_frame(ends, ref_kmers = ref_kmers)

  expect_equal(corrected$reference_motif, c("A", "C", "G", "A", "C", "G"))
  expect_equal(corrected$reference_frequency, c(1 / 3, 0, 2 / 3, 0, 0, 0))
  expect_equal(corrected$correction_motif_count, c(2L, 2L, 2L, 0L, 0L, 0L))
  expect_equal(corrected$reference_scale, c(2 / 3, 0, 4 / 3, 0, 0, 0))
  expect_equal(corrected$reference_corrected_count, c(1.5, 0, 3.75, 0, 0, 0))
})

test_that("reference correction maps selected grouped rows by group name", {
  end_path <- make_dense_grouped_end_motif_zarr_fixture()
  ref_path <- make_sparse_grouped_ref_kmer_zarr_fixture(
    motifs = c("A", "C", "G"),
    sparse_row = c(1L, 1L),
    sparse_motif = c(0L, 2L),
    sparse_frequency = c(0.5, 0.5),
    group_labels = list("beta", "alpha"),
    sparse_shape = c(2L, 3L),
    row_scaling_factor = c(0, 2),
    eligible_windows = c(0L, 2L),
    blacklisted_fraction = c(0, 0.125)
  )
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  corrected <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    group_idxs = 1L
  )

  expect_equal(corrected$group_name, c("alpha", "alpha", "alpha"))
  expect_equal(corrected$reference_frequency, c(0.5, 0, 0.5))
  expect_equal(corrected$correction_motif_count, c(2L, 2L, 2L))
  expect_equal(corrected$reference_scale, c(1, 0, 1))
  expect_equal(corrected$reference_corrected_count, c(1, 0, 5))
})

test_that("reference correction rejects missing reference motifs by default", {
  end_path <- make_dense_windowed_end_motif_zarr_fixture()
  ref_path <- make_dense_windowed_ref_kmer_zarr_fixture(
    motifs = c("A", "G"),
    frequencies = matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2L, byrow = TRUE),
    row_scaling_factor = c(4, 4),
    chromosome_names = c("chr1", "chr2"),
    row_chromosome = c(0L, 1L),
    row_start_bp = c(10L, 30L),
    row_end_bp = c(12L, 36L),
    blacklisted_fraction = c(0, 0.125)
  )
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  expect_error(
    end_motif_data_frame(ends, ref_kmers = ref_kmers),
    "unsupported_motifs = \"drop\"",
    fixed = TRUE
  )
})

test_that("reference correction can drop unsupported motifs", {
  end_path <- make_dense_windowed_end_motif_zarr_fixture()
  ref_path <- make_dense_windowed_ref_kmer_zarr_fixture(
    motifs = c("A", "G"),
    frequencies = matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2L, byrow = TRUE),
    row_scaling_factor = c(4, 4),
    chromosome_names = c("chr1", "chr2"),
    row_chromosome = c(0L, 1L),
    row_start_bp = c(10L, 30L),
    row_end_bp = c(12L, 36L),
    blacklisted_fraction = c(0, 0.125)
  )
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  corrected <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    unsupported_motifs = "drop"
  )

  expect_equal(corrected$reference_motif, c("A", "G", "A", "G"))
  expect_equal(corrected$correction_motif_count, c(2L, 2L, 2L, 2L))
  expect_equal(corrected$reference_scale, c(1, 1, 1, 1))
  expect_equal(corrected$reference_corrected_count, c(0, 2, 1.5, 0))
})

test_that("reference-corrected matrix extractors reject drop policy", {
  end_path <- make_dense_windowed_end_motif_zarr_fixture()
  ref_path <- make_dense_windowed_ref_kmer_zarr_fixture(
    motifs = c("A", "G"),
    frequencies = matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2L, byrow = TRUE),
    row_scaling_factor = c(4, 4),
    chromosome_names = c("chr1", "chr2"),
    row_chromosome = c(0L, 1L),
    row_start_bp = c(10L, 30L),
    row_end_bp = c(12L, 36L),
    blacklisted_fraction = c(0, 0.125)
  )
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  expect_error(
    dense_corrected_counts_matrix(ends, ref_kmers, unsupported_motifs = "drop"),
    "fixed-shape dense_corrected_counts_matrix",
    fixed = TRUE
  )
  expect_error(
    sparse_corrected_counts_matrix(ends, ref_kmers, unsupported_motifs = "drop"),
    "fixed-shape sparse_corrected_counts_matrix",
    fixed = TRUE
  )
})

test_that("reference correction can keep unsupported motifs as NA", {
  end_path <- make_dense_windowed_end_motif_zarr_fixture()
  ref_path <- make_dense_windowed_ref_kmer_zarr_fixture(
    motifs = c("A", "G"),
    frequencies = matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2L, byrow = TRUE),
    row_scaling_factor = c(4, 4),
    chromosome_names = c("chr1", "chr2"),
    row_chromosome = c(0L, 1L),
    row_start_bp = c(10L, 30L),
    row_end_bp = c(12L, 36L),
    blacklisted_fraction = c(0, 0.125)
  )
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  corrected <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    unsupported_motifs = "keep_na"
  )

  expect_equal(corrected$correction_motif_count, c(2L, 2L, 2L, 2L, 2L, 2L))
  expect_equal(corrected$reference_scale, c(1, 1, 0, 1, 1, 0))
  expect_equal(corrected$reference_corrected_count, c(0, 2, 0, 1.5, 0, NA_real_))
})

test_that("reference correction requires opt-in for global reference bias", {
  end_path <- make_dense_windowed_end_motif_zarr_fixture()
  ref_path <- make_dense_global_ref_kmer_zarr_fixture()
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  expect_error(
    end_motif_data_frame(ends, ref_kmers = ref_kmers),
    "use_global_bias = TRUE",
    fixed = TRUE
  )
})

test_that("reference correction rejects global-bias flag for matched reference rows", {
  end_path <- make_dense_windowed_end_motif_zarr_fixture()
  ref_path <- make_dense_windowed_ref_kmer_zarr_fixture()
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  expect_error(
    end_motif_data_frame(ends, ref_kmers = ref_kmers, use_global_bias = TRUE),
    "use_global_bias = TRUE requires a global reference k-mer output",
    fixed = TRUE
  )
})

test_that("reference correction can use global reference bias", {
  end_path <- make_dense_windowed_end_motif_zarr_fixture()
  ref_path <- make_dense_global_ref_kmer_zarr_fixture()
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  corrected <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    use_global_bias = TRUE
  )
  selected <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    window_idxs = 2L,
    motifs = "_G",
    use_global_bias = TRUE
  )
  filtered <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    max_blacklisted_fraction = 0.1,
    use_global_bias = TRUE
  )

  expect_equal(corrected$correction_motif_count, c(3L, 3L, 3L, 3L, 3L, 3L))
  expect_equal(corrected$reference_scale, c(1, 0.5, 1.5, 1, 0.5, 1.5))
  expect_equal(
    corrected$reference_corrected_count,
    c(0, 4, 0, 1.5, 0, 8 / 3)
  )
  expected <- corrected[
    corrected$window_idx == 2L & corrected$motif == "_G",
    ,
    drop = FALSE
  ]
  row.names(expected) <- NULL
  expect_equal(selected, expected)
  expected_filtered <- corrected[corrected$window_idx == 1L, , drop = FALSE]
  row.names(expected_filtered) <- NULL
  expect_equal(filtered, expected_filtered)
  expect_error(
    end_motif_data_frame(
      ends,
      ref_kmers = ref_kmers,
      groups = "alpha",
      use_global_bias = TRUE
    ),
    "Unused argument(s): groups",
    fixed = TRUE
  )
})
