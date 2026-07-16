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

  expect_true("corrected_count" %in% names(corrected))
  expect_true("corrected_frequency" %in% names(corrected))
  expect_equal(
    corrected$corrected_count,
    c(0, 4, 0, 1, 0, 16 / 3)
  )
})

test_that("reference correction rejects non-finite corrected counts", {
  smallest_positive_double <- 2^-1074
  end_path <- make_dense_windowed_end_motif_zarr_fixture()
  ref_path <- make_dense_windowed_ref_kmer_zarr_fixture(
    motifs = c("A", "G", "T"),
    frequencies = matrix(
      c(
        0.5,
        smallest_positive_double,
        0.5,
        0.5,
        smallest_positive_double,
        0.5
      ),
      nrow = 2L,
      byrow = TRUE
    ),
    chromosome_names = c("chr1", "chr2"),
    row_chromosome = c(0L, 1L),
    row_start_bp = c(10L, 30L),
    row_end_bp = c(12L, 36L),
    blacklisted_fraction = c(0, 0.125),
    row_scaling_factor = c(3, 3)
  )
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  expect_error(
    end_motif_data_frame(ends, ref_kmers = ref_kmers),
    "non-finite corrected counts for motifs: _G",
    fixed = TRUE
  )
})

test_that("empty sparse selections still validate use_global_bias", {
  end_path <- make_dense_windowed_end_motif_zarr_fixture()
  ref_path <- make_dense_windowed_ref_kmer_zarr_fixture(
    motifs = c("A", "G", "T"),
    chromosome_names = c("chr1", "chr2"),
    row_chromosome = c(0L, 1L),
    row_start_bp = c(10L, 30L),
    row_end_bp = c(12L, 36L),
    blacklisted_fraction = c(0, 0.125)
  )
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  expect_error(
    sparse_corrected_counts_matrix(
      ends,
      ref_kmers,
      motifs = character(),
      use_global_bias = "yes"
    ),
    "use_global_bias must be TRUE or FALSE",
    fixed = TRUE
  )
})

test_that("empty sparse selections still validate matched rows", {
  end_path <- make_dense_windowed_end_motif_zarr_fixture()
  ref_path <- make_dense_windowed_ref_kmer_zarr_fixture(
    motifs = c("A", "G", "T"),
    chromosome_names = c("chr1", "chr2"),
    row_chromosome = c(0L, 1L),
    row_start_bp = c(11L, 30L),
    row_end_bp = c(12L, 36L),
    blacklisted_fraction = c(0, 0.125)
  )
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  expect_error(
    sparse_corrected_counts_matrix(
      ends,
      ref_kmers,
      motifs = character()
    ),
    "End-motif and reference k-mer rows do not match",
    fixed = TRUE
  )
})

test_that("two-sided correction requires reference k-mers", {
  end_path <- make_dense_global_end_motif_zarr_fixture()
  ends <- read_end_motifs(end_path)

  expect_error(
    end_motif_data_frame(ends, two_sided_correction = "joint"),
    "two_sided_correction requires ref_kmers",
    fixed = TRUE
  )
})

test_that("motif-group reference correction rejects two-sided mode", {
  end_path <- make_dense_global_end_motif_group_zarr_fixture()
  ref_path <- make_dense_global_ref_kmer_group_zarr_fixture()
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  for (mode in c("joint", "split", "outside", "inside")) {
    expect_error(
      end_motif_data_frame(
        ends,
        ref_kmers = ref_kmers,
        two_sided_correction = mode
      ),
      "Motif-group",
      fixed = TRUE
    )
  }
})

test_that("reference correction rejects canonical reference", {
  end_path <- make_dense_windowed_end_motif_zarr_fixture()
  ref_path <- make_dense_windowed_ref_kmer_zarr_fixture(
    motifs = c("A", "C"),
    canonical = TRUE,
    frequencies = matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2L, byrow = TRUE),
    row_scaling_factor = c(2, 2)
  )
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  expect_error(
    end_motif_data_frame(ends, ref_kmers = ref_kmers),
    "non-canonical",
    fixed = TRUE
  )
})

test_that("one-sided reference correction rejects explicit mode", {
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

  expect_error(
    end_motif_data_frame(
      ends,
      ref_kmers = ref_kmers,
      two_sided_correction = "joint"
    ),
    "One-sided",
    fixed = TRUE
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

  expect_equal(corrected$corrected_count, c(1, 0, 1.25, 0))
})

test_that("two-sided reference correction modes use the mode axis", {
  # This fixture is shared with the Rust and Python correction tests.
  # Every joint and side reference frequency is positive and differs from its
  # uniform baseline, so every correction mode has a visible effect.
  end_path <- make_dense_global_end_motif_zarr_fixture(
    motifs = c("A_C", "A_G", "T_C", "T_G"),
    counts = matrix(c(2, 4, 6, 8), nrow = 1L)
  )
  ref_path <- make_dense_global_ref_kmer_zarr_fixture(
    motifs = c("AC", "AG", "TC", "TG"),
    frequencies = matrix(c(1 / 8, 1 / 8, 1 / 4, 1 / 2), nrow = 1L),
    row_scaling_factor = 4
  )
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  expect_error(
    end_motif_data_frame(ends, ref_kmers = ref_kmers),
    "two-sided",
    fixed = TRUE
  )

  joint <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    two_sided_correction = "joint"
  )
  # Four positive reference motifs make the uniform frequency 1/4. Relative
  # to uniform, frequencies [1/8, 1/8, 1/4, 1/2] give correction factors
  # [1/2, 1/2, 1, 2] for [AC, AG, TC, TG]. Dividing original counts
  # [2, 4, 6, 8] by those factors gives [4, 8, 6, 4]. Their total is 22, so
  # dividing each corrected count by 22 gives [2/11, 4/11, 3/11, 2/11].
  expect_identical(joint$motif_idx, seq_len(4L))
  expect_equal(joint$motif, c("A_C", "A_G", "T_C", "T_G"))
  expect_equal(joint$count, c(2, 4, 6, 8))
  expect_equal(joint$corrected_count, c(4, 8, 6, 4))
  expect_equal(joint$corrected_frequency, c(2 / 11, 4 / 11, 3 / 11, 2 / 11))

  split <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    two_sided_correction = "split"
  )
  expect_identical(split$motif_idx, seq_len(4L))
  expect_equal(split$motif, c("A_C", "A_G", "T_C", "T_G"))
  expect_true("corrected_count" %in% names(split))
  expect_true("corrected_frequency" %in% names(split))
  expect_false("reference_frequency" %in% names(split))
  expect_equal(split$count, c(2, 4, 6, 8))
  # Two positive labels on each side make each side's uniform frequency 1/2.
  # Outside frequencies A=1/4 and T=3/4 give factors 1/2 and 3/2. Inside
  # frequencies C=3/8 and G=5/8 give factors 3/4 and 5/4. Multiplying matching
  # side factors gives [3/8, 5/8, 9/8, 15/8] for [A_C, A_G, T_C, T_G].
  # Dividing original counts [2, 4, 6, 8] by those factors gives
  # [16/3, 32/5, 16/3, 64/15]. The corrected counts total
  # 64/3, so normalization gives [1/4, 3/10, 1/4, 1/5].
  expect_equal(split$corrected_count, c(16 / 3, 32 / 5, 16 / 3, 64 / 15))
  expect_equal(
    split$corrected_frequency,
    c(1 / 4, 3 / 10, 1 / 4, 1 / 5)
  )

  outside <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    two_sided_correction = "outside"
  )
  expect_identical(outside$motif_idx, seq_len(2L))
  expect_equal(outside$motif, c("A_", "T_"))
  # Counts aggregate to [6, 14]. Relative to the uniform outside frequency
  # 1/2, reference frequencies [1/4, 3/4] give factors [1/2, 3/2]. Dividing
  # the aggregated counts by them gives [12, 28/3]. These total 64/3, so
  # normalization gives frequencies [9/16, 7/16].
  expect_equal(outside$count, c(6, 14))
  expect_equal(outside$corrected_count, c(12, 28 / 3))
  expect_equal(outside$corrected_frequency, c(9 / 16, 7 / 16))
  outside_dense <- dense_corrected_counts_matrix(
    ends,
    ref_kmers,
    two_sided_correction = "outside"
  )
  expect_equal(unname(outside_dense), matrix(c(12, 28 / 3), nrow = 1L))
  expect_identical(colnames(outside_dense), c("A_", "T_"))
  outside_sparse <- sparse_corrected_counts_matrix(
    ends,
    ref_kmers,
    two_sided_correction = "outside"
  )
  expect_equal(unname(as.matrix(outside_sparse)), matrix(c(12, 28 / 3), nrow = 1L))
  expect_identical(colnames(outside_sparse), c("A_", "T_"))

  inside <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    two_sided_correction = "inside"
  )
  # Counts aggregate to [8, 12]. Relative to the uniform inside frequency
  # 1/2, reference frequencies [3/8, 5/8] give factors [3/4, 5/4]. Dividing
  # the aggregated counts by them gives [32/3, 48/5]. These total 304/15, so
  # normalization gives frequencies [10/19, 9/19].
  expect_identical(inside$motif_idx, seq_len(2L))
  expect_equal(inside$motif, c("_C", "_G"))
  expect_equal(inside$count, c(8, 12))
  expect_equal(inside$corrected_count, c(32 / 3, 48 / 5))
  expect_equal(inside$corrected_frequency, c(10 / 19, 9 / 19))
  selected_inside <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    two_sided_correction = "inside",
    motifs = "_G"
  )
  expected_inside <- inside[inside$motif == "_G", , drop = FALSE]
  row.names(expected_inside) <- NULL
  expect_equal(selected_inside, expected_inside, ignore_attr = TRUE)
  expect_equal(selected_inside$count, 12)
  expect_equal(selected_inside$corrected_count, 48 / 5)
  expect_equal(selected_inside$corrected_frequency, 9 / 19)
  expect_error(
    end_motif_data_frame(
      ends,
      ref_kmers = ref_kmers,
      two_sided_correction = "outside",
      motif_idxs = 1L
    ),
    "motif index selectors",
    fixed = TRUE
  )
  expect_error(
    end_motif_data_frame(
      ends,
      ref_kmers = ref_kmers,
      two_sided_correction = "outside",
      motifs = "A_C"
    ),
    "Side-mode motif axis",
    fixed = TRUE
  )
})

test_that("side correction applies unsupported policy after aggregation", {
  end_path <- make_dense_global_end_motif_zarr_fixture(
    motifs = c("A_C", "A_G", "T_C", "T_G"),
    counts = matrix(c(2, 4, 6, 8), nrow = 1L)
  )
  ref_path <- make_dense_global_ref_kmer_zarr_fixture(
    motifs = "AC",
    frequencies = matrix(1, nrow = 1L),
    row_scaling_factor = 1
  )
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  expect_error(
    end_motif_data_frame(
      ends,
      ref_kmers = ref_kmers,
      two_sided_correction = "inside"
    ),
    "_G",
    fixed = TRUE
  )

  dropped <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    two_sided_correction = "inside",
    unsupported_motifs = "drop"
  )
  expect_identical(dropped$motif, "_C")
  expect_equal(dropped$corrected_count, 8)
  expect_equal(dropped$corrected_frequency, 1)

  kept <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    two_sided_correction = "inside",
    unsupported_motifs = "keep_na"
  )
  expect_identical(kept$motif, c("_C", "_G"))
  expect_equal(kept$corrected_count, c(8, NA_real_))
  expect_true(all(is.na(kept$corrected_frequency)))
})

test_that("corrected frequencies are zero when corrected total is zero", {
  end_path <- make_dense_global_end_motif_zarr_fixture(
    motifs = c("A_C", "A_G", "T_C", "T_G"),
    counts = matrix(0, nrow = 1L, ncol = 4L)
  )
  ref_path <- make_dense_global_ref_kmer_zarr_fixture(
    motifs = c("AC", "AG", "TC", "TG"),
    frequencies = matrix(0.25, nrow = 1L, ncol = 4L),
    row_scaling_factor = 4
  )
  ends <- read_end_motifs(end_path)
  ref_kmers <- read_ref_kmers(ref_path)

  corrected <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    two_sided_correction = "split"
  )

  expect_equal(corrected$corrected_count, rep(0, 4L))
  expect_equal(corrected$corrected_frequency, rep(0, 4L))
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
  expect_equal(selected$corrected_count, 1)
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
  expect_equal(corrected$corrected_count, c(0, 4, 0))
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
    "Positive-count end motifs have no positive reference-based correction factor",
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

  expect_equal(corrected$corrected_count, c(1.5, 0, 3.75, 0, 0, 0))
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
  expect_equal(corrected$corrected_count, c(1, 0, 5))
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

  expect_equal(corrected$corrected_count, c(0, 2, 1.5, 0))
  expect_equal(corrected$corrected_frequency, c(0, 1, 1, 0))
  expect_false("_C" %in% corrected$motif)
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

  expect_equal(corrected$corrected_count, c(0, 2, 0, 1.5, 0, NA_real_))
  expect_equal(
    corrected$corrected_frequency,
    c(0, 1, 0, NA_real_, NA_real_, NA_real_)
  )
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

  expect_equal(
    corrected$corrected_count,
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
