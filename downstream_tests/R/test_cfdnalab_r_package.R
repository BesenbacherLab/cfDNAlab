source(file.path("downstream_tests", "R", "common.R"))

library(cfdnalab)
library(testthat)

test_that("R helper package reads midpoint profiles", {
  midpoints <- read_midpoints(midpoint_zarr_path())

  expect_identical(schema_version(midpoints), 1L)
  expect_identical(
    group_metadata(midpoints)$group_name,
    c("LYL1", "beta-site", "gamma_long")
  )
  expect_identical(
    length_bins(midpoints)$length_start_bp,
    c(30L, 50L, 70L)
  )
  expect_identical(
    positions(midpoints)$position_bin_start_bp,
    c(0L, 2L, 4L, 6L, 8L)
  )

  beta_profile <- midpoint_data_frame(midpoints, groups = "beta-site", length_bin_idxs = 2L)
  lyl1_range <- midpoint_data_frame(midpoints, groups = "LYL1", with_length_range = c(50L, 100L))
  expect_identical(beta_profile$position_idx, c(1L, 2L, 3L, 4L, 5L))
  expect_equal(beta_profile$count, c(0, 0, 1.5, 0, 0.5))
  expect_equal(lyl1_range$length_bin_idx, c(rep(2L, 5), rep(3L, 5)))
  expect_equal(lyl1_range$count, c(0, 0, 1.5, 0.5, 0, 0, 0, 0, 0, 2))
  expect_equal(
    profile_array(midpoints, group_idx = 3L, length_bin_idx = 1L),
    c(0, 0, 2.5, 0, 0)
  )
})

test_that("R helper package reads dense global end motifs", {
  dense_global <- read_end_motifs(dense_global_end_zarr_path())

  expect_identical(schema_version(dense_global), 2L)
  expect_identical(storage_mode(dense_global), "dense")
  expect_identical(row_mode(dense_global), "global")
  expect_identical(motifs(dense_global)$motif, c("_A", "_C", "_G", "_T"))
  expect_identical(has_motif(dense_global, "_A"), TRUE)
  expect_identical(has_motif(dense_global, "_AA"), FALSE)
  expect_equal(
    unname(dense_counts_vector(dense_global)),
    c(1, 0, 1, 0)
  )
  expect_equal(
    end_motif_data_frame(dense_global)$count,
    c(1, 0, 1, 0)
  )
  expect_error(
    end_motif_data_frame(dense_global, two_sided_correction = "joint"),
    "two_sided_correction requires ref_kmers",
    fixed = TRUE
  )
})

test_that("R helper package reads sparse windowed end motifs", {
  sparse_windowed <- read_end_motifs(sparse_windowed_end_zarr_path())

  expect_identical(storage_mode(sparse_windowed), "sparse_coo")
  expect_identical(row_mode(sparse_windowed), "bed")
  expect_identical(window_metadata(sparse_windowed)$window_idx, c(1L, 2L, 3L))
  expect_identical(window_metadata(sparse_windowed)$chrom, c("chr1", "chr1", "chr2"))
  expect_identical(window_metadata(sparse_windowed)$start, c(10L, 19L, 10L))
  expect_identical(window_metadata(sparse_windowed)$end, c(11L, 20L, 11L))
  expect_equal(
    as.matrix(sparse_counts_matrix(sparse_windowed)),
    matrix(c(0, 1, 1, 0, 0, 1), nrow = 3, byrow = TRUE)
  )
  expect_equal(
    end_motif_data_frame(
      sparse_windowed,
      motifs = "_A",
      densify = TRUE,
      max_blacklisted_fraction = 0
    )$count,
    c(0, 1, 0)
  )
  expect_equal(end_motif_data_frame(sparse_windowed, window_idxs = 1L)$count, 1)
  expect_equal(end_motif_data_frame(sparse_windowed, motifs = "_G")$count, c(1, 1))

  ordered_dense <- end_motif_data_frame(
    sparse_windowed,
    window_idxs = c(2L, 1L),
    motifs = c("_G", "_A"),
    densify = TRUE
  )
  expect_equal(ordered_dense$window_idx, c(2L, 2L, 1L, 1L))
  expect_equal(ordered_dense$motif, c("_G", "_A", "_G", "_A"))
  expect_equal(ordered_dense$count, c(0, 1, 1, 0))
})

test_that("R helper package reads sparse windowed selected motif-file end motifs", {
  selected <- read_end_motifs(sparse_windowed_selected_motifs_end_zarr_path())

  expect_identical(storage_mode(selected), "sparse_coo")
  expect_identical(row_mode(selected), "bed")
  expect_identical(motifs(selected)$motif, c("GT_AC", "AC_GT"))
  expect_false(has_motif(selected, "TT_TT"))
  expect_equal(
    as.matrix(sparse_counts_matrix(selected)),
    matrix(c(0, 1, 1, 0, 0, 1), nrow = 3, byrow = TRUE)
  )
  expect_equal(
    as.matrix(sparse_counts_matrix(selected, motifs = c("AC_GT", "GT_AC"))),
    matrix(c(1, 0, 0, 1, 1, 0), nrow = 3, byrow = TRUE)
  )
  expect_error(
    end_motif_data_frame(selected, motifs = "TT_TT", densify = TRUE),
    "Unknown end-motif label"
  )
})

test_that("R helper package corrects two-sided end motifs without same motifs file", {
  ends <- read_end_motifs(sparse_windowed_two_sided_end_zarr_path())
  ref_kmers <- read_ref_kmers(sparse_windowed_end_motif_ref_kmer_zarr_path())

  expect_identical(motifs(ends)$motif, c("AC_GT", "GT_AC"))
  expect_identical(kmer_size(ref_kmers), 4L)
  expect_identical(
    motifs(ref_kmers)$motif,
    c("AAAA", "AAAC", "AACG", "ACGT", "CGTA", "CGTT", "GTAC", "GTTT", "TACG", "TTTT")
  )
  expected_counts <- matrix(c(1, 0, 0, 1, 1, 0), nrow = 3L, byrow = TRUE)
  expected_count_vector <- c(1, 0, 0, 1, 1, 0)
  # The stored columns are [AC_GT, GT_AC], and the three sample rows contain
  # [[1, 0], [0, 1], [1, 0]]. Ten positive reference k-mers make the joint
  # uniform frequency 1/10. Relative to uniform, ACGT frequency 1/4 gives
  # correction factor (1/4)/(1/10) = 5/2, while GTAC frequency 3/20 gives
  # (3/20)/(1/10) = 3/2. Dividing each observed count by its factor gives
  # [[2/5, 0], [0, 2/3], [2/5, 0]].
  expected_joint <- matrix(c(2 / 5, 0, 0, 2 / 3, 2 / 5, 0), nrow = 3L, byrow = TRUE)
  # Six positive labels on each side make the side-wise uniform frequency
  # 1/6. Outside AC has frequency 1/4, giving factor 3/2, while outside GT
  # has frequency 1/5, giving factor 6/5. Inside GT and AC have the same
  # respective frequencies and factors. Split therefore divides AC_GT by
  # (3/2)*(3/2)=9/4 and GT_AC by (6/5)*(6/5)=36/25. Because each observed
  # count is 1, the corrected values are 4/9 and 25/36. Outside correction
  # divides AC_GT by 3/2 and GT_AC by 6/5. Inside correction uses the same
  # factors. In stored row and column order, the split, outside, and inside
  # matrices are therefore [[4/9, 0], [0, 25/36], [4/9, 0]],
  # [[2/3, 0], [0, 5/6], [2/3, 0]], and
  # [[2/3, 0], [0, 5/6], [2/3, 0]], respectively.
  expected_split <- matrix(c(4 / 9, 0, 0, 25 / 36, 4 / 9, 0), nrow = 3L, byrow = TRUE)
  expected_outside <- matrix(c(2 / 3, 0, 0, 5 / 6, 2 / 3, 0), nrow = 3L, byrow = TRUE)
  expected_inside <- matrix(c(2 / 3, 0, 0, 5 / 6, 2 / 3, 0), nrow = 3L, byrow = TRUE)

  expect_error(
    end_motif_data_frame(ends, ref_kmers = ref_kmers, densify = TRUE),
    "two-sided",
    fixed = TRUE
  )

  joint <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    densify = TRUE,
    two_sided_correction = "joint"
  )
  expect_equal(joint$window_idx, c(1L, 1L, 2L, 2L, 3L, 3L))
  expect_equal(joint$motif, rep(c("AC_GT", "GT_AC"), 3L))
  expect_equal(joint$count, expected_count_vector)
  expect_equal(joint$corrected_count, as.vector(t(expected_joint)))
  expect_equal(joint$corrected_frequency, expected_count_vector)

  expect_equal(
    dense_corrected_counts_matrix(
      ends,
      ref_kmers,
      allow_densify = TRUE,
      two_sided_correction = "split"
    ),
    expected_split
  )
  expect_equal(
    as.matrix(sparse_corrected_counts_matrix(
      ends,
      ref_kmers,
      two_sided_correction = "split"
    )),
    expected_split
  )

  outside <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    densify = TRUE,
    two_sided_correction = "outside"
  )
  expect_equal(outside$motif, rep(c("AC_", "GT_"), 3L))
  expect_equal(outside$count, expected_count_vector)
  expect_equal(outside$corrected_count, as.vector(t(expected_outside)))
  expect_equal(
    unname(dense_corrected_counts_matrix(
      ends,
      ref_kmers,
      allow_densify = TRUE,
      two_sided_correction = "outside"
    )),
    expected_outside
  )

  inside <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    densify = TRUE,
    two_sided_correction = "inside"
  )
  expect_equal(inside$motif, rep(c("_GT", "_AC"), 3L))
  expect_equal(inside$count, expected_count_vector)
  expect_equal(inside$corrected_count, as.vector(t(expected_inside)))
  expect_equal(
    unname(as.matrix(sparse_corrected_counts_matrix(
      ends,
      ref_kmers,
      two_sided_correction = "inside"
    ))),
    expected_inside
  )
})

test_that("R helper package corrects two-sided end motifs with same motifs file", {
  ends <- read_end_motifs(sparse_windowed_selected_motifs_end_zarr_path())
  ref_kmers <- read_ref_kmers(sparse_windowed_selected_end_motifs_ref_kmer_zarr_path())

  expect_identical(motifs(ends)$motif, c("GT_AC", "AC_GT"))
  expect_identical(motifs(ref_kmers)$motif, c("GTAC", "ACGT", "GTTT", "TTTT"))
  expected_counts <- matrix(c(0, 1, 1, 0, 0, 1), nrow = 3L, byrow = TRUE)
  expected_count_vector <- c(0, 1, 1, 0, 0, 1)
  # The stored columns are [GT_AC, AC_GT], and the three sample rows contain
  # [[0, 1], [1, 0], [0, 1]]. Four positive reference k-mers make the joint
  # uniform frequency 1/4. GTAC and ACGT frequencies 6/19 and 10/19 give
  # factors (6/19)/(1/4) = 24/19 and (10/19)/(1/4) = 40/19. Dividing each
  # observed count by its factor gives
  # [[0, 19/40], [19/24, 0], [0, 19/40]].
  expected_joint <- matrix(c(0, 19 / 40, 19 / 24, 0, 0, 19 / 40), nrow = 3L, byrow = TRUE)
  # Three positive labels on each side make each side's uniform frequency
  # 1/3. Outside AC and GT have frequencies 10/19 and 8/19, giving factors
  # 30/19 and 24/19. Inside GT and AC have frequencies 10/19 and 6/19,
  # giving factors 30/19 and 18/19. Split therefore divides AC_GT by
  # (30/19)*(30/19)=900/361 and GT_AC by
  # (24/19)*(18/19)=432/361. Because each observed count is 1, the corrected
  # values are 361/900 and 361/432. Outside correction divides AC_GT by
  # 30/19 and GT_AC by 24/19. Inside correction divides AC_GT by 30/19 and
  # GT_AC by 18/19. In stored order, the split, outside, and inside matrices
  # are therefore [[0, 361/900], [361/432, 0], [0, 361/900]],
  # [[0, 19/30], [19/24, 0], [0, 19/30]], and
  # [[0, 19/30], [19/18, 0], [0, 19/30]], respectively.
  expected_split <- matrix(c(0, 361 / 900, 361 / 432, 0, 0, 361 / 900), nrow = 3L, byrow = TRUE)
  expected_outside <- matrix(c(0, 19 / 30, 19 / 24, 0, 0, 19 / 30), nrow = 3L, byrow = TRUE)
  expected_inside <- matrix(c(0, 19 / 30, 19 / 18, 0, 0, 19 / 30), nrow = 3L, byrow = TRUE)

  joint <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    densify = TRUE,
    two_sided_correction = "joint"
  )
  expect_equal(joint$motif, rep(c("GT_AC", "AC_GT"), 3L))
  expect_equal(joint$count, expected_count_vector)
  expect_equal(joint$corrected_count, as.vector(t(expected_joint)))
  expect_equal(joint$corrected_frequency, expected_count_vector)

  expect_equal(
    dense_corrected_counts_matrix(
      ends,
      ref_kmers,
      allow_densify = TRUE,
      two_sided_correction = "split"
    ),
    expected_split
  )
  expect_equal(
    as.matrix(sparse_corrected_counts_matrix(
      ends,
      ref_kmers,
      two_sided_correction = "split"
    )),
    expected_split
  )

  outside <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    densify = TRUE,
    two_sided_correction = "outside"
  )
  expect_equal(outside$motif, rep(c("GT_", "AC_"), 3L))
  expect_equal(outside$count, expected_count_vector)
  expect_equal(outside$corrected_count, as.vector(t(expected_outside)))
  expect_equal(
    unname(dense_corrected_counts_matrix(
      ends,
      ref_kmers,
      allow_densify = TRUE,
      two_sided_correction = "outside"
    )),
    expected_outside
  )

  selected_outside <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    motifs = "AC_",
    densify = TRUE,
    two_sided_correction = "outside"
  )
  expect_equal(selected_outside$motif, c("AC_", "AC_", "AC_"))
  expect_equal(selected_outside$corrected_count, c(19 / 30, 0, 19 / 30))

  inside <- end_motif_data_frame(
    ends,
    ref_kmers = ref_kmers,
    densify = TRUE,
    two_sided_correction = "inside"
  )
  expect_equal(inside$motif, rep(c("_AC", "_GT"), 3L))
  expect_equal(inside$count, expected_count_vector)
  expect_equal(inside$corrected_count, as.vector(t(expected_inside)))
  expect_equal(
    unname(as.matrix(sparse_corrected_counts_matrix(
      ends,
      ref_kmers,
      two_sided_correction = "inside"
    ))),
    expected_inside
  )
})

test_that("R helper package reads sparse grouped end motifs", {
  sparse_grouped <- read_end_motifs(sparse_grouped_end_zarr_path())

  expect_identical(storage_mode(sparse_grouped), "sparse_coo")
  expect_identical(row_mode(sparse_grouped), "grouped_bed")
  expect_identical(group_idx(sparse_grouped, "alpha"), 2L)
  expect_identical(
    group_metadata(sparse_grouped)$group_name,
    c("beta", "alpha", "gamma")
  )
  expect_equal(
    as.matrix(sparse_counts_matrix(sparse_grouped)),
    matrix(c(1, 2, 1, 0, 0, 0), nrow = 3, byrow = TRUE)
  )
  expect_equal(end_motif_data_frame(sparse_grouped, groups = "beta")$count, c(1, 2))
  expect_equal(
    end_motif_data_frame(
      sparse_grouped,
      groups = "beta",
      densify = TRUE,
      max_blacklisted_fraction = 0
    )$count,
    c(1, 2)
  )

  expect_equal(nrow(end_motif_data_frame(sparse_grouped, groups = "gamma")), 0L)
  gamma_dense <- end_motif_data_frame(sparse_grouped, groups = "gamma", densify = TRUE)
  expect_equal(gamma_dense$group_name, c("gamma", "gamma"))
  expect_equal(gamma_dense$motif, c("_A", "_G"))
  expect_equal(gamma_dense$count, c(0, 0))
})

test_that("R helper package reads sparse grouped motif-group end motifs", {
  motif_grouped <- read_end_motifs(sparse_grouped_motif_group_end_zarr_path())

  expect_identical(schema_version(motif_grouped), 2L)
  expect_identical(storage_mode(motif_grouped), "sparse_coo")
  expect_identical(row_mode(motif_grouped), "grouped_bed")
  expect_equal(
    motifs(motif_grouped),
    data.frame(
      motif_idx = c(1L, 2L),
      motif = c("left-hit", "right-hit"),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_equal(motif_idx(motif_grouped, "right-hit"), 2L)
  expect_true(has_motif(motif_grouped, "left-hit"))
  expect_equal(
    as.matrix(sparse_counts_matrix(motif_grouped)),
    matrix(c(2, 1, 0, 1, 0, 0), nrow = 3, byrow = TRUE)
  )
  expect_equal(
    as.matrix(sparse_counts_matrix(
      motif_grouped,
      groups = c("alpha", "beta"),
      motifs = c("right-hit", "left-hit")
    )),
    matrix(c(1, 0, 1, 2), nrow = 2, byrow = TRUE)
  )

  rows <- end_motif_data_frame(motif_grouped)
  expect_equal(names(rows), c(
    "group_idx",
    "group_name",
    "eligible_windows",
    "blacklisted_fraction",
    "motif_idx",
    "motif",
    "count"
  ))
  expect_equal(rows$group_name, c("beta", "beta", "alpha"))
  expect_equal(rows$motif_idx, c(1L, 2L, 2L))
  expect_equal(rows$motif, c("left-hit", "right-hit", "right-hit"))
  expect_equal(rows$count, c(2, 1, 1))

  alpha_dense <- end_motif_data_frame(motif_grouped, groups = "alpha", densify = TRUE)
  expect_equal(alpha_dense$motif, c("left-hit", "right-hit"))
  expect_equal(alpha_dense$count, c(0, 1))
  expect_error(
    end_motif_data_frame(motif_grouped, motifs = "_A"),
    "Unknown end-motif label",
    fixed = TRUE
  )
})

test_that("R helper package reads sparse grouped wide motif-group end motifs", {
  motif_grouped <- read_end_motifs(sparse_grouped_wide_motif_group_end_zarr_path())

  expect_identical(schema_version(motif_grouped), 2L)
  expect_identical(storage_mode(motif_grouped), "sparse_coo")
  expect_identical(row_mode(motif_grouped), "grouped_bed")
  expect_equal(
    motifs(motif_grouped),
    data.frame(
      motif_idx = c(1L, 2L),
      motif = c("left-hit-wide", "right-hit-wide"),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_equal(motif_idx(motif_grouped, "left-hit-wide"), 1L)
  expect_equal(
    as.matrix(sparse_counts_matrix(motif_grouped)),
    matrix(c(2, 1, 0, 1, 0, 0), nrow = 3, byrow = TRUE)
  )
  expect_equal(
    as.matrix(sparse_counts_matrix(
      motif_grouped,
      groups = c("alpha", "beta"),
      motifs = c("left-hit-wide", "right-hit-wide")
    )),
    matrix(c(0, 1, 2, 1), nrow = 2, byrow = TRUE)
  )

  beta_dense <- end_motif_data_frame(motif_grouped, groups = "beta", densify = TRUE)
  expect_equal(beta_dense$motif, c("left-hit-wide", "right-hit-wide"))
  expect_equal(beta_dense$count, c(2, 1))
  expect_error(
    end_motif_data_frame(motif_grouped, motifs = "GT_AC"),
    "Unknown end-motif label",
    fixed = TRUE
  )
})

test_that("R helper package reads dense global reference k-mers", {
  ref_kmers <- read_ref_kmers(dense_global_ref_kmer_zarr_path())

  expect_s3_class(ref_kmers, "cfdnalab_global_ref_kmer_frequencies")
  expect_identical(schema_version(ref_kmers), 1L)
  expect_identical(storage_mode(ref_kmers), "dense")
  expect_identical(row_mode(ref_kmers), "global")
  expect_identical(motif_axis_kind(ref_kmers), "motif")
  expect_identical(kmer_size(ref_kmers), 3L)
  expect_true(canonical(ref_kmers))
  expect_true(all_motifs(ref_kmers))
  expect_identical(assign_by(ref_kmers), "count-overlap")
  expect_identical(length(motifs(ref_kmers)$motif), 32L)
  expect_identical(head(motifs(ref_kmers)$motif, 4L), c("AAA", "AAC", "AAG", "AAT"))
  expect_identical(tail(motifs(ref_kmers)$motif, 4L), c("TCA", "TCC", "TCG", "TCT"))
  expect_equal(
    dense_frequencies_matrix(ref_kmers, motifs = c("AAA", "ACG", "TAC", "CAT")),
    matrix(c(4 / 36, 7 / 36, 6 / 36, 0), nrow = 1L),
    tolerance = 1e-8
  )
  expect_equal(
    dense_counts_matrix(ref_kmers, motifs = c("AAA", "ACG", "TAC", "CAT")),
    matrix(c(4, 7, 6, 0), nrow = 1L),
    tolerance = 1e-8
  )
  expect_equal(row_scaling_factors(ref_kmers)$row_scaling_factor, 36, tolerance = 1e-8)
  expect_equal(
    ref_kmer_data_frame(ref_kmers, motifs = c("TAC", "AAA")),
    data.frame(
      row_label = c("global", "global"),
      motif_idx = c(26L, 1L),
      motif = c("TAC", "AAA"),
      frequency = c(6 / 36, 4 / 36),
      count = c(6, 4),
      stringsAsFactors = FALSE
    ),
    tolerance = 1e-8,
    ignore_attr = TRUE
  )
})

test_that("R helper package reads sparse windowed reference k-mers", {
  ref_kmers <- read_ref_kmers(sparse_windowed_ref_kmer_zarr_path())

  expect_s3_class(ref_kmers, "cfdnalab_windowed_ref_kmer_frequencies")
  expect_identical(storage_mode(ref_kmers), "sparse_coo")
  expect_identical(row_mode(ref_kmers), "bed")
  expect_identical(motifs(ref_kmers)$motif, c("CGT", "AAA", "TAC", "CCC", "GGG", "ACG", "GTA"))
  expect_false(has_motif(ref_kmers, "TTT"))
  expect_error(dense_frequencies_matrix(ref_kmers), "Use sparse_frequencies_matrix")
  expect_equal(
    dense_counts_matrix(ref_kmers, allow_densify = TRUE),
    matrix(
      c(
        0, 1, 0, 1, 1, 0, 0,
        1, 0, 4 / 3, 0, 1 / 3, 1, 2 / 3,
        2 / 3, 0, 0, 1, 1 / 3, 0, 1 / 3,
        5 / 3, 0, 1, 0, 0, 1, 2
      ),
      nrow = 4L,
      byrow = TRUE
    ),
    tolerance = 1e-8
  )
  expect_equal(
    row_scaling_factors(ref_kmers)$row_scaling_factor,
    c(3, 13 / 3, 7 / 3, 17 / 3),
    tolerance = 1e-8
  )
  expect_equal(
    window_metadata(ref_kmers),
    data.frame(
      window_idx = 1:4,
      chrom = c("chr1", "chr1", "chr2", "chr2"),
      start = c(0L, 8L, 2L, 12L),
      end = c(9L, 16L, 13L, 20L),
      blacklisted_fraction = c(0, 1 / 8, 2 / 11, 0),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_equal(
    ref_kmer_data_frame(
      ref_kmers,
      window_idxs = c(4L, 1L),
      motifs = c("GTA", "AAA"),
      densify = TRUE
    ),
    data.frame(
      window_idx = c(4L, 4L, 1L, 1L),
      chrom = c("chr2", "chr2", "chr1", "chr1"),
      start = c(12L, 12L, 0L, 0L),
      end = c(20L, 20L, 9L, 9L),
      blacklisted_fraction = c(0, 0, 0, 0),
      motif_idx = c(7L, 2L, 7L, 2L),
      motif = c("GTA", "AAA", "GTA", "AAA"),
      frequency = c(6 / 17, 0, 0, 1 / 3),
      count = c(2, 0, 0, 1),
      stringsAsFactors = FALSE
    ),
    tolerance = 1e-8,
    ignore_attr = TRUE
  )
})

test_that("R helper package reads sparse grouped reference k-mers", {
  ref_kmers <- read_ref_kmers(sparse_grouped_ref_kmer_zarr_path())

  expect_s3_class(ref_kmers, "cfdnalab_grouped_ref_kmer_frequencies")
  expect_identical(storage_mode(ref_kmers), "sparse_coo")
  expect_identical(row_mode(ref_kmers), "grouped_bed")
  expect_identical(group_idx(ref_kmers, "alpha"), 2L)
  expect_identical(motifs(ref_kmers)$motif, c("CGT", "AAA", "TAC", "CCC", "GGG", "ACG", "GTA"))
  expect_equal(
    group_metadata(ref_kmers),
    data.frame(
      group_idx = c(1L, 2L),
      group_name = c("beta", "alpha"),
      eligible_windows = c(2L, 2L),
      blacklisted_fraction = c(0.1, 1 / 16),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_equal(
    as.matrix(sparse_counts_matrix(
      ref_kmers,
      groups = c("alpha", "beta"),
      motifs = c("GTA", "AAA")
    )),
    matrix(c(8 / 3, 0, 1 / 3, 1), nrow = 2L, byrow = TRUE),
    tolerance = 1e-8
  )
  expect_equal(
    ref_kmer_data_frame(ref_kmers),
    data.frame(
      group_idx = c(rep(1L, 5L), rep(2L, 5L)),
      group_name = c(rep("beta", 5L), rep("alpha", 5L)),
      eligible_windows = rep(2L, 10L),
      blacklisted_fraction = c(rep(0.1, 5L), rep(1 / 16, 5L)),
      motif_idx = c(1L, 2L, 4L, 5L, 7L, 1L, 3L, 5L, 6L, 7L),
      motif = c("CGT", "AAA", "CCC", "GGG", "GTA", "CGT", "TAC", "GGG", "ACG", "GTA"),
      frequency = c(1 / 8, 3 / 16, 3 / 8, 1 / 4, 1 / 16, 4 / 15, 7 / 30, 1 / 30, 1 / 5, 4 / 15),
      count = c(2 / 3, 1, 2, 4 / 3, 1 / 3, 8 / 3, 7 / 3, 1 / 3, 2, 8 / 3),
      stringsAsFactors = FALSE
    ),
    tolerance = 1e-8,
    ignore_attr = TRUE
  )
})

test_that("R helper package reads dense grouped motif-group reference k-mers", {
  ref_kmers <- read_ref_kmers(dense_grouped_motif_group_ref_kmer_zarr_path())

  expect_s3_class(ref_kmers, "cfdnalab_grouped_ref_kmer_frequencies")
  expect_identical(storage_mode(ref_kmers), "dense")
  expect_identical(row_mode(ref_kmers), "grouped_bed")
  expect_identical(motif_axis_kind(ref_kmers), "motif_group")
  expect_true(all_motifs(ref_kmers))
  expect_identical(motifs(ref_kmers)$motif, c("absent", "edge", "gc_rich", "homopolymer", "transition"))
  expect_equal(
    dense_counts_matrix(ref_kmers),
    matrix(c(0, 1 / 3, 10 / 3, 1, 2 / 3, 0, 5, 1 / 3, 0, 14 / 3), nrow = 2L, byrow = TRUE),
    tolerance = 1e-8
  )
  expect_equal(
    dense_counts_matrix(ref_kmers, groups = c("alpha", "beta"), motifs = c("edge", "absent")),
    matrix(c(5, 0, 1 / 3, 0), nrow = 2L, byrow = TRUE),
    tolerance = 1e-8
  )
  expect_equal(
    ref_kmer_data_frame(ref_kmers, groups = c("alpha", "beta"), motifs = c("edge", "absent")),
    data.frame(
      group_idx = c(2L, 2L, 1L, 1L),
      group_name = c("alpha", "alpha", "beta", "beta"),
      eligible_windows = c(2L, 2L, 2L, 2L),
      blacklisted_fraction = c(1 / 16, 1 / 16, 0.1, 0.1),
      motif_idx = c(2L, 1L, 2L, 1L),
      motif = c("edge", "absent", "edge", "absent"),
      frequency = c(1 / 2, 0, 1 / 16, 0),
      count = c(5, 0, 1 / 3, 0),
      stringsAsFactors = FALSE
    ),
    tolerance = 1e-8,
    ignore_attr = TRUE
  )
})

test_that("R helper package reads global length counts", {
  lengths <- read_lengths(global_length_counts_path())

  expect_s3_class(lengths, "cfdnalab_global_length_counts")
  expect_equal(
    length_bins(lengths),
    data.frame(
      length_bin_idx = c(1L, 2L, 3L),
      length_start_bp = c(30L, 50L, 70L),
      length_end_bp = c(50L, 70L, 100L),
      length_midpoint_bp = c(40, 60, 85),
      length_width_bp = c(20L, 20L, 30L),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_equal(length_counts_vector(lengths), c(count_30_50 = 3, count_50_70 = 2, count_70_100 = 1))

  fractions <- length_data_frame(lengths, value = "fraction")
  selected_fractions <- length_data_frame(
    lengths,
    with_length_range = c(50L, 100L),
    value = "fraction",
    denominator = "selected_bins"
  )
  expect_equal(fractions$fraction, c(0.5, 1 / 3, 1 / 6), tolerance = 1e-8)
  expect_equal(selected_fractions$length_bin_idx, c(2L, 3L))
  expect_equal(selected_fractions$fraction, c(2 / 3, 1 / 3), tolerance = 1e-8)
  expect_equal(length_data_frame(lengths, value = "density")$density, c(0.025, 1 / 60, 1 / 180), tolerance = 1e-8)
})

test_that("R helper package reads windowed length counts", {
  lengths <- read_lengths(windowed_length_counts_path())

  expect_s3_class(lengths, "cfdnalab_windowed_length_counts")
  expect_equal(
    window_metadata(lengths),
    data.frame(
      window_idx = 1:4,
      chrom = rep("chr1", 4),
      start = c(0L, 100L, 200L, 300L),
      end = c(100L, 200L, 300L, 360L),
      blacklisted_fraction = c(0.04, 0.05, 0.1, 0.25),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )

  expect_equal(
    length_counts_matrix(lengths),
    matrix(
      c(2, 0, 0, 0, 2, 0, 0, 0, 1, 1, 0, 0),
      nrow = 4,
      byrow = TRUE,
      dimnames = list(NULL, c("count_30_50", "count_50_70", "count_70_100"))
    ),
    tolerance = 1e-8
  )

  selected <- length_data_frame(lengths, window_idxs = c(2L, 4L), value = "fraction", keep_wide = TRUE)
  expect_equal(names(selected), c(
    "window_idx",
    "chrom",
    "start",
    "end",
    "blacklisted_fraction",
    "fraction_30_50",
    "fraction_50_70",
    "fraction_70_100"
  ))
  expect_equal(selected$fraction_30_50, c(0, 1), tolerance = 1e-8)
  expect_equal(selected$fraction_50_70, c(1, 0), tolerance = 1e-8)
  expect_equal(selected$fraction_70_100, c(0, 0), tolerance = 1e-8)

  filtered <- length_data_frame(lengths, max_blacklisted_fraction = 0.05)
  expect_identical(unique(filtered$window_idx), c(1L, 2L))

  range_fraction <- length_data_frame(
    lengths,
    window_idxs = c(2L, 4L),
    with_length_range = c(50L, 100L),
    value = "fraction",
    denominator = "selected_bins",
    keep_wide = TRUE
  )
  expect_equal(range_fraction$window_idx, c(2L, 4L))
  expect_equal(range_fraction$fraction_50_70, c(1, NA_real_), tolerance = 1e-8)
  expect_equal(range_fraction$fraction_70_100, c(0, NA_real_), tolerance = 1e-8)
})

test_that("R helper package reads grouped length counts", {
  lengths <- read_lengths(grouped_length_counts_path())

  expect_s3_class(lengths, "cfdnalab_grouped_length_counts")
  expect_equal(
    group_metadata(lengths),
    data.frame(
      group_idx = 1:4,
      group_name = c("beta", "alpha", "gamma", "zero"),
      eligible_windows = c(2L, 1L, 1L, 1L),
      blacklisted_fraction = c(0.07, 0.05, 0.25, 0.333),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_identical(group_idx(lengths, "gamma"), 3L)

  beta <- length_data_frame(lengths, groups = "beta")
  expect_equal(beta$count, c(2, 0, 1), tolerance = 1e-8)

  wide_density <- length_data_frame(lengths, groups = c("alpha", "zero"), value = "density", keep_wide = TRUE)
  expect_equal(names(wide_density), c(
    "group_idx",
    "group_name",
    "eligible_windows",
    "blacklisted_fraction",
    "density_30_50",
    "density_50_70",
    "density_70_100"
  ))
  expect_equal(wide_density$density_30_50, c(0, NA_real_), tolerance = 1e-8)
  expect_equal(wide_density$density_50_70, c(1 / 20, NA_real_), tolerance = 1e-8)
  expect_equal(wide_density$density_70_100, c(0, NA_real_), tolerance = 1e-8)

  selected_range <- length_data_frame(
    lengths,
    groups = c("beta", "zero"),
    with_length_range = c(50L, 100L),
    value = "fraction",
    denominator = "selected_bins"
  )
  expect_equal(selected_range$group_name, c("beta", "beta", "zero", "zero"))
  expect_equal(selected_range$length_bin_idx, c(2L, 3L, 2L, 3L))
  expect_equal(selected_range$fraction, c(0, 1, NA_real_, NA_real_), tolerance = 1e-8)
})

test_that("R helper package reads no-blacklist windowed length counts", {
  lengths <- read_lengths(windowed_length_counts_no_blacklist_path())

  expect_s3_class(lengths, "cfdnalab_windowed_length_counts")
  expect_equal(
    window_metadata(lengths),
    data.frame(
      window_idx = 1:4,
      chrom = rep("chr1", 4),
      start = c(0L, 100L, 200L, 300L),
      end = c(100L, 200L, 300L, 360L),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_equal(
    length_counts_matrix(lengths),
    matrix(
      c(2, 0, 0, 0, 2, 0, 0, 0, 1, 1, 0, 0),
      nrow = 4,
      byrow = TRUE,
      dimnames = list(NULL, c("count_30_50", "count_50_70", "count_70_100"))
    ),
    tolerance = 1e-8
  )
  expect_equal(
    length_data_frame(lengths, max_blacklisted_fraction = 1)$count,
    c(2, 0, 0, 0, 2, 0, 0, 0, 1, 1, 0, 0),
    tolerance = 1e-8
  )
  expect_error(
    length_data_frame(lengths, max_blacklisted_fraction = 0.5),
    "has no blacklisted_fraction column"
  )
})

test_that("R helper package reads no-blacklist grouped length counts", {
  lengths <- read_lengths(grouped_length_counts_no_blacklist_path())

  expect_s3_class(lengths, "cfdnalab_grouped_length_counts")
  expect_equal(
    group_metadata(lengths),
    data.frame(
      group_idx = 1:4,
      group_name = c("beta", "alpha", "gamma", "zero"),
      eligible_windows = c(2L, 1L, 1L, 1L),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_equal(length_data_frame(lengths, groups = "beta")$count, c(2, 0, 1), tolerance = 1e-8)
  expect_equal(
    length_data_frame(lengths, group_idxs = 4L, value = "fraction")$fraction,
    c(NA_real_, NA_real_, NA_real_)
  )
  expect_error(
    length_data_frame(lengths, max_blacklisted_fraction = 0.5),
    "has no blacklisted_fraction column"
  )
})
