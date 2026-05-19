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

  expect_identical(schema_version(dense_global), 1L)
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
})

test_that("R helper package reads sparse windowed end motifs", {
  sparse_windowed <- read_end_motifs(sparse_windowed_end_zarr_path())

  expect_identical(storage_mode(sparse_windowed), "sparse_coo")
  expect_identical(row_mode(sparse_windowed), "bed")
  expect_identical(window_metadata(sparse_windowed)$window_idx, c(1L, 2L))
  expect_identical(window_metadata(sparse_windowed)$chrom, c("chr1", "chr1"))
  expect_identical(window_metadata(sparse_windowed)$start, c(10L, 19L))
  expect_identical(window_metadata(sparse_windowed)$end, c(11L, 20L))
  expect_equal(
    as.matrix(sparse_counts_matrix(sparse_windowed)),
    matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
  )
  expect_equal(
    end_motif_data_frame(
      sparse_windowed,
      motifs = "_A",
      densify = TRUE,
      max_blacklisted_fraction = 0
    )$count,
    c(0, 1)
  )
  expect_equal(end_motif_data_frame(sparse_windowed, window_idxs = 1L)$count, 1)
  expect_equal(end_motif_data_frame(sparse_windowed, motifs = "_G")$count, 1)

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
  expect_equal(range_fraction$fraction_50_70, c(1, NA), tolerance = 1e-8)
  expect_equal(range_fraction$fraction_70_100, c(0, NA), tolerance = 1e-8)
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
  expect_equal(wide_density$density_30_50, c(0, NA), tolerance = 1e-8)
  expect_equal(wide_density$density_50_70, c(1 / 20, NA), tolerance = 1e-8)
  expect_equal(wide_density$density_70_100, c(0, NA), tolerance = 1e-8)

  selected_range <- length_data_frame(
    lengths,
    groups = c("beta", "zero"),
    with_length_range = c(50L, 100L),
    value = "fraction",
    denominator = "selected_bins"
  )
  expect_equal(selected_range$group_name, c("beta", "beta", "zero", "zero"))
  expect_equal(selected_range$length_bin_idx, c(2L, 3L, 2L, 3L))
  expect_equal(selected_range$fraction, c(0, 1, NA, NA), tolerance = 1e-8)
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
  expect_equal(length_data_frame(lengths, group_idxs = 4L, value = "fraction")$fraction, c(NA, NA, NA))
  expect_error(
    length_data_frame(lengths, max_blacklisted_fraction = 0.5),
    "has no blacklisted_fraction column"
  )
})
