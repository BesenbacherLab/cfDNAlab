source(file.path("downstream_tests", "R", "common.R"))

library(cfdnalab)
library(testthat)

test_that("R helper package reads midpoint profiles", {
  midpoints <- read_midpoints(midpoint_zarr_path())

  expect_identical(schema_version(midpoints), 1L)
  expect_identical(
    groups(midpoints)$group_name,
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

  beta_profile <- profile_data_frame(midpoints, group = "beta-site", length_bin_idx = 2L)
  expect_identical(beta_profile$position_idx, c(1L, 2L, 3L, 4L, 5L))
  expect_equal(beta_profile$count, c(0, 0, 1.5, 0, 0.5))
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
    dense_data_frame(dense_global)$count,
    c(1, 0, 1, 0)
  )
})

test_that("R helper package reads sparse windowed end motifs", {
  sparse_windowed <- read_end_motifs(sparse_windowed_end_zarr_path())

  expect_identical(storage_mode(sparse_windowed), "sparse_coo")
  expect_identical(row_mode(sparse_windowed), "bed")
  expect_identical(windows(sparse_windowed)$window_idx, c(1L, 2L))
  expect_equal(
    as.matrix(sparse_counts_matrix(sparse_windowed)),
    matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
  )
  expect_equal(sparse_data_frame_for_window(sparse_windowed, 1L)$count, 1)
  expect_equal(sparse_data_frame_for_motif(sparse_windowed, "_G")$count, 1)
})

test_that("R helper package reads sparse grouped end motifs", {
  sparse_grouped <- read_end_motifs(sparse_grouped_end_zarr_path())

  expect_identical(storage_mode(sparse_grouped), "sparse_coo")
  expect_identical(row_mode(sparse_grouped), "grouped_bed")
  expect_identical(group_idx(sparse_grouped, "alpha"), 2L)
  expect_identical(
    groups(sparse_grouped)$group_name,
    c("beta", "alpha", "gamma")
  )
  expect_equal(
    as.matrix(sparse_counts_matrix(sparse_grouped)),
    matrix(c(1, 2, 1, 0, 0, 0), nrow = 3, byrow = TRUE)
  )
  expect_equal(sparse_data_frame_for_group(sparse_grouped, "beta")$count, c(1, 2))
})
