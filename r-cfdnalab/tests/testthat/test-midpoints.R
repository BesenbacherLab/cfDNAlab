test_that("midpoint loader reads locally generated schema fixture", {
  midpoints <- read_midpoints(make_midpoint_zarr_fixture())

  expect_s3_class(midpoints, "cfdnalab_midpoint_profiles")
  expect_equal(schema_version(midpoints), 1L)
  expect_equal(group_idx(midpoints, "long_group"), 2L)
  expect_equal(length_bin_idx(midpoints, 60), 2L)
  expect_output(print(midpoints), "<cfDNAlab midpoint profiles>", fixed = TRUE)
  expect_equal(
    groups(midpoints),
    data.frame(
      group_idx = c(1L, 2L),
      group_name = c("A", "long_group"),
      eligible_intervals = c(1L, 3L),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_equal(
    length_bins(midpoints),
    data.frame(
      length_bin_idx = c(1L, 2L),
      length_start_bp = c(30L, 60L),
      length_end_bp = c(60L, 90L),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_equal(
    positions(midpoints)$position_bin_start_bp,
    c(0L, 5L, 10L, 15L)
  )

  profile <- profile_data_frame(midpoints, group = "long_group", length = 60)
  expect_equal(profile$count, c(5, 0, 0, 6.5), tolerance = 1e-8)
  profile_by_index <- profile_data_frame(midpoints, group_idx = 1L, length_bin_idx = 1L)
  expect_equal(profile_by_index$count, c(0, 1.5, 0, 2.25), tolerance = 1e-8)
  expect_false(any(grepl("idx0|index0", names(profile_by_index))))
  expect_equal(
    profile_array(midpoints, group_idx = 1L, length_bin_idx = 2L),
    c(3, 0, 4.5, 0),
    tolerance = 1e-8
  )
  expect_equal(dim(midpoint_array(midpoints)), c(2L, 2L, 4L))
  expect_error(
    profile_data_frame(midpoints, group_idx = 0L, length_bin_idx = 1L),
    "group_idx 0 is outside 1..2",
    fixed = TRUE
  )
  expect_error(
    profile_array(midpoints, group_idx = 1L, length_bin_idx = 0L),
    "length_bin_idx 0 is outside 1..2",
    fixed = TRUE
  )
  expect_error(
    profile_array(midpoints, group_idx = 1L, group = "A", length_bin_idx = 1L),
    "Use either group_idx or group, not both",
    fixed = TRUE
  )
  expect_error(
    profile_array(midpoints, group_idx = 1L, length_bin_idx = 1L, group_index = 0L),
    "Unused argument(s): group_index",
    fixed = TRUE
  )
})

test_that("midpoint loader reads metadata and one profile", {
  testthat::skip_if_not_installed("zarr")

  midpoints <- read_midpoints(midpoint_fixture_path())

  expect_s3_class(midpoints, "cfdnalab_midpoint_profiles")
  expect_equal(schema_version(midpoints), 1L)

  expect_equal(
    groups(midpoints),
    data.frame(
      group_idx = c(1L, 2L, 3L),
      group_name = c("alpha", "beta-site", "gamma_long"),
      eligible_intervals = c(2, 2, 2),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_equal(
    length_bins(midpoints),
    data.frame(
      length_bin_idx = c(1L, 2L, 3L),
      length_start_bp = c(30, 50, 70),
      length_end_bp = c(50, 70, 100),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_equal(
    positions(midpoints),
    data.frame(
      position_idx = c(1L, 2L, 3L, 4L, 5L),
      position_bin_start_bp = c(0L, 2L, 4L, 6L, 8L),
      position_bin_end_bp = c(2L, 4L, 6L, 8L, 10L),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )

  expect_equal(group_idx(midpoints, "beta-site"), 2L)
  expect_equal(length_bin_idx(midpoints, 50), 2L)
  expect_error(
    length_bin_idx(midpoints, 50.5),
    "Fragment length must be a single non-negative integer",
    fixed = TRUE
  )

  profile <- profile_data_frame(midpoints, group = "beta-site", length_bin_idx = 2L)
  expect_equal(profile$group_idx, rep(2L, 5L))
  expect_equal(profile$group_name, rep("beta-site", 5L))
  expect_equal(profile$length_bin_idx, rep(2L, 5L))
  expect_equal(profile$position_idx, c(1L, 2L, 3L, 4L, 5L))
  expect_equal(profile$count, c(0, 0, 1.5, 0, 0.5), tolerance = 1e-8)
})

test_that("midpoint arrays preserve declared shape", {
  testthat::skip_if_not_installed("zarr")

  midpoints <- read_midpoints(midpoint_fixture_path())
  counts <- midpoint_array(midpoints)

  expect_equal(dim(counts), c(3L, 3L, 5L))
  expect_equal(counts[1, 1, 1], 1)
  expect_equal(counts[1, 2, 3], 1.5)
  expect_equal(counts[3, 1, 3], 2.5)
  expect_equal(
    profile_array(midpoints, group_idx = 3L, length_bin_idx = 1L),
    c(0, 0, 2.5, 0, 0),
    tolerance = 1e-8
  )
})
