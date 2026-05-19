write_length_tsv_fixture <- function(lines) {
  path <- tempfile(fileext = ".tsv")
  writeLines(lines, path, useBytes = TRUE)
  path
}

test_that("global length counts expose bins, matrix, vector, and long values", {
  path <- write_length_tsv_fixture(c(
    "count_30\tcount_31_40",
    "12\t3"
  ))

  lengths <- read_lengths(path)

  expect_s3_class(lengths, "cfdnalab_global_length_counts")
  expect_output(print(lengths), "<cfDNAlab length counts>", fixed = TRUE)
  expect_equal(
    length_bins(lengths),
    data.frame(
      length_bin_idx = c(1L, 2L),
      length_start_bp = c(30L, 31L),
      length_end_bp = c(31L, 40L),
      length_midpoint_bp = c(30.5, 35.5),
      length_width_bp = c(1L, 9L),
      count_column = c("count_30", "count_31_40"),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_equal(length_bin_idx(lengths, 30L), 1L)
  expect_equal(length_bin_idx(lengths, 39L), 2L)
  expect_error(length_bin_idx(lengths, 40L), "No length-count bin contains length 40")
  expect_equal(
    length_counts_matrix(lengths),
    matrix(c(12, 3), nrow = 1L, dimnames = list(NULL, c("count_30", "count_31_40"))),
    tolerance = 1e-8
  )
  expect_equal(
    length_counts_vector(lengths),
    c(count_30 = 12, count_31_40 = 3),
    tolerance = 1e-8
  )

  fractions <- length_data_frame(lengths, value = "fraction")
  expect_equal(fractions$fraction, c(0.8, 0.2), tolerance = 1e-8)
  densities <- length_data_frame(lengths, value = "density")
  expect_equal(densities$density, c(0.8, 0.2 / 9), tolerance = 1e-8)
  expect_false(any(grepl("idx0|index0", names(fractions))))
})

test_that("global length counts support wide count, fraction, and density frames", {
  path <- write_length_tsv_fixture(c(
    "count_30_50\tcount_50_70\tcount_70_100",
    "3\t2\t1"
  ))

  lengths <- read_lengths(path)

  expect_equal(
    length_data_frame(lengths, keep_wide = TRUE),
    data.frame(
      count_30_50 = 3,
      count_50_70 = 2,
      count_70_100 = 1,
      check.names = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_equal(
    length_data_frame(lengths, value = "fraction", keep_wide = TRUE),
    data.frame(
      fraction_30_50 = 0.5,
      fraction_50_70 = 1 / 3,
      fraction_70_100 = 1 / 6,
      check.names = FALSE
    ),
    tolerance = 1e-8,
    ignore_attr = TRUE
  )
  expect_equal(
    length_data_frame(lengths, value = "density", keep_wide = TRUE),
    data.frame(
      density_30_50 = 0.025,
      density_50_70 = 1 / 60,
      density_70_100 = 1 / 180,
      check.names = FALSE
    ),
    tolerance = 1e-8,
    ignore_attr = TRUE
  )
  expect_error(length_data_frame(lengths, keep_wide = NA), "keep_wide must be TRUE or FALSE")
})

test_that("windowed length counts expose window metadata and selected frames", {
  path <- write_length_tsv_fixture(c(
    "chrom\tstart\tend\tblacklisted_fraction\tcount_30\tcount_31_40",
    "chr1\t10\t20\t0.25\t12\t3",
    "chr2\t30\t45\t0\t0\t5"
  ))

  lengths <- read_lengths(path)

  expect_s3_class(lengths, "cfdnalab_windowed_length_counts")
  expect_equal(
    window_metadata(lengths),
    data.frame(
      window_idx = c(1L, 2L),
      chrom = c("chr1", "chr2"),
      start = c(10L, 30L),
      end = c(20L, 45L),
      blacklisted_fraction = c(0.25, 0),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )

  selected <- length_data_frame(lengths, window_idx = 2L)
  expect_equal(selected$window_idx, c(2L, 2L))
  expect_equal(selected$count, c(0, 5), tolerance = 1e-8)

  filtered <- length_data_frame(lengths, max_blacklisted_fraction = 0.1)
  expect_equal(unique(filtered$window_idx), 2L)
  expect_error(
    length_data_frame(lengths, max_blacklisted_fraction = 1.1),
    "max_blacklisted_fraction must be a single finite fraction in 0..1"
  )

  wide_fraction <- length_data_frame(lengths, value = "fraction", keep_wide = TRUE)
  expect_equal(names(wide_fraction), c(
    "window_idx",
    "chrom",
    "start",
    "end",
    "blacklisted_fraction",
    "fraction_30",
    "fraction_31_40"
  ))
  expect_equal(wide_fraction$fraction_30, c(0.8, 0), tolerance = 1e-8)
  expect_equal(wide_fraction$fraction_31_40, c(0.2, 1), tolerance = 1e-8)
})

test_that("windowed length selectors validate one-based indices", {
  path <- write_length_tsv_fixture(c(
    "chrom\tstart\tend\tcount_30\tcount_31_40",
    "chr1\t10\t20\t1\t2",
    "chr1\t30\t40\t3\t4"
  ))

  lengths <- read_lengths(path)

  expect_error(length_data_frame(lengths, window_idx = 0L), "window_idx contains values outside 1..2")
  expect_error(length_data_frame(lengths, window_idx = 3L), "window_idx contains values outside 1..2")
  expect_error(length_data_frame(lengths, window_idx = 1.5), "window_idx must contain integer values")
  expect_equal(length_data_frame(lengths, window_idx = integer(0))$count, numeric(0))
})

test_that("grouped length counts expose group metadata and group selectors", {
  path <- write_length_tsv_fixture(c(
    "group_name\teligible_windows\tcount_30\tcount_31_40",
    "alpha\t2\t12\t3",
    "beta\t0\t0\t0"
  ))

  lengths <- read_lengths(path)

  expect_s3_class(lengths, "cfdnalab_grouped_length_counts")
  expect_equal(
    group_metadata(lengths),
    data.frame(
      group_idx = c(1L, 2L),
      group_name = c("alpha", "beta"),
      eligible_windows = c(2L, 0L),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_equal(group_idx(lengths, "beta"), 2L)
  expect_error(group_idx(lengths, "gamma"), "Unknown length-count group name")

  alpha <- length_data_frame(lengths, group = "alpha")
  expect_equal(alpha$group_idx, c(1L, 1L))
  expect_equal(alpha$count, c(12, 3), tolerance = 1e-8)

  beta_wide <- length_data_frame(lengths, group_idx = 2L, value = "density", keep_wide = TRUE)
  expect_equal(names(beta_wide), c("group_idx", "group_name", "eligible_windows", "density_30", "density_31_40"))
  expect_true(all(is.na(beta_wide[c("density_30", "density_31_40")])))
})

test_that("grouped length selectors validate names and one-based indices", {
  path <- write_length_tsv_fixture(c(
    "group_name\teligible_windows\tblacklisted_fraction\tcount_30\tcount_31_40",
    "alpha\t2\t0.25\t12\t3",
    "beta\t1\t0\t0\t5",
    "gamma\t0\t1\t0\t0"
  ))

  lengths <- read_lengths(path)

  expect_error(length_data_frame(lengths, group = "alpha", group_idx = 1L), "Use either group or group_idx")
  expect_error(length_data_frame(lengths, group = "missing"), "Unknown length-count group name")
  expect_error(length_data_frame(lengths, group_idx = 4L), "group_idx contains values outside 1..3")
  expect_equal(unique(length_data_frame(lengths, max_blacklisted_fraction = 0.25)$group_name), c("alpha", "beta"))
  expect_equal(length_data_frame(lengths, group = c("beta", "alpha"))$group_idx, c(2L, 2L, 1L, 1L))
})

test_that("length TSV validation rejects ambiguous or unsupported shapes", {
  missing_zstd <- write_length_tsv_fixture(c(
    "chrom\tstart\tend\tcount_30",
    "chr1\t10\t20\t1"
  ))
  zst_path <- sub("\\.tsv$", ".tsv.zst", missing_zstd)
  file.copy(missing_zstd, zst_path)
  if (!nzchar(Sys.which("zstd"))) {
    expect_error(read_lengths(zst_path), "requires the zstd command-line tool")
  }

  bad_count <- write_length_tsv_fixture(c(
    "count_30\tother",
    "1\t2"
  ))
  expect_error(read_lengths(bad_count), "count columns must be contiguous")

  no_blacklist <- write_length_tsv_fixture(c(
    "chrom\tstart\tend\tcount_30",
    "chr1\t10\t20\t1"
  ))
  lengths <- read_lengths(no_blacklist)
  expect_equal(length_data_frame(lengths)$count, 1)
  expect_equal(length_data_frame(lengths, max_blacklisted_fraction = 1)$count, 1)
  expect_error(
    length_data_frame(lengths, max_blacklisted_fraction = 0.1),
    "has no blacklisted_fraction column"
  )

  multiple_global_rows <- write_length_tsv_fixture(c(
    "count_30",
    "1",
    "2"
  ))
  expect_error(read_lengths(multiple_global_rows), "Global length-count output must contain exactly one row")

  duplicate_bins <- write_length_tsv_fixture(c(
    "count_30\tcount_30_31",
    "1\t2"
  ))
  expect_error(read_lengths(duplicate_bins), "duplicate length bins")

  unsupported_metadata <- write_length_tsv_fixture(c(
    "window_idx\tchrom\tstart\tend\tcount_30",
    "1\tchr1\t10\t20\t1"
  ))
  expect_error(read_lengths(unsupported_metadata), "Could not infer length-count output mode")
})

test_that("length-bin lookup rejects gaps and overlapping bins", {
  gapped <- read_lengths(write_length_tsv_fixture(c(
    "count_30_40\tcount_50_60",
    "1\t2"
  )))
  expect_error(length_bin_idx(gapped, 45L), "No length-count bin contains length 45")

  overlapping <- read_lengths(write_length_tsv_fixture(c(
    "count_30_50\tcount_40_60",
    "1\t2"
  )))
  expect_error(length_bin_idx(overlapping, 45L), "Multiple length-count bins contain length 45")
})
