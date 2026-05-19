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
    length_counts_matrix(lengths, with_length_range = c(31L, 35L)),
    matrix(3, nrow = 1L, dimnames = list(NULL, "count_31_40")),
    tolerance = 1e-8
  )
  expect_equal(
    length_counts_vector(lengths),
    c(count_30 = 12, count_31_40 = 3),
    tolerance = 1e-8
  )
  expect_equal(
    length_counts_vector(lengths, with_lengths = 39L),
    c(count_31_40 = 3),
    tolerance = 1e-8
  )

  fractions <- length_data_frame(lengths, value = "fraction")
  expect_equal(fractions$fraction, c(0.8, 0.2), tolerance = 1e-8)
  fraction_range <- length_data_frame(lengths, with_length_range = c(31L, 35L), value = "fraction")
  expect_equal(fraction_range$length_bin_idx, 2L)
  expect_equal(fraction_range$fraction, 0.2, tolerance = 1e-8)
  fraction_range_selected <- length_data_frame(
    lengths,
    with_length_range = c(31L, 35L),
    value = "fraction",
    denominator = "selected_bins"
  )
  expect_equal(fraction_range_selected$length_bin_idx, 2L)
  expect_equal(fraction_range_selected$fraction, 1, tolerance = 1e-8)
  densities <- length_data_frame(lengths, value = "density")
  expect_equal(densities$density, c(0.8, 0.2 / 9), tolerance = 1e-8)
  density_range_selected <- length_data_frame(
    lengths,
    with_length_range = c(31L, 35L),
    value = "density",
    denominator = "selected_bins",
    keep_wide = TRUE
  )
  expect_equal(density_range_selected$density_31_40, 1 / 9, tolerance = 1e-8)
  expect_false(any(grepl("idx0|index0", names(fractions))))
  expect_error(
    length_data_frame(lengths, max_blacklisted_fraction = 1),
    "Unused argument\\(s\\): max_blacklisted_fraction"
  )
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

test_that("wide length frames preserve TSV count column labels internally", {
  path <- write_length_tsv_fixture(c(
    "count_030\tcount_031_040",
    "12\t3"
  ))

  lengths <- read_lengths(path)

  expect_false("count_column" %in% names(length_bins(lengths)))
  expect_equal(
    names(length_data_frame(lengths, keep_wide = TRUE)),
    c("count_030", "count_031_040")
  )
  expect_equal(
    names(length_data_frame(lengths, value = "fraction", keep_wide = TRUE)),
    c("fraction_030", "fraction_031_040")
  )
})

test_that("length-bin selectors preserve order, boundaries, and row totals", {
  path <- write_length_tsv_fixture(c(
    "count_30_50\tcount_50_70\tcount_70_100",
    "3\t2\t1"
  ))

  lengths <- read_lengths(path)

  expect_equal(
    length_counts_matrix(lengths, with_lengths = c(75L, 30L)),
    matrix(c(1, 3), nrow = 1L, dimnames = list(NULL, c("count_70_100", "count_30_50"))),
    tolerance = 1e-8
  )
  expect_equal(
    length_counts_matrix(lengths, length_bin_idxs = c(3L, 1L)),
    matrix(c(1, 3), nrow = 1L, dimnames = list(NULL, c("count_70_100", "count_30_50"))),
    tolerance = 1e-8
  )

  empty_counts <- length_counts_matrix(lengths, with_lengths = integer(0))
  expect_equal(dim(empty_counts), c(1L, 0L))
  expect_null(colnames(empty_counts))

  exact_range <- length_data_frame(lengths, with_length_range = c(50L, 70L), value = "fraction")
  expect_equal(exact_range$length_bin_idx, 2L)
  expect_equal(exact_range$fraction, 2 / 6, tolerance = 1e-8)
  exact_range_selected <- length_data_frame(
    lengths,
    with_length_range = c(50L, 70L),
    value = "fraction",
    denominator = "selected_bins"
  )
  expect_equal(exact_range_selected$length_bin_idx, 2L)
  expect_equal(exact_range_selected$fraction, 1, tolerance = 1e-8)

  left_edge_range <- length_data_frame(lengths, with_length_range = c(49L, 50L))
  expect_equal(left_edge_range$length_bin_idx, 1L)
  expect_equal(left_edge_range$count, 3, tolerance = 1e-8)

  right_edge_range <- length_data_frame(lengths, with_length_range = c(70L, 100L))
  expect_equal(right_edge_range$length_bin_idx, 3L)
  expect_equal(right_edge_range$count, 1, tolerance = 1e-8)

  overlap_range <- length_data_frame(lengths, with_length_range = c(49L, 71L))
  expect_equal(overlap_range$length_bin_idx, c(1L, 2L, 3L))
  expect_equal(overlap_range$count, c(3, 2, 1), tolerance = 1e-8)

  wide_selected <- length_data_frame(lengths, length_bin_idxs = c(3L, 1L), keep_wide = TRUE)
  expect_equal(names(wide_selected), c("count_70_100", "count_30_50"))
  expect_equal(wide_selected$count_70_100, 1, tolerance = 1e-8)
  expect_equal(wide_selected$count_30_50, 3, tolerance = 1e-8)
  wide_selected_fraction <- length_data_frame(
    lengths,
    length_bin_idxs = c(3L, 1L),
    value = "fraction",
    denominator = "selected_bins",
    keep_wide = TRUE
  )
  expect_equal(names(wide_selected_fraction), c("fraction_70_100", "fraction_30_50"))
  expect_equal(wide_selected_fraction$fraction_70_100, 0.25, tolerance = 1e-8)
  expect_equal(wide_selected_fraction$fraction_30_50, 0.75, tolerance = 1e-8)
  wide_empty <- length_data_frame(lengths, length_bin_idxs = integer(0), keep_wide = TRUE)
  expect_equal(dim(wide_empty), c(1L, 0L))
  expect_equal(names(wide_empty), character(0))
  expect_equal(length_data_frame(lengths, with_lengths = integer(0))$count, numeric(0))
  expect_equal(length_data_frame(lengths, length_bin_idxs = integer(0), value = "density")$density, numeric(0))
})

test_that("length-bin selectors reject invalid values", {
  path <- write_length_tsv_fixture(c(
    "count_30_50\tcount_50_70",
    "3\t2"
  ))

  lengths <- read_lengths(path)

  expect_error(length_counts_matrix(lengths, with_lengths = TRUE), "Fragment length must be")
  expect_error(length_counts_matrix(lengths, with_lengths = c(35L, 35L)), "with_lengths contains duplicate values")
  expect_error(length_counts_matrix(lengths, length_bin_idxs = TRUE), "length_bin_idxs must contain integer values")
  expect_error(length_counts_matrix(lengths, length_bin_idxs = c(1L, 1L)), "length_bin_idxs contains duplicate values")
  expect_error(
    length_data_frame(lengths, with_length_range = c(30L)),
    "with_length_range must be a numeric vector of two non-negative integer bp bounds"
  )
  expect_error(
    length_data_frame(lengths, with_length_range = c(30L, 50L, 70L)),
    "with_length_range must be a numeric vector of two non-negative integer bp bounds"
  )
  expect_error(
    length_data_frame(lengths, with_length_range = c(30.5, 50)),
    "with_length_range must be a numeric vector of two non-negative integer bp bounds"
  )
  expect_error(
    length_data_frame(lengths, with_length_range = c(-1L, 50L)),
    "with_length_range must be a numeric vector of two non-negative integer bp bounds"
  )
  expect_error(
    length_data_frame(lengths, with_length_range = c(50L, 50L)),
    "with_length_range start must be smaller than end"
  )
  expect_error(
    length_data_frame(lengths, value = "fraction", denominator = "selected"),
    "denominator must be one of"
  )
  expect_error(
    length_data_frame(lengths, value = "frac"),
    "value must be one of"
  )
})

test_that("windowed length counts expose window metadata and selected frames", {
  path <- write_length_tsv_fixture(c(
    "chrom\tstart\tend\tblacklisted_fraction\tcount_30\tcount_31_40",
    "chr1\t10\t20\t0.25\t12\t3",
    "chr2\t30\t45\t0\t0\t5"
  ))

  lengths <- read_lengths(path)

  expect_s3_class(lengths, "cfdnalab_windowed_length_counts")
  expect_output(print(lengths), "Mode: windowed", fixed = TRUE)
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

  selected <- length_data_frame(lengths, window_idxs = 2L)
  expect_equal(selected$window_idx, c(2L, 2L))
  expect_equal(selected$count, c(0, 5), tolerance = 1e-8)
  selected_range <- length_data_frame(lengths, window_idxs = 2L, with_length_range = c(31L, 35L))
  expect_equal(selected_range$window_idx, 2L)
  expect_equal(selected_range$length_bin_idx, 2L)
  expect_equal(selected_range$count, 5, tolerance = 1e-8)
  expect_equal(
    length_counts_matrix(lengths, window_idxs = 2L),
    matrix(c(0, 5), nrow = 1L, dimnames = list(NULL, c("count_30", "count_31_40"))),
    tolerance = 1e-8
  )
  expect_equal(
    length_counts_matrix(lengths, window_idxs = 2L, with_length_range = c(31L, 35L)),
    matrix(5, nrow = 1L, dimnames = list(NULL, "count_31_40")),
    tolerance = 1e-8
  )

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
  metadata_only <- length_data_frame(lengths, window_idxs = 2L, length_bin_idxs = integer(0), keep_wide = TRUE)
  expect_equal(names(metadata_only), c("window_idx", "chrom", "start", "end", "blacklisted_fraction"))
  expect_equal(metadata_only$window_idx, 2L)
})

test_that("windowed length counts keep numeric chromosome labels as character", {
  path <- write_length_tsv_fixture(c(
    "chrom\tstart\tend\tcount_30",
    "1\t10\t20\t5"
  ))

  lengths <- read_lengths(path)

  expect_identical(window_metadata(lengths)$chrom, "1")
})

test_that("windowed length selectors validate one-based indices", {
  path <- write_length_tsv_fixture(c(
    "chrom\tstart\tend\tcount_30\tcount_31_40",
    "chr1\t10\t20\t1\t2",
    "chr1\t30\t40\t3\t4"
  ))

  lengths <- read_lengths(path)

  expect_error(length_data_frame(lengths, window_idxs = 0L), "window_idxs contains values outside 1..2")
  expect_error(length_data_frame(lengths, window_idxs = 3L), "window_idxs contains values outside 1..2")
  expect_error(length_data_frame(lengths, window_idxs = 1.5), "window_idxs must contain integer values")
  expect_error(length_data_frame(lengths, window_idxs = c(1L, 1L)), "window_idxs contains duplicate values")
  expect_error(length_counts_matrix(lengths, window_idxs = c(1L, 1L)), "window_idxs contains duplicate values")
  expect_equal(length_data_frame(lengths, window_idxs = integer(0))$count, numeric(0))
})

test_that("grouped length counts expose group metadata and group selectors", {
  path <- write_length_tsv_fixture(c(
    "group_name\teligible_windows\tcount_30\tcount_31_40",
    "alpha\t2\t12\t3",
    "beta\t0\t0\t0"
  ))

  lengths <- read_lengths(path)

  expect_s3_class(lengths, "cfdnalab_grouped_length_counts")
  expect_output(print(lengths), "Mode: grouped", fixed = TRUE)
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

  alpha <- length_data_frame(lengths, groups = "alpha")
  expect_equal(alpha$group_idx, c(1L, 1L))
  expect_equal(alpha$count, c(12, 3), tolerance = 1e-8)

  beta_wide <- length_data_frame(lengths, group_idxs = 2L, value = "density", keep_wide = TRUE)
  expect_equal(names(beta_wide), c("group_idx", "group_name", "eligible_windows", "density_30", "density_31_40"))
  expect_true(all(is.na(beta_wide[c("density_30", "density_31_40")])))
  expect_equal(
    length_counts_matrix(lengths, group_idxs = 2L),
    matrix(c(0, 0), nrow = 1L, dimnames = list(NULL, c("count_30", "count_31_40"))),
    tolerance = 1e-8
  )
  alpha_length <- length_data_frame(lengths, groups = "alpha", with_lengths = 39L)
  expect_equal(alpha_length$length_bin_idx, 2L)
  expect_equal(alpha_length$count, 3, tolerance = 1e-8)
  selected_fraction <- length_data_frame(
    lengths,
    groups = c("alpha", "beta"),
    with_length_range = c(31L, 35L),
    value = "fraction",
    denominator = "selected_bins"
  )
  expect_equal(selected_fraction$group_idx, c(1L, 2L))
  expect_equal(selected_fraction$fraction, c(1, NA_real_), tolerance = 1e-8)
  metadata_only <- length_data_frame(lengths, groups = "alpha", length_bin_idxs = integer(0), keep_wide = TRUE)
  expect_equal(names(metadata_only), c("group_idx", "group_name", "eligible_windows"))
  expect_equal(metadata_only$group_name, "alpha")
})

test_that("grouped length selectors preserve requested row and bin order", {
  path <- write_length_tsv_fixture(c(
    "group_name\teligible_windows\tcount_30\tcount_31_40",
    "alpha\t2\t1\t2",
    "beta\t3\t3\t4",
    "gamma\t5\t5\t6"
  ))

  lengths <- read_lengths(path)

  selected_wide <- length_data_frame(
    lengths,
    groups = c("gamma", "alpha"),
    length_bin_idxs = c(2L, 1L),
    keep_wide = TRUE
  )
  expect_equal(selected_wide$group_idx, c(3L, 1L))
  expect_equal(selected_wide$group_name, c("gamma", "alpha"))
  expect_equal(names(selected_wide), c(
    "group_idx",
    "group_name",
    "eligible_windows",
    "count_31_40",
    "count_30"
  ))
  expect_equal(selected_wide$count_31_40, c(6, 2), tolerance = 1e-8)
  expect_equal(selected_wide$count_30, c(5, 1), tolerance = 1e-8)
  expect_equal(
    length_counts_matrix(lengths, groups = c("gamma", "alpha"), length_bin_idxs = c(2L, 1L)),
    matrix(c(6, 2, 5, 1), nrow = 2L, dimnames = list(NULL, c("count_31_40", "count_30"))),
    tolerance = 1e-8
  )
})

test_that("grouped length duplicate names only affect requested names", {
  path <- write_length_tsv_fixture(c(
    "group_name\teligible_windows\tcount_30",
    "alpha\t2\t12",
    "beta\t1\t5",
    "alpha\t3\t7"
  ))

  lengths <- read_lengths(path)

  expect_equal(group_idx(lengths, "beta"), 2L)
  expect_equal(length_data_frame(lengths, groups = "beta")$count, 5)
  expect_error(length_data_frame(lengths, groups = "alpha"), "Length-count group name is not unique")
})

test_that("grouped length selectors validate names and one-based indices", {
  path <- write_length_tsv_fixture(c(
    "group_name\teligible_windows\tblacklisted_fraction\tcount_30\tcount_31_40",
    "alpha\t2\t0.25\t12\t3",
    "beta\t1\t0\t0\t5",
    "gamma\t0\t1\t0\t0"
  ))

  lengths <- read_lengths(path)

  expect_error(length_data_frame(lengths, groups = "alpha", group_idxs = 1L), "Use either groups or group_idxs")
  expect_error(length_counts_matrix(lengths, groups = "alpha", group_idxs = 1L), "Use either groups or group_idxs")
  expect_error(length_data_frame(lengths, groups = "missing"), "Unknown length-count group name")
  expect_error(length_data_frame(lengths, group_idxs = 4L), "group_idxs contains values outside 1..3")
  expect_error(length_data_frame(lengths, groups = list("alpha", 1L)), "groups must contain character strings")
  expect_error(length_data_frame(lengths, groups = c("alpha", "alpha")), "groups contains duplicate values")
  expect_error(length_data_frame(lengths, group_idxs = c(1L, 1L)), "group_idxs contains duplicate values")
  expect_error(length_counts_matrix(lengths, group_idxs = c(1L, 1L)), "group_idxs contains duplicate values")
  expect_error(
    length_data_frame(lengths, with_lengths = 39L, with_length_range = c(31L, 35L)),
    "Use only one of with_lengths, with_length_range, or length_bin_idxs"
  )
  expect_error(
    length_data_frame(lengths, with_lengths = c(31L, 39L)),
    "with_lengths values must resolve to distinct length bins; 31 and 39 both resolve to length_bin_idx 2"
  )
  expect_error(
    length_data_frame(lengths, with_length_range = c(40L, 45L)),
    "with_length_range does not overlap any length bins"
  )
  expect_error(length_data_frame(lengths, value = "invalid"), "value must be one of")
  expect_equal(length_data_frame(lengths, groups = character(0))$count, numeric(0))
  expect_equal(length_data_frame(lengths, length_bin_idxs = integer(0))$count, numeric(0))
  expect_equal(length_data_frame(lengths, group_idxs = integer(0))$count, numeric(0))
  expect_equal(unique(length_data_frame(lengths, max_blacklisted_fraction = 0.25)$group_name), c("alpha", "beta"))
  expect_equal(length_data_frame(lengths, groups = c("beta", "alpha"))$group_idx, c(2L, 2L, 1L, 1L))
})

test_that("length TSV validation rejects ambiguous or unsupported shapes", {
  missing_path <- tempfile(fileext = ".tsv")
  expect_error(read_lengths(missing_path), "Length-count TSV does not exist")

  directory_path <- tempfile(fileext = ".tsv")
  dir.create(directory_path)
  expect_error(read_lengths(directory_path), "exists but is a directory")

  wrong_extension <- write_length_tsv_fixture(c(
    "count_30",
    "1"
  ))
  wrong_extension <- sub("\\.tsv$", ".csv", wrong_extension)
  file.create(wrong_extension)
  expect_error(read_lengths(wrong_extension), "must end in '.tsv' or '.tsv.zst'")

  wrong_gzip_extension <- sub("\\.csv$", ".tsv.gz", wrong_extension)
  file.create(wrong_gzip_extension)
  expect_error(read_lengths(wrong_gzip_extension), "must end in '.tsv' or '.tsv.zst'")

  expect_error(read_lengths(123), "must be a single path string")

  missing_zstd <- write_length_tsv_fixture(c(
    "chrom\tstart\tend\tcount_30",
    "chr1\t10\t20\t1"
  ))
  zst_path <- sub("\\.tsv$", ".tsv.zst", missing_zstd)
  file.copy(missing_zstd, zst_path)
  if (!nzchar(Sys.which("zstd"))) {
    expect_error(read_lengths(zst_path), "requires the zstd command-line tool")
  }

  if (nzchar(Sys.which("zstd"))) {
    source_tsv <- write_length_tsv_fixture(c(
      "count_30",
      "5"
    ))
    compressed_tsv <- sub("\\.tsv$", ".tsv.zst", source_tsv)
    system2(Sys.which("zstd"), c("-q", "-f", "-o", compressed_tsv, source_tsv))
    expect_equal(length_counts_vector(read_lengths(compressed_tsv)), c(count_30 = 5))
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

  duplicate_columns <- write_length_tsv_fixture(c(
    "count_30\tcount_30",
    "1\t2"
  ))
  expect_error(read_lengths(duplicate_columns), "column names must be unique")

  unsupported_metadata <- write_length_tsv_fixture(c(
    "window_idx\tchrom\tstart\tend\tcount_30",
    "1\tchr1\t10\t20\t1"
  ))
  expect_error(read_lengths(unsupported_metadata), "Could not infer length-count output mode")

  negative_count <- write_length_tsv_fixture(c(
    "count_30",
    "-1"
  ))
  expect_error(read_lengths(negative_count), "count_30 must contain finite non-negative values")

  nonfinite_count <- write_length_tsv_fixture(c(
    "count_30",
    "Inf"
  ))
  expect_error(read_lengths(nonfinite_count), "count_30 must contain finite non-negative values")

  negative_window_start <- write_length_tsv_fixture(c(
    "chrom\tstart\tend\tcount_30",
    "chr1\t-1\t20\t1"
  ))
  expect_error(read_lengths(negative_window_start), "start must contain non-negative integer values")
})

test_that("length-bin lookup rejects gaps and overlapping bins", {
  gapped <- read_lengths(write_length_tsv_fixture(c(
    "count_30_40\tcount_50_60",
    "1\t2"
  )))
  expect_error(length_bin_idx(gapped, 45L), "No length-count bin contains length 45")
  expect_error(length_bin_idx(gapped, -1L), "Fragment length must be a single non-negative integer")
  expect_error(length_bin_idx(gapped, 45.5), "Fragment length must be a single non-negative integer")
  expect_error(length_bin_idx(gapped, "45"), "Fragment length must be a single non-negative integer")

  overlapping <- read_lengths(write_length_tsv_fixture(c(
    "count_30_50\tcount_40_60",
    "1\t2"
  )))
  expect_error(length_bin_idx(overlapping, 45L), "Multiple length-count bins contain length 45")
})
