source(file.path("downstream_tests", "R", "common.R"))

library(testthat)

zstd_magic <- as.raw(c(0x28, 0xb5, 0x2f, 0xfd))

expect_zstd_frame <- function(path) {
  expect_identical(readBin(path, what = "raw", n = 4L), zstd_magic)
}

read_length_counts_table <- function(path) {
  zstd <- Sys.which("zstd")
  if (!nzchar(zstd)) {
    stop("Reading compressed length-count fixtures requires the zstd command-line tool")
  }
  data.table::fread(
    cmd = paste(shQuote(zstd), "-dc", shQuote(path)),
    data.table = FALSE,
    check.names = FALSE
  )
}

test_that("length-count TSV fixtures are zstd frames", {
  expect_zstd_frame(global_length_counts_path())
  expect_zstd_frame(windowed_length_counts_path())
  expect_zstd_frame(grouped_length_counts_path())
  expect_zstd_frame(windowed_length_counts_no_blacklist_path())
  expect_zstd_frame(grouped_length_counts_no_blacklist_path())
})

test_that("data.table reads global length-count fixtures", {
  global <- read_length_counts_table(global_length_counts_path())

  expect_identical(names(global), c("count_30_50", "count_50_70", "count_70_100"))
  expect_identical(as.integer(unlist(global[1, ], use.names = FALSE)), c(3L, 2L, 1L))
})

test_that("data.table reads windowed length-count fixtures", {
  windowed <- read_length_counts_table(windowed_length_counts_path())
  windowed_no_blacklist <- read_length_counts_table(windowed_length_counts_no_blacklist_path())

  expect_identical(names(windowed), c(
    "chrom",
    "start",
    "end",
    "blacklisted_fraction",
    "count_30_50",
    "count_50_70",
    "count_70_100"
  ))
  expect_identical(names(windowed_no_blacklist), c(
    "chrom",
    "start",
    "end",
    "count_30_50",
    "count_50_70",
    "count_70_100"
  ))
  expect_identical(windowed$chrom, rep("chr1", 4))
  expect_identical(windowed$start, c(0L, 100L, 200L, 300L))
  expect_identical(windowed$end, c(100L, 200L, 300L, 360L))
  expect_equal(windowed$blacklisted_fraction, c(0.04, 0.05, 0.1, 0.25))
  expect_identical(
    unname(as.matrix(windowed[c("count_30_50", "count_50_70", "count_70_100")])),
    matrix(
      c(2L, 0L, 0L, 0L, 2L, 0L, 0L, 0L, 1L, 1L, 0L, 0L),
      nrow = 4L,
      byrow = TRUE
    )
  )
})

test_that("data.table reads grouped length-count fixtures", {
  grouped <- read_length_counts_table(grouped_length_counts_path())
  grouped_no_blacklist <- read_length_counts_table(grouped_length_counts_no_blacklist_path())

  expect_identical(names(grouped), c(
    "group_name",
    "eligible_windows",
    "blacklisted_fraction",
    "count_30_50",
    "count_50_70",
    "count_70_100"
  ))
  expect_identical(names(grouped_no_blacklist), c(
    "group_name",
    "eligible_windows",
    "count_30_50",
    "count_50_70",
    "count_70_100"
  ))
  expect_identical(grouped$group_name, c("beta", "alpha", "gamma", "zero"))
  expect_identical(grouped$eligible_windows, c(2L, 1L, 1L, 1L))
  expect_equal(grouped$blacklisted_fraction, c(0.07, 0.05, 0.25, 0.333))
  expect_identical(
    unname(as.matrix(grouped[c("count_30_50", "count_50_70", "count_70_100")])),
    matrix(
      c(2L, 0L, 1L, 0L, 2L, 0L, 1L, 0L, 0L, 0L, 0L, 0L),
      nrow = 4L,
      byrow = TRUE
    )
  )
})
