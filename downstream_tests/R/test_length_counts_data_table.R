source(file.path("downstream_tests", "R", "common.R"))

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

global <- read_length_counts_table(global_length_counts_path())
stopifnot(identical(names(global), c("count_30_50", "count_50_70", "count_70_100")))
stopifnot(identical(as.integer(unlist(global[1, ], use.names = FALSE)), c(3L, 2L, 1L)))

windowed <- read_length_counts_table(windowed_length_counts_path())
windowed_no_blacklist <- read_length_counts_table(windowed_length_counts_no_blacklist_path())
stopifnot(identical(names(windowed), c(
  "chrom",
  "start",
  "end",
  "blacklisted_fraction",
  "count_30_50",
  "count_50_70",
  "count_70_100"
)))
stopifnot(identical(names(windowed_no_blacklist), c(
  "chrom",
  "start",
  "end",
  "count_30_50",
  "count_50_70",
  "count_70_100"
)))
stopifnot(identical(windowed$chrom, rep("chr1", 4)))
stopifnot(identical(windowed$start, c(0L, 100L, 200L, 300L)))
stopifnot(identical(windowed$end, c(100L, 200L, 300L, 360L)))
stopifnot(isTRUE(all.equal(windowed$blacklisted_fraction, c(0.04, 0.05, 0.1, 0.25))))
stopifnot(identical(as.matrix(windowed[c("count_30_50", "count_50_70", "count_70_100")]), matrix(
  c(2L, 0L, 0L, 0L, 2L, 0L, 0L, 0L, 1L, 1L, 0L, 0L),
  nrow = 4L,
  byrow = TRUE
)))

grouped <- read_length_counts_table(grouped_length_counts_path())
grouped_no_blacklist <- read_length_counts_table(grouped_length_counts_no_blacklist_path())
stopifnot(identical(names(grouped), c(
  "group_name",
  "eligible_windows",
  "blacklisted_fraction",
  "count_30_50",
  "count_50_70",
  "count_70_100"
)))
stopifnot(identical(names(grouped_no_blacklist), c(
  "group_name",
  "eligible_windows",
  "count_30_50",
  "count_50_70",
  "count_70_100"
)))
stopifnot(identical(grouped$group_name, c("beta", "alpha", "gamma", "zero")))
stopifnot(identical(grouped$eligible_windows, c(2L, 1L, 1L, 1L)))
stopifnot(isTRUE(all.equal(grouped$blacklisted_fraction, c(0.07, 0.05, 0.25, 0.333))))
stopifnot(identical(as.matrix(grouped[c("count_30_50", "count_50_70", "count_70_100")]), matrix(
  c(2L, 0L, 1L, 0L, 2L, 0L, 1L, 0L, 0L, 0L, 0L, 0L),
  nrow = 4L,
  byrow = TRUE
)))
