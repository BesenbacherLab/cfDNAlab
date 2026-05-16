source(file.path("downstream_tests", "R", "common.R"))
store_path <- midpoint_zarr_path()

read_array <- function(name, index) {
  Rarr::read_zarr_array(file.path(store_path, name), index = index)
}

read_group_names_from_fallback <- function() {
  bytes <- read_array("group_name_utf8", index = list(NULL, NULL))
  nbytes <- read_array("group_name_nbytes", index = list(NULL))
  decode_group_names(bytes, nbytes)
}

counts <- read_array("counts", index = list(NULL, NULL, NULL))
stopifnot(identical(dim(counts), c(3L, 3L, 5L)))
stopifnot(isTRUE(all.equal(counts[1, 1, 1], 0)))
stopifnot(isTRUE(all.equal(counts[1, 2, 3], 1)))
stopifnot(isTRUE(all.equal(counts[1, 2, 4], 1)))
stopifnot(isTRUE(all.equal(counts[1, 2, 5], 0.5)))
stopifnot(isTRUE(all.equal(counts[1, 3, 3], 1)))
stopifnot(isTRUE(all.equal(counts[2, 2, 1], 1.5)))
stopifnot(isTRUE(all.equal(counts[2, 2, 2], 0.5)))
stopifnot(isTRUE(all.equal(counts[2, 3, 3], 0.5)))
stopifnot(isTRUE(all.equal(counts[3, 1, 3], 0.5)))
stopifnot(isTRUE(all.equal(counts[3, 1, 4], 0.5)))

group <- read_array("group", index = list(NULL))
length_bin <- read_array("length_bin", index = list(NULL))
position <- read_array("position", index = list(NULL))
stopifnot(identical(as.integer(group), c(0L, 1L, 2L)))
stopifnot(identical(as.integer(length_bin), c(0L, 1L, 2L)))
stopifnot(identical(as.integer(position), c(0L, 1L, 2L, 3L, 4L)))

group_names <- tryCatch(
  read_array("group_name", index = list(NULL)),
  error = function(error) read_group_names_from_fallback()
)
stopifnot(identical(as.character(group_names), c("alpha", "beta-site", "gamma_long")))
stopifnot(identical(read_group_names_from_fallback(), c("alpha", "beta-site", "gamma_long")))

group_index <- 1L
length_start_bp <- read_array("length_start_bp", index = list(NULL))
length_end_bp <- read_array("length_end_bp", index = list(NULL))
position_bin_start_bp <- read_array("position_bin_start_bp", index = list(NULL))
position_bin_end_bp <- read_array("position_bin_end_bp", index = list(NULL))
eligible_intervals <- read_array("eligible_intervals", index = list(NULL))
group_name_nbytes <- read_array("group_name_nbytes", index = list(NULL))
group_name_byte <- read_array("group_name_byte", index = list(NULL))

stopifnot(identical(as.integer(length_start_bp), c(30L, 50L, 70L)))
stopifnot(identical(as.integer(length_end_bp), c(50L, 70L, 100L)))
stopifnot(identical(as.integer(position_bin_start_bp), c(0L, 2L, 4L, 6L, 8L)))
stopifnot(identical(as.integer(position_bin_end_bp), c(2L, 4L, 6L, 8L, 10L)))
stopifnot(identical(as.integer(eligible_intervals), c(2L, 2L, 2L)))
stopifnot(identical(as.integer(group_name_nbytes), c(5L, 9L, 10L)))
stopifnot(identical(as.integer(group_name_byte), 0:9))

profile <- counts[group_index, , ]
rows <- do.call(rbind, lapply(seq_along(length_start_bp), function(length_index) {
  do.call(rbind, lapply(seq_along(position_bin_start_bp), function(position_index) {
    data.frame(
      group_name = as.character(group_names[group_index]),
      eligible_intervals = eligible_intervals[group_index],
      length_start_bp = length_start_bp[length_index],
      length_end_bp = length_end_bp[length_index],
      position_bin_start_bp = position_bin_start_bp[position_index],
      position_bin_end_bp = position_bin_end_bp[position_index],
      count = profile[length_index, position_index]
    )
  }))
}))

stopifnot(identical(nrow(rows), 15L))
stopifnot(isTRUE(all.equal(rows$count, c(
  0, 0, 0, 0, 0,
  0, 0, 1, 1, 0.5,
  0, 0, 1, 0, 0
))))
