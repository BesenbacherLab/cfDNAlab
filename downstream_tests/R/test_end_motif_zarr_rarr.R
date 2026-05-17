source(file.path("downstream_tests", "R", "common.R"))

read_array <- function(store_path, name, index) {
  Rarr::read_zarr_array(file.path(store_path, name), index = index)
}

read_sparse_as_dense <- function(store_path) {
  shape <- as.integer(read_array(store_path, file.path("sparse", "shape"), index = list(NULL)))
  row <- as.integer(read_array(store_path, file.path("sparse", "row"), index = list(NULL)))
  motif <- as.integer(read_array(store_path, file.path("sparse", "motif"), index = list(NULL)))
  count <- read_array(store_path, file.path("sparse", "count"), index = list(NULL))
  dense <- matrix(0, nrow = shape[1], ncol = shape[2])
  for (entry_index in seq_along(count)) {
    dense[row[entry_index] + 1L, motif[entry_index] + 1L] <- count[entry_index]
  }
  dense
}

dense_path <- dense_global_end_zarr_path()
dense_motifs <- decode_motif_ascii(
  read_array(dense_path, "motif_ascii", index = list(NULL, NULL))
)
dense_counts <- read_array(dense_path, "counts", index = list(NULL, NULL))

stopifnot(identical(dense_motifs, c("_A", "_C", "_G", "_T")))
stopifnot(identical(dim(dense_counts), c(1L, 4L)))
stopifnot(isTRUE(all.equal(dense_counts[1, ], c(1, 0, 1, 0))))
stopifnot(identical(
  labels_from_array_attributes(dense_path, "row", "row_label"),
  "global"
))

windowed_path <- sparse_windowed_end_zarr_path()
windowed_motifs <- decode_motif_ascii(
  read_array(windowed_path, "motif_ascii", index = list(NULL, NULL))
)
windowed_counts <- read_sparse_as_dense(windowed_path)
window_start_bp <- read_array(windowed_path, "row_start_bp", index = list(NULL))
window_end_bp <- read_array(windowed_path, "row_end_bp", index = list(NULL))

stopifnot(identical(windowed_motifs, c("_A", "_G")))
stopifnot(identical(dim(windowed_counts), c(2L, 2L)))
stopifnot(isTRUE(all.equal(windowed_counts, matrix(c(0, 1, 1, 0), nrow = 2L))))
stopifnot(identical(as.integer(window_start_bp), c(10L, 19L)))
stopifnot(identical(as.integer(window_end_bp), c(11L, 20L)))

grouped_path <- sparse_grouped_end_zarr_path()
grouped_motifs <- decode_motif_ascii(
  read_array(grouped_path, "motif_ascii", index = list(NULL, NULL))
)
grouped_counts <- read_sparse_as_dense(grouped_path)
group_names <- labels_from_array_attributes(grouped_path, "group", "group_name")
eligible_windows <- read_array(grouped_path, "eligible_windows", index = list(NULL))

stopifnot(identical(grouped_motifs, c("_A", "_G")))
stopifnot(identical(group_names, c("beta", "alpha", "gamma")))
stopifnot(identical(as.integer(eligible_windows), c(2L, 1L, 1L)))
stopifnot(identical(dim(grouped_counts), c(3L, 2L)))
stopifnot(isTRUE(all.equal(
  grouped_counts,
  matrix(c(1, 1, 0, 2, 0, 0), nrow = 3L)
)))
