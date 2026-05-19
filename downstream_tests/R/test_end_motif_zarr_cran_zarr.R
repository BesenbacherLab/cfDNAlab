source(file.path("downstream_tests", "R", "common.R"))

dense_path <- dense_global_end_zarr_path()
root <- zarr::open_zarr(dense_path, read_only = TRUE)

motifs <- decode_motif_ascii(read_cran_zarr_array(root, "motif_ascii"))
counts <- read_cran_zarr_array(root, "counts")

stopifnot(identical(motifs, c("_A", "_C", "_G", "_T")))
stopifnot(identical(dim(counts), c(1L, 4L)))
stopifnot(isTRUE(all.equal(counts[1, ], c(1, 0, 1, 0))))
stopifnot(identical(
  labels_from_array_attributes(dense_path, "row", "row_label"),
  "global"
))

windowed_path <- sparse_windowed_end_zarr_path()
windowed_root <- zarr::open_zarr(windowed_path, read_only = TRUE)

windowed_motifs <- decode_motif_ascii(read_cran_zarr_array(windowed_root, "motif_ascii"))
stopifnot(identical(windowed_motifs, c("_A", "_G")))
stopifnot(identical(
  labels_from_array_attributes(windowed_path, "chromosome", "chromosome_name"),
  c("chr1", "chr2")
))
stopifnot(identical(as.integer(read_cran_zarr_array(windowed_root, "row")), c(0L, 1L, 2L)))
stopifnot(identical(as.integer(read_cran_zarr_array(windowed_root, "row_chromosome")), c(0L, 0L, 1L)))
stopifnot(identical(as.integer(read_cran_zarr_array(windowed_root, "row_start_bp")), c(10L, 19L, 10L)))
stopifnot(identical(as.integer(read_cran_zarr_array(windowed_root, "row_end_bp")), c(11L, 20L, 11L)))
stopifnot(isTRUE(all.equal(read_cran_zarr_array(windowed_root, "blacklisted_fraction"), c(0, 0, 0))))
stopifnot(identical(as.integer(read_cran_zarr_array(windowed_root, "sparse/shape")), c(3L, 2L)))
stopifnot(identical(as.integer(read_cran_zarr_array(windowed_root, "sparse/row")), c(0L, 1L, 2L)))
stopifnot(identical(as.integer(read_cran_zarr_array(windowed_root, "sparse/motif")), c(1L, 0L, 1L)))
stopifnot(isTRUE(all.equal(read_cran_zarr_array(windowed_root, "sparse/count"), c(1, 1, 1))))

grouped_path <- sparse_grouped_end_zarr_path()
grouped_root <- zarr::open_zarr(grouped_path, read_only = TRUE)

grouped_motifs <- decode_motif_ascii(read_cran_zarr_array(grouped_root, "motif_ascii"))
stopifnot(identical(grouped_motifs, c("_A", "_G")))
stopifnot(identical(
  labels_from_array_attributes(grouped_path, "group", "group_name"),
  c("beta", "alpha", "gamma")
))
stopifnot(identical(
  labels_from_array_attributes(grouped_path, "sparse/sparse_dimension", "sparse_dimension_name"),
  c("row", "motif")
))
stopifnot(identical(as.integer(read_cran_zarr_array(grouped_root, "group")), c(0L, 1L, 2L)))
stopifnot(identical(as.integer(read_cran_zarr_array(grouped_root, "eligible_windows")), c(2L, 1L, 1L)))
stopifnot(isTRUE(all.equal(read_cran_zarr_array(grouped_root, "blacklisted_fraction"), c(0, 0, 0))))
stopifnot(identical(as.integer(read_cran_zarr_array(grouped_root, "sparse/shape")), c(3L, 2L)))
stopifnot(identical(as.integer(read_cran_zarr_array(grouped_root, "sparse/row")), c(0L, 0L, 1L)))
stopifnot(identical(as.integer(read_cran_zarr_array(grouped_root, "sparse/motif")), c(0L, 1L, 0L)))
stopifnot(isTRUE(all.equal(read_cran_zarr_array(grouped_root, "sparse/count"), c(1, 2, 1))))
