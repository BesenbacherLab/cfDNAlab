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
