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

selected_path <- sparse_windowed_selected_motifs_end_zarr_path()
selected_root <- zarr::open_zarr(selected_path, read_only = TRUE)
selected_metadata <- jsonlite::fromJSON(file.path(selected_path, "zarr.json"), simplifyVector = FALSE)

selected_motifs <- decode_motif_ascii(read_cran_zarr_array(selected_root, "motif_ascii"))
stopifnot(identical(selected_metadata$attributes$cfdnalab_schema_version, 2L))
stopifnot(identical(selected_metadata$attributes$motif_axis_kind, "motif"))
stopifnot(identical(selected_motifs, c("GT_AC", "AC_GT")))
stopifnot(identical(
  labels_from_array_attributes(selected_path, "chromosome", "chromosome_name"),
  c("chr1", "chr2")
))
stopifnot(identical(as.integer(read_cran_zarr_array(selected_root, "sparse/shape")), c(3L, 2L)))
stopifnot(identical(as.integer(read_cran_zarr_array(selected_root, "sparse/row")), c(0L, 1L, 2L)))
stopifnot(identical(as.integer(read_cran_zarr_array(selected_root, "sparse/motif")), c(1L, 0L, 1L)))
stopifnot(isTRUE(all.equal(read_cran_zarr_array(selected_root, "sparse/count"), c(1, 1, 1))))

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

motif_grouped_path <- sparse_grouped_motif_group_end_zarr_path()
motif_grouped_root <- zarr::open_zarr(motif_grouped_path, read_only = TRUE)
motif_grouped_metadata <- jsonlite::fromJSON(file.path(motif_grouped_path, "zarr.json"), simplifyVector = FALSE)

stopifnot(identical(motif_grouped_metadata$attributes$cfdnalab_schema_version, 2L))
stopifnot(identical(motif_grouped_metadata$attributes$motif_axis_kind, "motif_group"))
stopifnot(!dir.exists(file.path(motif_grouped_path, "motif_ascii")))
stopifnot(identical(
  labels_from_array_attributes(motif_grouped_path, "motif_index", "motif_group"),
  c("left-hit", "right-hit")
))
stopifnot(identical(
  labels_from_array_attributes(motif_grouped_path, "group", "group_name"),
  c("beta", "alpha", "gamma")
))
stopifnot(identical(as.integer(read_cran_zarr_array(motif_grouped_root, "motif_index")), c(0L, 1L)))
stopifnot(identical(as.integer(read_cran_zarr_array(motif_grouped_root, "sparse/shape")), c(3L, 2L)))
stopifnot(identical(as.integer(read_cran_zarr_array(motif_grouped_root, "sparse/row")), c(0L, 0L, 1L)))
stopifnot(identical(as.integer(read_cran_zarr_array(motif_grouped_root, "sparse/motif")), c(0L, 1L, 1L)))
stopifnot(isTRUE(all.equal(read_cran_zarr_array(motif_grouped_root, "sparse/count"), c(2, 1, 1))))

wide_motif_grouped_path <- sparse_grouped_wide_motif_group_end_zarr_path()
wide_motif_grouped_root <- zarr::open_zarr(wide_motif_grouped_path, read_only = TRUE)
wide_motif_grouped_metadata <- jsonlite::fromJSON(file.path(wide_motif_grouped_path, "zarr.json"), simplifyVector = FALSE)

stopifnot(identical(wide_motif_grouped_metadata$attributes$cfdnalab_schema_version, 2L))
stopifnot(identical(wide_motif_grouped_metadata$attributes$motif_axis_kind, "motif_group"))
stopifnot(!dir.exists(file.path(wide_motif_grouped_path, "motif_ascii")))
stopifnot(identical(
  labels_from_array_attributes(wide_motif_grouped_path, "motif_index", "motif_group"),
  c("left-hit-wide", "right-hit-wide")
))
stopifnot(identical(as.integer(read_cran_zarr_array(wide_motif_grouped_root, "sparse/shape")), c(3L, 2L)))
stopifnot(identical(as.integer(read_cran_zarr_array(wide_motif_grouped_root, "sparse/row")), c(0L, 0L, 1L)))
stopifnot(identical(as.integer(read_cran_zarr_array(wide_motif_grouped_root, "sparse/motif")), c(0L, 1L, 1L)))
stopifnot(isTRUE(all.equal(read_cran_zarr_array(wide_motif_grouped_root, "sparse/count"), c(2, 1, 1))))
