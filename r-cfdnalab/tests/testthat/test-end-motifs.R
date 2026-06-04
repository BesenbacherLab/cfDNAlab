test_that("dense global end motifs read from locally generated schema fixture", {
  ends <- read_end_motifs(make_dense_global_end_motif_zarr_fixture())

  expect_s3_class(ends, "cfdnalab_global_end_motif_counts")
  expect_equal(schema_version(ends), 1L)
  expect_equal(storage_mode(ends), "dense")
  expect_equal(row_mode(ends), "global")
  expect_output(print(ends), "<cfDNAlab end-motif counts>", fixed = TRUE)
  expect_equal(motifs(ends)$motif, c("_A", "_C", "_G", "_T"))
  expect_true(has_motif(ends, "_G"))
  expect_false(has_motif(ends, "_AA"))
  expect_equal(motif_idx(ends, "_G"), 3L)
  expect_equal(motifs(ends)$motif_idx, 1:4)
  expect_equal(
    dense_counts_vector(ends),
    c("_A" = 1, "_C" = 0, "_G" = 2.5, "_T" = 0),
    tolerance = 1e-8
  )
  expect_equal(
    dense_counts_matrix(ends),
    matrix(c(1, 0, 2.5, 0), nrow = 1L),
    tolerance = 1e-8
  )
  expect_equal(
    dense_counts_matrix(ends, motifs = "_G"),
    matrix(2.5, nrow = 1L),
    tolerance = 1e-8
  )
  expect_equal(
    dense_counts_matrix(ends, motif_idxs = c(3L, 1L)),
    matrix(c(2.5, 1), nrow = 1L),
    tolerance = 1e-8
  )
  expect_equal(
    as.matrix(sparse_counts_matrix(ends)),
    matrix(c(1, 0, 2.5, 0), nrow = 1L),
    tolerance = 1e-8
  )
  expect_equal(
    as.matrix(sparse_counts_matrix(ends, motifs = "_G")),
    matrix(2.5, nrow = 1L),
    tolerance = 1e-8
  )
  expect_equal(end_motif_data_frame(ends)$count, c(1, 0, 2.5, 0), tolerance = 1e-8)
  expect_equal(
    end_motif_data_frame(ends, motifs = "_G")$count,
    2.5,
    tolerance = 1e-8
  )
  expect_false(any(grepl("idx0|index0", names(end_motif_data_frame(ends)))))
})

test_that("sparse global end motifs can be densified explicitly", {
  ends <- read_end_motifs(make_sparse_global_end_motif_zarr_fixture())

  expect_s3_class(ends, "cfdnalab_global_end_motif_counts")
  expect_equal(storage_mode(ends), "sparse_coo")
  expect_equal(row_mode(ends), "global")
  expect_equal(
    as.matrix(sparse_counts_matrix(ends)),
    matrix(c(1.25, 0, 3.5, 0), nrow = 1L),
    tolerance = 1e-8
  )
  expect_equal(
    as.matrix(sparse_counts_matrix(ends, motifs = "_G")),
    matrix(3.5, nrow = 1L),
    tolerance = 1e-8
  )
  expect_error(dense_counts_vector(ends), "Use sparse_counts_matrix")
  expect_equal(
    dense_counts_vector(ends, allow_densify = TRUE),
    c("_A" = 1.25, "_C" = 0, "_G" = 3.5, "_T" = 0),
    tolerance = 1e-8
  )
  expect_equal(
    dense_counts_matrix(ends, motifs = c("_C", "_G"), allow_densify = TRUE),
    matrix(c(0, 3.5), nrow = 1L),
    tolerance = 1e-8
  )
  expect_equal(
    end_motif_data_frame(ends)$count,
    c(1.25, 3.5),
    tolerance = 1e-8
  )
  expect_equal(end_motif_data_frame(ends)$row_label, c("global", "global"))
  expect_equal(end_motif_data_frame(ends)$motif_idx, c(1L, 3L))
  expect_equal(
    end_motif_data_frame(ends, motifs = "_C")$count,
    numeric(0),
    tolerance = 1e-8
  )
  expect_equal(
    end_motif_data_frame(ends, densify = TRUE)$count,
    c(1.25, 0, 3.5, 0),
    tolerance = 1e-8
  )
  expect_equal(
    end_motif_data_frame(ends, motifs = "_G", densify = TRUE)$count,
    3.5,
    tolerance = 1e-8
  )
})

test_that("sparse end motifs without stored counts report a loader error", {
  expect_error(
    read_end_motifs(make_empty_sparse_end_motif_metadata_fixture()),
    "No end-motif counts are available",
    fixed = TRUE
  )
})

test_that("dense global motif-group end motifs load through motif-axis selectors", {
  ends <- read_end_motifs(make_dense_global_end_motif_group_zarr_fixture())

  expect_s3_class(ends, "cfdnalab_global_end_motif_counts")
  expect_equal(schema_version(ends), 2L)
  expect_equal(storage_mode(ends), "dense")
  expect_equal(row_mode(ends), "global")
  expect_equal(
    motifs(ends),
    data.frame(
      motif_idx = c(1L, 2L),
      motif = c("short", "group-two"),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_equal(motif_idx(ends, "group-two"), 2L)
  expect_true(has_motif(ends, "short"))
  expect_equal(
    dense_counts_vector(ends),
    stats::setNames(c(1.5, 3), c("short", "group-two")),
    tolerance = 1e-8
  )
  expect_equal(
    dense_counts_matrix(ends, motifs = "group-two"),
    matrix(3, nrow = 1L),
    tolerance = 1e-8
  )
  expect_equal(
    dense_counts_matrix(ends, motif_idxs = c(2L, 1L)),
    matrix(c(3, 1.5), nrow = 1L),
    tolerance = 1e-8
  )
  expect_equal(
    as.matrix(sparse_counts_matrix(ends, motifs = c("group-two", "short"))),
    matrix(c(3, 1.5), nrow = 1L),
    tolerance = 1e-8
  )
  group_frame <- end_motif_data_frame(ends, motifs = "group-two")
  expect_equal(names(group_frame), c("row_label", "motif_idx", "motif", "count"))
  expect_equal(group_frame$motif_idx, 2L)
  expect_equal(group_frame$motif, "group-two")
  expect_equal(group_frame$count, 3, tolerance = 1e-8)
  expect_false(any(grepl("idx0|index0", names(group_frame))))
  expect_error(
    end_motif_data_frame(ends, motifs = "_A"),
    "Unknown end-motif label",
    fixed = TRUE
  )
})

test_that("sparse global motif-group end motifs densify and keep motif columns", {
  ends <- read_end_motifs(make_sparse_global_end_motif_group_zarr_fixture())

  expect_s3_class(ends, "cfdnalab_global_end_motif_counts")
  expect_equal(schema_version(ends), 2L)
  expect_equal(storage_mode(ends), "sparse_coo")
  expect_equal(row_mode(ends), "global")
  expect_equal(
    motifs(ends),
    data.frame(
      motif_idx = c(1L, 2L),
      motif = c("short", "group-two"),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_equal(
    as.matrix(sparse_counts_matrix(ends)),
    matrix(c(1.5, 0), nrow = 1L),
    tolerance = 1e-8
  )
  expect_equal(
    as.matrix(sparse_counts_matrix(ends, motif_idxs = c(2L, 1L))),
    matrix(c(0, 1.5), nrow = 1L),
    tolerance = 1e-8
  )
  expect_error(dense_counts_matrix(ends), "Use sparse_counts_matrix")
  expect_equal(
    dense_counts_matrix(ends, allow_densify = TRUE),
    matrix(c(1.5, 0), nrow = 1L),
    tolerance = 1e-8
  )
  expect_equal(
    dense_counts_matrix(ends, motifs = "group-two", allow_densify = TRUE),
    matrix(0, nrow = 1L),
    tolerance = 1e-8
  )

  stored_frame <- end_motif_data_frame(ends)
  expect_equal(names(stored_frame), c("row_label", "motif_idx", "motif", "count"))
  expect_equal(stored_frame$motif_idx, 1L)
  expect_equal(stored_frame$motif, "short")
  expect_equal(stored_frame$count, 1.5, tolerance = 1e-8)

  dense_frame <- end_motif_data_frame(ends, densify = TRUE)
  expect_equal(names(dense_frame), c("row_label", "motif_idx", "motif", "count"))
  expect_equal(dense_frame$motif_idx, c(1L, 2L))
  expect_equal(dense_frame$motif, c("short", "group-two"))
  expect_equal(dense_frame$count, c(1.5, 0), tolerance = 1e-8)

  missing_group_frame <- end_motif_data_frame(ends, motifs = "group-two", densify = TRUE)
  expect_equal(missing_group_frame$motif_idx, 2L)
  expect_equal(missing_group_frame$motif, "group-two")
  expect_equal(missing_group_frame$count, 0, tolerance = 1e-8)
  expect_error(
    end_motif_data_frame(ends, motifs = "_A"),
    "Unknown end-motif label",
    fixed = TRUE
  )
})

test_that("dense windowed end motifs expose dense and sparse views", {
  ends <- read_end_motifs(make_dense_windowed_end_motif_zarr_fixture())

  expect_s3_class(ends, "cfdnalab_windowed_end_motif_counts")
  expect_equal(storage_mode(ends), "dense")
  expect_equal(row_mode(ends), "bed")
  expect_equal(
    dense_counts_matrix(ends),
    matrix(c(0, 2, 0, 1.5, 0, 4), nrow = 2L, byrow = TRUE),
    tolerance = 1e-8
  )
  expect_equal(
    dense_counts_matrix(ends, window_idxs = 2L, motifs = c("_T", "_A")),
    matrix(c(4, 1.5), nrow = 1L),
    tolerance = 1e-8
  )
  expect_equal(
    as.matrix(sparse_counts_matrix(ends)),
    matrix(c(0, 2, 0, 1.5, 0, 4), nrow = 2L, byrow = TRUE),
    tolerance = 1e-8
  )
  expect_equal(
    as.matrix(sparse_counts_matrix(ends, window_idxs = 2L, motifs = c("_T", "_A"))),
    matrix(c(4, 1.5), nrow = 1L),
    tolerance = 1e-8
  )
  expect_equal(
    end_motif_data_frame(ends, window_idxs = 1L)$count,
    c(0, 2, 0),
    tolerance = 1e-8
  )
  expect_equal(
    end_motif_data_frame(ends, window_idxs = 2L)$count,
    c(1.5, 0, 4),
    tolerance = 1e-8
  )
  expect_equal(nrow(end_motif_data_frame(ends, window_idxs = 2L, max_blacklisted_fraction = 0.1)), 0L)
  filtered_motif <- end_motif_data_frame(ends, motifs = "_A", max_blacklisted_fraction = 0.1)
  expect_equal(filtered_motif$window_idx, 1L)
  expect_equal(filtered_motif$count, 0, tolerance = 1e-8)
  expect_error(
    end_motif_data_frame(ends, window_idxs = 1L, max_blacklisted_fraction = 1.1),
    "max_blacklisted_fraction must be a single finite fraction in 0..1"
  )
  expect_error(
    end_motif_data_frame(ends, window_idxs = 0L),
    "window_idxs contains values outside 1..2",
    fixed = TRUE
  )
  expect_error(
    end_motif_data_frame(ends, window_idxs = 1L, window_index = 0L),
    "Unused argument(s): window_index",
    fixed = TRUE
  )
  expect_equal(
    end_motif_data_frame(ends, motifs = "_G")$count,
    c(2, 0),
    tolerance = 1e-8
  )
})

test_that("sparse windowed end motifs read from locally generated schema fixture", {
  ends <- read_end_motifs(make_sparse_windowed_end_motif_zarr_fixture())

  expect_s3_class(ends, "cfdnalab_windowed_end_motif_counts")
  expect_equal(storage_mode(ends), "sparse_coo")
  expect_equal(row_mode(ends), "bed")
  expect_error(
    window_metadata(read_end_motifs(make_dense_global_end_motif_zarr_fixture())),
    "no applicable method"
  )
  expect_equal(
    window_metadata(ends),
    data.frame(
      window_idx = c(1L, 2L, 3L),
      chrom = c("chr1", "chr1", "chr2"),
      start = c(10L, 20L, 30L),
      end = c(12L, 25L, 36L),
      blacklisted_fraction = c(0, 0.25, 0),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_equal(
    as.matrix(sparse_counts_matrix(ends)),
    matrix(c(0, 2, 0, 1.5, 0, 4, 0, 0, 3), nrow = 3L, byrow = TRUE),
    tolerance = 1e-8
  )
  expect_equal(
    as.matrix(sparse_counts_matrix(ends, window_idxs = c(3L, 1L), motifs = c("_T", "_G"))),
    matrix(c(3, 0, 0, 2), nrow = 2L, byrow = TRUE),
    tolerance = 1e-8
  )
  expect_error(dense_counts_matrix(ends), "Use sparse_counts_matrix")
  expect_equal(
    dense_counts_matrix(ends, allow_densify = TRUE),
    matrix(c(0, 2, 0, 1.5, 0, 4, 0, 0, 3), nrow = 3L, byrow = TRUE),
    tolerance = 1e-8
  )
  expect_equal(
    dense_counts_matrix(
      ends,
      window_idxs = c(3L, 1L),
      motifs = c("_T", "_G"),
      allow_densify = TRUE
    ),
    matrix(c(3, 0, 0, 2), nrow = 2L, byrow = TRUE),
    tolerance = 1e-8
  )

  motif_frame <- end_motif_data_frame(ends, motifs = "_T")
  expect_equal(motif_frame$window_idx, c(2L, 3L))
  expect_equal(motif_frame$count, c(4, 3), tolerance = 1e-8)
  expect_equal(
    end_motif_data_frame(ends)$count,
    c(2, 1.5, 4, 3),
    tolerance = 1e-8
  )
  expect_equal(end_motif_data_frame(ends)$window_idx, c(1L, 2L, 2L, 3L))
  expect_equal(end_motif_data_frame(ends)$motif_idx, c(2L, 1L, 3L, 3L))
  expect_equal(
    end_motif_data_frame(ends, window_idxs = 2L, densify = TRUE)$count,
    c(1.5, 0, 4),
    tolerance = 1e-8
  )
  expect_equal(nrow(end_motif_data_frame(ends, window_idxs = 2L, max_blacklisted_fraction = 0.1)), 0L)
  filtered_sparse_motif <- end_motif_data_frame(ends, motifs = "_T", max_blacklisted_fraction = 0.1)
  expect_equal(filtered_sparse_motif$window_idx, 3L)
  expect_equal(filtered_sparse_motif$count, 3, tolerance = 1e-8)

  window_frame <- end_motif_data_frame(ends, window_idxs = 1L)
  expect_equal(window_frame$motif, "_G")
  expect_equal(window_frame$count, 2, tolerance = 1e-8)
  expect_equal(
    end_motif_data_frame(ends, window_idxs = 2L)$count,
    c(1.5, 4),
    tolerance = 1e-8
  )
  expect_error(
    end_motif_data_frame(ends, window_idxs = 0L),
    "window_idxs contains values outside 1..3",
    fixed = TRUE
  )
})

test_that("size-mode end motifs use the windowed interface", {
  ends <- read_end_motifs(make_sparse_windowed_end_motif_zarr_fixture(row_mode = "size"))

  expect_s3_class(ends, "cfdnalab_windowed_end_motif_counts")
  expect_equal(row_mode(ends), "size")
  expect_equal(nrow(window_metadata(ends)), 3L)
  expect_equal(
    end_motif_data_frame(ends, window_idxs = 1L)$motif,
    "_G"
  )
})

test_that("dense grouped end motifs expose group helpers without densification flags", {
  ends <- read_end_motifs(make_dense_grouped_end_motif_zarr_fixture())

  expect_s3_class(ends, "cfdnalab_grouped_end_motif_counts")
  expect_equal(storage_mode(ends), "dense")
  expect_equal(row_mode(ends), "grouped_bed")
  expect_equal(group_idx(ends, "alpha"), 1L)
  expect_equal(
    dense_counts_matrix(ends, groups = "alpha", motifs = c("_G", "_A")),
    matrix(c(5, 1), nrow = 1L),
    tolerance = 1e-8
  )
  expect_equal(
    as.matrix(sparse_counts_matrix(ends, groups = "alpha", motifs = c("_G", "_A"))),
    matrix(c(5, 1), nrow = 1L),
    tolerance = 1e-8
  )
  expect_equal(
    end_motif_data_frame(ends, group_idxs = 1L)$count,
    c(1, 0, 5),
    tolerance = 1e-8
  )
  expect_error(
    end_motif_data_frame(ends, group_idxs = 0L),
    "group_idxs contains values outside 1..2",
    fixed = TRUE
  )
  expect_error(
    end_motif_data_frame(ends, group_idxs = 1L, group_index = 0L),
    "Unused argument(s): group_index",
    fixed = TRUE
  )
  expect_equal(
    end_motif_data_frame(ends, groups = "beta")$count,
    c(0, 0, 0),
    tolerance = 1e-8
  )
  expect_equal(nrow(end_motif_data_frame(ends, groups = "alpha", max_blacklisted_fraction = 0.1)), 0L)
  expect_equal(
    end_motif_data_frame(ends, motifs = "_G")$count,
    c(5, 0),
    tolerance = 1e-8
  )
})

test_that("sparse grouped end motifs read from locally generated schema fixture", {
  ends <- read_end_motifs(make_sparse_grouped_end_motif_zarr_fixture())

  expect_s3_class(ends, "cfdnalab_grouped_end_motif_counts")
  expect_equal(storage_mode(ends), "sparse_coo")
  expect_equal(row_mode(ends), "grouped_bed")
  expect_equal(group_idx(ends, "beta"), 2L)
  expect_equal(
    group_metadata(ends),
    data.frame(
      group_idx = c(1L, 2L),
      group_name = c("alpha", "beta"),
      eligible_windows = c(2L, 0L),
      blacklisted_fraction = c(0.125, 0),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_equal(
    as.matrix(sparse_counts_matrix(ends)),
    matrix(c(1, 0, 5, 0, 0, 0), nrow = 2L, byrow = TRUE),
    tolerance = 1e-8
  )
  expect_equal(
    as.matrix(sparse_counts_matrix(ends, groups = c("beta", "alpha"), motifs = c("_G", "_A"))),
    matrix(c(0, 0, 5, 1), nrow = 2L, byrow = TRUE),
    tolerance = 1e-8
  )
  expect_equal(
    dense_counts_matrix(
      ends,
      groups = c("beta", "alpha"),
      motifs = c("_G", "_A"),
      allow_densify = TRUE
    ),
    matrix(c(0, 0, 5, 1), nrow = 2L, byrow = TRUE),
    tolerance = 1e-8
  )

  alpha <- end_motif_data_frame(ends, groups = "alpha")
  expect_equal(alpha$motif, c("_A", "_G"))
  expect_equal(alpha$count, c(1, 5), tolerance = 1e-8)

  beta <- end_motif_data_frame(ends, groups = "beta")
  expect_equal(nrow(beta), 0L)
  expect_equal(
    end_motif_data_frame(ends, group_idxs = 1L)$count,
    c(1, 5),
    tolerance = 1e-8
  )
  expect_equal(nrow(end_motif_data_frame(ends, groups = "alpha", max_blacklisted_fraction = 0.1)), 0L)
  expect_equal(
    end_motif_data_frame(ends, motifs = "_G")$count,
    5,
    tolerance = 1e-8
  )
  expect_equal(nrow(end_motif_data_frame(ends, motifs = "_G", max_blacklisted_fraction = 0.1)), 0L)
  expect_equal(
    end_motif_data_frame(ends)$count,
    c(1, 5),
    tolerance = 1e-8
  )
  expect_false(any(grepl("idx0|index0", names(end_motif_data_frame(ends)))))
  expect_equal(
    end_motif_data_frame(ends, group_idxs = 1L, densify = TRUE)$count,
    c(1, 0, 5),
    tolerance = 1e-8
  )
  expect_equal(end_motif_data_frame(ends, groups = "beta", densify = TRUE)$count, c(0, 0, 0))
})

test_that("dense global end motifs load as data frames and vectors", {
  testthat::skip_if_not_installed("zarr")

  ends <- read_end_motifs(dense_global_end_fixture_path())

  expect_s3_class(ends, "cfdnalab_global_end_motif_counts")
  expect_equal(schema_version(ends), 1L)
  expect_equal(storage_mode(ends), "dense")
  expect_equal(row_mode(ends), "global")
  expect_equal(
    motifs(ends),
    data.frame(
      motif_idx = c(1L, 2L, 3L, 4L),
      motif = c("_A", "_C", "_G", "_T"),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_true(has_motif(ends, "_A"))
  expect_false(has_motif(ends, "_AA"))
  expect_equal(motif_idx(ends, "_G"), 3L)
  expect_equal(
    dense_counts_vector(ends),
    c("_A" = 1, "_C" = 0, "_G" = 1, "_T" = 0),
    tolerance = 1e-8
  )
  expect_equal(
    end_motif_data_frame(ends)$count,
    c(1, 0, 1, 0),
    tolerance = 1e-8
  )
  expect_equal(
    dense_counts_matrix(ends),
    matrix(c(1, 0, 1, 0), nrow = 1),
    tolerance = 1e-8
  )
  expect_equal(
    end_motif_data_frame(ends, motifs = "_A")$count,
    1,
    tolerance = 1e-8
  )
})

test_that("sparse windowed end motifs stay sparse unless densification is requested", {
  testthat::skip_if_not_installed("zarr")

  ends <- read_end_motifs(sparse_windowed_end_fixture_path())

  expect_s3_class(ends, "cfdnalab_windowed_end_motif_counts")
  expect_equal(schema_version(ends), 1L)
  expect_equal(storage_mode(ends), "sparse_coo")
  expect_equal(row_mode(ends), "bed")
  expect_equal(motifs(ends)$motif, c("_A", "_G"))

  expect_equal(
    window_metadata(ends),
    data.frame(
      window_idx = c(1L, 2L, 3L),
      chrom = c("chr1", "chr1", "chr2"),
      start = c(10L, 19L, 10L),
      end = c(11L, 20L, 11L),
      blacklisted_fraction = c(0, 0, 0),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_identical(window_metadata(ends)$start, c(10L, 19L, 10L))
  expect_identical(window_metadata(ends)$end, c(11L, 20L, 11L))

  sparse_counts <- sparse_counts_matrix(ends)
  expect_s4_class(sparse_counts, "sparseMatrix")
  expect_equal(as.matrix(sparse_counts), matrix(c(0, 1, 1, 0, 0, 1), nrow = 3, byrow = TRUE))
  expect_error(dense_counts_matrix(ends), "Use sparse_counts_matrix")
  expect_equal(
    dense_counts_matrix(ends, allow_densify = TRUE),
    matrix(c(0, 1, 1, 0, 0, 1), nrow = 3, byrow = TRUE),
    tolerance = 1e-8
  )

  motif_frame <- end_motif_data_frame(ends, motifs = "_G")
  expect_equal(motif_frame$window_idx, c(1L, 3L))
  expect_equal(motif_frame$motif, c("_G", "_G"))
  expect_equal(motif_frame$count, c(1, 1))

  window_frame <- end_motif_data_frame(ends, window_idxs = 1L)
  expect_equal(window_frame$window_idx, 1L)
  expect_equal(window_frame$motif, "_G")
  expect_equal(window_frame$count, 1)
})

test_that("sparse grouped end motifs expose group metadata and stored rows", {
  testthat::skip_if_not_installed("zarr")

  ends <- read_end_motifs(sparse_grouped_end_fixture_path())

  expect_s3_class(ends, "cfdnalab_grouped_end_motif_counts")
  expect_equal(schema_version(ends), 1L)
  expect_equal(storage_mode(ends), "sparse_coo")
  expect_equal(row_mode(ends), "grouped_bed")
  expect_equal(group_idx(ends, "alpha"), 2L)
  expect_equal(
    group_metadata(ends),
    data.frame(
      group_idx = c(1L, 2L, 3L),
      group_name = c("beta", "alpha", "gamma"),
      eligible_windows = c(2L, 1L, 1L),
      blacklisted_fraction = c(0, 0, 0),
      stringsAsFactors = FALSE
    ),
    ignore_attr = TRUE
  )
  expect_identical(group_metadata(ends)$eligible_windows, c(2L, 1L, 1L))

  sparse_counts <- sparse_counts_matrix(ends)
  expect_equal(
    as.matrix(sparse_counts),
    matrix(c(1, 2, 1, 0, 0, 0), nrow = 3, byrow = TRUE),
    tolerance = 1e-8
  )

  rows <- end_motif_data_frame(ends)
  expect_equal(rows$group_idx, c(1L, 1L, 2L))
  expect_equal(rows$motif_idx, c(1L, 2L, 1L))
  expect_equal(rows$motif, c("_A", "_G", "_A"))
  expect_equal(rows$count, c(1, 2, 1), tolerance = 1e-8)

  beta <- end_motif_data_frame(ends, groups = "beta")
  expect_equal(beta$group_name, c("beta", "beta"))
  expect_equal(beta$motif, c("_A", "_G"))
  expect_equal(beta$count, c(1, 2), tolerance = 1e-8)

  dense_beta <- end_motif_data_frame(ends, groups = "beta", densify = TRUE)
  expect_equal(dense_beta$count, c(1, 2), tolerance = 1e-8)
})
