fixture_path <- function(default_name, env_var) {
  from_env <- Sys.getenv(env_var)
  candidates <- if (nzchar(from_env)) {
    from_env
  } else {
    c(
      file.path("..", "downstream_tests", "tmp", default_name),
      file.path("downstream_tests", "tmp", default_name)
    )
  }

  for (candidate in candidates) {
    if (dir.exists(candidate)) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }

  testthat::skip(
    paste0(
      "Missing cfDNAlab-generated Zarr fixture ",
      default_name,
      ". Generate downstream_tests/tmp fixtures first."
    )
  )
}

midpoint_fixture_path <- function() {
  fixture_path("tiny.midpoint_profiles.zarr", "CFDNALAB_MIDPOINT_ZARR")
}

dense_global_end_fixture_path <- function() {
  fixture_path("tiny_dense_global.end_motifs.zarr", "CFDNALAB_ENDS_DENSE_GLOBAL_ZARR")
}

sparse_windowed_end_fixture_path <- function() {
  fixture_path("tiny_sparse_windowed.end_motifs.zarr", "CFDNALAB_ENDS_SPARSE_WINDOWED_ZARR")
}

sparse_grouped_end_fixture_path <- function() {
  fixture_path("tiny_sparse_grouped.end_motifs.zarr", "CFDNALAB_ENDS_SPARSE_GROUPED_ZARR")
}

dense_global_ref_kmer_fixture_path <- function() {
  fixture_path("tiny_ref_kmers_dense_global.ref_kmers.zarr", "CFDNALAB_REF_KMERS_DENSE_GLOBAL_ZARR")
}

sparse_windowed_ref_kmer_fixture_path <- function() {
  fixture_path("tiny_ref_kmers_sparse_windowed.ref_kmers.zarr", "CFDNALAB_REF_KMERS_SPARSE_WINDOWED_ZARR")
}

sparse_grouped_ref_kmer_fixture_path <- function() {
  fixture_path("tiny_ref_kmers_sparse_grouped.ref_kmers.zarr", "CFDNALAB_REF_KMERS_SPARSE_GROUPED_ZARR")
}

local_zarr_store_path <- function(name) {
  file.path(tempdir(), paste0(name, "-", basename(tempfile()), ".zarr"))
}

write_fixture_json <- function(path, value) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(value, path, auto_unbox = TRUE, pretty = TRUE)
}

patch_zarr_metadata <- function(
  store_path,
  array_name = NULL,
  dimension_names = NULL,
  attributes = NULL,
  fill_value = NULL
) {
  metadata_path <- if (is.null(array_name)) {
    file.path(store_path, "zarr.json")
  } else {
    do.call(
      file.path,
      as.list(c(store_path, strsplit(array_name, "/", fixed = TRUE)[[1L]], "zarr.json"))
    )
  }

  metadata <- jsonlite::fromJSON(metadata_path, simplifyVector = FALSE)
  if (!is.null(dimension_names)) {
    metadata$dimension_names <- as.list(dimension_names)
  }
  if (!is.null(attributes)) {
    metadata$attributes <- attributes
  }
  if (!is.null(fill_value)) {
    metadata$fill_value <- fill_value
  }
  write_fixture_json(metadata_path, metadata)
}

write_minimal_array_metadata <- function(store_path, array_name, shape) {
  metadata_path <- do.call(
    file.path,
    as.list(c(store_path, strsplit(array_name, "/", fixed = TRUE)[[1L]], "zarr.json"))
  )
  write_fixture_json(
    metadata_path,
    list(
      zarr_format = 3L,
      node_type = "array",
      shape = as.list(shape)
    )
  )
}

add_fixture_array <- function(
  zarr_store,
  store_path,
  group_path,
  name,
  values,
  data_type,
  dimension_names,
  attributes = list()
) {
  shape <- dim(values)
  if (is.null(shape)) {
    shape <- length(values)
  }
  definition <- zarr::define_array(data_type, shape)
  array <- zarr_store$add_array(group_path, name, definition)
  array$write(values)
  array_path <- if (identical(group_path, "/")) {
    name
  } else {
    paste0(sub("^/", "", group_path), "/", name)
  }
  patch_zarr_metadata(
    store_path = store_path,
    array_name = array_path,
    dimension_names = dimension_names,
    attributes = attributes,
    fill_value = fixture_fill_value(data_type)
  )
  array
}

fixture_fill_value <- function(data_type) {
  # Match the Rust writers: fill values sit outside each public array's valid
  # data domain so readers do not convert real zeroes into missing values.
  switch(
    data_type,
    int32 = -1L,
    float64 = -1,
    uint8 = 255L,
    NULL
  )
}

motif_ascii_matrix <- function(labels) {
  width <- max(nchar(labels, type = "bytes"))
  bytes <- matrix(0L, nrow = length(labels), ncol = width)
  for (label_idx in seq_along(labels)) {
    label_bytes <- as.integer(charToRaw(labels[[label_idx]]))
    bytes[label_idx, seq_along(label_bytes)] <- label_bytes
  }
  bytes
}

ref_kmer_root_attributes <- function(
  storage_mode,
  row_mode,
  motif_axis_kind = "motif",
  kmer_size = 2L,
  canonical = FALSE,
  orientation = "both",
  all_motifs = FALSE,
  assign_by = "count-overlap",
  value_units = "reference_kmer_frequency",
  count_units = "reference_kmer_count",
  row_scaling_factor_array = "row_scaling_factor",
  count_reconstruction = "reference_kmer_count = frequency * row_scaling_factor[row]",
  schema = "ref_kmer_frequencies",
  schema_version = 2L
) {
  attributes <- list(
    cfdnalab_schema = schema,
    cfdnalab_schema_version = schema_version,
    storage_mode = storage_mode,
    row_mode = row_mode,
    motif_axis_kind = motif_axis_kind,
    value_units = value_units,
    count_units = count_units,
    row_scaling_factor_array = row_scaling_factor_array,
    count_reconstruction = count_reconstruction,
    kmer_size = kmer_size,
    canonical = canonical,
    orientation = orientation,
    all_motifs = all_motifs,
    assign_by = assign_by
  )
  if (identical(storage_mode, "dense")) {
    attributes$primary_array <- "frequencies"
  } else {
    attributes$primary_group <- "sparse"
    attributes$sparse_format <- "coo"
    attributes$sparse_indices_base <- 0L
  }
  attributes
}

add_reference_contig_footprint_array <- function(zarr_store, store_path) {
  footprint_json <- jsonlite::toJSON(
    list(
      list(name = "chr2", length = 100L),
      list(name = "chr10", length = 120L)
    ),
    auto_unbox = TRUE
  )
  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "reference_contig_footprint_json",
    as.integer(charToRaw(footprint_json)),
    "uint8",
    "json_byte"
  )
}

add_ref_kmer_motif_axis <- function(zarr_store, store_path, motifs) {
  add_fixture_array(zarr_store, store_path, "/", "motif_index", seq_along(motifs) - 1L, "int32", "motif")
  add_fixture_array(zarr_store, store_path, "/", "motif_byte", seq_len(nchar(motifs[[1L]], type = "bytes")) - 1L, "int32", "motif_byte")
  add_fixture_array(zarr_store, store_path, "/", "motif_ascii", motif_ascii_matrix(motifs), "uint8", c("motif", "motif_byte"))
}

add_ref_kmer_window_metadata <- function(
  zarr_store,
  store_path,
  chromosome_names = c("chr2", "chr10"),
  row_chromosome = c(0L, 1L),
  row_start_bp = c(10L, 40L),
  row_end_bp = c(20L, 60L),
  blacklisted_fraction = c(0.25, 0)
) {
  row_idx0 <- seq_along(row_start_bp) - 1L
  chromosome_idx0 <- seq_along(chromosome_names) - 1L
  add_fixture_array(zarr_store, store_path, "/", "row", row_idx0, "int32", "row")
  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "chromosome",
    chromosome_idx0,
    "int32",
    "chromosome",
    list(label_field = "chromosome_name", labels = as.list(chromosome_names))
  )
  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "row_chromosome",
    row_chromosome,
    "int32",
    "row"
  )
  add_fixture_array(zarr_store, store_path, "/", "row_start_bp", row_start_bp, "int32", "row")
  add_fixture_array(zarr_store, store_path, "/", "row_end_bp", row_end_bp, "int32", "row")
  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "blacklisted_fraction",
    blacklisted_fraction,
    "float64",
    "row"
  )
}

make_midpoint_zarr_fixture <- function() {
  testthat::skip_if_not_installed("zarr")

  store_path <- local_zarr_store_path("r-midpoint-fixture")
  zarr_store <- zarr::create_zarr(store_path)
  patch_zarr_metadata(
    store_path,
    attributes = list(
      cfdnalab_schema = "midpoint_profiles",
      cfdnalab_schema_version = 1L
    )
  )

  counts <- array(0, dim = c(2L, 2L, 4L))
  counts[1L, 1L, ] <- c(0, 1.5, 0, 2.25)
  counts[1L, 2L, ] <- c(3, 0, 4.5, 0)
  counts[2L, 1L, ] <- c(0.25, 0, 0.75, 1)
  counts[2L, 2L, ] <- c(5, 0, 0, 6.5)

  add_fixture_array(zarr_store, store_path, "/", "counts", counts, "float64", c("group", "length_bin", "position"))
  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "group",
    c(0L, 1L),
    "int32",
    "group",
    list(label_field = "group_name", labels = list("A", "long_group"))
  )
  add_fixture_array(zarr_store, store_path, "/", "eligible_intervals", c(1L, 3L), "int32", "group")
  add_fixture_array(zarr_store, store_path, "/", "length_bin", c(0L, 1L), "int32", "length_bin")
  add_fixture_array(zarr_store, store_path, "/", "length_start_bp", c(30L, 60L), "int32", "length_bin")
  add_fixture_array(zarr_store, store_path, "/", "length_end_bp", c(60L, 90L), "int32", "length_bin")
  add_fixture_array(zarr_store, store_path, "/", "position", c(0L, 1L, 2L, 3L), "int32", "position")
  add_fixture_array(zarr_store, store_path, "/", "position_bin_start_bp", c(0L, 5L, 10L, 15L), "int32", "position")
  add_fixture_array(zarr_store, store_path, "/", "position_bin_end_bp", c(5L, 10L, 15L, 20L), "int32", "position")

  store_path
}

make_dense_global_end_motif_zarr_fixture <- function(
  motifs = c("_A", "_C", "_G", "_T"),
  counts = matrix(c(1, 0, 2.5, 0), nrow = 1L)
) {
  testthat::skip_if_not_installed("zarr")

  store_path <- local_zarr_store_path("r-dense-global-end-fixture")
  zarr_store <- zarr::create_zarr(store_path)
  patch_zarr_metadata(
    store_path,
    attributes = list(
      cfdnalab_schema = "end_motif_counts",
      cfdnalab_schema_version = 1L,
      storage_mode = "dense",
      row_mode = "global"
    )
  )

  add_fixture_array(zarr_store, store_path, "/", "motif_index", seq_along(motifs) - 1L, "int32", "motif")
  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "motif_byte",
    seq_len(nchar(motifs[[1L]], type = "bytes")) - 1L,
    "int32",
    "motif_byte"
  )
  add_fixture_array(zarr_store, store_path, "/", "motif_ascii", motif_ascii_matrix(motifs), "uint8", c("motif", "motif_byte"))
  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "row",
    0L,
    "int32",
    "row",
    list(label_field = "row_label", labels = list("global"))
  )
  add_fixture_array(zarr_store, store_path, "/", "counts", counts, "float64", c("row", "motif"))

  store_path
}

make_dense_global_end_motif_group_zarr_fixture <- function() {
  testthat::skip_if_not_installed("zarr")

  store_path <- local_zarr_store_path("r-dense-global-end-motif-group-fixture")
  zarr_store <- zarr::create_zarr(store_path)
  patch_zarr_metadata(
    store_path,
    attributes = list(
      cfdnalab_schema = "end_motif_counts",
      cfdnalab_schema_version = 2L,
      storage_mode = "dense",
      row_mode = "global",
      motif_axis_kind = "motif_group"
    )
  )

  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "motif_index",
    0:1,
    "int32",
    "motif",
    list(label_field = "motif_group", labels = list("short", "group-two"))
  )
  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "row",
    0L,
    "int32",
    "row",
    list(label_field = "row_label", labels = list("global"))
  )
  add_fixture_array(zarr_store, store_path, "/", "counts", matrix(c(1.5, 3), nrow = 1L), "float64", c("row", "motif"))

  store_path
}

make_sparse_windowed_end_motif_zarr_fixture <- function(row_mode = "bed") {
  testthat::skip_if_not_installed("zarr")

  store_path <- local_zarr_store_path("r-sparse-windowed-end-fixture")
  zarr_store <- zarr::create_zarr(store_path)
  zarr_store$add_group("/", "sparse")
  patch_zarr_metadata(
    store_path,
    attributes = list(
      cfdnalab_schema = "end_motif_counts",
      cfdnalab_schema_version = 1L,
      storage_mode = "sparse_coo",
      row_mode = row_mode
    )
  )

  motifs <- c("_A", "_G", "_T")
  add_fixture_array(zarr_store, store_path, "/", "motif_index", 0:2, "int32", "motif")
  add_fixture_array(zarr_store, store_path, "/", "motif_byte", c(0L, 1L), "int32", "motif_byte")
  add_fixture_array(zarr_store, store_path, "/", "motif_ascii", motif_ascii_matrix(motifs), "uint8", c("motif", "motif_byte"))
  add_fixture_array(zarr_store, store_path, "/", "row", 0:2, "int32", "row")
  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "chromosome",
    0:1,
    "int32",
    "chromosome",
    list(label_field = "chromosome_name", labels = list("chr1", "chr2"))
  )
  add_fixture_array(zarr_store, store_path, "/", "row_chromosome", c(0L, 0L, 1L), "int32", "row")
  add_fixture_array(zarr_store, store_path, "/", "row_start_bp", c(10L, 20L, 30L), "int32", "row")
  add_fixture_array(zarr_store, store_path, "/", "row_end_bp", c(12L, 25L, 36L), "int32", "row")
  add_fixture_array(zarr_store, store_path, "/", "blacklisted_fraction", c(0, 0.25, 0), "float64", "row")
  add_fixture_array(zarr_store, store_path, "/sparse", "row", c(0L, 1L, 1L, 2L), "int32", "nnz")
  add_fixture_array(zarr_store, store_path, "/sparse", "motif", c(1L, 0L, 2L, 2L), "int32", "nnz")
  add_fixture_array(zarr_store, store_path, "/sparse", "count", c(2, 1.5, 4, 3), "float64", "nnz")
  add_fixture_array(zarr_store, store_path, "/sparse", "shape", c(3L, 3L), "int32", "sparse_dimension")
  add_fixture_array(
    zarr_store,
    store_path,
    "/sparse",
    "sparse_dimension",
    0:1,
    "int32",
    "sparse_dimension",
    list(label_field = "sparse_dimension_name", labels = list("row", "motif"))
  )

  store_path
}

make_sparse_global_end_motif_zarr_fixture <- function() {
  testthat::skip_if_not_installed("zarr")

  store_path <- local_zarr_store_path("r-sparse-global-end-fixture")
  zarr_store <- zarr::create_zarr(store_path)
  zarr_store$add_group("/", "sparse")
  patch_zarr_metadata(
    store_path,
    attributes = list(
      cfdnalab_schema = "end_motif_counts",
      cfdnalab_schema_version = 1L,
      storage_mode = "sparse_coo",
      row_mode = "global"
    )
  )

  motifs <- c("_A", "_C", "_G", "_T")
  add_fixture_array(zarr_store, store_path, "/", "motif_index", 0:3, "int32", "motif")
  add_fixture_array(zarr_store, store_path, "/", "motif_byte", c(0L, 1L), "int32", "motif_byte")
  add_fixture_array(zarr_store, store_path, "/", "motif_ascii", motif_ascii_matrix(motifs), "uint8", c("motif", "motif_byte"))
  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "row",
    0L,
    "int32",
    "row",
    list(label_field = "row_label", labels = list("global"))
  )
  add_fixture_array(zarr_store, store_path, "/sparse", "row", c(0L, 0L), "int32", "nnz")
  add_fixture_array(zarr_store, store_path, "/sparse", "motif", c(0L, 2L), "int32", "nnz")
  add_fixture_array(zarr_store, store_path, "/sparse", "count", c(1.25, 3.5), "float64", "nnz")
  add_fixture_array(zarr_store, store_path, "/sparse", "shape", c(1L, 4L), "int32", "sparse_dimension")
  add_fixture_array(
    zarr_store,
    store_path,
    "/sparse",
    "sparse_dimension",
    0:1,
    "int32",
    "sparse_dimension",
    list(label_field = "sparse_dimension_name", labels = list("row", "motif"))
  )

  store_path
}

make_empty_sparse_end_motif_metadata_fixture <- function() {
  store_path <- local_zarr_store_path("r-empty-sparse-end-fixture")
  write_fixture_json(
    file.path(store_path, "zarr.json"),
    list(
      zarr_format = 3L,
      node_type = "group",
      attributes = list(
        cfdnalab_schema = "end_motif_counts",
        cfdnalab_schema_version = 2L,
        storage_mode = "sparse_coo",
        row_mode = "global",
        motif_axis_kind = "motif"
      )
    )
  )
  write_minimal_array_metadata(store_path, "motif_index", c(3L))
  write_minimal_array_metadata(store_path, "sparse/count", c(0L))
  store_path
}

make_sparse_global_end_motif_group_zarr_fixture <- function() {
  testthat::skip_if_not_installed("zarr")

  store_path <- local_zarr_store_path("r-sparse-global-end-motif-group-fixture")
  zarr_store <- zarr::create_zarr(store_path)
  zarr_store$add_group("/", "sparse")
  patch_zarr_metadata(
    store_path,
    attributes = list(
      cfdnalab_schema = "end_motif_counts",
      cfdnalab_schema_version = 2L,
      storage_mode = "sparse_coo",
      row_mode = "global",
      motif_axis_kind = "motif_group"
    )
  )

  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "motif_index",
    0:1,
    "int32",
    "motif",
    list(label_field = "motif_group", labels = list("short", "group-two"))
  )
  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "row",
    0L,
    "int32",
    "row",
    list(label_field = "row_label", labels = list("global"))
  )
  add_fixture_array(zarr_store, store_path, "/sparse", "row", 0L, "int32", "nnz")
  add_fixture_array(zarr_store, store_path, "/sparse", "motif", 0L, "int32", "nnz")
  add_fixture_array(zarr_store, store_path, "/sparse", "count", 1.5, "float64", "nnz")
  add_fixture_array(zarr_store, store_path, "/sparse", "shape", c(1L, 2L), "int32", "sparse_dimension")
  add_fixture_array(
    zarr_store,
    store_path,
    "/sparse",
    "sparse_dimension",
    0:1,
    "int32",
    "sparse_dimension",
    list(label_field = "sparse_dimension_name", labels = list("row", "motif"))
  )

  store_path
}

make_dense_windowed_end_motif_zarr_fixture <- function() {
  testthat::skip_if_not_installed("zarr")

  store_path <- local_zarr_store_path("r-dense-windowed-end-fixture")
  zarr_store <- zarr::create_zarr(store_path)
  patch_zarr_metadata(
    store_path,
    attributes = list(
      cfdnalab_schema = "end_motif_counts",
      cfdnalab_schema_version = 1L,
      storage_mode = "dense",
      row_mode = "bed"
    )
  )

  motifs <- c("_A", "_G", "_T")
  add_fixture_array(zarr_store, store_path, "/", "motif_index", 0:2, "int32", "motif")
  add_fixture_array(zarr_store, store_path, "/", "motif_byte", c(0L, 1L), "int32", "motif_byte")
  add_fixture_array(zarr_store, store_path, "/", "motif_ascii", motif_ascii_matrix(motifs), "uint8", c("motif", "motif_byte"))
  add_fixture_array(zarr_store, store_path, "/", "row", 0:1, "int32", "row")
  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "chromosome",
    0:1,
    "int32",
    "chromosome",
    list(label_field = "chromosome_name", labels = list("chr1", "chr2"))
  )
  add_fixture_array(zarr_store, store_path, "/", "row_chromosome", c(0L, 1L), "int32", "row")
  add_fixture_array(zarr_store, store_path, "/", "row_start_bp", c(10L, 30L), "int32", "row")
  add_fixture_array(zarr_store, store_path, "/", "row_end_bp", c(12L, 36L), "int32", "row")
  add_fixture_array(zarr_store, store_path, "/", "blacklisted_fraction", c(0, 0.125), "float64", "row")
  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "counts",
    matrix(c(0, 2, 0, 1.5, 0, 4), nrow = 2L, byrow = TRUE),
    "float64",
    c("row", "motif")
  )

  store_path
}

make_dense_grouped_end_motif_zarr_fixture <- function() {
  testthat::skip_if_not_installed("zarr")

  store_path <- local_zarr_store_path("r-dense-grouped-end-fixture")
  zarr_store <- zarr::create_zarr(store_path)
  patch_zarr_metadata(
    store_path,
    attributes = list(
      cfdnalab_schema = "end_motif_counts",
      cfdnalab_schema_version = 1L,
      storage_mode = "dense",
      row_mode = "grouped_bed"
    )
  )

  motifs <- c("_A", "_C", "_G")
  add_fixture_array(zarr_store, store_path, "/", "motif_index", 0:2, "int32", "motif")
  add_fixture_array(zarr_store, store_path, "/", "motif_byte", c(0L, 1L), "int32", "motif_byte")
  add_fixture_array(zarr_store, store_path, "/", "motif_ascii", motif_ascii_matrix(motifs), "uint8", c("motif", "motif_byte"))
  add_fixture_array(zarr_store, store_path, "/", "row", 0:1, "int32", "row")
  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "group",
    0:1,
    "int32",
    "row",
    list(label_field = "group_name", labels = list("alpha", "beta"))
  )
  add_fixture_array(zarr_store, store_path, "/", "eligible_windows", c(2L, 0L), "int32", "row")
  add_fixture_array(zarr_store, store_path, "/", "blacklisted_fraction", c(0.125, 0), "float64", "row")
  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "counts",
    matrix(c(1, 0, 5, 0, 0, 0), nrow = 2L, byrow = TRUE),
    "float64",
    c("row", "motif")
  )

  store_path
}

make_sparse_grouped_end_motif_zarr_fixture <- function() {
  testthat::skip_if_not_installed("zarr")

  store_path <- local_zarr_store_path("r-sparse-grouped-end-fixture")
  zarr_store <- zarr::create_zarr(store_path)
  zarr_store$add_group("/", "sparse")
  patch_zarr_metadata(
    store_path,
    attributes = list(
      cfdnalab_schema = "end_motif_counts",
      cfdnalab_schema_version = 1L,
      storage_mode = "sparse_coo",
      row_mode = "grouped_bed"
    )
  )

  motifs <- c("_A", "_C", "_G")
  add_fixture_array(zarr_store, store_path, "/", "motif_index", 0:2, "int32", "motif")
  add_fixture_array(zarr_store, store_path, "/", "motif_byte", c(0L, 1L), "int32", "motif_byte")
  add_fixture_array(zarr_store, store_path, "/", "motif_ascii", motif_ascii_matrix(motifs), "uint8", c("motif", "motif_byte"))
  add_fixture_array(zarr_store, store_path, "/", "row", 0:1, "int32", "row")
  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "group",
    0:1,
    "int32",
    "row",
    list(label_field = "group_name", labels = list("alpha", "beta"))
  )
  add_fixture_array(zarr_store, store_path, "/", "eligible_windows", c(2L, 0L), "int32", "row")
  add_fixture_array(zarr_store, store_path, "/", "blacklisted_fraction", c(0.125, 0), "float64", "row")
  add_fixture_array(zarr_store, store_path, "/sparse", "row", c(0L, 0L), "int32", "nnz")
  add_fixture_array(zarr_store, store_path, "/sparse", "motif", c(0L, 2L), "int32", "nnz")
  add_fixture_array(zarr_store, store_path, "/sparse", "count", c(1, 5), "float64", "nnz")
  add_fixture_array(zarr_store, store_path, "/sparse", "shape", c(2L, 3L), "int32", "sparse_dimension")
  add_fixture_array(
    zarr_store,
    store_path,
    "/sparse",
    "sparse_dimension",
    0:1,
    "int32",
    "sparse_dimension",
    list(label_field = "sparse_dimension_name", labels = list("row", "motif"))
  )

  store_path
}

make_dense_windowed_ref_kmer_zarr_fixture <- function(
  motifs = c("AA", "AC", "GT"),
  canonical = FALSE,
  frequencies = matrix(c(0.25, 0, 0.75, 0.5, 0.5, 0), nrow = 2L, byrow = TRUE),
  row_scaling_factor = c(4, 2),
  chromosome_names = c("chr2", "chr10"),
  row_chromosome = c(0L, 1L),
  row_start_bp = c(10L, 40L),
  row_end_bp = c(20L, 60L),
  blacklisted_fraction = c(0.25, 0),
  root_attributes = NULL
) {
  testthat::skip_if_not_installed("zarr")

  store_path <- local_zarr_store_path("r-dense-windowed-ref-kmer-fixture")
  zarr_store <- zarr::create_zarr(store_path)
  if (is.null(root_attributes)) {
    root_attributes <- ref_kmer_root_attributes(
      storage_mode = "dense",
      row_mode = "bed",
      kmer_size = nchar(motifs[[1L]], type = "bytes"),
      canonical = canonical
    )
  }
  patch_zarr_metadata(store_path, attributes = root_attributes)

  add_ref_kmer_motif_axis(zarr_store, store_path, motifs)
  add_ref_kmer_window_metadata(
    zarr_store,
    store_path,
    chromosome_names = chromosome_names,
    row_chromosome = row_chromosome,
    row_start_bp = row_start_bp,
    row_end_bp = row_end_bp,
    blacklisted_fraction = blacklisted_fraction
  )
  add_fixture_array(zarr_store, store_path, "/", "row_scaling_factor", row_scaling_factor, "float64", "row")
  add_reference_contig_footprint_array(zarr_store, store_path)
  add_fixture_array(zarr_store, store_path, "/", "frequencies", frequencies, "float64", c("row", "motif"))

  store_path
}

make_sparse_grouped_ref_kmer_zarr_fixture <- function(
  motifs = c("AA", "AC", "GT"),
  sparse_row = c(0L, 0L, 1L),
  sparse_motif = c(0L, 2L, 1L),
  sparse_frequency = c(0.25, 0.75, 1),
  group_labels = list("A", "long_group", "empty"),
  sparse_shape = c(length(group_labels), length(motifs)),
  row_scaling_factor = c(4, 2, 0),
  eligible_windows = c(1L, 2L, 0L),
  blacklisted_fraction = c(0, 0.125, 0),
  sparse_dimension_labels = list("row", "motif")
) {
  testthat::skip_if_not_installed("zarr")

  row_idx0 <- seq_along(group_labels) - 1L
  store_path <- local_zarr_store_path("r-sparse-grouped-ref-kmer-fixture")
  zarr_store <- zarr::create_zarr(store_path)
  zarr_store$add_group("/", "sparse")
  patch_zarr_metadata(
    store_path,
    attributes = ref_kmer_root_attributes(
      storage_mode = "sparse_coo",
      row_mode = "grouped_bed",
      kmer_size = nchar(motifs[[1L]], type = "bytes")
    )
  )

  add_ref_kmer_motif_axis(zarr_store, store_path, motifs)
  add_fixture_array(zarr_store, store_path, "/", "row", row_idx0, "int32", "row")
  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "group",
    row_idx0,
    "int32",
    "row",
    list(label_field = "group_name", labels = group_labels)
  )
  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "eligible_windows",
    eligible_windows,
    "int32",
    "row"
  )
  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "blacklisted_fraction",
    blacklisted_fraction,
    "float64",
    "row"
  )
  add_fixture_array(zarr_store, store_path, "/", "row_scaling_factor", row_scaling_factor, "float64", "row")
  add_reference_contig_footprint_array(zarr_store, store_path)
  add_fixture_array(zarr_store, store_path, "/sparse", "row", sparse_row, "int32", "nnz")
  add_fixture_array(zarr_store, store_path, "/sparse", "motif", sparse_motif, "int32", "nnz")
  add_fixture_array(zarr_store, store_path, "/sparse", "frequency", sparse_frequency, "float64", "nnz")
  add_fixture_array(zarr_store, store_path, "/sparse", "shape", sparse_shape, "int32", "sparse_dimension")
  add_fixture_array(
    zarr_store,
    store_path,
    "/sparse",
    "sparse_dimension",
    0:1,
    "int32",
    "sparse_dimension",
    list(label_field = "sparse_dimension_name", labels = sparse_dimension_labels)
  )

  store_path
}

make_dense_global_ref_kmer_zarr_fixture <- function(
  motifs = c("A", "G", "T"),
  frequencies = matrix(c(1 / 3, 1 / 6, 1 / 2), nrow = 1L),
  row_scaling_factor = 6
) {
  testthat::skip_if_not_installed("zarr")

  store_path <- local_zarr_store_path("r-dense-global-ref-kmer-fixture")
  zarr_store <- zarr::create_zarr(store_path)
  patch_zarr_metadata(
    store_path,
    attributes = ref_kmer_root_attributes(
      storage_mode = "dense",
      row_mode = "global",
      kmer_size = nchar(motifs[[1L]], type = "bytes")
    )
  )

  add_ref_kmer_motif_axis(zarr_store, store_path, motifs)
  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "row",
    0L,
    "int32",
    "row",
    list(label_field = "row_label", labels = list("global"))
  )
  add_fixture_array(zarr_store, store_path, "/", "row_scaling_factor", row_scaling_factor, "float64", "row")
  add_reference_contig_footprint_array(zarr_store, store_path)
  add_fixture_array(zarr_store, store_path, "/", "frequencies", frequencies, "float64", c("row", "motif"))

  store_path
}

make_dense_global_ref_kmer_group_zarr_fixture <- function() {
  testthat::skip_if_not_installed("zarr")

  store_path <- local_zarr_store_path("r-dense-global-ref-kmer-group-fixture")
  zarr_store <- zarr::create_zarr(store_path)
  patch_zarr_metadata(
    store_path,
    attributes = ref_kmer_root_attributes(
      storage_mode = "dense",
      row_mode = "global",
      motif_axis_kind = "motif_group"
    )
  )

  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "motif_index",
    0:1,
    "int32",
    "motif",
    list(label_field = "motif_group", labels = list("left", "right"))
  )
  add_fixture_array(
    zarr_store,
    store_path,
    "/",
    "row",
    0L,
    "int32",
    "row",
    list(label_field = "row_label", labels = list("global"))
  )
  add_fixture_array(zarr_store, store_path, "/", "row_scaling_factor", 4, "float64", "row")
  add_reference_contig_footprint_array(zarr_store, store_path)
  add_fixture_array(zarr_store, store_path, "/", "frequencies", matrix(c(0.25, 0.75), nrow = 1L), "float64", c("row", "motif"))

  store_path
}
