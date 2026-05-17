midpoint_zarr_path <- function() {
  fixture_path(
    "CFDNALAB_MIDPOINT_ZARR",
    "tiny.midpoint_profiles.zarr",
    "midpoint"
  )
}

dense_global_end_zarr_path <- function() {
  fixture_path(
    "CFDNALAB_ENDS_DENSE_GLOBAL_ZARR",
    "tiny_dense_global.end_motifs.zarr",
    "dense global end-motif"
  )
}

sparse_windowed_end_zarr_path <- function() {
  fixture_path(
    "CFDNALAB_ENDS_SPARSE_WINDOWED_ZARR",
    "tiny_sparse_windowed.end_motifs.zarr",
    "sparse windowed end-motif"
  )
}

sparse_grouped_end_zarr_path <- function() {
  fixture_path(
    "CFDNALAB_ENDS_SPARSE_GROUPED_ZARR",
    "tiny_sparse_grouped.end_motifs.zarr",
    "sparse grouped end-motif"
  )
}

fixture_path <- function(env_var, default_name, label) {
  store_path <- Sys.getenv(env_var)
  if (!nzchar(store_path)) {
    store_path <- file.path("downstream_tests", "tmp", default_name)
  }
  if (!dir.exists(store_path)) {
    stop(
      "Missing cfDNAlab-generated ", label, " Zarr fixture: ",
      store_path,
      ". Generate it with the ignored downstream fixture integration tests."
    )
  }
  store_path
}

zarr_array_attributes <- function(store_path, array_name) {
  metadata_path <- file.path(store_path, array_name, "zarr.json")
  if (!file.exists(metadata_path)) {
    stop("Missing Zarr array metadata file: ", metadata_path)
  }
  metadata <- jsonlite::fromJSON(metadata_path, simplifyVector = FALSE)
  metadata$attributes
}

labels_from_array_attributes <- function(store_path, array_name, label_field) {
  attributes <- zarr_array_attributes(store_path, array_name)
  if (!identical(attributes$label_field, label_field)) {
    stop(
      "Expected label_field = ", label_field,
      " for array ", array_name,
      ", found ", attributes$label_field
    )
  }
  unlist(attributes$labels, use.names = FALSE)
}

decode_motif_ascii <- function(bytes) {
  apply(bytes, 1, function(row) {
    rawToChar(as.raw(as.integer(row)))
  })
}
