midpoint_zarr_path <- function() {
  store_path <- Sys.getenv("CFDNALAB_MIDPOINT_ZARR")
  if (!nzchar(store_path)) {
    store_path <- file.path("downstream_tests", "tmp", "tiny.midpoint_profiles.zarr")
  }
  if (!dir.exists(store_path)) {
    stop(
      "Missing cfDNAlab-generated midpoint Zarr fixture: ",
      store_path,
      ". Generate it with the ignored generate_midpoint_zarr_fixture_with_cfdnalab integration test."
    )
  }
  store_path
}

decode_group_names <- function(bytes, nbytes) {
  vapply(seq_along(nbytes), function(group_index) {
    rawToChar(as.raw(bytes[group_index, seq_len(nbytes[group_index])]))
  }, character(1))
}
