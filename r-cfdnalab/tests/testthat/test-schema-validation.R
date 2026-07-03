write_json <- function(path, value) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(jsonlite::toJSON(value, auto_unbox = TRUE), path)
}

test_that("root attribute reader reports missing attributes", {
  store_path <- tempfile(fileext = ".zarr")
  dir.create(store_path)
  write_json(file.path(store_path, "zarr.json"), list(node_type = "group"))

  expect_error(
    cf_root_attributes(store_path),
    "root metadata is missing attributes",
    fixed = TRUE
  )
})

test_that("schema validation accepts supported schema version range", {
  attrs <- list(
    cfdnalab_schema = "midpoint_profiles",
    cfdnalab_schema_version = 1L
  )

  expect_true(cf_validate_schema(attrs, "midpoint_profiles", "midpoint profile"))
})

test_that("schema validation rejects wrong schema and unsupported versions", {
  expect_error(
    cf_validate_schema(
      list(cfdnalab_schema = "other", cfdnalab_schema_version = 1L),
      "midpoint_profiles",
      "midpoint profile"
    ),
    "Expected cfdnalab_schema",
    fixed = TRUE
  )
  expect_error(
    cf_validate_schema(
      list(cfdnalab_schema = "midpoint_profiles", cfdnalab_schema_version = 99L),
      "midpoint_profiles",
      "midpoint profile"
    ),
    "Unsupported midpoint profile schema version",
    fixed = TRUE
  )
  expect_error(
    cf_validate_schema(
      list(cfdnalab_schema = "midpoint_profiles", cfdnalab_schema_version = "1"),
      "midpoint_profiles",
      "midpoint profile"
    ),
    "Unsupported midpoint profile schema version",
    fixed = TRUE
  )
  expect_error(
    cf_validate_schema(
      list(cfdnalab_schema = "midpoint_profiles", cfdnalab_schema_version = 1.5),
      "midpoint_profiles",
      "midpoint profile"
    ),
    "Unsupported midpoint profile schema version",
    fixed = TRUE
  )
})

test_that("schema validation requires a per-schema supported version range", {
  expect_error(
    cf_validate_schema(
      list(cfdnalab_schema = "unregistered_schema", cfdnalab_schema_version = 1L),
      "unregistered_schema",
      "unregistered"
    ),
    "No supported schema-version range is registered",
    fixed = TRUE
  )
})

test_that("dimension-name validation reads nested array metadata", {
  store_path <- tempfile(fileext = ".zarr")
  write_json(
    file.path(store_path, "sparse", "row", "zarr.json"),
    list(
      node_type = "array",
      dimension_names = list("nnz"),
      attributes = list()
    )
  )

  expect_true(cf_validate_dimension_names(store_path, "sparse/row", "nnz"))
  expect_error(
    cf_validate_dimension_names(store_path, "sparse/row", "row"),
    "sparse/row dimensions must be row",
    fixed = TRUE
  )
})

test_that("dimension-name validation rejects non-character dimensions", {
  store_path <- tempfile(fileext = ".zarr")
  write_json(
    file.path(store_path, "row", "zarr.json"),
    list(
      node_type = "array",
      dimension_names = list(1),
      attributes = list()
    )
  )

  expect_error(
    cf_validate_dimension_names(store_path, "row", "row"),
    "row dimensions must be character strings",
    fixed = TRUE
  )
})

test_that("label reader validates label field and length", {
  store_path <- tempfile(fileext = ".zarr")
  write_json(
    file.path(store_path, "group", "zarr.json"),
    list(
      node_type = "array",
      dimension_names = list("group"),
      attributes = list(
        label_field = "group_name",
        labels = list("alpha", "beta")
      )
    )
  )

  expect_equal(cf_read_labels(store_path, "group", "group_name", 2L), c("alpha", "beta"))
  expect_error(
    cf_read_labels(store_path, "group", "chromosome_name", 2L),
    "metadata must declare label_field = chromosome_name",
    fixed = TRUE
  )
  expect_error(
    cf_read_labels(store_path, "group", "group_name", 3L),
    "labels length (2) does not match axis length (3)",
    fixed = TRUE
  )
})

test_that("label reader rejects non-character labels", {
  store_path <- tempfile(fileext = ".zarr")
  write_json(
    file.path(store_path, "group", "zarr.json"),
    list(
      node_type = "array",
      dimension_names = list("group"),
      attributes = list(
        label_field = "group_name",
        labels = list(1, 2)
      )
    )
  )

  expect_error(
    cf_read_labels(store_path, "group", "group_name", 2L),
    "labels must be character strings",
    fixed = TRUE
  )
})

test_that("label reader rejects control characters", {
  store_path <- tempfile(fileext = ".zarr")
  write_json(
    file.path(store_path, "group", "zarr.json"),
    list(
      node_type = "array",
      dimension_names = list("group"),
      attributes = list(
        label_field = "group_name",
        labels = list("alpha", "bad\nlabel")
      )
    )
  )

  expect_error(
    cf_read_labels(store_path, "group", "group_name", 2L),
    "labels must not contain control characters",
    fixed = TRUE
  )
})
