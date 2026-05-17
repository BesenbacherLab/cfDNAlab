test_that("zero-based axis validation accepts empty and contiguous axes", {
  expect_true(cf_validate_axis(integer(), "empty_axis"))
  expect_true(cf_validate_axis(c(0L, 1L, 2L), "axis"))
})

test_that("zero-based axis validation rejects gaps and wrong starts", {
  expect_error(
    cf_validate_axis(c(1L, 2L), "axis"),
    "axis must contain contiguous zero-based indices",
    fixed = TRUE
  )
  expect_error(
    cf_validate_axis(c(0L, 2L), "axis"),
    "axis must contain contiguous zero-based indices",
    fixed = TRUE
  )
})

test_that("zero-based axis validation rejects coerced values", {
  expect_error(
    cf_validate_axis(c("0", "1"), "axis"),
    "axis must contain integer values",
    fixed = TRUE
  )
  expect_error(
    cf_validate_axis(c(0.1, 1.1), "axis"),
    "axis must contain integer values",
    fixed = TRUE
  )
})

test_that("index vector validation rejects non-integers and out-of-range indices", {
  expect_true(cf_validate_index_vector(c(0L, 2L), 3L, "sparse/row"))
  expect_error(
    cf_validate_index_vector(c(0L, 3L), 3L, "sparse/row"),
    "sparse/row contains an index outside 0..2",
    fixed = TRUE
  )
  expect_error(
    cf_validate_index_vector(c(0L, -1L), 3L, "sparse/row"),
    "sparse/row must contain non-negative integer values",
    fixed = TRUE
  )
  expect_error(
    cf_validate_index_vector(c(0L, 1.5), 3L, "sparse/row"),
    "sparse/row must contain integer values",
    fixed = TRUE
  )
  expect_error(
    cf_validate_index_vector(
      .Machine$integer.max + 1,
      .Machine$integer.max + 2,
      "sparse/row"
    ),
    "sparse/row values must fit in R integer range",
    fixed = TRUE
  )
})

test_that("interval and fraction validation reject invalid metadata", {
  expect_true(cf_validate_half_open_intervals(c(0L, 10L), c(5L, 12L), "start", "end"))
  expect_true(cf_validate_half_open_intervals(
    bit64::as.integer64(c(0, 10)),
    bit64::as.integer64(c(5, 12)),
    "start",
    "end"
  ))
  expect_error(
    cf_validate_half_open_intervals(c(0L, 10L), c(5L, 10L), "start", "end"),
    "start must be smaller than end",
    fixed = TRUE
  )
  expect_true(cf_validate_fraction_vector(c(0, 0.5, 1), "blacklisted_fraction"))
  expect_error(
    cf_validate_fraction_vector(c(0, 1.5), "blacklisted_fraction"),
    "blacklisted_fraction must contain finite fractions in 0..1",
    fixed = TRUE
  )
})

test_that("rank-1 vector reads preserve integer64 coordinates", {
  store <- list(
    get_node = function(path) {
      expect_equal(path, "/row_start_bp")
      list(
        read = function() {
          values <- bit64::as.integer64(c(10, 19))
          dim(values) <- 2L
          values
        }
      )
    }
  )

  values <- cf_read_vector(store, "row_start_bp", "test")

  expect_null(dim(values))
  expect_s3_class(values, "integer64")
  expect_equal(as.character(values), c("10", "19"))
  expect_true(cf_validate_half_open_intervals(
    values,
    bit64::as.integer64(c(11, 20)),
    "row_start_bp",
    "row_end_bp"
  ))
})

test_that("non-negative numeric validation rejects invalid counts", {
  expect_true(cf_validate_nonnegative_numeric_vector(c(0, 1.5), "sparse/count"))
  expect_error(
    cf_validate_nonnegative_numeric_vector(c(0, -1), "sparse/count"),
    "sparse/count must contain finite non-negative values",
    fixed = TRUE
  )
  expect_error(
    cf_validate_nonnegative_numeric_vector(c(0, NaN), "sparse/count"),
    "sparse/count must contain finite non-negative values",
    fixed = TRUE
  )
})

test_that("public R index validation rejects zero-based and coerced values", {
  expect_equal(cf_validate_r_index(1L, 3L, "group_idx"), 1L)
  expect_error(
    cf_validate_r_index(0L, 3L, "group_idx"),
    "group_idx 0 is outside 1..3",
    fixed = TRUE
  )
  expect_error(
    cf_validate_r_index(4L, 3L, "group_idx"),
    "group_idx 4 is outside 1..3",
    fixed = TRUE
  )
  expect_error(
    cf_validate_r_index("1", 3L, "group_idx"),
    "group_idx must be a single integer",
    fixed = TRUE
  )
  expect_error(
    cf_validate_r_index(Inf, 3L, "group_idx"),
    "group_idx must be a single integer",
    fixed = TRUE
  )
})

test_that("internal zero-based index conversion is explicit", {
  expect_equal(cf_validate_index0(0L, 3L, "group_idx0"), 0L)
  expect_equal(cf_r_index_to_index0(1L), 0L)
  expect_equal(cf_index0_to_r_index(0L), 1L)
  expect_error(
    cf_validate_index0(3L, 3L, "group_idx0"),
    "group_idx0 3 is outside 0..2",
    fixed = TRUE
  )
})

test_that("scalar string validation distinguishes missing and non-string values", {
  expect_equal(cf_validate_scalar_string("alpha", "group_name"), "alpha")
  expect_error(
    cf_validate_scalar_string(1L, "group_name"),
    "group_name must be a single character string",
    fixed = TRUE
  )
  expect_error(
    cf_validate_scalar_string(c("a", "b"), "group_name"),
    "group_name must be a single character string",
    fixed = TRUE
  )
})

test_that("motif ASCII decoding validates shape and byte width", {
  bytes <- matrix(
    as.integer(charToRaw("_A_C")),
    nrow = 2L,
    byrow = TRUE
  )

  expect_equal(cf_decode_motif_ascii(bytes, n_motifs = 2L, motif_width = 2L), c("_A", "_C"))
  expect_error(
    cf_decode_motif_ascii(bytes, n_motifs = 3L, motif_width = 2L),
    "motif_ascii row count",
    fixed = TRUE
  )
  expect_error(
    cf_decode_motif_ascii(bytes, n_motifs = 2L, motif_width = 3L),
    "motif_ascii column count",
    fixed = TRUE
  )
  bytes_with_invalid_value <- bytes
  bytes_with_invalid_value[[1L, 1L]] <- 128L
  expect_error(
    cf_decode_motif_ascii(bytes_with_invalid_value, n_motifs = 2L, motif_width = 2L),
    "ASCII byte values in 0..127",
    fixed = TRUE
  )
  bytes_with_infinite_value <- bytes
  bytes_with_infinite_value[[1L, 1L]] <- Inf
  expect_error(
    cf_decode_motif_ascii(bytes_with_infinite_value, n_motifs = 2L, motif_width = 2L),
    "ASCII byte values in 0..127",
    fixed = TRUE
  )
})
