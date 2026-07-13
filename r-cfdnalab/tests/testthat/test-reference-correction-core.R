# Keep the counts, frequencies, and rational expectations identical to the Rust
# and Python core correction tests. R-facing motif indices remain one-based.
shared_reference_correction_end_rows <- function() {
  data.frame(
    row_label = rep("global", 4L),
    motif_idx = seq_len(4L),
    motif = c("A_C", "A_G", "T_C", "T_G"),
    count = c(2, 4, 6, 8),
    stringsAsFactors = FALSE
  )
}

shared_reference_correction_ref_rows <- function() {
  data.frame(
    row_label = rep("global", 4L),
    reference_motif = c("AC", "AG", "TC", "TG"),
    reference_frequency = c(1 / 8, 1 / 8, 1 / 4, 1 / 2),
    stringsAsFactors = FALSE
  )
}

shared_reference_correction_mode <- function(mode, side_labels = character()) {
  list(
    mode = mode,
    outside_width = 1L,
    inside_width = 1L,
    side_labels = side_labels
  )
}

add_shared_corrected_frequencies <- function(corrected) {
  cf_add_corrected_frequency(corrected, "row_label", "error")
}

test_that("joint core uses full motif frequencies", {
  corrected <- cf_exact_reference_corrected_rows(
    list(motif_axis_kind = "motif"),
    shared_reference_correction_end_rows(),
    shared_reference_correction_ref_rows(),
    "row_label",
    "error"
  )
  corrected <- add_shared_corrected_frequencies(corrected)

  # Four positive reference motifs make the uniform frequency 1/4. Relative
  # to uniform, frequencies [1/8, 1/8, 1/4, 1/2] give correction factors
  # [1/2, 1/2, 1, 2] for [AC, AG, TC, TG]. Dividing original counts
  # [2, 4, 6, 8] by those factors gives [4, 8, 6, 4]. Their total is 22, so
  # dividing each corrected count by 22 gives [2/11, 4/11, 3/11, 2/11].
  expect_identical(corrected$motif_idx, seq_len(4L))
  expect_equal(corrected$reference_denominator, c(1 / 2, 1 / 2, 1, 2))
  expect_equal(corrected$corrected_count, c(4, 8, 6, 4))
  expect_equal(corrected$corrected_frequency, c(2 / 11, 4 / 11, 3 / 11, 2 / 11))
})

test_that("split core multiplies outside and inside denominators", {
  corrected <- cf_split_reference_corrected_rows(
    shared_reference_correction_end_rows(),
    shared_reference_correction_ref_rows(),
    "row_label",
    shared_reference_correction_mode("split"),
    "error"
  )
  corrected <- add_shared_corrected_frequencies(corrected)

  # Two positive labels on each side make each side's uniform frequency 1/2.
  # Outside frequencies A=1/4 and T=3/4 give factors 1/2 and 3/2. Inside
  # frequencies C=3/8 and G=5/8 give factors 3/4 and 5/4. Multiplying matching
  # side factors gives [3/8, 5/8, 9/8, 15/8] for [A_C, A_G, T_C, T_G].
  # Dividing original counts [2, 4, 6, 8] by those factors gives
  # [16/3, 32/5, 16/3, 64/15]. Their total is 64/3, so normalization gives
  # frequencies [1/4, 3/10, 1/4, 1/5].
  expect_identical(corrected$motif_idx, seq_len(4L))
  expect_equal(corrected$reference_denominator, c(3 / 8, 5 / 8, 9 / 8, 15 / 8))
  expect_equal(corrected$corrected_count, c(16 / 3, 32 / 5, 16 / 3, 64 / 15))
  expect_equal(corrected$corrected_frequency, c(1 / 4, 3 / 10, 1 / 4, 1 / 5))
})

test_that("outside core aggregates counts before correction", {
  end_rows <- shared_reference_correction_end_rows()
  end_rows$optional_metadata <- NA_character_
  output_columns <- names(end_rows)
  end_rows$.cfdnalab_row_order <- 1L
  corrected <- cf_side_reference_corrected_rows(
    end_rows,
    shared_reference_correction_ref_rows(),
    "row_label",
    shared_reference_correction_mode("outside", c("A_", "T_")),
    output_columns,
    "error"
  )
  corrected <- add_shared_corrected_frequencies(corrected)

  # Counts aggregate to A_=2+4=6 and T_=6+8=14. Two positive outside labels
  # make the uniform frequency 1/2. Relative to uniform, frequencies A=1/4
  # and T=3/4 give factors 1/2 and 3/2. Dividing the aggregated counts by
  # those factors gives [12, 28/3]. Their total is 64/3, so normalization
  # gives frequencies [9/16, 7/16].
  expect_identical(corrected$motif, c("A_", "T_"))
  expect_identical(corrected$motif_idx, seq_len(2L))
  expect_true(all(is.na(corrected$optional_metadata)))
  expect_equal(corrected$count, c(6, 14))
  expect_equal(corrected$reference_denominator, c(1 / 2, 3 / 2))
  expect_equal(corrected$corrected_count, c(12, 28 / 3))
  expect_equal(corrected$corrected_frequency, c(9 / 16, 7 / 16))
})

test_that("inside core aggregates counts before correction", {
  end_rows <- shared_reference_correction_end_rows()
  output_columns <- names(end_rows)
  end_rows$.cfdnalab_row_order <- 1L
  corrected <- cf_side_reference_corrected_rows(
    end_rows,
    shared_reference_correction_ref_rows(),
    "row_label",
    shared_reference_correction_mode("inside", c("_C", "_G")),
    output_columns,
    "error"
  )
  corrected <- add_shared_corrected_frequencies(corrected)

  # Counts aggregate to _C=2+6=8 and _G=4+8=12. Two positive inside labels
  # make the uniform frequency 1/2. Relative to uniform, frequencies C=3/8
  # and G=5/8 give factors 3/4 and 5/4. Dividing the aggregated counts by
  # those factors gives [32/3, 48/5]. Their total is 304/15, so normalization
  # gives frequencies [10/19, 9/19].
  expect_identical(corrected$motif, c("_C", "_G"))
  expect_identical(corrected$motif_idx, seq_len(2L))
  expect_equal(corrected$count, c(8, 12))
  expect_equal(corrected$reference_denominator, c(3 / 4, 5 / 4))
  expect_equal(corrected$corrected_count, c(32 / 3, 48 / 5))
  expect_equal(corrected$corrected_frequency, c(10 / 19, 9 / 19))
})

test_that("corrected frequencies remain finite when direct totals would overflow", {
  corrected <- data.frame(
    row_label = c("global", "global"),
    corrected_count = c(.Machine$double.xmax, .Machine$double.xmax),
    stringsAsFactors = FALSE
  )

  # Scaling both counts by their row maximum gives [1, 1], whose total is 2.
  # The normalized frequencies are therefore exactly [1/2, 1/2] without
  # summing the original values to infinity.
  corrected <- cf_add_corrected_frequency(corrected, "row_label", "error")

  expect_equal(corrected$corrected_frequency, c(1 / 2, 1 / 2))
  expect_true(all(is.finite(corrected$corrected_frequency)))
})
