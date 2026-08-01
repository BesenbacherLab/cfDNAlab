use super::*;

fn interval(start: u32, end: u32) -> Interval<u32> {
    Interval::new(start, end).expect("test interval should be valid")
}

fn fragment(
    start: u32,
    end: u32,
    segments: Option<&[(u32, u32)]>,
) -> FragmentWithSegments {
    FragmentWithSegments {
        tid: 0,
        interval: interval(start, end),
        segments: segments.map(|entries| {
            entries
                .iter()
                .map(|(segment_start, segment_end)| interval(*segment_start, *segment_end))
                .collect()
        }),
        gc_tag: Default::default(),
        #[cfg(feature = "cmd_overlapping_lengths_correction")]
        overlap_length_weight: 1.0,
    }
}

fn segment_tuples(segments: &[Interval<u32>]) -> Vec<(u32, u32)> {
    segments.iter().map(Interval::as_tuple).collect()
}

#[test]
fn fragment_span_trim_roundtrips_supported_modes() {
    // Arrange
    let values = ["at-most=165", "exactly=165"];

    for value in values {
        // Act
        let parsed = value
            .parse::<FragmentSpanTrim>()
            .expect("supported trim mode should parse");

        // Assert
        assert_eq!(parsed.to_string(), value);
        assert_eq!(
            parsed
                .to_string()
                .parse::<FragmentSpanTrim>()
                .expect("displayed trim mode should parse"),
            parsed
        );
    }
}

#[test]
fn fragment_span_trim_accepts_minimum_target_and_canonicalizes_mode() {
    // Arrange
    let value = "AT-MOST=1";

    // Act
    let parsed = value
        .parse::<FragmentSpanTrim>()
        .expect("minimum supported target should parse");

    // Assert
    assert_eq!(parsed.target_length(), 1);
    assert_eq!(parsed.to_string(), "at-most=1");
}

#[test]
fn fragment_span_trim_rejects_even_target() {
    // Act
    let error = "at-most=166"
        .parse::<FragmentSpanTrim>()
        .expect_err("even target should fail");

    // Assert
    assert!(error.contains("must be odd"), "unexpected error: {error}");
}

#[test]
fn fragment_span_trim_rejects_target_below_minimum_trim_length() {
    // Act
    let error = "exactly=0"
        .parse::<FragmentSpanTrim>()
        .expect_err("zero-length target should fail");

    // Assert
    assert!(
        error.contains("at least 1 bp"),
        "unexpected error: {error}"
    );
}

#[test]
fn fragment_span_trim_rejects_target_above_supported_fragment_length() {
    // Act
    let error = "at-most=50001"
        .parse::<FragmentSpanTrim>()
        .expect_err("target above supported fragment length should fail");

    // Assert
    assert!(error.contains("<= 50000 bp"), "unexpected error: {error}");
}

#[test]
fn fragment_span_trim_rejects_unknown_or_incomplete_modes() {
    // Act
    let unknown_error = "resize=165"
        .parse::<FragmentSpanTrim>()
        .expect_err("unknown mode should fail");
    let incomplete_error = "at-most"
        .parse::<FragmentSpanTrim>()
        .expect_err("missing target should fail");

    // Assert
    assert!(
        unknown_error.contains("unsupported trim mode"),
        "unexpected error: {unknown_error}"
    );
    assert!(
        incomplete_error.contains("invalid trim specification"),
        "unexpected error: {incomplete_error}"
    );
}

#[test]
fn fragment_span_trim_rejects_non_integer_and_overflowing_targets() {
    // Act
    let non_integer_error = "at-most=61.5"
        .parse::<FragmentSpanTrim>()
        .expect_err("non-integer target should fail");
    let overflowing_error = "exactly=4294967296"
        .parse::<FragmentSpanTrim>()
        .expect_err("target above u32 should fail");

    // Assert
    assert!(
        non_integer_error.contains("invalid trim target"),
        "unexpected error: {non_integer_error}"
    );
    assert!(
        overflowing_error.contains("invalid trim target"),
        "unexpected error: {overflowing_error}"
    );
}

#[test]
fn at_most_centers_odd_target_on_existing_even_fragment_midpoint() {
    // Arrange
    let original = fragment(100, 300, None);
    let target_length = 165;
    let midpoint = midpoint_random_even_for_fragment("chr1", 100, 200);
    let expected = interval(midpoint - 82, midpoint + 83);

    // Act
    let counting_segments = trimmed_counting_segments(
        &original,
        "chr1",
        1_000,
        FragmentSpanTrim::AtMost { target_length },
    )
    .expect("trim should succeed");

    // Assert
    assert_eq!(counting_segments, vec![expected]);
    assert_eq!(counting_segments[0].len(), target_length);
}

#[test]
fn at_most_centers_target_on_existing_odd_fragment_midpoint() {
    // Arrange
    // The unique midpoint of [100, 301) is base 200. A centered 165 bp span therefore starts
    // 82 bases before it and ends 83 bases after it.
    let original = fragment(100, 301, None);

    // Act
    let counting_segments = trimmed_counting_segments(
        &original,
        "chr1",
        1_000,
        FragmentSpanTrim::AtMost { target_length: 165 },
    )
    .expect("trim should succeed");

    // Assert
    assert_eq!(counting_segments, vec![interval(118, 283)]);
}

#[test]
fn at_most_one_base_counts_only_the_midpoint_base() {
    // Arrange
    let original = fragment(100, 301, None);

    // Act
    let counting_segments = trimmed_counting_segments(
        &original,
        "chr1",
        1_000,
        FragmentSpanTrim::AtMost { target_length: 1 },
    )
    .expect("one-base trim should succeed");

    // Assert
    assert_eq!(counting_segments, vec![interval(200, 201)]);
}

#[test]
fn at_most_leaves_shorter_fragment_unchanged() {
    // Arrange
    let original = fragment(100, 200, None);

    // Act
    let counting_segments = trimmed_counting_segments(
        &original,
        "chr1",
        1_000,
        FragmentSpanTrim::AtMost { target_length: 165 },
    )
    .expect("at-most projection should succeed");

    // Assert
    assert_eq!(counting_segments, vec![original.interval]);
    assert_eq!(counting_segments[0].len(), 100);
}

#[test]
fn at_most_leaves_fragment_at_requested_span_unchanged() {
    // Arrange
    let original = fragment(100, 265, None);

    // Act
    let counting_segments = trimmed_counting_segments(
        &original,
        "chr1",
        1_000,
        FragmentSpanTrim::AtMost { target_length: 165 },
    )
    .expect("equal-span projection should succeed");

    // Assert
    assert_eq!(counting_segments, vec![original.interval]);
}

#[test]
fn exactly_leaves_fragment_at_requested_span_unchanged() {
    // Arrange
    let original = fragment(100, 265, None);

    // Act
    let counting_segments = trimmed_counting_segments(
        &original,
        "chr1",
        1_000,
        FragmentSpanTrim::Exactly { target_length: 165 },
    )
    .expect("equal-span projection should succeed");

    // Assert
    assert_eq!(counting_segments, vec![original.interval]);
}

#[test]
fn exactly_trims_longer_fragment_to_requested_span() {
    // Arrange
    let original = fragment(100, 301, None);

    // Act
    let counting_segments = trimmed_counting_segments(
        &original,
        "chr1",
        1_000,
        FragmentSpanTrim::Exactly { target_length: 165 },
    )
    .expect("exact projection should succeed");

    // Assert
    assert_eq!(counting_segments, vec![interval(118, 283)]);
}

#[test]
fn exactly_extends_shorter_fragment_to_requested_span() {
    // Arrange
    let original = fragment(100, 201, None);

    // Act
    let counting_segments = trimmed_counting_segments(
        &original,
        "chr1",
        1_000,
        FragmentSpanTrim::Exactly { target_length: 165 },
    )
    .expect("exact projection should succeed");

    // Assert
    assert_eq!(counting_segments, vec![interval(68, 233)]);
    assert_eq!(counting_segments[0].len(), 165);
}

#[test]
fn exactly_clips_extension_at_chromosome_start_without_shifting() {
    // Arrange
    let original = fragment(0, 101, None);

    // Act
    let counting_segments = trimmed_counting_segments(
        &original,
        "chr1",
        1_000,
        FragmentSpanTrim::Exactly { target_length: 165 },
    )
    .expect("boundary projection should succeed");

    // Assert
    // The selected midpoint is 50. The intended span is [-32, 133), so clipping retains
    // [0, 133) instead of shifting the span right to recover 165 bp.
    assert_eq!(counting_segments, vec![interval(0, 133)]);
}

#[test]
fn exactly_clips_extension_at_chromosome_end_without_shifting() {
    // Arrange
    let original = fragment(899, 1_000, None);

    // Act
    let counting_segments = trimmed_counting_segments(
        &original,
        "chr1",
        1_000,
        FragmentSpanTrim::Exactly { target_length: 165 },
    )
    .expect("boundary projection should succeed");

    // Assert
    // The midpoint is 949. The intended span is [867, 1032), so clipping retains [867, 1000)
    // instead of shifting the span left to recover 165 bp.
    assert_eq!(counting_segments, vec![interval(867, 1_000)]);
}

#[test]
fn at_most_preserves_segments_when_fragment_is_shorter_than_target() {
    // Arrange
    let original = fragment(100, 201, Some(&[(100, 120), (180, 201)]));

    // Act
    let counting_segments = trimmed_counting_segments(
        &original,
        "chr1",
        1_000,
        FragmentSpanTrim::AtMost { target_length: 165 },
    )
    .expect("short segmented fragment should remain unchanged");

    // Assert
    assert_eq!(
        segment_tuples(&counting_segments),
        vec![(100, 120), (180, 201)]
    );
}

#[test]
fn trimming_intersects_existing_segments_with_centered_span() {
    // Arrange
    let original = fragment(
        100,
        301,
        Some(&[(100, 140), (180, 220), (260, 301)]),
    );

    // Act
    let counting_segments = trimmed_counting_segments(
        &original,
        "chr1",
        1_000,
        FragmentSpanTrim::AtMost { target_length: 165 },
    )
    .expect("segmented trim should succeed");

    // Assert
    assert_eq!(
        segment_tuples(&counting_segments),
        vec![(118, 140), (180, 220), (260, 283)]
    );
}

#[test]
fn exactly_adds_outer_flanks_without_filling_internal_segment_gap() {
    // Arrange
    let original = fragment(100, 201, Some(&[(100, 120), (180, 201)]));

    // Act
    let counting_segments = trimmed_counting_segments(
        &original,
        "chr1",
        1_000,
        FragmentSpanTrim::Exactly { target_length: 165 },
    )
    .expect("segmented extension should succeed");

    // Assert
    assert_eq!(
        segment_tuples(&counting_segments),
        vec![(68, 120), (180, 233)]
    );
    assert_eq!(
        counting_segments
            .iter()
            .map(Interval::len)
            .sum::<u32>(),
        105,
        "the 60 bp internal gap must remain excluded from the 165 bp outer span"
    );
}

#[test]
fn trimming_to_internal_gap_returns_no_counting_segments() {
    // Arrange
    // The selected midpoint is 200, and the 1 bp trimmed span is [200, 201). Both retained
    // segments lie outside that span, so restoring any bases would violate segment semantics.
    let original = fragment(100, 301, Some(&[(100, 150), (250, 301)]));

    // Act
    let counting_segments = trimmed_counting_segments(
        &original,
        "chr1",
        1_000,
        FragmentSpanTrim::AtMost { target_length: 1 },
    )
    .expect("segmented trim should succeed");

    // Assert
    assert!(counting_segments.is_empty());
}
