use super::{
    ClassifiedGCTagWeight, GCTagValue, SanitizedGCWeight, combine_gc_tag_values,
    sanitize_gc_weight,
};
use crate::shared::base::ZEROISH_F32_TOLERANCE;

#[test]
fn sanitize_gc_weight_snaps_small_negative_values_to_zero() {
    // Arrange
    // A small negative residue should be treated the same as a small positive one
    let small_negative = -(ZEROISH_F32_TOLERANCE as f64) / 2.0;

    // Act
    let sanitized = sanitize_gc_weight(small_negative);

    // Assert
    assert_eq!(sanitized, SanitizedGCWeight::Usable(0.0));
}

#[test]
fn sanitize_gc_weight_marks_large_negative_values_as_out_of_range() {
    // Arrange
    // A meaningfully negative GC weight is not recoverable numeric noise
    let large_negative = -3.0_f64;

    // Act
    let sanitized = sanitize_gc_weight(large_negative);

    // Assert
    assert_eq!(
        sanitized,
        SanitizedGCWeight::Unusable {
            out_of_range: true
        }
    );
}

#[test]
fn gc_tag_value_from_number_marks_large_negative_values_as_invalid_and_out_of_range() {
    // Arrange
    // Record parsing should preserve the same invalid/out-of-range classification
    // as direct sanitizer use
    let large_negative = -3.0_f32;

    // Act
    let parsed = GCTagValue::from_number(large_negative);

    // Assert
    assert_eq!(
        parsed.classify().expect("parsed tag should classify"),
        ClassifiedGCTagWeight::Invalid {
            out_of_range: true
        }
    );
}

#[test]
fn combine_gc_tag_values_reuses_single_usable_weight_when_other_mate_is_missing() {
    // Arrange
    // A single usable mate should still provide the fragment weight when the other mate
    // is simply missing rather than invalid
    let usable_tag = GCTagValue {
        weight: Some(1.75),
        was_missing: false,
        had_invalid: false,
        was_out_of_range: false,
    };
    let missing_tag = GCTagValue::missing();

    // Act
    let combined_left = combine_gc_tag_values(&usable_tag, &missing_tag);
    let combined_right = combine_gc_tag_values(&missing_tag, &usable_tag);

    // Assert
    assert_eq!(
        combined_left.classify().expect("combined left tag should classify"),
        ClassifiedGCTagWeight::Usable(1.75)
    );
    assert_eq!(
        combined_right
            .classify()
            .expect("combined right tag should classify"),
        ClassifiedGCTagWeight::Usable(1.75)
    );
}

#[test]
fn combine_gc_tag_values_keeps_missing_when_both_mates_are_missing() {
    // Arrange
    // The fragment is only truly missing when neither mate contributes a usable weight
    let left_missing_tag = GCTagValue::missing();
    let right_missing_tag = GCTagValue::missing();

    // Act
    let combined = combine_gc_tag_values(&left_missing_tag, &right_missing_tag);

    // Assert
    assert_eq!(
        combined.classify().expect("combined tag should classify"),
        ClassifiedGCTagWeight::Missing
    );
}

#[test]
fn combine_gc_tag_values_keeps_zero_weight_even_when_other_mate_is_missing() {
    // Arrange
    // An explicit zero keeps top priority, even if the mate is missing
    let zero_tag = GCTagValue {
        weight: Some(0.0),
        was_missing: false,
        had_invalid: false,
        was_out_of_range: false,
    };
    let missing_tag = GCTagValue::missing();

    // Act
    let combined = combine_gc_tag_values(&zero_tag, &missing_tag);

    // Assert
    assert_eq!(
        combined.classify().expect("combined tag should classify"),
        ClassifiedGCTagWeight::Usable(0.0)
    );
}

#[test]
fn combine_gc_tag_values_treats_positive_zero_snap_threshold_as_zero_precedence() {
    // Arrange
    // The inclusive positive snap boundary should collapse to explicit zero
    // before any averaging happens
    let zeroish_tag = GCTagValue {
        weight: Some(ZEROISH_F32_TOLERANCE),
        was_missing: false,
        had_invalid: false,
        was_out_of_range: false,
    };
    let usable_tag = GCTagValue {
        weight: Some(2.0),
        was_missing: false,
        had_invalid: false,
        was_out_of_range: false,
    };

    // Act
    let combined = combine_gc_tag_values(&zeroish_tag, &usable_tag);

    // Assert
    assert_eq!(
        combined.classify().expect("combined tag should classify"),
        ClassifiedGCTagWeight::Usable(0.0)
    );
}

#[test]
fn combine_gc_tag_values_treats_negative_zero_snap_threshold_as_zero_precedence() {
    // Arrange
    // The snap window is symmetric around zero. More negative values are covered
    // by the invalid/out-of-range tests above.
    let zeroish_tag = GCTagValue {
        weight: Some(-ZEROISH_F32_TOLERANCE),
        was_missing: false,
        had_invalid: false,
        was_out_of_range: false,
    };
    let usable_tag = GCTagValue {
        weight: Some(2.0),
        was_missing: false,
        had_invalid: false,
        was_out_of_range: false,
    };

    // Act
    let combined = combine_gc_tag_values(&zeroish_tag, &usable_tag);

    // Assert
    assert_eq!(
        combined.classify().expect("combined tag should classify"),
        ClassifiedGCTagWeight::Usable(0.0)
    );
}

#[test]
fn combine_gc_tag_values_keeps_invalid_above_single_usable_or_missing_mate() {
    // Arrange
    // Invalid metadata should poison the fragment instead of being hidden by
    // a usable mate or a missing mate
    let usable_tag = GCTagValue {
        weight: Some(2.0),
        was_missing: false,
        had_invalid: false,
        was_out_of_range: false,
    };
    let invalid_tag = GCTagValue {
        weight: None,
        was_missing: false,
        had_invalid: true,
        was_out_of_range: false,
    };
    let missing_tag = GCTagValue::missing();

    // Act
    let combined_with_usable = combine_gc_tag_values(&usable_tag, &invalid_tag);
    let combined_with_missing = combine_gc_tag_values(&missing_tag, &invalid_tag);

    // Assert
    assert_eq!(
        combined_with_usable
            .classify()
            .expect("combined usable+invalid tag should classify"),
        ClassifiedGCTagWeight::Invalid {
            out_of_range: false
        }
    );
    assert_eq!(
        combined_with_missing
            .classify()
            .expect("combined missing+invalid tag should classify"),
        ClassifiedGCTagWeight::Invalid {
            out_of_range: false
        }
    );
}
