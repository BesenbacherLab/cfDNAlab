use super::*;

#[test]
fn gc_length_trim_rare_validation_accepts_fractional_values_below_one() {
    // Arrange / Act / Assert
    validate_gc_length_trim_rare(0.0).expect("zero trim should preserve existing behavior");
    validate_gc_length_trim_rare(0.05).expect("small trim fraction should be accepted");
    validate_gc_length_trim_rare(0.999_999).expect("values below one should be accepted");
}

#[test]
fn gc_length_trim_rare_validation_rejects_invalid_values() {
    // Arrange
    let invalid_values = [-0.1, 1.0, 1.1, f64::NAN, f64::INFINITY, f64::NEG_INFINITY];

    for invalid_value in invalid_values {
        // Act
        let error = validate_gc_length_trim_rare(invalid_value)
            .expect_err("invalid trim fraction should fail");

        // Assert
        assert!(
            error
                .to_string()
                .contains("--gc-length-trim-rare must be finite and within [0, 1)"),
            "unexpected error for {invalid_value}: {error}"
        );
    }
}
