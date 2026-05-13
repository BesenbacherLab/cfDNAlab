use super::{
    DEFAULT_SAVGOL_WINDOW_BP, MidpointSmoothing, order3_coefficients, parse_midpoint_smoothing,
};

const TOLERANCE: f64 = 1e-8;
const SCIPY_REGRESSION_TOLERANCE: f64 = 1e-12;

#[test]
fn default_smoothing_is_explicit_165bp_savgol() {
    assert_eq!(
        MidpointSmoothing::default(),
        MidpointSmoothing::SavGol {
            window_bp: DEFAULT_SAVGOL_WINDOW_BP
        }
    );
    assert_eq!(MidpointSmoothing::default().to_string(), "savgol=165");
}

#[test]
fn parse_smoothing_accepts_none_and_explicit_savgol_window() {
    assert_eq!(
        parse_midpoint_smoothing("none").expect("none smoothing mode should parse"),
        MidpointSmoothing::None
    );
    assert_eq!(MidpointSmoothing::None.to_string(), "none");
    assert_eq!(
        parse_midpoint_smoothing("savgol=101").expect("explicit odd SavGol window should parse"),
        MidpointSmoothing::SavGol { window_bp: 101 }
    );
}

#[test]
fn parse_smoothing_rejects_even_and_too_short_savgol_windows() {
    let even_error =
        parse_midpoint_smoothing("savgol=100").expect_err("even SavGol windows should fail");
    assert!(
        even_error.contains("must be odd"),
        "unexpected even-window error: {even_error}"
    );

    let short_error =
        parse_midpoint_smoothing("savgol=3").expect_err("order-3 SavGol needs enough points");
    assert!(
        short_error.contains("at least 5 bp"),
        "unexpected short-window error: {short_error}"
    );
}

#[test]
fn parse_smoothing_rejects_unknown_methods_and_malformed_windows() {
    let bare_savgol_error =
        parse_midpoint_smoothing("savgol").expect_err("bare SavGol should not parse");
    assert!(
        bare_savgol_error.contains("Use 'none' or 'savgol=<odd_bp>'"),
        "unexpected bare-savgol error: {bare_savgol_error}"
    );

    let raw_error = parse_midpoint_smoothing("raw").expect_err("raw was renamed to none");
    assert!(
        raw_error.contains("Use 'none' or 'savgol=<odd_bp>'"),
        "unexpected raw-mode error: {raw_error}"
    );

    let unknown_error =
        parse_midpoint_smoothing("gaussian=21").expect_err("unsupported methods should fail");
    assert!(
        unknown_error.contains("unsupported smoothing method"),
        "unexpected unsupported-method error: {unknown_error}"
    );

    let malformed_error =
        parse_midpoint_smoothing("savgol=abc").expect_err("non-integer windows should fail");
    assert!(
        malformed_error.contains("invalid Savitzky-Golay window"),
        "unexpected malformed-window error: {malformed_error}"
    );
}

#[test]
fn order3_coefficients_match_scipy_regression_values() {
    // These constants were generated once from the local SciPy checkout with
    // scipy.signal.savgol_coeffs(window, 3, use="conv") and then hard-coded here
    // This guards against implementation drift. The tests below document the hand derivation
    let scipy_expected = [
        (
            5_u32,
            vec![
                -0.08571428571428572,
                0.34285714285714286,
                0.48571428571428570,
                0.34285714285714286,
                -0.08571428571428572,
            ],
        ),
        (
            9_u32,
            vec![
                -0.09090909090909091,
                0.06060606060606061,
                0.16883116883116883,
                0.23376623376623376,
                0.25541125541125540,
                0.23376623376623376,
                0.16883116883116883,
                0.06060606060606061,
                -0.09090909090909091,
            ],
        ),
    ];

    for (window_size, expected_coefficients) in scipy_expected {
        let coefficients =
            order3_coefficients(window_size).expect("valid regression window should resolve");
        for (observed, expected) in coefficients.iter().zip(expected_coefficients) {
            assert!(
                (observed - expected).abs() <= SCIPY_REGRESSION_TOLERANCE,
                "window {window_size}: expected SciPy coefficient {expected}, got {observed}"
            );
        }
    }
}

#[test]
fn order3_coefficients_match_hand_derived_7bp_window() {
    // For offsets -3..3, the order-3 center smoother coefficients reduce to:
    // [-2, 3, 6, 7, 6, 3, -2] / 21.
    // This is derived from the closed-form coefficient formula in the module docs:
    // m = 3, denominator = 7 * 45 = 315, baseline = 35
    // offset 0: 3 * 35 / 315 = 7 / 21
    // offset 1: 3 * 30 / 315 = 6 / 21
    // offset 2: 3 * 15 / 315 = 3 / 21
    // offset 3: 3 * -10 / 315 = -2 / 21
    let coefficients = order3_coefficients(7).expect("7 bp order-3 coefficients should exist");
    let expected = [-2.0, 3.0, 6.0, 7.0, 6.0, 3.0, -2.0].map(|value| value / 21.0);

    for (observed, expected) in coefficients.iter().zip(expected) {
        assert!(
            (observed - expected).abs() <= TOLERANCE,
            "expected coefficient {expected}, got {observed}"
        );
    }
}

#[test]
fn order3_coefficients_satisfy_moment_constraints_for_165bp_window() {
    let coefficients =
        order3_coefficients(165).expect("165 bp order-3 coefficients should exist");
    let half_window = 82_i32;

    let mut moment0 = 0.0;
    let mut moment1 = 0.0;
    let mut moment2 = 0.0;
    let mut moment3 = 0.0;

    for (coefficient, offset) in coefficients.iter().zip(-half_window..=half_window) {
        let offset = f64::from(offset);
        moment0 += coefficient;
        moment1 += coefficient * offset;
        moment2 += coefficient * offset.powi(2);
        moment3 += coefficient * offset.powi(3);
    }

    assert!(
        (moment0 - 1.0).abs() <= TOLERANCE,
        "constant moment should be 1, got {moment0}"
    );
    assert!(
        moment1.abs() <= TOLERANCE,
        "linear moment should be 0, got {moment1}"
    );
    assert!(
        moment2.abs() <= TOLERANCE,
        "quadratic moment should be 0, got {moment2}"
    );
    assert!(
        moment3.abs() <= TOLERANCE,
        "cubic moment should be 0, got {moment3}"
    );
}

#[test]
fn order3_coefficients_are_symmetric() {
    let coefficients = order3_coefficients(21).expect("21 bp coefficients should exist");

    for coefficient_idx in 0..coefficients.len() {
        let mirror_idx = coefficients.len() - 1 - coefficient_idx;
        assert!(
            (coefficients[coefficient_idx] - coefficients[mirror_idx]).abs() <= TOLERANCE,
            "coefficients should be symmetric at {coefficient_idx} and {mirror_idx}"
        );
    }
}

#[test]
fn order3_coefficients_preserve_polynomials_at_center() {
    let coefficients = order3_coefficients(21).expect("21 bp coefficients should exist");
    let half_window = 10_i32;

    let mut constant = 0.0;
    let mut linear = 0.0;
    let mut quadratic = 0.0;
    let mut cubic = 0.0;
    for (coefficient, offset) in coefficients.iter().zip(-half_window..=half_window) {
        let offset = f64::from(offset);
        constant += coefficient * 5.0;
        linear += coefficient * (2.0 * offset + 5.0);
        quadratic += coefficient * (offset.powi(2) + 2.0 * offset + 5.0);
        cubic += coefficient * (offset.powi(3) + offset.powi(2) + 2.0 * offset + 5.0);
    }

    // At the center offset 0, all four polynomials above evaluate to 5
    for (name, value) in [
        ("constant", constant),
        ("linear", linear),
        ("quadratic", quadratic),
        ("cubic", cubic),
    ] {
        assert!(
            (value - 5.0).abs() <= TOLERANCE,
            "{name} polynomial should be preserved at the center, got {value}"
        );
    }
}
