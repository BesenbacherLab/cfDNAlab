// Additional cfDNAlab validation for the fixed-size SciPy 1.13.0 optimizer translation
//
// Faithful adaptations of upstream SciPy regressions live in
// `scipy_bfgs_scipy_1_13_tests.rs`. The tests here cover source-level defaults, hand-derived helper
// behavior, the numerical-gradient path used by LIONHEART, and cfDNAlab's explicit error contract.

use super::*;

#[test]
fn scipy_1_13_bfgs_defaults_are_locked_for_three_parameter_model() {
    // LIONHEART passes no method, Jacobian, or optimizer options. These are therefore scientific
    // compatibility settings rather than implementation preferences
    assert_eq!(PARAMETER_COUNT, 3);
    assert_eq!(DEFAULT_MAXIMUM_ITERATIONS, 600);
    assert_eq!(DEFAULT_GRADIENT_TOLERANCE, 1.0e-5);
    assert_eq!(FORWARD_DIFFERENCE_STEP, f64::EPSILON.sqrt());
    assert_eq!(WOLFE_SUFFICIENT_DECREASE, 1.0e-4);
    assert_eq!(WOLFE_CURVATURE, 0.9);
    assert_eq!(
        BFGS_WOLFE_ONE_SETTINGS.relative_step_tolerance,
        1.0e-14
    );
    assert_eq!(BFGS_WOLFE_ONE_SETTINGS.minimum_step, 1.0e-100);
    assert_eq!(BFGS_WOLFE_ONE_SETTINGS.maximum_step, 1.0e100);
    assert_eq!(BFGS_WOLFE_ONE_SETTINGS.maximum_iterations, 100);
    assert_eq!(BFGS_WOLFE_TWO_SETTINGS.maximum_step, Some(1.0e100));
    assert_eq!(BFGS_WOLFE_TWO_SETTINGS.maximum_iterations, 10);
    assert_eq!(BFGS_WOLFE_TWO_SETTINGS.zoom_maximum_iterations, 10);
    assert_eq!(infinity_norm([3.0, -4.0, 12.0]), 12.0);
    assert!(infinity_norm([f64::NAN, 1.0, 2.0]).is_nan());
}

#[test]
fn initial_step_uses_scipy_objective_history_formula() {
    // SciPy calculates min(1, 1.01 * 2 * (9 - 12) / -36) = 0.16833333333333333
    let step = initial_line_search_step(9.0, 12.0, -36.0);

    assert!((step - 0.168_333_333_333_333_33).abs() < 1.0e-15);
}

#[test]
fn initial_step_keeps_scipy_displacement_scaling_for_large_gradients() {
    // SciPy seeds previous_value = value + ||gradient|| / 2. Along the negative gradient,
    // derivative = -||gradient||^2, so the proposed parameter displacement is exactly 1.01 for
    // every gradient norm above 1.01, even when the objective scale is extremely large
    for gradient_norm in [2.0_f64, 2.0e50] {
        let value = 7.0;
        let previous_value = value + gradient_norm / 2.0;
        let directional_derivative = -gradient_norm.powi(2);

        let step = initial_line_search_step(value, previous_value, directional_derivative);
        let displacement = step * gradient_norm;

        assert!((displacement - 1.01).abs() <= 1.01e-14);
    }
}

#[test]
fn quadratic_interpolation_finds_hand_derived_minimum() {
    // f(x) = (x - 2)^2 has f(0) = 4, f'(0) = -4, and f(4) = 4
    let minimum = quadratic_interpolant_minimum(0.0, 4.0, -4.0, 4.0, 4.0);

    assert_eq!(minimum, Some(2.0));
}

#[test]
fn cubic_interpolation_finds_hand_derived_minimum() {
    // f(x) = x^3 - 3x^2 has f(0) = 0, f'(0) = 0, f(1) = -2, and f(3) = 0
    // The stationary points are 0 and 2. SciPy's selected cubic root is the minimum at 2
    let minimum = cubic_interpolant_minimum(0.0, 0.0, 0.0, 1.0, -2.0, 3.0, 0.0);

    assert_eq!(minimum, Some(2.0));
}

#[test]
fn wolfe_one_accepts_step_satisfying_both_strong_wolfe_conditions() -> Result<()> {
    let point = [0.0, 0.0, 0.0];
    let evaluate = |point: [f64; PARAMETER_COUNT]| {
        Ok((point[0] - 3.0).powi(2) + point[1].powi(2) + point[2].powi(2))
    };
    let value = evaluate(point)?;
    let gradient = forward_difference_gradient(point, value, &evaluate)?;
    let direction = gradient.map(|component| -component);
    let derivative = dot(gradient, direction);
    let previous_value = value + euclidean_norm(gradient) / 2.0;

    let accepted =
        line_search_wolfe1(point, direction, gradient, value, previous_value, &evaluate)?
            .expect("the convex quadratic should have a Wolfe-1 step");

    let accepted_derivative = dot(accepted.gradient, direction);
    assert!(accepted.value <= value + WOLFE_SUFFICIENT_DECREASE * accepted.step * derivative);
    assert!(accepted_derivative.abs() <= WOLFE_CURVATURE * derivative.abs());
    Ok(())
}

#[test]
fn bfgs_converges_on_anisotropic_three_parameter_quadratic() -> Result<()> {
    let optimum = [2.0, -3.0, 7.0];
    let fitted = minimize_bfgs([8.0, -0.5, 4.0], &|point| {
        Ok(2.0 * (point[0] - optimum[0]).powi(2)
            + 5.0 * (point[1] - optimum[1]).powi(2)
            + 11.0 * (point[2] - optimum[2]).powi(2))
    })?;

    for (value, target) in fitted.iter().zip(optimum) {
        assert!((value - target).abs() < 1.0e-4);
    }
    Ok(())
}

#[test]
fn lionheart_numerical_gradient_path_recovers_entropy_parameters() -> Result<()> {
    // The dedicated upstream test checks SciPy's objective-value expectation. This stronger
    // cfDNAlab test also checks the fitted parameters because LIONHEART uses this exact default
    // numerical-gradient path. The first feature is constant, so its parameter remains zero.
    let features = [
        [1.0, 1.0, 1.0],
        [1.0, 1.0, 0.0],
        [1.0, 0.0, 1.0],
        [1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
    ];
    let target_moments = [1.0, 0.3, 0.5];
    let objective = |parameters: [f64; PARAMETER_COUNT]| {
        let log_probabilities = features.map(|row| dot(row, parameters));
        let log_normalizer = log_probabilities
            .into_iter()
            .map(f64::exp)
            .sum::<f64>()
            .ln();
        Ok(log_normalizer - dot(target_moments, parameters))
    };
    let known_solution = [0.0, -0.524_869_316, 0.487_525_860];

    let fitted = minimize_bfgs([0.0, 0.0, 0.0], &objective)?;
    let fitted_value = objective(fitted)?;
    let known_value = objective(known_solution)?;

    assert!((fitted_value - known_value).abs() <= 1.0e-6);
    for (parameter, known_parameter) in fitted.into_iter().zip(known_solution) {
        assert!((parameter - known_parameter).abs() <= 1.0e-4);
    }
    Ok(())
}

#[test]
fn bfgs_reports_non_finite_initial_objective() {
    let error = minimize_bfgs([0.0, 0.0, 0.0], &|_point| Ok(f64::NAN))
        .expect_err("a non-finite initial objective must fail fitting");

    assert!(
        error
            .to_string()
            .contains("initial objective is non-finite")
    );
}

#[test]
fn bfgs_reports_non_finite_forward_difference_gradient() {
    // The objective is finite only at the starting point. Each forward-difference probe therefore
    // becomes NaN. Unlike SciPy's unsuccessful NaN result, cfDNAlab reports a fitting error because
    // the command cannot produce a scientifically meaningful model from a non-finite gradient.
    let error = minimize_bfgs([0.0, 0.0, 0.0], &|point| {
        Ok(if point == [0.0, 0.0, 0.0] {
            0.0
        } else {
            f64::NAN
        })
    })
    .expect_err("a non-finite numerical gradient must fail fitting");

    assert!(
        error
            .to_string()
            .contains("initial finite-difference gradient is non-finite")
    );
}

#[test]
fn bfgs_propagates_objective_error_from_forward_difference_probe() {
    let initial_point = [0.0, 0.0, 0.0];
    let error = minimize_bfgs(initial_point, &|point| {
        if point == initial_point {
            Ok(0.0)
        } else {
            bail!("deliberate finite-difference objective failure")
        }
    })
    .expect_err("an objective error from a finite-difference probe must be preserved");

    assert!(
        error
            .to_string()
            .contains("deliberate finite-difference objective failure")
    );
}

#[test]
fn bfgs_propagates_objective_error_from_line_search_probe() {
    let error = minimize_bfgs([0.0, 0.0, 0.0], &|point| {
        if point[0] > 0.5 {
            bail!("deliberate line-search objective failure")
        }
        Ok((point[0] - 1.0).powi(2) + point[1].powi(2) + point[2].powi(2))
    })
    .expect_err("an objective error from a line-search probe must be preserved");

    assert!(
        error
            .to_string()
            .contains("deliberate line-search objective failure")
    );
}

#[test]
fn forward_difference_uses_scaled_representable_step_for_large_parameter() {
    let point = [1.0e20, 0.0, 0.0];
    let evaluate = |candidate: [f64; PARAMETER_COUNT]| Ok(candidate[0] / 1.0e20);
    let value = evaluate(point).expect("linear objective should be finite");

    let gradient = forward_difference_gradient(point, value, &evaluate)
        .expect("linear objective finite differences should succeed");

    assert!(gradient[0].is_finite() && gradient[0] > 0.0);
    // The expected derivative is exactly 1e-20. The tolerance allows rounding in the large
    // representable displacement while remaining eight orders of magnitude below the derivative
    assert!((gradient[0] - 1.0e-20).abs() <= 1.0e-27);
    assert_eq!(gradient[1], 0.0);
    assert_eq!(gradient[2], 0.0);
}

#[test]
fn higher_value_dcstep_brackets_minimum_and_keeps_finite_step() {
    // The current point has a worse value than the best point, which is DCSRCH's first case
    let updated = safeguarded_dcstep(
        0.0, 1.0, -1.0, 0.0, 1.0, -1.0, 1.0, 2.0, 1.0, false, 0.0, 4.0,
    );

    assert!(updated.bracketed);
    assert_eq!(updated.other_step, 1.0);
    assert!(updated.next_step.is_finite());
    assert!(updated.next_step > 0.0 && updated.next_step < 1.0);
}

#[test]
fn opposite_derivative_dcstep_brackets_between_best_and_current_steps() {
    // The current point improves the value but changes the derivative sign, which is DCSRCH's
    // second interpolation case
    let updated = safeguarded_dcstep(
        0.0, 1.0, -1.0, 0.0, 1.0, -1.0, 1.0, 0.5, 1.0, false, 0.0, 4.0,
    );

    assert!(updated.bracketed);
    assert_eq!(updated.best_step, 1.0);
    assert_eq!(updated.other_step, 0.0);
    assert!(updated.next_step > 0.0 && updated.next_step < 1.0);
}

#[test]
fn reduced_derivative_dcstep_extends_unbracketed_search() {
    // The value improves and the negative derivative becomes smaller in magnitude, which is
    // DCSRCH's third interpolation case
    let updated = safeguarded_dcstep(
        0.0, 1.0, -2.0, 0.0, 1.0, -2.0, 1.0, 0.0, -0.5, false, 0.0, 4.0,
    );

    assert!(!updated.bracketed);
    assert_eq!(updated.best_step, 1.0);
    assert!(updated.next_step > 1.0 && updated.next_step <= 4.0);
}

#[test]
fn non_reducing_derivative_dcstep_uses_unbracketed_step_limit() {
    // The value improves but the derivative magnitude increases without a bracket. DCSRCH's fourth
    // case therefore advances to the permitted upper step limit
    let updated = safeguarded_dcstep(
        0.0, 1.0, -1.0, 0.0, 1.0, -1.0, 1.0, 0.0, -2.0, false, 0.0, 4.0,
    );

    assert!(!updated.bracketed);
    assert_eq!(updated.best_step, 1.0);
    assert_eq!(updated.next_step, 4.0);
}
