// Faithful Rust adaptations of the SciPy 1.13.0 regressions that exercise the optimizer path
// translated in this module
//
// The source tests are `scipy/optimize/tests/test_linesearch.py` and
// `scipy/optimize/tests/test_optimize.py` at tag v1.13.0. Objective functions, derivatives,
// seeded historical objective values, algorithm settings, and behavioral assertions are retained.
// Rust-specific assertions only adapt SciPy return tuples to the smaller internal API used here.

use super::*;
use std::cell::Cell;

/// Defaults of SciPy's standalone `scalar_search_wolfe1` API.
///
/// BFGS overrides the minimum and maximum with `1e-100` and `1e100`. The upstream scalar tests do
/// not, so they must use the public scalar defaults instead of the production BFGS settings.
const SCIPY_SCALAR_WOLFE_ONE_SETTINGS: WolfeOneSettings = WolfeOneSettings {
    sufficient_decrease: WOLFE_SUFFICIENT_DECREASE,
    curvature: WOLFE_CURVATURE,
    relative_step_tolerance: WOLFE_ONE_RELATIVE_STEP_TOLERANCE,
    minimum_step: 1.0e-8,
    maximum_step: 50.0,
    maximum_iterations: WOLFE_ONE_MAXIMUM_ITERATIONS,
};

/// Defaults of SciPy's standalone `scalar_search_wolfe2` API.
const SCIPY_SCALAR_WOLFE_TWO_SETTINGS: WolfeTwoSettings = WolfeTwoSettings {
    sufficient_decrease: WOLFE_SUFFICIENT_DECREASE,
    curvature: WOLFE_CURVATURE,
    maximum_step: None,
    maximum_iterations: WOLFE_TWO_MAXIMUM_ITERATIONS,
    zoom_maximum_iterations: WOLFE_TWO_ZOOM_MAXIMUM_ITERATIONS,
};

#[derive(Clone, Copy)]
struct ScalarLineSearchCase {
    name: &'static str,
    objective: fn(f64) -> f64,
    derivative: fn(f64) -> f64,
    historical_values: [f64; 3],
}

/// SciPy's `_scalar_func_1`: a polynomial with competing cubic and quartic terms.
fn scalar_objective_polynomial(step: f64) -> f64 {
    -step - step.powi(3) + step.powi(4)
}

fn scalar_derivative_polynomial(step: f64) -> f64 {
    -1.0 - 3.0 * step.powi(2) + 4.0 * step.powi(3)
}

/// SciPy's `_scalar_func_2`: an exponential-to-quadratic transition.
fn scalar_objective_exponential(step: f64) -> f64 {
    (-4.0 * step).exp() + step.powi(2)
}

fn scalar_derivative_exponential(step: f64) -> f64 {
    -4.0 * (-4.0 * step).exp() + 2.0 * step
}

/// SciPy's `_scalar_func_3`: a non-convex line with several local extrema.
fn scalar_objective_oscillating(step: f64) -> f64 {
    -(10.0 * step).sin()
}

fn scalar_derivative_oscillating(step: f64) -> f64 {
    -10.0 * (10.0 * step).cos()
}

/// Return the seeded inputs used by SciPy after its test setup creates a 20 by 20 random matrix.
///
/// SciPy calls `np.random.seed(1234)`, consumes 400 values for `self.A`, and then draws three values
/// for each scalar function. The values are test inputs, not expected optimizer outputs. Making them
/// explicit avoids adding a NumPy-compatible random-number generator for a fixed upstream fixture.
fn scalar_line_search_cases() -> [ScalarLineSearchCase; 3] {
    [
        ScalarLineSearchCase {
            name: "polynomial",
            objective: scalar_objective_polynomial,
            derivative: scalar_derivative_polynomial,
            historical_values: [
                -0.226_632_293_676_835_93,
                -0.923_830_720_176_923,
                0.355_838_997_472_820_33,
            ],
        },
        ScalarLineSearchCase {
            name: "exponential",
            objective: scalar_objective_exponential,
            derivative: scalar_derivative_exponential,
            historical_values: [
                -1.270_063_478_386_288_5,
                -0.195_472_412_696_758_52,
                -0.463_419_399_217_462_95,
            ],
        },
        ScalarLineSearchCase {
            name: "oscillating",
            objective: scalar_objective_oscillating,
            derivative: scalar_derivative_oscillating,
            historical_values: [
                0.989_414_847_285_850_7,
                1.388_647_381_631_717,
                1.087_713_736_524_141_2,
            ],
        },
    ]
}

/// Reproduce SciPy's 50-ULP comparison for finite scalar line-search results.
fn assert_float_within_50_ulps(actual: f64, expected: f64, context: &str) {
    assert!(
        actual.is_finite() && expected.is_finite(),
        "{context} compared non-finite values: actual={actual}, expected={expected}"
    );
    if actual == expected {
        return;
    }
    let ordered_bits = |value: f64| {
        let bits = value.to_bits();
        if bits & (1_u64 << 63) == 0 {
            bits | (1_u64 << 63)
        } else {
            !bits
        }
    };
    let ulp_distance = ordered_bits(actual).abs_diff(ordered_bits(expected));
    assert!(
        ulp_distance <= 50,
        "{context} differed by {ulp_distance} ULPs: actual={actual}, expected={expected}"
    );
}

/// Reproduce NumPy's default `assert_allclose` scalar tolerance with an explicit absolute tolerance.
fn assert_allclose_with_absolute_tolerance(actual: f64, expected: f64, absolute_tolerance: f64) {
    const DEFAULT_RELATIVE_TOLERANCE: f64 = 1.0e-7;
    let allowed_difference =
        absolute_tolerance + DEFAULT_RELATIVE_TOLERANCE * expected.abs();
    assert!(
        (actual - expected).abs() <= allowed_difference,
        "actual={actual}, expected={expected}, allowed difference={allowed_difference}"
    );
}

/// Apply the same strong-Wolfe assertions as SciPy's `assert_wolfe` helper.
fn assert_strong_wolfe_conditions(
    step: f64,
    objective: &impl Fn(f64) -> f64,
    derivative: &impl Fn(f64) -> f64,
    sufficient_decrease: f64,
    curvature: f64,
    case_name: &str,
) {
    let initial_value = objective(0.0);
    let initial_derivative = derivative(0.0);
    let accepted_value = objective(step);
    let accepted_derivative = derivative(step);

    assert!(
        accepted_value <= initial_value + sufficient_decrease * step * initial_derivative,
        "{case_name} failed sufficient decrease at step {step}"
    );
    assert!(
        accepted_derivative.abs() <= curvature * initial_derivative.abs(),
        "{case_name} failed the curvature condition at step {step}"
    );
}

/// Drive the translated DCSRCH state with SciPy's analytic scalar derivative.
fn scalar_search_wolfe1(
    objective: &impl Fn(f64) -> f64,
    derivative: &impl Fn(f64) -> f64,
    initial_value: f64,
    historical_value: Option<f64>,
    initial_derivative: f64,
) -> Option<LineSearchPoint> {
    let mut step = historical_value
        .map(|previous| initial_line_search_step(initial_value, previous, initial_derivative))
        .unwrap_or(1.0);
    let mut candidate = LineSearchPoint {
        step,
        value: initial_value,
        gradient: [initial_derivative, 0.0, 0.0],
    };
    let mut state = DcsrchState::new(
        initial_value,
        initial_derivative,
        step,
        SCIPY_SCALAR_WOLFE_ONE_SETTINGS,
    )?;

    for iteration in 0..SCIPY_SCALAR_WOLFE_ONE_SETTINGS.maximum_iterations {
        let candidate_derivative = if iteration == 0 {
            initial_derivative
        } else {
            candidate.gradient[0]
        };
        match state.iterate(
            step,
            candidate.value,
            candidate_derivative,
            iteration == 0,
        ) {
            DcsrchAction::Evaluate(next_step) => {
                step = next_step;
                candidate = LineSearchPoint {
                    step,
                    value: objective(step),
                    gradient: [derivative(step), 0.0, 0.0],
                };
            }
            DcsrchAction::Converged => return Some(candidate),
            DcsrchAction::Failed => return None,
        }
    }
    None
}

/// Drive the translated Wolfe-2 core with SciPy's analytic scalar derivative.
fn scalar_search_wolfe2_with_analytic_derivative(
    objective: &impl Fn(f64) -> f64,
    derivative: &impl Fn(f64) -> f64,
    initial_value: f64,
    historical_value: Option<f64>,
    initial_derivative: f64,
    settings: WolfeTwoSettings,
) -> Option<WolfeTwoPoint> {
    let initial_step = historical_value
        .map(|previous| initial_line_search_step(initial_value, previous, initial_derivative))
        .unwrap_or(1.0);
    let evaluate_value = |step| Ok(objective(step));
    let evaluate_gradient = |step, _value| Ok([derivative(step), 0.0, 0.0]);

    scalar_search_wolfe2(
        [1.0, 0.0, 0.0],
        initial_value,
        initial_derivative,
        initial_step,
        settings,
        &evaluate_value,
        &evaluate_gradient,
    )
    .expect("analytic scalar objective evaluation cannot fail")
}

/// Adapt SciPy's vector Wolfe-2 wrapper to the fixed three-parameter implementation.
///
/// The upstream boundary regression is two-dimensional. A trailing zero coordinate preserves every
/// objective value, gradient, directional derivative, and accepted step while exercising the vector
/// wrapper logic that the scalar regression does not cover.
fn line_search_wolfe2_with_analytic_gradient(
    objective: &impl Fn([f64; PARAMETER_COUNT]) -> f64,
    gradient: &impl Fn([f64; PARAMETER_COUNT]) -> [f64; PARAMETER_COUNT],
    point: [f64; PARAMETER_COUNT],
    direction: [f64; PARAMETER_COUNT],
    settings: WolfeTwoSettings,
) -> Option<WolfeTwoPoint> {
    let initial_value = objective(point);
    let initial_derivative = dot(gradient(point), direction);
    let evaluate_value = |step| Ok(objective(add_scaled(point, direction, step)));
    let evaluate_gradient = |step, _value| Ok(gradient(add_scaled(point, direction, step)));

    scalar_search_wolfe2(
        direction,
        initial_value,
        initial_derivative,
        1.0,
        settings,
        &evaluate_value,
        &evaluate_gradient,
    )
    .expect("analytic vector objective evaluation cannot fail")
}

#[test]
fn scipy_test_scalar_search_wolfe1() {
    let mut case_count = 0_usize;
    for case in scalar_line_search_cases() {
        for historical_value in case.historical_values {
            case_count += 1;
            let initial_value = (case.objective)(0.0);
            let initial_derivative = (case.derivative)(0.0);

            let accepted = scalar_search_wolfe1(
                &case.objective,
                &case.derivative,
                initial_value,
                Some(historical_value),
                initial_derivative,
            )
            .unwrap_or_else(|| {
                panic!(
                    "Wolfe-1 did not converge for {} with historical value {}",
                    case.name, historical_value
                )
            });

            assert_float_within_50_ulps(
                accepted.value,
                (case.objective)(accepted.step),
                case.name,
            );
            assert_strong_wolfe_conditions(
                accepted.step,
                &case.objective,
                &case.derivative,
                WOLFE_SUFFICIENT_DECREASE,
                WOLFE_CURVATURE,
                case.name,
            );
        }
    }
    // Match SciPy's guard that the seeded iterator exercised more than a single fixture
    assert!(case_count > 3);
}

#[test]
fn scipy_test_scalar_search_wolfe2() {
    for case in scalar_line_search_cases() {
        for historical_value in case.historical_values {
            let initial_value = (case.objective)(0.0);
            let initial_derivative = (case.derivative)(0.0);

            let accepted = scalar_search_wolfe2_with_analytic_derivative(
                &case.objective,
                &case.derivative,
                initial_value,
                Some(historical_value),
                initial_derivative,
                SCIPY_SCALAR_WOLFE_TWO_SETTINGS,
            )
            .unwrap_or_else(|| {
                panic!(
                    "Wolfe-2 did not converge for {} with historical value {}",
                    case.name, historical_value
                )
            });

            assert_float_within_50_ulps(
                accepted.value,
                (case.objective)(accepted.step),
                case.name,
            );
            if let Some(accepted_gradient) = accepted.gradient {
                assert_float_within_50_ulps(
                    accepted_gradient[0],
                    (case.derivative)(accepted.step),
                    case.name,
                );
            }
            assert_strong_wolfe_conditions(
                accepted.step,
                &case.objective,
                &case.derivative,
                WOLFE_SUFFICIENT_DECREASE,
                WOLFE_CURVATURE,
                case.name,
            );
        }
    }
}

#[test]
fn scipy_test_scalar_search_wolfe2_with_low_amax() {
    let objective = |step: f64| (step - 5.0).powi(2);
    let derivative = |step: f64| 2.0 * (step - 5.0);
    let settings = WolfeTwoSettings {
        maximum_step: Some(0.001),
        ..SCIPY_SCALAR_WOLFE_TWO_SETTINGS
    };

    let accepted = scalar_search_wolfe2_with_analytic_derivative(
        &objective,
        &derivative,
        objective(0.0),
        None,
        derivative(0.0),
        settings,
    );

    assert!(accepted.is_none());
}

#[test]
fn scipy_test_scalar_search_wolfe2_wrong_basin_regression() {
    // SciPy gh-12157 and gh-13073: the minimum is at 4 / 3, while the old implementation returned
    // step 2.0 in a different basin.
    let objective = |step: f64| {
        if step < 1.0 {
            -1.5 * std::f64::consts::PI * (step - 1.0)
        } else {
            (1.5 * std::f64::consts::PI * step - std::f64::consts::PI).cos()
        }
    };
    let derivative = |step: f64| {
        if step < 1.0 {
            -1.5 * std::f64::consts::PI
        } else {
            -1.5 * std::f64::consts::PI
                * (1.5 * std::f64::consts::PI * step - std::f64::consts::PI).sin()
        }
    };

    let accepted = scalar_search_wolfe2_with_analytic_derivative(
        &objective,
        &derivative,
        objective(0.0),
        None,
        derivative(0.0),
        SCIPY_SCALAR_WOLFE_TWO_SETTINGS,
    )
    .expect("SciPy's wrong-basin regression objective should converge");

    assert!(accepted.step < 1.5, "accepted step was {}", accepted.step);
}

#[test]
fn scipy_test_line_search_wolfe2_bounds() {
    // This is SciPy gh-7475 with one trailing zero coordinate for the fixed three-parameter port
    let objective = |point: [f64; PARAMETER_COUNT]| dot(point, point);
    let gradient = |point: [f64; PARAMETER_COUNT]| point.map(|component| 2.0 * component);
    let point = [-60.0, 0.0, 0.0];
    let direction = [1.0, 0.0, 0.0];
    let scalar_objective = |step: f64| objective(add_scaled(point, direction, step));
    let scalar_derivative = |step: f64| dot(gradient(add_scaled(point, direction, step)), direction);
    let settings_at_boundary = WolfeTwoSettings {
        curvature: 0.5,
        maximum_step: Some(30.0),
        ..SCIPY_SCALAR_WOLFE_TWO_SETTINGS
    };
    let accepted = line_search_wolfe2_with_analytic_gradient(
        &objective,
        &gradient,
        point,
        direction,
        settings_at_boundary,
    )
    .expect("Wolfe-2 should accept the exact maximum-step boundary");

    assert_strong_wolfe_conditions(
        accepted.step,
        &scalar_objective,
        &scalar_derivative,
        WOLFE_SUFFICIENT_DECREASE,
        0.5,
        "gh-7475 exact boundary",
    );

    let settings_below_boundary = WolfeTwoSettings {
        maximum_step: Some(29.0),
        ..settings_at_boundary
    };
    let rejected = line_search_wolfe2_with_analytic_gradient(
        &objective,
        &gradient,
        point,
        direction,
        settings_below_boundary,
    );
    assert!(rejected.is_none());

    let settings_too_few_iterations = WolfeTwoSettings {
        maximum_step: None,
        maximum_iterations: 5,
        ..settings_at_boundary
    };
    let exhausted = line_search_wolfe2_with_analytic_gradient(
        &objective,
        &gradient,
        point,
        direction,
        settings_too_few_iterations,
    )
    .expect("SciPy returns its last trial when maxiter is exhausted");
    // With no `amax`, SciPy doubles the fifth evaluated step of 16 to the unverified return step 32
    assert_eq!(exhausted.step, 32.0);
    assert!(exhausted.gradient.is_none());
}

#[test]
fn scipy_test_wolfe_terminate() {
    // SciPy requires both searches to accept this immediately suitable first trial using at most
    // the two initial and two trial evaluations.
    for search_name in ["Wolfe-1", "Wolfe-2"] {
        let evaluation_count = Cell::new(0_usize);
        let objective = |step: f64| {
            evaluation_count.set(evaluation_count.get() + 1);
            -step + 0.05 * step.powi(2)
        };
        let derivative = |step: f64| {
            evaluation_count.set(evaluation_count.get() + 1);
            -1.0 + 0.1 * step
        };
        let initial_value = objective(0.0);
        let initial_derivative = derivative(0.0);

        let accepted_step = if search_name == "Wolfe-1" {
            scalar_search_wolfe1(
                &objective,
                &derivative,
                initial_value,
                None,
                initial_derivative,
            )
            .map(|accepted| accepted.step)
        } else {
            scalar_search_wolfe2_with_analytic_derivative(
                &objective,
                &derivative,
                initial_value,
                None,
                initial_derivative,
                SCIPY_SCALAR_WOLFE_TWO_SETTINGS,
            )
            .map(|accepted| accepted.step)
        }
        .unwrap_or_else(|| panic!("{search_name} should accept the suitable first trial"));

        assert!(
            evaluation_count.get() <= 4,
            "{search_name} used {} objective/derivative evaluations",
            evaluation_count.get()
        );
        assert_strong_wolfe_conditions(
            accepted_step,
            &objective,
            &derivative,
            WOLFE_SUFFICIENT_DECREASE,
            WOLFE_CURVATURE,
            search_name,
        );
    }
}

#[test]
fn scipy_test_bfgs_with_default_numerical_jacobian() -> Result<()> {
    // This is the BFGS, jac=None branch of SciPy's `test_finite_differences_jac`, using its
    // constrained-entropy fixture and objective-value expectation.
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
    let scipy_solution = [0.0, -0.524_869_316, 0.487_525_860];

    let fitted = minimize_bfgs([0.0, 0.0, 0.0], &objective)?;

    assert_allclose_with_absolute_tolerance(
        objective(fitted)?,
        objective(scipy_solution)?,
        1.0e-6,
    );
    Ok(())
}
