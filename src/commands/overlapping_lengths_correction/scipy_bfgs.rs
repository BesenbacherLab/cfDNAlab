// This module is a Rust translation of parts of SciPy 1.13.0. The SciPy copyright, license,
// disclaimer, and source-specific notices are reproduced in `THIRD_PARTY_LICENSES.md`
// Any defects or behavioral differences introduced by this Rust translation are the responsibility
// of the cfDNAlab developers, not the SciPy or MINPACK contributors
//
// The translated sources are `scipy/optimize/_optimize.py`, `_linesearch.py`, `_dcsrch.py`, and
// the numerical-differentiation path used when `_optimize.py` receives no Jacobian
//
// ******NOTICE***************
// SciPy's scipy/optimize/_optimize.py module was originally written by Travis E. Oliphant
//
// You may copy and use that SciPy module as you see fit with no
// guarantee implied provided you keep this notice in all copies.
// *****END NOTICE************
//
// DCSRCH and DCSTEP provenance retained from SciPy's `_dcsrch.py`:
// MINPACK-1 Project. June 1983. Argonne National Laboratory.
// Jorge J. More' and David J. Thuente.
// MINPACK-2 Project. November 1993. Argonne National Laboratory and University of Minnesota.
// Brett M. Averick, Richard G. Carter, and Jorge J. More'.

//! Fixed-size Rust port of the SciPy 1.13.0 BFGS path used by LIONHEART.
//!
//! LIONHEART calls `scipy.optimize.minimize` without a method, bounds, or Jacobian. SciPy therefore
//! selects BFGS, estimates gradients with absolute forward differences, tries its MINPACK-derived
//! DCSRCH Wolfe-1 line search, and falls back to its pure Python Wolfe-2 line search. This module
//! keeps that sequence and its defaults together so model code does not obscure optimizer details.
//! cfDNAlab has a stricter error contract for invalid optimizer states. It rejects a Wolfe-2 step
//! when the iteration limit is reached before the strong-Wolfe conditions are established, and it
//! rejects a non-finite or non-positive BFGS update curvature instead of substituting a value.

use anyhow::{Result, bail, ensure};

const PARAMETER_COUNT: usize = 3;
const DEFAULT_MAXIMUM_ITERATIONS: usize = PARAMETER_COUNT * 200;
const DEFAULT_GRADIENT_TOLERANCE: f64 = 1.0e-5;
const FORWARD_DIFFERENCE_STEP: f64 = 1.490_116_119_384_765_6e-8;
const WOLFE_SUFFICIENT_DECREASE: f64 = 1.0e-4;
const WOLFE_CURVATURE: f64 = 0.9;
const WOLFE_ONE_RELATIVE_STEP_TOLERANCE: f64 = 1.0e-14;
const WOLFE_ONE_MAXIMUM_ITERATIONS: usize = 100;
const WOLFE_TWO_MAXIMUM_ITERATIONS: usize = 10;
const WOLFE_TWO_ZOOM_MAXIMUM_ITERATIONS: usize = 10;
const MINIMUM_LINE_SEARCH_STEP: f64 = 1.0e-100;
const MAXIMUM_LINE_SEARCH_STEP: f64 = 1.0e100;

/// Minimize a three-parameter objective with SciPy 1.13.0's default BFGS algorithm.
///
/// The implementation follows `_minimize_bfgs` for the exact three-parameter case needed by the
/// overlapping fragment length model. It uses an identity inverse Hessian, the infinity norm for
/// convergence, an absolute `sqrt(f64::EPSILON)` forward-difference step, 600 iterations, and
/// SciPy's Wolfe-1 then Wolfe-2 line-search sequence. LIONHEART does not supply any optimizer
/// options, so these are the settings its fit receives.
pub(super) fn minimize_bfgs(
    mut point: [f64; PARAMETER_COUNT],
    evaluate: &impl Fn([f64; PARAMETER_COUNT]) -> Result<f64>,
) -> Result<[f64; PARAMETER_COUNT]> {
    let mut value = evaluate(point)?;
    ensure!(value.is_finite(), "BFGS initial objective is non-finite");
    let mut gradient = forward_difference_gradient(point, value, evaluate)?;
    ensure!(
        gradient.iter().all(|component| component.is_finite()),
        "BFGS initial finite-difference gradient is non-finite"
    );

    let mut inverse_hessian = identity_matrix();
    // SciPy seeds line-search history so the first proposed parameter displacement is near one
    let mut previous_value = value + euclidean_norm(gradient) / 2.0;
    let mut gradient_norm = infinity_norm(gradient);

    for _iteration in 0..DEFAULT_MAXIMUM_ITERATIONS {
        if gradient_norm <= DEFAULT_GRADIENT_TOLERANCE {
            return Ok(point);
        }

        let direction =
            matrix_vector_product(inverse_hessian, gradient).map(|component| -component);
        let accepted =
            line_search_wolfe12(point, direction, gradient, value, previous_value, evaluate)?;
        ensure!(
            accepted.value.is_finite(),
            "BFGS accepted a non-finite objective value"
        );
        ensure!(
            accepted
                .gradient
                .iter()
                .all(|component| component.is_finite()),
            "BFGS accepted a non-finite finite-difference gradient"
        );
        let step_vector = direction.map(|component| accepted.step * component);
        let next_point = add(point, step_vector);
        let gradient_change = subtract(accepted.gradient, gradient);

        point = next_point;
        previous_value = value;
        value = accepted.value;
        gradient = accepted.gradient;
        gradient_norm = infinity_norm(gradient);

        if gradient_norm <= DEFAULT_GRADIENT_TOLERANCE {
            return Ok(point);
        }
        // SciPy's xrtol default is zero, so only an exactly zero displacement satisfies this check
        if accepted.step * euclidean_norm(direction) <= 0.0 {
            return Ok(point);
        }
        let curvature = dot(gradient_change, step_vector);
        inverse_hessian =
            bfgs_inverse_hessian_update(inverse_hessian, step_vector, gradient_change, curvature)?;
    }

    bail!(
        "overlapping fragment length mixture fit did not converge after {} BFGS iterations",
        DEFAULT_MAXIMUM_ITERATIONS
    )
}

/// Result shared by the primary and fallback line searches.
#[derive(Clone, Copy, Debug)]
struct LineSearchPoint {
    step: f64,
    value: f64,
    gradient: [f64; PARAMETER_COUNT],
}

/// Wolfe-2 result before its optional derivative has been resolved.
///
/// A converged search contains the gradient used to verify the curvature condition. SciPy instead
/// returns its last step without a derivative when the outer Wolfe-2 loop reaches `maxiter`.
/// SciPy's BFGS caller then evaluates that gradient before updating its inverse Hessian. Retaining
/// the optional gradient reproduces that scalar-search result for upstream regression tests, while
/// cfDNAlab's production wrapper rejects it.
#[derive(Clone, Copy, Debug)]
struct WolfeTwoPoint {
    step: f64,
    value: f64,
    gradient: Option<[f64; PARAMETER_COUNT]>,
}

/// Settings accepted by SciPy's scalar Wolfe-2 implementation.
///
/// LIONHEART's BFGS call uses the production constant below. Keeping the settings explicit allows
/// upstream SciPy boundary regressions to exercise their smaller `amax`, altered `c2`, and reduced
/// iteration-count cases without changing production behavior.
#[derive(Clone, Copy, Debug)]
struct WolfeTwoSettings {
    sufficient_decrease: f64,
    curvature: f64,
    maximum_step: Option<f64>,
    maximum_iterations: usize,
    zoom_maximum_iterations: usize,
}

const BFGS_WOLFE_TWO_SETTINGS: WolfeTwoSettings = WolfeTwoSettings {
    sufficient_decrease: WOLFE_SUFFICIENT_DECREASE,
    curvature: WOLFE_CURVATURE,
    maximum_step: Some(MAXIMUM_LINE_SEARCH_STEP),
    maximum_iterations: WOLFE_TWO_MAXIMUM_ITERATIONS,
    zoom_maximum_iterations: WOLFE_TWO_ZOOM_MAXIMUM_ITERATIONS,
};

/// Settings for SciPy's MINPACK-derived scalar Wolfe-1 search.
#[derive(Clone, Copy, Debug)]
struct WolfeOneSettings {
    sufficient_decrease: f64,
    curvature: f64,
    relative_step_tolerance: f64,
    minimum_step: f64,
    maximum_step: f64,
    maximum_iterations: usize,
}

const BFGS_WOLFE_ONE_SETTINGS: WolfeOneSettings = WolfeOneSettings {
    sufficient_decrease: WOLFE_SUFFICIENT_DECREASE,
    curvature: WOLFE_CURVATURE,
    relative_step_tolerance: WOLFE_ONE_RELATIVE_STEP_TOLERANCE,
    minimum_step: MINIMUM_LINE_SEARCH_STEP,
    maximum_step: MAXIMUM_LINE_SEARCH_STEP,
    maximum_iterations: WOLFE_ONE_MAXIMUM_ITERATIONS,
};

/// Run SciPy's Wolfe-1 search and use Wolfe-2 only when Wolfe-1 cannot find a step.
fn line_search_wolfe12(
    point: [f64; PARAMETER_COUNT],
    direction: [f64; PARAMETER_COUNT],
    gradient: [f64; PARAMETER_COUNT],
    value: f64,
    previous_value: f64,
    evaluate: &impl Fn([f64; PARAMETER_COUNT]) -> Result<f64>,
) -> Result<LineSearchPoint> {
    if let Some(accepted) =
        line_search_wolfe1(point, direction, gradient, value, previous_value, evaluate)?
    {
        return Ok(accepted);
    }
    if let Some(accepted) = line_search_wolfe2(
        point,
        direction,
        gradient,
        value,
        Some(previous_value),
        BFGS_WOLFE_TWO_SETTINGS,
        evaluate,
    )? {
        return Ok(accepted);
    }
    bail!("SciPy BFGS Wolfe-1 and Wolfe-2 line searches both failed")
}

/// Estimate the gradient exactly as SciPy does when BFGS receives `jac=None`.
fn forward_difference_gradient(
    point: [f64; PARAMETER_COUNT],
    value: f64,
    evaluate: &impl Fn([f64; PARAMETER_COUNT]) -> Result<f64>,
) -> Result<[f64; PARAMETER_COUNT]> {
    let mut gradient = [0.0; PARAMETER_COUNT];
    for parameter_index in 0..PARAMETER_COUNT {
        let mut shifted = point;
        shifted[parameter_index] += FORWARD_DIFFERENCE_STEP;
        // SciPy divides by the representable displacement after floating-point addition
        let mut actual_step = shifted[parameter_index] - point[parameter_index];
        if actual_step == 0.0 {
            // SciPy substitutes an automatically scaled step when the requested absolute step is
            // too small to change a large floating-point parameter
            let direction = if point[parameter_index] >= 0.0 {
                1.0
            } else {
                -1.0
            };
            actual_step =
                FORWARD_DIFFERENCE_STEP * direction * point[parameter_index].abs().max(1.0);
            shifted[parameter_index] = point[parameter_index] + actual_step;
            actual_step = shifted[parameter_index] - point[parameter_index];
        }
        let shifted_value = evaluate(shifted)?;
        gradient[parameter_index] = (shifted_value - value) / actual_step;
    }
    Ok(gradient)
}

/// Evaluate the objective and its forward-difference gradient at a line-search step.
fn evaluate_step(
    point: [f64; PARAMETER_COUNT],
    direction: [f64; PARAMETER_COUNT],
    step: f64,
    evaluate: &impl Fn([f64; PARAMETER_COUNT]) -> Result<f64>,
) -> Result<LineSearchPoint> {
    let candidate = add_scaled(point, direction, step);
    let value = evaluate(candidate)?;
    let gradient = forward_difference_gradient(candidate, value, evaluate)?;
    Ok(LineSearchPoint {
        step,
        value,
        gradient,
    })
}

/// Calculate SciPy's history-dependent initial line-search step.
fn initial_line_search_step(value: f64, previous_value: f64, derivative: f64) -> f64 {
    if derivative != 0.0 {
        let estimated_step = (1.01 * 2.0 * (value - previous_value) / derivative).min(1.0);
        if estimated_step >= 0.0 {
            return estimated_step;
        }
    }
    1.0
}

/// Run SciPy's primary MINPACK-derived DCSRCH Wolfe-1 line search.
fn line_search_wolfe1(
    point: [f64; PARAMETER_COUNT],
    direction: [f64; PARAMETER_COUNT],
    gradient: [f64; PARAMETER_COUNT],
    value: f64,
    previous_value: f64,
    evaluate: &impl Fn([f64; PARAMETER_COUNT]) -> Result<f64>,
) -> Result<Option<LineSearchPoint>> {
    let initial_derivative = dot(gradient, direction);
    let mut step = initial_line_search_step(value, previous_value, initial_derivative);
    let Some(mut state) =
        DcsrchState::new(value, initial_derivative, step, BFGS_WOLFE_ONE_SETTINGS)
    else {
        return Ok(None);
    };
    let mut candidate = LineSearchPoint {
        step,
        value,
        gradient,
    };

    for iteration in 0..BFGS_WOLFE_ONE_SETTINGS.maximum_iterations {
        let derivative = if iteration == 0 {
            initial_derivative
        } else {
            dot(candidate.gradient, direction)
        };
        match state.iterate(step, candidate.value, derivative, iteration == 0) {
            DcsrchAction::Evaluate(next_step) => {
                if !next_step.is_finite() {
                    return Ok(None);
                }
                step = next_step;
                candidate = evaluate_step(point, direction, step, evaluate)?;
            }
            DcsrchAction::Converged => return Ok(Some(candidate)),
            DcsrchAction::Failed => return Ok(None),
        }
    }
    Ok(None)
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum DcsrchAction {
    Evaluate(f64),
    Converged,
    Failed,
}

/// Mutable search interval used by SciPy's port of MINPACK DCSRCH.
///
/// `best_*` corresponds to SciPy's `stx`, `fx`, and `gx`. `other_*` corresponds to `sty`, `fy`,
/// and `gy`. Keeping descriptive names makes the four safeguarded interpolation cases readable
/// without changing their arithmetic.
#[derive(Clone, Copy, Debug)]
struct DcsrchState {
    settings: WolfeOneSettings,
    stage: u8,
    bracketed: bool,
    initial_value: f64,
    initial_derivative: f64,
    sufficient_decrease_test: f64,
    best_step: f64,
    best_value: f64,
    best_derivative: f64,
    other_step: f64,
    other_value: f64,
    other_derivative: f64,
    allowed_minimum: f64,
    allowed_maximum: f64,
    interval_width: f64,
    previous_interval_width: f64,
}

impl DcsrchState {
    fn new(
        initial_value: f64,
        initial_derivative: f64,
        initial_step: f64,
        settings: WolfeOneSettings,
    ) -> Option<Self> {
        if initial_step < settings.minimum_step
            || initial_step > settings.maximum_step
            || initial_derivative >= 0.0
        {
            return None;
        }
        let interval_width = settings.maximum_step - settings.minimum_step;
        Some(Self {
            settings,
            stage: 1,
            bracketed: false,
            initial_value,
            initial_derivative,
            sufficient_decrease_test: settings.sufficient_decrease * initial_derivative,
            best_step: 0.0,
            best_value: initial_value,
            best_derivative: initial_derivative,
            other_step: 0.0,
            other_value: initial_value,
            other_derivative: initial_derivative,
            allowed_minimum: 0.0,
            allowed_maximum: initial_step + 4.0 * initial_step,
            interval_width,
            previous_interval_width: interval_width / 0.5,
        })
    }

    /// Advance DCSRCH after receiving the function value and directional derivative it requested.
    fn iterate(
        &mut self,
        mut step: f64,
        value: f64,
        derivative: f64,
        starting: bool,
    ) -> DcsrchAction {
        if starting {
            return DcsrchAction::Evaluate(step);
        }

        let armijo_limit = self.initial_value + step * self.sufficient_decrease_test;
        if self.stage == 1 && value <= armijo_limit && derivative >= 0.0 {
            self.stage = 2;
        }

        let rounding_stalled =
            self.bracketed && (step <= self.allowed_minimum || step >= self.allowed_maximum);
        let interval_too_small = self.bracketed
            && self.allowed_maximum - self.allowed_minimum
                <= self.settings.relative_step_tolerance * self.allowed_maximum;
        let maximum_step_stalled = step == self.settings.maximum_step
            && value <= armijo_limit
            && derivative <= self.sufficient_decrease_test;
        let minimum_step_stalled = step == self.settings.minimum_step
            && (value > armijo_limit || derivative >= self.sufficient_decrease_test);
        let converged = value <= armijo_limit
            && derivative.abs() <= self.settings.curvature * -self.initial_derivative;
        if converged {
            return DcsrchAction::Converged;
        }
        if rounding_stalled || interval_too_small || maximum_step_stalled || minimum_step_stalled {
            return DcsrchAction::Failed;
        }

        if self.stage == 1 && value <= self.best_value && value > armijo_limit {
            // DCSRCH first optimizes a modified function until it has sufficient decrease
            let modified_value = value - step * self.sufficient_decrease_test;
            let mut modified_best_value =
                self.best_value - self.best_step * self.sufficient_decrease_test;
            let mut modified_other_value =
                self.other_value - self.other_step * self.sufficient_decrease_test;
            let modified_derivative = derivative - self.sufficient_decrease_test;
            let mut modified_best_derivative = self.best_derivative - self.sufficient_decrease_test;
            let mut modified_other_derivative =
                self.other_derivative - self.sufficient_decrease_test;
            let updated = safeguarded_dcstep(
                self.best_step,
                modified_best_value,
                modified_best_derivative,
                self.other_step,
                modified_other_value,
                modified_other_derivative,
                step,
                modified_value,
                modified_derivative,
                self.bracketed,
                self.allowed_minimum,
                self.allowed_maximum,
            );
            self.best_step = updated.best_step;
            modified_best_value = updated.best_value;
            modified_best_derivative = updated.best_derivative;
            self.other_step = updated.other_step;
            modified_other_value = updated.other_value;
            modified_other_derivative = updated.other_derivative;
            step = updated.next_step;
            self.bracketed = updated.bracketed;
            self.best_value = modified_best_value + self.best_step * self.sufficient_decrease_test;
            self.best_derivative = modified_best_derivative + self.sufficient_decrease_test;
            self.other_value =
                modified_other_value + self.other_step * self.sufficient_decrease_test;
            self.other_derivative = modified_other_derivative + self.sufficient_decrease_test;
        } else {
            let updated = safeguarded_dcstep(
                self.best_step,
                self.best_value,
                self.best_derivative,
                self.other_step,
                self.other_value,
                self.other_derivative,
                step,
                value,
                derivative,
                self.bracketed,
                self.allowed_minimum,
                self.allowed_maximum,
            );
            self.best_step = updated.best_step;
            self.best_value = updated.best_value;
            self.best_derivative = updated.best_derivative;
            self.other_step = updated.other_step;
            self.other_value = updated.other_value;
            self.other_derivative = updated.other_derivative;
            step = updated.next_step;
            self.bracketed = updated.bracketed;
        }

        if self.bracketed {
            if (self.other_step - self.best_step).abs() >= 0.66 * self.previous_interval_width {
                step = self.best_step + 0.5 * (self.other_step - self.best_step);
            }
            self.previous_interval_width = self.interval_width;
            self.interval_width = (self.other_step - self.best_step).abs();
        }

        if self.bracketed {
            self.allowed_minimum = self.best_step.min(self.other_step);
            self.allowed_maximum = self.best_step.max(self.other_step);
        } else {
            self.allowed_minimum = step + 1.1 * (step - self.best_step);
            self.allowed_maximum = step + 4.0 * (step - self.best_step);
        }
        step = step.clamp(self.settings.minimum_step, self.settings.maximum_step);
        if self.bracketed
            && (step <= self.allowed_minimum
                || step >= self.allowed_maximum
                || self.allowed_maximum - self.allowed_minimum
                    <= self.settings.relative_step_tolerance * self.allowed_maximum)
        {
            step = self.best_step;
        }
        DcsrchAction::Evaluate(step)
    }
}

#[derive(Clone, Copy, Debug)]
struct DcstepResult {
    best_step: f64,
    best_value: f64,
    best_derivative: f64,
    other_step: f64,
    other_value: f64,
    other_derivative: f64,
    next_step: f64,
    bracketed: bool,
}

/// Compute the safeguarded interpolation step used by DCSRCH.
///
/// The four branches are the four cases in SciPy 1.13.0's `dcstep`: a higher function value,
/// opposite derivative signs, a smaller derivative magnitude with the same sign, and a derivative
/// magnitude that does not decrease. The formulas are kept in the same evaluation order.
#[allow(clippy::too_many_arguments)]
fn safeguarded_dcstep(
    mut best_step: f64,
    mut best_value: f64,
    mut best_derivative: f64,
    mut other_step: f64,
    mut other_value: f64,
    mut other_derivative: f64,
    current_step: f64,
    current_value: f64,
    current_derivative: f64,
    mut bracketed: bool,
    step_minimum: f64,
    step_maximum: f64,
) -> DcstepResult {
    let derivative_sign_changed =
        numeric_sign(current_derivative) * numeric_sign(best_derivative) < 0.0;
    let next_step;

    if current_value > best_value {
        let theta = 3.0 * (best_value - current_value) / (current_step - best_step)
            + best_derivative
            + current_derivative;
        let scale = max_three(theta.abs(), best_derivative.abs(), current_derivative.abs());
        let mut gamma = scale
            * ((theta / scale).powi(2) - (best_derivative / scale) * (current_derivative / scale))
                .sqrt();
        if current_step < best_step {
            gamma *= -1.0;
        }
        let cubic_ratio = ((gamma - best_derivative) + theta)
            / (((gamma - best_derivative) + gamma) + current_derivative);
        let cubic_step = best_step + cubic_ratio * (current_step - best_step);
        let quadratic_step = best_step
            + ((best_derivative
                / ((best_value - current_value) / (current_step - best_step) + best_derivative))
                / 2.0)
                * (current_step - best_step);
        next_step = if (cubic_step - best_step).abs() <= (quadratic_step - best_step).abs() {
            cubic_step
        } else {
            cubic_step + (quadratic_step - cubic_step) / 2.0
        };
        bracketed = true;
    } else if derivative_sign_changed {
        let theta = 3.0 * (best_value - current_value) / (current_step - best_step)
            + best_derivative
            + current_derivative;
        let scale = max_three(theta.abs(), best_derivative.abs(), current_derivative.abs());
        let mut gamma = scale
            * ((theta / scale).powi(2) - (best_derivative / scale) * (current_derivative / scale))
                .sqrt();
        if current_step > best_step {
            gamma *= -1.0;
        }
        let cubic_ratio = ((gamma - current_derivative) + theta)
            / (((gamma - current_derivative) + gamma) + best_derivative);
        let cubic_step = current_step + cubic_ratio * (best_step - current_step);
        let secant_step = current_step
            + (current_derivative / (current_derivative - best_derivative))
                * (best_step - current_step);
        next_step = if (cubic_step - current_step).abs() > (secant_step - current_step).abs() {
            cubic_step
        } else {
            secant_step
        };
        bracketed = true;
    } else if current_derivative.abs() < best_derivative.abs() {
        let theta = 3.0 * (best_value - current_value) / (current_step - best_step)
            + best_derivative
            + current_derivative;
        let scale = max_three(theta.abs(), best_derivative.abs(), current_derivative.abs());
        let mut gamma = scale
            * 0.0_f64
                .max(
                    (theta / scale).powi(2)
                        - (best_derivative / scale) * (current_derivative / scale),
                )
                .sqrt();
        if current_step > best_step {
            gamma = -gamma;
        }
        let cubic_ratio = ((gamma - current_derivative) + theta)
            / ((gamma + (best_derivative - current_derivative)) + gamma);
        let cubic_step = if cubic_ratio < 0.0 && gamma != 0.0 {
            current_step + cubic_ratio * (best_step - current_step)
        } else if current_step > best_step {
            step_maximum
        } else {
            step_minimum
        };
        let secant_step = current_step
            + (current_derivative / (current_derivative - best_derivative))
                * (best_step - current_step);
        if bracketed {
            let closer_step =
                if (cubic_step - current_step).abs() < (secant_step - current_step).abs() {
                    cubic_step
                } else {
                    secant_step
                };
            next_step = if current_step > best_step {
                closer_step.min(current_step + 0.66 * (other_step - current_step))
            } else {
                closer_step.max(current_step + 0.66 * (other_step - current_step))
            };
        } else {
            let farther_step =
                if (cubic_step - current_step).abs() > (secant_step - current_step).abs() {
                    cubic_step
                } else {
                    secant_step
                };
            next_step = farther_step.clamp(step_minimum, step_maximum);
        }
    } else if bracketed {
        let theta = 3.0 * (current_value - other_value) / (other_step - current_step)
            + other_derivative
            + current_derivative;
        let scale = max_three(
            theta.abs(),
            other_derivative.abs(),
            current_derivative.abs(),
        );
        let mut gamma = scale
            * ((theta / scale).powi(2) - (other_derivative / scale) * (current_derivative / scale))
                .sqrt();
        if current_step > other_step {
            gamma = -gamma;
        }
        let cubic_ratio = ((gamma - current_derivative) + theta)
            / (((gamma - current_derivative) + gamma) + other_derivative);
        next_step = current_step + cubic_ratio * (other_step - current_step);
    } else if current_step > best_step {
        next_step = step_maximum;
    } else {
        next_step = step_minimum;
    }

    if current_value > best_value {
        other_step = current_step;
        other_value = current_value;
        other_derivative = current_derivative;
    } else {
        if derivative_sign_changed {
            other_step = best_step;
            other_value = best_value;
            other_derivative = best_derivative;
        }
        best_step = current_step;
        best_value = current_value;
        best_derivative = current_derivative;
    }

    DcstepResult {
        best_step,
        best_value,
        best_derivative,
        other_step,
        other_value,
        other_derivative,
        next_step,
        bracketed,
    }
}

fn max_three(first: f64, second: f64, third: f64) -> f64 {
    first.max(second).max(third)
}

/// Match `numpy.sign` for the finite values used by DCSRCH, including zero.
fn numeric_sign(value: f64) -> f64 {
    if value > 0.0 {
        1.0
    } else if value < 0.0 {
        -1.0
    } else {
        0.0
    }
}

/// Run SciPy's pure Rust equivalent of `scalar_search_wolfe2`.
///
/// SciPy uses this search only after DCSRCH fails. It doubles the trial step until the strong-Wolfe
/// conditions hold or a suitable step is bracketed. The zoom stage then uses safeguarded cubic,
/// quadratic, and bisection interpolation in that order.
fn line_search_wolfe2(
    point: [f64; PARAMETER_COUNT],
    direction: [f64; PARAMETER_COUNT],
    gradient: [f64; PARAMETER_COUNT],
    value: f64,
    previous_value: Option<f64>,
    settings: WolfeTwoSettings,
    evaluate: &impl Fn([f64; PARAMETER_COUNT]) -> Result<f64>,
) -> Result<Option<LineSearchPoint>> {
    let initial_derivative = dot(gradient, direction);
    let initial_step = previous_value
        .map(|previous| initial_line_search_step(value, previous, initial_derivative))
        .unwrap_or(1.0);
    let evaluate_value = |step| evaluate(add_scaled(point, direction, step));
    let evaluate_gradient = |step, candidate_value| {
        let candidate = add_scaled(point, direction, step);
        forward_difference_gradient(candidate, candidate_value, evaluate)
    };
    let accepted = scalar_search_wolfe2(
        direction,
        value,
        initial_derivative,
        initial_step,
        settings,
        &evaluate_value,
        &evaluate_gradient,
    )?;
    let Some(accepted) = accepted else {
        return Ok(None);
    };
    // SciPy returns the last trial after exhausting `maxiter`, even though the strong-Wolfe
    // conditions were not established. A fitted cfDNAlab model must not silently accept that step.
    let Some(gradient) = accepted.gradient else {
        bail!(
            "Wolfe-2 line search reached its limit of {} iterations without establishing the strong-Wolfe conditions",
            settings.maximum_iterations
        );
    };
    Ok(Some(LineSearchPoint {
        step: accepted.step,
        value: accepted.value,
        gradient,
    }))
}

/// Run the decision loop from SciPy 1.13.0's `scalar_search_wolfe2`.
///
/// Objective and gradient evaluation remain separate because SciPy checks sufficient decrease
/// before requesting the derivative. Production constructs the derivative from the numerical
/// parameter gradient, while the dedicated SciPy regression tests supply the analytic scalar
/// derivatives used upstream. Both paths therefore execute the same evaluation order, bracketing,
/// doubling, and zoom decisions.
fn scalar_search_wolfe2(
    direction: [f64; PARAMETER_COUNT],
    initial_value: f64,
    initial_derivative: f64,
    initial_step: f64,
    settings: WolfeTwoSettings,
    evaluate_value: &impl Fn(f64) -> Result<f64>,
    evaluate_gradient: &impl Fn(f64, f64) -> Result<[f64; PARAMETER_COUNT]>,
) -> Result<Option<WolfeTwoPoint>> {
    let mut lower_step = 0.0;
    let mut trial_step = settings
        .maximum_step
        .map_or(initial_step, |maximum_step| initial_step.min(maximum_step));
    let mut lower_value = initial_value;
    let mut lower_derivative = initial_derivative;
    let mut trial_value = evaluate_value(trial_step)?;

    for iteration in 0..settings.maximum_iterations {
        if trial_step == 0.0
            || settings
                .maximum_step
                .is_some_and(|maximum_step| lower_step > maximum_step)
        {
            return Ok(None);
        }
        if trial_value
            > initial_value + settings.sufficient_decrease * trial_step * initial_derivative
            || (iteration > 0 && trial_value >= lower_value)
        {
            return zoom_wolfe2(
                direction,
                lower_step,
                trial_step,
                lower_value,
                trial_value,
                lower_derivative,
                initial_value,
                initial_derivative,
                settings,
                evaluate_value,
                evaluate_gradient,
            );
        }

        let trial_gradient = evaluate_gradient(trial_step, trial_value)?;
        let trial_derivative = dot(trial_gradient, direction);
        if trial_derivative.abs() <= -settings.curvature * initial_derivative {
            return Ok(Some(WolfeTwoPoint {
                step: trial_step,
                value: trial_value,
                gradient: Some(trial_gradient),
            }));
        }
        if trial_derivative >= 0.0 {
            return zoom_wolfe2(
                direction,
                trial_step,
                lower_step,
                trial_value,
                lower_value,
                trial_derivative,
                initial_value,
                initial_derivative,
                settings,
                evaluate_value,
                evaluate_gradient,
            );
        }

        let next_step = settings
            .maximum_step
            .map_or(2.0 * trial_step, |maximum_step| {
                (2.0 * trial_step).min(maximum_step)
            });
        lower_step = trial_step;
        lower_value = trial_value;
        lower_derivative = trial_derivative;
        trial_step = next_step;
        trial_value = evaluate_value(trial_step)?;
    }
    // This is SciPy's non-converged `for ... else` return: the step and objective remain usable,
    // but the absent derivative tells the caller that the strong-Wolfe conditions were not proven.
    Ok(Some(WolfeTwoPoint {
        step: trial_step,
        value: trial_value,
        gradient: None,
    }))
}

/// Resolve a Wolfe-2 bracket using SciPy's interpolation safeguards.
#[allow(clippy::too_many_arguments)]
fn zoom_wolfe2(
    direction: [f64; PARAMETER_COUNT],
    mut lower_step: f64,
    mut upper_step: f64,
    mut lower_value: f64,
    mut upper_value: f64,
    mut lower_derivative: f64,
    initial_value: f64,
    initial_derivative: f64,
    settings: WolfeTwoSettings,
    evaluate_value: &impl Fn(f64) -> Result<f64>,
    evaluate_gradient: &impl Fn(f64, f64) -> Result<[f64; PARAMETER_COUNT]>,
) -> Result<Option<WolfeTwoPoint>> {
    const CUBIC_BOUNDARY_FRACTION: f64 = 0.2;
    const QUADRATIC_BOUNDARY_FRACTION: f64 = 0.1;
    let mut recent_step = 0.0;
    let mut recent_value = initial_value;

    for iteration in 0..=settings.zoom_maximum_iterations {
        let signed_width = upper_step - lower_step;
        let interval_start = lower_step.min(upper_step);
        let interval_end = lower_step.max(upper_step);

        let mut trial_step = if iteration > 0 {
            cubic_interpolant_minimum(
                lower_step,
                lower_value,
                lower_derivative,
                upper_step,
                upper_value,
                recent_step,
                recent_value,
            )
        } else {
            None
        };
        let cubic_margin = CUBIC_BOUNDARY_FRACTION * signed_width;
        if trial_step.is_none_or(|step| {
            step > interval_end - cubic_margin || step < interval_start + cubic_margin
        }) {
            trial_step = quadratic_interpolant_minimum(
                lower_step,
                lower_value,
                lower_derivative,
                upper_step,
                upper_value,
            );
            let quadratic_margin = QUADRATIC_BOUNDARY_FRACTION * signed_width;
            if trial_step.is_none_or(|step| {
                step > interval_end - quadratic_margin || step < interval_start + quadratic_margin
            }) {
                trial_step = Some(lower_step + 0.5 * signed_width);
            }
        }

        let Some(trial_step) = trial_step else {
            return Ok(None);
        };
        let trial_value = evaluate_value(trial_step)?;
        if trial_value
            > initial_value + settings.sufficient_decrease * trial_step * initial_derivative
            || trial_value >= lower_value
        {
            recent_step = upper_step;
            recent_value = upper_value;
            upper_step = trial_step;
            upper_value = trial_value;
            continue;
        }

        let trial_gradient = evaluate_gradient(trial_step, trial_value)?;
        let trial_derivative = dot(trial_gradient, direction);
        if trial_derivative.abs() <= -settings.curvature * initial_derivative {
            return Ok(Some(WolfeTwoPoint {
                step: trial_step,
                value: trial_value,
                gradient: Some(trial_gradient),
            }));
        }
        if trial_derivative * (upper_step - lower_step) >= 0.0 {
            recent_step = upper_step;
            recent_value = upper_value;
            upper_step = lower_step;
            upper_value = lower_value;
        } else {
            recent_step = lower_step;
            recent_value = lower_value;
        }
        lower_step = trial_step;
        lower_value = trial_value;
        lower_derivative = trial_derivative;
    }
    Ok(None)
}

/// Find the stationary point of SciPy's cubic interpolation polynomial.
fn cubic_interpolant_minimum(
    anchor_step: f64,
    anchor_value: f64,
    anchor_derivative: f64,
    second_step: f64,
    second_value: f64,
    third_step: f64,
    third_value: f64,
) -> Option<f64> {
    let second_distance = second_step - anchor_step;
    let third_distance = third_step - anchor_step;
    let denominator =
        (second_distance * third_distance).powi(2) * (second_distance - third_distance);
    let second_remainder = second_value - anchor_value - anchor_derivative * second_distance;
    let third_remainder = third_value - anchor_value - anchor_derivative * third_distance;
    let cubic = (third_distance.powi(2) * second_remainder
        - second_distance.powi(2) * third_remainder)
        / denominator;
    let quadratic = (-third_distance.powi(3) * second_remainder
        + second_distance.powi(3) * third_remainder)
        / denominator;
    let radical = quadratic * quadratic - 3.0 * cubic * anchor_derivative;
    let minimum = anchor_step + (-quadratic + radical.sqrt()) / (3.0 * cubic);
    minimum.is_finite().then_some(minimum)
}

/// Find the stationary point of SciPy's quadratic interpolation polynomial.
fn quadratic_interpolant_minimum(
    anchor_step: f64,
    anchor_value: f64,
    anchor_derivative: f64,
    second_step: f64,
    second_value: f64,
) -> Option<f64> {
    let distance = second_step - anchor_step;
    let quadratic = (second_value - anchor_value - anchor_derivative * distance) / distance.powi(2);
    let minimum = anchor_step - anchor_derivative / (2.0 * quadratic);
    minimum.is_finite().then_some(minimum)
}

/// Apply SciPy's rank-two inverse-Hessian BFGS update.
///
/// A valid strong-Wolfe step gives positive update curvature `y^T s`. SciPy substitutes an inverse
/// curvature of `1000` when the curvature is exactly zero. cfDNAlab instead rejects non-finite or
/// non-positive curvature because those values cannot produce a trustworthy positive-definite
/// inverse-Hessian update.
fn bfgs_inverse_hessian_update(
    inverse_hessian: [[f64; PARAMETER_COUNT]; PARAMETER_COUNT],
    step: [f64; PARAMETER_COUNT],
    gradient_change: [f64; PARAMETER_COUNT],
    curvature: f64,
) -> Result<[[f64; PARAMETER_COUNT]; PARAMETER_COUNT]> {
    ensure!(
        curvature.is_finite() && curvature > 0.0,
        "BFGS update curvature must be finite and positive, got {}",
        curvature
    );
    let inverse_curvature = 1.0 / curvature;
    ensure!(
        inverse_curvature.is_finite(),
        "BFGS inverse update curvature is non-finite for curvature {}",
        curvature
    );
    let mut left = identity_matrix();
    let mut right = identity_matrix();
    for row in 0..PARAMETER_COUNT {
        for column in 0..PARAMETER_COUNT {
            left[row][column] -= step[row] * gradient_change[column] * inverse_curvature;
            right[row][column] -= gradient_change[row] * step[column] * inverse_curvature;
        }
    }
    // Preserve NumPy's `A1 @ (H @ A2)` parenthesization because floating-point matrix
    // multiplication is not exactly associative
    let mut updated = matrix_multiply(left, matrix_multiply(inverse_hessian, right));
    for row in 0..PARAMETER_COUNT {
        for column in 0..PARAMETER_COUNT {
            updated[row][column] += inverse_curvature * step[row] * step[column];
        }
    }
    Ok(updated)
}

fn identity_matrix() -> [[f64; PARAMETER_COUNT]; PARAMETER_COUNT] {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

fn matrix_vector_product(
    matrix: [[f64; PARAMETER_COUNT]; PARAMETER_COUNT],
    vector: [f64; PARAMETER_COUNT],
) -> [f64; PARAMETER_COUNT] {
    matrix.map(|row| dot(row, vector))
}

fn matrix_multiply(
    left: [[f64; PARAMETER_COUNT]; PARAMETER_COUNT],
    right: [[f64; PARAMETER_COUNT]; PARAMETER_COUNT],
) -> [[f64; PARAMETER_COUNT]; PARAMETER_COUNT] {
    let mut output = [[0.0; PARAMETER_COUNT]; PARAMETER_COUNT];
    for row in 0..PARAMETER_COUNT {
        for column in 0..PARAMETER_COUNT {
            output[row][column] = (0..PARAMETER_COUNT)
                .map(|inner| left[row][inner] * right[inner][column])
                .sum();
        }
    }
    output
}

fn add(left: [f64; PARAMETER_COUNT], right: [f64; PARAMETER_COUNT]) -> [f64; PARAMETER_COUNT] {
    [left[0] + right[0], left[1] + right[1], left[2] + right[2]]
}

fn add_scaled(
    point: [f64; PARAMETER_COUNT],
    direction: [f64; PARAMETER_COUNT],
    scale: f64,
) -> [f64; PARAMETER_COUNT] {
    [
        point[0] + scale * direction[0],
        point[1] + scale * direction[1],
        point[2] + scale * direction[2],
    ]
}

fn subtract(left: [f64; PARAMETER_COUNT], right: [f64; PARAMETER_COUNT]) -> [f64; PARAMETER_COUNT] {
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
}

fn dot(left: [f64; PARAMETER_COUNT], right: [f64; PARAMETER_COUNT]) -> f64 {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}

fn euclidean_norm(vector: [f64; PARAMETER_COUNT]) -> f64 {
    dot(vector, vector).sqrt()
}

fn infinity_norm(vector: [f64; PARAMETER_COUNT]) -> f64 {
    if vector.iter().any(|component| component.is_nan()) {
        f64::NAN
    } else {
        vector.into_iter().map(f64::abs).fold(0.0, f64::max)
    }
}

#[cfg(test)]
mod scipy_1_13_tests {
    include!("scipy_bfgs_scipy_1_13_tests.rs");
}

#[cfg(test)]
mod tests {
    include!("scipy_bfgs_tests.rs");
}
