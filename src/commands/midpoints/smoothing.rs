use anyhow::{Result, ensure};
use std::{fmt, str::FromStr};

pub const DEFAULT_SAVGOL_WINDOW_BP: u32 = 165;
pub const SAVGOL_POLYNOMIAL_ORDER: u32 = 3;

/// Smoothing mode for final midpoint profiles.
///
/// `Raw` is the default and leaves grouped midpoint profiles unsmoothed. `SavGolDefault` selects
/// the documented Savitzky-Golay window, while `SavGol` stores an explicit odd window size in base
/// pairs. Polynomial order is fixed at 3 so the command exposes one signal-scale choice instead of
/// several coupled filter parameters.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MidpointSmoothing {
    Raw,
    SavGolDefault,
    SavGol { window_bp: u32 },
}

impl Default for MidpointSmoothing {
    fn default() -> Self {
        Self::Raw
    }
}

impl fmt::Display for MidpointSmoothing {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Raw => formatter.write_str("raw"),
            Self::SavGolDefault => formatter.write_str("savgol"),
            Self::SavGol { window_bp } => write!(formatter, "savgol={window_bp}"),
        }
    }
}

impl FromStr for MidpointSmoothing {
    type Err = String;

    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        parse_midpoint_smoothing(value)
    }
}

/// Parse the `--smooth` value accepted by the midpoint CLI.
///
/// The parser intentionally accepts `raw` so configuration round-trips and tests can name the raw
/// mode, even though users normally omit `--smooth` for raw output.
pub fn parse_midpoint_smoothing(value: &str) -> std::result::Result<MidpointSmoothing, String> {
    let trimmed = value.trim();
    if trimmed.eq_ignore_ascii_case("raw") {
        return Ok(MidpointSmoothing::Raw);
    }
    if trimmed.eq_ignore_ascii_case("savgol") {
        return Ok(MidpointSmoothing::SavGolDefault);
    }

    let Some((method, window_text)) = trimmed.split_once('=') else {
        return Err(format!(
            "unsupported smoothing mode '{value}'. Use 'savgol' or 'savgol=<odd_bp>'"
        ));
    };
    if !method.eq_ignore_ascii_case("savgol") {
        return Err(format!(
            "unsupported smoothing method '{method}'. Use 'savgol'"
        ));
    }

    let window_bp = window_text.parse::<u32>().map_err(|_| {
        format!("invalid Savitzky-Golay window '{window_text}'. Use an odd integer in bp")
    })?;
    if window_bp % 2 == 0 {
        return Err(format!(
            "Savitzky-Golay window must be odd, got {window_bp}"
        ));
    }
    if window_bp < 5 {
        return Err(format!(
            "order-{SAVGOL_POLYNOMIAL_ORDER} Savitzky-Golay smoothing requires a window of at least 5 bp, got {window_bp}"
        ));
    }

    Ok(MidpointSmoothing::SavGol { window_bp })
}

/// Build center-point Savitzky-Golay smoothing coefficients for polynomial order 3.
///
/// Savitzky-Golay smoothing fits a low-degree polynomial over a symmetric local window and
/// returns the fitted value at the center. For an odd window and polynomial order 3, symmetry
/// makes the center smoothing coefficients identical to the order-2 case.
///
/// The closed-form coefficient for offset `x` in `[-m, m]` is:
///
/// ```text
/// 3 * (3m^2 + 3m - 1 - 5x^2) / ((2m + 1) * (4m^2 + 4m - 3))
/// ```
///
/// These coefficients preserve constants, linear trends, quadratics, and cubics at the window
/// center. That gives small, auditable code without depending on an external Savitzky-Golay
/// package for runtime coefficient generation.
pub(crate) fn order3_coefficients(window_size: u32) -> Result<Vec<f64>> {
    ensure!(
        window_size % 2 == 1,
        "Savitzky-Golay window must be odd, got {window_size}"
    );
    ensure!(
        window_size >= 5,
        "order-3 Savitzky-Golay smoothing requires a window of at least 5 bp, got {window_size}"
    );

    let half_window = window_size / 2;
    let half_window_f64 = f64::from(half_window);
    let denominator = (2.0 * half_window_f64 + 1.0)
        * (4.0 * half_window_f64.powi(2) + 4.0 * half_window_f64 - 3.0);
    let baseline = 3.0 * half_window_f64.powi(2) + 3.0 * half_window_f64 - 1.0;

    let mut coefficients = Vec::with_capacity(window_size as usize);
    for offset in -(half_window as i32)..=(half_window as i32) {
        let offset_f64 = f64::from(offset);
        coefficients.push(3.0 * (baseline - 5.0 * offset_f64.powi(2)) / denominator);
    }

    Ok(coefficients)
}

#[cfg(test)]
mod tests {
    include!("smoothing_tests.rs");
}
