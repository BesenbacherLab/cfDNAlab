use std::fmt;

/// Rounds a floating-point value to a fixed number of decimal places.
///
/// When the requested precision is non-positive the value is rounded to the nearest integer.
///
/// # Parameters
/// - `x`: Value to round.
/// - `decimals`: Number of decimal places to preserve.
///
/// # Returns
/// The rounded value.
pub(crate) fn round_to(x: f64, decimals: i32) -> f64 {
    if decimals <= 0 {
        return x.round();
    }
    let f = 10f64.powi(decimals);
    (x * f).round() / f
}

/// Rounds using a precomputed scaling factor for repeated operations.
///
/// This variant avoids recomputing the power-of-ten factor when many values share the same
/// precision requirement.
///
/// # Parameters
/// - `x`: Value to round.
/// - `factor`: Precomputed `10f64.powi(decimals)` for the desired precision.
///
/// # Returns
/// The rounded value.
pub(crate) fn round_to_with_precomputed_factor(x: f64, factor: f64) -> f64 {
    if factor == 1.0 {
        return x.round();
    }
    (x * factor).round() / factor
}

/// Lightweight adapter for printing rounded numbers without heap allocations.
///
/// The stored value and decimal precision are used to format coverage values compactly in hot
/// loops.
pub(crate) struct CompactNumber {
    pub(crate) v: f64,
    pub(crate) decimals: i32,
}

impl fmt::Display for CompactNumber {
    /// Formats the number using the stored decimal precision without allocating.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Normalize negative zero (can apparently happen after rounding tiny negatives)
        // After this, values will be non-zero
        if self.v == 0.0 {
            return f.write_str("0");
        }
        if self.decimals <= 0 {
            // Integer path: Avoid float fmt; Round then print int
            let r = self.v.round();
            // After integer rounding, also normalize -0 just in case
            if r == 0.0 {
                return f.write_str("0");
            }
            // Write directly; no heap
            return write!(f, "{:.0}", r);
        }
        // Stack buffer; Write fixed decimals
        let mut buf = arrayvec::ArrayString::<64>::new();
        // Write with N decimals
        let _ = fmt::write(
            &mut buf,
            format_args!("{:.*}", self.decimals as usize, self.v),
        );
        // Trim zeros and trailing dot
        while buf.as_bytes().last() == Some(&b'0') {
            buf.pop();
        }
        if buf.as_bytes().last() == Some(&b'.') {
            buf.pop();
        }
        f.write_str(buf.as_str())
    }
}
