use anyhow::{Result, bail};
use rust_htslib::bam::record::{Aux, Record};
use std::sync::atomic::{AtomicUsize, Ordering};
use tracing::warn;

/// GC weight extracted from an AUX tag together with validity status.
#[derive(Debug, Clone, Copy, Default)]
pub struct GCTagValue {
    /// Parsed GC weight if present and valid.
    pub weight: Option<f32>,
    /// True when the tag was not present on the read.
    pub was_missing: bool,
    /// True when the tag was present but unusable (wrong type, NaN, or outside the supported
    /// positive range after zero-snapping).
    pub had_invalid: bool,
    /// True when the tag value was outside the supported positive range.
    pub was_out_of_range: bool,
}

/// Values at or below this threshold are treated as exact zero.
pub const ZEROISH_GC_WEIGHT_TOLERANCE: f32 = 2.0 * f32::EPSILON;

/// Smallest supported positive GC weight.
pub const MIN_REASONABLE_GC_WEIGHT: f32 = 1.0e-3;

/// Reject GC weights that are clearly out of range to avoid runaway coverage when tags are corrupt.
///
/// Correction weights are expected to hover around 0–a few hundred at most.
/// Values far beyond that are treated as invalid.
pub const MAX_REASONABLE_GC_WEIGHT: f32 = 1.0e3;
const MAX_GC_WARNINGS: usize = 5;
static EXTREME_GC_WARNINGS: AtomicUsize = AtomicUsize::new(0);

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SanitizedGCWeight {
    Usable(f64),
    Unusable { out_of_range: bool },
}

/// Explicit classification of fragment-level GC-tag state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ClassifiedGCTagWeight {
    /// Usable GC weight to apply to the fragment.
    Usable(f32),
    /// The requested GC tag was not present.
    Missing,
    /// The tag was present but unusable.
    Invalid { out_of_range: bool },
}

#[inline]
fn warn_extreme_gc_weight(v: f32) {
    let seen = EXTREME_GC_WARNINGS.fetch_add(1, Ordering::Relaxed);
    if seen < MAX_GC_WARNINGS {
        warn!(
            target: "gc-tag",
            "warning: GC tag weight {:.3e} is outside the supported positive range [{:.0e}, {:.0e}] (zero is valid); treating as invalid",
            v, MIN_REASONABLE_GC_WEIGHT, MAX_REASONABLE_GC_WEIGHT
        );
        if seen + 1 == MAX_GC_WARNINGS {
            warn!(target: "gc-tag", "warning: suppressing further GC tag weight warnings");
        }
    }
}

#[inline]
fn is_gc_weight_out_of_range(v: f32) -> bool {
    (!v.is_finite() && !v.is_nan())
        || (v < -ZEROISH_GC_WEIGHT_TOLERANCE)
        || (v > ZEROISH_GC_WEIGHT_TOLERANCE
            && (v < MIN_REASONABLE_GC_WEIGHT || v > MAX_REASONABLE_GC_WEIGHT))
}

#[inline]
fn is_zeroish_gc_weight(weight: Option<f32>) -> bool {
    matches!(weight, Some(value) if value.abs() <= ZEROISH_GC_WEIGHT_TOLERANCE)
}

/// Describe how unusable GC information is handled for user-facing logs.
pub fn gc_failure_action_description(neutralize_invalid_gc: bool) -> &'static str {
    if neutralize_invalid_gc {
        "fragment counted with weight 1.0"
    } else {
        "fragment skipped"
    }
}

/// Sanitize a GC weight from either a BAM aux tag or a GC correction package.
///
/// Rules:
/// - `NaN` -> unusable, not counted as out-of-range
/// - `abs(v) <= 2 * f32::EPSILON` -> exact zero
/// - negative values below `-2 * f32::EPSILON` -> unusable
/// - positive values in `(2 * eps, 1e-3)` -> unusable
/// - positive values in `[1e-3, 1e3]` -> usable
/// - positive values above `1e3` and infinities -> unusable
pub fn sanitize_gc_weight(v: f64) -> SanitizedGCWeight {
    if v.is_nan() {
        return SanitizedGCWeight::Unusable {
            out_of_range: false,
        };
    }
    if !v.is_finite() {
        return SanitizedGCWeight::Unusable { out_of_range: true };
    }
    if v.abs() <= ZEROISH_GC_WEIGHT_TOLERANCE as f64 {
        return SanitizedGCWeight::Usable(0.0);
    }
    if v < -(ZEROISH_GC_WEIGHT_TOLERANCE as f64) {
        return SanitizedGCWeight::Unusable { out_of_range: true };
    }
    if v < MIN_REASONABLE_GC_WEIGHT as f64 || v > MAX_REASONABLE_GC_WEIGHT as f64 {
        return SanitizedGCWeight::Unusable { out_of_range: true };
    }
    SanitizedGCWeight::Usable(v)
}

impl GCTagValue {
    #[inline]
    pub fn missing() -> Self {
        GCTagValue {
            weight: None,
            was_missing: true,
            had_invalid: false,
            was_out_of_range: false,
        }
    }

    #[inline]
    fn from_number(v: f32) -> Self {
        let sanitized = sanitize_gc_weight(v as f64);
        if is_gc_weight_out_of_range(v) {
            warn_extreme_gc_weight(v);
        }
        match sanitized {
            SanitizedGCWeight::Usable(weight) => GCTagValue {
                weight: Some(weight as f32),
                was_missing: false,
                had_invalid: false,
                was_out_of_range: false,
            },
            SanitizedGCWeight::Unusable { out_of_range } => GCTagValue {
                weight: None,
                was_missing: false,
                had_invalid: true,
                was_out_of_range: out_of_range,
            },
        }
    }

    /// Classify the fragment-level GC-tag state into one explicit branch.
    ///
    /// This keeps command code from having to manually decode the state from
    /// multiple flags and ensures impossible combinations fail loudly.
    pub fn classify(self) -> Result<ClassifiedGCTagWeight> {
        match (self.weight, self.was_missing, self.had_invalid) {
            (Some(weight), false, false) => Ok(ClassifiedGCTagWeight::Usable(weight)),
            (None, true, false) => Ok(ClassifiedGCTagWeight::Missing),
            (None, false, true) => Ok(ClassifiedGCTagWeight::Invalid {
                out_of_range: self.was_out_of_range,
            }),
            _ => bail!(
                "inconsistent GC tag state: weight={:?}, was_missing={}, had_invalid={}, was_out_of_range={}",
                self.weight,
                self.was_missing,
                self.had_invalid,
                self.was_out_of_range
            ),
        }
    }
}

/// Read a numeric GC weight from an AUX tag on a record.
///
/// Accepts integer and floating-point tag types. Returns `weight=None` when the
/// tag is missing. Marks `had_invalid=true` when the tag is present but has an
/// unsupported type, NaN, or a value outside the allowed range.
pub fn read_gc_tag_from_record(rec: &Record, tag: &[u8]) -> GCTagValue {
    match rec.aux(tag) {
        Ok(Aux::Float(v)) => GCTagValue::from_number(v),
        Ok(Aux::Double(v)) => GCTagValue::from_number(v as f32),
        Ok(Aux::U8(v)) => GCTagValue::from_number(v as f32),
        Ok(Aux::U16(v)) => GCTagValue::from_number(v as f32),
        Ok(Aux::U32(v)) => GCTagValue::from_number(v as f32),
        Ok(Aux::I8(v)) => GCTagValue::from_number(v as f32),
        Ok(Aux::I16(v)) => GCTagValue::from_number(v as f32),
        Ok(Aux::I32(v)) => GCTagValue::from_number(v as f32),
        Ok(_) => GCTagValue {
            weight: None,
            was_missing: false,
            had_invalid: true,
            was_out_of_range: false,
        },
        Err(_) => GCTagValue::missing(),
    }
}

/// Combine per-read GC weights into a fragment-level weight.
///
/// Rules:
/// - Zero on either mate -> fragment weight 0.0.
/// - Invalid tag on either mate -> invalid fragment weight.
/// - One usable mate plus one missing mate -> use the usable mate weight.
/// - Missing tag on both mates -> missing fragment weight.
/// - Otherwise average the two mate weights, then validate the fragment-level result.
pub fn combine_gc_tag_values(a: &GCTagValue, b: &GCTagValue) -> GCTagValue {
    if is_zeroish_gc_weight(a.weight) || is_zeroish_gc_weight(b.weight) {
        return GCTagValue {
            weight: Some(0.0),
            was_missing: false,
            had_invalid: false,
            was_out_of_range: false,
        };
    }
    if a.had_invalid || b.had_invalid {
        return GCTagValue {
            weight: None,
            was_missing: false,
            had_invalid: true,
            was_out_of_range: a.was_out_of_range || b.was_out_of_range,
        };
    }

    // From here both inputs are present and not explicitly invalid.
    match (a.weight, b.weight) {
        (Some(weight), None) | (None, Some(weight)) => GCTagValue {
            weight: Some(weight),
            was_missing: false,
            had_invalid: false,
            was_out_of_range: false,
        },
        (Some(wa), Some(wb)) => match sanitize_gc_weight(((wa as f64) + (wb as f64)) / 2.0) {
            SanitizedGCWeight::Usable(weight) => GCTagValue {
                weight: Some(weight as f32),
                was_missing: false,
                had_invalid: false,
                was_out_of_range: false,
            },
            SanitizedGCWeight::Unusable { out_of_range } => GCTagValue {
                weight: None,
                was_missing: false,
                had_invalid: true,
                was_out_of_range: out_of_range,
            },
        },
        (None, None) => GCTagValue::missing(),
    }
}

#[cfg(test)]
mod tests {
    include!("gc_tag_tests.rs");
}
