use rust_htslib::bam::record::{Aux, Record};
use std::sync::atomic::{AtomicUsize, Ordering};

/// GC weight extracted from an AUX tag together with validity status.
#[derive(Debug, Clone, Copy, Default)]
pub struct GcTagValue {
    /// Parsed GC weight if present and valid.
    pub weight: Option<f32>,
    /// True when the tag was present but invalid (wrong type, NaN, negative).
    pub had_invalid: bool,
    /// True when the tag value was non-NaN but outside the allowed numeric range (includes ±inf).
    pub was_out_of_range: bool,
}

/// Reject GC weights that are clearly out of range to avoid runaway coverage when tags are corrupt.
///
/// Correction weights are expected to hover around 0–a few hundred at most.
/// Values far beyond that are treated as invalid.
pub const MAX_REASONABLE_GC_WEIGHT: f32 = 1.0e3;
const MAX_GC_WARNINGS: usize = 5;
static EXTREME_GC_WARNINGS: AtomicUsize = AtomicUsize::new(0);

#[inline]
fn warn_extreme_gc_weight(v: f32) {
    let seen = EXTREME_GC_WARNINGS.fetch_add(1, Ordering::Relaxed);
    if seen < MAX_GC_WARNINGS {
        eprintln!(
            "warning: GC tag weight {:.3e} is outside [0, {:.0}]; treating as invalid",
            v, MAX_REASONABLE_GC_WEIGHT
        );
        if seen + 1 == MAX_GC_WARNINGS {
            eprintln!("warning: suppressing further GC tag weight warnings");
        }
    }
}

impl GcTagValue {
    #[inline]
    fn from_number(v: f32) -> Self {
        if v.is_nan() {
            // Invalid but not counted as out-of-range
            return GcTagValue {
                weight: None,
                had_invalid: true,
                was_out_of_range: false,
            };
        }

        if v.is_finite() && (0.0..=MAX_REASONABLE_GC_WEIGHT).contains(&v) {
            GcTagValue {
                weight: Some(v),
                had_invalid: false,
                was_out_of_range: false,
            }
        } else {
            warn_extreme_gc_weight(v);
            GcTagValue {
                weight: None,
                had_invalid: true,
                was_out_of_range: true,
            }
        }
    }
}

/// Read a numeric GC weight from an AUX tag on a record.
///
/// Accepts integer and floating-point tag types. Returns `weight=None` when the
/// tag is missing. Marks `had_invalid=true` when the tag is present but has an
/// unsupported type, NaN, or a value outside the allowed range.
pub fn read_gc_tag_from_record(rec: &Record, tag: &[u8]) -> GcTagValue {
    match rec.aux(tag) {
        Ok(Aux::Float(v)) => GcTagValue::from_number(v),
        Ok(Aux::Double(v)) => GcTagValue::from_number(v as f32),
        Ok(Aux::U8(v)) => GcTagValue::from_number(v as f32),
        Ok(Aux::U16(v)) => GcTagValue::from_number(v as f32),
        Ok(Aux::U32(v)) => GcTagValue::from_number(v as f32),
        Ok(Aux::I8(v)) => GcTagValue::from_number(v as f32),
        Ok(Aux::I16(v)) => GcTagValue::from_number(v as f32),
        Ok(Aux::I32(v)) => GcTagValue::from_number(v as f32),
        Ok(_) => GcTagValue {
            weight: None,
            had_invalid: true,
            was_out_of_range: false,
        },
        Err(_) => GcTagValue::default(), // Tag missing
    }
}

/// Combine per-read GC weights into a fragment-level weight.
///
/// Rules:
/// - Invalid tag on either mate -> invalid fragment weight.
/// - Zero on either mate -> fragment weight 0.0.
/// - Both present -> average them.
/// - One present -> use it.
/// - Both missing -> no weight.
pub fn combine_gc_tag_values(a: &GcTagValue, b: &GcTagValue) -> GcTagValue {
    if a.had_invalid || b.had_invalid {
        return GcTagValue {
            weight: None,
            had_invalid: true,
            was_out_of_range: a.was_out_of_range || b.was_out_of_range,
        };
    }

    // From here both inputs are valid and within range.
    match (a.weight, b.weight) {
        (Some(wa), Some(wb)) => {
            if wa == 0.0 || wb == 0.0 {
                GcTagValue {
                    weight: Some(0.0),
                    had_invalid: false,
                    was_out_of_range: false,
                }
            } else {
                GcTagValue {
                    weight: Some((wa + wb) / 2.0),
                    had_invalid: false,
                    was_out_of_range: false,
                }
            }
        }
        (Some(w), None) | (None, Some(w)) => GcTagValue {
            weight: Some(w),
            had_invalid: false,
            was_out_of_range: false,
        },
        (None, None) => GcTagValue::default(),
    }
}
