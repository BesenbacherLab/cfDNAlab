use rust_htslib::bam::record::{Aux, Record};

/// GC weight extracted from an AUX tag together with validity status.
#[derive(Debug, Clone, Copy, Default)]
pub struct GcTagValue {
    /// Parsed GC weight if present and valid.
    pub weight: Option<f32>,
    /// True when the tag was present but invalid (wrong type, NaN, negative).
    pub had_invalid: bool,
}

impl GcTagValue {
    #[inline]
    fn from_number(v: f32) -> Self {
        if v.is_finite() && v >= 0.0 {
            GcTagValue {
                weight: Some(v),
                had_invalid: false,
            }
        } else {
            GcTagValue {
                weight: None,
                had_invalid: true,
            }
        }
    }
}

/// Read a numeric GC weight from an AUX tag on a record.
///
/// Accepts integer and floating-point tag types. Returns `weight=None` when the
/// tag is missing. Marks `had_invalid=true` when the tag is present but has an
/// unsupported type or a non-finite/negative value.
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
        };
    }

    match (a.weight, b.weight) {
        (Some(wa), Some(wb)) => {
            if wa == 0.0 || wb == 0.0 {
                GcTagValue {
                    weight: Some(0.0),
                    had_invalid: false,
                }
            } else {
                GcTagValue {
                    weight: Some((wa + wb) / 2.0),
                    had_invalid: false,
                }
            }
        }
        (Some(w), None) | (None, Some(w)) => GcTagValue {
            weight: Some(w),
            had_invalid: false,
        },
        (None, None) => GcTagValue::default(),
    }
}
