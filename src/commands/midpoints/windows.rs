use crate::shared::bed::GroupedWindows;
use anyhow::{Context, Result};
use fxhash::FxHashMap;

/// Get window length and ensure it's the same for ALL windows.
pub fn ensure_uniform_window_len(
    windows_by_chr: &FxHashMap<String, GroupedWindows>,
) -> Result<usize> {
    let mut reference_len: Option<usize> = None;

    for (chr, gw) in windows_by_chr {
        for (start, end, _) in &gw.windows {
            let len = end.checked_sub(*start).with_context(|| {
                format!("Invalid window on {chr}: end ({end}) < start ({start})")
            })? as usize;

            match reference_len {
                None => reference_len = Some(len),
                Some(ref_len) if (len) != ref_len => {
                    anyhow::bail!(
                        "Non-uniform window length detected on {chr}: [{start},{end}) has len {}, expected {}",
                        len,
                        ref_len
                    );
                }
                _ => {}
            }
        }
    }

    reference_len.context("No windows found when checking uniform window length")
}
