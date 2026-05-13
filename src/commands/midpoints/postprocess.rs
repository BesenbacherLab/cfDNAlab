use crate::commands::midpoints::smoothing::{MidpointSmoothing, order3_coefficients};
use anyhow::{Context, Result, bail, ensure};
use ndarray::{Array3, ArrayView3, s};

/// Resolved layout for counted and written midpoint profiles.
///
/// Input BED intervals define `output_len`, the public number of positions. Smoothing may require
/// extra counted positions on both sides. Those positions are represented by `flanked_length` and
/// are removed before writing the final output. Binning always applies after smoothing and trimming.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct ProfileLayout {
    pub output_len: usize,
    pub flanked_length: usize,
    pub bin_size: u32,
    pub smoothing_window: Option<u32>,
    pub smoothing_flank: u32,
    pub output_positions: usize,
}

impl ProfileLayout {
    /// Resolve public output dimensions into the internal counted profile shape.
    ///
    /// This is the central validation point for smoothing and final binning. It keeps the output
    /// interval fixed, derives any computation-only smoothing flank, and fails early when the
    /// requested Savitzky-Golay window cannot fit the public interval.
    pub(crate) fn resolve(
        output_len: usize,
        bin_size: u32,
        smoothing: MidpointSmoothing,
    ) -> Result<Self> {
        ensure!(
            output_len > 0,
            "midpoint output interval length must be > 0"
        );
        ensure!(bin_size >= 1, "--bin-size must be at least 1");

        let smoothing_window = match smoothing {
            MidpointSmoothing::None => None,
            MidpointSmoothing::SavGol { window_bp } => Some(window_bp),
        };

        if let Some(window) = smoothing_window {
            if output_len < 7 {
                bail!(
                    "Cannot apply Savitzky-Golay smoothing to {output_len} bp output intervals. Use intervals of at least 7 bp or set --smoothing none."
                );
            }
            ensure!(
                window % 2 == 1,
                "Savitzky-Golay window must be odd, got {window}"
            );
            ensure!(
                window >= 5,
                "order-3 Savitzky-Golay smoothing requires a window of at least 5 bp, got {window}"
            );
            if window as usize > output_len {
                let suggested_window = largest_odd_window_that_fits(output_len)
                    .context("could not suggest a Savitzky-Golay window")?;
                bail!(
                    "Savitzky-Golay window {window} bp is longer than the {output_len} bp output interval. Use --smoothing savgol={suggested_window} or longer intervals."
                );
            }
        }

        let smoothing_flank = smoothing_window.map_or(0, |window| window / 2);
        let flanked_length = output_len
            .checked_add((smoothing_flank as usize).saturating_mul(2))
            .context("midpoint count interval length overflow")?;
        let output_positions = output_len.div_ceil(bin_size as usize);

        Ok(Self {
            output_len,
            flanked_length,
            bin_size,
            smoothing_window,
            smoothing_flank,
            output_positions,
        })
    }

    #[inline]
    pub(crate) fn is_identity(&self) -> bool {
        self.smoothing_window.is_none()
            && self.bin_size == 1
            && self.flanked_length == self.output_len
    }
}

fn largest_odd_window_that_fits(output_len: usize) -> Option<u32> {
    if output_len < 7 {
        return None;
    }
    let window = if output_len % 2 == 1 {
        output_len
    } else {
        output_len - 1
    };
    u32::try_from(window).ok()
}

/// Apply final smoothing, flank trimming, and binning to a merged midpoint tensor.
///
/// Counting can use computation-only flank positions so smoothing has real data at the output
/// edges. This step removes those flanks, optionally smooths each `(group, length_bin)` profile
/// along the position axis, and averages final bins when requested.
///
/// Returns `None` for the identity path so the caller can write the merged dense tensor without
/// allocating a second array for unsmoothed, unbinned profiles.
pub(crate) fn postprocess_profile(
    profile: ArrayView3<'_, f32>,
    layout: ProfileLayout,
) -> Result<Option<Array3<f32>>> {
    let (_, _, positions) = profile.dim();
    ensure!(
        positions == layout.flanked_length,
        "profile postprocessing expected {} counted positions, got {}",
        layout.flanked_length,
        positions
    );

    if layout.is_identity() {
        return Ok(None);
    }

    let smoothed_profile = if let Some(window) = layout.smoothing_window {
        Some(smooth_trimmed_profile(profile, window, layout.output_len)?)
    } else {
        ensure!(
            positions == layout.output_len,
            "unsmoothed midpoint profile expected {} positions, got {}",
            layout.output_len,
            positions
        );
        None
    };

    if layout.bin_size == 1 {
        return Ok(smoothed_profile);
    }

    let binned_profile = {
        let profile_to_bin = match smoothed_profile.as_ref() {
            Some(smoothed) => smoothed.view(),
            None => profile.view(),
        };
        bin_profile(profile_to_bin, layout.bin_size)
    };
    Ok(Some(binned_profile))
}

/// Smooth the expanded tensor and return only the original output positions.
///
/// No boundary mode is needed here. The counting step already expanded each interval by the
/// Savitzky-Golay support radius, so every retained output base has a complete filter window.
fn smooth_trimmed_profile(
    profile: ArrayView3<'_, f32>,
    window: u32,
    output_len: usize,
) -> Result<Array3<f32>> {
    let coefficients = order3_coefficients(window)?;
    let (num_groups, num_length_bins, positions) = profile.dim();
    let window_len = window as usize;
    let smoothing_flank = window_len / 2;
    let expected_positions = output_len
        .checked_add(smoothing_flank.saturating_mul(2))
        .context("smoothed profile position count overflow")?;
    ensure!(
        positions == expected_positions,
        "Savitzky-Golay smoothing expected {expected_positions} counted positions, got {positions}"
    );

    let mut out = Array3::<f32>::zeros((num_groups, num_length_bins, output_len));
    for group_idx in 0..num_groups {
        for length_idx in 0..num_length_bins {
            let input = profile.slice(s![group_idx, length_idx, ..]);
            for output_pos in 0..output_len {
                let mut value = 0.0_f64;
                // `output_pos` is the public retained coordinate, but the counted tensor starts
                // `smoothing_flank` bases earlier. Therefore `input[output_pos]` is the left edge
                // of the complete filter support, `input[output_pos + smoothing_flank]` is the
                // retained base itself, and the last coefficient sees the right flank.
                for (offset, coefficient) in coefficients.iter().enumerate() {
                    value += f64::from(input[output_pos + offset]) * coefficient;
                }
                out[[group_idx, length_idx, output_pos]] = value as f32;
            }
        }
    }

    Ok(out)
}

/// Average adjacent final profile positions.
///
/// The last bin can be shorter than `bin_size`, so it is divided by its real width. This keeps the
/// output on an average-count scale and makes sum-style downstream use a simple multiplication.
fn bin_profile(profile: ArrayView3<'_, f32>, bin_size: u32) -> Array3<f32> {
    let bin_width = bin_size as usize;
    let (num_groups, num_length_bins, positions) = profile.dim();
    let output_positions = positions.div_ceil(bin_width);
    let mut out = Array3::<f32>::zeros((num_groups, num_length_bins, output_positions));

    for group_idx in 0..num_groups {
        for length_idx in 0..num_length_bins {
            for output_pos in 0..output_positions {
                let start = output_pos * bin_width;
                let end = (start + bin_width).min(positions);
                let width = (end - start) as f32;
                let mut sum = 0.0_f32;
                for input_pos in start..end {
                    sum += profile[[group_idx, length_idx, input_pos]];
                }
                out[[group_idx, length_idx, output_pos]] = sum / width;
            }
        }
    }

    out
}

#[cfg(test)]
mod tests {
    include!("postprocess_tests.rs");
}
