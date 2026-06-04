use std::path::Path;

use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;
use ndarray::{Array2, ArrayView3, Axis, s};

use crate::shared::io::dot_join;
use crate::shared::plotters::{
    heatmap::{HeatmapFormat, HeatmapUpsample, write_heatmap},
    lineplot::write_line_plot_png,
};

/// Plot midpoint profiles and optional heatmaps for the requested groups.
///
/// Uses the merged midpoint counts to emit per-group line plots and length-binned heatmaps.
/// Assumes counts are shaped `(group, length_bin, position)`.
///
/// Parameters
/// ----------
/// - `prefix`:
///     File prefix for plot outputs.
/// - `output_dir`:
///     Directory where plots are written.
/// - `plot_groups`:
///     Group indices to plot.
/// - `length_bins`:
///     Length bin edges matching the counts.
/// - `group_idx_to_name`:
///     Mapping from group index to user-readable name.
/// - `counts`:
///     Midpoint counts shaped `(group, length_bin, position)`.
pub(crate) fn plot_midpoint_profiles(
    prefix: &str,
    output_dir: &Path,
    plot_groups: &[usize],
    length_bins: &[u32],
    group_idx_to_name: &FxHashMap<u64, String>,
    counts: ArrayView3<'_, f32>,
) -> Result<()> {
    if plot_groups.is_empty() {
        return Ok(());
    }

    let (num_groups, num_length_bins, window_size) = counts.dim();
    let y_edges: Vec<f64> = length_bins.iter().map(|&b| b as f64).collect();

    if let Some(bad_idx) = plot_groups
        .iter()
        .copied()
        .find(|&idx| !group_idx_to_name.contains_key(&(idx as u64)))
    {
        bail!(
            "plotting: group index {} does not exist. There are {} groups (0-based).",
            bad_idx,
            num_groups
        );
    }

    if window_size == 0 || num_groups == 0 {
        return Ok(());
    }

    let x_values: Vec<f64> = if window_size % 2 == 1 {
        let center = (window_size / 2) as i64;
        (0..window_size)
            .map(|idx| (idx as i64 - center) as f64)
            .collect()
    } else {
        (0..window_size).map(|idx| idx as f64).collect()
    };
    let x_edges: Vec<f64> = if window_size % 2 == 1 {
        let center = (window_size / 2) as f64;
        (0..=window_size)
            .map(|idx| idx as f64 - center - 0.5)
            .collect()
    } else {
        (0..=window_size).map(|idx| idx as f64 - 0.5).collect()
    };

    for &group_idx in plot_groups {
        let profile: Vec<f64> = counts
            .slice(s![group_idx, .., ..])
            .sum_axis(Axis(0))
            .iter()
            .map(|&v| v as f64)
            .collect();

        if profile.is_empty() {
            continue;
        }

        let plot_path = output_dir.join(dot_join(&[
            prefix,
            &format!("midpoint_profile.group_{}.png", group_idx),
        ]));
        let title = group_idx_to_name
            .get(&(group_idx as u64))
            .map(|name| format!("Midpoint profile ({})", name))
            .unwrap_or_else(|| format!("Midpoint profile (group {})", group_idx));
        write_line_plot_png(
            &plot_path, &title, "Position", "Count", &x_values, &profile, 1600, 900,
        )
        .with_context(|| format!("writing midpoint profile plot to {}", plot_path.display()))?;

        if num_length_bins > 1 {
            let mut heatmap_values: Array2<f64> = Array2::zeros((num_length_bins, window_size));
            heatmap_values.assign(&counts.slice(s![group_idx, .., ..]).mapv(f64::from));

            let title = group_idx_to_name
                .get(&(group_idx as u64))
                .map(|name| format!("Midpoint profile by length bin ({})", name))
                .unwrap_or_else(|| format!("Midpoint profile by length bin (group {})", group_idx));
            let heatmap_path = output_dir.join(dot_join(&[
                prefix,
                &format!("midpoint_profile.group_{}.heatmap.png", group_idx),
            ]));
            write_heatmap(
                &heatmap_path,
                &title,
                "Position",
                "Fragment length (bp)",
                &heatmap_values,
                Some(&x_edges),
                Some(&y_edges),
                Some(0f64),
                None,
                None,
                None,
                None,
                false,
                1,
                HeatmapUpsample::Nearest,
                1600,
                900,
                HeatmapFormat::Png,
            )
            .with_context(|| format!("writing midpoint heatmap to {}", heatmap_path.display()))?;
        }
    }

    Ok(())
}
