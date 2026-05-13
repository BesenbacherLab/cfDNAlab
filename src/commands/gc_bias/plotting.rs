use crate::commands::gc_bias::binning::BinnedAxis;
use crate::commands::gc_bias::binning::compute_bin_edges;
use crate::commands::gc_bias::load_reference_bias::ReferenceGCMetadata;
use crate::shared::io::dot_join;
use crate::shared::plotters::{
    heatmap::{HeatmapFormat, HeatmapUpsample, write_heatmap_with_histograms},
    histogram::HistogramSpec,
    lineplot::{LinePlotSeries, write_line_plot_png, write_multi_line_plot_png},
};
use anyhow::{Context, Result};
use ndarray::{Array1, Array2, Axis};
use std::path::{Path, PathBuf};

const SELECTED_LENGTH_BIAS_MIN_BP: usize = 80;
const SELECTED_LENGTH_BIAS_MAX_BP: usize = 220;
const SELECTED_LENGTH_BIAS_STEP_BP: usize = 20;

/// Plot GC bias diagnostics: average bias curves and heatmaps.
///
/// Generates line plots for weighted and unweighted average bias across GC bins, plus heatmaps
/// across length and GC bins for quick QC of correction matrices.
///
/// Parameters
/// ----------
/// - `output_dir`:
///     Destination for plot files.
/// - `prefix`:
///     Optional output-file prefix.
/// - `gc_bins`:
///     Mapping of GC-bin indices to contiguous GC ranges, used to build axis labels.
/// - `length_bins`:
///     Mapping of length-bin indices to contiguous fragment length ranges for axis labels.
/// - `correction_matrix`:
///     Correction factors shaped `(length_bin, gc_bin)`.
/// - `length_bin_frequencies`:
///     Weight for each length bin when computing weighted averages.
/// - `reference_metadata`:
///     Metadata describing fragment length bounds for labeling.
/// - `full_res_counts`:
///     Count grid on GC percent x fragment length used for unbinned histograms.
/// - `binned_counts`:
///     Count grid collapsed to GC/length bins for the binned heatmap histograms.
pub fn plot_gc_bias(
    output_dir: &Path,
    prefix: &str,
    gc_bins: &BinnedAxis,
    length_bins: &BinnedAxis,
    correction_matrix: &Array2<f64>,
    length_bin_frequencies: &Array1<f64>,
    reference_metadata: &ReferenceGCMetadata,
    full_res_counts: &Array2<f64>,
    binned_counts: &Array2<f64>,
) -> Result<Vec<PathBuf>> {
    let mut written_paths = Vec::new();

    let gc_edges = compute_bin_edges(gc_bins, 0, 100)?;
    let x_values: Vec<f64> = gc_edges
        .windows(2)
        .map(|window| {
            let start = window[0] as f64;
            let end = window[1] as f64;
            (start + end) / 2.0
        })
        .collect();
    let gc_edges_f: Vec<f64> = gc_edges.iter().map(|v| *v as f64).collect();
    let length_edges = compute_bin_edges(
        length_bins,
        reference_metadata.min_fragment_length as u32,
        reference_metadata.max_fragment_length as u32,
    )?;
    let length_edges_f: Vec<f64> = length_edges.iter().map(|v| *v as f64).collect();

    // Histograms tied to the full-resolution heatmap (GC% x fragment length)
    let gc_histogram = HistogramSpec::from_binned(
        (0..=full_res_counts.ncols()).map(|v| v as f64).collect(),
        full_res_counts.sum_axis(Axis(0)).to_vec(),
    )?;
    let length_histogram = HistogramSpec::from_binned(
        (0..=full_res_counts.nrows())
            .map(|i| reference_metadata.min_fragment_length as f64 + i as f64)
            .collect(),
        full_res_counts.sum_axis(Axis(1)).to_vec(),
    )?;

    // Histograms aligned to the binned heatmap (bin-index space)
    let gc_bin_histogram = HistogramSpec::from_binned(
        (0..=binned_counts.ncols())
            .map(|i| i as f64 - 0.5)
            .collect(),
        binned_counts.sum_axis(Axis(0)).to_vec(),
    )?;
    let length_bin_histogram = HistogramSpec::from_binned(
        (0..=binned_counts.nrows())
            .map(|i| i as f64 - 0.5)
            .collect(),
        binned_counts.sum_axis(Axis(1)).to_vec(),
    )?;

    // Bias matrix and per-GC averages:
    // - Convert corrections to bias (1/cf, keep zeros masked)
    // - Unweighted: simple mean across length bins ignoring zeroed cells
    // - Weighted: length-frequency weighted mean so common lengths influence more
    let num_gc_bins = correction_matrix.ncols();
    let bias_matrix = correction_matrix.mapv(|cf| if cf == 0.0 { 0.0 } else { 1.0 / cf });
    let mut unweighted_bias = vec![0.0; num_gc_bins];
    let mut unweighted_counts = vec![0usize; num_gc_bins];

    for length_biases in bias_matrix.outer_iter() {
        for (gc_idx, &bias) in length_biases.iter().enumerate() {
            if bias == 0.0 {
                continue;
            }
            unweighted_bias[gc_idx] += bias;
            unweighted_counts[gc_idx] += 1;
        }
    }

    for (bias, count) in unweighted_bias.iter_mut().zip(unweighted_counts.iter()) {
        if *count > 0 {
            *bias /= *count as f64;
        }
    }

    let mut weighted_bias = vec![0.0; num_gc_bins];
    let mut weight_per_gc = vec![0.0; num_gc_bins];

    for (length_biases, &length_weight) in bias_matrix.outer_iter().zip(length_bin_frequencies) {
        if length_weight == 0.0 {
            continue;
        }
        for (gc_idx, &bias) in length_biases.iter().enumerate() {
            if bias == 0.0 {
                continue;
            }
            weight_per_gc[gc_idx] += length_weight;
            weighted_bias[gc_idx] += length_weight * bias;
        }
    }

    for (bias, weight) in weighted_bias.iter_mut().zip(weight_per_gc.iter()) {
        if *weight > 0.0 {
            *bias /= *weight;
        }
    }

    // Line plots: average GC bias across lengths (unweighted and weighted)
    let plot_path_unweighted = output_dir.join(dot_join(&[
        prefix,
        "avg_gc_bias_across_lengths_unweighted.png",
    ]));
    write_line_plot_png(
        &plot_path_unweighted,
        "Average GC bias across fragment lengths (unweighted)",
        "GC bin (%)",
        "GC bias",
        &x_values,
        &unweighted_bias,
        1600,
        1000,
    )
    .with_context(|| format!("writing GC bias plot to {}", plot_path_unweighted.display()))?;
    written_paths.push(plot_path_unweighted);

    let plot_path_weighted = output_dir.join(dot_join(&[
        prefix,
        "avg_gc_bias_across_lengths_weighted.png",
    ]));
    write_line_plot_png(
        &plot_path_weighted,
        "Average GC bias across fragment lengths (weighted by length frequency)",
        "GC bin (%)",
        "GC bias",
        &x_values,
        &weighted_bias,
        1600,
        1000,
    )
    .with_context(|| format!("writing GC bias plot to {}", plot_path_weighted.display()))?;
    written_paths.push(plot_path_weighted);

    let selected_lengths_plot_path = output_dir.join(dot_join(&[
        prefix,
        "gc_bias_by_selected_lengths_80_220bp.png",
    ]));
    if write_selected_length_bias_plot(
        &selected_lengths_plot_path,
        &x_values,
        &bias_matrix,
        length_bins,
        reference_metadata,
    )
    .with_context(|| {
        format!(
            "writing selected-length GC bias plot to {}",
            selected_lengths_plot_path.display()
        )
    })? {
        written_paths.push(selected_lengths_plot_path);
    }

    let hm_width: u32 = 1000;
    let hm_height: u32 = 700;
    let scaling_factor = (hm_height as f32 / bias_matrix.nrows() as f32)
        .max(hm_width as f32 / bias_matrix.ncols() as f32)
        .ceil() as usize;

    // Full-resolution heatmap with GC% / length histograms
    let heatmap_path = output_dir.join(dot_join(&[prefix, "gc_bias_heatmap.png"]));
    write_heatmap_with_histograms(
        &heatmap_path,
        "GC bias per length and GC %",
        "GC (%)",
        "Fragment length (bp)",
        &bias_matrix,
        Some(&gc_edges_f),
        Some(&length_edges_f),
        Some(&gc_histogram),
        Some(&length_histogram),
        None,
        None,
        Some(1.0),
        None,
        None,
        true,
        scaling_factor,
        HeatmapUpsample::Nearest,
        hm_width,
        hm_height,
        HeatmapFormat::Png,
    )
    .with_context(|| format!("writing GC bias heatmap to {}", heatmap_path.display()))?;
    written_paths.push(heatmap_path);

    // Binned heatmap with bin-index histograms
    let heatmap_path = output_dir.join(dot_join(&[prefix, "gc_bias_heatmap.bins.png"]));
    write_heatmap_with_histograms(
        &heatmap_path,
        "GC bias per length bin and GC bin",
        "GC bin",
        "Fragment length bin",
        &bias_matrix,
        None,
        None,
        Some(&gc_bin_histogram),
        Some(&length_bin_histogram),
        None,
        None,
        Some(1.0),
        None,
        None,
        true,
        scaling_factor,
        HeatmapUpsample::Nearest,
        hm_width,
        hm_height,
        HeatmapFormat::Png,
    )
    .with_context(|| {
        format!(
            "writing GC bias heatmap (bins) with histograms to {}",
            heatmap_path.display()
        )
    })?;
    written_paths.push(heatmap_path);

    Ok(written_paths)
}

fn write_selected_length_bias_plot(
    output_path: &Path,
    gc_bin_midpoints_pct: &[f64],
    bias_matrix: &Array2<f64>,
    length_bins: &BinnedAxis,
    reference_metadata: &ReferenceGCMetadata,
) -> Result<bool> {
    let selected_lengths = selected_length_biases(bias_matrix, length_bins, reference_metadata);
    if selected_lengths.is_empty() {
        return Ok(false);
    }

    let gc_content: Vec<f64> = gc_bin_midpoints_pct
        .iter()
        .map(|gc_percent| gc_percent / 100.0)
        .collect();
    let labels: Vec<String> = selected_lengths
        .iter()
        .map(|(fragment_length, _)| format!("{fragment_length} bp"))
        .collect();
    let series: Vec<LinePlotSeries<'_>> = selected_lengths
        .iter()
        .zip(labels.iter())
        .map(|((_, bias_values), label)| LinePlotSeries {
            label,
            x_values: &gc_content,
            y_values: bias_values,
        })
        .collect();

    write_multi_line_plot_png(
        output_path,
        "GC bias by selected fragment length",
        "GC content",
        "GC bias",
        &series,
        1600,
        1000,
    )?;
    Ok(true)
}

fn selected_length_biases(
    bias_matrix: &Array2<f64>,
    length_bins: &BinnedAxis,
    reference_metadata: &ReferenceGCMetadata,
) -> Vec<(usize, Vec<f64>)> {
    let mut selected_lengths = Vec::new();
    for fragment_length in (SELECTED_LENGTH_BIAS_MIN_BP..=SELECTED_LENGTH_BIAS_MAX_BP)
        .step_by(SELECTED_LENGTH_BIAS_STEP_BP)
    {
        if fragment_length < reference_metadata.min_fragment_length
            || fragment_length > reference_metadata.max_fragment_length
        {
            continue;
        }

        let length_index = fragment_length - reference_metadata.min_fragment_length;
        let Some(&length_bin_index) = length_bins.index_to_bin.get(&length_index) else {
            continue;
        };
        if length_bin_index >= bias_matrix.nrows() {
            continue;
        }

        selected_lengths.push((fragment_length, bias_matrix.row(length_bin_index).to_vec()));
    }

    selected_lengths
}
