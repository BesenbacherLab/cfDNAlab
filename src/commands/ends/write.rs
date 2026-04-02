use crate::{
    commands::ends::{
        config::EndsConfig,
        config_structs::{ClipStrategy, KmerSource, WindowMotifAssigner},
    },
    shared::{
        io::{create_text_writer, dot_join},
        kmers::write::write_category_sparse_with_paths,
    },
};
use anyhow::{Context, Result};
use fxhash::FxHashMap;
use ndarray::Array2;
use ndarray_npy::write_npy;
use std::{io::Write, path::Path};

/// Write the final end-motif count outputs.
///
/// Dense output writes one `.npy` matrix with a shared motif order. Sparse
/// output writes a COO-style `.npz` plus the matching motif label file.
///
/// Parameters
/// ----------
/// - `output_dir`:
///   Directory where the final files should be written
/// - `prefix`:
///   Optional output-file prefix
/// - `bins`:
///   Final decoded per-window motif maps
/// - `motifs`:
///   Final motif column order
/// - `write_dense_output`:
///   Whether to write dense `.npy` output instead of sparse `.npz`
///
/// Returns
/// -------
/// - `Result<()>`:
///   `Ok(())` after the requested output files have been written
pub fn write_end_motif_outputs(
    output_dir: &Path,
    prefix: &str,
    bins: &[FxHashMap<String, f64>],
    motifs: &[String],
    write_dense_output: bool,
) -> Result<()> {
    let motifs_path = output_dir.join(dot_join(&[prefix, "end_motifs.txt"]));
    if write_dense_output {
        let counts_path = output_dir.join(dot_join(&[prefix, "end_motifs.npy"]));
        write_npy(&counts_path, &stack_end_motif_counts(bins, motifs)?)
            .with_context(|| format!("writing {}", counts_path.display()))?;

        let mut motifs_writer = create_text_writer(&motifs_path)
            .with_context(|| format!("create {}", motifs_path.display()))?;
        for motif in motifs {
            writeln!(motifs_writer, "{motif}")
                .with_context(|| format!("write {}", motifs_path.display()))?;
        }
        motifs_writer
            .finish()
            .with_context(|| format!("finalize {}", motifs_path.display()))?;
    } else {
        let counts_path = output_dir.join(dot_join(&[prefix, "end_motifs.sparse.npz"]));
        write_category_sparse_with_paths(bins, motifs, &counts_path, &motifs_path)?;
    }

    Ok(())
}

/// Write the small settings sidecar needed to interpret end-motif outputs.
///
/// This records the motif-definition settings and fragment-length filter basis,
/// but intentionally leaves out output-format details that are already obvious
/// from the files written next to it.
///
/// Parameters
/// ----------
/// - `output_dir`:
///   Directory where the sidecar should be written
/// - `prefix`:
///   Optional output-file prefix
/// - `opt`:
///   Full `ends` configuration used for the run
///
/// Returns
/// -------
/// - `Result<()>`:
///   `Ok(())` after the settings sidecar has been written
pub fn write_end_settings_json(output_dir: &Path, prefix: &str, opt: &EndsConfig) -> Result<()> {
    let settings_path = output_dir.join(dot_join(&[prefix, "end_motif_settings.json"]));
    let mut settings_writer = create_text_writer(&settings_path)
        .with_context(|| format!("create {}", settings_path.display()))?;
    writeln!(settings_writer, "{{")
        .with_context(|| format!("write {}", settings_path.display()))?;
    writeln!(
        settings_writer,
        "  \"source_inside\": \"{}\",",
        kmer_source_name(opt.source_inside)
    )
    .with_context(|| format!("write {}", settings_path.display()))?;
    writeln!(
        settings_writer,
        "  \"clip_strategy\": \"{}\",",
        clip_strategy_name(opt.clip.clip_strategy)
    )
    .with_context(|| format!("write {}", settings_path.display()))?;
    writeln!(
        settings_writer,
        "  \"window_assignment\": \"{}\",",
        window_assigner_name(opt.window_assignment.assign_by)
    )
    .with_context(|| format!("write {}", settings_path.display()))?;
    writeln!(
        settings_writer,
        "  \"collapse_complement\": {}",
        opt.collapse_complement
    )
    .with_context(|| format!("write {}", settings_path.display()))?;
    writeln!(settings_writer, "}}")
        .with_context(|| format!("write {}", settings_path.display()))?;
    settings_writer
        .finish()
        .with_context(|| format!("finalize {}", settings_path.display()))?;
    Ok(())
}

/// Stack sparse per-window motif maps into a dense matrix with a fixed column order.
///
/// Parameters
/// ----------
/// - `bins`:
///   Final decoded per-window motif maps
/// - `motifs`:
///   Fixed motif column order
///
/// Returns
/// -------
/// - `Result<Array2<f64>>`:
///   Dense matrix with one row per window and one column per motif
fn stack_end_motif_counts(
    bins: &[FxHashMap<String, f64>],
    motifs: &[String],
) -> Result<Array2<f64>> {
    let n_rows = bins.len();
    let n_cols = motifs.len();
    let mut mat = Array2::<f64>::zeros((n_rows, n_cols));
    let motif_columns: FxHashMap<&String, usize> = motifs
        .iter()
        .enumerate()
        .map(|(col, motif)| (motif, col))
        .collect();

    for (row, bin) in bins.iter().enumerate() {
        for (motif, &count) in bin {
            let col = motif_columns.get(motif).copied().with_context(|| {
                format!("missing dense output column for motif label '{motif}'")
            })?;
            mat[(row, col)] = count;
        }
    }

    Ok(mat)
}

/// Convert the inside-source enum to its JSON-sidecar string form.
///
/// Parameters
/// ----------
/// - `source`:
///   Inside-sequence source mode
///
/// Returns
/// -------
/// - `&'static str`:
///   Stable sidecar string for that setting
fn kmer_source_name(source: KmerSource) -> &'static str {
    match source {
        KmerSource::Read => "read",
        KmerSource::Reference => "reference",
    }
}

/// Convert the clip-strategy enum to its JSON-sidecar string form.
///
/// Parameters
/// ----------
/// - `strategy`:
///   Clip-handling mode
///
/// Returns
/// -------
/// - `&'static str`:
///   Stable sidecar string for that setting
fn clip_strategy_name(strategy: ClipStrategy) -> &'static str {
    match strategy {
        ClipStrategy::Aligned => "aligned",
        ClipStrategy::Raw => "raw",
        ClipStrategy::Drop => "drop",
    }
}

/// Convert the window-assignment mode to its JSON-sidecar string form.
///
/// Parameters
/// ----------
/// - `assigner`:
///   Window-assignment mode
///
/// Returns
/// -------
/// - `String`:
///   Stable sidecar string for that setting
fn window_assigner_name(assigner: WindowMotifAssigner) -> String {
    match assigner {
        WindowMotifAssigner::Endpoint => "endpoint".to_string(),
        WindowMotifAssigner::CountOverlap => "count-overlap".to_string(),
        WindowMotifAssigner::Any => "any".to_string(),
        WindowMotifAssigner::All => "all".to_string(),
        WindowMotifAssigner::Midpoint => "midpoint".to_string(),
        WindowMotifAssigner::Proportion(value) => {
            format!("proportion={}", format_proportion_threshold(value))
        }
    }
}

/// Format a proportion threshold in a stable human-readable form.
///
/// This avoids scientific notation and trims noisy trailing zeros so the
/// settings sidecar stays easy to read and stable across runs.
///
/// Parameters
/// ----------
/// - `value`:
///   Proportion threshold between 0.0 and 1.0
///
/// Returns
/// -------
/// - `String`:
///   Stable decimal representation for JSON-sidecar output
fn format_proportion_threshold(value: f64) -> String {
    let mut formatted = format!("{value:.15}");
    while formatted.contains('.') && formatted.ends_with('0') {
        formatted.pop();
    }
    if formatted.ends_with('.') {
        formatted.push('0');
    }
    formatted
}

#[cfg(test)]
mod tests {
    include!("write_tests.rs");
}
