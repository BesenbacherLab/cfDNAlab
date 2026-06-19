use crate::{
    commands::ends::{
        config::EndsConfig,
        config_structs::{BaseQualityFilter, ClipStrategy, KmerSource, WindowMotifAssigner},
        counting::EndMotifColumnKind,
    },
    shared::{
        indel_mode::IndelMotifFilterPolicy,
        io::{create_text_writer, dot_join},
    },
};
use anyhow::{Context, Result, ensure};
use std::{
    io::Write,
    path::{Path, PathBuf},
};

/// Write the small settings JSON needed to interpret end-motif outputs.
///
/// This records the motif-definition and target-axis settings needed to
/// interpret the count columns. The Zarr store remains the source of truth for
/// storage mode, row metadata, and count arrays.
///
/// Parameters
/// ----------
/// - `output_dir`:
///   Directory where the settings JSON should be written
/// - `prefix`:
///   Optional output-file prefix
/// - `opt`:
///   Full `ends` configuration used for the run
/// - `motifs_file_column_kind`:
///   Parsed motifs-file target mode when `--motifs-file` was used
///
/// Returns
/// -------
/// - `Result<PathBuf>`:
///   Path to the written settings JSON
pub(crate) fn write_end_settings_json(
    output_dir: &Path,
    prefix: &str,
    opt: &EndsConfig,
    motifs_file_column_kind: Option<EndMotifColumnKind>,
) -> Result<PathBuf> {
    ensure!(
        opt.motifs_file.is_some() == motifs_file_column_kind.is_some(),
        "internal error: motifs-file settings require both the path and parsed motifs-file mode"
    );
    let settings_path = output_dir.join(dot_join(&[prefix, "end_settings.json"]));
    let mut settings_writer = create_text_writer(&settings_path)
        .with_context(|| format!("create {}", settings_path.display()))?;
    let settings_entries: Vec<String> = [
        format!("  \"k_inside\": {}", opt.k_inside),
        format!("  \"k_outside\": {}", opt.k_outside),
        format!("  \"all_motifs\": {}", opt.all_motifs),
        format!(
            "  \"motifs_file\": {}",
            json_path_or_null(opt.motifs_file.as_deref())
        ),
        format!(
            "  \"motifs_file_mode\": {}",
            json_motifs_file_mode_or_null(motifs_file_column_kind)
        ),
        format!(
            "  \"source_inside\": \"{}\"",
            kmer_source_name(opt.source_inside)
        ),
        format!(
            "  \"clip_strategy\": \"{}\"",
            clip_strategy_name(opt.clip.clip_strategy)
        ),
        format!(
            "  \"window_assignment\": \"{}\"",
            window_assigner_name(opt.window_assignment.assign_by)
        ),
        format!(
            "  \"indel_filter\": \"{}\"",
            indel_filter_name(opt.indel_filter)
        ),
        format!(
            "  \"effective_indel_filter\": \"{}\"",
            effective_indel_filter_name(opt.indel_filter, opt.source_inside)
        ),
    ]
    .into_iter()
    .chain(base_quality_filter_settings_entry(&opt.bq_filter))
    .chain(collapse_complement_settings_entry(opt))
    .collect();

    writeln!(settings_writer, "{{")
        .with_context(|| format!("write {}", settings_path.display()))?;
    for (entry_index, entry) in settings_entries.iter().enumerate() {
        let comma = if entry_index + 1 == settings_entries.len() {
            ""
        } else {
            ","
        };
        writeln!(settings_writer, "{entry}{comma}")
            .with_context(|| format!("write {}", settings_path.display()))?;
    }
    writeln!(settings_writer, "}}")
        .with_context(|| format!("write {}", settings_path.display()))?;
    settings_writer
        .finish()
        .with_context(|| format!("finalize {}", settings_path.display()))?;
    Ok(settings_path)
}

fn json_path_or_null(path: Option<&Path>) -> String {
    path.map(|path| json_string(&path.to_string_lossy()))
        .unwrap_or_else(|| "null".to_string())
}

fn json_motifs_file_mode_or_null(column_kind: Option<EndMotifColumnKind>) -> &'static str {
    match column_kind {
        Some(EndMotifColumnKind::Motif) => "\"ungrouped\"",
        Some(EndMotifColumnKind::MotifGroup) => "\"grouped\"",
        None => "null",
    }
}

fn json_string(value: &str) -> String {
    serde_json::to_string(value).expect("JSON string serialization should not fail")
}

#[cfg_attr(not(feature = "ends_experimental"), allow(unused_variables))]
fn collapse_complement_settings_entry(opt: &EndsConfig) -> Option<String> {
    #[cfg(feature = "ends_experimental")]
    {
        return Some(format!(
            "  \"collapse_complement\": {}",
            opt.collapse_complement
        ));
    }

    #[cfg(not(feature = "ends_experimental"))]
    {
        None
    }
}

fn base_quality_filter_settings_entry(filters: &[BaseQualityFilter]) -> Option<String> {
    if filters.is_empty() {
        return None;
    }

    let joined = filters
        .iter()
        .map(|filter| format!("\"{}\"", filter.as_cli_expr()))
        .collect::<Vec<_>>()
        .join(", ");
    Some(format!("  \"bq_filters\": [{joined}]"))
}

/// Convert the inside-source enum to its settings JSON string form.
///
/// Parameters
/// ----------
/// - `source`:
///   Inside-sequence source mode
///
/// Returns
/// -------
/// - `&'static str`:
///   Stable settings string for that setting
fn kmer_source_name(source: KmerSource) -> &'static str {
    match source {
        KmerSource::Read => "read",
        KmerSource::Reference => "reference",
    }
}

/// Convert the clip-strategy enum to its settings JSON string form.
///
/// Parameters
/// ----------
/// - `strategy`:
///   Clip-handling mode
///
/// Returns
/// -------
/// - `&'static str`:
///   Stable settings string for that setting
pub(crate) fn clip_strategy_name(strategy: ClipStrategy) -> &'static str {
    match strategy {
        ClipStrategy::Aligned => "aligned",
        ClipStrategy::IncludeAtAlignedBoundary => "include-at-aligned-boundary",
        ClipStrategy::IncludeAtShiftedBoundary => "include-at-shifted-boundary",
        ClipStrategy::Skip => "skip",
    }
}

/// Convert the window-assignment mode to its settings JSON string form.
///
/// Parameters
/// ----------
/// - `assigner`:
///   Window-assignment mode
///
/// Returns
/// -------
/// - `String`:
///   Stable settings string for that setting
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

/// Convert the indel-filter policy to its settings JSON string form.
///
/// Parameters
/// ----------
/// - `policy`:
///   Indel-handling policy for end motifs
///
/// Returns
/// -------
/// - `&'static str`:
///   Stable settings string for that setting
fn indel_filter_name(policy: IndelMotifFilterPolicy) -> &'static str {
    match policy {
        IndelMotifFilterPolicy::Auto => "auto",
        IndelMotifFilterPolicy::SkipAffectedEnd => "skip-affected-end",
        IndelMotifFilterPolicy::SkipAffectedFragment => "skip-affected-fragment",
    }
}

/// Resolve the indel-filter policy that is actually applied during motif extraction.
///
/// The CLI-level `auto` value depends on where inside-fragment bases come from.
/// Read-backed motifs keep indel-affected ends, while reference-backed motifs
/// skip only the affected end.
///
/// Parameters
/// ----------
/// - `policy`:
///   Configured indel-handling policy for end motifs
/// - `source_inside`:
///   Source for inside-fragment motif bases
///
/// Returns
/// -------
/// - `&'static str`:
///   Effective settings string for that run
fn effective_indel_filter_name(
    policy: IndelMotifFilterPolicy,
    source_inside: KmerSource,
) -> &'static str {
    match policy {
        IndelMotifFilterPolicy::Auto => match source_inside {
            KmerSource::Read => "allow",
            KmerSource::Reference => "skip-affected-end",
        },
        IndelMotifFilterPolicy::SkipAffectedEnd => "skip-affected-end",
        IndelMotifFilterPolicy::SkipAffectedFragment => "skip-affected-fragment",
    }
}

/// Format a proportion threshold in a stable user-readable form.
///
/// This avoids scientific notation and trims noisy trailing zeros so the
/// settings JSON stays easy to read and stable across runs.
///
/// Parameters
/// ----------
/// - `value`:
///   Proportion threshold between 0.0 and 1.0
///
/// Returns
/// -------
/// - `String`:
///   Stable decimal representation for settings JSON output
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
