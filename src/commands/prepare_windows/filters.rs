use crate::commands::prepare_windows::{
    config::PrepareConfig,
    intermediate::{IntermediateWindow, parse_intermediate_line, write_intermediate_window},
    labels::{
        AtomicLabelPart, LabelKey, LabelSchema, LabelTuple, build_tuple_compositions,
        render_label_for_key,
    },
    writers::{ChromTempWriter, ensure_temp_writer_for_chrom, finalize_temp_writers},
};
use crate::shared::io::{TextWriter, create_text_writer, stdout_text_writer};
use crate::shared::tiled_run::make_temp_dir;
use anyhow::{Context, Result, bail};
use fxhash::{FxHashMap, FxHashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

/// Parsed `--exclude-label` rule.
#[derive(Clone, Debug)]
pub struct ExcludeRule {
    pub key: LabelKey,
    pub term: String,
}

/// Parsed `--min-per` rule.
#[derive(Clone, Debug)]
pub struct MinPerRule {
    pub key: LabelKey,
    pub min_count: u64,
}

/// Parse `--exclude-label` specs into resolved rules.
///
/// Use this to turn `KEY=TERM` strings into validated label keys that
/// can be applied during filtering. Terms are kept as raw strings so
/// near-side symbols and composition values can be matched exactly.
///
/// Parameters
/// ----------
/// - `specs`:
///     Raw `KEY=TERM` strings from the CLI.
/// - `label_schema`:
///     Resolved schema for validating label keys.
///
/// Returns
/// -------
/// - `rules`:
///     Parsed exclusion rules in the input order.
pub fn parse_exclude_rules(
    specs: &[String],
    label_schema: &LabelSchema,
) -> Result<Vec<ExcludeRule>> {
    if specs.is_empty() {
        return Ok(Vec::new());
    }

    let mut rules = Vec::with_capacity(specs.len());
    for spec in specs {
        let (key, term) = split_key_value(spec, "exclude label")?;
        let key = label_schema.resolve_key(&key)?;
        rules.push(ExcludeRule { key, term });
    }
    Ok(rules)
}

/// Parse `--min-per` specs into resolved rules.
///
/// This validates that each key is known and each count parses as a non-negative
/// integer so later filtering can use the compiled rules.
///
/// Parameters
/// ----------
/// - `specs`:
///     Raw `KEY=COUNT` strings from the CLI.
/// - `label_schema`:
///     Resolved schema for validating label keys.
///
/// Returns
/// -------
/// - `rules`:
///     Parsed minimum-per rules in the input order.
pub fn parse_min_per_rules(
    specs: &[String],
    label_schema: &LabelSchema,
) -> Result<Vec<MinPerRule>> {
    if specs.is_empty() {
        return Ok(Vec::new());
    }

    let mut rules = Vec::with_capacity(specs.len());
    for spec in specs {
        let (key, raw_count) = split_key_value(spec, "min-per")?;
        let count = raw_count
            .parse::<u64>()
            .with_context(|| format!("Invalid min-per count '{}'", raw_count))?;
        let key = label_schema.resolve_key(&key)?;
        rules.push(MinPerRule {
            key,
            min_count: count,
        });
    }
    Ok(rules)
}

fn split_key_value(input: &str, context: &str) -> Result<(String, String)> {
    let (raw_key, raw_value) = input
        .split_once('=')
        .with_context(|| format!("Invalid {} spec (missing '='): '{}'", context, input))?;

    let key = raw_key.trim();
    if key.is_empty() {
        bail!("Invalid {} spec '{}': key cannot be empty", context, input);
    }

    let value = raw_value.trim();
    Ok((key.to_string(), value.to_string()))
}

fn atomic_value_for_tuple(tuple: &LabelTuple, part: AtomicLabelPart) -> &str {
    match part {
        AtomicLabelPart::Input => tuple.input.as_str(),
        AtomicLabelPart::NearSide => tuple.near_side.as_deref().unwrap_or(""),
        AtomicLabelPart::NearName => tuple.near_name.as_deref().unwrap_or(""),
        AtomicLabelPart::Bin => tuple.bin.as_deref().unwrap_or(""),
        AtomicLabelPart::Cluster => tuple.cluster.as_deref().unwrap_or(""),
    }
}

#[derive(Clone, Debug)]
struct MinPerKeyState {
    key: LabelKey,
    min_count: u64,
    counts: FxHashMap<String, u64>,
    rejected_values: FxHashSet<String>,
}

impl MinPerKeyState {
    fn new(key: LabelKey, min_count: u64) -> Self {
        Self {
            key,
            min_count,
            counts: FxHashMap::default(),
            rejected_values: FxHashSet::default(),
        }
    }
}

/// Filter intermediate windows and write the final output.
///
/// This drops excluded windows before any min-per counting, applies min-per
/// rules with eager bucket removal until the result stabilizes, then writes
/// the requested label columns.
///
/// Parameters
/// ----------
/// - `cfg`:
///     Configuration with output settings and separators.
/// - `temp_entries`:
///     Intermediate temp files created by the main streaming pipeline.
/// - `label_schema`:
///     Resolved label compositions.
/// - `out_labels`:
///     Label keys to write in the final output.
/// - `min_per_rules`:
///     Minimum-per rules to enforce.
/// - `exclude_rules`:
///     Label exclusion rules applied before min-per filtering.
/// - `base_output_dir`:
///     Directory used for creating filtering temp files.
///
/// Returns
/// -------
/// `Ok(())` on success or an error if reading or writing fails.
pub fn filter_and_write_output(
    cfg: &PrepareConfig,
    temp_entries: &[(String, PathBuf)],
    label_schema: &LabelSchema,
    out_labels: &[LabelKey],
    min_per_rules: &[MinPerRule],
    exclude_rules: &[ExcludeRule],
    base_output_dir: &Path,
) -> Result<()> {
    let mut entries = sort_temp_entries(temp_entries);
    if entries.is_empty() {
        return write_output_from_entries(cfg, &entries, label_schema, out_labels, exclude_rules);
    }

    // Drop excluded windows before we do any min-per counting to avoid wasted work
    let mut exclude_rules_for_output: &[ExcludeRule] = exclude_rules;
    if !exclude_rules.is_empty() {
        let filtered_entries = filter_excluded_entries(
            &entries,
            cfg.sep,
            label_schema,
            exclude_rules,
            base_output_dir,
        )?;
        remove_temp_dir_for_entries(&entries)?;
        entries = sort_temp_entries(&filtered_entries);
        // Exclusions already applied to temp files so skip rechecking later
        exclude_rules_for_output = &[];
    }

    if entries.is_empty() {
        // No windows remain after exclusion and early filtering
        return write_output_from_entries(
            cfg,
            &entries,
            label_schema,
            out_labels,
            exclude_rules_for_output,
        );
    }

    let normalized = normalize_min_per_rules(min_per_rules);
    if normalized.is_empty() {
        let result = write_output_from_entries(
            cfg,
            &entries,
            label_schema,
            out_labels,
            exclude_rules_for_output,
        );
        remove_temp_dir_for_entries(&entries)?;
        return result;
    }

    if normalized.len() == 1 {
        let mut state = MinPerKeyState::new(normalized[0].key.clone(), normalized[0].min_count);
        compute_initial_counts(
            &entries,
            cfg.sep,
            label_schema,
            std::slice::from_mut(&mut state),
        )?;
        seed_rejected_values(std::slice::from_mut(&mut state));
        let result = filter_single_key_and_write_output(
            cfg,
            &entries,
            label_schema,
            out_labels,
            exclude_rules_for_output,
            &state,
        );
        remove_temp_dir_for_entries(&entries)?;
        return result;
    }

    let mut states: Vec<MinPerKeyState> = normalized
        .iter()
        .map(|rule| MinPerKeyState::new(rule.key.clone(), rule.min_count))
        .collect();
    compute_initial_counts(&entries, cfg.sep, label_schema, &mut states)?;
    seed_rejected_values(&mut states);

    loop {
        let (next_entries, rejections_added) = filter_min_per_pass(
            &entries,
            cfg.sep,
            label_schema,
            &mut states,
            base_output_dir,
        )?;
        remove_temp_dir_for_entries(&entries)?;
        entries = sort_temp_entries(&next_entries);
        if !rejections_added {
            break;
        }
    }

    let result = write_output_from_entries(
        cfg,
        &entries,
        label_schema,
        out_labels,
        exclude_rules_for_output,
    );
    remove_temp_dir_for_entries(&entries)?;
    result
}

fn filter_excluded_entries(
    entries: &[(String, PathBuf)],
    separator: char,
    label_schema: &LabelSchema,
    exclude_rules: &[ExcludeRule],
    base_output_dir: &Path,
) -> Result<Vec<(String, PathBuf)>> {
    let temp_dir = make_temp_dir(base_output_dir, "prepare_windows_exclude")?;
    let mut writers: FxHashMap<String, ChromTempWriter> = FxHashMap::default();

    let needs_compositions_for_exclude = exclude_rules
        .iter()
        .any(|rule| matches!(rule.key, LabelKey::Composition(_)));

    for (_chrom, path) in entries {
        let file = File::open(path)?;
        let reader = BufReader::with_capacity(1 << 20, file);
        for (line_idx, line_res) in reader.lines().enumerate() {
            let line = line_res?;
            if line.is_empty() {
                continue;
            }
            let window =
                parse_intermediate_line_with_context(&line, line_idx + 1, path, separator)?;
            let tuple_compositions = if needs_compositions_for_exclude {
                build_tuple_compositions(&window.label_tuples, label_schema)
            } else {
                Vec::new()
            };

            // Skip windows that match any exclusion rule before min-per counting
            if tuples_match_exclusion_with_compositions(
                &window.label_tuples,
                exclude_rules,
                &tuple_compositions,
            ) {
                continue;
            }

            let writer = ensure_temp_writer_for_chrom(&window.chrom, &temp_dir, &mut writers)?;
            write_intermediate_window(writer.writer(), &window, separator)?;
        }
    }

    let entries = finalize_temp_writers(&mut writers)?;
    if entries.is_empty() && temp_dir.exists() {
        fs::remove_dir_all(&temp_dir)
            .with_context(|| format!("Removing temp dir {}", temp_dir.display()))?;
    }
    Ok(entries)
}

fn normalize_min_per_rules(rules: &[MinPerRule]) -> Vec<MinPerRule> {
    let mut normalized: Vec<MinPerRule> = Vec::new();
    let mut positions: FxHashMap<LabelKey, usize> = FxHashMap::default();

    for rule in rules {
        if rule.min_count == 0 {
            continue;
        }
        if let Some(&idx) = positions.get(&rule.key) {
            // Keep the strictest minimum when the same key appears multiple times
            normalized[idx].min_count = normalized[idx].min_count.max(rule.min_count);
        } else {
            positions.insert(rule.key.clone(), normalized.len());
            normalized.push(rule.clone());
        }
    }
    normalized
}

fn compute_initial_counts(
    entries: &[(String, PathBuf)],
    separator: char,
    label_schema: &LabelSchema,
    states: &mut [MinPerKeyState],
) -> Result<()> {
    if states.is_empty() {
        return Ok(());
    }

    let needs_compositions_for_counts = states
        .iter()
        .any(|state| matches!(state.key, LabelKey::Composition(_)));

    for (_chrom, path) in entries {
        let file = File::open(path)?;
        let reader = BufReader::with_capacity(1 << 20, file);
        for (line_idx, line_res) in reader.lines().enumerate() {
            let line = line_res?;
            if line.is_empty() {
                continue;
            }
            let window =
                parse_intermediate_line_with_context(&line, line_idx + 1, path, separator)?;
            let tuple_compositions = if needs_compositions_for_counts {
                build_tuple_compositions(&window.label_tuples, label_schema)
            } else {
                Vec::new()
            };

            for state in states.iter_mut() {
                // Count each value once per window to avoid double counting tied labels
                let values = collect_unique_values_for_key(
                    &window.label_tuples,
                    &tuple_compositions,
                    &state.key,
                );
                for value in values {
                    *state.counts.entry(value).or_insert(0) += 1;
                }
            }
        }
    }
    Ok(())
}

fn seed_rejected_values(states: &mut [MinPerKeyState]) {
    // Seed rejected values from the initial counts
    for state in states {
        for (value, count) in state.counts.iter() {
            if *count < state.min_count {
                state.rejected_values.insert(value.clone());
            }
        }
    }
}

fn filter_single_key_and_write_output(
    cfg: &PrepareConfig,
    entries: &[(String, PathBuf)],
    label_schema: &LabelSchema,
    out_labels: &[LabelKey],
    exclude_rules: &[ExcludeRule],
    state: &MinPerKeyState,
) -> Result<()> {
    let mut out: TextWriter = if cfg.output.as_os_str() == "-" {
        stdout_text_writer()
    } else {
        create_text_writer(&cfg.output)?
    };

    // Allowed values are those meeting the minimum count
    let allowed_values: FxHashSet<&str> = state
        .counts
        .iter()
        .filter(|(_, count)| **count >= state.min_count)
        .map(|(value, _)| value.as_str())
        .collect();

    let needs_compositions_for_selection = matches!(state.key, LabelKey::Composition(_))
        || needs_compositions_for_output(out_labels, exclude_rules);

    for (_chrom, path) in entries {
        let file = File::open(path)?;
        let reader = BufReader::with_capacity(1 << 20, file);
        for (line_idx, line_res) in reader.lines().enumerate() {
            let line = line_res?;
            if line.is_empty() {
                continue;
            }
            let window = parse_intermediate_line_with_context(&line, line_idx + 1, path, cfg.sep)?;
            let tuple_compositions = if needs_compositions_for_selection {
                build_tuple_compositions(&window.label_tuples, label_schema)
            } else {
                Vec::new()
            };

            let mut kept_tuples: Vec<LabelTuple> = Vec::new();
            for (tuple_idx, tuple) in window.label_tuples.iter().enumerate() {
                let value = key_value_for_tuple(tuple, tuple_idx, &tuple_compositions, &state.key);
                if value.is_empty() || !allowed_values.contains(value) {
                    continue;
                }
                kept_tuples.push(tuple.clone());
            }

            if kept_tuples.is_empty() {
                continue;
            }

            let kept_tuple_compositions =
                if needs_compositions_for_output(out_labels, exclude_rules) {
                    build_tuple_compositions(&kept_tuples, label_schema)
                } else {
                    Vec::new()
                };

            if tuples_match_exclusion_with_compositions(
                &kept_tuples,
                exclude_rules,
                &kept_tuple_compositions,
            ) {
                continue;
            }

            // Emit only the tuples that survived the min-per filter
            write_output_line(
                &mut out,
                &window.chrom,
                window.start,
                window.end,
                &kept_tuples,
                &kept_tuple_compositions,
                label_schema,
                out_labels,
                cfg.sep,
            )?;
        }
    }

    out.finish()?;
    Ok(())
}

fn filter_min_per_pass(
    entries: &[(String, PathBuf)],
    separator: char,
    label_schema: &LabelSchema,
    states: &mut [MinPerKeyState],
    base_output_dir: &Path,
) -> Result<(Vec<(String, PathBuf)>, bool)> {
    let temp_dir = make_temp_dir(base_output_dir, "prepare_windows_filter")?;
    let mut writers: FxHashMap<String, ChromTempWriter> = FxHashMap::default();
    let mut rejections_added = false;

    let needs_compositions_for_min_per = states
        .iter()
        .any(|state| matches!(state.key, LabelKey::Composition(_)));

    for (_chrom, path) in entries {
        let file = File::open(path)?;
        let reader = BufReader::with_capacity(1 << 20, file);
        for (line_idx, line_res) in reader.lines().enumerate() {
            let line = line_res?;
            if line.is_empty() {
                continue;
            }
            let window =
                parse_intermediate_line_with_context(&line, line_idx + 1, path, separator)?;
            let tuple_compositions = if needs_compositions_for_min_per {
                build_tuple_compositions(&window.label_tuples, label_schema)
            } else {
                Vec::new()
            };

            // Track values before filtering so we can decrement counts once per window
            let mut values_before_filter: Vec<Vec<String>> = Vec::with_capacity(states.len());
            for state in states.iter() {
                values_before_filter.push(collect_unique_values_for_key(
                    &window.label_tuples,
                    &tuple_compositions,
                    &state.key,
                ));
            }

            let mut values_after_filter: Vec<Vec<String>> = vec![Vec::new(); states.len()];
            let mut kept_tuples: Vec<LabelTuple> = Vec::new();

            for (tuple_idx, tuple) in window.label_tuples.iter().enumerate() {
                if tuple_passes_min_per(tuple, tuple_idx, states, &tuple_compositions) {
                    kept_tuples.push(tuple.clone());
                    for (key_idx, state) in states.iter().enumerate() {
                        let value =
                            key_value_for_tuple(tuple, tuple_idx, &tuple_compositions, &state.key);
                        if !value.is_empty() {
                            values_after_filter[key_idx].push(value.to_string());
                        }
                    }
                }
            }

            for values in values_after_filter.iter_mut() {
                deduplicate_values(values);
            }

            // Decrement counts for values that disappeared after filtering this window
            for (idx, state) in states.iter_mut().enumerate() {
                for value in values_before_filter[idx].iter() {
                    if !values_after_filter[idx].contains(value) {
                        if decrement_value_count(state, value)? {
                            rejections_added = true;
                        }
                    }
                }
            }

            if kept_tuples.is_empty() {
                continue;
            }

            let writer = ensure_temp_writer_for_chrom(&window.chrom, &temp_dir, &mut writers)?;
            write_intermediate_window(
                writer.writer(),
                &IntermediateWindow {
                    chrom: window.chrom,
                    start: window.start,
                    end: window.end,
                    label_tuples: kept_tuples,
                },
                separator,
            )?;
        }
    }

    let entries = finalize_temp_writers(&mut writers)?;
    if entries.is_empty() && temp_dir.exists() {
        fs::remove_dir_all(&temp_dir)
            .with_context(|| format!("Removing temp dir {}", temp_dir.display()))?;
    }
    Ok((entries, rejections_added))
}

fn remove_temp_dir_for_entries(entries: &[(String, PathBuf)]) -> Result<()> {
    let Some((_chrom, path)) = entries.first() else {
        return Ok(());
    };
    // Use the first entry to locate the shared temp directory
    let dir = path.parent().context("Missing temp directory")?;
    fs::remove_dir_all(dir).with_context(|| format!("Removing temp dir {}", dir.display()))?;
    Ok(())
}

fn sort_temp_entries(entries: &[(String, PathBuf)]) -> Vec<(String, PathBuf)> {
    let mut sorted: Vec<(String, PathBuf)> = entries.iter().cloned().collect();
    sorted.sort_by(|a, b| a.0.cmp(&b.0));
    sorted
}

fn parse_intermediate_line_with_context(
    line: &str,
    line_no: usize,
    path: &Path,
    separator: char,
) -> Result<IntermediateWindow> {
    parse_intermediate_line(line, separator).with_context(|| {
        format!(
            "Failed to parse intermediate line {} in {}",
            line_no,
            path.display()
        )
    })
}

fn tuple_passes_min_per(
    tuple: &LabelTuple,
    tuple_idx: usize,
    states: &[MinPerKeyState],
    tuple_compositions: &[Vec<String>],
) -> bool {
    for state in states {
        let value = key_value_for_tuple(tuple, tuple_idx, tuple_compositions, &state.key);
        // Reject tuples that do not provide a value for any required key
        if value.is_empty() {
            return false;
        }
        // Reject tuples tied to a value already below the threshold
        if state.rejected_values.contains(value) {
            return false;
        }
    }
    true
}

fn key_value_for_tuple<'a>(
    tuple: &'a LabelTuple,
    tuple_idx: usize,
    tuple_compositions: &'a [Vec<String>],
    key: &'a LabelKey,
) -> &'a str {
    // Composition values are precomputed per tuple index
    match key {
        LabelKey::Atomic(part) => atomic_value_for_tuple(tuple, *part),
        LabelKey::Composition(idx) => tuple_compositions[tuple_idx][*idx].as_str(),
    }
}

fn collect_unique_values_for_key(
    tuples: &[LabelTuple],
    tuple_compositions: &[Vec<String>],
    key: &LabelKey,
) -> Vec<String> {
    let mut values: Vec<String> = Vec::new();
    for (idx, tuple) in tuples.iter().enumerate() {
        let value = key_value_for_tuple(tuple, idx, tuple_compositions, key);
        if !value.is_empty() {
            values.push(value.to_string());
        }
    }
    deduplicate_values(&mut values);
    values
}

fn deduplicate_values(values: &mut Vec<String>) {
    values.sort();
    values.dedup();
}

fn decrement_value_count(state: &mut MinPerKeyState, value: &str) -> Result<bool> {
    let entry = state
        .counts
        .get_mut(value)
        .with_context(|| format!("Missing count for value '{}'", value))?;
    if *entry == 0 {
        bail!("Count underflow for value '{}'", value);
    }
    *entry -= 1;
    // Mark new rejections when values fall below the threshold mid-pass
    if *entry < state.min_count && !state.rejected_values.contains(value) {
        state.rejected_values.insert(value.to_string());
        return Ok(true);
    }
    Ok(false)
}

fn needs_compositions_for_output(out_labels: &[LabelKey], exclude_rules: &[ExcludeRule]) -> bool {
    out_labels
        .iter()
        .any(|key| matches!(key, LabelKey::Composition(_)))
        || exclude_rules
            .iter()
            .any(|rule| matches!(rule.key, LabelKey::Composition(_)))
}

fn tuples_match_exclusion_with_compositions(
    tuples: &[LabelTuple],
    exclude_rules: &[ExcludeRule],
    tuple_compositions: &[Vec<String>],
) -> bool {
    if tuples.is_empty() || exclude_rules.is_empty() {
        return false;
    }

    for rule in exclude_rules {
        match &rule.key {
            LabelKey::Atomic(part) => {
                if tuples
                    .iter()
                    .any(|tuple| atomic_value_for_tuple(tuple, *part) == rule.term)
                {
                    return true;
                }
            }
            LabelKey::Composition(idx) => {
                for values in tuple_compositions {
                    if values
                        .get(*idx)
                        .map(|value| value == &rule.term)
                        .unwrap_or(false)
                    {
                        return true;
                    }
                }
            }
        }
    }
    false
}

fn write_output_from_entries(
    cfg: &PrepareConfig,
    entries: &[(String, PathBuf)],
    label_schema: &LabelSchema,
    out_labels: &[LabelKey],
    exclude_rules: &[ExcludeRule],
) -> Result<()> {
    let mut out: TextWriter = if cfg.output.as_os_str() == "-" {
        stdout_text_writer()
    } else {
        create_text_writer(&cfg.output)?
    };

    let needs_compositions_for_labels = needs_compositions_for_output(out_labels, exclude_rules);

    for (_chrom, path) in entries {
        let file = File::open(path)?;
        let reader = BufReader::with_capacity(1 << 20, file);
        for (line_idx, line_res) in reader.lines().enumerate() {
            let line = line_res?;
            if line.is_empty() {
                continue;
            }
            let window = parse_intermediate_line_with_context(&line, line_idx + 1, path, cfg.sep)?;
            if window.label_tuples.is_empty() {
                // If there are no label tuples, exclusions cannot match anything
                write_output_line(
                    &mut out,
                    &window.chrom,
                    window.start,
                    window.end,
                    &window.label_tuples,
                    &[],
                    label_schema,
                    out_labels,
                    cfg.sep,
                )?;
                continue;
            }

            let tuple_compositions = if needs_compositions_for_labels {
                build_tuple_compositions(&window.label_tuples, label_schema)
            } else {
                Vec::new()
            };
            if tuples_match_exclusion_with_compositions(
                &window.label_tuples,
                exclude_rules,
                &tuple_compositions,
            ) {
                continue;
            }

            write_output_line(
                &mut out,
                &window.chrom,
                window.start,
                window.end,
                &window.label_tuples,
                &tuple_compositions,
                label_schema,
                out_labels,
                cfg.sep,
            )?;
        }
    }

    out.finish()?;
    Ok(())
}

fn write_output_line<W: Write>(
    writer: &mut W,
    chrom: &str,
    start: u32,
    end: u32,
    tuples: &[LabelTuple],
    tuple_compositions: &[Vec<String>],
    label_schema: &LabelSchema,
    out_labels: &[LabelKey],
    separator: char,
) -> Result<()> {
    write!(
        writer,
        "{}{sep}{}{sep}{}",
        chrom,
        start,
        end,
        sep = separator
    )?;
    for key in out_labels {
        let label = render_label_for_key(tuples, tuple_compositions, key, label_schema);
        write!(writer, "{sep}{}", label, sep = separator)?;
    }
    writeln!(writer)?;
    Ok(())
}
