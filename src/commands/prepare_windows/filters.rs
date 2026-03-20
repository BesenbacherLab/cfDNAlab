use crate::commands::prepare_windows::{
    config::PrepareConfig,
    intermediate::{IntermediateWindow, parse_intermediate_line, write_intermediate_window},
    labels::{
        AtomicLabelPart, LabelKey, LabelPartRef, LabelSchema, LabelTuple, build_tuple_compositions,
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

/// Parsed `--exclude-labels` rule.
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

/// Parse `--exclude-labels` specs into resolved rules.
///
/// Use this to turn `KEY=TERM` strings into validated label keys that
/// can be applied during filtering. Rules that target unavailable labels
/// are rejected early. Terms are kept as raw strings so win-direction symbols
/// and composition values can be matched exactly.
///
/// Parameters
/// ----------
/// - `specs`:
///     Raw `KEY=TERM` strings from the CLI.
/// - `label_schema`:
///     Resolved schema for validating label keys.
/// - `available_parts`:
///     Atomic parts that are available for this run.
///
/// Returns
/// -------
/// - `rules`:
///     Parsed exclusion rules in the input order.
pub fn parse_exclude_rules(
    specs: &[String],
    label_schema: &LabelSchema,
    available_parts: &FxHashSet<AtomicLabelPart>,
) -> Result<Vec<ExcludeRule>> {
    if specs.is_empty() {
        return Ok(Vec::new());
    }

    let mut rules = Vec::with_capacity(specs.len());
    let mut composition_part_counts: Vec<Option<usize>> =
        vec![None; label_schema.compositions().len()];
    let mut composition_membership_cache: Vec<Option<Vec<AtomicLabelPart>>> =
        vec![None; label_schema.compositions().len()];
    for spec in specs {
        let (key, term) = split_key_value(spec, "exclude label")?;
        let key = label_schema.resolve_key(&key)?;
        ensure_key_membership_available(
            &key,
            label_schema,
            available_parts,
            &mut composition_membership_cache,
            "exclude label",
        )?;
        if let LabelKey::Composition(idx) = &key {
            let composition = &label_schema.compositions()[*idx];
            let expected_parts =
                composition_flattened_part_count(*idx, label_schema, &mut composition_part_counts);
            let provided_parts = term.split('.').count();
            if expected_parts != provided_parts {
                bail!(
                    "Exclude label '{}' expects {} dot-separated parts, got {}",
                    composition.name,
                    expected_parts,
                    provided_parts
                );
            }
        }
        rules.push(ExcludeRule { key, term });
    }
    Ok(rules)
}

fn composition_flattened_part_count(
    idx: usize,
    label_schema: &LabelSchema,
    composition_part_counts: &mut [Option<usize>],
) -> usize {
    if let Some(count) = composition_part_counts[idx] {
        // Cached counts avoid recomputing nested compositions
        return count;
    }
    let composition = &label_schema.compositions()[idx];
    let mut count = 0;
    for part in &composition.parts {
        match part {
            LabelPartRef::Atomic(_) => count += 1,
            LabelPartRef::Composition(inner_idx) => {
                count += composition_flattened_part_count(
                    *inner_idx,
                    label_schema,
                    composition_part_counts,
                );
            }
        }
    }
    // Memoize the flattened part count for reuse
    composition_part_counts[idx] = Some(count);
    count
}

fn ensure_key_membership_available(
    key: &LabelKey,
    label_schema: &LabelSchema,
    available_parts: &FxHashSet<AtomicLabelPart>,
    composition_membership_cache: &mut [Option<Vec<AtomicLabelPart>>],
    context: &str,
) -> Result<()> {
    let membership = key_membership_signature(key, label_schema, composition_membership_cache);
    let mut missing: Vec<&str> = Vec::new();
    for part in membership {
        if !available_parts.contains(&part) {
            missing.push(part.as_str());
        }
    }
    if missing.is_empty() {
        return Ok(());
    }
    missing.sort();
    missing.dedup();
    let key_name = label_key_name(key, label_schema);
    bail!(
        "Invalid {} key '{}': unavailable labels: {}",
        context,
        key_name,
        missing.join(", ")
    );
}

/// Validate that all composition definitions only use available atomic parts.
///
/// Parameters
/// ----------
/// - `label_schema`:
///     Schema providing composition definitions.
/// - `available_parts`:
///     Atomic parts that are available for this run.
///
/// Returns
/// -------
/// `Ok(())` when every composition only uses available parts.
pub fn validate_compositions_available(
    label_schema: &LabelSchema,
    available_parts: &FxHashSet<AtomicLabelPart>,
) -> Result<()> {
    let mut composition_membership_cache: Vec<Option<Vec<AtomicLabelPart>>> =
        vec![None; label_schema.compositions().len()];
    for idx in 0..label_schema.compositions().len() {
        let key = LabelKey::Composition(idx);
        ensure_key_membership_available(
            &key,
            label_schema,
            available_parts,
            &mut composition_membership_cache,
            "compose",
        )?;
    }
    Ok(())
}

/// Validate that resolved keys only use available atomic parts.
///
/// Parameters
/// ----------
/// - `keys`:
///     Resolved keys to validate.
/// - `label_schema`:
///     Schema providing composition definitions.
/// - `available_parts`:
///     Atomic parts that are available for this run.
/// - `context`:
///     Label source used in error messages.
///
/// Returns
/// -------
/// `Ok(())` when every key only uses available parts.
pub fn validate_available_keys(
    keys: &[LabelKey],
    label_schema: &LabelSchema,
    available_parts: &FxHashSet<AtomicLabelPart>,
    context: &str,
) -> Result<()> {
    let mut composition_membership_cache: Vec<Option<Vec<AtomicLabelPart>>> =
        vec![None; label_schema.compositions().len()];
    for key in keys {
        ensure_key_membership_available(
            key,
            label_schema,
            available_parts,
            &mut composition_membership_cache,
            context,
        )?;
    }
    Ok(())
}

/// Parse `--min-per` specs into resolved rules.
///
/// This validates that each key is known, available for this run, and each count
/// parses as a non-negative integer so later filtering can use the compiled rules.
///
/// Parameters
/// ----------
/// - `specs`:
///     Raw `KEY=COUNT` strings from the CLI.
/// - `label_schema`:
///     Resolved schema for validating label keys.
/// - `available_parts`:
///     Atomic parts that are available for this run.
///
/// Returns
/// -------
/// - `rules`:
///     Parsed minimum-per rules in the input order.
pub fn parse_min_per_rules(
    specs: &[String],
    label_schema: &LabelSchema,
    available_parts: &FxHashSet<AtomicLabelPart>,
) -> Result<Vec<MinPerRule>> {
    if specs.is_empty() {
        return Ok(Vec::new());
    }

    let mut rules = Vec::with_capacity(specs.len());
    let mut composition_membership_cache: Vec<Option<Vec<AtomicLabelPart>>> =
        vec![None; label_schema.compositions().len()];
    for spec in specs {
        let (key, raw_count) = split_key_value(spec, "min-per")?;
        let count = raw_count
            .parse::<u64>()
            .with_context(|| format!("Invalid min-per count '{}'", raw_count))?;
        let key = label_schema.resolve_key(&key)?;
        ensure_key_membership_available(
            &key,
            label_schema,
            available_parts,
            &mut composition_membership_cache,
            "min-per",
        )?;
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

#[inline]
fn atomic_value_for_tuple(tuple: &LabelTuple, part: AtomicLabelPart) -> &str {
    match part {
        AtomicLabelPart::Input => tuple.input.as_str(),
        AtomicLabelPart::NearWindowSide => tuple.near_side.as_deref().unwrap_or(""),
        AtomicLabelPart::NearName => tuple.near_name.as_deref().unwrap_or(""),
        AtomicLabelPart::Bin => tuple.bin.as_deref().unwrap_or(""),
        AtomicLabelPart::Cluster => tuple.cluster.as_deref().unwrap_or(""),
    }
}

/// Track counts and rejections for one min-per rule.
///
/// Holds the running count for each observed value, plus the set of values
/// that have already fallen below the minimum.
#[derive(Clone, Debug)]
pub struct MinPerKeyRuleState {
    key: LabelKey,
    min_count: u64,
    counts: FxHashMap<String, u64>,
    rejected_values: FxHashSet<String>,
}

impl MinPerKeyRuleState {
    /// Create a new min-per rule state.
    ///
    /// Parameters
    /// ----------
    /// - `key`:
    ///     Key to track for this rule.
    /// - `min_count`:
    ///     Minimum count required for a value to remain valid.
    ///
    /// Returns
    /// -------
    /// - `state`:
    ///     Initialized state with empty counts and rejections.
    pub fn new(key: LabelKey, min_count: u64) -> Self {
        Self {
            key,
            min_count,
            counts: FxHashMap::default(),
            rejected_values: FxHashSet::default(),
        }
    }

    /// Add a rejected value to the rule state.
    ///
    /// Parameters
    /// ----------
    /// - `value`:
    ///     Value to mark as rejected.
    pub fn add_rejected_value(&mut self, value: &str) {
        // Store owned values so later lookups can use borrowed str
        self.rejected_values.insert(value.to_string());
    }
}

/// Filter intermediate windows and write the final output.
///
/// This drops excluded windows and computes initial min-per
/// counts (when needed). It then applies min-per rules with bucket removal
/// until the result stabilizes, and writes the requested label columns.
///
/// Min-per rule states
/// -------------------
/// - Each `--min-per` rule builds a `MinPerKeyRuleState` that tracks counts for
///   every observed value of the rule key, for example `input=SampleA` or
///   `core=A.B`.
/// - Values that fall below the minimum are rejected and removed from all
///   affected windows before the next pass. The process repeats until no
///   new rejections are found.
///
/// Min-per Filtering Example
/// -------------------------
/// Given a window with two label tuples:
///
/// - Tuple 1: `input=A, near-name=TSS`
///
/// - Tuple 2: `input=B, near-name=TSS`
///
/// If the rule is `--min-per input=100` and value `A` is below the minimum,
/// Tuple 1 is removed. Tuple 2 can still remain, so the window is kept.
/// If all tuples are removed, the window is dropped.
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
/// - `base_temp_dir`:
///     Per-run temp directory used for creating filtering temp subdirectories.
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
    base_temp_dir: &Path,
    chrom_order: &[String],
) -> Result<()> {
    let chrom_positions = build_chrom_positions(chrom_order);
    let mut entries = sort_temp_entries(temp_entries, &chrom_positions)?;
    if entries.is_empty() {
        return write_empty_output(cfg);
    }

    // Remove duplicate rules (by membership) and zero-min rules
    let normalized_rules = normalize_min_per_rules(min_per_rules, label_schema);
    let mut rule_states: Vec<MinPerKeyRuleState> = normalized_rules
        .iter()
        .map(|rule| MinPerKeyRuleState::new(rule.key.clone(), rule.min_count))
        .collect();

    let do_any_filtering = !rule_states.is_empty() || !exclude_rules.is_empty();

    // Filter by exclude-rules when present and count post-filtering values per state (when min-per rules exist)
    if do_any_filtering {
        let filtered_entries = filter_excluded_entries_and_count(
            &entries,
            cfg.sep,
            label_schema,
            exclude_rules,
            base_temp_dir,
            if rule_states.is_empty() {
                None
            } else {
                Some(&mut rule_states)
            },
        )?;
        entries = sort_temp_entries(&filtered_entries, &chrom_positions)?;
    }

    if entries.is_empty() {
        // No windows remain after exclusion and early filtering
        return write_empty_output(cfg);
    }
    if normalized_rules.is_empty() {
        let result = write_output_from_entries(cfg, &entries, label_schema, out_labels);
        remove_temp_dir_for_entries(&entries)?;
        return result;
    }

    // Min-per rules exist so find the initial values to reject in the filtering passes
    initialize_rejected_values(&mut rule_states);

    // Single-rule is a simplified path (no iteration)
    if rule_states.len() == 1 {
        let result = filter_single_key_and_write_output(
            cfg,
            &entries,
            label_schema,
            out_labels,
            &rule_states[0],
        );
        remove_temp_dir_for_entries(&entries)?;
        return result;
    }

    // Multi-rule requires iterative filtering until stable
    loop {
        let (next_entries, rejections_added) = filter_min_per_pass(
            &entries,
            cfg.sep,
            label_schema,
            &mut rule_states,
            base_temp_dir,
        )?;
        remove_temp_dir_for_entries(&entries)?;
        entries = sort_temp_entries(&next_entries, &chrom_positions)?;
        if !rejections_added {
            break;
        }
    }

    let result = write_output_from_entries(cfg, &entries, label_schema, out_labels);
    remove_temp_dir_for_entries(&entries)?;
    result
}

/// Filter out excluded windows and optionally compute min-per counts.
///
/// This writes a new set of temp files that skip excluded windows. When
/// `rule_states` is provided, it also tallies per-window counts for each
/// min-per rule using only the kept windows.
///
/// Parameters
/// ----------
/// - `entries`:
///   Sorted temp files containing intermediate windows.
/// - `separator`:
///   Field separator used in the temp files.
/// - `label_schema`:
///   Schema for resolving composition values.
/// - `exclude_rules`:
///   Label exclusion rules applied before min-per filtering.
/// - `base_temp_dir`:
///   Parent directory for creating a new temp directory.
/// - `rule_states`:
///   Optional min-per rule states to update with counts.
/// Returns
/// -------
/// - `entries`:
///     Temp files containing only the kept windows.
fn filter_excluded_entries_and_count(
    entries: &[(String, PathBuf)],
    separator: char,
    label_schema: &LabelSchema,
    exclude_rules: &[ExcludeRule],
    base_temp_dir: &Path,
    rule_states: Option<&mut [MinPerKeyRuleState]>,
) -> Result<Vec<(String, PathBuf)>> {
    let has_exclusions = !exclude_rules.is_empty();
    let temp_dir = if has_exclusions {
        Some(make_temp_dir(base_temp_dir, "prepare_windows_exclude")?)
    } else {
        None
    };
    let mut writers: FxHashMap<String, ChromTempWriter> = FxHashMap::default();

    let mut rule_states = rule_states;
    let needs_compositions_for_exclude = has_exclusions
        && exclude_rules
            .iter()
            .any(|rule| matches!(rule.key, LabelKey::Composition(_)));
    let needs_compositions_for_counts = rule_states
        .as_ref()
        .map(|states| {
            states
                .iter()
                .any(|state| matches!(state.key, LabelKey::Composition(_)))
        })
        .unwrap_or(false);
    let needs_compositions = needs_compositions_for_exclude || needs_compositions_for_counts;

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
            let tuple_compositions = if needs_compositions {
                build_tuple_compositions(&window.label_tuples, label_schema)
            } else {
                Vec::new()
            };

            // Skip windows that match any exclusion rule before min-per counting
            if has_exclusions
                && tuples_match_exclusion_with_compositions(
                    &window.label_tuples,
                    exclude_rules,
                    &tuple_compositions,
                )
            {
                continue;
            }

            // Count initial values (post-filtering) for downstream min-per filtering
            if let Some(rule_states) = rule_states.as_deref_mut() {
                // Count each value once per window to avoid double counting tied labels
                for state in rule_states.iter_mut() {
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

            if has_exclusions {
                let temp_dir = temp_dir
                    .as_ref()
                    .expect("Temp dir should exist when exclusions are enabled");
                let writer = ensure_temp_writer_for_chrom(&window.chrom, temp_dir, &mut writers)?;
                write_intermediate_window(writer.writer(), &window, separator)?;
            }
        }
    }

    let entries = if has_exclusions {
        let entries = finalize_temp_writers(&mut writers)?;
        if entries.is_empty() {
            let temp_dir = temp_dir
                .as_ref()
                .expect("Temp dir should exist when exclusions are enabled");
            if temp_dir.exists() {
                fs::remove_dir_all(temp_dir)
                    .with_context(|| format!("Removing temp dir {}", temp_dir.display()))?;
            }
        }
        entries
    } else {
        entries.to_vec()
    };

    Ok(entries)
}

/// Normalize `--min-per` rules so downstream passes are deterministic.
///
/// Drops any rule with a zero minimum and merges duplicate keys by keeping
/// the strictest minimum. Duplicate detection compares the full membership
/// set for each key, ignoring order, so `a,b` matches `b,a`. The first
/// occurrence of each membership set defines the output order. Removed
/// duplicates are reported to the user.
///
/// Parameters
/// ----------
/// - `rules`:
///   Raw rules in the CLI order.
/// - `label_schema`:
///   Schema used to resolve composition membership sets.
///
/// Returns
/// -------
/// - `normalized`:
///   De-duplicated rules with zero-minimum entries removed.
///   Membership-based duplicates keep the strictest minimum.
pub fn normalize_min_per_rules(
    rules: &[MinPerRule],
    label_schema: &LabelSchema,
) -> Vec<MinPerRule> {
    let mut normalized: Vec<MinPerRule> = Vec::new();
    let mut positions: FxHashMap<LabelKey, usize> = FxHashMap::default();
    let mut membership_positions: FxHashMap<Vec<AtomicLabelPart>, usize> = FxHashMap::default();
    let mut composition_membership_cache: Vec<Option<Vec<AtomicLabelPart>>> =
        vec![None; label_schema.compositions().len()];
    let mut removed_duplicates: Vec<(String, String)> = Vec::new();

    for rule in rules {
        if rule.min_count == 0 {
            continue;
        }
        // Membership is the unique set of atomic parts the rule depends on
        let membership_parts =
            key_membership_signature(&rule.key, label_schema, &mut composition_membership_cache);
        if let Some(&idx) = membership_positions.get(&membership_parts) {
            // Keep the strictest minimum when identical memberships appear multiple times
            if normalized[idx].min_count < rule.min_count {
                normalized[idx].min_count = rule.min_count;
            }
            removed_duplicates.push((
                label_key_name(&rule.key, label_schema),
                label_key_name(&normalized[idx].key, label_schema),
            ));
            continue;
        }
        if let Some(&idx) = positions.get(&rule.key) {
            // Keep the strictest minimum when the same key appears multiple times
            normalized[idx].min_count = normalized[idx].min_count.max(rule.min_count);
            continue;
        }
        positions.insert(rule.key.clone(), normalized.len());
        membership_positions.insert(membership_parts, normalized.len());
        normalized.push(rule.clone());
    }

    if !removed_duplicates.is_empty() {
        eprintln!("Removed duplicate min-per rules with identical membership");
        for (removed, kept) in removed_duplicates {
            eprintln!("Removed '{}' because it matches '{}'", removed, kept);
        }
    }
    normalized
}

/// Resolve a key to its membership set of atomic parts.
///
/// This flattens compositions into their atomic parts so identical memberships
/// compare equal even when names differ.
///
/// Parameters
/// ----------
/// - `key`:
///   Key to resolve.
/// - `label_schema`:
///   Schema providing composition definitions.
/// - `composition_signatures`:
///   Cache for resolved composition memberships.
///
/// Returns
/// -------
/// - `membership`:
///   Sorted unique atomic parts that define the key membership.
fn key_membership_signature(
    key: &LabelKey,
    label_schema: &LabelSchema,
    composition_signatures: &mut [Option<Vec<AtomicLabelPart>>],
) -> Vec<AtomicLabelPart> {
    match key {
        LabelKey::Atomic(part) => vec![*part],
        LabelKey::Composition(idx) => {
            composition_membership_signature(*idx, label_schema, composition_signatures)
        }
    }
}

/// Resolve a composition to its atomic membership set.
///
/// This flattens nested compositions and deduplicates the resulting
/// atomic parts so membership comparisons ignore order.
///
/// Parameters
/// ----------
/// - `idx`:
///     Composition index in schema order.
/// - `label_schema`:
///     Schema providing composition definitions.
/// - `composition_signatures`:
///     Cache for resolved composition memberships.
///
/// Returns
/// -------
/// - `membership`:
///     Sorted unique atomic parts for the composition.
fn composition_membership_signature(
    idx: usize,
    label_schema: &LabelSchema,
    composition_signatures: &mut [Option<Vec<AtomicLabelPart>>],
) -> Vec<AtomicLabelPart> {
    if let Some(signature) = composition_signatures[idx].as_ref() {
        // Cached membership avoids recomputing nested compositions
        return signature.clone();
    }
    let composition = &label_schema.compositions()[idx];
    let mut membership: Vec<AtomicLabelPart> = Vec::new();
    for part in &composition.parts {
        match part {
            LabelPartRef::Atomic(atomic) => membership.push(*atomic),
            LabelPartRef::Composition(inner_idx) => membership.extend(
                composition_membership_signature(*inner_idx, label_schema, composition_signatures),
            ),
        }
    }
    membership.sort();
    membership.dedup();
    // Memoize the flattened membership so recursive lookups stay fast
    composition_signatures[idx] = Some(membership.clone());
    membership
}

/// Render a display name for a key.
///
/// Parameters
/// ----------
/// - `key`:
///     Key to label.
/// - `label_schema`:
///     Schema providing composition names.
///
/// Returns
/// -------
/// - `name`:
///     Human-readable key name.
fn label_key_name(key: &LabelKey, label_schema: &LabelSchema) -> String {
    match key {
        LabelKey::Atomic(part) => part.as_str().to_string(),
        LabelKey::Composition(idx) => label_schema.compositions()[*idx].name.clone(),
    }
}

/// Initialize rejections for values already below their minimums.
///
/// This avoids extra work in the first filtering pass by skipping values
/// that start out below the threshold.
///
/// Parameters
/// ----------
/// - `rule_states`:
///     Per-key rule state to update in place.
fn initialize_rejected_values(rule_states: &mut [MinPerKeyRuleState]) {
    // Initialize rejected values from the initial counts
    for state in rule_states {
        for (value, count) in state.counts.iter() {
            if *count < state.min_count {
                state.rejected_values.insert(value.clone());
            }
        }
    }
}

/// Filter by a single min-per key and stream results to the final output.
///
/// This is a fast path that avoids multi-key state churn when only one
/// min-per rule is active.
///
/// Parameters
/// ----------
/// - `cfg`:
///     Output configuration and delimiter settings.
/// - `entries`:
///     Sorted temp files containing intermediate windows.
/// - `label_schema`:
///     Schema for resolving composition values.
/// - `out_labels`:
///     Label keys to write in the final output.
/// - `state`:
///     Precomputed min-per state for the single key.
///
/// Returns
/// -------
/// `Ok(())` on success or an error if reading or writing fails.
fn filter_single_key_and_write_output(
    cfg: &PrepareConfig,
    entries: &[(String, PathBuf)],
    label_schema: &LabelSchema,
    out_labels: &[LabelKey],
    state: &MinPerKeyRuleState,
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
        || needs_compositions_for_output(out_labels, &[]);

    for (_chrom, path) in entries {
        let file = File::open(path)?;
        let reader = BufReader::with_capacity(1 << 20, file);
        for (line_idx, line_res) in reader.lines().enumerate() {
            let line = line_res?;
            if line.is_empty() {
                continue;
            }
            let window = parse_intermediate_line_with_context(&line, line_idx + 1, path, cfg.sep)?;
            // Build composition values once when min-per or output depends on them
            let tuple_compositions = if needs_compositions_for_selection {
                build_tuple_compositions(&window.label_tuples, label_schema)
            } else {
                Vec::new()
            };

            // Apply min-per by keeping only tuples whose key value is allowed
            let mut kept_tuples: Vec<LabelTuple> = Vec::new();
            for (tuple_idx, tuple) in window.label_tuples.iter().enumerate() {
                let value = key_value_for_tuple(tuple, tuple_idx, &tuple_compositions, &state.key);
                if value.is_empty() || !allowed_values.contains(value) {
                    continue;
                }
                kept_tuples.push(tuple.clone());
            }

            // Drop windows that lose all tuples after min-per filtering
            if kept_tuples.is_empty() {
                continue;
            }

            // Rebuild compositions for the kept tuples before exclusions and output
            let kept_tuple_compositions = if needs_compositions_for_output(out_labels, &[]) {
                build_tuple_compositions(&kept_tuples, label_schema)
            } else {
                Vec::new()
            };

            // Write only the tuples that survived the min-per filter
            write_output_line(
                &mut out,
                &window.chrom,
                window.start(),
                window.end(),
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

/// Run one filtering pass for multi-key min-per rules.
///
/// Writes a new set of temp files containing only tuples that pass all
/// current rejection sets. This pass may add new rejections as counts
/// fall below thresholds.
///
/// Parameters
/// ----------
/// - `entries`:
///     Sorted temp files containing intermediate windows.
/// - `separator`:
///     Field separator used in the temp files.
/// - `label_schema`:
///     Schema for resolving composition values.
/// - `rule_states`:
///     Per-key rule state to update.
/// - `base_temp_dir`:
///     Parent directory for creating a new pass directory.
///
/// Returns
/// -------
/// - `entries`:
///     Temp files for the next pass.
/// - `rejections_added`:
///     True when new rejections were discovered.
fn filter_min_per_pass(
    entries: &[(String, PathBuf)],
    separator: char,
    label_schema: &LabelSchema,
    rule_states: &mut [MinPerKeyRuleState],
    base_temp_dir: &Path,
) -> Result<(Vec<(String, PathBuf)>, bool)> {
    let temp_dir = make_temp_dir(base_temp_dir, "prepare_windows_filter")?;
    let mut writers: FxHashMap<String, ChromTempWriter> = FxHashMap::default();
    let mut rejections_added = false;

    let needs_compositions_for_min_per = rule_states
        .iter()
        .any(|state| matches!(state.key, LabelKey::Composition(_)));

    // One pass reads all windows, drops tuples that violate any min-per rule,
    // and decrements counts for values that disappear from a window.
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
            // Build composition values once when any rule targets a composition
            let tuple_compositions = if needs_compositions_for_min_per {
                build_tuple_compositions(&window.label_tuples, label_schema)
            } else {
                Vec::new()
            };

            let filter_data = collect_min_per_window_filter_data(
                &window.label_tuples,
                &tuple_compositions,
                rule_states,
            );
            let values_before_filter = filter_data.values_before_filter;
            let values_after_filter = filter_data.values_after_filter;
            let kept_tuples = filter_data.kept_tuples;

            // Decrement counts for values that disappeared after filtering this window
            for (idx, state) in rule_states.iter_mut().enumerate() {
                // Compare the before and after sets using the same per-state index
                for value in values_before_filter[idx].iter() {
                    // Decrement only when a value was present before but not after filtering
                    if !values_after_filter[idx].contains(value) {
                        if decrement_value_count(state, value)? {
                            rejections_added = true;
                        }
                    }
                }
            }

            // Drop windows that lose all tuples after min-per filtering
            if kept_tuples.is_empty() {
                continue;
            }

            let writer = ensure_temp_writer_for_chrom(&window.chrom, &temp_dir, &mut writers)?;
            let window_interval =
                crate::shared::interval::Interval::new(window.start(), window.end())?;
            let filtered_window = IntermediateWindow::new(
                window.chrom,
                window_interval,
                kept_tuples,
            );
            write_intermediate_window(writer.writer(), &filtered_window, separator)?;
        }
    }

    let entries = finalize_temp_writers(&mut writers)?;
    if entries.is_empty() && temp_dir.exists() {
        fs::remove_dir_all(&temp_dir)
            .with_context(|| format!("Removing temp dir {}", temp_dir.display()))?;
    }
    Ok((entries, rejections_added))
}

/// Values collected for one min-per filtering window.
///
/// Captures which tuples remain after applying rejection sets, plus the
/// per-rule values before and after filtering for count updates.
pub struct MinPerWindowFilterData {
    pub values_before_filter: Vec<Vec<String>>,
    pub values_after_filter: Vec<Vec<String>>,
    pub kept_tuples: Vec<LabelTuple>,
}

/// Collect per-rule values and kept tuples for a single window.
///
/// Builds the before and after value sets, applying each min-per rule in order.
/// The resulting vectors keep the same order as `rule_states` so index positions
/// match when counts are decremented.
///
/// Parameters
/// ----------
/// - `label_tuples`:
///     Tuples to evaluate for the window.
/// - `tuple_compositions`:
///     Precomputed composition values for each tuple.
/// - `rule_states`:
///     Min-per rule states providing keys and rejected values.
///
/// Returns
/// -------
/// - `data`:
///     Values before and after filtering, plus the kept tuples.
pub fn collect_min_per_window_filter_data(
    label_tuples: &[LabelTuple],
    tuple_compositions: &[Vec<String>],
    rule_states: &[MinPerKeyRuleState],
) -> MinPerWindowFilterData {
    // Track values before filtering so we can decrement counts once per window
    let mut values_before_filter: Vec<Vec<String>> = Vec::with_capacity(rule_states.len());
    for state in rule_states.iter() {
        values_before_filter.push(collect_unique_values_for_key(
            label_tuples,
            tuple_compositions,
            &state.key,
        ));
    }

    // Keep tuples that pass every min-per rule and record their values
    let mut values_after_filter: Vec<Vec<String>> = vec![Vec::new(); rule_states.len()];
    let mut kept_tuples: Vec<LabelTuple> = Vec::new();

    // Index positions match rule_states, so values_before_filter[idx] and values_after_filter[idx]
    // always refer to the same min-per rule when we compare and decrement
    for (tuple_idx, tuple) in label_tuples.iter().enumerate() {
        // Collect values in rule_state order so we can reuse them for values_after_filter
        let mut tuple_values: Vec<&str> = Vec::with_capacity(rule_states.len());
        let mut passes = true;
        for state in rule_states.iter() {
            let value = key_value_for_tuple(tuple, tuple_idx, tuple_compositions, &state.key);
            // Reject tuples with missing values or previously rejected values
            if value.is_empty() || state.rejected_values.contains(value) {
                passes = false;
                break;
            }
            // Store the value so we can reuse it when building values_after_filter
            tuple_values.push(value);
        }
        if !passes {
            // Tuple values may be shorter when we break early, but we skip zip in that case
            continue;
        }
        kept_tuples.push(tuple.clone());
        // Reuse the collected values so we do not resolve the same keys twice
        // Preserve the same state ordering for the values_before_filter comparison
        for (value, values_for_key) in tuple_values.iter().zip(values_after_filter.iter_mut()) {
            // Push into the bucket for the matching state by index
            if !value.is_empty() {
                values_for_key.push((*value).to_string());
            }
        }
    }

    for values in values_after_filter.iter_mut() {
        deduplicate_values(values);
    }

    MinPerWindowFilterData {
        values_before_filter,
        values_after_filter,
        kept_tuples,
    }
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

fn build_chrom_positions(chrom_order: &[String]) -> FxHashMap<String, usize> {
    let mut positions = FxHashMap::default();
    // Preserve input chromosome order for deterministic temp replay
    for (idx, chrom) in chrom_order.iter().enumerate() {
        positions.insert(chrom.clone(), idx);
    }
    positions
}

fn sort_temp_entries(
    entries: &[(String, PathBuf)],
    chrom_positions: &FxHashMap<String, usize>,
) -> Result<Vec<(String, PathBuf)>> {
    let mut tagged_entries: Vec<(usize, (String, PathBuf))> = Vec::with_capacity(entries.len());
    for entry in entries {
        // entry.0 is the chromosome name paired with the temp path
        let position = chrom_positions
            .get(&entry.0)
            .copied()
            .with_context(|| format!("Missing chromosome '{}' in input ordering", entry.0))?;
        tagged_entries.push((position, entry.clone()));
    }
    // Sort by input chromosome order, each per-chrom temp file writes windows in streaming order
    // and streaming enforces non-decreasing starts within a chromosome
    tagged_entries.sort_unstable_by_key(|(position, _)| *position);
    Ok(tagged_entries
        .into_iter()
        .map(|(_position, entry)| entry)
        .collect())
}

#[inline]
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

/// Resolve the value for a key from a tuple.
///
/// Parameters
/// ----------
/// - `tuple`:
///     Label tuple to inspect.
/// - `tuple_idx`:
///     Index of the tuple in the current window.
/// - `tuple_compositions`:
///     Precomputed composition values for each tuple.
/// - `key`:
///     The key to resolve.
///
/// Returns
/// -------
/// - `value`:
///     Resolved value or an empty string when the key is missing.
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

/// Collect unique values for a key across tuples.
///
/// Parameters
/// ----------
/// - `tuples`:
///     Label tuples to inspect.
/// - `tuple_compositions`:
///     Precomputed composition values for each tuple.
/// - `key`:
///     Key to resolve.
///
/// Returns
/// -------
/// - `values`:
///     Sorted unique values with empty strings removed.
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

#[inline]
/// Sort and deduplicate a vector of values in place.
///
/// Parameters
/// ----------
/// - `values`:
///     Values to sort and deduplicate.
fn deduplicate_values(values: &mut Vec<String>) {
    values.sort();
    values.dedup();
}

/// Decrement a value count and update rejection state.
///
/// This detects when a count crosses below the minimum so the caller can
/// trigger another filtering pass.
///
/// Parameters
/// ----------
/// - `state`:
///     Per-key state holding counts and rejections.
/// - `value`:
///     Value to decrement.
///
/// Returns
/// -------
/// - `added_rejection`:
///     True when the value newly falls below the threshold.
fn decrement_value_count(state: &mut MinPerKeyRuleState, value: &str) -> Result<bool> {
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

#[inline]
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
) -> Result<()> {
    let mut out: TextWriter = if cfg.output.as_os_str() == "-" {
        stdout_text_writer()
    } else {
        create_text_writer(&cfg.output)?
    };

    let needs_compositions_for_labels = needs_compositions_for_output(out_labels, &[]);

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
                    window.start(),
                    window.end(),
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

            write_output_line(
                &mut out,
                &window.chrom,
                window.start(),
                window.end(),
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

fn write_empty_output(cfg: &PrepareConfig) -> Result<()> {
    // Ensure the output file exists even when no windows are produced
    let out: TextWriter = if cfg.output.as_os_str() == "-" {
        stdout_text_writer()
    } else {
        create_text_writer(&cfg.output)?
    };
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
