//! Parser and post-processing helpers for `ends --motifs-file`.
//!
//! A motifs file is parsed once before tile processing starts. The parser validates every public
//! label, builds the output target axis in file order, and precomputes encoded lookup keys for both
//! fragment ends. The counting code can then do one hash lookup per observed end motif and skip all
//! unselected motifs without allocating motif strings.
//!
//! The important orientation rule is that right-end observations are represented by their
//! reverse-complemented encoded state with `reverse_on_decode = true`. The lookup stores that state
//! explicitly. It does not canonicalize or collapse complements, because the motifs file may map a
//! motif and its reverse complement to different targets.

use crate::{
    commands::ends::{
        counting::{
            EncodedEndMotifKey, EndMotifColumnKind, EndMotifCounts, EndMotifHalfSpec,
            SelectedEndCountsByWindow, SelectedEndMotifLookup,
        },
        motifs::build_optional_kmer_spec,
        output::ensure_dense_end_motif_output_size,
    },
    shared::{
        base::rev_complement,
        kmers::kmer_codec::{MAX_RADIX5_KMER_SIZE, build_subspace_kmer_spec},
    },
};
use anyhow::{Context, Result, bail, ensure};
use fxhash::{FxHashMap, FxHashSet};
use std::{fs, path::Path, sync::Arc};

/// Shape of the motifs file after splitting each row on tabs.
///
/// The file is intentionally one mode at a time. A mixed file would make the count axis hard to
/// explain because some rows would target motifs and other rows would target groups.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MotifsFileMode {
    /// One-column file where every motif is its own output target
    Ungrouped,
    /// Two-column file where the second column names the output target
    Grouped,
}

/// Parsed and normalized representation of one motifs-file motif label.
///
/// The public label is always normalized to `<outside>_<inside>`, even when the input used the
/// short form allowed for one-sided motifs. `outside` and `inside` are stored separately because
/// the encoded lookup must preserve the two-half end motif model used by the counter.
struct ParsedEndMotifLabel {
    /// Normalized public motif label
    label: String,
    /// Uppercase outside-of-fragment sequence, empty when `k_outside = 0`
    outside: String,
    /// Uppercase inside-of-fragment sequence, empty when `k_inside = 0`
    inside: String,
}

struct ParsedMotifsFileRow {
    motif: ParsedEndMotifLabel,
    target_idx: u32,
    line_number: usize,
}

/// Parse a selected end-motifs file from disk.
///
/// This validates user-facing motif labels, builds the selected output target
/// axis, and precomputes the encoded lookup used by counting. Filtering must use
/// this encoded lookup before sparse count insertion.
///
/// Parameters
/// ----------
/// - `path`:
///   Motifs-file path supplied by the user
/// - `k_inside`:
///   Expected inside motif length from `--k-inside`
/// - `k_outside`:
///   Expected outside motif length from `--k-outside`
///
/// Returns
/// -------
/// - `Result<SelectedEndMotifLookup>`:
///   Parsed output axis and encoded lookup ready for tile-local counting
pub(crate) fn parse_selected_end_motifs_file(
    path: &Path,
    k_inside: usize,
    k_outside: usize,
) -> Result<SelectedEndMotifLookup> {
    let contents = fs::read_to_string(path)
        .with_context(|| format!("reading motifs file {}", path.display()))?;
    parse_selected_end_motifs(&contents, k_inside, k_outside)
        .with_context(|| format!("parsing motifs file {}", path.display()))
}

/// Convert reduced selected motif counts into compact target-indexed output bins.
///
/// Counts are stored by the target index assigned while parsing the motifs file.
/// This keeps observed targets in file order and remaps them to the compact
/// output axis used by the writer.
///
/// When `include_all_targets` is false, targets with no observed counts are dropped from the final
/// axis. When it is true, every motifs-file target is retained and the dense-output size check is
/// applied before the writer allocates the full dense matrix.
///
/// Parameters
/// ----------
/// - `counts_by_window`:
///   Reduced sparse counts keyed by global output row and original motifs-file target index
/// - `total_windows`:
///   Number of rows in the final output matrix
/// - `lookup`:
///   Parsed motifs-file lookup that owns the original target labels
/// - `include_all_targets`:
///   Whether to keep unobserved motifs-file targets in the final output axis
///
/// Returns
/// -------
/// - `Result<(Vec<FxHashMap<u32, f64>>, Vec<String>)>`:
///   Compact output bins keyed by final column index, plus final column labels in order
pub(crate) fn postprocess_selected_end_motif_counts(
    counts_by_window: SelectedEndCountsByWindow,
    total_windows: usize,
    lookup: &SelectedEndMotifLookup,
    include_all_targets: bool,
) -> Result<(Vec<FxHashMap<u32, f64>>, Vec<String>)> {
    let mut raw_bins = vec![FxHashMap::default(); total_windows];
    let mut observed_targets = vec![false; lookup.labels.len()];

    for (original_idx, counts) in counts_by_window {
        // Validate the global row index before touching the output vector
        let idx: usize = original_idx
            .try_into()
            .context("selected motif window index does not fit in usize")?;
        ensure!(
            idx < raw_bins.len(),
            "selected motif window index {} is out of bounds for {} output windows",
            idx,
            raw_bins.len()
        );

        for (target_idx, value) in counts {
            // Tile merge already applies this check, but keeping it here protects the writer
            // boundary if another caller is added later
            if !EndMotifCounts::should_store_weight(value)? {
                continue;
            }
            let target_position = target_idx as usize;
            ensure!(
                target_position < lookup.labels.len(),
                "selected motif target index {} is out of bounds for {} targets",
                target_position,
                lookup.labels.len()
            );
            observed_targets[target_position] = true;
            // Multiple tiles may contribute to the same final row and target
            *raw_bins[idx].entry(target_idx).or_insert(0.0) += value;
        }
    }

    // Build the public output axis in motifs-file order
    let mut output_index_by_target = vec![None; lookup.labels.len()];
    let mut motif_order = Vec::new();
    for (target_position, label) in lookup.labels.iter().enumerate() {
        if include_all_targets || observed_targets[target_position] {
            let output_index = u32::try_from(motif_order.len())
                .context("selected motif output index does not fit in u32")?;
            output_index_by_target[target_position] = Some(output_index);
            motif_order.push(label.clone());
        }
    }

    if include_all_targets {
        ensure_dense_end_motif_output_size(total_windows, motif_order.len())?;
    }

    // Remap original motifs-file target indices to compact output column indices
    let mut all_bins = vec![FxHashMap::default(); total_windows];
    for (row_idx, raw_bin) in raw_bins.into_iter().enumerate() {
        for (target_idx, value) in raw_bin {
            let target_position = target_idx as usize;
            let output_index = output_index_by_target
                .get(target_position)
                .copied()
                .flatten()
                .with_context(|| {
                    format!(
                        "selected motif target index {} is missing from the output axis",
                        target_position
                    )
                })?;
            all_bins[row_idx].insert(output_index, value);
        }
    }

    Ok((all_bins, motif_order))
}

/// Parse motifs-file text into a selected motif lookup.
///
/// This is the testable core behind [`parse_selected_end_motifs_file`]. It validates a consistent
/// one-column or two-column file, normalizes motif labels, assigns output target indices, and
/// inserts both left-end and right-end encoded states into the lookup.
///
/// Parameters
/// ----------
/// - `contents`:
///   Raw motifs-file text
/// - `k_inside`:
///   Expected inside motif length
/// - `k_outside`:
///   Expected outside motif length
///
/// Returns
/// -------
/// - `Result<SelectedEndMotifLookup>`:
///   Output labels and encoded target lookup
fn parse_selected_end_motifs(
    contents: &str,
    k_inside: usize,
    k_outside: usize,
) -> Result<SelectedEndMotifLookup> {
    ensure!(
        !contents.is_empty(),
        "--motifs-file must contain at least one motif"
    );

    let mut mode = None;
    let mut labels = Vec::new();
    let mut target_by_group = FxHashMap::default();
    let mut seen_motif_labels = FxHashSet::default();
    let mut rows = Vec::new();

    for (line_index, line) in contents.lines().enumerate() {
        let line_number = line_index + 1;
        // A motifs file is a simple TSV
        let line = line.trim_end_matches('\r');
        let columns: Vec<&str> = line.split('\t').collect();
        let row_mode = parse_row_mode(&columns, line_number)?;
        match mode {
            Some(expected_mode) => ensure!(
                expected_mode == row_mode,
                "line {line_number}: --motifs-file must use either one column for every row or two columns for every row"
            ),
            None => mode = Some(row_mode),
        }

        // Normalize once, before duplicate checks and lookup-key generation
        let motif = parse_motif_label(columns[0], k_inside, k_outside, line_number)?;
        ensure!(
            seen_motif_labels.insert(motif.label.clone()),
            "line {line_number}: duplicate motif `{}` in --motifs-file",
            motif.label
        );

        let target_idx = match row_mode {
            // In ungrouped mode each motif creates exactly one output column
            MotifsFileMode::Ungrouped => create_target(&mut labels, motif.label.clone())?,
            MotifsFileMode::Grouped => {
                // In grouped mode repeated group names intentionally share an output column
                let group_name = parse_group_name(columns[1], line_number)?;
                target_idx_for_group(&mut labels, &mut target_by_group, group_name)?
            }
        };

        rows.push(ParsedMotifsFileRow {
            motif,
            target_idx,
            line_number,
        });
    }

    let column_kind = match mode.context("--motifs-file must contain at least one motif")? {
        MotifsFileMode::Ungrouped => EndMotifColumnKind::Motif,
        MotifsFileMode::Grouped => EndMotifColumnKind::MotifGroup,
    };
    let (inside_spec, outside_spec) = build_selected_motif_half_specs(k_inside, k_outside, &rows)?;
    let mut lookup = FxHashMap::default();
    for row in &rows {
        // Store both observable end states for this public motif label
        for key in encoded_keys_for_motif(&row.motif, inside_spec.as_ref(), outside_spec.as_ref())?
        {
            insert_lookup_key(
                &mut lookup,
                key,
                row.target_idx,
                &row.motif.label,
                row.line_number,
            )?;
        }
    }

    Ok(SelectedEndMotifLookup {
        labels,
        column_kind,
        inside_spec,
        outside_spec,
        lookup,
    })
}

/// Determine whether one motifs-file row is grouped or ungrouped.
///
/// The parser requires the entire file to use one shape. Mixing one-column and two-column rows
/// would make the output axis ambiguous, so the caller compares each row with the first row's mode.
///
/// Parameters
/// ----------
/// - `columns`:
///   Tab-split row fields
/// - `line_number`:
///   One-based line number for diagnostics
///
/// Returns
/// -------
/// - `Result<MotifsFileMode>`:
///   Parsed row mode, or an error when the row has the wrong number of columns
fn parse_row_mode(columns: &[&str], line_number: usize) -> Result<MotifsFileMode> {
    match columns.len() {
        1 => Ok(MotifsFileMode::Ungrouped),
        2 => Ok(MotifsFileMode::Grouped),
        _ => bail!(
            "line {line_number}: expected one or two tab-separated columns in --motifs-file, got {}",
            columns.len()
        ),
    }
}

/// Parse and normalize one public motif label.
///
/// Input labels usually use `<outside>_<inside>`. The underscore may be omitted when exactly one
/// side has length zero. Output labels are always normalized back to `<outside>_<inside>` so the
/// writer and downstream readers see one stable representation.
///
/// Parameters
/// ----------
/// - `raw_label`:
///   Motif label from the first motifs-file column
/// - `k_inside`:
///   Expected inside length
/// - `k_outside`:
///   Expected outside length
/// - `line_number`:
///   One-based line number for diagnostics
///
/// Returns
/// -------
/// - `Result<ParsedEndMotifLabel>`:
///   Uppercase, side-split motif ready for encoding
fn parse_motif_label(
    raw_label: &str,
    k_inside: usize,
    k_outside: usize,
    line_number: usize,
) -> Result<ParsedEndMotifLabel> {
    ensure!(
        !raw_label.is_empty(),
        "line {line_number}: motif label is empty"
    );
    ensure!(
        raw_label.matches('_').count() <= 1,
        "line {line_number}: motif `{raw_label}` has more than one `_` separator"
    );

    // Require the explicit separator only when both sides contain bases
    let (outside_raw, inside_raw) = match raw_label.split_once('_') {
        Some((outside, inside)) => (outside, inside),
        None if k_outside > 0 && k_inside > 0 => bail!(
            "line {line_number}: motif `{raw_label}` must use `<outside>_<inside>` because both --k-outside and --k-inside are non-zero"
        ),
        None if k_outside > 0 => (raw_label, ""),
        None => ("", raw_label),
    };

    validate_motif_half(outside_raw, k_outside, "outside", raw_label, line_number)?;
    validate_motif_half(inside_raw, k_inside, "inside", raw_label, line_number)?;

    // Store the canonical uppercase spelling regardless of input casing
    let outside = outside_raw.to_ascii_uppercase();
    let inside = inside_raw.to_ascii_uppercase();
    Ok(ParsedEndMotifLabel {
        label: format!("{outside}_{inside}"),
        outside,
        inside,
    })
}

/// Validate one side of a parsed motif label.
///
/// This checks length and allowed bases before the k-mer codec is used. The later codec checks are
/// still kept because they protect this module against future parser changes.
///
/// Parameters
/// ----------
/// - `value`:
///   Raw motif side from the file
/// - `expected_len`:
///   Required length for this side
/// - `side_name`:
///   User-facing side name for diagnostics
/// - `raw_label`:
///   Original motif label from the file
/// - `line_number`:
///   One-based line number for diagnostics
///
/// Returns
/// -------
/// - `Result<()>`:
///   `Ok(())` when the side is valid DNA with the configured length
fn validate_motif_half(
    value: &str,
    expected_len: usize,
    side_name: &str,
    raw_label: &str,
    line_number: usize,
) -> Result<()> {
    ensure!(
        value.len() == expected_len,
        "line {line_number}: motif `{raw_label}` has {side_name} length {}, expected {expected_len}",
        value.len()
    );
    for base in value.bytes() {
        ensure!(
            matches!(base, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'),
            "line {line_number}: motif `{raw_label}` contains invalid base `{}`",
            base as char
        );
    }
    Ok(())
}

/// Parse a grouped motifs-file target label.
///
/// Group names are public labels written into JSON metadata and downstream data frames. They must
/// be non-empty ASCII labels without whitespace. Dots, underscores, and hyphens are accepted so
/// ordinary analysis labels do not need rewriting.
///
/// Parameters
/// ----------
/// - `raw_group_name`:
///   Group label from the second motifs-file column
/// - `line_number`:
///   One-based line number for diagnostics
///
/// Returns
/// -------
/// - `Result<String>`:
///   Accepted group name in its original spelling
fn parse_group_name(raw_group_name: &str, line_number: usize) -> Result<String> {
    ensure!(
        !raw_group_name.is_empty(),
        "line {line_number}: group name is empty"
    );
    for character in raw_group_name.chars() {
        ensure!(
            character.is_ascii_alphanumeric()
                || character == '.'
                || character == '_'
                || character == '-',
            "line {line_number}: group name `{raw_group_name}` contains invalid character `{character}`"
        );
    }
    Ok(raw_group_name.to_string())
}

/// Append one new output target label and return its target index.
///
/// Target indices are internal numeric coordinates used during counting and tile reduction. They
/// are assigned in motifs-file order and later compacted if unobserved targets are dropped.
///
/// Parameters
/// ----------
/// - `labels`:
///   Output target labels accumulated so far
/// - `label`:
///   Label for the new target
///
/// Returns
/// -------
/// - `Result<u32>`:
///   Internal target index for the appended label
fn create_target(labels: &mut Vec<String>, label: String) -> Result<u32> {
    let target_idx =
        u32::try_from(labels.len()).context("motifs-file target index does not fit in u32")?;
    labels.push(label);
    Ok(target_idx)
}

/// Return the existing target index for a group or create it in first-seen order.
///
/// Grouped motifs intentionally allow multiple motifs to count into the same output column. The
/// first row mentioning a group defines its position on the output axis.
///
/// Parameters
/// ----------
/// - `labels`:
///   Output target labels accumulated so far
/// - `target_by_group`:
///   Map from group label to target index
/// - `group_name`:
///   Parsed group label for the current motif
///
/// Returns
/// -------
/// - `Result<u32>`:
///   Target index for the group
fn target_idx_for_group(
    labels: &mut Vec<String>,
    target_by_group: &mut FxHashMap<String, u32>,
    group_name: String,
) -> Result<u32> {
    if let Some(target_idx) = target_by_group.get(&group_name) {
        return Ok(*target_idx);
    }

    let target_idx = create_target(labels, group_name.clone())?;
    target_by_group.insert(group_name, target_idx);
    Ok(target_idx)
}

fn build_selected_motif_half_specs(
    k_inside: usize,
    k_outside: usize,
    rows: &[ParsedMotifsFileRow],
) -> Result<(Option<EndMotifHalfSpec>, Option<EndMotifHalfSpec>)> {
    if k_inside > MAX_RADIX5_KMER_SIZE && k_inside == k_outside && k_inside > 0 {
        // One byte-backed selected subspace is enough when both large-k halves have the same
        // length. The full inside/outside pair is still checked by the encoded lookup, so sharing
        // this half-code universe only broadens the cheap prefilter and avoids duplicate per-tile
        // code arrays.
        let mut selected_halves = Vec::with_capacity(rows.len() * 4);
        for row in rows {
            collect_selected_half_states(&mut selected_halves, &row.motif.inside);
            collect_selected_half_states(&mut selected_halves, &row.motif.outside);
        }
        let shared_spec = Arc::new(
            build_subspace_kmer_spec(k_inside, &selected_halves)
                .with_context(|| "building shared selected motif-half subspace")?,
        );
        return Ok((
            Some(EndMotifHalfSpec::from_shared_subspace(shared_spec.clone())),
            Some(EndMotifHalfSpec::from_shared_subspace(shared_spec)),
        ));
    }

    Ok((
        build_optional_selected_motif_half_spec(k_inside, "--k-inside", rows, |motif| {
            motif.inside.as_str()
        })?,
        build_optional_selected_motif_half_spec(k_outside, "--k-outside", rows, |motif| {
            motif.outside.as_str()
        })?,
    ))
}

fn build_optional_selected_motif_half_spec(
    k: usize,
    label: &str,
    rows: &[ParsedMotifsFileRow],
    motif_half: impl Fn(&ParsedEndMotifLabel) -> &str,
) -> Result<Option<EndMotifHalfSpec>> {
    if k == 0 {
        return Ok(None);
    }
    if k <= MAX_RADIX5_KMER_SIZE {
        // Small motifs-file halves keep the ordinary full radix-5 code space. Subspaces are only
        // needed when the requested k cannot be represented by the full radix-5 codec.
        return build_optional_kmer_spec(k, label)
            .map(|spec| spec.map(EndMotifHalfSpec::from_radix5));
    }

    let mut selected_halves = Vec::with_capacity(rows.len() * 2);
    for row in rows {
        collect_selected_half_states(&mut selected_halves, motif_half(&row.motif));
    }
    let spec = build_subspace_kmer_spec(k, &selected_halves)
        .with_context(|| format!("building selected {label} subspace"))?;
    Ok(Some(EndMotifHalfSpec::from_subspace(spec)))
}

fn collect_selected_half_states(selected_halves: &mut Vec<String>, half: &str) {
    selected_halves.push(half.to_string());
    selected_halves.push(rev_complement(half));
}

/// Build the encoded lookup keys produced by one public motif label.
///
/// The left key is the motif as written. The right key uses reverse complements for both halves and
/// sets `reverse_on_decode = true`, matching the state produced by right-end counting. These two
/// states must remain separate because users may assign them to different groups.
///
/// Parameters
/// ----------
/// - `motif`:
///   Parsed motif label with separated inside and outside halves
/// - `inside_spec`:
///   Codec spec for the inside half, or `None` when `k_inside = 0`
/// - `outside_spec`:
///   Codec spec for the outside half, or `None` when `k_outside = 0`
///
/// Returns
/// -------
/// - `Result<[EncodedEndMotifKey, 2 entries]>`:
///   Encoded left-end and right-end states for the motif
fn encoded_keys_for_motif(
    motif: &ParsedEndMotifLabel,
    inside_spec: Option<&EndMotifHalfSpec>,
    outside_spec: Option<&EndMotifHalfSpec>,
) -> Result<[EncodedEndMotifKey; 2]> {
    let left_key = EncodedEndMotifKey {
        inside_code: encode_optional_motif_half(
            &motif.inside,
            inside_spec,
            "inside",
            &motif.label,
        )?,
        outside_code: encode_optional_motif_half(
            &motif.outside,
            outside_spec,
            "outside",
            &motif.label,
        )?,
        reverse_on_decode: false,
    };

    // Right-end counting stores the reverse-complement state and marks it for decode reversal
    let right_inside = rev_complement(&motif.inside);
    let right_outside = rev_complement(&motif.outside);
    let right_key = EncodedEndMotifKey {
        inside_code: encode_optional_motif_half(
            &right_inside,
            inside_spec,
            "inside",
            &motif.label,
        )?,
        outside_code: encode_optional_motif_half(
            &right_outside,
            outside_spec,
            "outside",
            &motif.label,
        )?,
        reverse_on_decode: true,
    };

    Ok([left_key, right_key])
}

/// Encode one optional motif half with the shared k-mer codec.
///
/// A missing spec means this side has `k=0`, so only an empty motif half is valid. A present spec
/// encodes the uppercase DNA string and rejects codec sentinel values before they reach the lookup.
///
/// Parameters
/// ----------
/// - `motif_half`:
///   Parsed motif half to encode
/// - `spec`:
///   Codec spec for this side, or `None` when the side is disabled
/// - `side_name`:
///   User-facing side name for diagnostics
/// - `motif_label`:
///   Normalized full motif label for diagnostics
///
/// Returns
/// -------
/// - `Result<u64>`:
///   Encoded motif-half code, or zero for a disabled side
fn encode_optional_motif_half(
    motif_half: &str,
    spec: Option<&EndMotifHalfSpec>,
    side_name: &str,
    motif_label: &str,
) -> Result<u64> {
    let Some(spec) = spec else {
        ensure!(
            motif_half.is_empty(),
            "motif `{motif_label}` has {side_name} bases even though that side has k=0"
        );
        return Ok(0);
    };

    // The parser already checked length and bases, but sentinel checks keep the codec boundary safe
    let code = spec.encode_kmer_bytes(motif_half.as_bytes());
    ensure!(
        !spec.code_is_invalid(code),
        "motif `{motif_label}` has {side_name} length {}, expected {}",
        motif_half.len(),
        spec.k()
    );
    Ok(code)
}

/// Insert one encoded motif state into the selected-target lookup.
///
/// Duplicate encoded states are allowed only when they resolve to the same target. A duplicate that
/// points at another target would make counting order-dependent, so it fails during parsing.
///
/// Parameters
/// ----------
/// - `lookup`:
///   Encoded state to target-index map under construction
/// - `key`:
///   Encoded left-end or right-end state for the motif
/// - `target_idx`:
///   Target index assigned to the current motif row
/// - `motif_label`:
///   Normalized motif label for diagnostics
/// - `line_number`:
///   One-based line number for diagnostics
///
/// Returns
/// -------
/// - `Result<()>`:
///   `Ok(())` when the state was inserted or was already assigned to the same target
fn insert_lookup_key(
    lookup: &mut FxHashMap<EncodedEndMotifKey, u32>,
    key: EncodedEndMotifKey,
    target_idx: u32,
    motif_label: &str,
    line_number: usize,
) -> Result<()> {
    if let Some(existing_target_idx) = lookup.get(&key) {
        ensure!(
            *existing_target_idx == target_idx,
            "line {line_number}: motif `{motif_label}` produces an encoded key already assigned to another target"
        );
        return Ok(());
    }

    lookup.insert(key, target_idx);
    Ok(())
}

#[cfg(test)]
mod tests {
    include!("motifs_file_tests.rs");
}
