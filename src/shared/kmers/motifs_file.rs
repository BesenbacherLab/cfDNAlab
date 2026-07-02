//! Parser and lookup builder for selected k-mer motif files.
//!
//! `ends --motifs-file` and `ref-kmers --motifs-file` use the same TSV shape:
//!
//! - one column: one concrete motif per output target
//! - two columns: motif plus a group label, where repeated groups share an output target
//!
//! The command-specific part is the motif label parser. `ends` labels are normalized to
//! `<outside>_<inside>`, while `ref-kmers` labels are normalized to a plain k-mer. A ref-kmers row
//! may also use one `_` separator, for example `AC_GT`, and that separator is removed before
//! validation so the counted k-mer is `ACGT`.
//!
//! The parser validates every public label, builds the output target axis, and precomputes encoded
//! lookup keys before tile processing starts. Ungrouped files keep motif-file order. Grouped files
//! sort group labels alphabetically. Counting can then do one hash lookup per observed motif and
//! skip all unselected motifs without allocating motif strings.

use crate::shared::{
    base::rev_complement,
    kmers::kmer_codec::{
        KmerCodes, KmerSpec, MAX_RADIX5_KMER_SIZE, SubspaceKmerSpec,
        build_left_aligned_codes_for_spec, build_optional_kmer_spec, build_subspace_kmer_spec,
    },
};
use anyhow::{Context, Result, bail, ensure};
use fxhash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, fs, path::Path, sync::Arc};

/// Encoded key for one selected motif state before final decoding.
///
/// The key can represent end motifs with separate inside and outside halves, or ref-kmers with the
/// full k-mer stored as `inside_code` and `outside_code = 0`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub(crate) struct EncodedMotifKey {
    pub(crate) inside_code: u64,
    pub(crate) outside_code: u64,
    pub(crate) reverse_on_decode: bool,
}

/// Meaning of the public `motif` axis for selected motif outputs.
///
/// Output matrices use numeric motif-column coordinates internally. This enum records whether
/// those numeric columns should be interpreted as concrete motif labels or motifs-file groups.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SelectedMotifColumnKind {
    /// Each column is one concrete motif
    Motif,
    /// Each column is one user-defined group from the second motifs-file column
    MotifGroup,
}

/// Precomputed motifs-file lookup used by selected counting paths.
///
/// A motifs file defines two related things:
/// - The output target axis, stored in `labels`
/// - The encoded motif states that should count into each target, stored in `lookup`
///
/// `lookup` is keyed by the same encoded state that tile-local counting builds.
#[derive(Debug, Clone)]
pub(crate) struct SelectedMotifLookup {
    /// Output labels in parser-assigned target order
    pub(crate) labels: Vec<String>,
    /// Whether `labels` are concrete motifs or user-defined motif groups
    pub(crate) column_kind: SelectedMotifColumnKind,
    /// Codec spec for the inside half, or the full ref-kmer, if present
    pub(crate) inside_spec: Option<SelectedMotifHalfSpec>,
    /// Codec spec for the outside half, if present
    pub(crate) outside_spec: Option<SelectedMotifHalfSpec>,
    /// Encoded motif key to original target index in `labels`
    pub(crate) lookup: FxHashMap<EncodedMotifKey, u32>,
}

impl SelectedMotifLookup {
    /// Return the motifs-file target for an encoded motif state.
    ///
    /// This is intentionally a thin map lookup. Counting without a motifs file never calls it, and
    /// motifs-file counting has already paid the parsing and validation cost before tile processing
    /// starts.
    #[inline]
    pub(crate) fn target_for(&self, key: EncodedMotifKey) -> Option<u32> {
        self.lookup.get(&key).copied()
    }
}

/// Codec used for one selected motif half during tile-local counting.
///
/// Full radix-5 specs are used up to the radix-5 limit. Larger selected motifs switch to a
/// byte-backed selected subspace because the full motif universe cannot be represented.
#[derive(Clone, Debug)]
pub(crate) enum SelectedMotifHalfSpec {
    /// Full radix-5 k-mer space
    Radix5(Arc<KmerSpec>),
    /// Byte-backed selected-k-mer subspace for motifs-file halves above the radix-5 limit
    Subspace(Arc<SubspaceKmerSpec>),
}

impl SelectedMotifHalfSpec {
    /// Wrap a full-space radix-5 spec in the selected motif codec enum.
    pub(crate) fn from_radix5(spec: KmerSpec) -> Self {
        SelectedMotifHalfSpec::Radix5(Arc::new(spec))
    }

    /// Wrap a byte-backed selected-subspace spec in the selected motif codec enum.
    pub(crate) fn from_subspace(spec: SubspaceKmerSpec) -> Self {
        SelectedMotifHalfSpec::Subspace(Arc::new(spec))
    }

    /// Wrap a shared byte-backed selected-subspace spec in the selected motif codec enum.
    pub(crate) fn from_shared_subspace(spec: Arc<SubspaceKmerSpec>) -> Self {
        SelectedMotifHalfSpec::Subspace(spec)
    }

    /// Return the motif-half length.
    #[inline]
    pub(crate) fn k(&self) -> usize {
        match self {
            SelectedMotifHalfSpec::Radix5(spec) => spec.k,
            SelectedMotifHalfSpec::Subspace(spec) => spec.k,
        }
    }

    /// Encode one exact motif-half byte slice.
    #[inline]
    pub(crate) fn encode_kmer_bytes(&self, seq: &[u8]) -> u64 {
        match self {
            SelectedMotifHalfSpec::Radix5(spec) => spec.encode_kmer_bytes(seq),
            SelectedMotifHalfSpec::Subspace(spec) => spec.encode_kmer_bytes(seq),
        }
    }

    /// Build per-position codes for a tile-local reference slice.
    #[inline]
    pub(crate) fn build_left_aligned_codes(&self, seq: &[u8]) -> KmerCodes {
        match self {
            SelectedMotifHalfSpec::Radix5(spec) => build_left_aligned_codes_for_spec(seq, spec),
            SelectedMotifHalfSpec::Subspace(spec) => spec.build_left_aligned_codes(seq),
        }
    }

    /// Return whether an encoded code is invalid for this motif half.
    #[inline]
    pub(crate) fn code_is_invalid(&self, code: u64) -> bool {
        match self {
            SelectedMotifHalfSpec::Radix5(spec) => {
                code == spec.sentinel_none() || code == spec.sentinel_n()
            }
            SelectedMotifHalfSpec::Subspace(spec) => code == spec.sentinel_missing(),
        }
    }

    /// Return the invalid code used when a reference-coordinate lookup has no full motif half.
    #[inline]
    pub(crate) fn missing_reference_code(&self) -> u64 {
        match self {
            SelectedMotifHalfSpec::Radix5(spec) => spec.sentinel_none(),
            SelectedMotifHalfSpec::Subspace(spec) => spec.sentinel_missing(),
        }
    }

    /// Return the invalid code used when blacklist masking touches a reference motif half.
    #[inline]
    pub(crate) fn masked_reference_code(&self) -> u64 {
        match self {
            SelectedMotifHalfSpec::Radix5(spec) => spec.sentinel_n(),
            SelectedMotifHalfSpec::Subspace(spec) => spec.sentinel_missing(),
        }
    }

    /// Return whether two specs can share one precomputed reference-code vector.
    #[inline]
    pub(crate) fn can_share_reference_codes_with(&self, other: &Self) -> bool {
        matches!(
            (self, other),
            (SelectedMotifHalfSpec::Radix5(left), SelectedMotifHalfSpec::Radix5(right)) if left.k == right.k
        ) || matches!(
            (self, other),
            (SelectedMotifHalfSpec::Subspace(left), SelectedMotifHalfSpec::Subspace(right)) if Arc::ptr_eq(left, right)
        )
    }
}

/// Command-specific interpretation of the first motifs-file column.
#[derive(Debug, Clone, Copy)]
enum SelectedMotifsFileKind {
    /// End-motif labels with separated outside and inside halves.
    EndMotifs {
        /// Expected inside motif length from `--k-inside`
        k_inside: usize,
        /// Expected outside motif length from `--k-outside`
        k_outside: usize,
    },
    /// Reference k-mer labels with one left-to-right k-mer.
    #[cfg(feature = "cmd_ref_kmers")]
    RefKmers {
        /// Expected k-mer length from `--kmer-size`
        kmer_size: usize,
    },
}

impl SelectedMotifsFileKind {
    /// Return the user-facing noun used in diagnostics.
    fn item_name(self) -> &'static str {
        match self {
            Self::EndMotifs { .. } => "motif",
            #[cfg(feature = "cmd_ref_kmers")]
            Self::RefKmers { .. } => "k-mer",
        }
    }
}

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

/// Parsed and normalized representation of one selected motif label.
///
/// The public `label` is command-specific. For `ends`, it is `<outside>_<inside>`. For
/// `ref-kmers`, it is the plain k-mer with any one separator removed. `outside` and `inside` are
/// stored separately because the encoded lookup uses the same two-half key type for both commands.
struct ParsedSelectedMotifLabel {
    /// Normalized public motif label
    label: String,
    /// Uppercase outside sequence, empty for ref-kmers and when `k_outside = 0`
    outside: String,
    /// Uppercase inside sequence, or the full k-mer for ref-kmers
    inside: String,
}

struct ParsedMotifsFileRow {
    motif: ParsedSelectedMotifLabel,
    target_label: String,
    line_number: usize,
}

/// Parse a selected end-motifs file from disk.
///
/// This validates user-facing motif labels, builds the selected output target axis, and precomputes
/// the encoded lookup used by counting. Filtering must use this encoded lookup before sparse count
/// insertion.
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
/// - `Result<SelectedMotifLookup>`:
///   Parsed output axis and encoded lookup ready for tile-local counting
pub(crate) fn parse_selected_end_motifs_file(
    path: &Path,
    k_inside: usize,
    k_outside: usize,
) -> Result<SelectedMotifLookup> {
    let contents = fs::read_to_string(path)
        .with_context(|| format!("reading motifs file {}", path.display()))?;
    parse_selected_motifs(
        &contents,
        SelectedMotifsFileKind::EndMotifs {
            k_inside,
            k_outside,
        },
    )
    .with_context(|| format!("parsing motifs file {}", path.display()))
}

/// Parse a selected reference k-mers file from disk.
///
/// A one-column file creates one output target per k-mer. A two-column file creates one output
/// target per distinct group label, ordered alphabetically by group name. Ref-kmer labels are
/// normalized to a plain uppercase k-mer. One `_` separator is accepted and removed, so `AC_GT`
/// and `ACGT` address the same selected k-mer.
///
/// Parameters
/// ----------
/// - `path`:
///   Motifs-file path supplied by the user
/// - `kmer_size`:
///   Required k-mer length from `--kmer-size`
///
/// Returns
/// -------
/// - `Result<SelectedMotifLookup>`:
///   Parsed output axis and encoded lookup ready for tile-local counting
#[cfg(feature = "cmd_ref_kmers")]
pub(crate) fn parse_selected_ref_kmers_file(
    path: &Path,
    kmer_size: usize,
) -> Result<SelectedMotifLookup> {
    let contents = fs::read_to_string(path)
        .with_context(|| format!("reading ref-kmers motifs file {}", path.display()))?;
    parse_selected_motifs(&contents, SelectedMotifsFileKind::RefKmers { kmer_size })
        .with_context(|| format!("parsing ref-kmers motifs file {}", path.display()))
}

fn parse_selected_motifs(
    contents: &str,
    kind: SelectedMotifsFileKind,
) -> Result<SelectedMotifLookup> {
    ensure!(
        !contents.is_empty(),
        "--motifs-file must contain at least one {}",
        kind.item_name()
    );
    #[cfg(feature = "cmd_ref_kmers")]
    if let SelectedMotifsFileKind::RefKmers { kmer_size } = kind {
        ensure!(kmer_size > 0, "`--kmer-size` must be positive");
    }

    let mut mode = None;
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
        let motif = parse_motif_label(columns[0], kind, line_number)?;
        ensure!(
            seen_motif_labels.insert(motif.label.clone()),
            "line {line_number}: duplicate {} `{}` in --motifs-file",
            kind.item_name(),
            motif.label
        );

        let target_label = match row_mode {
            // In ungrouped mode each motif creates exactly one output column
            MotifsFileMode::Ungrouped => motif.label.clone(),
            // In grouped mode repeated group names intentionally share an output column
            MotifsFileMode::Grouped => parse_group_name(columns[1], line_number)?,
        };

        rows.push(ParsedMotifsFileRow {
            motif,
            target_label,
            line_number,
        });
    }

    let file_mode = mode.context("--motifs-file must contain at least one motif")?;
    let column_kind = match file_mode {
        MotifsFileMode::Ungrouped => SelectedMotifColumnKind::Motif,
        MotifsFileMode::Grouped => SelectedMotifColumnKind::MotifGroup,
    };
    let labels = build_target_labels(file_mode, &rows);
    let target_idx_by_label = target_idx_by_label(&labels)?;
    let (inside_spec, outside_spec) = build_selected_motif_half_specs(kind, &rows)?;
    let mut lookup = FxHashMap::default();
    for row in &rows {
        let target_idx = target_idx_by_label
            .get(row.target_label.as_str())
            .copied()
            .with_context(|| {
                format!(
                    "line {}: motifs-file target `{}` was missing from the output axis",
                    row.line_number, row.target_label
                )
            })?;
        for key in encoded_keys_for_motif(
            &row.motif,
            kind,
            inside_spec.as_ref(),
            outside_spec.as_ref(),
        )? {
            insert_lookup_key(
                &mut lookup,
                key,
                target_idx,
                &row.motif.label,
                row.line_number,
            )?;
        }
    }

    Ok(SelectedMotifLookup {
        labels,
        column_kind,
        inside_spec,
        outside_spec,
        lookup,
    })
}

/// Build the public output labels for the selected motif targets.
///
/// Ungrouped files keep motif-file order. Grouped files sort by group name so the output axis does
/// not depend on where a group first appeared in the input.
fn build_target_labels(file_mode: MotifsFileMode, rows: &[ParsedMotifsFileRow]) -> Vec<String> {
    match file_mode {
        MotifsFileMode::Ungrouped => rows.iter().map(|row| row.target_label.clone()).collect(),
        MotifsFileMode::Grouped => {
            let labels = rows
                .iter()
                .map(|row| row.target_label.clone())
                .collect::<BTreeSet<_>>();
            labels.into_iter().collect()
        }
    }
}

fn target_idx_by_label(labels: &[String]) -> Result<FxHashMap<&str, u32>> {
    labels
        .iter()
        .enumerate()
        .map(|(target_idx, label)| {
            Ok((
                label.as_str(),
                u32::try_from(target_idx)
                    .context("motifs-file target index does not fit in u32")?,
            ))
        })
        .collect()
}

/// Determine whether one motifs-file row is grouped or ungrouped.
///
/// The parser requires the entire file to use one shape. Mixing one-column and two-column rows
/// would make the output axis ambiguous, so the caller compares each row with the first row's mode.
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

fn parse_motif_label(
    raw_label: &str,
    kind: SelectedMotifsFileKind,
    line_number: usize,
) -> Result<ParsedSelectedMotifLabel> {
    match kind {
        SelectedMotifsFileKind::EndMotifs {
            k_inside,
            k_outside,
        } => parse_end_motif_label(raw_label, k_inside, k_outside, line_number),
        #[cfg(feature = "cmd_ref_kmers")]
        SelectedMotifsFileKind::RefKmers { kmer_size } => {
            parse_ref_kmer_label(raw_label, kmer_size, line_number)
        }
    }
}

/// Parse and normalize one public end-motif label.
///
/// Input labels usually use `<outside>_<inside>`. The underscore may be omitted when exactly one
/// side has length zero. Output labels are always normalized back to `<outside>_<inside>` so the
/// writer and downstream readers see one stable representation.
fn parse_end_motif_label(
    raw_label: &str,
    k_inside: usize,
    k_outside: usize,
    line_number: usize,
) -> Result<ParsedSelectedMotifLabel> {
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
    Ok(ParsedSelectedMotifLabel {
        label: format!("{outside}_{inside}"),
        outside,
        inside,
    })
}

/// Parse and normalize one public reference k-mer label.
///
/// The public ref-kmers motif axis is a plain k-mer. One `_` separator is accepted so files can
/// reuse end-motif-style labels, but it is removed before length and base validation.
#[cfg(feature = "cmd_ref_kmers")]
fn parse_ref_kmer_label(
    raw_label: &str,
    kmer_size: usize,
    line_number: usize,
) -> Result<ParsedSelectedMotifLabel> {
    ensure!(
        !raw_label.is_empty(),
        "line {line_number}: k-mer label is empty"
    );
    ensure!(
        raw_label.matches('_').count() <= 1,
        "line {line_number}: k-mer `{raw_label}` has more than one `_` separator"
    );

    let full_kmer = match raw_label.split_once('_') {
        Some((outside, inside)) => format!("{outside}{inside}"),
        None => raw_label.to_string(),
    };
    validate_ref_kmer(&full_kmer, raw_label, kmer_size, line_number)?;

    let motif = full_kmer.to_ascii_uppercase();
    Ok(ParsedSelectedMotifLabel {
        label: motif.clone(),
        outside: String::new(),
        inside: motif,
    })
}

/// Validate one side of a parsed end-motif label.
///
/// This checks length and allowed bases before the k-mer codec is used. The later codec checks are
/// still kept because they protect this module against future parser changes.
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

#[cfg(feature = "cmd_ref_kmers")]
fn validate_ref_kmer(
    value: &str,
    raw_label: &str,
    expected_len: usize,
    line_number: usize,
) -> Result<()> {
    ensure!(
        value.len() == expected_len,
        "line {line_number}: k-mer `{raw_label}` has length {}, expected {expected_len}",
        value.len()
    );
    for base in value.bytes() {
        ensure!(
            matches!(base, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'),
            "line {line_number}: k-mer `{raw_label}` contains invalid base `{}`",
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

fn build_selected_motif_half_specs(
    kind: SelectedMotifsFileKind,
    rows: &[ParsedMotifsFileRow],
) -> Result<(Option<SelectedMotifHalfSpec>, Option<SelectedMotifHalfSpec>)> {
    match kind {
        SelectedMotifsFileKind::EndMotifs {
            k_inside,
            k_outside,
        } => build_selected_end_motif_half_specs(k_inside, k_outside, rows),
        #[cfg(feature = "cmd_ref_kmers")]
        SelectedMotifsFileKind::RefKmers { kmer_size } => {
            Ok((Some(build_selected_ref_kmer_spec(kmer_size, rows)?), None))
        }
    }
}

fn build_selected_end_motif_half_specs(
    k_inside: usize,
    k_outside: usize,
    rows: &[ParsedMotifsFileRow],
) -> Result<(Option<SelectedMotifHalfSpec>, Option<SelectedMotifHalfSpec>)> {
    if k_inside > MAX_RADIX5_KMER_SIZE && k_inside == k_outside && k_inside > 0 {
        // One byte-backed selected subspace is enough when both large-k halves have the same
        // length. The full inside/outside pair is still checked by the encoded lookup, so sharing
        // this half-code universe only broadens the cheap prefilter and avoids duplicate per-tile
        // code arrays
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
            Some(SelectedMotifHalfSpec::from_shared_subspace(
                shared_spec.clone(),
            )),
            Some(SelectedMotifHalfSpec::from_shared_subspace(shared_spec)),
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
    motif_half: impl Fn(&ParsedSelectedMotifLabel) -> &str,
) -> Result<Option<SelectedMotifHalfSpec>> {
    if k == 0 {
        return Ok(None);
    }
    if k <= MAX_RADIX5_KMER_SIZE {
        // Small motifs-file halves keep the ordinary full radix-5 code space. Subspaces are only
        // needed when the requested k cannot be represented by the full radix-5 codec
        return build_optional_kmer_spec(k, label)
            .map(|spec| spec.map(SelectedMotifHalfSpec::from_radix5));
    }

    let mut selected_halves = Vec::with_capacity(rows.len() * 2);
    for row in rows {
        collect_selected_half_states(&mut selected_halves, motif_half(&row.motif));
    }
    let spec = build_subspace_kmer_spec(k, &selected_halves)
        .with_context(|| format!("building selected {label} subspace"))?;
    Ok(Some(SelectedMotifHalfSpec::from_subspace(spec)))
}

#[cfg(feature = "cmd_ref_kmers")]
fn build_selected_ref_kmer_spec(
    kmer_size: usize,
    rows: &[ParsedMotifsFileRow],
) -> Result<SelectedMotifHalfSpec> {
    if kmer_size <= MAX_RADIX5_KMER_SIZE {
        return build_optional_kmer_spec(kmer_size, "kmer")?
            .map(SelectedMotifHalfSpec::from_radix5)
            .context("missing non-zero reference k-mer spec");
    }

    let selected_kmers: Vec<&str> = rows.iter().map(|row| row.motif.inside.as_str()).collect();
    let spec = build_subspace_kmer_spec(kmer_size, &selected_kmers)
        .with_context(|| "building selected reference k-mer subspace")?;
    Ok(SelectedMotifHalfSpec::from_subspace(spec))
}

fn collect_selected_half_states(selected_halves: &mut Vec<String>, half: &str) {
    selected_halves.push(half.to_string());
    selected_halves.push(rev_complement(half));
}

/// Build the encoded lookup keys produced by one public motif label.
///
/// End motifs get both observable end states. Reference k-mers get one left-to-right state with an
/// empty outside half.
fn encoded_keys_for_motif(
    motif: &ParsedSelectedMotifLabel,
    kind: SelectedMotifsFileKind,
    inside_spec: Option<&SelectedMotifHalfSpec>,
    outside_spec: Option<&SelectedMotifHalfSpec>,
) -> Result<Vec<EncodedMotifKey>> {
    match kind {
        SelectedMotifsFileKind::EndMotifs { .. } => {
            encoded_keys_for_end_motif(motif, inside_spec, outside_spec).map(Vec::from)
        }
        #[cfg(feature = "cmd_ref_kmers")]
        SelectedMotifsFileKind::RefKmers { .. } => Ok(vec![EncodedMotifKey {
            inside_code: encode_optional_motif_half(
                &motif.inside,
                inside_spec,
                "k-mer",
                &motif.label,
            )?,
            outside_code: 0,
            reverse_on_decode: false,
        }]),
    }
}

/// Build the encoded lookup keys produced by one public end-motif label.
///
/// The left key is the motif as written. The right key uses reverse complements for both halves and
/// sets `reverse_on_decode = true`, matching the state produced by right-end counting. These two
/// states must remain separate because users may assign them to different groups.
fn encoded_keys_for_end_motif(
    motif: &ParsedSelectedMotifLabel,
    inside_spec: Option<&SelectedMotifHalfSpec>,
    outside_spec: Option<&SelectedMotifHalfSpec>,
) -> Result<[EncodedMotifKey; 2]> {
    let left_key = EncodedMotifKey {
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
    let right_key = EncodedMotifKey {
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
fn encode_optional_motif_half(
    motif_half: &str,
    spec: Option<&SelectedMotifHalfSpec>,
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
fn insert_lookup_key(
    lookup: &mut FxHashMap<EncodedMotifKey, u32>,
    key: EncodedMotifKey,
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
