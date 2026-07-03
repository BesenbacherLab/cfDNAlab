//! Reference correction for end-motif output loaders.
//!
//! This module keeps the public API on `EndsOutput` and shares the correction
//! logic between dense and sparse selections. The correction divides end-motif
//! counts by the matching reference k-mer frequency multiplied by the number
//! of positive-frequency reference motifs in that correction row.

use crate::output_loaders::{
    DenseMatrix, EndMotifAxisKind, EndMotifCountSelection, EndMotifCountsData, EndMotifRowMetadata,
    EndMotifRowMode, EndMotifSparseCounts, EndMotifSparseEntry, EndsOutput, EndsSelector,
    OutputLoaderResult, RefKmerFrequencyData, RefKmerMotifAxisKind, RefKmerRowMetadata,
    RefKmerRowMode, RefKmerSparseFrequencyLookup, RefKmersOutput, WindowRow,
};
use anyhow::{Context, Result, bail, ensure};
use std::collections::{BTreeMap, BTreeSet};

/// Controls positive end-motif counts that cannot be corrected.
///
/// A motif is unsupported when the matching reference k-mer is absent from the
/// reference motif axis or has zero reference frequency in the matched
/// correction row. Missing sparse reference entries are treated as zero.
///
/// Fixed-shape Rust selections cannot drop unsupported cells without changing
/// row or motif axes. Use `KeepNaN` when downstream code should keep the shape
/// and mark unsupported positive counts as `NaN`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnsupportedReferencePolicy {
    /// Report an error if any positive count has no positive reference frequency.
    Error,
    /// Keep unsupported positive counts in the selected shape as `NaN`.
    KeepNaN,
}

/// Builder for selecting reference-corrected end-motif counts.
///
/// Start from `EndsOutput::select_corrected_counts(&ref_kmers)`, add the same
/// row and motif selectors used by `EndsOutput::select()`, and call `read()`.
/// The result is an `EndMotifCountSelection` whose dense or sparse storage mode
/// matches the end-motif source.
///
/// Correction uses motifs with positive reference frequency in each matched
/// reference row. It does not shrink that motif count to the motifs observed in
/// the end-motif sample or to motifs selected on this builder. This makes the
/// correction factor depend only on the reference row, so selecting a motif
/// gives the same value as selecting everything and filtering afterward.
///
/// Concrete end-motif labels are matched to reference k-mers by removing `_`.
/// For example, `AT_CG` is matched to `ATCG`. Motif-group outputs are matched
/// directly by group label. Both `cfdna ends` and `cfdna ref-kmers` write
/// forward-oriented motif labels.
#[derive(Debug, Clone)]
pub struct CorrectedEndMotifCountsSelector<'a> {
    ends: &'a EndsOutput,
    ref_kmers: &'a RefKmersOutput,
    ends_selector: EndsSelector<'a>,
    use_global_bias: bool,
    unsupported_reference_policy: UnsupportedReferencePolicy,
}

impl<'a> CorrectedEndMotifCountsSelector<'a> {
    /// Create a corrected selector that initially selects all rows and motifs.
    pub(crate) fn new(ends: &'a EndsOutput, ref_kmers: &'a RefKmersOutput) -> Self {
        Self {
            ends,
            ref_kmers,
            ends_selector: ends.select(),
            use_global_bias: false,
            unsupported_reference_policy: UnsupportedReferencePolicy::Error,
        }
    }

    /// Select generic output rows by zero-based row index.
    pub fn rows(mut self, row_indices: &[usize]) -> Self {
        self.ends_selector = self.ends_selector.rows(row_indices);
        self
    }

    /// Select window rows by zero-based window row index.
    ///
    /// `read()` returns an error if the end-motif output is not windowed.
    pub fn windows(mut self, window_indices: &[usize]) -> Self {
        self.ends_selector = self.ends_selector.windows(window_indices);
        self
    }

    /// Select grouped rows by zero-based group row index.
    ///
    /// `read()` returns an error if the end-motif output is not grouped.
    pub fn groups(mut self, group_indices: &[usize]) -> Self {
        self.ends_selector = self.ends_selector.groups(group_indices);
        self
    }

    /// Select grouped rows by group name.
    ///
    /// `read()` returns an error if the end-motif output is not grouped or any
    /// requested name is missing or duplicated.
    pub fn groups_by_name<S: AsRef<str>>(mut self, group_names: &[S]) -> Self {
        self.ends_selector = self.ends_selector.groups_by_name(group_names);
        self
    }

    /// Select motifs by zero-based motif index.
    pub fn motifs(mut self, motif_indices: &[usize]) -> Self {
        self.ends_selector = self.ends_selector.motifs(motif_indices);
        self
    }

    /// Select motifs by end-motif label or motif-group label.
    pub fn motifs_by_label<S: AsRef<str>>(mut self, motif_labels: &[S]) -> Self {
        self.ends_selector = self.ends_selector.motifs_by_label(motif_labels);
        self
    }

    /// Allow a global reference k-mer output to correct every end-motif row.
    ///
    /// By default, end-motif rows and reference k-mer rows must match exactly.
    /// Set this to `true` only when the reference output is global and that
    /// global composition should be broadcast to a windowed or grouped
    /// end-motif output. The option is accepted but unnecessary when both
    /// outputs are global, because the rows already match exactly.
    pub fn use_global_bias(mut self, use_global_bias: bool) -> Self {
        self.use_global_bias = use_global_bias;
        self
    }

    /// Set how unsupported positive end-motif counts are handled.
    pub fn unsupported_reference_policy(mut self, policy: UnsupportedReferencePolicy) -> Self {
        self.unsupported_reference_policy = policy;
        self
    }

    /// Read selected counts after applying reference correction.
    pub fn read(self) -> OutputLoaderResult<EndMotifCountSelection> {
        validate_reference_correction_motif_axes(self.ends, self.ref_kmers)?;
        validate_reference_correction_rows(self.ends, self.ref_kmers, self.use_global_bias)?;

        let selection = self.ends_selector.read()?;
        let reference_row_indices = selected_reference_row_indices(
            self.ends,
            &selection,
            self.ref_kmers,
            self.use_global_bias,
        )?;
        let reference_motifs = selected_reference_motifs(self.ends, self.ref_kmers, &selection)?;
        let support_counts_by_reference_row =
            reference_support_counts(self.ref_kmers, &reference_row_indices)?;
        let reference_frequencies = ReferenceFrequencies::new(self.ref_kmers);
        let corrected_data = corrected_counts_data(
            &selection,
            &reference_row_indices,
            &reference_motifs,
            &support_counts_by_reference_row,
            &reference_frequencies,
            self.unsupported_reference_policy,
        )?;
        Ok(selection.with_data(corrected_data)?)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum CorrectionRowKey {
    Global,
    Window {
        index: usize,
        chrom: String,
        start: u64,
        end: u64,
    },
    Group(String),
}

#[derive(Debug, Clone)]
struct ReferenceMotif {
    label: String,
    index: Option<usize>,
}

enum ReferenceFrequencies<'a> {
    Dense(&'a DenseMatrix<f64>),
    Sparse(RefKmerSparseFrequencyLookup),
}

impl<'a> ReferenceFrequencies<'a> {
    fn new(ref_kmers: &'a RefKmersOutput) -> Self {
        match ref_kmers.data() {
            RefKmerFrequencyData::Dense(frequencies) => Self::Dense(frequencies),
            RefKmerFrequencyData::Sparse(sparse) => Self::Sparse(sparse.to_lookup_index()),
        }
    }

    fn frequency(&self, row_index: usize, motif_index: usize) -> Option<f64> {
        match self {
            Self::Dense(frequencies) => frequencies.get(row_index, motif_index).copied(),
            Self::Sparse(lookup) => lookup.frequency(row_index, motif_index),
        }
    }
}

fn validate_reference_correction_motif_axes(
    ends: &EndsOutput,
    ref_kmers: &RefKmersOutput,
) -> Result<()> {
    match ends.motif_axis_kind() {
        EndMotifAxisKind::MotifGroup => {
            ensure!(
                ref_kmers.motif_axis_kind() == RefKmerMotifAxisKind::MotifGroup,
                "grouped end-motif output requires grouped reference k-mer output"
            );
        }
        EndMotifAxisKind::Motif => {
            ensure!(
                ref_kmers.motif_axis_kind() == RefKmerMotifAxisKind::Motif,
                "concrete end-motif output requires concrete reference k-mer output"
            );
            ensure!(
                !ref_kmers.canonical(),
                "reference correction requires non-canonical reference k-mer output"
            );
            let reference_kmer_size = usize::from(ref_kmers.kmer_size());
            if let Some(motif_label) = ends.motif_labels().iter().find(|motif_label| {
                reference_motif_label(ends, motif_label).len() != reference_kmer_size
            }) {
                bail!(
                    "end-motif width must match reference k-mer size ({reference_kmer_size}): {motif_label}"
                );
            }
        }
    }
    Ok(())
}

fn validate_reference_correction_rows(
    ends: &EndsOutput,
    ref_kmers: &RefKmersOutput,
    use_global_bias: bool,
) -> Result<()> {
    validate_global_bias_option(ref_kmers, use_global_bias)?;

    if use_global_reference_bias(ends, ref_kmers, use_global_bias) {
        return Ok(());
    }

    if !row_modes_match(ends.row_mode(), ref_kmers.row_mode()) {
        if ref_kmers.row_mode() == RefKmerRowMode::Global
            && ends.row_mode() != EndMotifRowMode::Global
        {
            bail!(
                "reference k-mer output is global but end-motif output is {}. Use use_global_bias(true) to apply the global reference bias to every end-motif row",
                describe_end_row_mode(ends.row_mode())
            );
        }
        bail!(
            "end-motif and reference k-mer row modes must match: {} != {}",
            describe_end_row_mode(ends.row_mode()),
            describe_ref_row_mode(ref_kmers.row_mode())
        );
    }

    let end_row_keys = all_end_row_keys(ends)?;
    let reference_row_keys = all_reference_row_keys(ref_kmers)?;
    ensure!(
        unique_sorted_keys(&end_row_keys).len() == end_row_keys.len(),
        "end-motif row labels are not unique enough for correction"
    );
    ensure!(
        unique_sorted_keys(&reference_row_keys).len() == reference_row_keys.len(),
        "reference k-mer row labels are not unique enough for correction"
    );
    ensure!(
        unique_sorted_keys(&end_row_keys) == unique_sorted_keys(&reference_row_keys),
        "end-motif and reference k-mer rows do not match. Run ref-kmers with the same windowing or grouping"
    );
    Ok(())
}

fn validate_global_bias_option(ref_kmers: &RefKmersOutput, use_global_bias: bool) -> Result<()> {
    ensure!(
        !use_global_bias || ref_kmers.row_mode() == RefKmerRowMode::Global,
        "use_global_bias(true) requires a global reference k-mer output"
    );
    Ok(())
}

fn row_modes_match(end_row_mode: EndMotifRowMode, ref_row_mode: RefKmerRowMode) -> bool {
    matches!(
        (end_row_mode, ref_row_mode),
        (EndMotifRowMode::Global, RefKmerRowMode::Global)
            | (EndMotifRowMode::SizeWindows, RefKmerRowMode::SizeWindows)
            | (EndMotifRowMode::BedWindows, RefKmerRowMode::BedWindows)
            | (EndMotifRowMode::Groups, RefKmerRowMode::Groups)
    )
}

fn use_global_reference_bias(
    ends: &EndsOutput,
    ref_kmers: &RefKmersOutput,
    use_global_bias: bool,
) -> bool {
    use_global_bias
        && ref_kmers.row_mode() == RefKmerRowMode::Global
        && ends.row_mode() != EndMotifRowMode::Global
}

fn selected_reference_row_indices(
    ends: &EndsOutput,
    selection: &EndMotifCountSelection,
    ref_kmers: &RefKmersOutput,
    use_global_bias: bool,
) -> Result<Vec<usize>> {
    if use_global_reference_bias(ends, ref_kmers, use_global_bias) {
        ensure!(
            ref_kmers.row_count() == 1,
            "global reference k-mer output must contain exactly one row"
        );
        return Ok(vec![0; selection.row_count()]);
    }

    let row_indices_by_reference_key = reference_indices_by_key(ref_kmers)?;
    selection
        .row_indices()
        .iter()
        .map(|&end_row_index| {
            let row_key = end_row_key(ends, end_row_index)?;
            row_indices_by_reference_key
                .get(&row_key)
                .copied()
                .with_context(|| "selected end-motif row has no matching reference k-mer row")
        })
        .collect()
}

fn selected_reference_motifs(
    ends: &EndsOutput,
    ref_kmers: &RefKmersOutput,
    selection: &EndMotifCountSelection,
) -> Result<Vec<ReferenceMotif>> {
    selection
        .motif_labels()
        .iter()
        .map(|motif_label| {
            let label = reference_motif_label(ends, motif_label);
            let index = match ref_kmers.motif_index(&label) {
                Ok(index) => Some(index),
                Err(_) => None,
            };
            Ok(ReferenceMotif { label, index })
        })
        .collect()
}

fn reference_support_counts(
    ref_kmers: &RefKmersOutput,
    reference_row_indices: &[usize],
) -> Result<BTreeMap<usize, usize>> {
    let requested_rows = reference_row_indices
        .iter()
        .copied()
        .collect::<BTreeSet<_>>();
    let mut support_counts = requested_rows
        .iter()
        .map(|&row_index| (row_index, 0usize))
        .collect::<BTreeMap<_, _>>();

    match ref_kmers.data() {
        RefKmerFrequencyData::Dense(frequencies) => {
            for &row_index in &requested_rows {
                let row_values = frequencies.row(row_index).with_context(|| {
                    format!("reference row index {row_index} is outside frequency matrix")
                })?;
                let positive_count = row_values
                    .iter()
                    .filter(|&&frequency| frequency > 0.0)
                    .count();
                support_counts.insert(row_index, positive_count);
            }
        }
        RefKmerFrequencyData::Sparse(sparse) => {
            for entry in sparse.entries() {
                if entry.frequency > 0.0 && requested_rows.contains(&entry.row_index) {
                    *support_counts.entry(entry.row_index).or_default() += 1;
                }
            }
        }
    }
    Ok(support_counts)
}

fn corrected_counts_data(
    selection: &EndMotifCountSelection,
    reference_row_indices: &[usize],
    reference_motifs: &[ReferenceMotif],
    reference_support_counts: &BTreeMap<usize, usize>,
    reference_frequencies: &ReferenceFrequencies<'_>,
    unsupported_reference_policy: UnsupportedReferencePolicy,
) -> Result<EndMotifCountsData> {
    match selection.data() {
        EndMotifCountsData::Dense(counts) => Ok(EndMotifCountsData::Dense(correct_dense_counts(
            counts,
            reference_row_indices,
            reference_motifs,
            reference_support_counts,
            reference_frequencies,
            unsupported_reference_policy,
        )?)),
        EndMotifCountsData::Sparse(sparse) => {
            Ok(EndMotifCountsData::Sparse(correct_sparse_counts(
                sparse,
                reference_row_indices,
                reference_motifs,
                reference_support_counts,
                reference_frequencies,
                unsupported_reference_policy,
            )?))
        }
    }
}

fn correct_dense_counts(
    counts: &DenseMatrix<f64>,
    reference_row_indices: &[usize],
    reference_motifs: &[ReferenceMotif],
    reference_support_counts: &BTreeMap<usize, usize>,
    reference_frequencies: &ReferenceFrequencies<'_>,
    unsupported_reference_policy: UnsupportedReferencePolicy,
) -> Result<DenseMatrix<f64>> {
    let (row_count, motif_count) = counts.shape();
    ensure!(
        reference_row_indices.len() == row_count,
        "reference row selection does not match end-motif row count"
    );
    ensure!(
        reference_motifs.len() == motif_count,
        "reference motif selection does not match end-motif motif count"
    );

    let mut unsupported_positive_labels = BTreeSet::new();
    let mut corrected_values = Vec::with_capacity(row_count.saturating_mul(motif_count));
    for (selected_row_index, &reference_row_index) in reference_row_indices.iter().enumerate() {
        let correction_motif_count = *reference_support_counts
            .get(&reference_row_index)
            .unwrap_or(&0);
        for (selected_motif_index, reference_motif) in reference_motifs.iter().enumerate() {
            let count = counts
                .get(selected_row_index, selected_motif_index)
                .copied()
                .with_context(|| {
                    format!(
                        "selected end-motif count coordinate ({selected_row_index}, {selected_motif_index}) is outside dense matrix"
                    )
                })?;
            let corrected_value = corrected_count(
                count,
                reference_row_index,
                correction_motif_count,
                reference_motif,
                reference_frequencies,
                unsupported_reference_policy,
                &mut unsupported_positive_labels,
            );
            corrected_values.push(corrected_value);
        }
    }
    ensure_no_unsupported_positive_counts(&unsupported_positive_labels)?;
    Ok(DenseMatrix::from_row_major(
        corrected_values,
        row_count,
        motif_count,
    )?)
}

fn correct_sparse_counts(
    counts: &EndMotifSparseCounts,
    reference_row_indices: &[usize],
    reference_motifs: &[ReferenceMotif],
    reference_support_counts: &BTreeMap<usize, usize>,
    reference_frequencies: &ReferenceFrequencies<'_>,
    unsupported_reference_policy: UnsupportedReferencePolicy,
) -> Result<EndMotifSparseCounts> {
    let (row_count, motif_count) = counts.shape();
    ensure!(
        reference_row_indices.len() == row_count,
        "reference row selection does not match sparse end-motif row count"
    );
    ensure!(
        reference_motifs.len() == motif_count,
        "reference motif selection does not match sparse end-motif motif count"
    );

    let mut unsupported_positive_labels = BTreeSet::new();
    let mut corrected_entries = Vec::new();
    for entry in counts.entries() {
        let reference_row_index =
            *reference_row_indices
                .get(entry.row_index)
                .with_context(|| {
                    format!(
                        "sparse end-motif row index {} is outside selected row count",
                        entry.row_index
                    )
                })?;
        let correction_motif_count = *reference_support_counts
            .get(&reference_row_index)
            .unwrap_or(&0);
        let reference_motif = reference_motifs.get(entry.motif_index).with_context(|| {
            format!(
                "sparse end-motif motif index {} is outside selected motif count",
                entry.motif_index
            )
        })?;
        let corrected_value = corrected_count(
            entry.count,
            reference_row_index,
            correction_motif_count,
            reference_motif,
            reference_frequencies,
            unsupported_reference_policy,
            &mut unsupported_positive_labels,
        );
        if corrected_value != 0.0 || corrected_value.is_nan() {
            corrected_entries.push(EndMotifSparseEntry {
                row_index: entry.row_index,
                motif_index: entry.motif_index,
                count: corrected_value,
            });
        }
    }
    ensure_no_unsupported_positive_counts(&unsupported_positive_labels)?;
    Ok(EndMotifSparseCounts::from_entries(
        row_count,
        motif_count,
        corrected_entries,
    ))
}

fn corrected_count(
    count: f64,
    reference_row_index: usize,
    correction_motif_count: usize,
    reference_motif: &ReferenceMotif,
    reference_frequencies: &ReferenceFrequencies<'_>,
    unsupported_reference_policy: UnsupportedReferencePolicy,
    unsupported_positive_labels: &mut BTreeSet<String>,
) -> f64 {
    let reference_frequency = reference_motif
        .index
        .and_then(|motif_index| reference_frequencies.frequency(reference_row_index, motif_index))
        .unwrap_or(0.0);
    let reference_scale = reference_frequency * correction_motif_count as f64;
    if reference_scale > 0.0 {
        return count / reference_scale;
    }
    if count > 0.0 {
        match unsupported_reference_policy {
            UnsupportedReferencePolicy::Error => {
                unsupported_positive_labels.insert(reference_motif.label.clone());
                0.0
            }
            UnsupportedReferencePolicy::KeepNaN => f64::NAN,
        }
    } else {
        0.0
    }
}

fn ensure_no_unsupported_positive_counts(
    unsupported_positive_labels: &BTreeSet<String>,
) -> Result<()> {
    if unsupported_positive_labels.is_empty() {
        return Ok(());
    }
    let labels = unsupported_positive_labels
        .iter()
        .map(String::as_str)
        .collect::<Vec<_>>()
        .join(", ");
    bail!(
        "positive-count end motifs have no positive reference frequency: {labels}. Use unsupported_reference_policy(UnsupportedReferencePolicy::KeepNaN) to keep them as NaN"
    );
}

fn reference_motif_label(ends: &EndsOutput, end_motif_label: &str) -> String {
    match ends.motif_axis_kind() {
        EndMotifAxisKind::Motif => end_motif_label.replace('_', ""),
        EndMotifAxisKind::MotifGroup => end_motif_label.to_string(),
    }
}

fn all_end_row_keys(ends: &EndsOutput) -> Result<Vec<CorrectionRowKey>> {
    match ends.row_metadata() {
        EndMotifRowMetadata::Global => Ok(vec![CorrectionRowKey::Global]),
        EndMotifRowMetadata::Windows { windows, .. } => windows
            .iter()
            .map(end_window_row_key)
            .collect::<Result<Vec<_>>>(),
        EndMotifRowMetadata::Groups(groups) => Ok(groups
            .iter()
            .map(|group| CorrectionRowKey::Group(group.name.clone()))
            .collect()),
    }
}

fn all_reference_row_keys(ref_kmers: &RefKmersOutput) -> Result<Vec<CorrectionRowKey>> {
    match ref_kmers.row_metadata() {
        RefKmerRowMetadata::Global => Ok(vec![CorrectionRowKey::Global]),
        RefKmerRowMetadata::Windows { windows, .. } => windows
            .iter()
            .map(reference_window_row_key)
            .collect::<Result<Vec<_>>>(),
        RefKmerRowMetadata::Groups(groups) => Ok(groups
            .iter()
            .map(|group| CorrectionRowKey::Group(group.name.clone()))
            .collect()),
    }
}

fn end_row_key(ends: &EndsOutput, row_index: usize) -> Result<CorrectionRowKey> {
    match ends.row_metadata() {
        EndMotifRowMetadata::Global => {
            ensure!(row_index == 0, "global end-motif row index must be 0");
            Ok(CorrectionRowKey::Global)
        }
        EndMotifRowMetadata::Windows { windows, .. } => {
            let window = windows.get(row_index).with_context(|| {
                format!("end-motif row index {row_index} is outside window metadata")
            })?;
            end_window_row_key(window)
        }
        EndMotifRowMetadata::Groups(groups) => {
            let group = groups.get(row_index).with_context(|| {
                format!("end-motif row index {row_index} is outside group metadata")
            })?;
            Ok(CorrectionRowKey::Group(group.name.clone()))
        }
    }
}

fn reference_indices_by_key(
    ref_kmers: &RefKmersOutput,
) -> Result<BTreeMap<CorrectionRowKey, usize>> {
    let mut indices_by_key = BTreeMap::new();
    for row_index in 0..ref_kmers.row_count() {
        let row_key = reference_row_key(ref_kmers, row_index)?;
        ensure!(
            indices_by_key.insert(row_key, row_index).is_none(),
            "reference k-mer row labels are not unique enough for correction"
        );
    }
    Ok(indices_by_key)
}

fn reference_row_key(ref_kmers: &RefKmersOutput, row_index: usize) -> Result<CorrectionRowKey> {
    match ref_kmers.row_metadata() {
        RefKmerRowMetadata::Global => {
            ensure!(row_index == 0, "global reference k-mer row index must be 0");
            Ok(CorrectionRowKey::Global)
        }
        RefKmerRowMetadata::Windows { windows, .. } => {
            let window = windows.get(row_index).with_context(|| {
                format!("reference row index {row_index} is outside window metadata")
            })?;
            reference_window_row_key(window)
        }
        RefKmerRowMetadata::Groups(groups) => {
            let group = groups.get(row_index).with_context(|| {
                format!("reference row index {row_index} is outside group metadata")
            })?;
            Ok(CorrectionRowKey::Group(group.name.clone()))
        }
    }
}

fn end_window_row_key(window: &WindowRow) -> Result<CorrectionRowKey> {
    window_row_key(window)
}

fn reference_window_row_key(window: &WindowRow) -> Result<CorrectionRowKey> {
    window_row_key(window)
}

fn window_row_key(window: &WindowRow) -> Result<CorrectionRowKey> {
    let (start, end) = window.interval.as_tuple();
    Ok(CorrectionRowKey::Window {
        index: window.index,
        chrom: window.chrom.clone(),
        start,
        end,
    })
}

fn unique_sorted_keys(row_keys: &[CorrectionRowKey]) -> Vec<CorrectionRowKey> {
    row_keys
        .iter()
        .cloned()
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

fn describe_end_row_mode(row_mode: EndMotifRowMode) -> &'static str {
    match row_mode {
        EndMotifRowMode::Global => "global",
        EndMotifRowMode::SizeWindows => "size windows",
        EndMotifRowMode::BedWindows => "BED windows",
        EndMotifRowMode::Groups => "groups",
    }
}

fn describe_ref_row_mode(row_mode: RefKmerRowMode) -> &'static str {
    match row_mode {
        RefKmerRowMode::Global => "global",
        RefKmerRowMode::SizeWindows => "size windows",
        RefKmerRowMode::BedWindows => "BED windows",
        RefKmerRowMode::Groups => "groups",
    }
}
