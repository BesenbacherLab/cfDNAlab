//! Reference correction for end-motif output loaders.
//!
//! This module keeps the public API on `EndsOutput` and shares the correction
//! logic between dense and sparse selections. For a reference row with `n`
//! positive motif frequencies, a uniform motif has frequency `1/n`. A motif's
//! correction factor is its reference frequency relative to that uniform
//! frequency, `reference_frequency / (1/n)`, or equivalently
//! `reference_frequency * n`. Corrected counts divide by this factor, so a
//! uniformly represented reference motif has factor 1 and leaves counts
//! unchanged.
//!
//! In `joint` mode, labels such as `AC_TG` match reference k-mer `ACTG` after
//! removing `_`. `split` keeps the full label but calculates outside and inside
//! factors separately and multiplies them. Side modes first sum sample counts
//! over the unused side, then correct and return labels such as `AC_` or `_TG`.

use crate::commands::ref_kmers::config::RefKmerOrientation;
use crate::output_loaders::{
    DenseMatrix, EndMotifAxisKind, EndMotifCountSelection, EndMotifCountsData, EndMotifRowMetadata,
    EndMotifRowMode, EndMotifSparseCounts, EndMotifSparseEntry, EndsOutput, EndsSelector,
    OutputLoaderResult, RefKmerFrequencyData, RefKmerMotifAxisKind, RefKmerRowMetadata,
    RefKmerRowMode, RefKmerSparseFrequencyLookup, RefKmersOutput, WindowRow,
};
use anyhow::{Context, Result, bail, ensure};
use std::collections::{BTreeMap, BTreeSet};

/// Chooses how two-sided end motifs are reference corrected.
///
/// The mode is required when loaded motif labels contain both outside and
/// inside bases, such as `AC_TG`. It is not accepted for one-sided outputs
/// such as `AC_` or `_TG`, where there is no two-sided choice to make.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TwoSidedCorrectionMode {
    /// Keep full labels such as `AC_TG` and correct each count with matching k-mer `ACTG`.
    Joint,
    /// Keep full labels, but calculate separate outside and inside correction factors.
    ///
    /// For `AC_TG`, separate correction factors are calculated for outside
    /// label `AC` and inside label `TG`, then multiplied.
    Split,
    /// Sum over inside bases and return corrected outside labels such as `AC_`.
    Outside,
    /// Sum over outside bases and return corrected inside labels such as `_TG`.
    Inside,
}

/// Controls positive end-motif counts that cannot be corrected.
///
/// An observed sample motif is unsupported when it has a positive count but no
/// positive correction factor under the selected mode. In `Joint` mode, the
/// factor is the matching reference k-mer frequency times the number of
/// positive reference k-mers in that row. In `Split`, the outside and inside
/// factors are calculated independently from aggregated side frequencies, and
/// the motif is supported only when their product is positive. `Outside` and
/// `Inside` use the corresponding side factor after sample counts are summed
/// by that side. Missing sparse reference entries have frequency zero.
///
/// Fixed-shape Rust selections cannot drop unsupported cells without changing
/// row or motif axes. Use `KeepNaN` when downstream code should keep the shape
/// and mark unsupported positive counts as `NaN`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnsupportedReferencePolicy {
    /// Report an error if a positive sample count has no positive correction factor.
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
/// Motif labels are matched to reference k-mers by removing `_`. For example,
/// `AT_CG` is matched to `ATCG` in `joint` mode. Motif-group outputs are
/// matched directly by group label and do not accept a two-sided mode. End-motif
/// labels are oriented from each fragment end inward. Reference correction therefore
/// requires `ref-kmers --orientation both`, whose frequencies average each reference
/// k-mer with its reverse complement. This assumes that left and right fragment ends contribute
/// equally within each output row, as expected for genomic windows of practical size.
///
/// For labels that contain both outside and inside bases, choose the mode that
/// matches the question:
///
/// - `Joint` keeps labels such as `AC_TG` and corrects each count with the
///   matching reference k-mer, such as `ACTG`.
/// - `Split` keeps labels such as `AC_TG`, but calculates the correction factor
///   from the two sides separately. For `AC_TG`, separate correction factors are
///   calculated for outside label `AC` and inside label `TG`. Those two
///   correction factors are multiplied, and the observed `AC_TG` count is
///   divided by that product. Use this when full two-sided motif labels remain in the
///   result, but the exact full reference k-mers are too sparse or the
///   correction should treat outside and inside sequence composition
///   separately.
/// - `Outside` returns labels such as `AC_`. For each outside label, all full
///   motif counts with that outside label are summed first. For example,
///   `AC_AA` and `AC_TG` both contribute to the `AC_` count. That summed count
///   is corrected using the outside label `AC`.
/// - `Inside` returns labels such as `_TG`. For each inside label, all full
///   motif counts with that inside label are summed first. For example,
///   `AA_TG` and `AC_TG` both contribute to the `_TG` count. That summed count
///   is corrected using the inside label `TG`.
///
/// For `Split`, `Outside`, and `Inside`, side-specific reference frequencies
/// are calculated from the loaded full-length reference k-mers. For example,
/// the outside frequency for `AC` is the sum of frequencies for loaded k-mers
/// with prefix `AC`, such as `ACTG` and `ACAA`. The inside frequency for `TG` is
/// the corresponding sum over loaded k-mers with suffix `TG`. Separate shorter
/// reference k-mer runs are not required.
///
/// A motifs file used for the reference output restricts these sums to the
/// k-mers in that file. Without a motifs file, all k-mers in the reference
/// output can contribute, including k-mers absent from the sample end-motif
/// output.
///
/// For `Outside` and `Inside`, repeated side labels are deduplicated in their
/// first loaded-motif occurrence order. The returned selection's motif labels
/// and motif indices describe that corrected side axis and its matrix columns.
#[derive(Debug, Clone)]
pub struct CorrectedEndMotifCountsSelector<'a> {
    ends: &'a EndsOutput,
    ref_kmers: &'a RefKmersOutput,
    ends_selector: EndsSelector<'a>,
    motif_selector: CorrectedMotifSelector,
    motif_selection_error: Option<String>,
    two_sided_correction_mode: Option<TwoSidedCorrectionMode>,
    use_global_bias: bool,
    unsupported_reference_policy: UnsupportedReferencePolicy,
}

impl<'a> CorrectedEndMotifCountsSelector<'a> {
    /// Create a corrected selector over an end-motif and reference output pair.
    ///
    /// All rows and motifs are initially selected. Exact row matching and errors
    /// for unsupported positive counts remain enabled until explicitly changed.
    pub(crate) fn new(ends: &'a EndsOutput, ref_kmers: &'a RefKmersOutput) -> Self {
        Self {
            ends,
            ref_kmers,
            ends_selector: ends.select(),
            motif_selector: CorrectedMotifSelector::All,
            motif_selection_error: None,
            two_sided_correction_mode: None,
            use_global_bias: false,
            unsupported_reference_policy: UnsupportedReferencePolicy::Error,
        }
    }

    /// Select output rows by zero-based index on the stored row axis.
    ///
    /// Row-mode validation and bounds checks are deferred to `read()`, matching
    /// the behavior of the underlying end-motif selector.
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

    /// Select motifs by zero-based index on the stored motif axis.
    ///
    /// Indices are valid for exact and split correction. Side-only correction
    /// rejects them because it constructs a different, aggregated motif axis.
    pub fn motifs(mut self, motif_indices: &[usize]) -> Self {
        self.set_motif_selector(
            CorrectedMotifSelector::Indices(motif_indices.to_vec()),
            "motifs",
        );
        self
    }

    /// Select motifs by label on the correction mode's output axis.
    ///
    /// Exact and split modes resolve stored end-motif or motif-group labels.
    /// Side-only modes resolve the derived labels such as `AC_` or `_TG`.
    pub fn motifs_by_label<S: AsRef<str>>(mut self, motif_labels: &[S]) -> Self {
        self.set_motif_selector(
            CorrectedMotifSelector::Labels(
                motif_labels
                    .iter()
                    .map(|motif_label| motif_label.as_ref().to_string())
                    .collect(),
            ),
            "motifs_by_label",
        );
        self
    }

    /// Record a motif selector, deferring the first selector conflict to `read()`.
    ///
    /// Builder setters return `Self`, not `Result`, so they cannot report a
    /// conflict immediately. Remembering the first conflict preserves normal
    /// method chaining while ensuring `read()` returns the useful error.
    fn set_motif_selector(
        &mut self,
        selector: CorrectedMotifSelector,
        selector_name: &'static str,
    ) {
        if let Some(previous_selector_name) = self.motif_selector.selector_name() {
            if self.motif_selection_error.is_none() {
                self.motif_selection_error = Some(format!(
                    "cannot combine {previous_selector_name}() and {selector_name}() on the motif axis"
                ));
            }
        } else {
            self.motif_selector = selector;
        }
    }

    /// Select the correction model for motifs with outside and inside bases.
    ///
    /// The mode is required when the loaded end-motif axis has both an outside
    /// and inside label, such as `AC_TG`. It is rejected for one-sided motifs
    /// such as `AC_` or `_TG`, and for motif-group outputs.
    pub fn two_sided_correction(mut self, mode: TwoSidedCorrectionMode) -> Self {
        self.two_sided_correction_mode = Some(mode);
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

    /// Set how positive counts without a positive reference denominator are handled.
    ///
    /// The default reports all affected motif labels as an error. The alternative
    /// keeps those coordinates as `NaN` while unsupported zero counts remain zero.
    pub fn unsupported_reference_policy(mut self, policy: UnsupportedReferencePolicy) -> Self {
        self.unsupported_reference_policy = policy;
        self
    }

    /// Validate the configured axes and return reference-corrected counts.
    ///
    /// Correction mode determines whether full labels are matched directly,
    /// corrected from separate side marginals, or aggregated onto a derived side
    /// axis before correction. Row compatibility is checked before counts are read.
    pub fn read(self) -> OutputLoaderResult<EndMotifCountSelection> {
        let correction_mode =
            resolve_correction_mode(self.ends, self.ref_kmers, self.two_sided_correction_mode)?;
        validate_reference_correction_rows(self.ends, self.ref_kmers, self.use_global_bias)?;

        match correction_mode {
            CorrectionMode::ExactLabel => self.read_exact_label_correction(),
            CorrectionMode::Split(shape) => self.read_split_correction(shape),
            CorrectionMode::Outside(shape) => self.read_side_correction(shape, SideMode::Outside),
            CorrectionMode::Inside(shape) => self.read_side_correction(shape, SideMode::Inside),
        }
    }

    /// Correct each selected motif with its matching full reference label.
    ///
    /// The correction factor is the matching reference frequency times the
    /// number of positive reference motifs in that row. Each selected count is
    /// divided by that factor.
    fn read_exact_label_correction(self) -> OutputLoaderResult<EndMotifCountSelection> {
        let selection = self.read_source_motif_selection()?;
        Ok(correct_exact_label_selection(
            self.ends,
            self.ref_kmers,
            selection,
            self.use_global_bias,
            self.unsupported_reference_policy,
        )?)
    }

    /// Correct full two-sided motifs with separate outside and inside factors.
    ///
    /// The returned axis remains the selected full-motif axis. Reference k-mer
    /// frequencies are summed by prefix and suffix separately, each side is
    /// compared with its own uniform frequency, and the matching outside and
    /// inside factors are multiplied for each full motif.
    fn read_split_correction(
        self,
        shape: TwoSidedMotifShape,
    ) -> OutputLoaderResult<EndMotifCountSelection> {
        let selection = self.read_source_motif_selection()?;
        let reference_row_indices = selected_reference_row_indices(
            self.ends,
            &selection,
            self.ref_kmers,
            self.use_global_bias,
        )?;
        let parsed_motifs = parse_end_motif_labels(selection.motif_labels(), shape)?;
        let reference_caches =
            side_reference_caches(self.ref_kmers, &reference_row_indices, shape)?;
        let corrected_data = split_corrected_counts_data(
            &selection,
            &reference_row_indices,
            &parsed_motifs,
            &reference_caches,
            self.unsupported_reference_policy,
        )?;
        Ok(selection.with_data(corrected_data)?)
    }

    /// Aggregate the full motif axis onto one side and correct those totals.
    ///
    /// Counts sharing a selected outside or inside label are summed per row. The
    /// returned selection replaces the source motif axis with derived side labels
    /// and uses matching marginal reference denominators.
    fn read_side_correction(
        self,
        shape: TwoSidedMotifShape,
        side_mode: SideMode,
    ) -> OutputLoaderResult<EndMotifCountSelection> {
        let requested_side_labels = self.side_mode_motif_labels()?;
        let selection = self.ends_selector.clone().read()?;
        let reference_row_indices = selected_reference_row_indices(
            self.ends,
            &selection,
            self.ref_kmers,
            self.use_global_bias,
        )?;
        let parsed_motifs = parse_end_motif_labels(selection.motif_labels(), shape)?;
        let side_axis = SideAxisSelection::new(&parsed_motifs, side_mode, requested_side_labels)?;
        let aggregated_data = aggregate_side_counts_data(&selection, &side_axis)?;
        let reference_caches =
            side_reference_caches(self.ref_kmers, &reference_row_indices, shape)?;
        let corrected_data = side_corrected_counts_data(
            &aggregated_data,
            &reference_row_indices,
            &side_axis.selected_labels,
            side_mode,
            &reference_caches,
            self.unsupported_reference_policy,
        )?;
        Ok(selection.with_derived_motif_axis(
            side_axis.selected_indices,
            side_axis.selected_labels,
            corrected_data,
        )?)
    }

    /// Read selected rows and motifs from the stored full-motif axis.
    ///
    /// Exact and split correction preserve this source axis. Outside and inside
    /// correction instead read the full axis first and apply label selection
    /// after constructing the derived side axis.
    fn read_source_motif_selection(&self) -> OutputLoaderResult<EndMotifCountSelection> {
        self.ensure_no_motif_selection_conflict()?;
        let selector = match &self.motif_selector {
            CorrectedMotifSelector::All => self.ends_selector.clone(),
            CorrectedMotifSelector::Indices(indices) => self.ends_selector.clone().motifs(indices),
            CorrectedMotifSelector::Labels(labels) => {
                self.ends_selector.clone().motifs_by_label(labels)
            }
        };
        selector.read()
    }

    /// Resolve the optional label selector used by a derived side axis.
    ///
    /// Label selection is copied for later side-axis construction. Stored motif
    /// indices are rejected because they do not identify aggregated side columns,
    /// and any earlier selector conflict is reported first.
    fn side_mode_motif_labels(&self) -> Result<Option<Vec<String>>> {
        self.ensure_no_motif_selection_conflict()?;
        match &self.motif_selector {
            CorrectedMotifSelector::All => Ok(None),
            CorrectedMotifSelector::Labels(labels) => Ok(Some(labels.clone())),
            CorrectedMotifSelector::Indices(_) => {
                bail!(
                    "motif index selectors are not supported for outside or inside reference correction. Use motif labels on the side-mode axis"
                );
            }
        }
    }

    /// Report a motif-selector conflict recorded during method chaining.
    ///
    /// Builder methods cannot return a `Result`, so the first conflict is stored
    /// and raised here before any output is read or corrected.
    fn ensure_no_motif_selection_conflict(&self) -> Result<()> {
        if let Some(selection_error) = &self.motif_selection_error {
            bail!("{selection_error}");
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum CorrectedMotifSelector {
    All,
    Indices(Vec<usize>),
    Labels(Vec<String>),
}

impl CorrectedMotifSelector {
    /// Return the builder method name represented by this selector.
    ///
    /// The name is used to explain conflicts between index- and label-based
    /// motif selection. Selecting the full axis has no method name.
    fn selector_name(&self) -> Option<&'static str> {
        match self {
            Self::All => None,
            Self::Indices(_) => Some("motifs"),
            Self::Labels(_) => Some("motifs_by_label"),
        }
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

enum ReferenceFrequencies<'a> {
    Dense(&'a DenseMatrix<f64>),
    Sparse(RefKmerSparseFrequencyLookup),
}

impl<'a> ReferenceFrequencies<'a> {
    /// Wrap dense or sparse reference frequencies behind a shared lookup API.
    ///
    /// Sparse output is indexed once here so correction code can use the same
    /// coordinate lookup regardless of the stored representation.
    fn new(ref_kmers: &'a RefKmersOutput) -> Self {
        match ref_kmers.data() {
            RefKmerFrequencyData::Dense(frequencies) => Self::Dense(frequencies),
            RefKmerFrequencyData::Sparse(sparse) => Self::Sparse(sparse.to_lookup_index()),
        }
    }

    /// Return the frequency at a reference row and motif coordinate.
    ///
    /// Missing sparse coordinates and out-of-range dense coordinates return
    /// `None`, allowing callers to treat absent reference support as zero.
    fn frequency(&self, row_index: usize, motif_index: usize) -> Option<f64> {
        match self {
            Self::Dense(frequencies) => frequencies.get(row_index, motif_index).copied(),
            Self::Sparse(lookup) => lookup.frequency(row_index, motif_index),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CorrectionMode {
    ExactLabel,
    Split(TwoSidedMotifShape),
    Outside(TwoSidedMotifShape),
    Inside(TwoSidedMotifShape),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct TwoSidedMotifShape {
    outside_width: usize,
    inside_width: usize,
}

impl TwoSidedMotifShape {
    /// Return the reference k-mer width covered by both motif sides.
    ///
    /// This is the required width of an underscore-free reference motif label.
    fn combined_width(self) -> usize {
        self.outside_width + self.inside_width
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SideMode {
    Outside,
    Inside,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ParsedEndMotif {
    label: String,
    outside: String,
    inside: String,
}

impl ParsedEndMotif {
    /// Format this motif for an outside- or inside-only output axis.
    ///
    /// The underscore is retained on the boundary so the derived label still
    /// states which side of the fragment end it represents.
    fn side_label(&self, side_mode: SideMode) -> String {
        match side_mode {
            SideMode::Outside => format!("{}_", self.outside),
            SideMode::Inside => format!("_{}", self.inside),
        }
    }
}

#[derive(Debug, Clone)]
struct SideAxisSelection {
    selected_indices: Vec<usize>,
    selected_labels: Vec<String>,
    source_motif_to_selected_index: Vec<Option<usize>>,
}

impl SideAxisSelection {
    /// Build a selected side axis and map source motifs onto its columns.
    ///
    /// Full-motif labels that share a side are mapped to the same output column.
    /// Optional label selection filters that derived axis without changing its
    /// original indices or first-occurrence order.
    fn new(
        parsed_motifs: &[ParsedEndMotif],
        side_mode: SideMode,
        requested_side_labels: Option<Vec<String>>,
    ) -> Result<Self> {
        let (full_side_labels, source_motif_to_full_side_index) =
            full_side_axis(parsed_motifs, side_mode);
        let selected_indices = match requested_side_labels.as_ref() {
            Some(labels) => resolve_side_label_indices(&full_side_labels, labels)?,
            None => (0..full_side_labels.len()).collect(),
        };
        let selected_labels = match requested_side_labels {
            Some(labels) => labels,
            None => full_side_labels.clone(),
        };
        let mut full_side_to_selected_index = vec![None; full_side_labels.len()];
        for (selected_index, &full_side_index) in selected_indices.iter().enumerate() {
            full_side_to_selected_index[full_side_index] = Some(selected_index);
        }
        let source_motif_to_selected_index = source_motif_to_full_side_index
            .into_iter()
            .map(|full_side_index| full_side_to_selected_index[full_side_index])
            .collect();
        Ok(Self {
            selected_indices,
            selected_labels,
            source_motif_to_selected_index,
        })
    }
}

#[derive(Debug, Clone)]
struct SideReferenceCache {
    outside_frequencies: BTreeMap<String, f64>,
    inside_frequencies: BTreeMap<String, f64>,
    outside_support_count: usize,
    inside_support_count: usize,
}

impl SideReferenceCache {
    /// Create empty side-frequency maps and support counts.
    ///
    /// Support counts remain zero until all full reference frequencies have
    /// been accumulated and `finalize_support_counts` is called.
    fn new() -> Self {
        Self {
            outside_frequencies: BTreeMap::new(),
            inside_frequencies: BTreeMap::new(),
            outside_support_count: 0,
            inside_support_count: 0,
        }
    }

    /// Add a full reference motif's frequency to both side marginals.
    ///
    /// The label must match the resolved combined width. Its prefix contributes
    /// to the outside total and its suffix contributes to the inside total.
    fn add_frequency(
        &mut self,
        reference_motif_label: &str,
        frequency: f64,
        shape: TwoSidedMotifShape,
    ) -> Result<()> {
        ensure!(
            reference_motif_label.len() == shape.combined_width(),
            "reference motif label must split into outside width {} and inside width {}: {reference_motif_label}",
            shape.outside_width,
            shape.inside_width
        );
        let outside = reference_motif_label[..shape.outside_width].to_string();
        let inside = reference_motif_label[shape.outside_width..].to_string();
        *self.outside_frequencies.entry(outside).or_insert(0.0) += frequency;
        *self.inside_frequencies.entry(inside).or_insert(0.0) += frequency;
        Ok(())
    }

    /// Count supported labels after all side frequencies are accumulated.
    ///
    /// Only strictly positive marginal frequencies define correction support,
    /// matching the exact-label normalization rule.
    fn finalize_support_counts(&mut self) {
        self.outside_support_count = self
            .outside_frequencies
            .values()
            .filter(|&&frequency| frequency > 0.0)
            .count();
        self.inside_support_count = self
            .inside_frequencies
            .values()
            .filter(|&&frequency| frequency > 0.0)
            .count();
    }

    /// Return the selected side's reference correction denominator.
    ///
    /// The marginal frequency is multiplied by the number of positive-frequency
    /// labels on that side, so uniform side composition has denominator one.
    /// Missing labels have a zero denominator.
    fn denominator(&self, side_label: &str, side_mode: SideMode) -> f64 {
        match side_mode {
            SideMode::Outside => {
                let outside = side_label.trim_end_matches('_');
                self.outside_frequencies
                    .get(outside)
                    .copied()
                    .unwrap_or(0.0)
                    * self.outside_support_count as f64
            }
            SideMode::Inside => {
                let inside = side_label.trim_start_matches('_');
                self.inside_frequencies.get(inside).copied().unwrap_or(0.0)
                    * self.inside_support_count as f64
            }
        }
    }

    /// Return the split-mode denominator for a full motif.
    ///
    /// Outside and inside marginal frequencies are normalized independently,
    /// then multiplied so the full motif axis can be retained.
    fn split_denominator(&self, parsed_motif: &ParsedEndMotif) -> f64 {
        let outside_denominator = self
            .outside_frequencies
            .get(&parsed_motif.outside)
            .copied()
            .unwrap_or(0.0)
            * self.outside_support_count as f64;
        let inside_denominator = self
            .inside_frequencies
            .get(&parsed_motif.inside)
            .copied()
            .unwrap_or(0.0)
            * self.inside_support_count as f64;
        outside_denominator * inside_denominator
    }
}

/// Validate the motif axes and resolve the requested correction strategy.
///
/// Grouped, empty, one-sided, and joint axes use exact-label correction.
/// Other two-sided modes retain the inferred side widths needed to calculate
/// separate or side-only reference denominators.
fn resolve_correction_mode(
    ends: &EndsOutput,
    ref_kmers: &RefKmersOutput,
    two_sided_correction_mode: Option<TwoSidedCorrectionMode>,
) -> Result<CorrectionMode> {
    validate_reference_correction_motif_axes(ends, ref_kmers, two_sided_correction_mode)?;

    if ends.motif_axis_kind() == EndMotifAxisKind::MotifGroup {
        return Ok(CorrectionMode::ExactLabel);
    }
    if ends.motif_labels().is_empty() {
        return Ok(CorrectionMode::ExactLabel);
    }

    let reference_kmer_size = usize::from(ref_kmers.kmer_size());
    let shape = infer_concrete_motif_shape(ends.motif_labels(), reference_kmer_size)?;
    validate_reference_labels_split_cleanly(ref_kmers, shape)?;
    if shape.outside_width == 0 || shape.inside_width == 0 {
        ensure!(
            two_sided_correction_mode.is_none(),
            "one-sided end-motif outputs do not accept two_sided_correction"
        );
        return Ok(CorrectionMode::ExactLabel);
    }

    match two_sided_correction_mode {
        Some(TwoSidedCorrectionMode::Joint) => Ok(CorrectionMode::ExactLabel),
        Some(TwoSidedCorrectionMode::Split) => Ok(CorrectionMode::Split(shape)),
        Some(TwoSidedCorrectionMode::Outside) => Ok(CorrectionMode::Outside(shape)),
        Some(TwoSidedCorrectionMode::Inside) => Ok(CorrectionMode::Inside(shape)),
        None => {
            bail!(
                "two-sided end-motif labels with both outside and inside bases require two_sided_correction(TwoSidedCorrectionMode::Joint, Split, Outside, or Inside)"
            );
        }
    }
}

/// Ensure sample and reference motif axes can be used for correction.
///
/// Reference output must include both sequence orientations and be non-canonical. Grouped axes
/// must match, and concrete end-motif labels must have the reference k-mer width. Two-sided mode
/// options are rejected for grouped axes before any counts are read.
fn validate_reference_correction_motif_axes(
    ends: &EndsOutput,
    ref_kmers: &RefKmersOutput,
    two_sided_correction_mode: Option<TwoSidedCorrectionMode>,
) -> Result<()> {
    ensure!(
        ref_kmers.orientation() == RefKmerOrientation::Both,
        "reference correction requires reference k-mer output generated with `--orientation both`"
    );
    ensure!(
        !ref_kmers.canonical(),
        "reference correction requires non-canonical reference k-mer output"
    );
    match ends.motif_axis_kind() {
        EndMotifAxisKind::MotifGroup => {
            ensure!(
                two_sided_correction_mode.is_none(),
                "motif-group end-motif outputs do not accept two_sided_correction"
            );
            ensure!(
                ref_kmers.motif_axis_kind() == RefKmerMotifAxisKind::MotifGroup,
                "grouped end-motif output requires grouped reference k-mer output"
            );
        }
        EndMotifAxisKind::Motif => {
            ensure!(
                ref_kmers.motif_axis_kind() == RefKmerMotifAxisKind::Motif,
                "end-motif output with motif labels requires reference k-mer output with motif labels"
            );
            let reference_kmer_size = usize::from(ref_kmers.kmer_size());
            if !ends.motif_labels().is_empty() {
                infer_concrete_motif_shape(ends.motif_labels(), reference_kmer_size)?;
            }
        }
    }
    Ok(())
}

/// Infer the outside and inside widths shared by concrete motif labels.
///
/// Every label must contain one separator, have the reference k-mer's total
/// width, and use the same split. An empty axis cannot define a shape.
fn infer_concrete_motif_shape(
    motif_labels: &[String],
    reference_kmer_size: usize,
) -> Result<TwoSidedMotifShape> {
    let mut inferred_shape = None;
    for motif_label in motif_labels {
        let (outside, inside) = split_end_motif_label(motif_label)?;
        let shape = TwoSidedMotifShape {
            outside_width: outside.len(),
            inside_width: inside.len(),
        };
        ensure!(
            shape.combined_width() == reference_kmer_size,
            "end-motif width must match reference k-mer size ({reference_kmer_size}): {motif_label}"
        );
        match inferred_shape {
            Some(previous_shape) => ensure!(
                previous_shape == shape,
                "all end-motif labels must use the same outside and inside widths"
            ),
            None => inferred_shape = Some(shape),
        }
    }
    inferred_shape.with_context(|| "cannot infer side widths from an empty motif axis")
}

/// Split an end-motif label at its outside/inside separator.
///
/// Exactly one underscore is required. Either returned side may be empty for a
/// valid one-sided motif label.
fn split_end_motif_label(motif_label: &str) -> Result<(&str, &str)> {
    let mut parts = motif_label.split('_');
    let outside = parts.next().with_context(|| {
        format!("end-motif label must contain exactly one '_' to separate outside and inside bases: {motif_label}")
    })?;
    let inside = parts.next().with_context(|| {
        format!("end-motif label must contain exactly one '_' to separate outside and inside bases: {motif_label}")
    })?;
    ensure!(
        parts.next().is_none(),
        "end-motif label must contain exactly one '_' to separate outside and inside bases: {motif_label}"
    );
    Ok((outside, inside))
}

/// Parse end-motif labels and verify that they match the inferred shape.
///
/// The returned values retain the original label and owned side strings so
/// later aggregation does not need to split or validate labels again.
fn parse_end_motif_labels(
    motif_labels: &[String],
    shape: TwoSidedMotifShape,
) -> Result<Vec<ParsedEndMotif>> {
    motif_labels
        .iter()
        .map(|motif_label| {
            let (outside, inside) = split_end_motif_label(motif_label)?;
            ensure!(
                outside.len() == shape.outside_width && inside.len() == shape.inside_width,
                "end-motif label does not match inferred outside and inside widths: {motif_label}"
            );
            Ok(ParsedEndMotif {
                label: motif_label.clone(),
                outside: outside.to_string(),
                inside: inside.to_string(),
            })
        })
        .collect()
}

/// Ensure reference labels can be split using the resolved side widths.
///
/// Reference labels contain no underscore, so their byte length must equal the
/// combined outside and inside width.
fn validate_reference_labels_split_cleanly(
    ref_kmers: &RefKmersOutput,
    shape: TwoSidedMotifShape,
) -> Result<()> {
    for reference_motif_label in ref_kmers.motif_labels() {
        ensure!(
            reference_motif_label.len() == shape.combined_width(),
            "reference motif label must split into outside width {} and inside width {}: {reference_motif_label}",
            shape.outside_width,
            shape.inside_width
        );
    }
    Ok(())
}

/// Build the full side-label axis and map each source motif onto it.
///
/// Labels are deduplicated in first-occurrence order. The parallel mapping says
/// which derived side column receives each full-motif count.
fn full_side_axis(
    parsed_motifs: &[ParsedEndMotif],
    side_mode: SideMode,
) -> (Vec<String>, Vec<usize>) {
    let mut seen_labels = BTreeMap::<String, usize>::new();
    let mut side_labels = Vec::new();
    let mut source_motif_to_full_side_index = Vec::with_capacity(parsed_motifs.len());
    for parsed_motif in parsed_motifs {
        let side_label = parsed_motif.side_label(side_mode);
        let side_index = match seen_labels.get(&side_label).copied() {
            Some(side_index) => side_index,
            None => {
                let side_index = side_labels.len();
                seen_labels.insert(side_label.clone(), side_index);
                side_labels.push(side_label);
                side_index
            }
        };
        source_motif_to_full_side_index.push(side_index);
    }
    (side_labels, source_motif_to_full_side_index)
}

/// Resolve requested side labels to full side-axis indices.
///
/// Requested order is preserved. Duplicate requests and labels absent from the
/// derived axis are reported as errors instead of being silently ignored.
fn resolve_side_label_indices(
    full_side_labels: &[String],
    requested_labels: &[String],
) -> Result<Vec<usize>> {
    let mut seen_requested_labels = BTreeSet::new();
    let full_side_indices = full_side_labels
        .iter()
        .enumerate()
        .map(|(side_index, label)| (label.as_str(), side_index))
        .collect::<BTreeMap<_, _>>();
    requested_labels
        .iter()
        .map(|requested_label| {
            ensure!(
                seen_requested_labels.insert(requested_label.as_str()),
                "motif label selector contains duplicate side-mode label '{requested_label}'"
            );
            full_side_indices
                .get(requested_label.as_str())
                .copied()
                .with_context(|| format!("side-mode motif axis has no label '{requested_label}'"))
        })
        .collect()
}

/// Validate that sample and reference rows can be matched for correction.
///
/// Without global broadcasting, row modes must agree and both complete key sets
/// must be unique and identical. A global reference may be broadcast only when
/// explicitly requested.
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
    let unique_end_row_keys = end_row_keys.iter().cloned().collect::<BTreeSet<_>>();
    let unique_reference_row_keys = reference_row_keys.iter().cloned().collect::<BTreeSet<_>>();
    ensure!(
        unique_end_row_keys.len() == end_row_keys.len(),
        "end-motif row labels are not unique enough for correction"
    );
    ensure!(
        unique_reference_row_keys.len() == reference_row_keys.len(),
        "reference k-mer row labels are not unique enough for correction"
    );
    ensure!(
        unique_end_row_keys == unique_reference_row_keys,
        "end-motif and reference k-mer rows do not match. Run ref-kmers with the same windowing or grouping"
    );
    Ok(())
}

/// Ensure global-bias broadcasting is used only with a global reference.
///
/// This catches an invalid option even when the sample and reference row modes
/// would otherwise happen to match.
fn validate_global_bias_option(ref_kmers: &RefKmersOutput, use_global_bias: bool) -> Result<()> {
    ensure!(
        !use_global_bias || ref_kmers.row_mode() == RefKmerRowMode::Global,
        "use_global_bias(true) requires a global reference k-mer output"
    );
    Ok(())
}

/// Return whether the end and reference row modes describe the same row axis.
///
/// Window kinds are matched explicitly, so size windows and BED windows are not
/// considered interchangeable even though both carry interval metadata.
fn row_modes_match(end_row_mode: EndMotifRowMode, ref_row_mode: RefKmerRowMode) -> bool {
    matches!(
        (end_row_mode, ref_row_mode),
        (EndMotifRowMode::Global, RefKmerRowMode::Global)
            | (EndMotifRowMode::SizeWindows, RefKmerRowMode::SizeWindows)
            | (EndMotifRowMode::BedWindows, RefKmerRowMode::BedWindows)
            | (EndMotifRowMode::Groups, RefKmerRowMode::Groups)
    )
}

/// Return whether one global reference row should be applied to every end row.
///
/// Broadcasting is active only for a non-global sample and an explicitly
/// enabled global reference.
fn use_global_reference_bias(
    ends: &EndsOutput,
    ref_kmers: &RefKmersOutput,
    use_global_bias: bool,
) -> bool {
    use_global_bias
        && ref_kmers.row_mode() == RefKmerRowMode::Global
        && ends.row_mode() != EndMotifRowMode::Global
}

/// Correct selected counts by matching full sample and reference motif labels.
///
/// The function resolves reference rows and motifs, counts positive reference
/// support per row, and divides counts by frequency times support. The selected
/// row and motif axes are preserved.
fn correct_exact_label_selection(
    ends: &EndsOutput,
    ref_kmers: &RefKmersOutput,
    selection: EndMotifCountSelection,
    use_global_bias: bool,
    unsupported_reference_policy: UnsupportedReferencePolicy,
) -> Result<EndMotifCountSelection> {
    let reference_row_indices =
        selected_reference_row_indices(ends, &selection, ref_kmers, use_global_bias)?;
    let reference_motif_indices = selected_reference_motif_indices(ends, ref_kmers, &selection);
    let support_counts_by_reference_row =
        reference_support_counts(ref_kmers, &reference_row_indices)?;
    let reference_frequencies = ReferenceFrequencies::new(ref_kmers);
    let corrected_data = corrected_counts_data(
        &selection,
        &reference_row_indices,
        &reference_motif_indices,
        &support_counts_by_reference_row,
        &reference_frequencies,
        unsupported_reference_policy,
    )?;
    Ok(selection.with_data(corrected_data)?)
}

/// Match selected end rows to corresponding reference row indices.
///
/// Global broadcasting repeats row zero. Keyed correction resolves every end
/// row through its full correction key and fails if any reference row is absent.
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

/// Build side-frequency caches for requested reference rows.
///
/// Dense and sparse reference representations contribute the same prefix and
/// suffix marginals. Positive support counts are finalized only after every
/// frequency in a requested row has been accumulated.
fn side_reference_caches(
    ref_kmers: &RefKmersOutput,
    reference_row_indices: &[usize],
    shape: TwoSidedMotifShape,
) -> Result<BTreeMap<usize, SideReferenceCache>> {
    let requested_rows = reference_row_indices
        .iter()
        .copied()
        .collect::<BTreeSet<_>>();
    let mut caches = requested_rows
        .iter()
        .map(|&row_index| (row_index, SideReferenceCache::new()))
        .collect::<BTreeMap<_, _>>();

    match ref_kmers.data() {
        RefKmerFrequencyData::Dense(frequencies) => {
            for &row_index in &requested_rows {
                let row_values = frequencies.row(row_index).with_context(|| {
                    format!("reference row index {row_index} is outside frequency matrix")
                })?;
                let cache = caches
                    .get_mut(&row_index)
                    .with_context(|| format!("missing side-reference cache for row {row_index}"))?;
                for (motif_index, &frequency) in row_values.iter().enumerate() {
                    let reference_motif_label =
                        ref_kmers.motif_labels().get(motif_index).with_context(|| {
                            format!("reference motif index {motif_index} is outside motif axis")
                        })?;
                    cache.add_frequency(reference_motif_label, frequency, shape)?;
                }
            }
        }
        RefKmerFrequencyData::Sparse(sparse) => {
            for entry in sparse.entries() {
                if !requested_rows.contains(&entry.row_index) {
                    continue;
                }
                let reference_motif_label = ref_kmers
                    .motif_labels()
                    .get(entry.motif_index)
                    .with_context(|| {
                        format!(
                            "reference sparse motif index {} is outside motif axis",
                            entry.motif_index
                        )
                    })?;
                let cache = caches.get_mut(&entry.row_index).with_context(|| {
                    format!("missing side-reference cache for row {}", entry.row_index)
                })?;
                cache.add_frequency(reference_motif_label, entry.frequency, shape)?;
            }
        }
    }

    for cache in caches.values_mut() {
        cache.finalize_support_counts();
    }
    Ok(caches)
}

/// Correct full motifs using independently calculated side denominators.
///
/// Each motif keeps its original column. Its outside and inside marginal
/// denominators are looked up from the matched reference row and multiplied
/// before the unsupported-reference policy is applied.
fn split_corrected_counts_data(
    selection: &EndMotifCountSelection,
    reference_row_indices: &[usize],
    parsed_motifs: &[ParsedEndMotif],
    reference_caches: &BTreeMap<usize, SideReferenceCache>,
    unsupported_reference_policy: UnsupportedReferencePolicy,
) -> Result<EndMotifCountsData> {
    ensure!(
        reference_row_indices.len() == selection.row_count(),
        "reference row selection does not match end-motif row count"
    );
    ensure!(
        parsed_motifs.len() == selection.motif_count(),
        "parsed motif count does not match selected motif count"
    );
    correct_counts_data(
        selection.data(),
        selection.motif_labels(),
        |selected_row_index, selected_motif_index| {
            let reference_row_index = reference_row_indices[selected_row_index];
            let reference_cache =
                reference_caches
                    .get(&reference_row_index)
                    .with_context(|| {
                        format!(
                            "missing side-reference cache for reference row {reference_row_index}"
                        )
                    })?;
            Ok(reference_cache.split_denominator(&parsed_motifs[selected_motif_index]))
        },
        unsupported_reference_policy,
    )
}

/// Aggregate full-motif counts onto the selected side axis.
///
/// The dense or sparse representation is preserved while all source motifs
/// mapped to the same selected side column are summed.
fn aggregate_side_counts_data(
    selection: &EndMotifCountSelection,
    side_axis: &SideAxisSelection,
) -> Result<EndMotifCountsData> {
    match selection.data() {
        EndMotifCountsData::Dense(counts) => Ok(EndMotifCountsData::Dense(
            aggregate_dense_side_counts(counts, side_axis)?,
        )),
        EndMotifCountsData::Sparse(sparse) => Ok(EndMotifCountsData::Sparse(
            aggregate_sparse_side_counts(sparse, side_axis)?,
        )),
    }
}

/// Sum dense full-motif counts into selected side-axis columns.
///
/// Every source coordinate contributes to its mapped side column unless that
/// side label was not selected. The function checks index arithmetic and rejects
/// non-finite aggregated counts.
fn aggregate_dense_side_counts(
    counts: &DenseMatrix<f64>,
    side_axis: &SideAxisSelection,
) -> Result<DenseMatrix<f64>> {
    let (row_count, source_motif_count) = counts.shape();
    ensure!(
        side_axis.source_motif_to_selected_index.len() == source_motif_count,
        "side-axis mapping does not match selected joint motif count"
    );
    let selected_side_count = side_axis.selected_labels.len();
    let mut aggregated_values = vec![0.0; row_count.saturating_mul(selected_side_count)];
    for selected_row_index in 0..row_count {
        for source_motif_index in 0..source_motif_count {
            let Some(selected_side_index) =
                side_axis.source_motif_to_selected_index[source_motif_index]
            else {
                continue;
            };
            let count = counts
                .get(selected_row_index, source_motif_index)
                .copied()
                .with_context(|| {
                    format!(
                        "selected end-motif count coordinate ({selected_row_index}, {source_motif_index}) is outside dense matrix"
                    )
                })?;
            let value_index = selected_row_index
                .checked_mul(selected_side_count)
                .and_then(|row_start| row_start.checked_add(selected_side_index))
                .with_context(|| "side-mode dense aggregation index overflowed")?;
            let aggregated_count = aggregated_values[value_index] + count;
            ensure!(
                aggregated_count.is_finite(),
                "side-mode aggregation produced non-finite count for motif '{}'",
                side_axis.selected_labels[selected_side_index]
            );
            aggregated_values[value_index] = aggregated_count;
        }
    }
    DenseMatrix::from_row_major(aggregated_values, row_count, selected_side_count)
}

/// Sum stored sparse counts into selected side-axis coordinates.
///
/// Coordinates mapped to the same row and side label are combined in a sorted
/// map. Zero totals are omitted, preserving sparse storage and derived shape.
fn aggregate_sparse_side_counts(
    counts: &EndMotifSparseCounts,
    side_axis: &SideAxisSelection,
) -> Result<EndMotifSparseCounts> {
    let (row_count, source_motif_count) = counts.shape();
    ensure!(
        side_axis.source_motif_to_selected_index.len() == source_motif_count,
        "side-axis mapping does not match selected joint motif count"
    );
    let selected_side_count = side_axis.selected_labels.len();
    let mut aggregated_counts = BTreeMap::<(usize, usize), f64>::new();
    for entry in counts.entries() {
        let Some(selected_side_index) = side_axis
            .source_motif_to_selected_index
            .get(entry.motif_index)
            .copied()
            .flatten()
        else {
            continue;
        };
        let aggregated_count = aggregated_counts
            .entry((entry.row_index, selected_side_index))
            .or_insert(0.0);
        *aggregated_count += entry.count;
        ensure!(
            aggregated_count.is_finite(),
            "side-mode aggregation produced non-finite count for motif '{}'",
            side_axis.selected_labels[selected_side_index]
        );
    }
    let entries = aggregated_counts
        .into_iter()
        .filter_map(|((row_index, motif_index), count)| {
            (count != 0.0).then_some(EndMotifSparseEntry {
                row_index,
                motif_index,
                count,
            })
        })
        .collect();
    Ok(EndMotifSparseCounts::from_entries(
        row_count,
        selected_side_count,
        entries,
    ))
}

/// Correct aggregated side counts with matching marginal denominators.
///
/// Each derived side label is looked up in the cache for its matched reference
/// row. The dense or sparse aggregated representation is retained.
fn side_corrected_counts_data(
    aggregated_data: &EndMotifCountsData,
    reference_row_indices: &[usize],
    side_labels: &[String],
    side_mode: SideMode,
    reference_caches: &BTreeMap<usize, SideReferenceCache>,
    unsupported_reference_policy: UnsupportedReferencePolicy,
) -> Result<EndMotifCountsData> {
    let (row_count, motif_count) = aggregated_data.shape();
    ensure!(
        reference_row_indices.len() == row_count,
        "reference row selection does not match end-motif row count"
    );
    ensure!(
        side_labels.len() == motif_count,
        "side label count does not match aggregated motif count"
    );
    correct_counts_data(
        aggregated_data,
        side_labels,
        |selected_row_index, selected_motif_index| {
            let reference_row_index = reference_row_indices[selected_row_index];
            let reference_cache =
                reference_caches
                    .get(&reference_row_index)
                    .with_context(|| {
                        format!(
                            "missing side-reference cache for reference row {reference_row_index}"
                        )
                    })?;
            Ok(reference_cache.denominator(&side_labels[selected_motif_index], side_mode))
        },
        unsupported_reference_policy,
    )
}

/// Match selected sample labels to optional reference motif indices.
///
/// Concrete labels lose their underscore before matching, while motif-group
/// labels match directly. Missing labels remain `None` and therefore have zero
/// reference support during correction.
fn selected_reference_motif_indices(
    ends: &EndsOutput,
    ref_kmers: &RefKmersOutput,
    selection: &EndMotifCountSelection,
) -> Vec<Option<usize>> {
    selection
        .motif_labels()
        .iter()
        .map(|motif_label| {
            let reference_label = reference_motif_label(ends, motif_label);
            ref_kmers.motif_index(&reference_label).ok()
        })
        .collect()
}

/// Count reference motifs used for uniform-support normalization per row.
///
/// Only strictly positive frequencies count as supported. Dense and sparse
/// reference stores are handled without densifying sparse data.
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

/// Correct exact-label counts using relative reference motif frequencies.
///
/// For each coordinate, the denominator is the matched frequency multiplied by
/// the number of positive-frequency motifs in that reference row. This makes a
/// uniform reference composition leave counts unchanged.
fn corrected_counts_data(
    selection: &EndMotifCountSelection,
    reference_row_indices: &[usize],
    reference_motif_indices: &[Option<usize>],
    reference_support_counts: &BTreeMap<usize, usize>,
    reference_frequencies: &ReferenceFrequencies<'_>,
    unsupported_reference_policy: UnsupportedReferencePolicy,
) -> Result<EndMotifCountsData> {
    ensure!(
        reference_row_indices.len() == selection.row_count(),
        "reference row selection does not match end-motif row count"
    );
    ensure!(
        reference_motif_indices.len() == selection.motif_count(),
        "reference motif selection does not match end-motif motif count"
    );
    correct_counts_data(
        selection.data(),
        selection.motif_labels(),
        |selected_row_index, selected_motif_index| {
            let reference_row_index = reference_row_indices[selected_row_index];
            let number_of_supported_motifs = reference_support_counts
                .get(&reference_row_index)
                .copied()
                .unwrap_or(0);
            let reference_frequency = reference_motif_indices[selected_motif_index]
                .and_then(|reference_motif_index| {
                    reference_frequencies.frequency(reference_row_index, reference_motif_index)
                })
                .unwrap_or(0.0);
            Ok(reference_frequency * number_of_supported_motifs as f64)
        },
        unsupported_reference_policy,
    )
}

/// Dispatch correction to the dense or sparse count representation.
///
/// Both paths use the same coordinate-specific denominator callback and
/// unsupported-reference policy, so only storage mechanics differ.
fn correct_counts_data(
    counts_data: &EndMotifCountsData,
    motif_labels: &[String],
    denominator_for_coordinate: impl Fn(usize, usize) -> Result<f64>,
    unsupported_reference_policy: UnsupportedReferencePolicy,
) -> Result<EndMotifCountsData> {
    match counts_data {
        EndMotifCountsData::Dense(counts) => Ok(EndMotifCountsData::Dense(correct_dense_counts(
            counts,
            motif_labels,
            &denominator_for_coordinate,
            unsupported_reference_policy,
        )?)),
        EndMotifCountsData::Sparse(counts) => {
            Ok(EndMotifCountsData::Sparse(correct_sparse_counts(
                counts,
                motif_labels,
                &denominator_for_coordinate,
                unsupported_reference_policy,
            )?))
        }
    }
}

/// Correct every coordinate in a dense count matrix.
///
/// The full matrix shape is preserved. Unsupported positive counts and
/// non-finite results are collected across coordinates and reported after the
/// complete pass, producing deterministic motif-level errors.
fn correct_dense_counts<F>(
    counts: &DenseMatrix<f64>,
    motif_labels: &[String],
    denominator_for_coordinate: &F,
    unsupported_reference_policy: UnsupportedReferencePolicy,
) -> Result<DenseMatrix<f64>>
where
    F: Fn(usize, usize) -> Result<f64>,
{
    let (row_count, motif_count) = counts.shape();
    ensure!(
        motif_labels.len() == motif_count,
        "motif label count does not match dense motif count"
    );

    let mut unsupported_positive_labels = BTreeSet::new();
    let mut non_finite_labels = BTreeSet::new();
    let mut corrected_values = Vec::with_capacity(row_count.saturating_mul(motif_count));
    for selected_row_index in 0..row_count {
        for (selected_motif_index, motif_label) in motif_labels.iter().enumerate() {
            let count = counts
                .get(selected_row_index, selected_motif_index)
                .copied()
                .with_context(|| {
                    format!(
                        "selected end-motif count coordinate ({selected_row_index}, {selected_motif_index}) is outside dense matrix"
                    )
                })?;
            let denominator = denominator_for_coordinate(selected_row_index, selected_motif_index)?;
            corrected_values.push(correct_with_denominator(
                count,
                denominator,
                motif_label,
                unsupported_reference_policy,
                &mut unsupported_positive_labels,
                &mut non_finite_labels,
            ));
        }
    }
    ensure_no_non_finite_corrected_counts(&non_finite_labels)?;
    ensure_no_unsupported_positive_counts(&unsupported_positive_labels)?;
    DenseMatrix::from_row_major(corrected_values, row_count, motif_count)
}

/// Correct each stored coordinate in a sparse count matrix.
///
/// Implicit zeros remain implicit. Corrected zeros are omitted, while `NaN`
/// values are stored so the keep-missing policy remains visible to callers.
fn correct_sparse_counts<F>(
    counts: &EndMotifSparseCounts,
    motif_labels: &[String],
    denominator_for_coordinate: &F,
    unsupported_reference_policy: UnsupportedReferencePolicy,
) -> Result<EndMotifSparseCounts>
where
    F: Fn(usize, usize) -> Result<f64>,
{
    let (row_count, motif_count) = counts.shape();
    ensure!(
        motif_labels.len() == motif_count,
        "motif label count does not match sparse motif count"
    );

    let mut unsupported_positive_labels = BTreeSet::new();
    let mut non_finite_labels = BTreeSet::new();
    let mut corrected_entries = Vec::new();
    for entry in counts.entries() {
        ensure!(
            entry.row_index < row_count,
            "sparse end-motif row index {} is outside selected row count",
            entry.row_index
        );
        let motif_label = motif_labels.get(entry.motif_index).with_context(|| {
            format!(
                "sparse end-motif motif index {} is outside selected motif count",
                entry.motif_index
            )
        })?;
        let denominator = denominator_for_coordinate(entry.row_index, entry.motif_index)?;
        let corrected_value = correct_with_denominator(
            entry.count,
            denominator,
            motif_label,
            unsupported_reference_policy,
            &mut unsupported_positive_labels,
            &mut non_finite_labels,
        );
        if corrected_value != 0.0 || corrected_value.is_nan() {
            corrected_entries.push(EndMotifSparseEntry {
                row_index: entry.row_index,
                motif_index: entry.motif_index,
                count: corrected_value,
            });
        }
    }
    ensure_no_non_finite_corrected_counts(&non_finite_labels)?;
    ensure_no_unsupported_positive_counts(&unsupported_positive_labels)?;
    Ok(EndMotifSparseCounts::from_entries(
        row_count,
        motif_count,
        corrected_entries,
    ))
}

/// Apply a correction denominator to one count.
///
/// Finite counts with positive denominators are divided normally. Invalid
/// arithmetic is recorded as non-finite, unsupported zero counts stay zero,
/// and unsupported positive counts follow the configured error or `NaN` policy.
fn correct_with_denominator(
    count: f64,
    denominator: f64,
    motif_label: &str,
    unsupported_reference_policy: UnsupportedReferencePolicy,
    unsupported_positive_labels: &mut BTreeSet<String>,
    non_finite_labels: &mut BTreeSet<String>,
) -> f64 {
    if !count.is_finite() || !denominator.is_finite() {
        non_finite_labels.insert(motif_label.to_string());
        return 0.0;
    }
    if denominator > 0.0 {
        let corrected_count = count / denominator;
        if corrected_count.is_finite() {
            return corrected_count;
        }
        non_finite_labels.insert(motif_label.to_string());
        return 0.0;
    }
    if count <= 0.0 {
        return 0.0;
    }
    match unsupported_reference_policy {
        UnsupportedReferencePolicy::Error => {
            unsupported_positive_labels.insert(motif_label.to_string());
            0.0
        }
        UnsupportedReferencePolicy::KeepNaN => f64::NAN,
    }
}

/// Fail when correction produced non-finite values for any motif.
///
/// Labels are collected in a sorted set during correction, so the final error
/// reports every affected motif in deterministic order.
fn ensure_no_non_finite_corrected_counts(non_finite_labels: &BTreeSet<String>) -> Result<()> {
    if non_finite_labels.is_empty() {
        return Ok(());
    }
    let labels = non_finite_labels
        .iter()
        .map(String::as_str)
        .collect::<Vec<_>>()
        .join(", ");
    bail!("reference correction produced non-finite corrected counts for motifs: {labels}");
}

/// Fail when positive counts lack a correction denominator under error policy.
///
/// All affected labels are reported together after dense or sparse correction,
/// rather than failing at the first unsupported coordinate.
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
        "positive-count end motifs have no positive correction denominator under the selected mode: {labels}. Use unsupported_reference_policy(UnsupportedReferencePolicy::KeepNaN) to keep them as NaN"
    );
}

/// Convert a sample motif label to its matching reference-axis label.
///
/// Concrete end motifs drop the outside/inside underscore. Motif-group labels
/// are already shared between outputs and therefore remain unchanged.
fn reference_motif_label(ends: &EndsOutput, end_motif_label: &str) -> String {
    match ends.motif_axis_kind() {
        EndMotifAxisKind::Motif => end_motif_label.replace('_', ""),
        EndMotifAxisKind::MotifGroup => end_motif_label.to_string(),
    }
}

/// Return correction keys for every end-motif row.
///
/// Global, windowed, and grouped metadata are converted to the same key type
/// used for exact sample-to-reference row comparison.
fn all_end_row_keys(ends: &EndsOutput) -> Result<Vec<CorrectionRowKey>> {
    match ends.row_metadata() {
        EndMotifRowMetadata::Global => Ok(vec![CorrectionRowKey::Global]),
        EndMotifRowMetadata::Windows { windows, .. } => windows
            .iter()
            .map(window_row_key)
            .collect::<Result<Vec<_>>>(),
        EndMotifRowMetadata::Groups(groups) => Ok(groups
            .iter()
            .map(|group| CorrectionRowKey::Group(group.name.clone()))
            .collect()),
    }
}

/// Return correction keys for every reference k-mer row.
///
/// Keys use the same representation as end-motif rows, allowing complete row
/// sets to be compared independently of their storage types.
fn all_reference_row_keys(ref_kmers: &RefKmersOutput) -> Result<Vec<CorrectionRowKey>> {
    match ref_kmers.row_metadata() {
        RefKmerRowMetadata::Global => Ok(vec![CorrectionRowKey::Global]),
        RefKmerRowMetadata::Windows { windows, .. } => windows
            .iter()
            .map(window_row_key)
            .collect::<Result<Vec<_>>>(),
        RefKmerRowMetadata::Groups(groups) => Ok(groups
            .iter()
            .map(|group| CorrectionRowKey::Group(group.name.clone()))
            .collect()),
    }
}

/// Resolve an end-motif row index to its correction key.
///
/// The index is checked against the row-mode-specific metadata. Window keys
/// include index, chromosome, start, and end, while grouped keys use the name.
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
            window_row_key(window)
        }
        EndMotifRowMetadata::Groups(groups) => {
            let group = groups.get(row_index).with_context(|| {
                format!("end-motif row index {row_index} is outside group metadata")
            })?;
            Ok(CorrectionRowKey::Group(group.name.clone()))
        }
    }
}

/// Map each unique reference correction key to its row index.
///
/// Duplicate keys are rejected because they would make sample-to-reference row
/// matching ambiguous.
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

/// Resolve a reference row index to its correction key.
///
/// The index is checked against the reference row metadata and converted to the
/// same global, window, or group key used by end-motif rows.
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
            window_row_key(window)
        }
        RefKmerRowMetadata::Groups(groups) => {
            let group = groups.get(row_index).with_context(|| {
                format!("reference row index {row_index} is outside group metadata")
            })?;
            Ok(CorrectionRowKey::Group(group.name.clone()))
        }
    }
}

/// Convert window metadata into a shared correction row key.
///
/// The checked interval is represented by its `start` and `end` coordinates
/// together with the stored window index and chromosome.
fn window_row_key(window: &WindowRow) -> Result<CorrectionRowKey> {
    let (start, end) = window.interval.as_tuple();
    Ok(CorrectionRowKey::Window {
        index: window.index,
        chrom: window.chrom.clone(),
        start,
        end,
    })
}

/// Return a human-readable end-motif row-mode name for errors.
///
/// Keeping this mapping centralized makes row-mode mismatch messages precise
/// without exposing enum debug names to users.
fn describe_end_row_mode(row_mode: EndMotifRowMode) -> &'static str {
    match row_mode {
        EndMotifRowMode::Global => "global",
        EndMotifRowMode::SizeWindows => "size windows",
        EndMotifRowMode::BedWindows => "BED windows",
        EndMotifRowMode::Groups => "groups",
    }
}

/// Return a human-readable reference row-mode name for errors.
///
/// The wording mirrors end-motif mode descriptions so mismatch errors compare
/// equivalent concepts.
fn describe_ref_row_mode(row_mode: RefKmerRowMode) -> &'static str {
    match row_mode {
        RefKmerRowMode::Global => "global",
        RefKmerRowMode::SizeWindows => "size windows",
        RefKmerRowMode::BedWindows => "BED windows",
        RefKmerRowMode::Groups => "groups",
    }
}

#[cfg(test)]
mod tests {
    include!("reference_correction_tests.rs");
}
