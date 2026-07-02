//! Loader for `cfdna ref-kmers` Zarr outputs.
//!
//! Reference k-mer outputs store row-wise frequencies plus one row scaling factor. Downstream code
//! can use the frequencies directly or reconstruct counts as:
//!
//! ```text
//! count = frequency * row_scaling_factor[row]
//! ```
//!
//! The loader reads and validates the store metadata eagerly. Dense frequency stores are read into
//! a `DenseMatrix<f64>`. Sparse stores are read as sorted COO vectors. Counts are exposed through
//! point lookups and dense reconstruction without changing the on-disk frequency contract.

use crate::{
    interval::Interval,
    output_loaders::{
        OutputLoaderError, OutputLoaderResult,
        common::{
            DenseMatrix, WindowRow, build_selection_index_map, contiguous_index_span,
            ensure_unique_indices, ensure_unique_labels, resolve_row_indices,
            validate_zarr_public_label,
        },
    },
    shared::reference::{ContigFootprintEntry, twobit_contig_footprint},
    shared::zarr::read_zarr_root_attributes,
};
use anyhow::{Context, Result, bail, ensure};
use fxhash::FxHashMap;
use serde_json::Value;
use std::{
    fs,
    path::{Path, PathBuf},
    sync::Arc,
};
use zarrs::{
    array::{Array, ElementOwned},
    filesystem::FilesystemStore,
};

const REF_KMER_SCHEMA_VERSION: u64 = 1;

/// Load a `cfdna ref-kmers` Zarr store.
///
/// The path must point to a Zarr directory written with the supported reference k-mer schema
/// version. The loader reads row metadata, motif labels, row scaling factors, and either dense
/// frequencies or sparse COO frequencies into owned Rust containers.
///
/// Parameters
/// ----------
/// - `path`:
///   Path to a `cfdna ref-kmers` Zarr output directory.
///
/// Returns
/// -------
/// - `RefKmersOutput`:
///   Loaded row metadata, motif labels, frequencies, and scaling factors.
pub fn load_ref_kmers_output(path: impl AsRef<Path>) -> OutputLoaderResult<RefKmersOutput> {
    RefKmersParser::new(path.as_ref())
        .load()
        .map_err(Into::into)
}

/// Loaded reference k-mer frequencies from `cfdna ref-kmers`.
#[derive(Debug, Clone, PartialEq)]
pub struct RefKmersOutput {
    row_metadata: RefKmerRowMetadata,
    motif_axis_kind: RefKmerMotifAxisKind,
    motif_labels: Vec<String>,
    motif_label_indices: FxHashMap<String, usize>,
    row_scaling_factors: Vec<f64>,
    data: RefKmerFrequencyData,
    kmer_size: u8,
    canonical: bool,
    all_motifs: bool,
    assign_by: String,
    reference_contig_footprint: Vec<ContigFootprintEntry>,
}

impl RefKmersOutput {
    /// Return how frequencies were stored on disk.
    pub fn storage_mode(&self) -> RefKmerStorageMode {
        self.data.storage_mode()
    }

    /// Return what each frequency row represents.
    pub fn row_mode(&self) -> RefKmerRowMode {
        self.row_metadata.mode()
    }

    /// Return whether motif-axis labels are concrete k-mers or motif groups.
    pub fn motif_axis_kind(&self) -> RefKmerMotifAxisKind {
        self.motif_axis_kind
    }

    /// Return the configured k-mer size.
    pub fn kmer_size(&self) -> u8 {
        self.kmer_size
    }

    /// Return whether the output collapsed reverse-complement-equivalent k-mers.
    pub fn canonical(&self) -> bool {
        self.canonical
    }

    /// Return whether the command requested all possible motifs on the output axis.
    pub fn all_motifs(&self) -> bool {
        self.all_motifs
    }

    /// Return the window-assignment mode recorded in the Zarr metadata.
    pub fn assign_by(&self) -> &str {
        &self.assign_by
    }

    /// Return the reference contig footprint stored in the Zarr package.
    ///
    /// This is the sorted `(contig name, contig length)` identity of the `.2bit` reference used
    /// when `ref-kmers` produced the package. It intentionally does not include file paths or
    /// sequence content.
    pub fn reference_contig_footprint(&self) -> &[ContigFootprintEntry] {
        &self.reference_contig_footprint
    }

    /// Return an error when a `.2bit` reference does not match this package.
    ///
    /// This computes the `.2bit` contig footprint and compares it to the footprint stored by
    /// `ref-kmers`, so downstream tools can fail before using a precomputed package with the wrong
    /// reference.
    pub fn ensure_reference_2bit_matches(
        &self,
        reference_2bit: impl AsRef<Path>,
    ) -> OutputLoaderResult<()> {
        let reference_2bit = reference_2bit.as_ref();
        let reference_contig_footprint =
            twobit_contig_footprint(reference_2bit).with_context(|| {
                format!(
                    "read reference contig footprint from {}",
                    reference_2bit.display()
                )
            })?;
        if self.reference_contig_footprint == reference_contig_footprint {
            return Ok(());
        }
        Err(OutputLoaderError::message(
            "reference k-mer package was built against a different reference contig footprint",
        ))
    }

    /// Return motif-axis labels in frequency-column order.
    pub fn motif_labels(&self) -> &[String] {
        &self.motif_labels
    }

    /// Return row scaling factors in frequency-row order.
    pub fn row_scaling_factors(&self) -> &[f64] {
        &self.row_scaling_factors
    }

    /// Return one row scaling factor by zero-based row index.
    pub fn row_scaling_factor(&self, row_index: usize) -> Option<f64> {
        self.row_scaling_factors.get(row_index).copied()
    }

    /// Return row metadata describing the frequency rows.
    pub fn row_metadata(&self) -> &RefKmerRowMetadata {
        &self.row_metadata
    }

    /// Return window metadata, or an error if this is not a windowed output.
    pub fn window_metadata(&self) -> OutputLoaderResult<&[WindowRow]> {
        match &self.row_metadata {
            RefKmerRowMetadata::Windows { windows, .. } => Ok(windows),
            _ => Err(OutputLoaderError::message(
                "reference k-mer output is not windowed",
            )),
        }
    }

    /// Return group metadata, or an error if this is not a grouped output.
    pub fn group_metadata(&self) -> OutputLoaderResult<&[RefKmerGroupRow]> {
        match &self.row_metadata {
            RefKmerRowMetadata::Groups(groups) => Ok(groups),
            _ => Err(OutputLoaderError::message(
                "reference k-mer output is not grouped",
            )),
        }
    }

    /// Return one window row by zero-based row index.
    pub fn window(&self, row_index: usize) -> OutputLoaderResult<Option<&WindowRow>> {
        Ok(self.window_metadata()?.get(row_index))
    }

    /// Return one group row by zero-based row index.
    pub fn group(&self, row_index: usize) -> OutputLoaderResult<Option<&RefKmerGroupRow>> {
        Ok(self.group_metadata()?.get(row_index))
    }

    /// Return the group row index for one group name.
    pub fn group_index(&self, group_name: &str) -> OutputLoaderResult<usize> {
        let groups = self.group_metadata()?;
        Ok(groups
            .iter()
            .find(|group| group.name == group_name)
            .map(|group| group.index)
            .with_context(|| format!("reference k-mer output has no group named '{group_name}'"))?)
    }

    /// Return whether one group name exists in a grouped output.
    ///
    /// This returns `false` for non-grouped outputs.
    pub fn has_group(&self, group_name: &str) -> bool {
        self.group_metadata()
            .is_ok_and(|groups| groups.iter().any(|group| group.name == group_name))
    }

    /// Return the motif index for one label.
    pub fn motif_index(&self, motif_label: &str) -> OutputLoaderResult<usize> {
        Ok(self
            .motif_label_indices
            .get(motif_label)
            .copied()
            .with_context(|| {
                format!("reference k-mer output has no motif label '{motif_label}'")
            })?)
    }

    /// Return whether one motif label exists.
    pub fn has_motif(&self, motif_label: &str) -> bool {
        self.motif_label_indices.contains_key(motif_label)
    }

    /// Return a compact description of the loaded reference k-mer output.
    pub fn output_metadata(&self) -> RefKmerOutputMetadata {
        RefKmerOutputMetadata {
            storage_mode: self.storage_mode(),
            row_mode: self.row_mode(),
            motif_axis_kind: self.motif_axis_kind(),
            row_count: self.row_count(),
            motif_count: self.motif_count(),
            kmer_size: self.kmer_size,
            canonical: self.canonical,
            all_motifs: self.all_motifs,
            assign_by: self.assign_by.clone(),
            reference_contig_footprint: self.reference_contig_footprint.clone(),
        }
    }

    /// Return the number of frequency rows.
    pub fn row_count(&self) -> usize {
        self.row_scaling_factors.len()
    }

    /// Return the number of motif columns.
    pub fn motif_count(&self) -> usize {
        self.motif_labels.len()
    }

    /// Return one frequency value, if both indices are in bounds.
    ///
    /// Sparse outputs return `0.0` for in-bounds entries that are not stored in the COO vectors.
    pub fn frequency(&self, row_index: usize, motif_index: usize) -> Option<f64> {
        match &self.data {
            RefKmerFrequencyData::Dense(frequencies) => {
                frequencies.get(row_index, motif_index).copied()
            }
            RefKmerFrequencyData::Sparse(sparse) => sparse.frequency(row_index, motif_index),
        }
    }

    /// Return one frequency value by motif label, if the row index is in bounds.
    pub fn frequency_for_motif(
        &self,
        row_index: usize,
        motif_label: &str,
    ) -> OutputLoaderResult<Option<f64>> {
        Ok(self.frequency(row_index, self.motif_index(motif_label)?))
    }

    /// Return one reconstructed count value, if both indices are in bounds.
    ///
    /// Ref-kmers stores frequencies on disk. Counts are reconstructed as
    /// `frequency(row, motif) * row_scaling_factor[row]`.
    pub fn count(&self, row_index: usize, motif_index: usize) -> Option<f64> {
        let scaling_factor = self.row_scaling_factor(row_index)?;
        self.frequency(row_index, motif_index)
            .map(|frequency| frequency * scaling_factor)
    }

    /// Return one reconstructed count by motif label, if the row index is in bounds.
    pub fn count_for_motif(
        &self,
        row_index: usize,
        motif_label: &str,
    ) -> OutputLoaderResult<Option<f64>> {
        Ok(self.count(row_index, self.motif_index(motif_label)?))
    }

    /// Return dense frequencies, or an error if this store is sparse.
    pub fn dense_frequencies(&self) -> OutputLoaderResult<&DenseMatrix<f64>> {
        match &self.data {
            RefKmerFrequencyData::Dense(frequencies) => Ok(frequencies),
            RefKmerFrequencyData::Sparse(_) => Err(OutputLoaderError::message(
                "reference k-mer output is sparse",
            )),
        }
    }

    /// Return sparse COO frequencies, or an error if this store is dense.
    pub fn sparse_frequencies(&self) -> OutputLoaderResult<&RefKmerSparseFrequencies> {
        match &self.data {
            RefKmerFrequencyData::Sparse(sparse) => Ok(sparse),
            RefKmerFrequencyData::Dense(_) => Err(OutputLoaderError::message(
                "reference k-mer output is dense",
            )),
        }
    }

    /// Reconstruct sparse count entries, or an error if this store is dense.
    ///
    /// Ref-kmers stores sparse frequencies on disk. This returns the stored sparse coordinates
    /// with each frequency multiplied by its row scaling factor.
    pub fn sparse_count_entries(&self) -> OutputLoaderResult<Vec<RefKmerSparseCountEntry>> {
        match &self.data {
            RefKmerFrequencyData::Sparse(sparse) => {
                Ok(sparse.count_entries(&self.row_scaling_factors)?)
            }
            RefKmerFrequencyData::Dense(_) => Err(OutputLoaderError::message(
                "reference k-mer output is dense",
            )),
        }
    }

    /// Return frequency storage as either a dense matrix or sparse COO entries.
    pub fn data(&self) -> &RefKmerFrequencyData {
        &self.data
    }

    /// Return frequencies as a dense row-major matrix.
    ///
    /// Dense outputs are cloned into the returned matrix. Sparse outputs are explicitly densified.
    pub fn to_dense_frequency_matrix(&self) -> OutputLoaderResult<DenseMatrix<f64>> {
        match &self.data {
            RefKmerFrequencyData::Dense(frequencies) => Ok(frequencies.clone()),
            RefKmerFrequencyData::Sparse(sparse) => sparse.to_dense_matrix(),
        }
    }

    /// Reconstruct counts as a dense row-major matrix.
    ///
    /// Each row is multiplied by its row scaling factor. This is the loader-side extraction path
    /// for downstream code that wants counts instead of frequencies.
    pub fn to_dense_count_matrix(&self) -> OutputLoaderResult<DenseMatrix<f64>> {
        self.data.to_dense_count_matrix(
            &self.row_scaling_factors,
            "reference k-mer dense count row offset overflow",
        )
    }

    /// Start a frequency selection.
    ///
    /// A new selector initially selects all rows and all motifs. Add row and motif constraints
    /// before calling `read()`.
    pub fn select(&self) -> RefKmersSelector<'_> {
        RefKmersSelector::new(self)
    }

    fn select_frequencies(
        &self,
        row_indices: Option<&[usize]>,
        motif_indices: Option<&[usize]>,
    ) -> Result<RefKmerFrequencySelection> {
        self.select_frequencies_with_label(row_indices, motif_indices, "row")
    }

    fn select_window_frequencies(
        &self,
        window_indices: Option<&[usize]>,
        motif_indices: Option<&[usize]>,
    ) -> Result<RefKmerFrequencySelection> {
        ensure!(
            matches!(
                self.row_mode(),
                RefKmerRowMode::SizeWindows | RefKmerRowMode::BedWindows
            ),
            "reference k-mer output is not windowed"
        );
        self.select_frequencies_with_label(window_indices, motif_indices, "window")
    }

    fn select_group_frequencies(
        &self,
        group_indices: Option<&[usize]>,
        motif_indices: Option<&[usize]>,
    ) -> Result<RefKmerFrequencySelection> {
        ensure!(
            self.row_mode() == RefKmerRowMode::Groups,
            "reference k-mer output is not grouped"
        );
        self.select_frequencies_with_label(group_indices, motif_indices, "group")
    }

    fn select_group_frequencies_by_name<S: AsRef<str>>(
        &self,
        group_names: Option<&[S]>,
        motif_indices: Option<&[usize]>,
    ) -> Result<RefKmerFrequencySelection> {
        let group_indices = match group_names {
            Some(group_names) => Some(resolve_ref_kmer_group_name_indices(self, group_names)?),
            None => {
                ensure!(
                    self.row_mode() == RefKmerRowMode::Groups,
                    "reference k-mer output is not grouped"
                );
                None
            }
        };
        self.select_frequencies_with_label(group_indices.as_deref(), motif_indices, "group")
    }

    fn select_frequencies_by_motif_label<S: AsRef<str>>(
        &self,
        row_indices: Option<&[usize]>,
        motif_labels: Option<&[S]>,
    ) -> Result<RefKmerFrequencySelection> {
        let motif_indices = match motif_labels {
            Some(motif_labels) => Some(resolve_ref_kmer_motif_label_indices(self, motif_labels)?),
            None => None,
        };
        self.select_frequencies_with_label(row_indices, motif_indices.as_deref(), "row")
    }

    fn select_frequencies_with_label(
        &self,
        row_indices: Option<&[usize]>,
        motif_indices: Option<&[usize]>,
        row_label: &str,
    ) -> Result<RefKmerFrequencySelection> {
        ensure!(
            row_indices.is_none() || self.row_mode() != RefKmerRowMode::Global,
            "global reference k-mer output has no selectable row axis"
        );
        let row_indices = resolve_row_indices(row_indices, self.row_count(), row_label)?;
        let motif_indices = resolve_ref_kmer_motif_indices(motif_indices, self.motif_count());
        ensure_unique_indices(&row_indices, row_label)?;
        ensure_unique_indices(&motif_indices, "motif")?;

        let row_metadata = self.selected_row_metadata(&row_indices, row_label)?;
        let row_scaling_factors = row_indices
            .iter()
            .map(|&row_index| {
                self.row_scaling_factors
                    .get(row_index)
                    .copied()
                    .with_context(|| {
                        format!(
                            "{row_label} index {row_index} is outside 0..{}",
                            self.row_count()
                        )
                    })
            })
            .collect::<Result<Vec<_>>>()?;
        let motif_labels = motif_indices
            .iter()
            .map(|&motif_index| {
                self.motif_labels
                    .get(motif_index)
                    .cloned()
                    .with_context(|| {
                        format!(
                            "motif index {motif_index} is outside 0..{}",
                            self.motif_count()
                        )
                    })
            })
            .collect::<Result<Vec<_>>>()?;

        let data = match &self.data {
            RefKmerFrequencyData::Dense(frequencies) => {
                RefKmerFrequencyData::Dense(select_dense_ref_kmer_frequencies(
                    frequencies,
                    &row_indices,
                    &motif_indices,
                    row_label,
                )?)
            }
            RefKmerFrequencyData::Sparse(sparse) => RefKmerFrequencyData::Sparse(
                sparse.select_frequencies(&row_indices, &motif_indices, row_label)?,
            ),
        };
        Ok(RefKmerFrequencySelection {
            row_metadata,
            motif_axis_kind: self.motif_axis_kind,
            row_indices,
            motif_indices,
            motif_labels,
            row_scaling_factors,
            data,
            kmer_size: self.kmer_size,
            canonical: self.canonical,
            source_all_motifs: self.all_motifs,
            assign_by: self.assign_by.clone(),
        })
    }

    fn selected_row_metadata(
        &self,
        row_indices: &[usize],
        row_label: &str,
    ) -> Result<RefKmerRowMetadata> {
        match &self.row_metadata {
            RefKmerRowMetadata::Global => Ok(RefKmerRowMetadata::Global),
            RefKmerRowMetadata::Windows {
                window_mode,
                windows,
            } => {
                let selected_windows = row_indices
                    .iter()
                    .map(|&row_index| {
                        windows.get(row_index).cloned().with_context(|| {
                            format!(
                                "{row_label} index {row_index} is outside 0..{}",
                                windows.len()
                            )
                        })
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok(RefKmerRowMetadata::Windows {
                    window_mode: *window_mode,
                    windows: selected_windows,
                })
            }
            RefKmerRowMetadata::Groups(groups) => {
                let selected_groups = row_indices
                    .iter()
                    .map(|&row_index| {
                        groups.get(row_index).cloned().with_context(|| {
                            format!(
                                "{row_label} index {row_index} is outside 0..{}",
                                groups.len()
                            )
                        })
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok(RefKmerRowMetadata::Groups(selected_groups))
            }
        }
    }
}

/// Builder for selecting rows and motifs from a `RefKmersOutput`.
///
/// The builder starts with all rows and all motifs selected. Set at most one selector per axis.
/// For example, use `motifs()` or `motifs_by_label()`, not both. Conflicting selector calls are
/// reported by `read()`.
#[derive(Debug, Clone)]
pub struct RefKmersSelector<'a> {
    output: &'a RefKmersOutput,
    rows: RefKmerRowSelector,
    motifs: RefKmerMotifSelector,
    selection_error: Option<String>,
}

impl<'a> RefKmersSelector<'a> {
    fn new(output: &'a RefKmersOutput) -> Self {
        Self {
            output,
            rows: RefKmerRowSelector::All,
            motifs: RefKmerMotifSelector::All,
            selection_error: None,
        }
    }

    /// Select generic output rows by zero-based row index.
    pub fn rows(self, row_indices: &[usize]) -> Self {
        self.set_rows(RefKmerRowSelector::Rows(row_indices.to_vec()), "rows")
    }

    /// Select window rows by zero-based window row index.
    ///
    /// `read()` returns an error if the loaded output is not windowed.
    pub fn windows(self, window_indices: &[usize]) -> Self {
        self.set_rows(
            RefKmerRowSelector::Windows(window_indices.to_vec()),
            "windows",
        )
    }

    /// Select grouped rows by zero-based group row index.
    ///
    /// `read()` returns an error if the loaded output is not grouped.
    pub fn groups(self, group_indices: &[usize]) -> Self {
        self.set_rows(RefKmerRowSelector::Groups(group_indices.to_vec()), "groups")
    }

    /// Select grouped rows by group name.
    ///
    /// `read()` returns an error if the loaded output is not grouped or any requested name is
    /// missing or duplicated.
    pub fn groups_by_name<S: AsRef<str>>(self, group_names: &[S]) -> Self {
        self.set_rows(
            RefKmerRowSelector::GroupNames(
                group_names
                    .iter()
                    .map(|group_name| group_name.as_ref().to_string())
                    .collect(),
            ),
            "groups_by_name",
        )
    }

    /// Select motifs by zero-based motif index.
    pub fn motifs(self, motif_indices: &[usize]) -> Self {
        self.set_motifs(
            RefKmerMotifSelector::Indices(motif_indices.to_vec()),
            "motifs",
        )
    }

    /// Select motifs by motif or motif-group label.
    pub fn motifs_by_label<S: AsRef<str>>(self, motif_labels: &[S]) -> Self {
        self.set_motifs(
            RefKmerMotifSelector::Labels(
                motif_labels
                    .iter()
                    .map(|motif_label| motif_label.as_ref().to_string())
                    .collect(),
            ),
            "motifs_by_label",
        )
    }

    fn set_rows(mut self, selector: RefKmerRowSelector, selector_name: &'static str) -> Self {
        if let Some(previous_selector_name) = self.rows.selector_name() {
            self.record_axis_conflict("row", previous_selector_name, selector_name);
        } else {
            self.rows = selector;
        }
        self
    }

    fn set_motifs(mut self, selector: RefKmerMotifSelector, selector_name: &'static str) -> Self {
        if let Some(previous_selector_name) = self.motifs.selector_name() {
            self.record_axis_conflict("motif", previous_selector_name, selector_name);
        } else {
            self.motifs = selector;
        }
        self
    }

    fn record_axis_conflict(
        &mut self,
        axis_name: &'static str,
        previous_selector_name: &'static str,
        selector_name: &'static str,
    ) {
        if self.selection_error.is_none() {
            self.selection_error = Some(format!(
                "cannot combine {previous_selector_name}() and {selector_name}() on the {axis_name} axis"
            ));
        }
    }

    fn ensure_no_selector_conflict(&self) -> Result<()> {
        if let Some(selection_error) = &self.selection_error {
            bail!("{selection_error}");
        }
        Ok(())
    }

    /// Read selected frequencies and scaling factors while preserving the loaded storage mode.
    pub fn read(self) -> OutputLoaderResult<RefKmerFrequencySelection> {
        self.ensure_no_selector_conflict()?;
        let (motif_indices, motif_labels) = match self.motifs {
            RefKmerMotifSelector::All => (None, None),
            RefKmerMotifSelector::Indices(indices) => (Some(indices), None),
            RefKmerMotifSelector::Labels(labels) => (None, Some(labels)),
        };
        let motif_indices = motif_indices.as_deref();
        let motif_labels = motif_labels.as_deref();

        let selection = match self.rows {
            RefKmerRowSelector::All => {
                if let Some(motif_labels) = motif_labels {
                    self.output
                        .select_frequencies_by_motif_label(None, Some(motif_labels))
                } else {
                    self.output.select_frequencies(None, motif_indices)
                }
            }
            RefKmerRowSelector::Rows(indices) => {
                if let Some(motif_labels) = motif_labels {
                    self.output.select_frequencies_by_motif_label(
                        Some(indices.as_slice()),
                        Some(motif_labels),
                    )
                } else {
                    self.output
                        .select_frequencies(Some(indices.as_slice()), motif_indices)
                }
            }
            RefKmerRowSelector::Windows(indices) => {
                let resolved_motif_indices =
                    resolve_ref_kmer_motif_labels(self.output, motif_labels)?;
                let motif_indices = motif_indices.or(resolved_motif_indices.as_deref());
                self.output
                    .select_window_frequencies(Some(indices.as_slice()), motif_indices)
            }
            RefKmerRowSelector::Groups(indices) => {
                let resolved_motif_indices =
                    resolve_ref_kmer_motif_labels(self.output, motif_labels)?;
                let motif_indices = motif_indices.or(resolved_motif_indices.as_deref());
                self.output
                    .select_group_frequencies(Some(indices.as_slice()), motif_indices)
            }
            RefKmerRowSelector::GroupNames(names) => {
                let resolved_motif_indices =
                    resolve_ref_kmer_motif_labels(self.output, motif_labels)?;
                let motif_indices = motif_indices.or(resolved_motif_indices.as_deref());
                self.output
                    .select_group_frequencies_by_name(Some(names.as_slice()), motif_indices)
            }
        }?;
        Ok(selection)
    }
}

#[derive(Debug, Clone)]
enum RefKmerRowSelector {
    All,
    Rows(Vec<usize>),
    Windows(Vec<usize>),
    Groups(Vec<usize>),
    GroupNames(Vec<String>),
}

impl RefKmerRowSelector {
    fn selector_name(&self) -> Option<&'static str> {
        match self {
            Self::All => None,
            Self::Rows(_) => Some("rows"),
            Self::Windows(_) => Some("windows"),
            Self::Groups(_) => Some("groups"),
            Self::GroupNames(_) => Some("groups_by_name"),
        }
    }
}

#[derive(Debug, Clone)]
enum RefKmerMotifSelector {
    All,
    Indices(Vec<usize>),
    Labels(Vec<String>),
}

impl RefKmerMotifSelector {
    fn selector_name(&self) -> Option<&'static str> {
        match self {
            Self::All => None,
            Self::Indices(_) => Some("motifs"),
            Self::Labels(_) => Some("motifs_by_label"),
        }
    }
}

/// Selected reference k-mer frequencies with row and motif-axis metadata.
///
/// Selections preserve the storage mode of the loaded output. Use frequency methods when
/// downstream code wants the stored values. Use count methods when it wants reconstructed counts.
#[derive(Debug, Clone, PartialEq)]
pub struct RefKmerFrequencySelection {
    row_metadata: RefKmerRowMetadata,
    motif_axis_kind: RefKmerMotifAxisKind,
    row_indices: Vec<usize>,
    motif_indices: Vec<usize>,
    motif_labels: Vec<String>,
    row_scaling_factors: Vec<f64>,
    data: RefKmerFrequencyData,
    kmer_size: u8,
    canonical: bool,
    source_all_motifs: bool,
    assign_by: String,
}

impl RefKmerFrequencySelection {
    /// Return how selected frequencies are stored.
    pub fn storage_mode(&self) -> RefKmerStorageMode {
        self.data.storage_mode()
    }

    /// Return what each selected frequency row represents.
    pub fn row_mode(&self) -> RefKmerRowMode {
        self.row_metadata.mode()
    }

    /// Return whether selected motif-axis labels are concrete k-mers or motif groups.
    pub fn motif_axis_kind(&self) -> RefKmerMotifAxisKind {
        self.motif_axis_kind
    }

    /// Return the configured k-mer size.
    pub fn kmer_size(&self) -> u8 {
        self.kmer_size
    }

    /// Return whether the source output collapsed reverse-complement-equivalent k-mers.
    pub fn canonical(&self) -> bool {
        self.canonical
    }

    /// Return whether the source output included every possible motif before this selection.
    pub fn source_all_motifs(&self) -> bool {
        self.source_all_motifs
    }

    /// Return the window-assignment mode recorded in the Zarr metadata.
    pub fn assign_by(&self) -> &str {
        &self.assign_by
    }

    /// Return selected row metadata in selection order.
    pub fn row_metadata(&self) -> &RefKmerRowMetadata {
        &self.row_metadata
    }

    /// Return selected window metadata, or an error if the selection is not windowed.
    pub fn window_metadata(&self) -> OutputLoaderResult<&[WindowRow]> {
        match &self.row_metadata {
            RefKmerRowMetadata::Windows { windows, .. } => Ok(windows),
            _ => Err(OutputLoaderError::message(
                "reference k-mer frequency selection is not windowed",
            )),
        }
    }

    /// Return selected group metadata, or an error if the selection is not grouped.
    pub fn group_metadata(&self) -> OutputLoaderResult<&[RefKmerGroupRow]> {
        match &self.row_metadata {
            RefKmerRowMetadata::Groups(groups) => Ok(groups),
            _ => Err(OutputLoaderError::message(
                "reference k-mer frequency selection is not grouped",
            )),
        }
    }

    /// Return selected source row indices in selection order.
    pub fn row_indices(&self) -> &[usize] {
        &self.row_indices
    }

    /// Return selected source motif indices in selection order.
    pub fn motif_indices(&self) -> &[usize] {
        &self.motif_indices
    }

    /// Return selected motif labels in selection order.
    pub fn motif_labels(&self) -> &[String] {
        &self.motif_labels
    }

    /// Return selected row scaling factors in selection order.
    pub fn row_scaling_factors(&self) -> &[f64] {
        &self.row_scaling_factors
    }

    /// Return one selected row scaling factor by zero-based selected row index.
    pub fn row_scaling_factor(&self, row_index: usize) -> Option<f64> {
        self.row_scaling_factors.get(row_index).copied()
    }

    /// Return the selected matrix shape as `(rows, motifs)`.
    pub fn shape(&self) -> (usize, usize) {
        match &self.data {
            RefKmerFrequencyData::Dense(frequencies) => frequencies.shape(),
            RefKmerFrequencyData::Sparse(sparse) => sparse.shape(),
        }
    }

    /// Return the number of selected rows.
    pub fn row_count(&self) -> usize {
        self.shape().0
    }

    /// Return the number of selected motifs.
    pub fn motif_count(&self) -> usize {
        self.shape().1
    }

    /// Return one selected frequency value, if both selection indices are in bounds.
    pub fn frequency(&self, row_index: usize, motif_index: usize) -> Option<f64> {
        match &self.data {
            RefKmerFrequencyData::Dense(frequencies) => {
                frequencies.get(row_index, motif_index).copied()
            }
            RefKmerFrequencyData::Sparse(sparse) => sparse.frequency(row_index, motif_index),
        }
    }

    /// Return one reconstructed selected count value, if both selection indices are in bounds.
    pub fn count(&self, row_index: usize, motif_index: usize) -> Option<f64> {
        let scaling_factor = self.row_scaling_factor(row_index)?;
        self.frequency(row_index, motif_index)
            .map(|frequency| frequency * scaling_factor)
    }

    /// Return dense selected frequencies, or an error if this selection is sparse.
    pub fn dense_frequencies(&self) -> OutputLoaderResult<&DenseMatrix<f64>> {
        match &self.data {
            RefKmerFrequencyData::Dense(frequencies) => Ok(frequencies),
            RefKmerFrequencyData::Sparse(_) => Err(OutputLoaderError::message(
                "reference k-mer frequency selection is sparse",
            )),
        }
    }

    /// Return sparse selected frequencies, or an error if this selection is dense.
    pub fn sparse_frequencies(&self) -> OutputLoaderResult<&RefKmerSparseFrequencies> {
        match &self.data {
            RefKmerFrequencyData::Sparse(sparse) => Ok(sparse),
            RefKmerFrequencyData::Dense(_) => Err(OutputLoaderError::message(
                "reference k-mer frequency selection is dense",
            )),
        }
    }

    /// Reconstruct sparse selected count entries, or an error if this selection is dense.
    pub fn sparse_count_entries(&self) -> OutputLoaderResult<Vec<RefKmerSparseCountEntry>> {
        match &self.data {
            RefKmerFrequencyData::Sparse(sparse) => {
                Ok(sparse.count_entries(&self.row_scaling_factors)?)
            }
            RefKmerFrequencyData::Dense(_) => Err(OutputLoaderError::message(
                "reference k-mer frequency selection is dense",
            )),
        }
    }

    /// Return selected frequency storage as either a dense matrix or sparse COO entries.
    pub fn data(&self) -> &RefKmerFrequencyData {
        &self.data
    }

    /// Return frequencies as a dense row-major matrix.
    pub fn to_dense_frequency_matrix(&self) -> OutputLoaderResult<DenseMatrix<f64>> {
        match &self.data {
            RefKmerFrequencyData::Dense(frequencies) => Ok(frequencies.clone()),
            RefKmerFrequencyData::Sparse(sparse) => sparse.to_dense_matrix(),
        }
    }

    /// Reconstruct counts as a dense row-major matrix.
    pub fn to_dense_count_matrix(&self) -> OutputLoaderResult<DenseMatrix<f64>> {
        self.data.to_dense_count_matrix(
            &self.row_scaling_factors,
            "reference k-mer selected dense count row offset overflow",
        )
    }
}

/// Compact metadata for loaded reference k-mer frequencies.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RefKmerOutputMetadata {
    /// Whether frequencies are dense or sparse COO.
    pub storage_mode: RefKmerStorageMode,
    /// Whether rows are global, windows, or grouped-BED groups.
    pub row_mode: RefKmerRowMode,
    /// Whether motif-axis labels are k-mers or motif groups.
    pub motif_axis_kind: RefKmerMotifAxisKind,
    /// Number of frequency rows.
    pub row_count: usize,
    /// Number of motif-axis labels.
    pub motif_count: usize,
    /// Configured k-mer size.
    pub kmer_size: u8,
    /// Whether reverse-complement-equivalent k-mers were collapsed.
    pub canonical: bool,
    /// Whether zero-frequency motifs were included.
    pub all_motifs: bool,
    /// Window-assignment mode recorded in the Zarr metadata.
    pub assign_by: String,
    /// Reference contig footprint stored in the package.
    pub reference_contig_footprint: Vec<ContigFootprintEntry>,
}

/// Reference k-mer frequency storage mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RefKmerStorageMode {
    /// Full dense `frequencies[row, motif]` array.
    Dense,
    /// Sparse coordinate arrays under `sparse/`.
    SparseCoo,
}

/// Meaning of reference k-mer rows.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RefKmerRowMode {
    /// One row covering the full selected reference.
    Global,
    /// Generated fixed-size genomic windows.
    SizeWindows,
    /// User-provided BED windows.
    BedWindows,
    /// Grouped-BED groups.
    Groups,
}

/// Meaning of the motif axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RefKmerMotifAxisKind {
    /// Concrete reference k-mer labels.
    Motif,
    /// Motif-group labels from `--motifs-file`.
    MotifGroup,
}

/// Source of window rows in a reference k-mer output.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RefKmerWindowMode {
    /// Generated fixed-size windows.
    Size,
    /// User-provided BED windows.
    Bed,
}

/// Row metadata for loaded reference k-mer frequencies.
#[derive(Debug, Clone, PartialEq)]
pub enum RefKmerRowMetadata {
    /// One global row.
    Global,
    /// One row per genomic window.
    Windows {
        /// Whether windows came from fixed-size generation or BED input.
        window_mode: RefKmerWindowMode,
        /// Window metadata in frequency-row order.
        windows: Vec<WindowRow>,
    },
    /// One row per grouped-BED group.
    Groups(Vec<RefKmerGroupRow>),
}

impl RefKmerRowMetadata {
    /// Return the row mode represented by this metadata.
    pub fn mode(&self) -> RefKmerRowMode {
        match self {
            Self::Global => RefKmerRowMode::Global,
            Self::Windows {
                window_mode: RefKmerWindowMode::Size,
                ..
            } => RefKmerRowMode::SizeWindows,
            Self::Windows {
                window_mode: RefKmerWindowMode::Bed,
                ..
            } => RefKmerRowMode::BedWindows,
            Self::Groups(_) => RefKmerRowMode::Groups,
        }
    }
}

/// Metadata for one grouped reference k-mer output row.
#[derive(Debug, Clone, PartialEq)]
pub struct RefKmerGroupRow {
    /// Zero-based row index in frequency-row order.
    pub index: usize,
    /// Public group name from the grouped BED input.
    pub name: String,
    /// Number of grouped BED windows contributing to the group.
    pub eligible_windows: u64,
    /// Length-weighted blacklist fraction across the group's windows.
    pub blacklisted_fraction: f64,
}

/// Native frequency storage read from a reference k-mer output.
#[derive(Debug, Clone, PartialEq)]
pub enum RefKmerFrequencyData {
    /// Dense frequency matrix with shape `(row, motif)`.
    Dense(DenseMatrix<f64>),
    /// Sparse COO frequency entries with dense shape metadata.
    Sparse(RefKmerSparseFrequencies),
}

impl RefKmerFrequencyData {
    /// Return the storage mode represented by this value.
    pub fn storage_mode(&self) -> RefKmerStorageMode {
        match self {
            Self::Dense(_) => RefKmerStorageMode::Dense,
            Self::Sparse(_) => RefKmerStorageMode::SparseCoo,
        }
    }

    /// Reconstruct counts directly from the stored frequency representation.
    ///
    /// Dense stores copy the dense frequency values and scale each row in place. Sparse stores
    /// allocate the final dense count matrix and write scaled stored entries directly, so missing
    /// sparse coordinates remain zero without a separate dense-frequency pass.
    fn to_dense_count_matrix(
        &self,
        row_scaling_factors: &[f64],
        row_offset_context: &'static str,
    ) -> OutputLoaderResult<DenseMatrix<f64>> {
        match self {
            Self::Dense(frequencies) => {
                let (row_count, motif_count) = frequencies.shape();
                dense_counts_from_row_major_frequencies(
                    frequencies.values_row_major().to_vec(),
                    row_count,
                    motif_count,
                    row_scaling_factors,
                    row_offset_context,
                )
            }
            Self::Sparse(sparse) => sparse.to_dense_count_matrix(row_scaling_factors),
        }
    }
}

/// Sparse COO reference k-mer frequency entries.
#[derive(Debug, Clone, PartialEq)]
pub struct RefKmerSparseFrequencies {
    row_count: usize,
    motif_count: usize,
    row_indices: Vec<usize>,
    motif_indices: Vec<usize>,
    frequencies: Vec<f64>,
}

impl RefKmerSparseFrequencies {
    /// Return the dense shape represented by the sparse entries.
    pub fn shape(&self) -> (usize, usize) {
        (self.row_count, self.motif_count)
    }

    /// Return zero-based row indices for stored COO entries.
    pub fn row_indices(&self) -> &[usize] {
        &self.row_indices
    }

    /// Return zero-based motif indices for stored COO entries.
    pub fn motif_indices(&self) -> &[usize] {
        &self.motif_indices
    }

    /// Return stored frequencies in COO entry order.
    pub fn frequencies(&self) -> &[f64] {
        &self.frequencies
    }

    /// Return stored frequencies as `(row, motif, frequency)` entries.
    pub fn entries(&self) -> impl ExactSizeIterator<Item = RefKmerSparseFrequencyEntry> + '_ {
        self.row_indices
            .iter()
            .copied()
            .zip(self.motif_indices.iter().copied())
            .zip(self.frequencies.iter().copied())
            .map(
                |((row_index, motif_index), frequency)| RefKmerSparseFrequencyEntry {
                    row_index,
                    motif_index,
                    frequency,
                },
            )
    }

    /// Return the number of stored COO entries.
    pub fn nnz(&self) -> usize {
        self.frequencies.len()
    }

    /// Return one frequency value, if both indices are in bounds.
    ///
    /// Missing in-bounds sparse coordinates are zero frequencies.
    pub fn frequency(&self, row_index: usize, motif_index: usize) -> Option<f64> {
        if row_index >= self.row_count || motif_index >= self.motif_count {
            return None;
        }
        let mut left_entry_index = 0;
        let mut right_entry_index = self.frequencies.len();
        while left_entry_index < right_entry_index {
            let middle_entry_index = left_entry_index + (right_entry_index - left_entry_index) / 2;
            let stored_coordinate = (
                self.row_indices[middle_entry_index],
                self.motif_indices[middle_entry_index],
            );
            match stored_coordinate.cmp(&(row_index, motif_index)) {
                std::cmp::Ordering::Less => left_entry_index = middle_entry_index + 1,
                std::cmp::Ordering::Equal => return Some(self.frequencies[middle_entry_index]),
                std::cmp::Ordering::Greater => right_entry_index = middle_entry_index,
            }
        }
        Some(0.0)
    }

    /// Reconstruct a dense frequency matrix from the stored COO entries.
    pub fn to_dense_matrix(&self) -> OutputLoaderResult<DenseMatrix<f64>> {
        let value_count = self
            .row_count
            .checked_mul(self.motif_count)
            .context("reference k-mer dense frequency shape overflow")?;
        let mut values = vec![0.0; value_count];
        for ((&row_index, &motif_index), &frequency) in self
            .row_indices
            .iter()
            .zip(&self.motif_indices)
            .zip(&self.frequencies)
        {
            let value_index = self.dense_value_index(row_index, motif_index)?;
            values[value_index] = frequency;
        }
        Ok(DenseMatrix::from_row_major(
            values,
            self.row_count,
            self.motif_count,
        )?)
    }

    /// Reconstruct a dense count matrix from the stored COO entries.
    ///
    /// This writes `frequency * row_scaling_factor[row]` directly into the returned dense matrix.
    /// Missing sparse coordinates stay zero, which avoids materializing and then scaling a dense
    /// frequency matrix first.
    fn to_dense_count_matrix(
        &self,
        row_scaling_factors: &[f64],
    ) -> OutputLoaderResult<DenseMatrix<f64>> {
        if row_scaling_factors.len() != self.row_count {
            return Err(OutputLoaderError::message(
                "reference k-mer row scaling factor count did not match sparse row count",
            ));
        }
        let value_count = self
            .row_count
            .checked_mul(self.motif_count)
            .context("reference k-mer dense count shape overflow")?;
        let mut values = vec![0.0; value_count];
        for ((&row_index, &motif_index), &frequency) in self
            .row_indices
            .iter()
            .zip(&self.motif_indices)
            .zip(&self.frequencies)
        {
            let value_index = self.dense_value_index(row_index, motif_index)?;
            values[value_index] = frequency * row_scaling_factors[row_index];
        }
        Ok(DenseMatrix::from_row_major(
            values,
            self.row_count,
            self.motif_count,
        )?)
    }

    /// Return the row-major dense index for a sparse coordinate.
    fn dense_value_index(&self, row_index: usize, motif_index: usize) -> OutputLoaderResult<usize> {
        if row_index >= self.row_count || motif_index >= self.motif_count {
            return Err(OutputLoaderError::message(format!(
                "reference k-mer sparse coordinate ({row_index}, {motif_index}) is outside dense shape ({}, {})",
                self.row_count, self.motif_count
            )));
        }
        Ok(row_index
            .checked_mul(self.motif_count)
            .and_then(|row_start| row_start.checked_add(motif_index))
            .context("reference k-mer sparse coordinate overflow")?)
    }

    fn select_frequencies(
        &self,
        row_indices: &[usize],
        motif_indices: &[usize],
        row_label: &str,
    ) -> Result<Self> {
        let row_index_map = build_selection_index_map(row_indices, self.row_count, row_label)?;
        let motif_index_map = build_selection_index_map(motif_indices, self.motif_count, "motif")?;

        let mut selected_entries = Vec::new();
        for entry_index in 0..self.frequencies.len() {
            let source_row_index = self.row_indices[entry_index];
            let source_motif_index = self.motif_indices[entry_index];
            let Some(selected_row_index) = row_index_map[source_row_index] else {
                continue;
            };
            let Some(selected_motif_index) = motif_index_map[source_motif_index] else {
                continue;
            };
            selected_entries.push((
                selected_row_index,
                selected_motif_index,
                self.frequencies[entry_index],
            ));
        }
        selected_entries.sort_by_key(|&(row_index, motif_index, _)| (row_index, motif_index));

        let mut selected_row_indices = Vec::with_capacity(selected_entries.len());
        let mut selected_motif_indices = Vec::with_capacity(selected_entries.len());
        let mut selected_frequencies = Vec::with_capacity(selected_entries.len());
        for (row_index, motif_index, frequency) in selected_entries {
            selected_row_indices.push(row_index);
            selected_motif_indices.push(motif_index);
            selected_frequencies.push(frequency);
        }

        Ok(Self {
            row_count: row_indices.len(),
            motif_count: motif_indices.len(),
            row_indices: selected_row_indices,
            motif_indices: selected_motif_indices,
            frequencies: selected_frequencies,
        })
    }

    fn count_entries(&self, row_scaling_factors: &[f64]) -> Result<Vec<RefKmerSparseCountEntry>> {
        ensure!(
            row_scaling_factors.len() == self.row_count,
            "reference k-mer row scaling factor count did not match sparse row count"
        );
        Ok(self
            .entries()
            .map(|entry| RefKmerSparseCountEntry {
                row_index: entry.row_index,
                motif_index: entry.motif_index,
                count: entry.frequency * row_scaling_factors[entry.row_index],
            })
            .collect())
    }
}

/// One stored entry from sparse reference k-mer COO frequencies.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RefKmerSparseFrequencyEntry {
    /// Zero-based row index in the dense matrix represented by the sparse frequencies.
    pub row_index: usize,
    /// Zero-based motif index in the dense matrix represented by the sparse frequencies.
    pub motif_index: usize,
    /// Stored frequency for `(row_index, motif_index)`.
    pub frequency: f64,
}

/// One reconstructed count entry from sparse reference k-mer COO frequencies.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RefKmerSparseCountEntry {
    /// Zero-based row index in the dense matrix represented by the sparse frequencies.
    pub row_index: usize,
    /// Zero-based motif index in the dense matrix represented by the sparse frequencies.
    pub motif_index: usize,
    /// Reconstructed count for `(row_index, motif_index)`.
    pub count: f64,
}

fn select_dense_ref_kmer_frequencies(
    frequencies: &DenseMatrix<f64>,
    row_indices: &[usize],
    motif_indices: &[usize],
    row_label: &str,
) -> Result<DenseMatrix<f64>> {
    let mut selected_values =
        Vec::with_capacity(row_indices.len().saturating_mul(motif_indices.len()));
    let contiguous_columns = contiguous_index_span(motif_indices);
    for &row_index in row_indices {
        ensure!(
            row_index < frequencies.row_count(),
            "{row_label} index {row_index} is outside 0..{}",
            frequencies.row_count()
        );
        if let Some((start_column, end_column)) = contiguous_columns {
            let contiguous_row_values = frequencies
                .row_values_for_column_range(row_index, start_column, end_column)
                .with_context(|| {
                    format!(
                        "motif range {start_column}..{end_column} is outside 0..{}",
                        frequencies.column_count()
                    )
                })?;
            selected_values.extend_from_slice(contiguous_row_values);
        } else {
            let row_values = frequencies.row(row_index).with_context(|| {
                format!(
                    "{row_label} index {row_index} is outside 0..{}",
                    frequencies.row_count()
                )
            })?;
            for &motif_index in motif_indices {
                let frequency = row_values.get(motif_index).copied().with_context(|| {
                    format!(
                        "motif index {motif_index} is outside 0..{}",
                        frequencies.column_count()
                    )
                })?;
                selected_values.push(frequency);
            }
        }
    }
    Ok(DenseMatrix::from_row_major(
        selected_values,
        row_indices.len(),
        motif_indices.len(),
    )?)
}

/// Scale dense row-major frequencies into dense row-major counts.
fn dense_counts_from_row_major_frequencies(
    mut values: Vec<f64>,
    row_count: usize,
    motif_count: usize,
    row_scaling_factors: &[f64],
    row_offset_context: &'static str,
) -> OutputLoaderResult<DenseMatrix<f64>> {
    if row_scaling_factors.len() != row_count {
        return Err(OutputLoaderError::message(
            "reference k-mer row scaling factor count did not match dense row count",
        ));
    }
    let expected_value_count = row_count
        .checked_mul(motif_count)
        .context("reference k-mer dense count shape overflow")?;
    if values.len() != expected_value_count {
        return Err(OutputLoaderError::message(format!(
            "reference k-mer dense count shape ({row_count}, {motif_count}) does not match {} row-major values",
            values.len()
        )));
    }
    for (row_index, &scaling_factor) in row_scaling_factors.iter().enumerate() {
        let row_start = row_index
            .checked_mul(motif_count)
            .context(row_offset_context)?;
        let row_end = row_start
            .checked_add(motif_count)
            .context(row_offset_context)?;
        for value in &mut values[row_start..row_end] {
            *value *= scaling_factor;
        }
    }
    Ok(DenseMatrix::from_row_major(values, row_count, motif_count)?)
}

fn resolve_ref_kmer_motif_indices(
    motif_indices: Option<&[usize]>,
    motif_count: usize,
) -> Vec<usize> {
    if let Some(motif_indices) = motif_indices {
        return motif_indices.to_vec();
    }
    (0..motif_count).collect()
}

fn resolve_ref_kmer_motif_label_indices<S: AsRef<str>>(
    output: &RefKmersOutput,
    motif_labels: &[S],
) -> Result<Vec<usize>> {
    ensure_unique_labels(motif_labels, "motif_labels")?;
    motif_labels
        .iter()
        .map(|motif_label| {
            output
                .motif_index(motif_label.as_ref())
                .map_err(anyhow::Error::from)
        })
        .collect()
}

fn resolve_ref_kmer_motif_labels<S: AsRef<str>>(
    output: &RefKmersOutput,
    motif_labels: Option<&[S]>,
) -> Result<Option<Vec<usize>>> {
    match motif_labels {
        Some(labels) => Ok(Some(resolve_ref_kmer_motif_label_indices(output, labels)?)),
        None => Ok(None),
    }
}

fn resolve_ref_kmer_group_name_indices<S: AsRef<str>>(
    output: &RefKmersOutput,
    group_names: &[S],
) -> Result<Vec<usize>> {
    ensure_unique_labels(group_names, "group_names")?;
    group_names
        .iter()
        .map(|group_name| {
            output
                .group_index(group_name.as_ref())
                .map_err(anyhow::Error::from)
        })
        .collect()
}

/// Parser for a `cfdna ref-kmers` Zarr store.
struct RefKmersParser {
    path: PathBuf,
}

impl RefKmersParser {
    /// Create a parser for one reference k-mer Zarr store path.
    fn new(path: &Path) -> Self {
        Self {
            path: path.to_path_buf(),
        }
    }

    /// Load and validate the reference k-mer Zarr store.
    fn load(&self) -> Result<RefKmersOutput> {
        validate_ref_kmer_zarr_store_path(&self.path)?;
        let root_attributes = read_zarr_root_attributes(&self.path).with_context(|| {
            format!("read reference k-mer Zarr metadata {}", self.path.display())
        })?;
        let root_metadata = RefKmerRootMetadata::from_attributes(&root_attributes)?;
        let store =
            Arc::new(FilesystemStore::new(&self.path).with_context(|| {
                format!("open reference k-mer Zarr store {}", self.path.display())
            })?);

        let motif_index = read_zarr_array1::<i32>(store.clone(), "/motif_index")?;
        validate_zero_based_axis(&motif_index, "motif_index")?;
        let motif_labels = read_motif_labels(
            &self.path,
            store.clone(),
            root_metadata.motif_axis_kind,
            motif_index.len(),
        )?;
        ensure_unique_labels(&motif_labels, "motif_labels")?;
        let mut motif_label_indices =
            FxHashMap::with_capacity_and_hasher(motif_labels.len(), Default::default());
        for (motif_index, motif_label) in motif_labels.iter().enumerate() {
            motif_label_indices.insert(motif_label.clone(), motif_index);
        }

        let row = read_zarr_array1::<i32>(store.clone(), "/row")?;
        validate_zero_based_axis(&row, "row")?;
        ensure!(!row.is_empty(), "reference k-mer output has no rows");
        let row_metadata =
            read_row_metadata(&self.path, store.clone(), root_metadata.row_mode, row.len())?;
        let row_scaling_factors = read_row_scaling_factors(store.clone(), row.len())?;
        let reference_contig_footprint = read_reference_contig_footprint(store.clone())?;

        let data =
            match root_metadata.storage_mode {
                RefKmerStorageMode::Dense => RefKmerFrequencyData::Dense(read_dense_frequencies(
                    store,
                    row.len(),
                    motif_labels.len(),
                )?),
                RefKmerStorageMode::SparseCoo => RefKmerFrequencyData::Sparse(
                    read_sparse_frequencies(&self.path, store, row.len(), motif_labels.len())?,
                ),
            };

        Ok(RefKmersOutput {
            row_metadata,
            motif_axis_kind: root_metadata.motif_axis_kind,
            motif_labels,
            motif_label_indices,
            row_scaling_factors,
            data,
            kmer_size: root_metadata.kmer_size,
            canonical: root_metadata.canonical,
            all_motifs: root_metadata.all_motifs,
            assign_by: root_metadata.assign_by,
            reference_contig_footprint,
        })
    }
}

/// Required root-level metadata from a reference k-mer Zarr store.
struct RefKmerRootMetadata {
    storage_mode: RefKmerStorageMode,
    row_mode: RefKmerRowMode,
    motif_axis_kind: RefKmerMotifAxisKind,
    kmer_size: u8,
    canonical: bool,
    all_motifs: bool,
    assign_by: String,
}

impl RefKmerRootMetadata {
    /// Parse and validate required root attributes from a Zarr metadata object.
    fn from_attributes(attributes: &Value) -> Result<Self> {
        ensure!(
            string_attr(attributes, "cfdnalab_schema")? == "ref_kmer_frequencies",
            "reference k-mer Zarr schema mismatch"
        );
        ensure!(
            u64_attr(attributes, "cfdnalab_schema_version")? == REF_KMER_SCHEMA_VERSION,
            "reference k-mer Zarr schema version mismatch: expected {}",
            REF_KMER_SCHEMA_VERSION
        );
        ensure!(
            string_attr(attributes, "value_units")? == "reference_kmer_frequency",
            "reference k-mer Zarr value_units must be reference_kmer_frequency"
        );
        ensure!(
            string_attr(attributes, "count_units")? == "reference_kmer_count",
            "reference k-mer Zarr count_units must be reference_kmer_count"
        );
        ensure!(
            string_attr(attributes, "row_scaling_factor_array")? == "row_scaling_factor",
            "reference k-mer Zarr row_scaling_factor_array must be row_scaling_factor"
        );

        let storage_mode = match string_attr(attributes, "storage_mode")? {
            "dense" => {
                ensure!(
                    string_attr(attributes, "primary_array")? == "frequencies",
                    "dense reference k-mer Zarr primary_array must be frequencies"
                );
                ensure!(
                    attributes.get("primary_group").is_some_and(Value::is_null),
                    "dense reference k-mer Zarr primary_group must be null"
                );
                RefKmerStorageMode::Dense
            }
            "sparse_coo" => {
                ensure!(
                    attributes.get("primary_array").is_some_and(Value::is_null),
                    "sparse reference k-mer Zarr primary_array must be null"
                );
                ensure!(
                    string_attr(attributes, "primary_group")? == "sparse",
                    "sparse reference k-mer Zarr primary_group must be sparse"
                );
                ensure!(
                    string_attr(attributes, "sparse_format")? == "coo",
                    "sparse reference k-mer Zarr sparse_format must be coo"
                );
                ensure!(
                    u64_attr(attributes, "sparse_indices_base")? == 0,
                    "sparse reference k-mer Zarr sparse_indices_base must be 0"
                );
                RefKmerStorageMode::SparseCoo
            }
            other => bail!("unsupported reference k-mer storage_mode '{other}'"),
        };
        let row_mode = match string_attr(attributes, "row_mode")? {
            "global" => RefKmerRowMode::Global,
            "size" => RefKmerRowMode::SizeWindows,
            "bed" => RefKmerRowMode::BedWindows,
            "grouped_bed" => RefKmerRowMode::Groups,
            other => bail!("unsupported reference k-mer row_mode '{other}'"),
        };
        let motif_axis_kind = match string_attr(attributes, "motif_axis_kind")? {
            "motif" => RefKmerMotifAxisKind::Motif,
            "motif_group" => RefKmerMotifAxisKind::MotifGroup,
            other => bail!("unsupported reference k-mer motif_axis_kind '{other}'"),
        };
        let kmer_size = u8_from_u64(u64_attr(attributes, "kmer_size")?, "kmer_size")?;
        ensure!(kmer_size > 0, "reference k-mer kmer_size must be positive");

        Ok(Self {
            storage_mode,
            row_mode,
            motif_axis_kind,
            kmer_size,
            canonical: bool_attr(attributes, "canonical")?,
            all_motifs: bool_attr(attributes, "all_motifs")?,
            assign_by: string_attr(attributes, "assign_by")?.to_string(),
        })
    }
}

fn validate_ref_kmer_zarr_store_path(path: &Path) -> Result<()> {
    ensure!(
        path.exists(),
        "reference k-mer Zarr store does not exist: {}",
        path.display()
    );
    ensure!(
        path.is_dir(),
        "reference k-mer Zarr store is not a directory: {}",
        path.display()
    );
    ensure!(
        path.file_name()
            .and_then(|name| name.to_str())
            .is_some_and(|name| name.ends_with(".zarr")),
        "reference k-mer Zarr store path must end with .zarr: {}",
        path.display()
    );
    ensure!(
        path.join("zarr.json").is_file(),
        "reference k-mer Zarr store is missing root zarr.json: {}",
        path.display()
    );
    Ok(())
}

fn read_reference_contig_footprint(
    store: Arc<FilesystemStore>,
) -> Result<Vec<ContigFootprintEntry>> {
    let reference_contig_footprint_json =
        read_zarr_array1::<u8>(store, "/reference_contig_footprint_json")?;
    serde_json::from_slice(&reference_contig_footprint_json)
        .context("invalid reference_contig_footprint_json in reference k-mer package")
}

fn read_motif_labels(
    root_path: &Path,
    store: Arc<FilesystemStore>,
    motif_axis_kind: RefKmerMotifAxisKind,
    motif_count: usize,
) -> Result<Vec<String>> {
    match motif_axis_kind {
        RefKmerMotifAxisKind::Motif => read_motif_ascii_labels(store, motif_count),
        RefKmerMotifAxisKind::MotifGroup => {
            read_zarr_labels(root_path, "motif_index", "motif_group", motif_count)
        }
    }
}

fn read_motif_ascii_labels(store: Arc<FilesystemStore>, motif_count: usize) -> Result<Vec<String>> {
    let motif_byte = read_zarr_array1::<i32>(store.clone(), "/motif_byte")?;
    validate_zero_based_axis(&motif_byte, "motif_byte")?;
    let motif_width = motif_byte.len();
    ensure!(
        motif_width > 0 || motif_count == 0,
        "motif_ascii cannot decode non-empty motif axis with zero motif_byte width"
    );

    let (bytes, shape) = read_zarr_array_values::<u8>(store, "/motif_ascii")?;
    ensure!(
        shape == [motif_count, motif_width],
        "motif_ascii shape {:?} did not match expected [{}, {}]",
        shape,
        motif_count,
        motif_width
    );
    if motif_width == 0 {
        return Ok(Vec::new());
    }
    bytes
        .chunks_exact(motif_width)
        .enumerate()
        .map(|(motif_index, motif_bytes)| {
            ensure!(
                motif_bytes.is_ascii(),
                "motif_ascii row {motif_index} contains non-ASCII motif bytes"
            );
            String::from_utf8(motif_bytes.to_vec()).context("motif_ascii contains invalid UTF-8")
        })
        .collect()
}

fn read_row_metadata(
    root_path: &Path,
    store: Arc<FilesystemStore>,
    row_mode: RefKmerRowMode,
    row_count: usize,
) -> Result<RefKmerRowMetadata> {
    match row_mode {
        RefKmerRowMode::Global => {
            let labels = read_zarr_labels(root_path, "row", "row_label", row_count)?;
            ensure!(
                labels == ["global"],
                "global reference k-mer row labels must be [\"global\"]"
            );
            Ok(RefKmerRowMetadata::Global)
        }
        RefKmerRowMode::SizeWindows | RefKmerRowMode::BedWindows => {
            read_window_row_metadata(root_path, store, row_mode, row_count)
        }
        RefKmerRowMode::Groups => read_group_row_metadata(root_path, store, row_count),
    }
}

fn read_window_row_metadata(
    root_path: &Path,
    store: Arc<FilesystemStore>,
    row_mode: RefKmerRowMode,
    row_count: usize,
) -> Result<RefKmerRowMetadata> {
    let chromosome = read_zarr_array1::<i32>(store.clone(), "/chromosome")?;
    validate_zero_based_axis(&chromosome, "chromosome")?;
    let chromosome_names =
        read_zarr_labels(root_path, "chromosome", "chromosome_name", chromosome.len())?;

    let row_chromosome = read_zarr_array1::<i32>(store.clone(), "/row_chromosome")?;
    let row_start_bp = read_zarr_array1::<i64>(store.clone(), "/row_start_bp")?;
    let row_end_bp = read_zarr_array1::<i64>(store.clone(), "/row_end_bp")?;
    let blacklisted_fraction = read_zarr_array1::<f64>(store, "/blacklisted_fraction")?;
    ensure_same_len(&row_chromosome, row_count, "row_chromosome")?;
    ensure_same_len(&row_start_bp, row_count, "row_start_bp")?;
    ensure_same_len(&row_end_bp, row_count, "row_end_bp")?;
    ensure_same_len(&blacklisted_fraction, row_count, "blacklisted_fraction")?;

    let mut windows = Vec::with_capacity(row_count);
    for row_index in 0..row_count {
        let chromosome_index = usize_from_i32(row_chromosome[row_index], "row_chromosome")?;
        let chrom = chromosome_names
            .get(chromosome_index)
            .with_context(|| {
                format!("row_chromosome {chromosome_index} is outside the chromosome axis")
            })?
            .clone();
        let start = u64_from_i64(row_start_bp[row_index], "row_start_bp")?;
        let end = u64_from_i64(row_end_bp[row_index], "row_end_bp")?;
        let interval = Interval::new(start, end).map_err(|error| {
            anyhow::anyhow!("reference k-mer row {row_index} has invalid interval: {error}")
        })?;
        let blacklisted_fraction = blacklisted_fraction[row_index];
        ensure_valid_fraction(blacklisted_fraction, "blacklisted_fraction")?;
        windows.push(WindowRow {
            index: row_index,
            chrom,
            interval,
            blacklisted_fraction: Some(blacklisted_fraction),
        });
    }

    let window_mode = match row_mode {
        RefKmerRowMode::SizeWindows => RefKmerWindowMode::Size,
        RefKmerRowMode::BedWindows => RefKmerWindowMode::Bed,
        _ => bail!("internal reference k-mer loader row-mode mismatch"),
    };
    Ok(RefKmerRowMetadata::Windows {
        window_mode,
        windows,
    })
}

fn read_group_row_metadata(
    root_path: &Path,
    store: Arc<FilesystemStore>,
    row_count: usize,
) -> Result<RefKmerRowMetadata> {
    let group = read_zarr_array1::<i32>(store.clone(), "/group")?;
    validate_zero_based_axis(&group, "group")?;
    ensure_same_len(&group, row_count, "group")?;

    let group_names = read_zarr_labels(root_path, "group", "group_name", row_count)?;
    ensure_unique_labels(&group_names, "group_names")?;

    let eligible_windows = read_zarr_array1::<i32>(store.clone(), "/eligible_windows")?;
    let blacklisted_fraction = read_zarr_array1::<f64>(store, "/blacklisted_fraction")?;
    ensure_same_len(&eligible_windows, row_count, "eligible_windows")?;
    ensure_same_len(&blacklisted_fraction, row_count, "blacklisted_fraction")?;

    let groups = (0..row_count)
        .map(|row_index| {
            let eligible_windows = u64_from_i32(eligible_windows[row_index], "eligible_windows")?;
            let blacklisted_fraction = blacklisted_fraction[row_index];
            ensure_valid_fraction(blacklisted_fraction, "blacklisted_fraction")?;
            Ok(RefKmerGroupRow {
                index: row_index,
                name: group_names[row_index].clone(),
                eligible_windows,
                blacklisted_fraction,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(RefKmerRowMetadata::Groups(groups))
}

fn read_row_scaling_factors(store: Arc<FilesystemStore>, row_count: usize) -> Result<Vec<f64>> {
    let row_scaling_factors = read_zarr_array1::<f64>(store, "/row_scaling_factor")?;
    ensure_same_len(&row_scaling_factors, row_count, "row_scaling_factor")?;
    for (row_index, &scaling_factor) in row_scaling_factors.iter().enumerate() {
        ensure!(
            scaling_factor.is_finite() && scaling_factor >= 0.0,
            "row_scaling_factor contains value outside finite and non-negative range at row {row_index}: {scaling_factor}"
        );
    }
    Ok(row_scaling_factors)
}

fn read_dense_frequencies(
    store: Arc<FilesystemStore>,
    row_count: usize,
    motif_count: usize,
) -> Result<DenseMatrix<f64>> {
    let (values, shape) = read_zarr_array_values::<f64>(store, "/frequencies")?;
    ensure!(
        shape == [row_count, motif_count],
        "dense reference k-mer frequencies shape {:?} did not match expected [{}, {}]",
        shape,
        row_count,
        motif_count
    );
    ensure_valid_frequencies(&values, "dense reference k-mer frequencies")?;
    DenseMatrix::from_row_major(values, row_count, motif_count)
}

fn read_sparse_frequencies(
    root_path: &Path,
    store: Arc<FilesystemStore>,
    row_count: usize,
    motif_count: usize,
) -> Result<RefKmerSparseFrequencies> {
    let sparse_dimension = read_zarr_array1::<i32>(store.clone(), "/sparse/sparse_dimension")?;
    ensure!(
        sparse_dimension == [0, 1],
        "sparse/sparse_dimension must be [0, 1]"
    );
    let sparse_dimension_labels = read_zarr_labels(
        root_path,
        "sparse/sparse_dimension",
        "sparse_dimension_name",
        2,
    )?;
    ensure!(
        sparse_dimension_labels == ["row", "motif"],
        "sparse dimension labels must be [\"row\", \"motif\"]"
    );

    let sparse_shape = read_zarr_array1::<i32>(store.clone(), "/sparse/shape")?;
    ensure!(
        sparse_shape.len() == 2,
        "sparse/shape must have length 2, found {}",
        sparse_shape.len()
    );
    let sparse_row_count = usize_from_i32(sparse_shape[0], "sparse row count")?;
    let sparse_motif_count = usize_from_i32(sparse_shape[1], "sparse motif count")?;
    ensure!(
        (sparse_row_count, sparse_motif_count) == (row_count, motif_count),
        "sparse shape ({sparse_row_count}, {sparse_motif_count}) did not match metadata ({row_count}, {motif_count})"
    );

    let row = read_zarr_array1::<i32>(store.clone(), "/sparse/row")?;
    let motif = read_zarr_array1::<i32>(store.clone(), "/sparse/motif")?;
    let frequencies = read_zarr_array1::<f64>(store, "/sparse/frequency")?;
    ensure_same_len(&motif, row.len(), "sparse/motif")?;
    ensure_same_len(&frequencies, row.len(), "sparse/frequency")?;
    ensure_valid_frequencies(&frequencies, "sparse reference k-mer frequencies")?;

    let mut row_indices = Vec::with_capacity(row.len());
    let mut motif_indices = Vec::with_capacity(motif.len());
    let mut previous_coordinate = None;
    for entry_index in 0..row.len() {
        let row_index = usize_from_i32(row[entry_index], "sparse/row")?;
        let motif_index = usize_from_i32(motif[entry_index], "sparse/motif")?;
        ensure!(
            row_index < row_count,
            "sparse row index {row_index} is outside 0..{row_count}"
        );
        ensure!(
            motif_index < motif_count,
            "sparse motif index {motif_index} is outside 0..{motif_count}"
        );
        let coordinate = (row_index, motif_index);
        if let Some(previous_coordinate) = previous_coordinate {
            ensure!(
                previous_coordinate < coordinate,
                "sparse reference k-mer COO entries must be sorted and unique"
            );
        }
        previous_coordinate = Some(coordinate);
        row_indices.push(row_index);
        motif_indices.push(motif_index);
    }

    Ok(RefKmerSparseFrequencies {
        row_count,
        motif_count,
        row_indices,
        motif_indices,
        frequencies,
    })
}

fn ensure_valid_frequencies(frequencies: &[f64], label: &str) -> Result<()> {
    for (frequency_index, &frequency) in frequencies.iter().enumerate() {
        ensure!(
            frequency.is_finite() && (0.0..=1.0).contains(&frequency),
            "{label} contain value outside [0, 1] at index {frequency_index}: {frequency}"
        );
    }
    Ok(())
}

fn read_zarr_array_values<T>(
    store: Arc<FilesystemStore>,
    array_path: &str,
) -> Result<(Vec<T>, Vec<usize>)>
where
    T: ElementOwned,
{
    let array = Array::open(store, array_path)?;
    let shape = array
        .shape()
        .iter()
        .map(|dimension| usize::try_from(*dimension).context("Zarr dimension exceeds usize"))
        .collect::<Result<Vec<_>>>()?;
    let values = array
        .retrieve_array_subset(&array.subset_all())
        .with_context(|| format!("read Zarr array {array_path}"))?;
    Ok((values, shape))
}

fn read_zarr_array1<T>(store: Arc<FilesystemStore>, array_path: &str) -> Result<Vec<T>>
where
    T: ElementOwned,
{
    let (values, shape) = read_zarr_array_values(store, array_path)?;
    ensure!(
        shape.len() == 1,
        "Zarr array {array_path} must be rank 1, found rank {}",
        shape.len()
    );
    Ok(values)
}

fn read_zarr_labels(
    root_path: &Path,
    array_path: &str,
    expected_label_field: &str,
    expected_len: usize,
) -> Result<Vec<String>> {
    let attributes = read_zarr_array_attributes(root_path, array_path)?;
    ensure!(
        string_attr(&attributes, "label_field")? == expected_label_field,
        "Zarr array {array_path} label_field must be {expected_label_field}"
    );
    let labels = attributes
        .get("labels")
        .and_then(Value::as_array)
        .with_context(|| format!("Zarr array {array_path} is missing labels"))?;
    ensure!(
        labels.len() == expected_len,
        "Zarr array {array_path} has {} labels, expected {expected_len}",
        labels.len()
    );
    labels
        .iter()
        .map(|label| {
            let label = label
                .as_str()
                .with_context(|| format!("Zarr array {array_path} label is not a string"))?;
            validate_zarr_public_label(label, expected_label_field)?;
            Ok(label.to_string())
        })
        .collect()
}

fn read_zarr_array_attributes(root_path: &Path, array_path: &str) -> Result<Value> {
    let metadata_path = zarr_metadata_path(root_path, array_path);
    let metadata: Value = serde_json::from_str(
        &fs::read_to_string(&metadata_path)
            .with_context(|| format!("read Zarr metadata {}", metadata_path.display()))?,
    )
    .with_context(|| format!("parse Zarr metadata {}", metadata_path.display()))?;
    metadata.get("attributes").cloned().with_context(|| {
        format!(
            "Zarr metadata {} is missing attributes",
            metadata_path.display()
        )
    })
}

fn zarr_metadata_path(root_path: &Path, array_path: &str) -> PathBuf {
    let mut path = root_path.to_path_buf();
    for component in array_path.trim_start_matches('/').split('/') {
        path.push(component);
    }
    path.join("zarr.json")
}

fn validate_zero_based_axis(values: &[i32], axis_name: &str) -> Result<()> {
    for (index, &value) in values.iter().enumerate() {
        ensure!(
            value == i32::try_from(index).context("axis index exceeds i32")?,
            "{axis_name} must be a zero-based coordinate axis"
        );
    }
    Ok(())
}

fn ensure_same_len<T>(values: &[T], expected_len: usize, array_name: &str) -> Result<()> {
    ensure!(
        values.len() == expected_len,
        "{array_name} has {} entries, expected {expected_len}",
        values.len()
    );
    Ok(())
}

fn ensure_valid_fraction(value: f64, field_name: &str) -> Result<()> {
    ensure!(
        value.is_finite() && (0.0..=1.0).contains(&value),
        "{field_name} must be a finite fraction in [0, 1], got {value}"
    );
    Ok(())
}

fn string_attr<'a>(attributes: &'a Value, name: &str) -> Result<&'a str> {
    attributes
        .get(name)
        .and_then(Value::as_str)
        .with_context(|| {
            format!("reference k-mer Zarr metadata is missing string attribute {name}")
        })
}

fn u64_attr(attributes: &Value, name: &str) -> Result<u64> {
    attributes
        .get(name)
        .and_then(Value::as_u64)
        .with_context(|| {
            format!("reference k-mer Zarr metadata is missing integer attribute {name}")
        })
}

fn bool_attr(attributes: &Value, name: &str) -> Result<bool> {
    attributes
        .get(name)
        .and_then(Value::as_bool)
        .with_context(|| format!("reference k-mer Zarr metadata is missing bool attribute {name}"))
}

fn usize_from_i32(value: i32, field_name: &str) -> Result<usize> {
    ensure!(value >= 0, "{field_name} must be non-negative, got {value}");
    Ok(usize::try_from(value).expect("non-negative i32 always fits usize"))
}

fn u64_from_i32(value: i32, field_name: &str) -> Result<u64> {
    ensure!(value >= 0, "{field_name} must be non-negative, got {value}");
    Ok(u64::try_from(value).expect("non-negative i32 always fits u64"))
}

fn u64_from_i64(value: i64, field_name: &str) -> Result<u64> {
    ensure!(value >= 0, "{field_name} must be non-negative, got {value}");
    Ok(u64::try_from(value).expect("non-negative i64 always fits u64"))
}

fn u8_from_u64(value: u64, field_name: &str) -> Result<u8> {
    value
        .try_into()
        .with_context(|| format!("{field_name} value {value} exceeds u8"))
}
