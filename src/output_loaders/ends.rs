//! Loader for `cfdna ends` end-motif Zarr outputs.
//!
//! End-motif outputs are Zarr stores with a row axis and a motif axis. Rows are
//! global, genomic windows, or grouped-BED groups. The motif axis contains
//! either concrete end motifs or motif-group labels from `--motifs-file`.
//!
//! The loader reads and validates the store metadata eagerly. Dense count
//! stores are read into a `DenseMatrix<f64>`. Sparse stores are read as COO
//! vectors. `select()` returns a selector builder for rows and motifs.
//!
//! Type overview:
//!
//! ```text
//! load_ends_output(path)
//!     -> EndsOutput
//!         row_metadata: EndMotifRowMetadata
//!             Global
//!             Windows { window_mode, windows: Vec<WindowRow> }
//!             Groups(Vec<EndMotifGroupRow>)
//!         motif_axis_kind: EndMotifAxisKind
//!         motif_labels: Vec<String>
//!         data: EndMotifCountsData
//!             Dense(DenseMatrix<f64>)
//!             Sparse(EndMotifSparseCounts)
//!
//! EndsOutput::select()
//!     -> EndsSelector
//!         -> read()
//!             -> EndMotifCountSelection
//!                 row_metadata: EndMotifRowMetadata
//!                 row_indices: Vec<usize>
//!                 motif_indices: Vec<usize>
//!                 motif_labels: Vec<String>
//!                 data: EndMotifCountsData
//! ```
//!
//! Selections preserve the original storage mode, and sparse selections are
//! only densified when `to_dense_matrix()` is called explicitly. Selections
//! include selected row and motif metadata next to selected counts, so count
//! rows can be paired with their windows or groups without re-reading the Zarr
//! metadata.

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
    shared::zarr::read_zarr_root_attributes,
};
use anyhow::{Context, Result, bail, ensure};
use fxhash::FxHashMap;
use serde_json::Value;
use std::{
    fmt, fs,
    path::{Path, PathBuf},
    sync::Arc,
};
use zarrs::{
    array::{Array, ElementOwned},
    filesystem::FilesystemStore,
};

const END_MOTIF_SCHEMA_VERSION: u64 = 2;

/// Load a `cfdna ends` end-motif Zarr store.
///
/// The path must point to a Zarr directory written with the supported
/// end-motif schema version. The loader reads row metadata, motif labels, and
/// either dense counts or sparse COO counts into owned Rust containers.
///
/// Parameters
/// ----------
/// - `path`:
///   Path to a `cfdna ends` Zarr output directory.
///
/// Returns
/// -------
/// - `EndsOutput`:
///   Loaded row metadata, motif labels, and dense or sparse count storage.
///
/// ```no_run
/// use cfdnalab::output_loaders::load_ends_output;
///
/// let ends = load_ends_output("sample.end_motifs.zarr")?;
/// let selected = ends
///     .select()
///     .groups_by_name(&["promoter", "enhancer"])
///     .motifs_by_label(&["_AA", "_TT"])
///     .read()?;
/// let dense_counts = selected.to_dense_matrix()?;
///
/// for (group, motif_counts) in selected.group_metadata()?.iter().zip(dense_counts.rows()) {
///     let motif_total = motif_counts.iter().copied().sum::<f64>();
///     println!("{}: {motif_total}", group.name);
/// }
///
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn load_ends_output(path: impl AsRef<Path>) -> OutputLoaderResult<EndsOutput> {
    EndsParser::new(path.as_ref()).load().map_err(Into::into)
}

/// Loaded end-motif counts from `cfdna ends`.
#[derive(Debug, Clone, PartialEq)]
pub struct EndsOutput {
    row_metadata: EndMotifRowMetadata,
    motif_axis_kind: EndMotifAxisKind,
    motif_labels: Vec<String>,
    motif_label_indices: FxHashMap<String, usize>,
    data: EndMotifCountsData,
}

impl EndsOutput {
    /// Return how counts were stored on disk.
    pub fn storage_mode(&self) -> EndMotifStorageMode {
        self.data.storage_mode()
    }

    /// Return what each count row represents.
    pub fn row_mode(&self) -> EndMotifRowMode {
        self.row_metadata.mode()
    }

    /// Return whether motif-axis labels are concrete motifs or motif groups.
    pub fn motif_axis_kind(&self) -> EndMotifAxisKind {
        self.motif_axis_kind
    }

    /// Return motif-axis labels in count-column order.
    pub fn motif_labels(&self) -> &[String] {
        &self.motif_labels
    }

    /// Return a compact description of the loaded end-motif output.
    ///
    /// This combines storage mode, row mode, motif-axis kind, row count, and
    /// motif count in one value for logging or quick checks.
    pub fn output_metadata(&self) -> EndMotifOutputMetadata {
        EndMotifOutputMetadata {
            storage_mode: self.storage_mode(),
            row_mode: self.row_mode(),
            motif_axis_kind: self.motif_axis_kind(),
            row_count: self.row_count(),
            motif_count: self.motif_count(),
        }
    }

    /// Return the number of count rows.
    pub fn row_count(&self) -> usize {
        match &self.data {
            EndMotifCountsData::Dense(counts) => counts.shape().0,
            EndMotifCountsData::Sparse(sparse) => sparse.shape().0,
        }
    }

    /// Return the number of motif columns.
    pub fn motif_count(&self) -> usize {
        self.motif_labels.len()
    }

    /// Return row metadata describing the count rows.
    pub fn row_metadata(&self) -> &EndMotifRowMetadata {
        &self.row_metadata
    }

    /// Return window metadata, or an error if this is not a windowed output.
    pub fn window_metadata(&self) -> OutputLoaderResult<&[WindowRow]> {
        match &self.row_metadata {
            EndMotifRowMetadata::Windows { windows, .. } => Ok(windows),
            _ => Err(OutputLoaderError::message(
                "end-motif output is not windowed",
            )),
        }
    }

    /// Return group metadata, or an error if this is not a grouped output.
    pub fn group_metadata(&self) -> OutputLoaderResult<&[EndMotifGroupRow]> {
        match &self.row_metadata {
            EndMotifRowMetadata::Groups(groups) => Ok(groups),
            _ => Err(OutputLoaderError::message(
                "end-motif output is not grouped",
            )),
        }
    }

    /// Return one window row by zero-based row index.
    ///
    /// This returns an error if the loaded output is not windowed.
    ///
    /// Parameters
    /// ----------
    /// - `row_index`:
    ///   Zero-based row index in the window metadata.
    pub fn window(&self, row_index: usize) -> OutputLoaderResult<Option<&WindowRow>> {
        Ok(self.window_metadata()?.get(row_index))
    }

    /// Return one group row by zero-based row index.
    ///
    /// This returns an error if the loaded output is not grouped.
    ///
    /// Parameters
    /// ----------
    /// - `row_index`:
    ///   Zero-based row index in the group metadata.
    pub fn group(&self, row_index: usize) -> OutputLoaderResult<Option<&EndMotifGroupRow>> {
        Ok(self.group_metadata()?.get(row_index))
    }

    /// Return the group row index for one group name.
    ///
    /// Group names are expected to identify one row. The loader reports an
    /// error if the output is not grouped, the name is missing, or the file
    /// contains duplicate group names.
    ///
    /// Parameters
    /// ----------
    /// - `group_name`:
    ///   Group label to resolve to a zero-based row index.
    pub fn group_index(&self, group_name: &str) -> OutputLoaderResult<usize> {
        let groups = self.group_metadata()?;
        Ok(groups
            .iter()
            .filter(|group| group.name == group_name)
            .map(|group| group.index)
            .next()
            .with_context(|| format!("end-motif output has no group named '{group_name}'"))?)
    }

    /// Return whether one group name exists in a grouped output.
    ///
    /// This returns `false` for non-grouped outputs.
    ///
    /// Parameters
    /// ----------
    /// - `group_name`:
    ///   Group label to look up.
    pub fn has_group(&self, group_name: &str) -> bool {
        self.group_metadata()
            .is_ok_and(|groups| groups.iter().any(|group| group.name == group_name))
    }

    /// Return the motif index for one label.
    ///
    /// Parameters
    /// ----------
    /// - `motif_label`:
    ///   Motif or motif-group label to resolve to a zero-based motif index.
    pub fn motif_index(&self, motif_label: &str) -> OutputLoaderResult<usize> {
        Ok(self
            .motif_label_indices
            .get(motif_label)
            .copied()
            .with_context(|| format!("end-motif output has no motif label '{motif_label}'"))?)
    }

    /// Return whether one motif label exists.
    ///
    /// Parameters
    /// ----------
    /// - `motif_label`:
    ///   Motif or motif-group label to look up.
    pub fn has_motif(&self, motif_label: &str) -> bool {
        self.motif_label_indices.contains_key(motif_label)
    }

    /// Return one count value, if both indices are in bounds.
    ///
    /// Sparse outputs return `0.0` for in-bounds entries that are not stored in
    /// the COO vectors. Each sparse call does one binary search over sorted COO
    /// entries. For repeated arbitrary cell lookups, build a reusable lookup
    /// once with `sparse_counts()?.to_lookup_index()`. For block-style access,
    /// use `select()` or `sparse_counts()?.entries()`.
    ///
    /// Parameters
    /// ----------
    /// - `row_index`:
    ///   Zero-based count row index.
    /// - `motif_index`:
    ///   Zero-based motif axis index.
    pub fn count(&self, row_index: usize, motif_index: usize) -> Option<f64> {
        match &self.data {
            EndMotifCountsData::Dense(counts) => counts.get(row_index, motif_index).copied(),
            EndMotifCountsData::Sparse(sparse) => sparse.count(row_index, motif_index),
        }
    }

    /// Return one count value by motif label, if the row index is in bounds.
    ///
    /// Parameters
    /// ----------
    /// - `row_index`:
    ///   Zero-based count row index.
    /// - `motif_label`:
    ///   Motif or motif-group label to resolve before reading the count.
    pub fn count_for_motif(
        &self,
        row_index: usize,
        motif_label: &str,
    ) -> OutputLoaderResult<Option<f64>> {
        Ok(self.count(row_index, self.motif_index(motif_label)?))
    }

    /// Return dense counts, or an error if this store is sparse.
    pub fn dense_counts(&self) -> OutputLoaderResult<&DenseMatrix<f64>> {
        match &self.data {
            EndMotifCountsData::Dense(counts) => Ok(counts),
            EndMotifCountsData::Sparse(_) => {
                Err(OutputLoaderError::message("end-motif output is sparse"))
            }
        }
    }

    /// Return sparse COO counts, or an error if this store is dense.
    pub fn sparse_counts(&self) -> OutputLoaderResult<&EndMotifSparseCounts> {
        match &self.data {
            EndMotifCountsData::Sparse(sparse) => Ok(sparse),
            EndMotifCountsData::Dense(_) => {
                Err(OutputLoaderError::message("end-motif output is dense"))
            }
        }
    }

    /// Return count storage as either a dense matrix or sparse COO entries.
    ///
    /// The variant matches the storage mode of the loaded Zarr store.
    pub fn data(&self) -> &EndMotifCountsData {
        &self.data
    }

    /// Start a count selection.
    ///
    /// A new selector initially selects all rows and all motifs. Add row and
    /// motif constraints before calling `read()`.
    pub fn select(&self) -> EndsSelector<'_> {
        EndsSelector::new(self)
    }

    /// Start a reference-corrected count selection.
    ///
    /// This is available when both `cmd_ends` and `cmd_ref_kmers` features are
    /// enabled. The selector mirrors `select()` and returns corrected counts in
    /// the same dense or sparse shape as an ordinary end-motif count selection.
    #[cfg(feature = "cmd_ref_kmers")]
    pub fn select_corrected_counts<'a>(
        &'a self,
        ref_kmers: &'a crate::output_loaders::RefKmersOutput,
    ) -> crate::output_loaders::CorrectedEndMotifCountsSelector<'a> {
        crate::output_loaders::CorrectedEndMotifCountsSelector::new(self, ref_kmers)
    }

    /// Return selected rows and motifs while preserving the output storage mode.
    ///
    /// Passing `None` for `row_indices` selects all rows. Passing `None` for
    /// `motif_indices` selects all motifs. Dense outputs return dense
    /// selections, and sparse outputs return sparse COO selections.
    ///
    /// For sparse output, selection copies only stored COO entries whose source
    /// row and motif are both selected. Stored entries outside the selected
    /// axes are skipped. Missing in-bounds cells remain implicit zero counts.
    pub(crate) fn select_counts(
        &self,
        row_indices: Option<&[usize]>,
        motif_indices: Option<&[usize]>,
    ) -> Result<EndMotifCountSelection> {
        self.select_counts_with_label(row_indices, motif_indices, "row")
    }

    /// Return selected window rows and motifs.
    ///
    /// This returns an error if the loaded output is not windowed.
    pub(crate) fn select_window_counts(
        &self,
        window_indices: Option<&[usize]>,
        motif_indices: Option<&[usize]>,
    ) -> Result<EndMotifCountSelection> {
        ensure!(
            matches!(
                self.row_mode(),
                EndMotifRowMode::SizeWindows | EndMotifRowMode::BedWindows
            ),
            "end-motif output is not windowed"
        );
        self.select_counts_with_label(window_indices, motif_indices, "window")
    }

    /// Return selected group rows and motifs.
    ///
    /// This returns an error if the loaded output is not grouped.
    pub(crate) fn select_group_counts(
        &self,
        group_indices: Option<&[usize]>,
        motif_indices: Option<&[usize]>,
    ) -> Result<EndMotifCountSelection> {
        ensure!(
            self.row_mode() == EndMotifRowMode::Groups,
            "end-motif output is not grouped"
        );
        self.select_counts_with_label(group_indices, motif_indices, "group")
    }

    /// Return selected group names and motifs.
    ///
    /// Passing `None` for `group_names` selects all groups.
    pub(crate) fn select_group_counts_by_name<S: AsRef<str>>(
        &self,
        group_names: Option<&[S]>,
        motif_indices: Option<&[usize]>,
    ) -> Result<EndMotifCountSelection> {
        let group_indices = match group_names {
            Some(group_names) => Some(resolve_group_name_indices(self, group_names)?),
            None => {
                ensure!(
                    self.row_mode() == EndMotifRowMode::Groups,
                    "end-motif output is not grouped"
                );
                None
            }
        };
        self.select_counts_with_label(group_indices.as_deref(), motif_indices, "group")
    }

    /// Return selected rows and motif labels.
    ///
    /// Passing `None` for `motif_labels` selects all motifs.
    pub(crate) fn select_counts_by_motif_label<S: AsRef<str>>(
        &self,
        row_indices: Option<&[usize]>,
        motif_labels: Option<&[S]>,
    ) -> Result<EndMotifCountSelection> {
        let motif_indices = match motif_labels {
            Some(motif_labels) => Some(resolve_motif_label_indices(self, motif_labels)?),
            None => None,
        };
        self.select_counts_with_label(row_indices, motif_indices.as_deref(), "row")
    }

    /// Select rows and motifs after resolving axis defaults and labels.
    fn select_counts_with_label(
        &self,
        row_indices: Option<&[usize]>,
        motif_indices: Option<&[usize]>,
        row_label: &str,
    ) -> Result<EndMotifCountSelection> {
        ensure!(
            row_indices.is_none() || self.row_mode() != EndMotifRowMode::Global,
            "global end-motif output has no selectable row axis"
        );
        let row_indices = resolve_row_indices(row_indices, self.row_count(), row_label)?;
        let motif_indices = resolve_motif_indices(motif_indices, self.motif_count())?;
        ensure_unique_indices(&row_indices, row_label)?;
        ensure_unique_indices(&motif_indices, "motif")?;
        let selected_row_metadata = self.selected_row_metadata(&row_indices, row_label)?;

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

        let data =
            match &self.data {
                EndMotifCountsData::Dense(counts) => EndMotifCountsData::Dense(
                    select_dense_counts(counts, &row_indices, &motif_indices, row_label)?,
                ),
                EndMotifCountsData::Sparse(sparse) => EndMotifCountsData::Sparse(
                    sparse.select_counts(&row_indices, &motif_indices, row_label)?,
                ),
            };
        Ok(EndMotifCountSelection {
            row_metadata: selected_row_metadata,
            row_indices,
            motif_indices,
            motif_labels,
            data,
        })
    }

    /// Build row metadata for selected source row indices.
    ///
    /// The returned metadata keeps the selector order and preserves whether the
    /// loaded output is global, windowed, or grouped.
    fn selected_row_metadata(
        &self,
        row_indices: &[usize],
        row_label: &str,
    ) -> Result<EndMotifRowMetadata> {
        match &self.row_metadata {
            EndMotifRowMetadata::Global => Ok(EndMotifRowMetadata::Global),
            EndMotifRowMetadata::Windows {
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
                Ok(EndMotifRowMetadata::Windows {
                    window_mode: *window_mode,
                    windows: selected_windows,
                })
            }
            EndMotifRowMetadata::Groups(groups) => {
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
                Ok(EndMotifRowMetadata::Groups(selected_groups))
            }
        }
    }
}

/// Builder for selecting rows and motifs from an `EndsOutput`.
///
/// The builder starts with all rows and all motifs selected. Set at most one
/// selector per axis; for example, use `motifs()` or `motifs_by_label()`, not
/// both. Conflicting selector calls are reported by `read()` together with
/// bounds, row-mode, and dense/sparse selection errors.
#[derive(Debug, Clone)]
pub struct EndsSelector<'a> {
    output: &'a EndsOutput,
    rows: EndMotifRowSelector,
    motifs: MotifAxisSelector,
    selection_error: Option<String>,
}

impl<'a> EndsSelector<'a> {
    /// Create a selector that initially includes all rows and motifs.
    fn new(output: &'a EndsOutput) -> Self {
        Self {
            output,
            rows: EndMotifRowSelector::All,
            motifs: MotifAxisSelector::All,
            selection_error: None,
        }
    }

    /// Select generic output rows by zero-based row index.
    ///
    /// Parameters
    /// ----------
    /// - `row_indices`:
    ///   Source row indices in output order. The returned selection keeps
    ///   this order and rejects duplicates.
    pub fn rows(self, row_indices: &[usize]) -> Self {
        self.set_rows(EndMotifRowSelector::Rows(row_indices.to_vec()), "rows")
    }

    /// Select window rows by zero-based window row index.
    ///
    /// `read()` returns an error if the loaded output is not windowed.
    ///
    /// Parameters
    /// ----------
    /// - `window_indices`:
    ///   Window row indices in output order. The returned selection keeps
    ///   this order and rejects duplicates.
    pub fn windows(self, window_indices: &[usize]) -> Self {
        self.set_rows(
            EndMotifRowSelector::Windows(window_indices.to_vec()),
            "windows",
        )
    }

    /// Select grouped rows by zero-based group row index.
    ///
    /// `read()` returns an error if the loaded output is not grouped.
    ///
    /// Parameters
    /// ----------
    /// - `group_indices`:
    ///   Group row indices in output order. The returned selection keeps this
    ///   order and rejects duplicates.
    pub fn groups(self, group_indices: &[usize]) -> Self {
        self.set_rows(
            EndMotifRowSelector::Groups(group_indices.to_vec()),
            "groups",
        )
    }

    /// Select grouped rows by group name.
    ///
    /// `read()` returns an error if the loaded output is not grouped or any
    /// requested name is missing or duplicated.
    ///
    /// Parameters
    /// ----------
    /// - `group_names`:
    ///   Group labels from grouped output metadata. The returned selection
    ///   follows this order and rejects duplicates.
    pub fn groups_by_name<S: AsRef<str>>(self, group_names: &[S]) -> Self {
        self.set_rows(
            EndMotifRowSelector::GroupNames(
                group_names
                    .iter()
                    .map(|group_name| group_name.as_ref().to_string())
                    .collect(),
            ),
            "groups_by_name",
        )
    }

    /// Select motifs by zero-based motif index.
    ///
    /// Parameters
    /// ----------
    /// - `motif_indices`:
    ///   Motif axis indices in output order. The returned selection keeps
    ///   this order and rejects duplicates.
    pub fn motifs(self, motif_indices: &[usize]) -> Self {
        self.set_motifs(MotifAxisSelector::Indices(motif_indices.to_vec()), "motifs")
    }

    /// Select motifs by motif or motif-group label.
    ///
    /// Parameters
    /// ----------
    /// - `motif_labels`:
    ///   Motif labels or motif-group labels from the output motif axis. The
    ///   returned selection follows this order and rejects duplicates.
    pub fn motifs_by_label<S: AsRef<str>>(self, motif_labels: &[S]) -> Self {
        self.set_motifs(
            MotifAxisSelector::Labels(
                motif_labels
                    .iter()
                    .map(|motif_label| motif_label.as_ref().to_string())
                    .collect(),
            ),
            "motifs_by_label",
        )
    }

    /// Record a row-axis selector or remember the first row-axis conflict.
    fn set_rows(mut self, selector: EndMotifRowSelector, selector_name: &'static str) -> Self {
        if let Some(previous_selector_name) = self.rows.selector_name() {
            self.record_axis_conflict("row", previous_selector_name, selector_name);
        } else {
            self.rows = selector;
        }
        self
    }

    /// Record a motif-axis selector or remember the first motif-axis conflict.
    fn set_motifs(mut self, selector: MotifAxisSelector, selector_name: &'static str) -> Self {
        if let Some(previous_selector_name) = self.motifs.selector_name() {
            self.record_axis_conflict("motif", previous_selector_name, selector_name);
        } else {
            self.motifs = selector;
        }
        self
    }

    /// Store the first selector conflict for later reporting by `read()`.
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

    /// Return an error if the builder recorded a selector-axis conflict.
    fn ensure_no_selector_conflict(&self) -> Result<()> {
        if let Some(selection_error) = &self.selection_error {
            bail!("{selection_error}");
        }
        Ok(())
    }

    /// Read the selected counts while preserving the loaded storage mode.
    ///
    /// Sparse selections copy only stored COO entries whose source row and
    /// motif are both selected. Stored entries outside the selected axes are
    /// skipped. Missing in-bounds cells remain implicit zero counts.
    pub fn read(self) -> OutputLoaderResult<EndMotifCountSelection> {
        self.ensure_no_selector_conflict()?;
        let (motif_indices, motif_labels) = match self.motifs {
            MotifAxisSelector::All => (None, None),
            MotifAxisSelector::Indices(indices) => (Some(indices), None),
            MotifAxisSelector::Labels(labels) => (None, Some(labels)),
        };
        let motif_indices = motif_indices.as_deref();
        let motif_labels = motif_labels.as_deref();

        let selection = match self.rows {
            EndMotifRowSelector::All => {
                if let Some(motif_labels) = motif_labels {
                    self.output
                        .select_counts_by_motif_label(None, Some(motif_labels))
                } else {
                    self.output.select_counts(None, motif_indices)
                }
            }
            EndMotifRowSelector::Rows(indices) => {
                if let Some(motif_labels) = motif_labels {
                    self.output
                        .select_counts_by_motif_label(Some(indices.as_slice()), Some(motif_labels))
                } else {
                    self.output
                        .select_counts(Some(indices.as_slice()), motif_indices)
                }
            }
            EndMotifRowSelector::Windows(indices) => {
                let resolved_motif_indices = resolve_motif_labels(self.output, motif_labels)?;
                let motif_indices = motif_indices.or(resolved_motif_indices.as_deref());
                self.output
                    .select_window_counts(Some(indices.as_slice()), motif_indices)
            }
            EndMotifRowSelector::Groups(indices) => {
                let resolved_motif_indices = resolve_motif_labels(self.output, motif_labels)?;
                let motif_indices = motif_indices.or(resolved_motif_indices.as_deref());
                self.output
                    .select_group_counts(Some(indices.as_slice()), motif_indices)
            }
            EndMotifRowSelector::GroupNames(names) => {
                let resolved_motif_indices = resolve_motif_labels(self.output, motif_labels)?;
                let motif_indices = motif_indices.or(resolved_motif_indices.as_deref());
                self.output
                    .select_group_counts_by_name(Some(names.as_slice()), motif_indices)
            }
        }?;
        Ok(selection)
    }
}

/// Row-axis selector state recorded by `EndsSelector`.
#[derive(Debug, Clone)]
enum EndMotifRowSelector {
    /// Select all rows.
    All,
    /// Select generic output rows by index.
    Rows(Vec<usize>),
    /// Select window rows by index and require windowed output.
    Windows(Vec<usize>),
    /// Select grouped rows by index and require grouped output.
    Groups(Vec<usize>),
    /// Select grouped rows by group label and require grouped output.
    GroupNames(Vec<String>),
}

impl EndMotifRowSelector {
    /// Return the public selector method name for conflict messages.
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

/// Motif-axis selector state recorded by `EndsSelector`.
#[derive(Debug, Clone)]
enum MotifAxisSelector {
    /// Select all motifs.
    All,
    /// Select motifs by motif-axis index.
    Indices(Vec<usize>),
    /// Select motifs by motif or motif-group label.
    Labels(Vec<String>),
}

impl MotifAxisSelector {
    /// Return the public selector method name for conflict messages.
    fn selector_name(&self) -> Option<&'static str> {
        match self {
            Self::All => None,
            Self::Indices(_) => Some("motifs"),
            Self::Labels(_) => Some("motifs_by_label"),
        }
    }
}

/// Compact metadata for loaded end-motif counts.
///
/// This is intended for quick inspection and logging. It collects the output
/// settings that otherwise live behind separate accessors on `EndsOutput`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EndMotifOutputMetadata {
    /// Whether counts are dense or sparse COO.
    pub storage_mode: EndMotifStorageMode,
    /// Whether rows are global, windows, or grouped-BED groups.
    pub row_mode: EndMotifRowMode,
    /// Whether motif-axis labels are motifs or motif groups.
    pub motif_axis_kind: EndMotifAxisKind,
    /// Number of count rows.
    pub row_count: usize,
    /// Number of motif-axis labels.
    pub motif_count: usize,
}

impl fmt::Display for EndMotifOutputMetadata {
    /// Render one-line output context for logs or interactive inspection.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "storage_mode={}, row_mode={}, motif_axis={}, row_count={}, motif_count={}",
            describe_end_motif_storage_mode(self.storage_mode),
            describe_end_motif_row_mode(self.row_mode),
            describe_end_motif_axis_kind(self.motif_axis_kind),
            self.row_count,
            self.motif_count
        )
    }
}

/// End-motif count storage mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EndMotifStorageMode {
    /// Full dense `counts[row, motif]` array.
    Dense,
    /// Sparse coordinate arrays under `sparse/`.
    SparseCoo,
}

/// Meaning of end-motif count rows.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EndMotifRowMode {
    /// One row covering the full selected input.
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
pub enum EndMotifAxisKind {
    /// Concrete end-motif labels.
    Motif,
    /// Motif-group labels from `--motifs-file`.
    MotifGroup,
}

fn describe_end_motif_storage_mode(storage_mode: EndMotifStorageMode) -> &'static str {
    match storage_mode {
        EndMotifStorageMode::Dense => "dense",
        EndMotifStorageMode::SparseCoo => "sparse COO",
    }
}

fn describe_end_motif_row_mode(row_mode: EndMotifRowMode) -> &'static str {
    match row_mode {
        EndMotifRowMode::Global => "global",
        EndMotifRowMode::SizeWindows => "size windows",
        EndMotifRowMode::BedWindows => "BED windows",
        EndMotifRowMode::Groups => "groups",
    }
}

fn describe_end_motif_axis_kind(axis_kind: EndMotifAxisKind) -> &'static str {
    match axis_kind {
        EndMotifAxisKind::Motif => "motifs",
        EndMotifAxisKind::MotifGroup => "motif groups",
    }
}

/// Row metadata for loaded end-motif counts.
#[derive(Debug, Clone, PartialEq)]
pub enum EndMotifRowMetadata {
    /// One global row.
    Global,
    /// One row per genomic window.
    Windows {
        /// Whether windows came from fixed-size generation or BED input.
        window_mode: EndMotifWindowMode,
        /// Window metadata in count-row order.
        windows: Vec<WindowRow>,
    },
    /// One row per grouped-BED group.
    Groups(Vec<EndMotifGroupRow>),
}

impl EndMotifRowMetadata {
    /// Return the row mode represented by this metadata.
    pub fn mode(&self) -> EndMotifRowMode {
        match self {
            Self::Global => EndMotifRowMode::Global,
            Self::Windows {
                window_mode: EndMotifWindowMode::Size,
                ..
            } => EndMotifRowMode::SizeWindows,
            Self::Windows {
                window_mode: EndMotifWindowMode::Bed,
                ..
            } => EndMotifRowMode::BedWindows,
            Self::Groups(_) => EndMotifRowMode::Groups,
        }
    }
}

/// Source of window rows in an end-motif output.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EndMotifWindowMode {
    /// Generated fixed-size windows.
    Size,
    /// User-provided BED windows.
    Bed,
}

/// Metadata for one grouped end-motif output row.
#[derive(Debug, Clone, PartialEq)]
pub struct EndMotifGroupRow {
    /// Zero-based row index in count-row order.
    pub index: usize,
    /// Public group name from the grouped BED input.
    pub name: String,
    /// Number of grouped BED windows contributing to the group.
    pub eligible_windows: u64,
    /// Length-weighted blacklist fraction across the group's windows.
    pub blacklisted_fraction: f64,
}

/// Native count storage read from an end-motif output.
#[derive(Debug, Clone, PartialEq)]
pub enum EndMotifCountsData {
    /// Dense count matrix with shape `(row, motif)`.
    Dense(DenseMatrix<f64>),
    /// Sparse COO count entries with dense shape metadata.
    Sparse(EndMotifSparseCounts),
}

impl EndMotifCountsData {
    /// Return the storage mode represented by this value.
    pub fn storage_mode(&self) -> EndMotifStorageMode {
        match self {
            Self::Dense(_) => EndMotifStorageMode::Dense,
            Self::Sparse(_) => EndMotifStorageMode::SparseCoo,
        }
    }

    /// Return the dense shape represented by this storage.
    #[cfg(feature = "cmd_ref_kmers")]
    pub(crate) fn shape(&self) -> (usize, usize) {
        match self {
            Self::Dense(counts) => counts.shape(),
            Self::Sparse(sparse) => sparse.shape(),
        }
    }
}

/// Selected end-motif counts with row and motif-axis metadata.
///
/// Selections preserve the storage mode of the loaded output. Use
/// `dense_counts()` or `sparse_counts()` when the storage mode is known. Use
/// `to_dense_matrix()` only when a dense matrix is explicitly wanted.
#[derive(Debug, Clone, PartialEq)]
pub struct EndMotifCountSelection {
    row_metadata: EndMotifRowMetadata,
    row_indices: Vec<usize>,
    motif_indices: Vec<usize>,
    motif_labels: Vec<String>,
    data: EndMotifCountsData,
}

impl EndMotifCountSelection {
    /// Replace only the count storage while keeping selected metadata.
    ///
    /// Reference correction first uses the ordinary selector, so row metadata,
    /// selected source row indices, motif indices, and motif labels all come
    /// from the same path as uncorrected counts. This helper swaps in corrected
    /// dense or sparse counts after checking that the replacement storage has
    /// the same `(row, motif)` shape as the selection metadata.
    #[cfg(feature = "cmd_ref_kmers")]
    pub(crate) fn with_data(self, data: EndMotifCountsData) -> Result<Self> {
        ensure!(
            data.shape() == self.shape(),
            "replacement end-motif count storage shape {:?} does not match selection shape {:?}",
            data.shape(),
            self.shape()
        );
        Ok(Self { data, ..self })
    }

    /// Return how selected counts are stored.
    pub fn storage_mode(&self) -> EndMotifStorageMode {
        self.data.storage_mode()
    }

    /// Return selected row metadata in selection order.
    pub fn row_metadata(&self) -> &EndMotifRowMetadata {
        &self.row_metadata
    }

    /// Return selected window metadata, or an error if the selection is not windowed.
    pub fn window_metadata(&self) -> OutputLoaderResult<&[WindowRow]> {
        match &self.row_metadata {
            EndMotifRowMetadata::Windows { windows, .. } => Ok(windows),
            _ => Err(OutputLoaderError::message(
                "end-motif count selection is not windowed",
            )),
        }
    }

    /// Return selected group metadata, or an error if the selection is not grouped.
    pub fn group_metadata(&self) -> OutputLoaderResult<&[EndMotifGroupRow]> {
        match &self.row_metadata {
            EndMotifRowMetadata::Groups(groups) => Ok(groups),
            _ => Err(OutputLoaderError::message(
                "end-motif count selection is not grouped",
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

    /// Return the selected matrix shape as `(rows, motifs)`.
    pub fn shape(&self) -> (usize, usize) {
        match &self.data {
            EndMotifCountsData::Dense(counts) => counts.shape(),
            EndMotifCountsData::Sparse(sparse) => sparse.shape(),
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

    /// Return one selected count value, if both selection indices are in bounds.
    ///
    /// Sparse selections return `0.0` for in-bounds entries that are not stored
    /// in the COO vectors.
    ///
    /// **NOTE**: Each sparse call does one binary search over sorted COO
    /// entries. For repeated arbitrary cell lookups, build a reusable lookup
    /// once with `sparse_counts()?.to_lookup_index()`. For block-style access,
    /// use `sparse_counts()?.entries()` or `to_dense_matrix()`.
    ///
    /// Parameters
    /// ----------
    /// - `row_index`:
    ///   Zero-based row index within the selected matrix.
    /// - `motif_index`:
    ///   Zero-based motif index within the selected matrix.
    pub fn count(&self, row_index: usize, motif_index: usize) -> Option<f64> {
        match &self.data {
            EndMotifCountsData::Dense(counts) => counts.get(row_index, motif_index).copied(),
            EndMotifCountsData::Sparse(sparse) => sparse.count(row_index, motif_index),
        }
    }

    /// Return dense selected counts, or an error if this selection is sparse.
    pub fn dense_counts(&self) -> OutputLoaderResult<&DenseMatrix<f64>> {
        match &self.data {
            EndMotifCountsData::Dense(counts) => Ok(counts),
            EndMotifCountsData::Sparse(_) => Err(OutputLoaderError::message(
                "end-motif count selection is sparse",
            )),
        }
    }

    /// Return sparse selected counts, or an error if this selection is dense.
    pub fn sparse_counts(&self) -> OutputLoaderResult<&EndMotifSparseCounts> {
        match &self.data {
            EndMotifCountsData::Sparse(sparse) => Ok(sparse),
            EndMotifCountsData::Dense(_) => Err(OutputLoaderError::message(
                "end-motif count selection is dense",
            )),
        }
    }

    /// Return selected count storage as either a dense matrix or sparse COO entries.
    ///
    /// The variant matches this selection's storage mode.
    pub fn data(&self) -> &EndMotifCountsData {
        &self.data
    }

    /// Return a dense matrix for the selection.
    ///
    /// Dense selections are cloned into the returned matrix. Sparse selections
    /// are explicitly densified into a row-major matrix.
    pub fn to_dense_matrix(&self) -> OutputLoaderResult<DenseMatrix<f64>> {
        match &self.data {
            EndMotifCountsData::Dense(counts) => Ok(counts.clone()),
            EndMotifCountsData::Sparse(sparse) => sparse.to_dense_matrix(),
        }
    }
}

/// Sparse COO end-motif count entries.
///
/// COO stores explicit coordinates and their counts. The returned slices are
/// index-matched: `counts()[entry_index]` belongs at
/// `(row_indices()[entry_index], motif_indices()[entry_index])`.
/// For entry `entry_index`, the row, motif, and count are:
///
/// ```text
/// row_indices()[entry_index], motif_indices()[entry_index], counts()[entry_index]
/// ```
///
/// Coordinates are sorted by `(row, motif)`. Missing in-bounds coordinates are
/// implicit zero counts. Use this representation when iterating stored entries
/// or keeping sparse memory use. Use `to_dense_matrix()` when downstream code
/// needs dense row and column operations. Use `to_lookup_index()` when downstream
/// code needs many arbitrary `count(row, motif)` lookups without densifying.
#[derive(Debug, Clone, PartialEq)]
pub struct EndMotifSparseCounts {
    row_count: usize,
    motif_count: usize,
    row_indices: Vec<usize>,
    motif_indices: Vec<usize>,
    counts: Vec<f64>,
}

impl EndMotifSparseCounts {
    /// Build sparse counts from corrected entries and a dense shape.
    ///
    /// This is an internal constructor for reference-correction results. The
    /// entries are produced from an already loaded sparse selection, so this is
    /// not user-input validation. The caller supplies the selected matrix shape
    /// and entries whose coordinates belong to that shape. This function sorts
    /// the entries to preserve the sparse COO ordering contract.
    #[cfg(feature = "cmd_ref_kmers")]
    pub(crate) fn from_entries(
        row_count: usize,
        motif_count: usize,
        mut entries: Vec<EndMotifSparseEntry>,
    ) -> Self {
        entries.sort_by_key(|entry| (entry.row_index, entry.motif_index));

        let mut row_indices = Vec::with_capacity(entries.len());
        let mut motif_indices = Vec::with_capacity(entries.len());
        let mut counts = Vec::with_capacity(entries.len());
        for entry in entries {
            row_indices.push(entry.row_index);
            motif_indices.push(entry.motif_index);
            counts.push(entry.count);
        }

        Self {
            row_count,
            motif_count,
            row_indices,
            motif_indices,
            counts,
        }
    }

    /// Return the dense shape represented by the sparse entries.
    pub fn shape(&self) -> (usize, usize) {
        (self.row_count, self.motif_count)
    }

    /// Return zero-based row indices for stored COO entries.
    ///
    /// `row_indices()[entry_index]` gives the row for `counts()[entry_index]`
    /// and `motif_indices()[entry_index]`.
    pub fn row_indices(&self) -> &[usize] {
        &self.row_indices
    }

    /// Return zero-based motif indices for stored COO entries.
    ///
    /// `motif_indices()[entry_index]` gives the motif for
    /// `counts()[entry_index]` and `row_indices()[entry_index]`.
    pub fn motif_indices(&self) -> &[usize] {
        &self.motif_indices
    }

    /// Return stored counts in COO entry order.
    ///
    /// `counts()[entry_index]` belongs at
    /// `(row_indices()[entry_index], motif_indices()[entry_index])`.
    /// Coordinates not present in these slices are implicit zero counts, as
    /// long as they are inside `shape()`.
    ///
    /// Prefer `entries()` when iterating stored entries in ordinary Rust code.
    /// Calling `count()` repeatedly performs one binary search per point lookup.
    /// Use `to_lookup_index()` for many arbitrary point lookups without dense
    /// matrix allocation. Use `to_dense_matrix()` when the downstream work needs
    /// dense row or column access.
    pub fn counts(&self) -> &[f64] {
        &self.counts
    }

    /// Return stored counts as `(row, motif, count)` entries.
    ///
    /// Use this when downstream code only needs stored entries. Missing
    /// in-bounds coordinates are not yielded, because they are implicit zero
    /// counts.
    pub fn entries(&self) -> impl ExactSizeIterator<Item = EndMotifSparseEntry> + '_ {
        self.row_indices
            .iter()
            .copied()
            .zip(self.motif_indices.iter().copied())
            .zip(self.counts.iter().copied())
            .map(|((row_index, motif_index), count)| EndMotifSparseEntry {
                row_index,
                motif_index,
                count,
            })
    }

    /// Return the number of stored COO entries.
    pub fn nnz(&self) -> usize {
        self.counts.len()
    }

    /// Return one count value, if both indices are in bounds.
    ///
    /// Missing in-bounds sparse coordinates are zero counts.
    ///
    /// **NOTE**: This method does one binary search per call over the sorted
    /// COO entries. For repeated arbitrary cell lookups, build
    /// `to_lookup_index()` once. Use `entries()` for stored-entry iteration, or
    /// `to_dense_matrix()` when dense access is needed.
    ///
    /// Parameters
    /// ----------
    /// - `row_index`:
    ///   Zero-based row index in the dense matrix represented by this sparse
    ///   object.
    /// - `motif_index`:
    ///   Zero-based motif index in the dense matrix represented by this
    ///   sparse object.
    pub fn count(&self, row_index: usize, motif_index: usize) -> Option<f64> {
        if row_index >= self.row_count || motif_index >= self.motif_count {
            return None;
        }
        // Sparse COO entries are validated as sorted by `(row, motif)` when
        // the store is loaded. Binary search avoids scanning every stored entry
        // for each lookup. The left and right bounds form a half-open range of
        // candidate COO entry indices, and the middle entry is the coordinate
        // used to decide which half can still contain the requested count.
        let mut left_entry_index = 0;
        let mut right_entry_index = self.counts.len();
        while left_entry_index < right_entry_index {
            let middle_entry_index = left_entry_index + (right_entry_index - left_entry_index) / 2;
            let stored_coordinate = (
                self.row_indices[middle_entry_index],
                self.motif_indices[middle_entry_index],
            );
            match stored_coordinate.cmp(&(row_index, motif_index)) {
                std::cmp::Ordering::Less => left_entry_index = middle_entry_index + 1,
                std::cmp::Ordering::Equal => return Some(self.counts[middle_entry_index]),
                std::cmp::Ordering::Greater => right_entry_index = middle_entry_index,
            }
        }
        Some(0.0)
    }

    /// Build an `FxHashMap`-backed lookup index for repeated point queries.
    ///
    /// This copies the stored entries into an `FxHashMap` keyed by
    /// `(row_index, motif_index)`. It uses more memory than sorted COO, but
    /// `EndMotifSparseCountLookup::count()` avoids the per-lookup binary search
    /// done by `EndMotifSparseCounts::count()`. Use this for repeated
    /// arbitrary cell lookups when a dense matrix would be too large.
    pub fn to_lookup_index(&self) -> EndMotifSparseCountLookup {
        let mut counts_by_coordinate =
            FxHashMap::with_capacity_and_hasher(self.counts.len(), Default::default());
        for entry in self.entries() {
            counts_by_coordinate.insert((entry.row_index, entry.motif_index), entry.count);
        }
        EndMotifSparseCountLookup {
            row_count: self.row_count,
            motif_count: self.motif_count,
            counts_by_coordinate,
        }
    }

    /// Reconstruct a dense count matrix from the stored COO entries.
    pub fn to_dense_matrix(&self) -> OutputLoaderResult<DenseMatrix<f64>> {
        let value_count = self
            .row_count
            .checked_mul(self.motif_count)
            .context("end-motif dense shape overflow")?;
        let mut values = vec![0.0; value_count];
        for ((&row_index, &motif_index), &count) in self
            .row_indices
            .iter()
            .zip(&self.motif_indices)
            .zip(&self.counts)
        {
            let value_index = row_index
                .checked_mul(self.motif_count)
                .and_then(|row_start| row_start.checked_add(motif_index))
                .context("end-motif sparse coordinate overflow")?;
            values[value_index] = count;
        }
        Ok(DenseMatrix::from_row_major(
            values,
            self.row_count,
            self.motif_count,
        )?)
    }

    /// Build a sparse selection from selected source row and motif indices.
    fn select_counts(
        &self,
        row_indices: &[usize],
        motif_indices: &[usize],
        row_label: &str,
    ) -> Result<Self> {
        let row_index_map = build_selection_index_map(row_indices, self.row_count, row_label)?;
        let motif_index_map = build_selection_index_map(motif_indices, self.motif_count, "motif")?;

        // Keep only stored entries whose source row and motif are both
        // selected. Coordinates not stored in the source sparse output are
        // implicit zero counts and are not materialized in the selection. The
        // final sort restores COO order when selectors request rows or motifs
        // in a different order than the source file
        let mut selected_entries = Vec::new();
        for entry_index in 0..self.counts.len() {
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
                self.counts[entry_index],
            ));
        }
        selected_entries.sort_by_key(|&(row_index, motif_index, _)| (row_index, motif_index));

        let mut selected_row_indices = Vec::with_capacity(selected_entries.len());
        let mut selected_motif_indices = Vec::with_capacity(selected_entries.len());
        let mut selected_counts = Vec::with_capacity(selected_entries.len());
        for (row_index, motif_index, count) in selected_entries {
            selected_row_indices.push(row_index);
            selected_motif_indices.push(motif_index);
            selected_counts.push(count);
        }

        Ok(Self {
            row_count: row_indices.len(),
            motif_count: motif_indices.len(),
            row_indices: selected_row_indices,
            motif_indices: selected_motif_indices,
            counts: selected_counts,
        })
    }
}

/// One stored entry from sparse end-motif COO counts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EndMotifSparseEntry {
    /// Zero-based row index in the dense matrix represented by the sparse counts.
    pub row_index: usize,
    /// Zero-based motif index in the dense matrix represented by the sparse counts.
    pub motif_index: usize,
    /// Stored count for `(row_index, motif_index)`.
    pub count: f64,
}

/// `FxHashMap`-backed lookup index for sparse end-motif counts.
///
/// Build this with `EndMotifSparseCounts::to_lookup_index()` when downstream
/// code needs many arbitrary point lookups. This stores one `FxHashMap` entry
/// per stored sparse count, so it uses more memory than sorted COO but keeps
/// missing in-bounds coordinates as implicit zero counts.
#[derive(Debug, Clone, PartialEq)]
pub struct EndMotifSparseCountLookup {
    row_count: usize,
    motif_count: usize,
    counts_by_coordinate: FxHashMap<(usize, usize), f64>,
}

impl EndMotifSparseCountLookup {
    /// Return the dense shape represented by this lookup index.
    pub fn shape(&self) -> (usize, usize) {
        (self.row_count, self.motif_count)
    }

    /// Return one count value, if both indices are in bounds.
    ///
    /// Missing in-bounds coordinates are zero counts. This method uses a hash
    /// map lookup instead of binary search over the COO entries.
    ///
    /// Parameters
    /// ----------
    /// - `row_index`:
    ///   Zero-based row index in the represented dense matrix.
    /// - `motif_index`:
    ///   Zero-based motif index in the represented dense matrix.
    pub fn count(&self, row_index: usize, motif_index: usize) -> Option<f64> {
        if row_index >= self.row_count || motif_index >= self.motif_count {
            return None;
        }
        Some(
            self.counts_by_coordinate
                .get(&(row_index, motif_index))
                .copied()
                .unwrap_or(0.0),
        )
    }
}

/// Build a dense selection from selected source row and motif indices.
fn select_dense_counts(
    counts: &DenseMatrix<f64>,
    row_indices: &[usize],
    motif_indices: &[usize],
    row_label: &str,
) -> Result<DenseMatrix<f64>> {
    let mut selected_values =
        Vec::with_capacity(row_indices.len().saturating_mul(motif_indices.len()));
    let contiguous_columns = contiguous_index_span(motif_indices);
    for &row_index in row_indices {
        ensure!(
            row_index < counts.row_count(),
            "{row_label} index {row_index} is outside 0..{}",
            counts.row_count()
        );
        if let Some((start_column, end_column)) = contiguous_columns {
            let contiguous_row_values = counts
                .row_values_for_column_range(row_index, start_column, end_column)
                .with_context(|| {
                    format!(
                        "motif range {start_column}..{end_column} is outside 0..{}",
                        counts.column_count()
                    )
                })?;
            selected_values.extend_from_slice(contiguous_row_values);
        } else {
            let row_values = counts.row(row_index).with_context(|| {
                format!(
                    "{row_label} index {row_index} is outside 0..{}",
                    counts.row_count()
                )
            })?;
            for &motif_index in motif_indices {
                let count = row_values.get(motif_index).copied().with_context(|| {
                    format!(
                        "motif index {motif_index} is outside 0..{}",
                        counts.column_count()
                    )
                })?;
                selected_values.push(count);
            }
        }
    }
    DenseMatrix::from_row_major(selected_values, row_indices.len(), motif_indices.len())
}

/// Resolve optional motif indices to an explicit motif-axis selection.
fn resolve_motif_indices(
    motif_indices: Option<&[usize]>,
    motif_count: usize,
) -> Result<Vec<usize>> {
    if let Some(motif_indices) = motif_indices {
        return Ok(motif_indices.to_vec());
    }
    Ok((0..motif_count).collect())
}

/// Resolve motif labels to motif-axis indices.
fn resolve_motif_label_indices<S: AsRef<str>>(
    output: &EndsOutput,
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

/// Resolve optional motif labels to optional motif-axis indices.
fn resolve_motif_labels<S: AsRef<str>>(
    output: &EndsOutput,
    motif_labels: Option<&[S]>,
) -> Result<Option<Vec<usize>>> {
    match motif_labels {
        Some(labels) => Ok(Some(resolve_motif_label_indices(output, labels)?)),
        None => Ok(None),
    }
}

/// Resolve group names to row indices in grouped output metadata.
fn resolve_group_name_indices<S: AsRef<str>>(
    output: &EndsOutput,
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

/// Parser for a `cfdna ends` Zarr store.
///
/// The parser validates root metadata, reads row and motif axes, then reads the
/// count data in the storage mode declared by root metadata. Dense counts become
/// a `DenseMatrix`, while sparse counts stay in COO form.
struct EndsParser {
    path: PathBuf,
}

impl EndsParser {
    /// Create a parser for one end-motif Zarr store path.
    fn new(path: &Path) -> Self {
        Self {
            path: path.to_path_buf(),
        }
    }

    /// Load and validate the end-motif Zarr store into owned metadata and counts.
    fn load(&self) -> Result<EndsOutput> {
        validate_zarr_store_path(&self.path)?;
        let root_attributes = read_zarr_root_attributes(&self.path)
            .with_context(|| format!("read end-motif Zarr metadata {}", self.path.display()))?;
        let root_metadata = EndMotifRootMetadata::from_attributes(&root_attributes)?;
        let store = Arc::new(
            FilesystemStore::new(&self.path)
                .with_context(|| format!("open end-motif Zarr store {}", self.path.display()))?,
        );

        // Motif labels define the count matrix columns
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

        // Each count row represents one global output row, genomic window, or
        // grouped-BED group, depending on row_mode. The /row array is the
        // zero-based coordinate axis for those row records.
        let row = read_zarr_array1::<i32>(store.clone(), "/row")?;
        validate_zero_based_axis(&row, "row")?;
        ensure!(!row.is_empty(), "end-motif output has no rows");

        let row_metadata =
            read_row_metadata(&self.path, store.clone(), root_metadata.row_mode, row.len())?;
        let data = match root_metadata.storage_mode {
            EndMotifStorageMode::Dense => {
                EndMotifCountsData::Dense(read_dense_counts(store, row.len(), motif_labels.len())?)
            }
            EndMotifStorageMode::SparseCoo => EndMotifCountsData::Sparse(read_sparse_counts(
                &self.path,
                store,
                row.len(),
                motif_labels.len(),
            )?),
        };

        Ok(EndsOutput {
            row_metadata,
            motif_axis_kind: root_metadata.motif_axis_kind,
            motif_labels,
            motif_label_indices,
            data,
        })
    }
}

/// Required root-level metadata from an end-motif Zarr store.
///
/// The parsed values decide which row metadata arrays and count arrays are
/// valid for the store, so this is checked before reading row or count data.
struct EndMotifRootMetadata {
    storage_mode: EndMotifStorageMode,
    row_mode: EndMotifRowMode,
    motif_axis_kind: EndMotifAxisKind,
}

impl EndMotifRootMetadata {
    /// Parse and validate required root attributes from a Zarr metadata object.
    fn from_attributes(attributes: &Value) -> Result<Self> {
        ensure!(
            string_attr(attributes, "cfdnalab_schema")? == "end_motif_counts",
            "end-motif Zarr schema mismatch"
        );
        ensure!(
            u64_attr(attributes, "cfdnalab_schema_version")? == END_MOTIF_SCHEMA_VERSION,
            "end-motif Zarr schema version mismatch: expected {}",
            END_MOTIF_SCHEMA_VERSION
        );
        ensure!(
            string_attr(attributes, "count_units")? == "weighted_end_motif_count",
            "end-motif Zarr count_units must be weighted_end_motif_count"
        );

        let storage_mode = match string_attr(attributes, "storage_mode")? {
            "dense" => {
                ensure!(
                    string_attr(attributes, "primary_array")? == "counts",
                    "dense end-motif Zarr primary_array must be counts"
                );
                ensure!(
                    attributes.get("primary_group").is_some_and(Value::is_null),
                    "dense end-motif Zarr primary_group must be null"
                );
                EndMotifStorageMode::Dense
            }
            "sparse_coo" => {
                ensure!(
                    attributes.get("primary_array").is_some_and(Value::is_null),
                    "sparse end-motif Zarr primary_array must be null"
                );
                ensure!(
                    string_attr(attributes, "primary_group")? == "sparse",
                    "sparse end-motif Zarr primary_group must be sparse"
                );
                ensure!(
                    string_attr(attributes, "sparse_format")? == "coo",
                    "sparse end-motif Zarr sparse_format must be coo"
                );
                ensure!(
                    u64_attr(attributes, "sparse_indices_base")? == 0,
                    "sparse end-motif Zarr sparse_indices_base must be 0"
                );
                EndMotifStorageMode::SparseCoo
            }
            other => bail!("unsupported end-motif storage_mode '{other}'"),
        };
        let row_mode = match string_attr(attributes, "row_mode")? {
            "global" => EndMotifRowMode::Global,
            "size" => EndMotifRowMode::SizeWindows,
            "bed" => EndMotifRowMode::BedWindows,
            "grouped_bed" => EndMotifRowMode::Groups,
            other => bail!("unsupported end-motif row_mode '{other}'"),
        };
        let motif_axis_kind = match string_attr(attributes, "motif_axis_kind")? {
            "motif" => EndMotifAxisKind::Motif,
            "motif_group" => EndMotifAxisKind::MotifGroup,
            other => bail!("unsupported end-motif motif_axis_kind '{other}'"),
        };

        Ok(Self {
            storage_mode,
            row_mode,
            motif_axis_kind,
        })
    }
}

/// Validate that a path points to an end-motif Zarr store directory.
fn validate_zarr_store_path(path: &Path) -> Result<()> {
    ensure!(
        path.exists(),
        "end-motif Zarr store does not exist: {}",
        path.display()
    );
    ensure!(
        path.is_dir(),
        "end-motif Zarr store is not a directory: {}",
        path.display()
    );
    ensure!(
        path.file_name()
            .and_then(|name| name.to_str())
            .is_some_and(|name| name.ends_with(".zarr")),
        "end-motif Zarr store path must end with .zarr: {}",
        path.display()
    );
    ensure!(
        path.join("zarr.json").is_file(),
        "end-motif Zarr store is missing root zarr.json: {}",
        path.display()
    );
    Ok(())
}

/// Read motif-axis labels for motif or motif-group outputs.
fn read_motif_labels(
    root_path: &Path,
    store: Arc<FilesystemStore>,
    motif_axis_kind: EndMotifAxisKind,
    motif_count: usize,
) -> Result<Vec<String>> {
    match motif_axis_kind {
        EndMotifAxisKind::Motif => read_motif_ascii_labels(store, motif_count),
        EndMotifAxisKind::MotifGroup => {
            read_zarr_labels(root_path, "motif_index", "motif_group", motif_count)
        }
    }
}

/// Decode fixed-width motif labels from the `motif_ascii` byte matrix.
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
            let motif = String::from_utf8(motif_bytes.to_vec())
                .context("motif_ascii contains invalid UTF-8")?;
            validate_zarr_public_label(&motif, "motif")?;
            Ok(motif)
        })
        .collect()
}

/// Dispatch to the row metadata reader required by the root `row_mode`.
fn read_row_metadata(
    root_path: &Path,
    store: Arc<FilesystemStore>,
    row_mode: EndMotifRowMode,
    row_count: usize,
) -> Result<EndMotifRowMetadata> {
    match row_mode {
        EndMotifRowMode::Global => {
            let labels = read_zarr_labels(root_path, "row", "row_label", row_count)?;
            ensure!(
                labels == ["global"],
                "global end-motif row labels must be [\"global\"]"
            );
            Ok(EndMotifRowMetadata::Global)
        }
        EndMotifRowMode::SizeWindows | EndMotifRowMode::BedWindows => {
            read_window_row_metadata(root_path, store, row_mode, row_count)
        }
        EndMotifRowMode::Groups => read_group_row_metadata(root_path, store, row_count),
    }
}

/// Read genomic-window row metadata for size-window and BED-window outputs.
fn read_window_row_metadata(
    root_path: &Path,
    store: Arc<FilesystemStore>,
    row_mode: EndMotifRowMode,
    row_count: usize,
) -> Result<EndMotifRowMetadata> {
    // Read chromosome dictionary first because row_chromosome stores
    // chromosome-axis indices
    let chromosome = read_zarr_array1::<i32>(store.clone(), "/chromosome")?;
    validate_zero_based_axis(&chromosome, "chromosome")?;
    let chromosome_names =
        read_zarr_labels(root_path, "chromosome", "chromosome_name", chromosome.len())?;

    // Per-row arrays all use the /row axis length
    let row_chromosome = read_zarr_array1::<i32>(store.clone(), "/row_chromosome")?;
    let row_start_bp = read_zarr_array1::<i64>(store.clone(), "/row_start_bp")?;
    let row_end_bp = read_zarr_array1::<i64>(store.clone(), "/row_end_bp")?;
    let blacklisted_fraction = read_zarr_array1::<f64>(store, "/blacklisted_fraction")?;
    ensure_same_len(&row_chromosome, row_count, "row_chromosome")?;
    ensure_same_len(&row_start_bp, row_count, "row_start_bp")?;
    ensure_same_len(&row_end_bp, row_count, "row_end_bp")?;
    ensure_same_len(&blacklisted_fraction, row_count, "blacklisted_fraction")?;

    // Build WindowRow values in row-axis order so indices match count rows
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
            anyhow::anyhow!("end-motif row {row_index} has invalid interval: {error}")
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
        EndMotifRowMode::SizeWindows => EndMotifWindowMode::Size,
        EndMotifRowMode::BedWindows => EndMotifWindowMode::Bed,
        _ => bail!("internal end-motif loader row-mode mismatch"),
    };
    Ok(EndMotifRowMetadata::Windows {
        window_mode,
        windows,
    })
}

/// Read grouped-BED row metadata for group-mode outputs.
fn read_group_row_metadata(
    root_path: &Path,
    store: Arc<FilesystemStore>,
    row_count: usize,
) -> Result<EndMotifRowMetadata> {
    // Group is the row coordinate axis for grouped-BED output
    let group = read_zarr_array1::<i32>(store.clone(), "/group")?;
    validate_zero_based_axis(&group, "group")?;
    ensure_same_len(&group, row_count, "group")?;

    // Group labels live on the row axis, one group name per count row
    let group_names = read_zarr_labels(root_path, "group", "group_name", row_count)?;
    ensure_unique_labels(&group_names, "group_names")?;

    // Group-level arrays must have one value per count row
    let eligible_windows = read_zarr_array1::<i32>(store.clone(), "/eligible_windows")?;
    let blacklisted_fraction = read_zarr_array1::<f64>(store, "/blacklisted_fraction")?;
    ensure_same_len(&eligible_windows, row_count, "eligible_windows")?;
    ensure_same_len(&blacklisted_fraction, row_count, "blacklisted_fraction")?;

    // Build group rows in row-axis order so indices match count rows
    let groups = (0..row_count)
        .map(|row_index| {
            let eligible_windows = u64_from_i32(eligible_windows[row_index], "eligible_windows")?;
            let blacklisted_fraction = blacklisted_fraction[row_index];
            ensure_valid_fraction(blacklisted_fraction, "blacklisted_fraction")?;
            Ok(EndMotifGroupRow {
                index: row_index,
                name: group_names[row_index].clone(),
                eligible_windows,
                blacklisted_fraction,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(EndMotifRowMetadata::Groups(groups))
}

/// Read and validate the dense count matrix.
fn read_dense_counts(
    store: Arc<FilesystemStore>,
    row_count: usize,
    motif_count: usize,
) -> Result<DenseMatrix<f64>> {
    let (values, shape) = read_zarr_array_values::<f64>(store, "/counts")?;
    ensure!(
        shape == [row_count, motif_count],
        "dense end-motif counts shape {:?} did not match expected [{}, {}]",
        shape,
        row_count,
        motif_count
    );
    ensure_non_negative_finite_counts(&values, "dense end-motif counts")?;
    DenseMatrix::from_row_major(values, row_count, motif_count)
}

/// Read sparse COO counts and validate coordinate arrays against metadata shape.
fn read_sparse_counts(
    root_path: &Path,
    store: Arc<FilesystemStore>,
    row_count: usize,
    motif_count: usize,
) -> Result<EndMotifSparseCounts> {
    // The sparse group declares which dense axes the COO coordinates address
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

    // The sparse shape must match the row and motif axes already read from metadata
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

    // COO arrays store one row coordinate, motif coordinate, and count per stored entry
    let row = read_zarr_array1::<i32>(store.clone(), "/sparse/row")?;
    let motif = read_zarr_array1::<i32>(store.clone(), "/sparse/motif")?;
    let counts = read_zarr_array1::<f64>(store, "/sparse/count")?;
    ensure_same_len(&motif, row.len(), "sparse/motif")?;
    ensure_same_len(&counts, row.len(), "sparse/count")?;
    ensure_non_negative_finite_counts(&counts, "sparse end-motif counts")?;

    // Validate coordinates before storing them as usize values used for lookup
    // and densification
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
            // Sorted unique coordinates make binary-search lookup possible
            ensure!(
                previous_coordinate < coordinate,
                "sparse end-motif COO entries must be sorted and unique"
            );
        }
        previous_coordinate = Some(coordinate);
        row_indices.push(row_index);
        motif_indices.push(motif_index);
    }

    Ok(EndMotifSparseCounts {
        row_count,
        motif_count,
        row_indices,
        motif_indices,
        counts,
    })
}

/// Validate count values that should come from non-negative weighted counts.
fn ensure_non_negative_finite_counts(counts: &[f64], label: &str) -> Result<()> {
    for (count_index, &count) in counts.iter().enumerate() {
        ensure!(
            count.is_finite() && count >= 0.0,
            "{label} contain value outside finite and non-negative range at index {count_index}: {count}"
        );
    }
    Ok(())
}

/// Read all values and shape metadata from one Zarr array.
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

/// Read a rank-one Zarr array.
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

/// Read string labels from a Zarr array's attributes.
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

/// Read the `attributes` object from a Zarr array metadata file.
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

/// Build the local `zarr.json` path for an array inside a Zarr store.
fn zarr_metadata_path(root_path: &Path, array_path: &str) -> PathBuf {
    let mut path = root_path.to_path_buf();
    for component in array_path.trim_start_matches('/').split('/') {
        path.push(component);
    }
    path.join("zarr.json")
}

/// Validate that an integer coordinate array is exactly `0..len`.
fn validate_zero_based_axis(values: &[i32], axis_name: &str) -> Result<()> {
    for (index, &value) in values.iter().enumerate() {
        ensure!(
            value == i32::try_from(index).context("axis index exceeds i32")?,
            "{axis_name} must be a zero-based coordinate axis"
        );
    }
    Ok(())
}

/// Ensure a metadata array has the expected number of entries.
fn ensure_same_len<T>(values: &[T], expected_len: usize, array_name: &str) -> Result<()> {
    ensure!(
        values.len() == expected_len,
        "{array_name} has {} entries, expected {expected_len}",
        values.len()
    );
    Ok(())
}

/// Ensure a value is a finite fraction in the inclusive range `[0, 1]`.
fn ensure_valid_fraction(value: f64, field_name: &str) -> Result<()> {
    ensure!(
        value.is_finite() && (0.0..=1.0).contains(&value),
        "{field_name} must be a finite fraction in [0, 1], got {value}"
    );
    Ok(())
}

/// Read a required string attribute from a Zarr metadata object.
fn string_attr<'a>(attributes: &'a Value, name: &str) -> Result<&'a str> {
    attributes
        .get(name)
        .and_then(Value::as_str)
        .with_context(|| format!("end-motif Zarr metadata is missing string attribute {name}"))
}

/// Read a required unsigned integer attribute from a Zarr metadata object.
fn u64_attr(attributes: &Value, name: &str) -> Result<u64> {
    attributes
        .get(name)
        .and_then(Value::as_u64)
        .with_context(|| format!("end-motif Zarr metadata is missing integer attribute {name}"))
}

/// Convert a non-negative `i32` field to `usize`.
fn usize_from_i32(value: i32, field_name: &str) -> Result<usize> {
    ensure!(value >= 0, "{field_name} must be non-negative, got {value}");
    Ok(usize::try_from(value).expect("non-negative i32 always fits usize"))
}

/// Convert a non-negative `i32` field to `u64`.
fn u64_from_i32(value: i32, field_name: &str) -> Result<u64> {
    ensure!(value >= 0, "{field_name} must be non-negative, got {value}");
    Ok(u64::try_from(value).expect("non-negative i32 always fits u64"))
}

/// Convert a non-negative `i64` field to `u64`.
fn u64_from_i64(value: i64, field_name: &str) -> Result<u64> {
    ensure!(value >= 0, "{field_name} must be non-negative, got {value}");
    Ok(u64::try_from(value).expect("non-negative i64 always fits u64"))
}
