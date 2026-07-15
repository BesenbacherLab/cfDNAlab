//! Loader for `cfdna lengths` output tables.
//!
//! A loaded lengths output is a count table. Columns are fragment length bins,
//! and rows are either one global row, genomic windows, or grouped-BED groups.
//! The loader keeps the count matrix, the length-bin axis, and the row metadata
//! together so code can work with the output without knowing the TSV header
//! layout.
//!
//! Type overview:
//!
//! ```text
//! load_lengths_output(path)
//!     -> LengthsOutput
//!         row_metadata: LengthRowMetadata
//!             Global
//!             Windows(Vec<WindowRow>)
//!             Groups(Vec<LengthGroupRow>)
//!         length_bins: Vec<LengthBin>
//!         counts: DenseMatrix<f64>
//!
//! LengthsOutput::select()
//!     -> LengthsSelector
//!         -> read()
//!             -> LengthCountSelection
//!                 row_metadata: LengthRowMetadata
//!                 row_indices: Vec<usize>
//!                 length_bins: Vec<LengthBin>
//!                 counts: DenseMatrix<f64>
//! ```
//!
//! The main object is `LengthsOutput`. Use it directly for the common count
//! operations:
//!
//! - `length_bins()` returns the fragment length bins in output-column order.
//!
//! - `counts()` returns the full dense matrix with shape
//!   `(number_of_rows, number_of_length_bins)`.
//!
//! - `count(row_index, length_bin_index)` returns a single value using the same
//!   zero-based indices stored in row metadata and length bins.
//!
//! - `length_bin_for_length()` and `length_bins_overlapping_range()` resolve
//!   biological length queries to output columns.
//!
//! Row metadata is separate from count values. `row_metadata()` tells you what
//! each matrix row represents. `window_metadata()` and `group_metadata()` are
//! row-mode-specific accessors for code that expects a windowed or grouped
//! output, and they return an error if the loaded file has a different row mode.
//!
//! `select()` returns a selector builder for extracting an owned
//! `LengthCountSelection`. Rows can be addressed as generic output rows,
//! windows, group indices, or group names. Fragment length bins can be addressed
//! by output-column index or by overlap with a half-open fragment length range.
//! Omitted axes select everything on that axis. The selected matrix preserves
//! requested order, rejects duplicate selectors, and copies whole row spans when
//! the selected length bins are contiguous.

use crate::{
    interval::Interval,
    output_loaders::{
        OutputLoaderError, OutputLoaderResult,
        common::{
            DenseMatrix, LengthBin, WindowRow, contiguous_index_span, ensure_unique_indices,
            ensure_unique_labels, resolve_row_indices,
        },
    },
    shared::{
        constants::{MAX_SUPPORTED_FRAGMENT_LENGTH, MIN_ACGT_BASES_FOR_GC_FRACTION},
        io::{open_text_reader, open_text_reader_in_background},
    },
};
use anyhow::{Context, Result, bail, ensure};
use std::{
    fmt,
    io::BufRead,
    path::{Path, PathBuf},
};

/// Load a `cfdna lengths` count table from a plain, gzip, or zstd-compressed TSV.
///
/// The loader reads the full count table eagerly. Length outputs are wide TSVs
/// where count columns are fragment length bins and rows are either one global
/// row, genomic windows, or grouped-BED groups. Public indices are zero-based,
/// and length bins use checked half-open interval semantics.
///
/// Most workflows can use the count methods directly. Use `row_metadata()`,
/// `window_metadata()`, or `group_metadata()` when row coordinates or group
/// labels are needed.
/// `select().read()` returns a `LengthCountSelection`, which is an owned count
/// matrix plus the selected source row indices and selected length-bin metadata.
///
/// Parameters
/// ----------
/// - `path`:
///   Path to a `cfdna lengths` TSV output. Plain text, gzip, and zstd
///   compressed files are supported.
///
/// Returns
/// -------
/// - `LengthsOutput`:
///   Loaded row metadata, fragment length bins, and the full count matrix.
///
/// ```no_run
/// use cfdnalab::{
///     interval::Interval,
///     output_loaders::load_lengths_output,
/// };
///
/// let lengths = load_lengths_output("sample.length_counts.tsv.zst")?;
///
/// let count_30bp = lengths
///     .length_bin_for_length(30)
///     .and_then(|bin_index| lengths.count(0, bin_index));
/// println!("first-row 30 bp count: {count_30bp:?}");
///
/// let selected = lengths
///     .select()
///     .groups_by_name(&["promoter", "enhancer"])
///     .length_range(Interval::new(100, 151)?)
///     .read()?;
///
/// // Build a downstream summary without parsing the TSV header yourself.
/// let selected_length_total_by_group = selected
///     .group_metadata()?
///     .iter()
///     .zip(selected.counts().rows())
///     .map(|(group, selected_length_counts)| {
///         let count_sum = selected_length_counts.iter().copied().sum::<f64>();
///         (group.name.as_str(), count_sum)
///     })
///     .collect::<Vec<_>>();
///
/// for (group_name, count_sum) in selected_length_total_by_group {
///     println!("{group_name}: {count_sum}");
/// }
///
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn load_lengths_output(path: impl AsRef<Path>) -> OutputLoaderResult<LengthsOutput> {
    let read_in_background =
        std::thread::available_parallelism().is_ok_and(|thread_count| thread_count.get() > 1);
    LengthsParser::new(path.as_ref(), read_in_background)
        .load()
        .map_err(Into::into)
}

#[cfg(all(test, feature = "testing"))]
include!("lengths_background_reading_benchmark.rs");

/// Row aggregation mode detected from the length-count table schema.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LengthOutputMode {
    /// One global row with only count columns.
    Global,
    /// One row per genomic window.
    Windows,
    /// One row per grouped-BED group.
    Groups,
}

impl LengthOutputMode {
    /// Return a short label used in error messages.
    fn description(self) -> &'static str {
        match self {
            Self::Global => "global",
            Self::Windows => "windowed",
            Self::Groups => "grouped",
        }
    }
}

/// Loaded output from `cfdna lengths`.
#[derive(Debug, Clone, PartialEq)]
pub struct LengthsOutput {
    row_metadata: LengthRowMetadata,
    length_bins: Vec<LengthBin>,
    counts: DenseMatrix<f64>,
}

impl LengthsOutput {
    /// Return the detected row aggregation mode.
    pub fn row_mode(&self) -> LengthOutputMode {
        self.row_metadata.mode()
    }

    /// Return row metadata describing what each count-matrix row represents.
    pub fn row_metadata(&self) -> &LengthRowMetadata {
        &self.row_metadata
    }

    /// Return a compact description of the loaded lengths output.
    ///
    /// This combines row mode, row count, fragment length-bin count, and the
    /// covered fragment length range in one value for logging or quick checks.
    pub fn output_metadata(&self) -> LengthOutputMetadata {
        LengthOutputMetadata {
            row_mode: self.row_mode(),
            row_count: self.row_count(),
            length_bin_count: self.length_bin_count(),
            min_fragment_length: self.length_bins.first().map(|bin| bin.start()),
            max_fragment_length_exclusive: self.length_bins.last().map(|bin| bin.end()),
        }
    }

    /// Return window metadata, or an error if this is not a windowed output.
    pub fn window_metadata(&self) -> OutputLoaderResult<&[WindowRow]> {
        match &self.row_metadata {
            LengthRowMetadata::Windows(windows) => Ok(windows),
            _ => Err(OutputLoaderError::message("lengths output is not windowed")),
        }
    }

    /// Return group metadata, or an error if this is not a grouped output.
    pub fn group_metadata(&self) -> OutputLoaderResult<&[LengthGroupRow]> {
        match &self.row_metadata {
            LengthRowMetadata::Groups(groups) => Ok(groups),
            _ => Err(OutputLoaderError::message("lengths output is not grouped")),
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
    pub fn group(&self, row_index: usize) -> OutputLoaderResult<Option<&LengthGroupRow>> {
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
            .with_context(|| format!("lengths output has no group named '{group_name}'"))?)
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

    /// Return fragment length bins in output-column order.
    pub fn length_bins(&self) -> &[LengthBin] {
        &self.length_bins
    }

    /// Return the number of fragment length bins.
    pub fn length_bin_count(&self) -> usize {
        self.length_bins.len()
    }

    /// Return the number of output rows.
    pub fn row_count(&self) -> usize {
        self.counts().shape().0
    }

    /// Return the count matrix.
    pub fn counts(&self) -> &DenseMatrix<f64> {
        &self.counts
    }

    /// Return all count values in row-major order.
    pub fn counts_row_major(&self) -> &[f64] {
        self.counts().values_row_major()
    }

    /// Return one count value, if both indices are in bounds.
    ///
    /// Parameters
    /// ----------
    /// - `row_index`:
    ///   Zero-based row index in the count matrix.
    /// - `length_bin_index`:
    ///   Zero-based fragment length-bin column index.
    pub fn count(&self, row_index: usize, length_bin_index: usize) -> Option<f64> {
        self.counts().get(row_index, length_bin_index).copied()
    }

    /// Return the index of the bin containing `fragment_length_bp`.
    ///
    /// Parameters
    /// ----------
    /// - `fragment_length_bp`:
    ///   Fragment length in bp to locate on the length-bin axis.
    pub fn length_bin_for_length(&self, fragment_length_bp: u32) -> Option<usize> {
        self.length_bins()
            .iter()
            .find(|bin| bin.start() <= fragment_length_bp && fragment_length_bp < bin.end())
            .map(|bin| bin.idx())
    }

    /// Return length bins overlapping a half-open fragment length range.
    ///
    /// A bin is selected when `[bin_start, bin_end)` intersects
    /// `[range_start, range_end)`. Bins that only touch the query boundary are
    /// not selected.
    ///
    /// Parameters
    /// ----------
    /// - `range`:
    ///   Half-open fragment length interval `[start, end)` in bp.
    pub fn length_bins_overlapping_range(
        &self,
        range: Interval<u32>,
    ) -> OutputLoaderResult<Vec<LengthBin>> {
        let selected = self
            .length_bins()
            .iter()
            .copied()
            .filter(|bin| bin.interval.intersects(range))
            .collect::<Vec<_>>();
        if selected.is_empty() {
            return Err(OutputLoaderError::message(format!(
                "length range [{}, {}) bp does not overlap any length bins",
                range.start(),
                range.end()
            )));
        }
        Ok(selected)
    }

    /// Start a count selection.
    ///
    /// A new selector initially selects all rows and all fragment length bins.
    /// Add row and length-bin constraints before calling `read()`.
    pub fn select(&self) -> LengthsSelector<'_> {
        LengthsSelector::new(self)
    }

    /// Return an owned count matrix for selected rows and length bins.
    ///
    /// Passing `None` for `row_indices` selects all rows. Passing `None` for
    /// both length-bin selectors selects all length bins. Use either
    /// `length_bin_indices` or `length_range`, not both.
    pub(crate) fn select_counts(
        &self,
        row_indices: Option<&[usize]>,
        length_bin_indices: Option<&[usize]>,
        length_range: Option<Interval<u32>>,
    ) -> Result<LengthCountSelection> {
        self.select_counts_with_label(row_indices, length_bin_indices, length_range, "row")
    }

    /// Return counts for selected window row indices and length bins.
    ///
    /// Passing `None` for `window_indices` selects all windows.
    pub(crate) fn select_window_counts(
        &self,
        window_indices: Option<&[usize]>,
        length_bin_indices: Option<&[usize]>,
        length_range: Option<Interval<u32>>,
    ) -> Result<LengthCountSelection> {
        self.ensure_row_mode(LengthOutputMode::Windows)?;
        self.select_counts_with_label(window_indices, length_bin_indices, length_range, "window")
    }

    /// Return counts for selected group row indices and length bins.
    ///
    /// Passing `None` for `group_indices` selects all groups.
    pub(crate) fn select_group_counts(
        &self,
        group_indices: Option<&[usize]>,
        length_bin_indices: Option<&[usize]>,
        length_range: Option<Interval<u32>>,
    ) -> Result<LengthCountSelection> {
        self.ensure_row_mode(LengthOutputMode::Groups)?;
        self.select_counts_with_label(group_indices, length_bin_indices, length_range, "group")
    }

    /// Return counts for selected group names and length bins.
    ///
    /// Passing `None` for `group_names` selects all groups.
    pub(crate) fn select_group_counts_by_name<S: AsRef<str>>(
        &self,
        group_names: Option<&[S]>,
        length_bin_indices: Option<&[usize]>,
        length_range: Option<Interval<u32>>,
    ) -> Result<LengthCountSelection> {
        let group_indices = match group_names {
            Some(group_names) => Some(self.resolve_group_name_indices(group_names)?),
            None => {
                self.ensure_row_mode(LengthOutputMode::Groups)?;
                None
            }
        };
        self.select_counts_with_label(
            group_indices.as_deref(),
            length_bin_indices,
            length_range,
            "group",
        )
    }

    /// Select rows and length bins after resolving row labels for error messages.
    fn select_counts_with_label(
        &self,
        row_indices: Option<&[usize]>,
        length_bin_indices: Option<&[usize]>,
        length_range: Option<Interval<u32>>,
        row_label: &str,
    ) -> Result<LengthCountSelection> {
        ensure!(
            row_indices.is_none() || self.row_mode() != LengthOutputMode::Global,
            "global lengths output has no selectable row axis"
        );
        // Resolve optional selectors first so the copy path works with
        // concrete matrix row and column indices
        let row_indices = resolve_row_indices(row_indices, self.row_count(), row_label)?;
        let length_bin_indices =
            self.resolve_length_bin_indices(length_bin_indices, length_range)?;

        // Reject repeated selectors rather than returning duplicated rows or
        // columns with duplicated metadata entries
        ensure_unique_indices(&row_indices, row_label)?;
        ensure_unique_indices(&length_bin_indices, "length bin")?;

        // Carry selected row and length-axis metadata alongside the copied
        // matrix values so the returned selection can be interpreted on its own
        let selected_row_metadata = self.selected_row_metadata(&row_indices, row_label)?;
        let selected_length_bins = length_bin_indices
            .iter()
            .map(|&length_bin_index| {
                self.length_bins
                    .get(length_bin_index)
                    .copied()
                    .with_context(|| {
                        format!(
                            "length bin index {length_bin_index} is outside 0..{}",
                            self.length_bins.len()
                        )
                    })
            })
            .collect::<Result<Vec<_>>>()?;

        // Build a new row-major matrix in requested row order. When selected
        // columns form one contiguous range, copy that range directly from
        // matrix storage
        let mut selected_values =
            Vec::with_capacity(row_indices.len().saturating_mul(length_bin_indices.len()));
        let contiguous_columns = contiguous_index_span(&length_bin_indices);
        for &row_index in &row_indices {
            ensure!(
                row_index < self.row_count(),
                "{row_label} index {row_index} is outside 0..{}",
                self.row_count()
            );
            if let Some((start_column, end_column)) = contiguous_columns {
                let contiguous_row_values = self
                    .counts
                    .row_values_for_column_range(row_index, start_column, end_column)
                    .with_context(|| {
                        format!(
                            "length bin range {start_column}..{end_column} is outside 0..{}",
                            self.length_bins.len()
                        )
                    })?;
                selected_values.extend_from_slice(contiguous_row_values);
            } else {
                let row_values = self.counts.row(row_index).with_context(|| {
                    format!(
                        "{row_label} index {row_index} is outside 0..{}",
                        self.row_count()
                    )
                })?;
                for &length_bin_index in &length_bin_indices {
                    let count = row_values.get(length_bin_index).copied().with_context(|| {
                        format!(
                            "length bin index {length_bin_index} is outside 0..{}",
                            self.length_bins.len()
                        )
                    })?;
                    selected_values.push(count);
                }
            }
        }

        // Re-wrap the copied values with their selection shape so callers can
        // use the same DenseMatrix API as on the full output
        let counts = DenseMatrix::from_row_major(
            selected_values,
            row_indices.len(),
            length_bin_indices.len(),
        )?;
        Ok(LengthCountSelection {
            row_metadata: selected_row_metadata,
            row_indices,
            length_bins: selected_length_bins,
            counts,
        })
    }

    /// Copy row metadata for the selected source row indices.
    fn selected_row_metadata(
        &self,
        row_indices: &[usize],
        row_label: &str,
    ) -> Result<LengthRowMetadata> {
        match &self.row_metadata {
            LengthRowMetadata::Global => Ok(LengthRowMetadata::Global),
            LengthRowMetadata::Windows(windows) => {
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
                Ok(LengthRowMetadata::Windows(selected_windows))
            }
            LengthRowMetadata::Groups(groups) => {
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
                Ok(LengthRowMetadata::Groups(selected_groups))
            }
        }
    }

    /// Resolve optional length-bin selectors to concrete source column indices.
    fn resolve_length_bin_indices(
        &self,
        length_bin_indices: Option<&[usize]>,
        length_range: Option<Interval<u32>>,
    ) -> Result<Vec<usize>> {
        ensure!(
            length_bin_indices.is_none() || length_range.is_none(),
            "use either length_bin_indices or length_range, not both"
        );
        if let Some(length_bin_indices) = length_bin_indices {
            return Ok(length_bin_indices.to_vec());
        }
        if let Some(length_range) = length_range {
            return Ok(self
                .length_bins_overlapping_range(length_range)?
                .iter()
                .map(|length_bin| length_bin.idx())
                .collect());
        }
        Ok((0..self.length_bins.len()).collect())
    }

    /// Resolve group labels to source row indices.
    fn resolve_group_name_indices<S: AsRef<str>>(&self, group_names: &[S]) -> Result<Vec<usize>> {
        ensure_unique_labels(group_names, "group_names")?;
        group_names
            .iter()
            .map(|group_name| {
                self.group_index(group_name.as_ref())
                    .map_err(anyhow::Error::from)
            })
            .collect()
    }

    /// Return an error when the loaded row mode does not match the typed selector.
    fn ensure_row_mode(&self, expected: LengthOutputMode) -> Result<()> {
        ensure!(
            self.row_mode() == expected,
            "lengths output is not {}",
            expected.description()
        );
        Ok(())
    }
}

/// Builder for selecting rows and fragment length bins from a `LengthsOutput`.
///
/// The builder starts with all rows and all length bins selected. Set at most
/// one selector per axis. For example, use `length_bins()` or `length_range()`,
/// not both. Conflicting selector calls are reported by `read()` together with
/// bounds validation and data-copy errors.
#[derive(Debug, Clone)]
pub struct LengthsSelector<'a> {
    output: &'a LengthsOutput,
    rows: LengthRowSelector,
    lengths: LengthAxisSelector,
    selection_error: Option<String>,
}

impl<'a> LengthsSelector<'a> {
    /// Start a selector with all rows and length bins selected.
    fn new(output: &'a LengthsOutput) -> Self {
        Self {
            output,
            rows: LengthRowSelector::All,
            lengths: LengthAxisSelector::All,
            selection_error: None,
        }
    }

    /// Select generic output rows by zero-based row index.
    ///
    /// Parameters
    /// ----------
    /// - `row_indices`:
    ///   Source row indices in output-file order. The returned selection keeps
    ///   this order and rejects duplicates.
    pub fn rows(self, row_indices: &[usize]) -> Self {
        self.set_rows(LengthRowSelector::Rows(row_indices.to_vec()), "rows")
    }

    /// Select window rows by zero-based window row index.
    ///
    /// `read()` returns an error if the loaded output is not windowed.
    ///
    /// Parameters
    /// ----------
    /// - `window_indices`:
    ///   Window row indices in output-file order. The returned selection keeps
    ///   this order and rejects duplicates.
    pub fn windows(self, window_indices: &[usize]) -> Self {
        self.set_rows(
            LengthRowSelector::Windows(window_indices.to_vec()),
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
    ///   Group row indices in output-file order. The returned selection keeps
    ///   this order and rejects duplicates.
    pub fn groups(self, group_indices: &[usize]) -> Self {
        self.set_rows(LengthRowSelector::Groups(group_indices.to_vec()), "groups")
    }

    /// Select grouped rows by group name.
    ///
    /// `read()` returns an error if the loaded output is not grouped or any
    /// requested name is missing or duplicated.
    ///
    /// Parameters
    /// ----------
    /// - `group_names`:
    ///   Group labels from the grouped output metadata. The returned selection
    ///   follows this order and rejects duplicates.
    pub fn groups_by_name<S: AsRef<str>>(self, group_names: &[S]) -> Self {
        self.set_rows(
            LengthRowSelector::GroupNames(
                group_names
                    .iter()
                    .map(|group_name| group_name.as_ref().to_string())
                    .collect(),
            ),
            "groups_by_name",
        )
    }

    /// Select fragment length bins by zero-based length-bin index.
    ///
    /// Parameters
    /// ----------
    /// - `length_bin_indices`:
    ///   Length-bin column indices in output order. The returned selection
    ///   keeps this order and rejects duplicates.
    pub fn length_bins(self, length_bin_indices: &[usize]) -> Self {
        self.set_lengths(
            LengthAxisSelector::Indices(length_bin_indices.to_vec()),
            "length_bins",
        )
    }

    /// Select fragment length bins overlapping a half-open length range.
    ///
    /// Parameters
    /// ----------
    /// - `range`:
    ///   Half-open fragment length interval `[start, end)` in bp. Every output
    ///   length bin that intersects the interval is selected.
    pub fn length_range(self, range: Interval<u32>) -> Self {
        self.set_lengths(LengthAxisSelector::Range(range), "length_range")
    }

    /// Set the row selector or record a row-axis selector conflict.
    fn set_rows(mut self, selector: LengthRowSelector, selector_name: &'static str) -> Self {
        if let Some(previous_selector_name) = self.rows.selector_name() {
            self.record_axis_conflict("row", previous_selector_name, selector_name);
        } else {
            self.rows = selector;
        }
        self
    }

    /// Set the length selector or record a fragment length-axis selector conflict.
    fn set_lengths(mut self, selector: LengthAxisSelector, selector_name: &'static str) -> Self {
        if let Some(previous_selector_name) = self.lengths.selector_name() {
            self.record_axis_conflict("fragment length", previous_selector_name, selector_name);
        } else {
            self.lengths = selector;
        }
        self
    }

    /// Store the first selector conflict so `read()` reports it as a normal error.
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

    /// Return any selector conflict recorded while building the selector.
    fn ensure_no_selector_conflict(&self) -> Result<()> {
        if let Some(selection_error) = &self.selection_error {
            bail!("{selection_error}");
        }
        Ok(())
    }

    /// Read the selected counts into an owned matrix with axis metadata.
    pub fn read(self) -> OutputLoaderResult<LengthCountSelection> {
        self.ensure_no_selector_conflict()?;
        let (length_bin_indices, length_range) = match self.lengths {
            LengthAxisSelector::All => (None, None),
            LengthAxisSelector::Indices(indices) => (Some(indices), None),
            LengthAxisSelector::Range(range) => (None, Some(range)),
        };
        let length_bin_indices = length_bin_indices.as_deref();

        let selection = match self.rows {
            LengthRowSelector::All => {
                self.output
                    .select_counts(None, length_bin_indices, length_range)
            }
            LengthRowSelector::Rows(indices) => self.output.select_counts(
                Some(indices.as_slice()),
                length_bin_indices,
                length_range,
            ),
            LengthRowSelector::Windows(indices) => self.output.select_window_counts(
                Some(indices.as_slice()),
                length_bin_indices,
                length_range,
            ),
            LengthRowSelector::Groups(indices) => self.output.select_group_counts(
                Some(indices.as_slice()),
                length_bin_indices,
                length_range,
            ),
            LengthRowSelector::GroupNames(names) => self.output.select_group_counts_by_name(
                Some(names.as_slice()),
                length_bin_indices,
                length_range,
            ),
        }?;
        Ok(selection)
    }
}

/// Row-axis selector state recorded by `LengthsSelector`.
#[derive(Debug, Clone)]
enum LengthRowSelector {
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

impl LengthRowSelector {
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

/// Fragment length-axis selector state recorded by `LengthsSelector`.
#[derive(Debug, Clone)]
enum LengthAxisSelector {
    /// Select all fragment length bins.
    All,
    /// Select fragment length bins by index.
    Indices(Vec<usize>),
    /// Select fragment length bins overlapping a half-open length range.
    Range(Interval<u32>),
}

impl LengthAxisSelector {
    /// Return the public selector method name for conflict messages.
    fn selector_name(&self) -> Option<&'static str> {
        match self {
            Self::All => None,
            Self::Indices(_) => Some("length_bins"),
            Self::Range(_) => Some("length_range"),
        }
    }
}

/// Compact metadata for a loaded length-count table.
///
/// This is intended for quick inspection and logging. It collects the output
/// settings that otherwise live behind separate accessors on `LengthsOutput`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LengthOutputMetadata {
    /// Whether count rows are global, genomic windows, or grouped-BED groups.
    pub row_mode: LengthOutputMode,
    /// Number of count rows.
    pub row_count: usize,
    /// Number of fragment length bins.
    pub length_bin_count: usize,
    /// Inclusive lower bound of the first fragment length bin.
    pub min_fragment_length: Option<u32>,
    /// Exclusive upper bound of the last fragment length bin.
    pub max_fragment_length_exclusive: Option<u32>,
}

impl fmt::Display for LengthOutputMetadata {
    /// Render one-line output context for logs or interactive inspection.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "row_mode={}, row_count={}, length_bin_count={}, fragment_length_range={}",
            describe_length_output_mode(self.row_mode),
            self.row_count,
            self.length_bin_count,
            describe_fragment_length_range(
                self.min_fragment_length,
                self.max_fragment_length_exclusive
            )
        )
    }
}

fn describe_fragment_length_range(
    min_fragment_length: Option<u32>,
    max_fragment_length_exclusive: Option<u32>,
) -> String {
    match (min_fragment_length, max_fragment_length_exclusive) {
        (Some(min_fragment_length), Some(max_fragment_length_exclusive)) => {
            format!("[{min_fragment_length}, {max_fragment_length_exclusive}) bp")
        }
        _ => "none".to_string(),
    }
}

fn describe_length_output_mode(row_mode: LengthOutputMode) -> &'static str {
    match row_mode {
        LengthOutputMode::Global => "global",
        LengthOutputMode::Windows => "windows",
        LengthOutputMode::Groups => "groups",
    }
}

/// Row metadata for a loaded length-count table.
#[derive(Debug, Clone, PartialEq)]
pub enum LengthRowMetadata {
    /// One global row with no additional metadata.
    Global,
    /// One row per genomic window.
    Windows(Vec<WindowRow>),
    /// One row per grouped-BED group.
    Groups(Vec<LengthGroupRow>),
}

impl LengthRowMetadata {
    /// Return the row aggregation mode.
    pub fn mode(&self) -> LengthOutputMode {
        match self {
            Self::Global => LengthOutputMode::Global,
            Self::Windows(_) => LengthOutputMode::Windows,
            Self::Groups(_) => LengthOutputMode::Groups,
        }
    }
}

/// Metadata for one grouped length-count output row.
#[derive(Debug, Clone, PartialEq)]
pub struct LengthGroupRow {
    /// Zero-based row index in output-file order.
    pub index: usize,
    /// Public group name from the grouped BED input.
    pub name: String,
    /// Number of eligible windows contributing to the group row.
    pub eligible_windows: u64,
    /// Fraction of grouped bases covered by blacklisted intervals, when written.
    pub blacklisted_fraction: Option<f64>,
}

/// Owned count matrix extracted from a `LengthsOutput`.
///
/// Rows are in the order requested by the selector. For `groups_by_name()`, this
/// is the same order as the requested group names. `row_metadata()` stores the
/// selected row metadata in that order, and `row_indices()` maps each selected
/// row back to its zero-based row in the source output. Columns are in
/// `length_bins()` order. Use the matrix or row-major values to build derived
/// summaries, tables, or model inputs from the selected subset.
#[derive(Debug, Clone, PartialEq)]
pub struct LengthCountSelection {
    row_metadata: LengthRowMetadata,
    row_indices: Vec<usize>,
    length_bins: Vec<LengthBin>,
    counts: DenseMatrix<f64>,
}

impl LengthCountSelection {
    /// Return selected row metadata in selection order.
    pub fn row_metadata(&self) -> &LengthRowMetadata {
        &self.row_metadata
    }

    /// Return selected window metadata, or an error if the selection is not windowed.
    pub fn window_metadata(&self) -> OutputLoaderResult<&[WindowRow]> {
        match &self.row_metadata {
            LengthRowMetadata::Windows(windows) => Ok(windows),
            _ => Err(OutputLoaderError::message(
                "length count selection is not windowed",
            )),
        }
    }

    /// Return selected group metadata, or an error if the selection is not grouped.
    pub fn group_metadata(&self) -> OutputLoaderResult<&[LengthGroupRow]> {
        match &self.row_metadata {
            LengthRowMetadata::Groups(groups) => Ok(groups),
            _ => Err(OutputLoaderError::message(
                "length count selection is not grouped",
            )),
        }
    }

    /// Return selected source row indices in selection order.
    pub fn row_indices(&self) -> &[usize] {
        &self.row_indices
    }

    /// Return selected length bins in selection order.
    pub fn length_bins(&self) -> &[LengthBin] {
        &self.length_bins
    }

    /// Return the selected matrix shape as `(rows, length_bins)`.
    pub fn shape(&self) -> (usize, usize) {
        self.counts.shape()
    }

    /// Return the number of selected rows.
    pub fn row_count(&self) -> usize {
        self.row_indices.len()
    }

    /// Return the number of selected fragment length bins.
    pub fn length_bin_count(&self) -> usize {
        self.length_bins.len()
    }

    /// Return selected counts with shape `(n_rows, n_length_bins)`.
    pub fn counts(&self) -> &DenseMatrix<f64> {
        &self.counts
    }

    /// Return selected count values in row-major order.
    pub fn counts_row_major(&self) -> &[f64] {
        self.counts.values_row_major()
    }

    /// Return one selected count value, if both selection indices are in bounds.
    ///
    /// Parameters
    /// ----------
    /// - `row_index`:
    ///   Zero-based row index within the selected matrix.
    /// - `length_bin_index`:
    ///   Zero-based length-bin index within the selected matrix.
    pub fn count(&self, row_index: usize, length_bin_index: usize) -> Option<f64> {
        self.counts.get(row_index, length_bin_index).copied()
    }
}

/// Parser for one `cfdna lengths` output file.
///
/// The parser owns only the input path. `load()` opens the file, detects the
/// row schema from the header, parses all rows, and then hands off to
/// `LengthsSchema::finish()` to build the public `LengthsOutput`.
struct LengthsParser {
    path: PathBuf,
    read_in_background: bool,
}

impl LengthsParser {
    /// Store the input path until `load()` opens it.
    fn new(path: &Path, read_in_background: bool) -> Self {
        Self {
            path: path.to_path_buf(),
            read_in_background,
        }
    }

    /// Read the TSV header, parse all data rows, and build a `LengthsOutput`.
    fn load(&self) -> Result<LengthsOutput> {
        let mut reader = if self.read_in_background {
            open_text_reader_in_background(&self.path)
        } else {
            open_text_reader(&self.path)
        }
        .with_context(|| format!("open lengths output {}", self.path.display()))?;
        let mut header_line = String::new();
        ensure!(
            reader
                .read_line(&mut header_line)
                .with_context(|| format!("read header from {}", self.path.display()))?
                > 0,
            "lengths output {} is empty; header required",
            self.path.display()
        );

        let header = split_header(&header_line);
        ensure!(
            !header.is_empty(),
            "lengths output {} has an empty header",
            self.path.display()
        );
        let schema = LengthsSchema::from_header(&self.path, &header)?;

        let mut rows = Vec::new();
        for (line_offset, line_result) in reader.lines().enumerate() {
            // User-facing error messages use 1-based file lines, and line 1 is the header
            let line_number = line_offset + 2;
            let line = line_result
                .with_context(|| format!("read line {line_number} from {}", self.path.display()))?;
            if line.is_empty() {
                bail!(
                    "lengths output {} line {line_number} is empty",
                    self.path.display()
                );
            }
            rows.push(schema.parse_row(&self.path, line_number, &line)?);
        }

        schema.finish(&self.path, rows)
    }
}

/// Internal row layout detected from the leading TSV columns.
///
/// The flags record whether optional blacklist metadata is present, which
/// controls how many metadata columns are consumed before count columns begin.
#[derive(Debug, Clone, Copy)]
enum RowMode {
    Global,
    Windows { has_blacklisted_fraction: bool },
    Groups { has_blacklisted_fraction: bool },
}

/// Parsed TSV schema for a `cfdna lengths` table.
///
/// This keeps the row mode, parsed fragment length bins, and count-column
/// offset together so row parsing cannot drift away from header interpretation.
struct LengthsSchema {
    mode: RowMode,
    length_bins: Vec<LengthBin>,
    metadata_columns: usize,
}

impl LengthsSchema {
    /// Detect the row schema and length-bin columns from the TSV header.
    fn from_header(path: &Path, header: &[&str]) -> Result<Self> {
        // The leading columns identify the row mode. Everything after those
        // metadata columns must be count columns on the fragment length axis
        let (mode, metadata_columns) = if header[0].starts_with("count_") {
            (RowMode::Global, 0)
        } else if header.starts_with(&["chrom", "start", "end", "blacklisted_fraction"]) {
            (
                RowMode::Windows {
                    has_blacklisted_fraction: true,
                },
                4,
            )
        } else if header.starts_with(&["chrom", "start", "end"]) {
            (
                RowMode::Windows {
                    has_blacklisted_fraction: false,
                },
                3,
            )
        } else if header.starts_with(&["group_name", "eligible_windows", "blacklisted_fraction"]) {
            (
                RowMode::Groups {
                    has_blacklisted_fraction: true,
                },
                3,
            )
        } else if header.starts_with(&["group_name", "eligible_windows"]) {
            (
                RowMode::Groups {
                    has_blacklisted_fraction: false,
                },
                2,
            )
        } else {
            bail!(
                "lengths output {} has unsupported header; expected count columns, window metadata, or group metadata",
                path.display()
            );
        };

        ensure!(
            metadata_columns < header.len(),
            "lengths output {} has no count columns",
            path.display()
        );
        let length_bins = parse_length_bins(path, &header[metadata_columns..])?;
        Ok(Self {
            mode,
            length_bins,
            metadata_columns,
        })
    }

    /// Parse one TSV data row according to the header-derived schema.
    fn parse_row(&self, path: &Path, line_number: usize, line: &str) -> Result<ParsedLengthRow> {
        let fields: Vec<&str> = line.split('\t').collect();
        ensure!(
            fields.len() == self.metadata_columns + self.length_bins.len(),
            "lengths output {} line {line_number} has {} columns, expected {}",
            path.display(),
            fields.len(),
            self.metadata_columns + self.length_bins.len()
        );

        // Count columns are always the suffix after the metadata columns, so
        // the same parsing path works for global, windowed, and grouped rows
        let counts = fields[self.metadata_columns..]
            .iter()
            .enumerate()
            .map(|(count_index, value)| {
                parse_count_field(
                    path,
                    line_number,
                    &format!("count column {count_index}"),
                    value,
                )
            })
            .collect::<Result<Vec<_>>>()?;

        // Parse only the metadata shape promised by the header-derived schema
        let row_metadata = match self.mode {
            RowMode::Global => ParsedRowMetadata::Global,
            RowMode::Windows {
                has_blacklisted_fraction,
            } => {
                let start = parse_u64_field(path, line_number, "start", fields[1])?;
                let end = parse_u64_field(path, line_number, "end", fields[2])?;
                let blacklisted_fraction = if has_blacklisted_fraction {
                    Some(parse_fraction_field(
                        path,
                        line_number,
                        "blacklisted_fraction",
                        fields[3],
                    )?)
                } else {
                    None
                };
                ParsedRowMetadata::Window {
                    chrom: fields[0].to_string(),
                    interval: crate::interval::Interval::new(start, end).map_err(|error| {
                        anyhow::anyhow!(
                            "lengths output {} line {line_number} has invalid window interval: {error}",
                            path.display()
                        )
                    })?,
                    blacklisted_fraction,
                }
            }
            RowMode::Groups {
                has_blacklisted_fraction,
            } => {
                let blacklisted_fraction = if has_blacklisted_fraction {
                    Some(parse_fraction_field(
                        path,
                        line_number,
                        "blacklisted_fraction",
                        fields[2],
                    )?)
                } else {
                    None
                };
                ParsedRowMetadata::Group {
                    name: fields[0].to_string(),
                    eligible_windows: parse_u64_field(
                        path,
                        line_number,
                        "eligible_windows",
                        fields[1],
                    )?,
                    blacklisted_fraction,
                }
            }
        };

        Ok(ParsedLengthRow {
            metadata: row_metadata,
            counts,
        })
    }

    /// Convert parsed TSV rows into the public lengths output object.
    ///
    /// `parse_row()` keeps each input line as parsed metadata plus count
    /// values. This method validates final table-level invariants, such as
    /// requiring exactly one row for global outputs, assigns public row indices,
    /// moves all count values into one row-major `DenseMatrix`, and attaches the
    /// row metadata variant that matches the header-derived row mode.
    ///
    /// Parameters
    /// ----------
    /// - `path`:
    ///   Input path used only for user-facing error messages.
    /// - `rows`:
    ///   Parsed data rows in file order.
    ///
    /// Returns
    /// -------
    /// - `LengthsOutput`:
    ///   Final loaded output with row metadata, length bins, and dense counts.
    fn finish(self, path: &Path, rows: Vec<ParsedLengthRow>) -> Result<LengthsOutput> {
        let row_count = rows.len();
        ensure!(
            row_count > 0,
            "lengths output {} has no data rows",
            path.display()
        );
        let length_bin_count = self.length_bins.len();

        // Convert private parsed metadata into the public metadata enum that
        // matches the row mode detected from the header. Counts are moved into
        // the final dense storage in the same pass so wide tables do not keep a
        // second full copy of count values alive during finishing.
        match self.mode {
            RowMode::Global => {
                ensure!(
                    row_count == 1,
                    "global lengths output {} has {row_count} data rows, expected 1",
                    path.display()
                );
                let mut values = Vec::with_capacity(row_count.saturating_mul(length_bin_count));
                for row in rows {
                    let ParsedLengthRow { metadata, counts } = row;
                    match metadata {
                        ParsedRowMetadata::Global => values.extend(counts),
                        _ => bail!("internal lengths loader row-mode mismatch"),
                    }
                }
                let counts = DenseMatrix::from_row_major(values, row_count, length_bin_count)?;
                Ok(LengthsOutput {
                    row_metadata: LengthRowMetadata::Global,
                    length_bins: self.length_bins,
                    counts,
                })
            }
            RowMode::Windows { .. } => {
                let mut values = Vec::with_capacity(row_count.saturating_mul(length_bin_count));
                let windows = rows
                    .into_iter()
                    .enumerate()
                    .map(|(index, row)| {
                        let ParsedLengthRow { metadata, counts } = row;
                        values.extend(counts);
                        match metadata {
                            ParsedRowMetadata::Window {
                                chrom,
                                interval,
                                blacklisted_fraction,
                            } => Ok(WindowRow {
                                index,
                                chrom,
                                interval,
                                blacklisted_fraction,
                            }),
                            _ => bail!("internal lengths loader row-mode mismatch"),
                        }
                    })
                    .collect::<Result<Vec<_>>>()?;
                let counts = DenseMatrix::from_row_major(values, row_count, length_bin_count)?;
                Ok(LengthsOutput {
                    row_metadata: LengthRowMetadata::Windows(windows),
                    length_bins: self.length_bins,
                    counts,
                })
            }
            RowMode::Groups { .. } => {
                let mut values = Vec::with_capacity(row_count.saturating_mul(length_bin_count));
                let groups = rows
                    .into_iter()
                    .enumerate()
                    .map(|(index, row)| {
                        let ParsedLengthRow { metadata, counts } = row;
                        values.extend(counts);
                        match metadata {
                            ParsedRowMetadata::Group {
                                name,
                                eligible_windows,
                                blacklisted_fraction,
                            } => Ok(LengthGroupRow {
                                index,
                                name,
                                eligible_windows,
                                blacklisted_fraction,
                            }),
                            _ => bail!("internal lengths loader row-mode mismatch"),
                        }
                    })
                    .collect::<Result<Vec<_>>>()?;
                let group_names = groups
                    .iter()
                    .map(|group| group.name.as_str())
                    .collect::<Vec<_>>();
                ensure_unique_labels(&group_names, "group_names")?;
                let counts = DenseMatrix::from_row_major(values, row_count, length_bin_count)?;
                Ok(LengthsOutput {
                    row_metadata: LengthRowMetadata::Groups(groups),
                    length_bins: self.length_bins,
                    counts,
                })
            }
        }
    }
}

/// One parsed data line before rows are collected into the final output matrix.
///
/// Metadata is kept separate from counts because the final public row metadata
/// type depends on the table mode, while counts are always stored in the same
/// row-major matrix layout.
struct ParsedLengthRow {
    metadata: ParsedRowMetadata,
    counts: Vec<f64>,
}

/// Row identity parsed from one non-header line.
enum ParsedRowMetadata {
    Global,
    Window {
        chrom: String,
        interval: crate::interval::Interval<u64>,
        blacklisted_fraction: Option<f64>,
    },
    Group {
        name: String,
        eligible_windows: u64,
        blacklisted_fraction: Option<f64>,
    },
}

/// Split a TSV header line into tab-delimited column names.
fn split_header(header_line: &str) -> Vec<&str> {
    header_line
        .trim_end_matches(['\r', '\n'])
        .split('\t')
        .collect()
}

/// Parse all count column headers into fragment length bins.
fn parse_length_bins(path: &Path, headers: &[&str]) -> Result<Vec<LengthBin>> {
    let mut previous_end = None;
    headers
        .iter()
        .enumerate()
        .map(|(index, header)| {
            let bin = parse_length_bin_header(path, index, header)?;
            if let Some(previous_end_bp) = previous_end {
                ensure!(
                    bin.start() == previous_end_bp,
                    "lengths output {} has non-contiguous length count columns at {}: previous end {}, next start {}",
                    path.display(),
                    header,
                    previous_end_bp,
                    bin.start()
                );
            }
            previous_end = Some(bin.end());
            Ok(bin)
        })
        .collect()
}

/// Parse one `count_*` column header into a length-bin interval.
fn parse_length_bin_header(path: &Path, index: usize, header: &str) -> Result<LengthBin> {
    let suffix = header.strip_prefix("count_").with_context(|| {
        format!(
            "lengths output {} count header '{}' must start with count_",
            path.display(),
            header
        )
    })?;
    let parts = suffix.split('_').collect::<Vec<_>>();
    let (start_bp, end_bp) = match parts.as_slice() {
        [length] => {
            let start_bp = parse_u32_header(path, header, length)?;
            let end_bp = start_bp.checked_add(1).with_context(|| {
                format!(
                    "lengths output {} count header '{}' overflows single-bp bin end",
                    path.display(),
                    header
                )
            })?;
            (start_bp, end_bp)
        }
        [start, end] => (
            parse_u32_header(path, header, start)?,
            parse_u32_header(path, header, end)?,
        ),
        _ => {
            bail!(
                "lengths output {} count header '{}' must be count_<length> or count_<start>_<end>",
                path.display(),
                header
            );
        }
    };
    ensure!(
        end_bp > start_bp,
        "lengths output {} count header '{}' has end <= start",
        path.display(),
        header
    );
    ensure!(
        start_bp >= MIN_ACGT_BASES_FOR_GC_FRACTION,
        "lengths output {} count header '{}' starts below minimum supported fragment length {}",
        path.display(),
        header,
        MIN_ACGT_BASES_FOR_GC_FRACTION
    );
    ensure!(
        end_bp <= MAX_SUPPORTED_FRAGMENT_LENGTH + 1,
        "lengths output {} count header '{}' ends above maximum supported fragment length edge {}",
        path.display(),
        header,
        MAX_SUPPORTED_FRAGMENT_LENGTH + 1
    );
    LengthBin::new(start_bp, end_bp, index).map_err(|error| {
        anyhow::anyhow!(
            "lengths output {} count header '{}' has invalid length interval: {error}",
            path.display(),
            header
        )
    })
}

/// Parse one unsigned integer embedded in a count column header.
fn parse_u32_header(path: &Path, header: &str, value: &str) -> Result<u32> {
    value.parse::<u32>().with_context(|| {
        format!(
            "lengths output {} count header '{}' contains invalid fragment length '{}'",
            path.display(),
            header,
            value
        )
    })
}

/// Parse one non-negative integer data field with file and line context.
fn parse_u64_field(path: &Path, line_number: usize, field_name: &str, value: &str) -> Result<u64> {
    value.parse::<u64>().with_context(|| {
        format!(
            "lengths output {} line {line_number} has invalid {field_name} '{}'",
            path.display(),
            value
        )
    })
}

/// Parse one floating-point data field with file and line context.
fn parse_f64_field(path: &Path, line_number: usize, field_name: &str, value: &str) -> Result<f64> {
    value.parse::<f64>().with_context(|| {
        format!(
            "lengths output {} line {line_number} has invalid {field_name} '{}'",
            path.display(),
            value
        )
    })
}

/// Parse one finite non-negative count value.
fn parse_count_field(
    path: &Path,
    line_number: usize,
    field_name: &str,
    value: &str,
) -> Result<f64> {
    let count = parse_f64_field(path, line_number, field_name, value)?;
    ensure!(
        count.is_finite() && count >= 0.0,
        "lengths output {} line {line_number} has {field_name} outside finite and non-negative range: {count}",
        path.display()
    );
    Ok(count)
}

/// Parse one finite fraction data field and reject values outside `[0, 1]`.
fn parse_fraction_field(
    path: &Path,
    line_number: usize,
    field_name: &str,
    value: &str,
) -> Result<f64> {
    let fraction = parse_f64_field(path, line_number, field_name, value)?;
    ensure!(
        fraction.is_finite() && (0.0..=1.0).contains(&fraction),
        "lengths output {} line {line_number} has {field_name} outside [0, 1]: {fraction}",
        path.display()
    );
    Ok(fraction)
}
