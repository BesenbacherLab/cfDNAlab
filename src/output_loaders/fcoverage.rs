//! Loader for non-positional `cfdna fcoverage` aggregate tables.
//!
//! A loaded fcoverage output has three pieces of information: row metadata,
//! a signal label, and row values. Row metadata says whether rows are genomic
//! windows or grouped-BED groups. The data are either scalar aggregate values
//! for `average` or `total` outputs, or summary statistics for `summary_stats`
//! outputs.
//!
//! Type overview:
//!
//! ```text
//! load_fcoverage_output(path)
//!     -> FCoverageOutput
//!         row_metadata: FCoverageRowMetadata
//!             Windows(Vec<FCoverageWindowRow>)
//!             Groups(Vec<FCoverageGroupRow>)
//!         signal: FCoverageSignal
//!         data: FCoverageData
//!             Values { value_mode, values: Vec<f64> }
//!             SummaryStats(Vec<FCoverageSummaryStats>)
//!
//! FCoverageOutput::select()
//!     -> FCoverageSelector
//!         -> read()
//!             -> FCoverageSelection
//!                 Values(FCoverageValueSelection)
//!                 SummaryStats(FCoverageSummaryStatsSelection)
//!                 row_metadata: FCoverageRowMetadata
//! ```
//!
//! Use `FCoverageOutput` directly for shared information such as `row_mode()`,
//! `signal()`, and `row_count()`. Use `values()` and `value_mode()` for scalar
//! aggregate outputs, or `summary_stats()` for summary-stat outputs. Use
//! `select()` to extract rows from either data mode with the same row-selector
//! API. Selections include selected row metadata next to selected values, so
//! each value or summary-stat row can be paired with its window or group.

use crate::{
    interval::Interval,
    output_loaders::{
        OutputLoaderError, OutputLoaderResult,
        common::{ensure_unique_indices, resolve_row_indices},
    },
    shared::io::open_text_reader,
};
use anyhow::{Context, Result, bail, ensure};
use fxhash::{FxHashMap, FxHashSet};
use std::{
    io::BufRead,
    path::{Path, PathBuf},
};

/// Load a non-positional `cfdna fcoverage` aggregate table.
///
/// Supported inputs are the aggregate TSV outputs for `average`, `total`, and
/// `summary_stats`, including grouped outputs and length-normalized
/// `fragment_mass` headers. Positional bedGraph and per-window positional TSV
/// outputs are intentionally out of scope and are rejected before parsing.
///
/// Parameters
/// ----------
/// - `path`:
///     Path to a non-positional aggregate `cfdna fcoverage` TSV output. Plain
///     text, gzip, and zstd compressed files are supported.
///
/// Returns
/// -------
/// - `FCoverageOutput`:
///     Loaded row metadata, signal label, and scalar values or summary stats.
///
/// ```no_run
/// use cfdnalab::output_loaders::load_fcoverage_output;
///
/// let fcoverage = load_fcoverage_output("sample.fcoverage.average.tsv.zst")?;
/// let selected = fcoverage.select().windows(&[0, 2, 4]).read()?;
///
/// for (window, value) in selected
///     .window_metadata()?
///     .iter()
///     .zip(selected.values()?.iter().copied())
/// {
///     println!(
///         "{}:{}-{} {value}",
///         window.chrom,
///         window.interval.start(),
///         window.interval.end()
///     );
/// }
///
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn load_fcoverage_output(path: impl AsRef<Path>) -> OutputLoaderResult<FCoverageOutput> {
    let path = path.as_ref();
    ensure_non_positional_path(path)?;
    FCoverageParser::new(path, None).load().map_err(Into::into)
}

/// Load a non-positional grouped `cfdna fcoverage` aggregate table with group names.
///
/// The aggregate TSV stores numeric `group_idx` values. Passing the matching
/// `group_index.tsv` file attaches group names to grouped row metadata and
/// enables `group_index()`, `has_group()`, and `groups_by_name()` lookups.
///
/// Parameters
/// ----------
/// - `path`:
///     Path to a grouped non-positional aggregate `cfdna fcoverage` TSV output.
///     Plain text, gzip, and zstd compressed files are supported.
/// - `group_index_path`:
///     Path to the matching group-index file with `group_idx` and `group_name`
///     columns.
///
/// Returns
/// -------
/// - `FCoverageOutput`:
///     Loaded row metadata with group names, signal label, and scalar values or
///     summary stats.
///
/// ```no_run
/// use cfdnalab::output_loaders::load_fcoverage_output_with_group_index;
///
/// let fcoverage = load_fcoverage_output_with_group_index(
///     "sample.fcoverage.total_on_unique_bases.tsv.zst",
///     "sample.group_index.tsv",
/// )?;
/// let selected = fcoverage
///     .select()
///     .groups_by_name(&["promoters", "enhancers"])
///     .read()?;
///
/// for (group, value) in selected
///     .group_metadata()?
///     .iter()
///     .zip(selected.values()?.iter().copied())
/// {
///     let group_name = group.name.as_deref().unwrap_or("<unnamed>");
///     println!("{group_name}\t{value}");
/// }
///
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn load_fcoverage_output_with_group_index(
    path: impl AsRef<Path>,
    group_index_path: impl AsRef<Path>,
) -> OutputLoaderResult<FCoverageOutput> {
    let path = path.as_ref();
    ensure_non_positional_path(path)?;
    FCoverageParser::new(path, Some(group_index_path.as_ref()))
        .load()
        .map_err(Into::into)
}

/// Loaded non-positional output from `cfdna fcoverage`.
#[derive(Debug, Clone, PartialEq)]
pub struct FCoverageOutput {
    row_metadata: FCoverageRowMetadata,
    signal: FCoverageSignal,
    filename_metadata: FCoverageFilenameMetadata,
    data: FCoverageData,
}

impl FCoverageOutput {
    /// Return whether rows are genomic windows or grouped rows.
    pub fn row_mode(&self) -> FCoverageRowMode {
        self.row_metadata.mode()
    }

    /// Return row metadata describing what each output row represents.
    pub fn row_metadata(&self) -> &FCoverageRowMetadata {
        &self.row_metadata
    }

    /// Return command-mode hints parsed from the output filename.
    ///
    /// These values come only from canonical cfDNAlab filename parts such as
    /// `fcoverage.total_on_unique_bases.tsv.zst` or
    /// `length_normalized.restored_mean`. Renamed files report `Unknown` for
    /// fields that are not present in the filename.
    pub fn filename_metadata(&self) -> &FCoverageFilenameMetadata {
        &self.filename_metadata
    }

    /// Return genomic window metadata, or an error if this is not a windowed output.
    pub fn window_metadata(&self) -> OutputLoaderResult<&[FCoverageWindowRow]> {
        match &self.row_metadata {
            FCoverageRowMetadata::Windows(windows) => Ok(windows),
            _ => Err(OutputLoaderError::message(
                "fcoverage output is not windowed",
            )),
        }
    }

    /// Return group metadata, or an error if this is not a grouped output.
    pub fn group_metadata(&self) -> OutputLoaderResult<&[FCoverageGroupRow]> {
        match &self.row_metadata {
            FCoverageRowMetadata::Groups(groups) => Ok(groups),
            _ => Err(OutputLoaderError::message(
                "fcoverage output is not grouped",
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
    ///     Zero-based row index in the window metadata.
    pub fn window(&self, row_index: usize) -> OutputLoaderResult<Option<&FCoverageWindowRow>> {
        Ok(self.window_metadata()?.get(row_index))
    }

    /// Return one group row by zero-based row index.
    ///
    /// This returns an error if the loaded output is not grouped.
    ///
    /// Parameters
    /// ----------
    /// - `row_index`:
    ///     Zero-based row index in the group metadata.
    pub fn group(&self, row_index: usize) -> OutputLoaderResult<Option<&FCoverageGroupRow>> {
        Ok(self.group_metadata()?.get(row_index))
    }

    /// Return the row index for a loaded group name.
    ///
    /// Group names are only available when the output was loaded with
    /// `load_fcoverage_output_with_group_index()`.
    ///
    /// Parameters
    /// ----------
    /// - `group_name`:
    ///     Exact group name from the group-index file.
    pub fn group_index(&self, group_name: &str) -> OutputLoaderResult<usize> {
        let groups = self.group_metadata()?;
        if !groups.iter().any(|group| group.name.is_some()) {
            return Err(OutputLoaderError::message(
                "fcoverage output has no group names loaded; use load_fcoverage_output_with_group_index() with a group-index file",
            ));
        }
        Ok(groups
            .iter()
            .find(|group| group.name.as_deref() == Some(group_name))
            .map(|group| group.index)
            .with_context(|| format!("fcoverage output has no group named '{group_name}'"))?)
    }

    /// Return whether a loaded grouped output contains a group name.
    ///
    /// This returns `false` for windowed outputs and grouped outputs loaded
    /// without a group-index file.
    pub fn has_group(&self, group_name: &str) -> bool {
        self.group_metadata()
            .map(|groups| {
                groups
                    .iter()
                    .any(|group| group.name.as_deref() == Some(group_name))
            })
            .unwrap_or(false)
    }

    /// Return the signal label used in value columns, such as `coverage` or `fragment_mass`.
    pub fn signal(&self) -> &FCoverageSignal {
        &self.signal
    }

    /// Return whether the output contains scalar values or summary stats.
    pub fn data(&self) -> &FCoverageData {
        &self.data
    }

    /// Return the scalar aggregate mode, or an error for summary-stat outputs.
    pub fn value_mode(&self) -> OutputLoaderResult<FCoverageValueMode> {
        match &self.data {
            FCoverageData::Values { value_mode, .. } => Ok(*value_mode),
            FCoverageData::SummaryStats(_) => Err(OutputLoaderError::message(
                "fcoverage output contains summary stats",
            )),
        }
    }

    /// Return the number of output rows.
    pub fn row_count(&self) -> usize {
        match &self.data {
            FCoverageData::Values { values, .. } => values.len(),
            FCoverageData::SummaryStats(stats) => stats.len(),
        }
    }

    /// Return scalar aggregate values in output-row order.
    pub fn values(&self) -> OutputLoaderResult<&[f64]> {
        match &self.data {
            FCoverageData::Values { values, .. } => Ok(values),
            FCoverageData::SummaryStats(_) => Err(OutputLoaderError::message(
                "fcoverage output contains summary stats",
            )),
        }
    }

    /// Return one scalar aggregate value, if `row_index` is in bounds.
    ///
    /// Parameters
    /// ----------
    /// - `row_index`:
    ///     Zero-based output row index.
    pub fn value(&self, row_index: usize) -> OutputLoaderResult<Option<f64>> {
        Ok(self.values()?.get(row_index).copied())
    }

    /// Return summary statistics in output-row order.
    pub fn summary_stats(&self) -> OutputLoaderResult<&[FCoverageSummaryStats]> {
        match &self.data {
            FCoverageData::SummaryStats(stats) => Ok(stats),
            FCoverageData::Values { .. } => Err(OutputLoaderError::message(
                "fcoverage output contains scalar values",
            )),
        }
    }

    /// Return one summary-stat row, if `row_index` is in bounds.
    ///
    /// Parameters
    /// ----------
    /// - `row_index`:
    ///     Zero-based output row index.
    pub fn summary_stat(
        &self,
        row_index: usize,
    ) -> OutputLoaderResult<Option<&FCoverageSummaryStats>> {
        Ok(self.summary_stats()?.get(row_index))
    }

    /// Start a row selection.
    ///
    /// A new selector initially selects all rows. Add a row constraint before
    /// calling `read()` when only some rows are needed.
    pub fn select(&self) -> FCoverageSelector<'_> {
        FCoverageSelector::new(self)
    }

    /// Return selected scalar aggregate values.
    ///
    /// Passing `None` for `row_indices` selects all rows.
    pub(crate) fn select_values(
        &self,
        row_indices: Option<&[usize]>,
    ) -> Result<FCoverageValueSelection> {
        self.select_values_with_label(row_indices, "row")
    }

    /// Return scalar aggregate values for selected window rows.
    ///
    /// Passing `None` for `window_indices` selects all windows.
    pub(crate) fn select_window_values(
        &self,
        window_indices: Option<&[usize]>,
    ) -> Result<FCoverageValueSelection> {
        self.ensure_row_mode(FCoverageRowMode::Windows)?;
        self.select_values_with_label(window_indices, "window")
    }

    /// Return scalar aggregate values for selected group rows.
    ///
    /// Passing `None` for `group_indices` selects all groups.
    pub(crate) fn select_group_values(
        &self,
        group_indices: Option<&[usize]>,
    ) -> Result<FCoverageValueSelection> {
        self.ensure_row_mode(FCoverageRowMode::Groups)?;
        self.select_values_with_label(group_indices, "group")
    }

    /// Return selected summary-stat rows.
    ///
    /// Passing `None` for `row_indices` selects all rows.
    pub(crate) fn select_summary_stats(
        &self,
        row_indices: Option<&[usize]>,
    ) -> Result<FCoverageSummaryStatsSelection> {
        self.select_summary_stats_with_label(row_indices, "row")
    }

    /// Return summary stats for selected window rows.
    ///
    /// Passing `None` for `window_indices` selects all windows.
    pub(crate) fn select_window_summary_stats(
        &self,
        window_indices: Option<&[usize]>,
    ) -> Result<FCoverageSummaryStatsSelection> {
        self.ensure_row_mode(FCoverageRowMode::Windows)?;
        self.select_summary_stats_with_label(window_indices, "window")
    }

    /// Return summary stats for selected group rows.
    ///
    /// Passing `None` for `group_indices` selects all groups.
    pub(crate) fn select_group_summary_stats(
        &self,
        group_indices: Option<&[usize]>,
    ) -> Result<FCoverageSummaryStatsSelection> {
        self.ensure_row_mode(FCoverageRowMode::Groups)?;
        self.select_summary_stats_with_label(group_indices, "group")
    }

    /// Resolve group names to row indices in requested order.
    fn resolve_group_name_indices<S: AsRef<str>>(&self, group_names: &[S]) -> Result<Vec<usize>> {
        self.ensure_row_mode(FCoverageRowMode::Groups)?;
        let groups = self.group_metadata()?;
        ensure!(
            groups.iter().any(|group| group.name.is_some()),
            "fcoverage output has no group names loaded; use load_fcoverage_output_with_group_index() with a group-index file"
        );
        let mut seen = FxHashSet::default();
        for group_name in group_names {
            let group_name = group_name.as_ref();
            ensure!(
                seen.insert(group_name.to_string()),
                "duplicate group name '{}'",
                group_name
            );
        }
        group_names
            .iter()
            .map(|group_name| {
                self.group_index(group_name.as_ref())
                    .map_err(anyhow::Error::from)
            })
            .collect()
    }

    /// Select scalar value rows after resolving a row label for errors.
    fn select_values_with_label(
        &self,
        row_indices: Option<&[usize]>,
        row_label: &str,
    ) -> Result<FCoverageValueSelection> {
        let row_indices = resolve_row_indices(row_indices, self.row_count(), row_label)?;
        ensure_unique_indices(&row_indices, row_label)?;
        let selected_row_metadata = self.selected_row_metadata(&row_indices, row_label)?;
        let values = self.values()?;
        let selected_values = row_indices
            .iter()
            .map(|&row_index| {
                values.get(row_index).copied().with_context(|| {
                    format!(
                        "{row_label} index {row_index} is outside 0..{}",
                        self.row_count()
                    )
                })
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(FCoverageValueSelection {
            row_metadata: selected_row_metadata,
            row_indices,
            value_mode: self.value_mode()?,
            signal: self.signal.clone(),
            values: selected_values,
        })
    }

    /// Select summary-stat rows after resolving a row label for errors.
    fn select_summary_stats_with_label(
        &self,
        row_indices: Option<&[usize]>,
        row_label: &str,
    ) -> Result<FCoverageSummaryStatsSelection> {
        let row_indices = resolve_row_indices(row_indices, self.row_count(), row_label)?;
        ensure_unique_indices(&row_indices, row_label)?;
        let selected_row_metadata = self.selected_row_metadata(&row_indices, row_label)?;
        let stats = self.summary_stats()?;
        let selected_stats = row_indices
            .iter()
            .map(|&row_index| {
                stats.get(row_index).cloned().with_context(|| {
                    format!(
                        "{row_label} index {row_index} is outside 0..{}",
                        self.row_count()
                    )
                })
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(FCoverageSummaryStatsSelection {
            row_metadata: selected_row_metadata,
            row_indices,
            signal: self.signal.clone(),
            stats: selected_stats,
        })
    }

    /// Build row metadata for selected source row indices.
    ///
    /// The returned metadata keeps the selector order. Generic row selectors
    /// preserve the loaded row mode, while window and group selectors have
    /// already checked the expected mode before this helper is called.
    fn selected_row_metadata(
        &self,
        row_indices: &[usize],
        row_label: &str,
    ) -> Result<FCoverageRowMetadata> {
        match &self.row_metadata {
            FCoverageRowMetadata::Windows(windows) => {
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
                Ok(FCoverageRowMetadata::Windows(selected_windows))
            }
            FCoverageRowMetadata::Groups(groups) => {
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
                Ok(FCoverageRowMetadata::Groups(selected_groups))
            }
        }
    }

    /// Return an error when the loaded row mode does not match the typed selector.
    fn ensure_row_mode(&self, expected: FCoverageRowMode) -> Result<()> {
        ensure!(
            self.row_mode() == expected,
            "fcoverage output is not {}",
            expected.description()
        );
        Ok(())
    }
}

/// Builder for selecting rows from an `FCoverageOutput`.
///
/// The builder starts with all rows selected. Set at most one row selector;
/// conflicting selector calls are reported by `read()` together with row-mode
/// validation and data-copy errors.
#[derive(Debug, Clone)]
pub struct FCoverageSelector<'a> {
    output: &'a FCoverageOutput,
    rows: FCoverageRowSelector,
    selection_error: Option<String>,
}

impl<'a> FCoverageSelector<'a> {
    /// Start a selector with all rows selected.
    fn new(output: &'a FCoverageOutput) -> Self {
        Self {
            output,
            rows: FCoverageRowSelector::All,
            selection_error: None,
        }
    }

    /// Select generic output rows by zero-based row index.
    ///
    /// Parameters
    /// ----------
    /// - `row_indices`:
    ///     Source row indices in output-file order. The returned selection keeps
    ///     this order and rejects duplicates.
    pub fn rows(self, row_indices: &[usize]) -> Self {
        self.set_rows(FCoverageRowSelector::Rows(row_indices.to_vec()), "rows")
    }

    /// Select window rows by zero-based window row index.
    ///
    /// `read()` returns an error if the loaded output is not windowed.
    ///
    /// Parameters
    /// ----------
    /// - `window_indices`:
    ///     Window row indices in output-file order. The returned selection keeps
    ///     this order and rejects duplicates.
    pub fn windows(self, window_indices: &[usize]) -> Self {
        self.set_rows(
            FCoverageRowSelector::Windows(window_indices.to_vec()),
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
    ///     Group row indices in output-file order. The returned selection keeps
    ///     this order and rejects duplicates.
    pub fn groups(self, group_indices: &[usize]) -> Self {
        self.set_rows(
            FCoverageRowSelector::Groups(group_indices.to_vec()),
            "groups",
        )
    }

    /// Select grouped rows by group name.
    ///
    /// `read()` returns an error if the loaded output is not grouped or if it
    /// was loaded without `load_fcoverage_output_with_group_index()`.
    ///
    /// Parameters
    /// ----------
    /// - `group_names`:
    ///     Group names in requested output order. The returned selection keeps
    ///     this order and rejects duplicates.
    pub fn groups_by_name<S: AsRef<str>>(self, group_names: &[S]) -> Self {
        self.set_rows(
            FCoverageRowSelector::GroupNames(
                group_names
                    .iter()
                    .map(|group_name| group_name.as_ref().to_string())
                    .collect(),
            ),
            "groups_by_name",
        )
    }

    /// Set the row selector or record a row-axis selector conflict.
    fn set_rows(mut self, selector: FCoverageRowSelector, selector_name: &'static str) -> Self {
        if let Some(previous_selector_name) = self.rows.selector_name() {
            self.record_row_conflict(previous_selector_name, selector_name);
        } else {
            self.rows = selector;
        }
        self
    }

    /// Store the first row selector conflict so `read()` reports it as an error.
    fn record_row_conflict(
        &mut self,
        previous_selector_name: &'static str,
        selector_name: &'static str,
    ) {
        if self.selection_error.is_none() {
            self.selection_error = Some(format!(
                "cannot combine {previous_selector_name}() and {selector_name}() on the row axis"
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

    /// Read selected rows into the selection type that matches the loaded file.
    pub fn read(self) -> OutputLoaderResult<FCoverageSelection> {
        self.ensure_no_selector_conflict()?;
        let selection = match self.output.data() {
            FCoverageData::Values { .. } => match self.rows {
                FCoverageRowSelector::All => self
                    .output
                    .select_values(None)
                    .map(FCoverageSelection::Values),
                FCoverageRowSelector::Rows(indices) => self
                    .output
                    .select_values(Some(indices.as_slice()))
                    .map(FCoverageSelection::Values),
                FCoverageRowSelector::Windows(indices) => self
                    .output
                    .select_window_values(Some(indices.as_slice()))
                    .map(FCoverageSelection::Values),
                FCoverageRowSelector::Groups(indices) => self
                    .output
                    .select_group_values(Some(indices.as_slice()))
                    .map(FCoverageSelection::Values),
                FCoverageRowSelector::GroupNames(names) => {
                    let indices = self.output.resolve_group_name_indices(&names)?;
                    self.output
                        .select_group_values(Some(indices.as_slice()))
                        .map(FCoverageSelection::Values)
                }
            },
            FCoverageData::SummaryStats(_) => match self.rows {
                FCoverageRowSelector::All => self
                    .output
                    .select_summary_stats(None)
                    .map(FCoverageSelection::SummaryStats),
                FCoverageRowSelector::Rows(indices) => self
                    .output
                    .select_summary_stats(Some(indices.as_slice()))
                    .map(FCoverageSelection::SummaryStats),
                FCoverageRowSelector::Windows(indices) => self
                    .output
                    .select_window_summary_stats(Some(indices.as_slice()))
                    .map(FCoverageSelection::SummaryStats),
                FCoverageRowSelector::Groups(indices) => self
                    .output
                    .select_group_summary_stats(Some(indices.as_slice()))
                    .map(FCoverageSelection::SummaryStats),
                FCoverageRowSelector::GroupNames(names) => {
                    let indices = self.output.resolve_group_name_indices(&names)?;
                    self.output
                        .select_group_summary_stats(Some(indices.as_slice()))
                        .map(FCoverageSelection::SummaryStats)
                }
            },
        }?;
        Ok(selection)
    }
}

/// Row-axis selector state recorded by `FCoverageSelector`.
#[derive(Debug, Clone)]
enum FCoverageRowSelector {
    /// Select all rows.
    All,
    /// Select generic output rows by index.
    Rows(Vec<usize>),
    /// Select window rows by index and require windowed output.
    Windows(Vec<usize>),
    /// Select grouped rows by index and require grouped output.
    Groups(Vec<usize>),
    /// Select grouped rows by loaded group name and require grouped output.
    GroupNames(Vec<String>),
}

impl FCoverageRowSelector {
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

/// Row aggregation mode detected from the aggregate table schema.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FCoverageRowMode {
    /// One row per genomic window.
    Windows,
    /// One row per grouped BED group.
    Groups,
}

impl FCoverageRowMode {
    /// Return a short label used in error messages.
    fn description(self) -> &'static str {
        match self {
            Self::Windows => "windowed",
            Self::Groups => "grouped",
        }
    }
}

/// Scalar aggregate mode for non-summary `fcoverage` tables.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FCoverageValueMode {
    /// Average value per eligible position.
    Average,
    /// Total value across eligible positions.
    Total,
}

/// Whether a fcoverage aggregate filename names ordinary or unique-base grouping.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FCoverageAggregationBasis {
    /// The filename uses `average`, `total`, or `summary_stats`.
    Ordinary,
    /// The filename uses `*_on_unique_bases`.
    UniqueBases,
    /// The filename does not contain a recognized aggregate action.
    Unknown,
}

/// Fragment length-normalization mode parsed from a fcoverage filename.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FCoverageLengthNormalization {
    /// The filename has no length-normalization marker.
    Off,
    /// The filename contains `length_normalized`.
    UnitMass,
    /// The filename contains `length_normalized.restored_mean`.
    RestoredMean,
    /// The filename does not contain enough recognized parts to decide.
    Unknown,
}

/// Command-mode hints parsed from a fcoverage output filename.
///
/// This metadata is intentionally lightweight. It reflects canonical filename
/// parts when they are present and uses `Unknown` when a user renamed the file
/// or supplied another valid aggregate table name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FCoverageFilenameMetadata {
    aggregation_basis: FCoverageAggregationBasis,
    length_normalization: FCoverageLengthNormalization,
}

impl FCoverageFilenameMetadata {
    /// Return whether the filename names ordinary or unique-base aggregation.
    pub fn aggregation_basis(&self) -> FCoverageAggregationBasis {
        self.aggregation_basis
    }

    /// Return length-normalization mode parsed from the filename.
    pub fn length_normalization(&self) -> FCoverageLengthNormalization {
        self.length_normalization
    }
}

/// Row metadata for a loaded fcoverage aggregate table.
#[derive(Debug, Clone, PartialEq)]
pub enum FCoverageRowMetadata {
    /// One row per genomic window.
    Windows(Vec<FCoverageWindowRow>),
    /// One row per grouped-BED group.
    Groups(Vec<FCoverageGroupRow>),
}

impl FCoverageRowMetadata {
    /// Return the row aggregation mode.
    pub fn mode(&self) -> FCoverageRowMode {
        match self {
            Self::Windows(_) => FCoverageRowMode::Windows,
            Self::Groups(_) => FCoverageRowMode::Groups,
        }
    }
}

/// Public signal label used by `fcoverage` aggregate value columns.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FCoverageSignal {
    label: String,
}

impl FCoverageSignal {
    /// Store the signal label parsed from value or summary-stat headers.
    fn new(label: impl Into<String>) -> Self {
        Self {
            label: label.into(),
        }
    }

    /// Return the exact suffix used in the output value columns.
    pub fn label(&self) -> &str {
        &self.label
    }
}

/// Loaded fcoverage data values.
#[derive(Debug, Clone, PartialEq)]
pub enum FCoverageData {
    /// Scalar aggregate values from `average` or `total` outputs.
    Values {
        /// Whether the values are averages or totals.
        value_mode: FCoverageValueMode,
        /// Values in output-row order.
        values: Vec<f64>,
    },
    /// Raw and derived summary statistics in output-row order.
    SummaryStats(Vec<FCoverageSummaryStats>),
}

/// Selected rows from an `FCoverageOutput`.
#[derive(Debug, Clone, PartialEq)]
pub enum FCoverageSelection {
    /// Selected scalar aggregate values.
    Values(FCoverageValueSelection),
    /// Selected summary-stat rows.
    SummaryStats(FCoverageSummaryStatsSelection),
}

impl FCoverageSelection {
    /// Return selected row metadata in selection order.
    pub fn row_metadata(&self) -> &FCoverageRowMetadata {
        match self {
            Self::Values(selection) => selection.row_metadata(),
            Self::SummaryStats(selection) => selection.row_metadata(),
        }
    }

    /// Return selected window metadata, or an error if the selection is not windowed.
    pub fn window_metadata(&self) -> OutputLoaderResult<&[FCoverageWindowRow]> {
        match self {
            Self::Values(selection) => selection.window_metadata(),
            Self::SummaryStats(selection) => selection.window_metadata(),
        }
    }

    /// Return selected group metadata, or an error if the selection is not grouped.
    pub fn group_metadata(&self) -> OutputLoaderResult<&[FCoverageGroupRow]> {
        match self {
            Self::Values(selection) => selection.group_metadata(),
            Self::SummaryStats(selection) => selection.group_metadata(),
        }
    }

    /// Return selected source row indices in requested order.
    pub fn row_indices(&self) -> &[usize] {
        match self {
            Self::Values(selection) => selection.row_indices(),
            Self::SummaryStats(selection) => selection.row_indices(),
        }
    }

    /// Return the number of selected rows.
    pub fn row_count(&self) -> usize {
        self.row_indices().len()
    }

    /// Return the signal label used by the selected data.
    pub fn signal(&self) -> &FCoverageSignal {
        match self {
            Self::Values(selection) => selection.signal(),
            Self::SummaryStats(selection) => selection.signal(),
        }
    }

    /// Return the scalar aggregate mode, or an error for summary-stat selections.
    pub fn value_mode(&self) -> OutputLoaderResult<FCoverageValueMode> {
        match self {
            Self::Values(selection) => Ok(selection.value_mode()),
            Self::SummaryStats(_) => Err(OutputLoaderError::message(
                "fcoverage selection contains summary stats",
            )),
        }
    }

    /// Return selected scalar values, or an error for summary-stat selections.
    pub fn values(&self) -> OutputLoaderResult<&[f64]> {
        match self {
            Self::Values(selection) => Ok(selection.values()),
            Self::SummaryStats(_) => Err(OutputLoaderError::message(
                "fcoverage selection contains summary stats",
            )),
        }
    }

    /// Return one selected scalar value, or an error for summary-stat selections.
    ///
    /// Parameters
    /// ----------
    /// - `selected_row_index`:
    ///     Zero-based row index within the selected rows.
    pub fn value(&self, selected_row_index: usize) -> OutputLoaderResult<Option<f64>> {
        Ok(self.values()?.get(selected_row_index).copied())
    }

    /// Return selected summary stats, or an error for scalar-value selections.
    pub fn summary_stats(&self) -> OutputLoaderResult<&[FCoverageSummaryStats]> {
        match self {
            Self::SummaryStats(selection) => Ok(selection.stats()),
            Self::Values(_) => Err(OutputLoaderError::message(
                "fcoverage selection contains scalar values",
            )),
        }
    }

    /// Return one selected summary-stat row, or an error for scalar-value selections.
    ///
    /// Parameters
    /// ----------
    /// - `selected_row_index`:
    ///     Zero-based row index within the selected rows.
    pub fn summary_stat(
        &self,
        selected_row_index: usize,
    ) -> OutputLoaderResult<Option<&FCoverageSummaryStats>> {
        Ok(self.summary_stats()?.get(selected_row_index))
    }
}

/// Selected scalar aggregate fcoverage rows.
#[derive(Debug, Clone, PartialEq)]
pub struct FCoverageValueSelection {
    row_metadata: FCoverageRowMetadata,
    row_indices: Vec<usize>,
    value_mode: FCoverageValueMode,
    signal: FCoverageSignal,
    values: Vec<f64>,
}

impl FCoverageValueSelection {
    /// Return selected row metadata in selection order.
    pub fn row_metadata(&self) -> &FCoverageRowMetadata {
        &self.row_metadata
    }

    /// Return selected window metadata, or an error if the selection is not windowed.
    pub fn window_metadata(&self) -> OutputLoaderResult<&[FCoverageWindowRow]> {
        match &self.row_metadata {
            FCoverageRowMetadata::Windows(windows) => Ok(windows),
            _ => Err(OutputLoaderError::message(
                "fcoverage value selection is not windowed",
            )),
        }
    }

    /// Return selected group metadata, or an error if the selection is not grouped.
    pub fn group_metadata(&self) -> OutputLoaderResult<&[FCoverageGroupRow]> {
        match &self.row_metadata {
            FCoverageRowMetadata::Groups(groups) => Ok(groups),
            _ => Err(OutputLoaderError::message(
                "fcoverage value selection is not grouped",
            )),
        }
    }

    /// Return selected source row indices in selection order.
    pub fn row_indices(&self) -> &[usize] {
        &self.row_indices
    }

    /// Return the number of selected rows.
    pub fn row_count(&self) -> usize {
        self.row_indices.len()
    }

    /// Return whether selected values are averages or totals.
    pub fn value_mode(&self) -> FCoverageValueMode {
        self.value_mode
    }

    /// Return the signal label used in selected value columns.
    pub fn signal(&self) -> &FCoverageSignal {
        &self.signal
    }

    /// Return selected values in requested row order.
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    /// Return one selected value, if `selected_row_index` is in bounds.
    ///
    /// Parameters
    /// ----------
    /// - `selected_row_index`:
    ///     Zero-based row index within the selected scalar values.
    pub fn value(&self, selected_row_index: usize) -> Option<f64> {
        self.values.get(selected_row_index).copied()
    }
}

/// Selected summary-stat fcoverage rows.
#[derive(Debug, Clone, PartialEq)]
pub struct FCoverageSummaryStatsSelection {
    row_metadata: FCoverageRowMetadata,
    row_indices: Vec<usize>,
    signal: FCoverageSignal,
    stats: Vec<FCoverageSummaryStats>,
}

impl FCoverageSummaryStatsSelection {
    /// Return selected row metadata in selection order.
    pub fn row_metadata(&self) -> &FCoverageRowMetadata {
        &self.row_metadata
    }

    /// Return selected window metadata, or an error if the selection is not windowed.
    pub fn window_metadata(&self) -> OutputLoaderResult<&[FCoverageWindowRow]> {
        match &self.row_metadata {
            FCoverageRowMetadata::Windows(windows) => Ok(windows),
            _ => Err(OutputLoaderError::message(
                "fcoverage summary-stat selection is not windowed",
            )),
        }
    }

    /// Return selected group metadata, or an error if the selection is not grouped.
    pub fn group_metadata(&self) -> OutputLoaderResult<&[FCoverageGroupRow]> {
        match &self.row_metadata {
            FCoverageRowMetadata::Groups(groups) => Ok(groups),
            _ => Err(OutputLoaderError::message(
                "fcoverage summary-stat selection is not grouped",
            )),
        }
    }

    /// Return selected source row indices in selection order.
    pub fn row_indices(&self) -> &[usize] {
        &self.row_indices
    }

    /// Return the number of selected rows.
    pub fn row_count(&self) -> usize {
        self.row_indices.len()
    }

    /// Return the signal label used in selected summary-stat columns.
    pub fn signal(&self) -> &FCoverageSignal {
        &self.signal
    }

    /// Return selected summary statistics in requested row order.
    pub fn stats(&self) -> &[FCoverageSummaryStats] {
        &self.stats
    }

    /// Return one selected summary-stat row, if `selected_row_index` is in bounds.
    ///
    /// Parameters
    /// ----------
    /// - `selected_row_index`:
    ///     Zero-based row index within the selected summary-stat rows.
    pub fn stat(&self, selected_row_index: usize) -> Option<&FCoverageSummaryStats> {
        self.stats.get(selected_row_index)
    }
}

/// Metadata for one windowed `fcoverage` aggregate row.
#[derive(Debug, Clone, PartialEq)]
pub struct FCoverageWindowRow {
    /// Zero-based row index in output-file order.
    pub index: usize,
    /// Chromosome or contig label from the output file.
    pub chrom: String,
    /// Checked half-open genomic interval for this row.
    pub interval: Interval<u64>,
    /// Number of positions in the row span that were blacklisted.
    pub blacklisted_positions: u64,
    /// Number of positions eligible for coverage after masking.
    pub eligible_positions: Option<u64>,
}

impl FCoverageWindowRow {
    /// Return `blacklisted_positions / interval.len()`.
    pub fn blacklisted_fraction(&self) -> f64 {
        self.blacklisted_positions as f64 / self.interval.len() as f64
    }
}

/// Metadata for one grouped `fcoverage` aggregate row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FCoverageGroupRow {
    /// Zero-based row index in output-file order.
    pub index: usize,
    /// Group index written by `cfdna fcoverage`.
    pub group_idx: u64,
    /// Group name loaded from a group-index file, when one was provided.
    pub name: Option<String>,
    /// Total represented group span in positions.
    pub span_positions: u64,
    /// Number of represented positions that were blacklisted.
    pub blacklisted_positions: u64,
    /// Number of positions eligible for coverage after masking.
    pub eligible_positions: u64,
}

impl FCoverageGroupRow {
    /// Return `blacklisted_positions / span_positions`, or `NaN` for an empty span.
    pub fn blacklisted_fraction(&self) -> f64 {
        if self.span_positions == 0 {
            f64::NAN
        } else {
            self.blacklisted_positions as f64 / self.span_positions as f64
        }
    }
}

/// Raw and derived summary statistics for one `fcoverage` output row.
#[derive(Debug, Clone, PartialEq)]
pub struct FCoverageSummaryStats {
    /// Positions with non-zero coverage or fragment mass.
    pub nonzero_positions: u64,
    /// Fraction of eligible positions with non-zero signal.
    pub covered_fraction: f64,
    /// Sum of the signal over eligible positions.
    pub total: f64,
    /// Sum of squared signal over eligible positions.
    pub total_squared: f64,
    /// Mean signal over eligible positions.
    pub average: f64,
    /// Signal variance over eligible positions.
    pub variance: f64,
    /// Signal standard deviation over eligible positions.
    pub sd: f64,
    /// Signal coefficient of variation.
    pub coefficient_of_variation: FCoverageCoefficientOfVariation,
}

/// Coefficient-of-variation value from a summary-stat row.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FCoverageCoefficientOfVariation {
    /// Ordinary numeric value, including `NaN`.
    Value(f64),
    /// Writer-side display cap, for example `>1e6`.
    GreaterThan(f64),
}

/// Parser for one non-positional `cfdna fcoverage` TSV output.
///
/// The parser owns the input path, detects the table schema from the header,
/// parses each data line into row metadata and data values, and then collects
/// them into the public `FCoverageOutput`.
struct FCoverageParser {
    path: PathBuf,
    group_index_path: Option<PathBuf>,
}

impl FCoverageParser {
    /// Store the input path until `load()` opens it.
    fn new(path: &Path, group_index_path: Option<&Path>) -> Self {
        Self {
            path: path.to_path_buf(),
            group_index_path: group_index_path.map(Path::to_path_buf),
        }
    }

    /// Read the TSV header, parse all data rows, and build an `FCoverageOutput`.
    fn load(&self) -> Result<FCoverageOutput> {
        let mut reader = open_text_reader(&self.path)
            .with_context(|| format!("open fcoverage output {}", self.path.display()))?;
        let mut header_line = String::new();
        ensure!(
            reader
                .read_line(&mut header_line)
                .with_context(|| format!("read header from {}", self.path.display()))?
                > 0,
            "fcoverage output {} is empty; header required",
            self.path.display()
        );

        let header = split_header(&header_line);
        let schema = FCoverageSchema::from_header(&self.path, &header)?;
        let group_names_by_idx = self
            .group_index_path
            .as_deref()
            .map(read_group_index)
            .transpose()?;
        ensure!(
            group_names_by_idx.is_none() || schema.row_mode() == FCoverageRowMode::Groups,
            "fcoverage group-index file can only be used with grouped output: {}",
            self.path.display()
        );

        let mut rows = Vec::new();
        for (line_offset, line_result) in reader.lines().enumerate() {
            let line_number = line_offset + 2;
            let line = line_result
                .with_context(|| format!("read line {line_number} from {}", self.path.display()))?;
            if line.is_empty() {
                bail!(
                    "fcoverage output {} line {line_number} is empty",
                    self.path.display()
                );
            }
            rows.push(schema.parse_row(
                &self.path,
                line_number,
                &line,
                group_names_by_idx.as_ref(),
            )?);
        }

        if let Some(group_names_by_idx) = &group_names_by_idx {
            let group_index_path = self
                .group_index_path
                .as_deref()
                .expect("group-index path exists when group names were loaded");
            ensure_group_index_matches_rows(
                &self.path,
                group_index_path,
                group_names_by_idx,
                &rows,
            )?;
        }

        schema.finish(rows, parse_filename_metadata(&self.path))
    }
}

/// Parsed header schema for a non-positional fcoverage aggregate table.
///
/// The schema records both row mode and value mode so data lines can be parsed
/// without re-inspecting header strings.
#[derive(Debug, Clone)]
enum FCoverageSchema {
    Value {
        row_mode: FCoverageRowMode,
        value_mode: FCoverageValueMode,
        signal: FCoverageSignal,
    },
    SummaryStats {
        row_mode: FCoverageRowMode,
        signal: FCoverageSignal,
    },
}

impl FCoverageSchema {
    /// Detect row mode, value mode, and signal label from the TSV header.
    fn from_header(path: &Path, header: &[&str]) -> Result<Self> {
        match header {
            [
                "chromosome",
                "start",
                "end",
                value_header,
                "blacklisted_positions",
            ] => {
                let (value_mode, signal) = parse_value_header(path, value_header)?;
                Ok(Self::Value {
                    row_mode: FCoverageRowMode::Windows,
                    value_mode,
                    signal,
                })
            }
            [
                "group_idx",
                "span_positions",
                "blacklisted_positions",
                "eligible_positions",
                value_header,
            ] => {
                let (value_mode, signal) = parse_value_header(path, value_header)?;
                Ok(Self::Value {
                    row_mode: FCoverageRowMode::Groups,
                    value_mode,
                    signal,
                })
            }
            [
                "chromosome",
                "start",
                "end",
                "span_positions",
                "blacklisted_positions",
                "eligible_positions",
                "nonzero_positions",
                "covered_fraction",
                total_header,
                total_squared_header,
                average_header,
                variance_header,
                sd_header,
                coefficient_of_variation_header,
            ] => {
                let signal = parse_summary_signal(
                    path,
                    total_header,
                    total_squared_header,
                    average_header,
                    variance_header,
                    sd_header,
                    coefficient_of_variation_header,
                )?;
                Ok(Self::SummaryStats {
                    row_mode: FCoverageRowMode::Windows,
                    signal,
                })
            }
            [
                "group_idx",
                "span_positions",
                "blacklisted_positions",
                "eligible_positions",
                "nonzero_positions",
                "covered_fraction",
                total_header,
                total_squared_header,
                average_header,
                variance_header,
                sd_header,
                coefficient_of_variation_header,
            ] => {
                let signal = parse_summary_signal(
                    path,
                    total_header,
                    total_squared_header,
                    average_header,
                    variance_header,
                    sd_header,
                    coefficient_of_variation_header,
                )?;
                Ok(Self::SummaryStats {
                    row_mode: FCoverageRowMode::Groups,
                    signal,
                })
            }
            _ => bail!(
                "fcoverage output {} has unsupported aggregate header",
                path.display()
            ),
        }
    }

    /// Return whether this schema describes windowed or grouped rows.
    fn row_mode(&self) -> FCoverageRowMode {
        match self {
            Self::Value { row_mode, .. } | Self::SummaryStats { row_mode, .. } => *row_mode,
        }
    }

    /// Parse one TSV data row according to the header-derived schema.
    fn parse_row(
        &self,
        path: &Path,
        line_number: usize,
        line: &str,
        group_names_by_idx: Option<&FxHashMap<u64, String>>,
    ) -> Result<ParsedRow> {
        let fields = line.split('\t').collect::<Vec<_>>();
        match self {
            Self::Value {
                row_mode: FCoverageRowMode::Windows,
                value_mode,
                ..
            } => {
                ensure_column_count(path, line_number, &fields, 5)?;
                let row_index = line_number - 2;
                let window = parse_window_row(path, line_number, row_index, &fields, None)?;
                let value = parse_scalar_value_field(
                    path,
                    line_number,
                    *value_mode,
                    fields[3],
                    Some(window.interval.len() - window.blacklisted_positions),
                )?;
                Ok(ParsedRow::Value {
                    row_metadata: ParsedRowMetadata::Window(window),
                    value,
                })
            }
            Self::Value {
                row_mode: FCoverageRowMode::Groups,
                value_mode,
                ..
            } => {
                ensure_column_count(path, line_number, &fields, 5)?;
                let row_index = line_number - 2;
                let group =
                    parse_group_row(path, line_number, row_index, &fields, group_names_by_idx)?;
                let value = parse_scalar_value_field(
                    path,
                    line_number,
                    *value_mode,
                    fields[4],
                    Some(group.eligible_positions),
                )?;
                Ok(ParsedRow::Value {
                    row_metadata: ParsedRowMetadata::Group(group),
                    value,
                })
            }
            Self::SummaryStats {
                row_mode: FCoverageRowMode::Windows,
                ..
            } => {
                ensure_column_count(path, line_number, &fields, 14)?;
                let row_index = line_number - 2;
                let span_positions =
                    parse_u64_field(path, line_number, "span_positions", fields[3])?;
                let window =
                    parse_window_row(path, line_number, row_index, &fields, Some(span_positions))?;
                let stats = parse_window_summary_stats(
                    path,
                    line_number,
                    &fields,
                    window
                        .eligible_positions
                        .expect("summary-stat window rows always store eligible_positions"),
                )?;
                Ok(ParsedRow::SummaryStats {
                    row_metadata: ParsedRowMetadata::Window(window),
                    stats,
                })
            }
            Self::SummaryStats {
                row_mode: FCoverageRowMode::Groups,
                ..
            } => {
                ensure_column_count(path, line_number, &fields, 12)?;
                let row_index = line_number - 2;
                let group =
                    parse_group_row(path, line_number, row_index, &fields, group_names_by_idx)?;
                let stats = parse_group_summary_stats(
                    path,
                    line_number,
                    &fields,
                    group.eligible_positions,
                )?;
                Ok(ParsedRow::SummaryStats {
                    row_metadata: ParsedRowMetadata::Group(group),
                    stats,
                })
            }
        }
    }

    /// Convert parsed TSV rows into the public fcoverage output object.
    ///
    /// `parse_row()` stores each line as a parsed value row or summary-stat row
    /// with row metadata. This method uses the header-derived schema variant to
    /// collect only the matching row kind, keeps row metadata in file order, and
    /// attaches the signal label and value mode needed by downstream callers.
    ///
    /// Parameters
    /// ----------
    /// - `rows`:
    ///     Parsed data rows in file order.
    ///
    /// Returns
    /// -------
    /// - `FCoverageOutput`:
    ///     Final loaded output with row metadata, signal label, and either
    ///     scalar values or summary statistics.
    fn finish(
        self,
        rows: Vec<ParsedRow>,
        filename_metadata: FCoverageFilenameMetadata,
    ) -> Result<FCoverageOutput> {
        match self {
            Self::Value {
                value_mode, signal, ..
            } => {
                let (row_metadata, values) = collect_value_rows(rows)?;
                Ok(FCoverageOutput {
                    row_metadata,
                    signal,
                    filename_metadata,
                    data: FCoverageData::Values { value_mode, values },
                })
            }
            Self::SummaryStats { signal, .. } => {
                let (row_metadata, stats) = collect_summary_rows(rows)?;
                Ok(FCoverageOutput {
                    row_metadata,
                    signal,
                    filename_metadata,
                    data: FCoverageData::SummaryStats(stats),
                })
            }
        }
    }
}

/// One parsed data line before rows are collected into the public output.
#[derive(Debug, Clone, PartialEq)]
enum ParsedRow {
    Value {
        row_metadata: ParsedRowMetadata,
        value: f64,
    },
    SummaryStats {
        row_metadata: ParsedRowMetadata,
        stats: FCoverageSummaryStats,
    },
}

/// Row identity parsed from one fcoverage data line.
#[derive(Debug, Clone, PartialEq)]
enum ParsedRowMetadata {
    Window(FCoverageWindowRow),
    Group(FCoverageGroupRow),
}

/// Require every group-index sidecar row to correspond to a grouped TSV row.
fn ensure_group_index_matches_rows(
    path: &Path,
    group_index_path: &Path,
    group_names_by_idx: &FxHashMap<u64, String>,
    rows: &[ParsedRow],
) -> Result<()> {
    let mut seen_group_indices = FxHashSet::default();
    for row in rows {
        match row {
            ParsedRow::Value {
                row_metadata: ParsedRowMetadata::Group(group),
                ..
            }
            | ParsedRow::SummaryStats {
                row_metadata: ParsedRowMetadata::Group(group),
                ..
            } => {
                seen_group_indices.insert(group.group_idx);
            }
            ParsedRow::Value {
                row_metadata: ParsedRowMetadata::Window(_),
                ..
            }
            | ParsedRow::SummaryStats {
                row_metadata: ParsedRowMetadata::Window(_),
                ..
            } => {}
        }
    }

    let mut missing_group_indices = group_names_by_idx
        .keys()
        .copied()
        .filter(|group_idx| !seen_group_indices.contains(group_idx))
        .collect::<Vec<_>>();
    missing_group_indices.sort_unstable();
    if let Some(missing_group_idx) = missing_group_indices.first() {
        bail!(
            "fcoverage group-index file {} contains group_idx {} with no matching row in {}",
            group_index_path.display(),
            missing_group_idx,
            path.display()
        );
    }
    Ok(())
}

/// Split parsed scalar-value rows into row metadata and value vectors.
fn collect_value_rows(rows: Vec<ParsedRow>) -> Result<(FCoverageRowMetadata, Vec<f64>)> {
    let mut row_metadata_entries = Vec::new();
    let mut values = Vec::new();
    for row in rows {
        let ParsedRow::Value {
            row_metadata,
            value,
        } = row
        else {
            bail!("internal fcoverage loader row-mode mismatch");
        };
        row_metadata_entries.push(row_metadata);
        values.push(value);
    }
    Ok((collect_row_metadata(row_metadata_entries)?, values))
}

/// Split parsed summary-stat rows into row metadata and summary-stat vectors.
fn collect_summary_rows(
    rows: Vec<ParsedRow>,
) -> Result<(FCoverageRowMetadata, Vec<FCoverageSummaryStats>)> {
    let mut row_metadata_entries = Vec::new();
    let mut stats = Vec::new();
    for row in rows {
        let ParsedRow::SummaryStats {
            row_metadata,
            stats: row_stats,
        } = row
        else {
            bail!("internal fcoverage loader row-mode mismatch");
        };
        row_metadata_entries.push(row_metadata);
        stats.push(row_stats);
    }
    Ok((collect_row_metadata(row_metadata_entries)?, stats))
}

/// Convert parsed row metadata into the windowed or grouped public variant.
fn collect_row_metadata(rows: Vec<ParsedRowMetadata>) -> Result<FCoverageRowMetadata> {
    let Some(first_row) = rows.first() else {
        bail!("fcoverage output has no data rows");
    };
    match first_row {
        ParsedRowMetadata::Window(_) => rows
            .into_iter()
            .map(|row| match row {
                ParsedRowMetadata::Window(window) => Ok(window),
                ParsedRowMetadata::Group(_) => bail!("internal fcoverage loader row-mode mismatch"),
            })
            .collect::<Result<Vec<_>>>()
            .map(FCoverageRowMetadata::Windows),
        ParsedRowMetadata::Group(_) => rows
            .into_iter()
            .try_fold(
                (Vec::new(), FxHashSet::default()),
                |(mut groups, mut seen_group_indices), row| match row {
                    ParsedRowMetadata::Group(group) => {
                        ensure!(
                            seen_group_indices.insert(group.group_idx),
                            "fcoverage output has duplicate group_idx {}",
                            group.group_idx
                        );
                        groups.push(group);
                        Ok((groups, seen_group_indices))
                    }
                    ParsedRowMetadata::Window(_) => {
                        bail!("internal fcoverage loader row-mode mismatch")
                    }
                },
            )
            .map(|(groups, _)| FCoverageRowMetadata::Groups(groups)),
    }
}

/// Parse one window metadata row and validate span-related fields.
fn parse_window_row(
    path: &Path,
    line_number: usize,
    row_index: usize,
    fields: &[&str],
    span_positions: Option<u64>,
) -> Result<FCoverageWindowRow> {
    // Build the checked genomic interval first, because later position counts
    // are validated against its length
    let start = parse_u64_field(path, line_number, "start", fields[1])?;
    let end = parse_u64_field(path, line_number, "end", fields[2])?;
    let interval = Interval::new(start, end).map_err(|error| {
        anyhow::anyhow!(
            "fcoverage output {} line {line_number} has invalid window interval: {error}",
            path.display()
        )
    })?;
    if let Some(span_positions) = span_positions {
        ensure!(
            interval.len() == span_positions,
            "fcoverage output {} line {line_number} has span_positions {} but interval length {}",
            path.display(),
            span_positions,
            interval.len()
        );
    }
    // Blacklisted positions are present for all non-positional fcoverage rows
    let blacklisted_positions =
        parse_u64_field(path, line_number, "blacklisted_positions", fields[4])?;
    ensure!(
        blacklisted_positions <= interval.len(),
        "fcoverage output {} line {line_number} has blacklisted_positions {} greater than span {}",
        path.display(),
        blacklisted_positions,
        interval.len()
    );
    // Summary-stat outputs also carry eligible_positions, while average and
    // total outputs do not
    let eligible_positions = if span_positions.is_some() {
        let eligible_positions =
            parse_u64_field(path, line_number, "eligible_positions", fields[5])?;
        ensure!(
            eligible_positions <= interval.len(),
            "fcoverage output {} line {line_number} has eligible_positions {} greater than span {}",
            path.display(),
            eligible_positions,
            interval.len()
        );
        let expected_eligible_positions = interval.len() - blacklisted_positions;
        ensure!(
            eligible_positions == expected_eligible_positions,
            "fcoverage output {} line {line_number} has eligible_positions {} but span minus blacklisted_positions is {}",
            path.display(),
            eligible_positions,
            expected_eligible_positions
        );
        Some(eligible_positions)
    } else {
        None
    };
    Ok(FCoverageWindowRow {
        index: row_index,
        chrom: fields[0].to_string(),
        interval,
        blacklisted_positions,
        eligible_positions,
    })
}

/// Parse one grouped metadata row and validate position counts.
fn parse_group_row(
    path: &Path,
    line_number: usize,
    row_index: usize,
    fields: &[&str],
    group_names_by_idx: Option<&FxHashMap<u64, String>>,
) -> Result<FCoverageGroupRow> {
    // Grouped rows summarize all positions assigned to one grouped-BED label
    let group_idx = parse_u64_field(path, line_number, "group_idx", fields[0])?;
    let span_positions = parse_u64_field(path, line_number, "span_positions", fields[1])?;
    let blacklisted_positions =
        parse_u64_field(path, line_number, "blacklisted_positions", fields[2])?;
    let eligible_positions = parse_u64_field(path, line_number, "eligible_positions", fields[3])?;
    ensure!(
        blacklisted_positions <= span_positions,
        "fcoverage output {} line {line_number} has blacklisted_positions {} greater than span_positions {}",
        path.display(),
        blacklisted_positions,
        span_positions
    );
    ensure!(
        eligible_positions <= span_positions,
        "fcoverage output {} line {line_number} has eligible_positions {} greater than span_positions {}",
        path.display(),
        eligible_positions,
        span_positions
    );
    let expected_eligible_positions = span_positions - blacklisted_positions;
    ensure!(
        eligible_positions == expected_eligible_positions,
        "fcoverage output {} line {line_number} has eligible_positions {} but span_positions minus blacklisted_positions is {}",
        path.display(),
        eligible_positions,
        expected_eligible_positions
    );
    let name = group_names_by_idx
        .map(|names_by_idx| {
            names_by_idx.get(&group_idx).cloned().with_context(|| {
                format!("fcoverage group-index file has no group_name for group_idx {group_idx}")
            })
        })
        .transpose()?;
    Ok(FCoverageGroupRow {
        index: row_index,
        group_idx,
        name,
        span_positions,
        blacklisted_positions,
        eligible_positions,
    })
}

/// Parse summary statistics for a windowed row.
fn parse_window_summary_stats(
    path: &Path,
    line_number: usize,
    fields: &[&str],
    eligible_positions: u64,
) -> Result<FCoverageSummaryStats> {
    parse_summary_stats_from_fields(path, line_number, fields, 6, eligible_positions)
}

/// Parse summary statistics for a grouped row.
fn parse_group_summary_stats(
    path: &Path,
    line_number: usize,
    fields: &[&str],
    eligible_positions: u64,
) -> Result<FCoverageSummaryStats> {
    parse_summary_stats_from_fields(path, line_number, fields, 4, eligible_positions)
}

/// Parse summary-stat value columns once the row-metadata prefix length is known.
fn parse_summary_stats_from_fields(
    path: &Path,
    line_number: usize,
    fields: &[&str],
    nonzero_positions_column: usize,
    eligible_positions: u64,
) -> Result<FCoverageSummaryStats> {
    let covered_fraction_column = nonzero_positions_column + 1;
    let nonzero_positions = parse_u64_field(
        path,
        line_number,
        "nonzero_positions",
        fields[nonzero_positions_column],
    )?;
    ensure!(
        nonzero_positions <= eligible_positions,
        "fcoverage output {} line {line_number} has nonzero_positions {} greater than eligible_positions {}",
        path.display(),
        nonzero_positions,
        eligible_positions
    );
    let covered_fraction = parse_fraction_or_zero_support_nan(
        path,
        line_number,
        "covered_fraction",
        fields[covered_fraction_column],
        eligible_positions,
    )?;
    let allow_zero_support_nan = eligible_positions == 0;
    let total = parse_non_negative_summary_value(
        path,
        line_number,
        "total",
        fields[covered_fraction_column + 1],
        allow_zero_support_nan,
    )?;
    let total_squared = parse_non_negative_summary_value(
        path,
        line_number,
        "total_squared",
        fields[covered_fraction_column + 2],
        allow_zero_support_nan,
    )?;
    let average = parse_non_negative_summary_value(
        path,
        line_number,
        "average",
        fields[covered_fraction_column + 3],
        allow_zero_support_nan,
    )?;
    let variance = parse_non_negative_summary_value(
        path,
        line_number,
        "variance",
        fields[covered_fraction_column + 4],
        allow_zero_support_nan,
    )?;
    let sd = parse_non_negative_summary_value(
        path,
        line_number,
        "sd",
        fields[covered_fraction_column + 5],
        allow_zero_support_nan,
    )?;
    let coefficient_of_variation =
        parse_coefficient_of_variation(path, line_number, fields[covered_fraction_column + 6])?;
    validate_coefficient_of_variation(path, line_number, coefficient_of_variation)?;
    Ok(FCoverageSummaryStats {
        nonzero_positions,
        covered_fraction,
        total,
        total_squared,
        average,
        variance,
        sd,
        coefficient_of_variation,
    })
}

/// Parse one scalar value header into value mode and signal label.
fn parse_value_header(path: &Path, header: &str) -> Result<(FCoverageValueMode, FCoverageSignal)> {
    if let Some(signal) = header.strip_prefix("average_") {
        ensure!(
            !signal.is_empty(),
            "fcoverage output {} has empty average signal label",
            path.display()
        );
        Ok((FCoverageValueMode::Average, FCoverageSignal::new(signal)))
    } else if let Some(signal) = header.strip_prefix("total_") {
        ensure!(
            !signal.is_empty(),
            "fcoverage output {} has empty total signal label",
            path.display()
        );
        Ok((FCoverageValueMode::Total, FCoverageSignal::new(signal)))
    } else {
        bail!(
            "fcoverage output {} has unsupported aggregate value column '{}'",
            path.display(),
            header
        );
    }
}

/// Parse and cross-check the shared signal label in summary-stat headers.
fn parse_summary_signal(
    path: &Path,
    total_header: &str,
    total_squared_header: &str,
    average_header: &str,
    variance_header: &str,
    sd_header: &str,
    coefficient_of_variation_header: &str,
) -> Result<FCoverageSignal> {
    let signal = strip_required_prefix(path, total_header, "total_")?;
    let expected = [
        (total_squared_header, "total_squared_"),
        (average_header, "average_"),
        (variance_header, "variance_"),
        (sd_header, "sd_"),
        (coefficient_of_variation_header, "coefficient_of_variation_"),
    ];
    for (header, prefix) in expected {
        let observed_signal = strip_required_prefix(path, header, prefix)?;
        ensure!(
            observed_signal == signal,
            "fcoverage output {} has inconsistent summary-stat signal labels: '{}' and '{}'",
            path.display(),
            signal,
            observed_signal
        );
    }
    Ok(FCoverageSignal::new(signal))
}

/// Strip a required column-name prefix and return the remaining signal label.
fn strip_required_prefix<'a>(path: &Path, header: &'a str, prefix: &str) -> Result<&'a str> {
    let signal = header.strip_prefix(prefix).with_context(|| {
        format!(
            "fcoverage output {} summary header '{}' must start with '{}'",
            path.display(),
            header,
            prefix
        )
    })?;
    ensure!(
        !signal.is_empty(),
        "fcoverage output {} summary header '{}' has empty signal label",
        path.display(),
        header
    );
    Ok(signal)
}

/// Parse a numeric coefficient of variation or writer-side display cap.
fn parse_coefficient_of_variation(
    path: &Path,
    line_number: usize,
    value: &str,
) -> Result<FCoverageCoefficientOfVariation> {
    if let Some(threshold) = value.strip_prefix('>') {
        let threshold = threshold.parse::<f64>().with_context(|| {
            format!(
                "fcoverage output {} line {line_number} has invalid coefficient_of_variation threshold '{}'",
                path.display(),
                value
            )
        })?;
        ensure!(
            threshold.is_finite() && threshold > 0.0,
            "fcoverage output {} line {line_number} has invalid coefficient_of_variation threshold '{}'",
            path.display(),
            value
        );
        Ok(FCoverageCoefficientOfVariation::GreaterThan(threshold))
    } else {
        Ok(FCoverageCoefficientOfVariation::Value(parse_f64_field(
            path,
            line_number,
            "coefficient_of_variation",
            value,
        )?))
    }
}

/// Parse one scalar average or total value with field-specific numeric rules.
fn parse_scalar_value_field(
    path: &Path,
    line_number: usize,
    value_mode: FCoverageValueMode,
    value: &str,
    eligible_positions: Option<u64>,
) -> Result<f64> {
    let parsed = parse_f64_field(path, line_number, "value", value)?;
    match value_mode {
        FCoverageValueMode::Average => {
            let allow_nan = eligible_positions == Some(0);
            validate_non_negative_f64(path, line_number, "average value", parsed, allow_nan)?;
        }
        FCoverageValueMode::Total => {
            validate_non_negative_f64(path, line_number, "total value", parsed, false)?;
        }
    }
    Ok(parsed)
}

/// Parse a fraction while allowing `NaN` only for zero-support rows.
fn parse_fraction_or_zero_support_nan(
    path: &Path,
    line_number: usize,
    field_name: &str,
    value: &str,
    eligible_positions: u64,
) -> Result<f64> {
    let fraction = parse_f64_field(path, line_number, field_name, value)?;
    if fraction.is_nan() && eligible_positions == 0 {
        return Ok(fraction);
    }
    ensure!(
        fraction.is_finite() && (0.0..=1.0).contains(&fraction),
        "fcoverage output {} line {line_number} has {field_name} outside [0, 1]: {fraction}",
        path.display()
    );
    Ok(fraction)
}

/// Parse a non-negative summary value while preserving zero-support `NaN`.
fn parse_non_negative_summary_value(
    path: &Path,
    line_number: usize,
    field_name: &str,
    value: &str,
    allow_nan: bool,
) -> Result<f64> {
    let parsed = parse_f64_field(path, line_number, field_name, value)?;
    validate_non_negative_f64(path, line_number, field_name, parsed, allow_nan)?;
    Ok(parsed)
}

/// Reject negative finite values and infinities in fcoverage numeric fields.
fn validate_non_negative_f64(
    path: &Path,
    line_number: usize,
    field_name: &str,
    value: f64,
    allow_nan: bool,
) -> Result<()> {
    if allow_nan && value.is_nan() {
        return Ok(());
    }
    ensure!(
        value.is_finite() && value >= 0.0,
        "fcoverage output {} line {line_number} has {field_name} outside finite and non-negative range: {value}",
        path.display()
    );
    Ok(())
}

/// Validate an already parsed coefficient of variation value.
fn validate_coefficient_of_variation(
    path: &Path,
    line_number: usize,
    value: FCoverageCoefficientOfVariation,
) -> Result<()> {
    match value {
        FCoverageCoefficientOfVariation::Value(value) if value.is_nan() => Ok(()),
        FCoverageCoefficientOfVariation::Value(value) => {
            validate_non_negative_f64(path, line_number, "coefficient_of_variation", value, false)
        }
        FCoverageCoefficientOfVariation::GreaterThan(_) => Ok(()),
    }
}

/// Reject command-generated positional fcoverage output suffixes.
fn ensure_non_positional_path(path: &Path) -> Result<()> {
    let file_name = path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or_default();
    let positional_suffixes = [
        ".fcoverage.per_position.bedgraph",
        ".fcoverage.per_position.bedgraph.gz",
        ".fcoverage.per_position.bedgraph.zst",
        ".fcoverage.per_position_per_window.tsv",
        ".fcoverage.per_position_per_window.tsv.gz",
        ".fcoverage.per_position_per_window.tsv.zst",
    ];
    ensure!(
        !positional_suffixes
            .iter()
            .any(|suffix| file_name.ends_with(suffix)),
        "positional fcoverage outputs are not supported by load_fcoverage_output: {}",
        path.display()
    );
    Ok(())
}

/// Parse command-mode hints from canonical fcoverage output filename parts.
fn parse_filename_metadata(path: &Path) -> FCoverageFilenameMetadata {
    let file_name = path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or_default();
    let aggregation_basis = if contains_filename_part(file_name, "average_on_unique_bases")
        || contains_filename_part(file_name, "total_on_unique_bases")
        || contains_filename_part(file_name, "summary_stats_on_unique_bases")
    {
        FCoverageAggregationBasis::UniqueBases
    } else if contains_filename_part(file_name, "average")
        || contains_filename_part(file_name, "total")
        || contains_filename_part(file_name, "summary_stats")
    {
        FCoverageAggregationBasis::Ordinary
    } else {
        FCoverageAggregationBasis::Unknown
    };
    let length_normalization = if file_name.contains(".length_normalized.restored_mean.fcoverage.")
    {
        FCoverageLengthNormalization::RestoredMean
    } else if file_name.contains(".length_normalized.fcoverage.") {
        FCoverageLengthNormalization::UnitMass
    } else if file_name.contains(".fcoverage.") {
        FCoverageLengthNormalization::Off
    } else {
        FCoverageLengthNormalization::Unknown
    };
    FCoverageFilenameMetadata {
        aggregation_basis,
        length_normalization,
    }
}

/// Return whether `part` appears as a dot-delimited filename component.
fn contains_filename_part(file_name: &str, part: &str) -> bool {
    let marker = format!(".{part}.");
    file_name.contains(&marker)
}

/// Read a grouped fcoverage group-index file with `group_idx` and `group_name` columns.
fn read_group_index(path: &Path) -> Result<FxHashMap<u64, String>> {
    let mut reader = open_text_reader(path)
        .with_context(|| format!("open fcoverage group-index file {}", path.display()))?;
    let mut header_line = String::new();
    ensure!(
        reader
            .read_line(&mut header_line)
            .with_context(|| format!("read header from {}", path.display()))?
            > 0,
        "fcoverage group-index file {} is empty; header required",
        path.display()
    );
    let header = split_header(&header_line);
    ensure!(
        header == ["group_idx", "group_name"],
        "fcoverage group-index file {} has unsupported header",
        path.display()
    );

    let mut names_by_idx = FxHashMap::default();
    let mut seen_names = FxHashSet::default();
    for (line_offset, line_result) in reader.lines().enumerate() {
        let line_number = line_offset + 2;
        let line = line_result
            .with_context(|| format!("read line {line_number} from {}", path.display()))?;
        if line.is_empty() {
            bail!(
                "fcoverage group-index file {} line {line_number} is empty",
                path.display()
            );
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        ensure!(
            fields.len() == 2,
            "fcoverage group-index file {} line {line_number} has {} columns, expected 2",
            path.display(),
            fields.len()
        );
        let group_idx = parse_u64_field(path, line_number, "group_idx", fields[0])?;
        let group_name = fields[1].to_string();
        ensure!(
            !group_name.is_empty(),
            "fcoverage group-index file {} line {line_number} has empty group_name",
            path.display()
        );
        ensure!(
            names_by_idx.insert(group_idx, group_name.clone()).is_none(),
            "fcoverage group-index file {} has duplicate group_idx {}",
            path.display(),
            group_idx
        );
        ensure!(
            seen_names.insert(group_name.clone()),
            "fcoverage group-index file {} has duplicate group_name '{}'",
            path.display(),
            group_name
        );
    }
    Ok(names_by_idx)
}

/// Split a TSV header line into tab-delimited column names.
fn split_header(header_line: &str) -> Vec<&str> {
    header_line
        .trim_end_matches(['\r', '\n'])
        .split('\t')
        .collect()
}

/// Require one parsed TSV line to have the expected number of columns.
fn ensure_column_count(
    path: &Path,
    line_number: usize,
    fields: &[&str],
    expected: usize,
) -> Result<()> {
    ensure!(
        fields.len() == expected,
        "fcoverage output {} line {line_number} has {} columns, expected {expected}",
        path.display(),
        fields.len()
    );
    Ok(())
}

/// Parse one non-negative integer data field with file and line context.
fn parse_u64_field(path: &Path, line_number: usize, field_name: &str, value: &str) -> Result<u64> {
    value.parse::<u64>().with_context(|| {
        format!(
            "fcoverage output {} line {line_number} has invalid {field_name} '{}'",
            path.display(),
            value
        )
    })
}

/// Parse one floating-point data field with file and line context.
fn parse_f64_field(path: &Path, line_number: usize, field_name: &str, value: &str) -> Result<f64> {
    value.parse::<f64>().with_context(|| {
        format!(
            "fcoverage output {} line {line_number} has invalid {field_name} '{}'",
            path.display(),
            value
        )
    })
}
