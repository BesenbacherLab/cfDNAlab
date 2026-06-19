//! Loader for `cfdna midpoints` profile Zarr outputs.
//!
//! A loaded midpoint output is a profile array with axes `(group, length_bin,
//! position)`. The loader reads and validates the Zarr root metadata and axis
//! metadata arrays when the store is loaded. The `counts` array is read only by
//! `read_all_counts()` or a selector returned by `select()`, so opening a large
//! profile store does not load all profile values.
//!
//! Type overview:
//!
//! ```text
//! load_midpoints_output(path)
//!     -> MidpointsOutput
//!         groups: Vec<MidpointGroupRow>
//!         length_bins: Vec<LengthBin>
//!         position_bins: Vec<MidpointPositionBin>
//!         count_shape: (groups, length_bins, positions)
//!         counts on disk: Zarr /counts array
//!
//! MidpointsOutput::select()
//!     -> MidpointsSelector
//!         -> read()
//!             -> MidpointCountSelection
//!                 groups: Vec<MidpointGroupRow>
//!                 length_bins: Vec<LengthBin>
//!                 position_bins: Vec<MidpointPositionBin>
//!                 counts: DenseArray3<f32>
//! ```
//!
//! Selections return owned group metadata, length bins, position bins, and a
//! `DenseArray3<f32>`. A selection preserves requested group, length-bin, and
//! position order, while the returned array is always stored in row-major
//! `(group, length_bin, position)` order.

use crate::{
    interval::{IndexedInterval, Interval},
    output_loaders::{
        OutputLoaderError, OutputLoaderResult,
        common::{
            DenseArray3, LengthBin, contiguous_index_span, ensure_unique_indices,
            ensure_unique_labels, resolve_row_indices, validate_zarr_public_label,
        },
    },
    shared::{
        constants::{MAX_SUPPORTED_FRAGMENT_LENGTH, MIN_ACGT_BASES_FOR_GC_FRACTION},
        zarr::read_zarr_root_attributes,
    },
};
use anyhow::{Context, Result, bail, ensure};
use fxhash::FxHashMap;
use serde_json::Value;
use std::{
    fs,
    ops::Range,
    path::{Path, PathBuf},
    sync::Arc,
};
use zarrs::{
    array::{Array, ArraySubset, ElementOwned},
    filesystem::FilesystemStore,
};

const MIDPOINT_SCHEMA_VERSION: u64 = 1;

/// One indexed half-open midpoint position bin in output-column order.
pub type MidpointPositionBin = IndexedInterval<u32, usize>;

/// Load a `cfdna midpoints` profile Zarr store.
///
/// The path must point to a midpoint profile Zarr directory. This function
/// validates the schema and reads group, fragment length, and position axes.
/// Count values are read later by `read_all_counts()` or `select().read()`.
/// Outputs with no groups, no length bins, or no position bins are rejected.
///
/// Parameters
/// ----------
/// - `path`:
///     Path to a `cfdna midpoints` Zarr output directory.
///
/// Returns
/// -------
/// - `MidpointsOutput`:
///     Loaded axis metadata and count shape. Count blocks are read from disk by
///     `read_all_counts()` or selection methods.
///
/// ```no_run
/// use cfdnalab::{
///     interval::Interval,
///     output_loaders::load_midpoints_output,
/// };
///
/// let midpoints = load_midpoints_output("sample.midpoint_profiles.zarr")?;
/// let group_index = midpoints.group_index("promoter")?;
/// let selected = midpoints
///     .select()
///     .groups(&[group_index])
///     .length_range(Interval::new(120, 181)?)
///     .read()?;
///
/// for (selected_group_index, group) in selected.groups().iter().enumerate() {
///     for (selected_length_index, length_bin) in selected.length_bins().iter().enumerate() {
///         let profile = selected
///             .profile(selected_group_index, selected_length_index)
///             .expect("selected indices are in bounds");
///         let profile_total = profile.iter().copied().sum::<f32>();
///         println!("{} {:?}: {profile_total}", group.name, length_bin.as_tuple());
///     }
/// }
///
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn load_midpoints_output(path: impl AsRef<Path>) -> OutputLoaderResult<MidpointsOutput> {
    MidpointsParser::new(path.as_ref())
        .load()
        .map_err(Into::into)
}

/// Loaded midpoint metadata and path-backed count access.
#[derive(Debug, Clone, PartialEq)]
pub struct MidpointsOutput {
    path: PathBuf,
    groups: Vec<MidpointGroupRow>,
    group_name_indices: FxHashMap<String, usize>,
    length_bins: Vec<LengthBin>,
    position_bins: Vec<MidpointPositionBin>,
    count_shape: (usize, usize, usize),
}

impl MidpointsOutput {
    /// Return group metadata in the same order as the first count-array dimension.
    pub fn group_metadata(&self) -> &[MidpointGroupRow] {
        &self.groups
    }

    /// Return the number of groups.
    pub fn group_count(&self) -> usize {
        self.groups.len()
    }

    /// Return fragment length bins in the same order as the second count-array dimension.
    pub fn length_bins(&self) -> &[LengthBin] {
        &self.length_bins
    }

    /// Return the number of fragment length bins.
    pub fn length_bin_count(&self) -> usize {
        self.length_bins.len()
    }

    /// Return interval-relative position bins in the same order as the third count-array dimension.
    pub fn position_bins(&self) -> &[MidpointPositionBin] {
        &self.position_bins
    }

    /// Return the number of interval-relative position bins.
    pub fn position_bin_count(&self) -> usize {
        self.position_bins.len()
    }

    /// Return the count array shape as `(groups, length_bins, positions)`.
    pub fn counts_shape(&self) -> (usize, usize, usize) {
        self.count_shape
    }

    /// Return one group row by zero-based group index.
    ///
    /// Parameters
    /// ----------
    /// - `group_index`:
    ///     Zero-based group axis index.
    pub fn group(&self, group_index: usize) -> Option<&MidpointGroupRow> {
        self.groups.get(group_index)
    }

    /// Return the group index for one group name.
    ///
    /// Group names are expected to identify one row. The loader reports an
    /// error if the name is missing or the file contains duplicate group names.
    ///
    /// Parameters
    /// ----------
    /// - `group_name`:
    ///     Group label to resolve to a zero-based group index.
    pub fn group_index(&self, group_name: &str) -> OutputLoaderResult<usize> {
        Ok(self
            .group_name_indices
            .get(group_name)
            .copied()
            .with_context(|| format!("midpoints output has no group named '{group_name}'"))?)
    }

    /// Return whether a group name exists.
    ///
    /// Parameters
    /// ----------
    /// - `group_name`:
    ///     Group label to look up.
    pub fn has_group(&self, group_name: &str) -> bool {
        self.group_name_indices.contains_key(group_name)
    }

    /// Return the index of the bin containing `fragment_length_bp`.
    ///
    /// Parameters
    /// ----------
    /// - `fragment_length_bp`:
    ///     Fragment length in bp to locate on the length-bin axis.
    pub fn length_bin_for_length(&self, fragment_length_bp: u32) -> Option<usize> {
        self.length_bins
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
    ///     Half-open fragment length interval `[start, end)` in bp.
    pub fn length_bins_overlapping_range(
        &self,
        range: Interval<u32>,
    ) -> OutputLoaderResult<Vec<LengthBin>> {
        let selected = self
            .length_bins
            .iter()
            .copied()
            .filter(|bin| bin.interval.intersects(range))
            .collect::<Vec<_>>();
        if selected.is_empty() {
            return Err(OutputLoaderError::message(format!(
                "length range [{}, {}) bp does not overlap any midpoint length bins",
                range.start(),
                range.end()
            )));
        }
        Ok(selected)
    }

    /// Return the index of the bin containing `position_bp`.
    ///
    /// Parameters
    /// ----------
    /// - `position_bp`:
    ///     Interval-relative position in bp to locate on the position axis.
    pub fn position_bin_for_position(&self, position_bp: u32) -> Option<usize> {
        self.position_bins
            .iter()
            .find(|bin| bin.start() <= position_bp && position_bp < bin.end())
            .map(|bin| bin.idx())
    }

    /// Return position bins overlapping a half-open interval-relative range.
    ///
    /// Parameters
    /// ----------
    /// - `range`:
    ///     Half-open interval-relative position range `[start, end)` in bp.
    pub fn position_bins_overlapping_range(
        &self,
        range: Interval<u32>,
    ) -> OutputLoaderResult<Vec<MidpointPositionBin>> {
        let selected = self
            .position_bins
            .iter()
            .copied()
            .filter(|bin| bin.interval.intersects(range))
            .collect::<Vec<_>>();
        if selected.is_empty() {
            return Err(OutputLoaderError::message(format!(
                "position range [{}, {}) bp does not overlap any midpoint position bins",
                range.start(),
                range.end()
            )));
        }
        Ok(selected)
    }

    /// Read the full count array.
    ///
    /// This can be large. Use `select()` when only some groups, length bins, or
    /// positions are needed.
    pub fn read_all_counts(&self) -> OutputLoaderResult<MidpointCountSelection> {
        self.select().read()
    }

    /// Start a count selection.
    ///
    /// A new selector initially selects all groups, all fragment length bins,
    /// and all positions. Add axis constraints before calling `read()`.
    pub fn select(&self) -> MidpointsSelector<'_> {
        MidpointsSelector::new(self)
    }

    /// Return selected groups, length bins, and positions.
    ///
    /// Passing `None` for a selector selects every value on that axis. Use
    /// either `length_bin_indices` or `length_range`, not both. Use either
    /// `position_indices` or `position_range`, not both.
    pub(crate) fn select_counts(
        &self,
        group_indices: Option<&[usize]>,
        length_bin_indices: Option<&[usize]>,
        length_range: Option<Interval<u32>>,
        position_indices: Option<&[usize]>,
        position_range: Option<Interval<u32>>,
    ) -> Result<MidpointCountSelection> {
        self.select_counts_with_group_indices(
            group_indices,
            length_bin_indices,
            length_range,
            position_indices,
            position_range,
        )
    }

    /// Return selected group names, length bins, and positions.
    ///
    /// Passing `None` for `group_names` selects all groups.
    pub(crate) fn select_group_counts_by_name<S: AsRef<str>>(
        &self,
        group_names: Option<&[S]>,
        length_bin_indices: Option<&[usize]>,
        length_range: Option<Interval<u32>>,
        position_indices: Option<&[usize]>,
        position_range: Option<Interval<u32>>,
    ) -> Result<MidpointCountSelection> {
        let group_indices = match group_names {
            Some(group_names) => Some(self.resolve_group_name_indices(group_names)?),
            None => None,
        };
        self.select_counts_with_group_indices(
            group_indices.as_deref(),
            length_bin_indices,
            length_range,
            position_indices,
            position_range,
        )
    }

    /// Select group, length-bin, and position axes after resolving group names.
    fn select_counts_with_group_indices(
        &self,
        group_indices: Option<&[usize]>,
        length_bin_indices: Option<&[usize]>,
        length_range: Option<Interval<u32>>,
        position_indices: Option<&[usize]>,
        position_range: Option<Interval<u32>>,
    ) -> Result<MidpointCountSelection> {
        let group_indices = resolve_row_indices(group_indices, self.groups.len(), "group")?;
        let length_bin_indices =
            self.resolve_length_bin_indices(length_bin_indices, length_range)?;
        let position_indices = self.resolve_position_indices(position_indices, position_range)?;
        ensure_unique_indices(&group_indices, "group")?;
        ensure_unique_indices(&length_bin_indices, "length bin")?;
        ensure_unique_indices(&position_indices, "position")?;

        let groups = group_indices
            .iter()
            .map(|&group_index| {
                self.groups.get(group_index).cloned().with_context(|| {
                    format!(
                        "group index {group_index} is outside 0..{}",
                        self.groups.len()
                    )
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let length_bins =
            selected_axis_values(&self.length_bins, &length_bin_indices, "length bin")?;
        let position_bins =
            selected_axis_values(&self.position_bins, &position_indices, "position")?;
        let counts = self.read_count_array_selection(
            &group_indices,
            &length_bin_indices,
            &position_indices,
        )?;

        Ok(MidpointCountSelection {
            group_indices,
            groups,
            length_bins,
            position_bins,
            counts,
        })
    }

    /// Resolve group labels to source group-axis indices.
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

    /// Resolve optional length-bin selectors to concrete source axis indices.
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

    /// Resolve optional position selectors to concrete source axis indices.
    fn resolve_position_indices(
        &self,
        position_indices: Option<&[usize]>,
        position_range: Option<Interval<u32>>,
    ) -> Result<Vec<usize>> {
        ensure!(
            position_indices.is_none() || position_range.is_none(),
            "use either position_indices or position_range, not both"
        );
        if let Some(position_indices) = position_indices {
            ensure!(
                !position_indices.is_empty(),
                "cannot select zero midpoint positions"
            );
            return Ok(position_indices.to_vec());
        }
        if let Some(position_range) = position_range {
            return Ok(self
                .position_bins_overlapping_range(position_range)?
                .iter()
                .map(|position_bin| position_bin.idx())
                .collect());
        }
        Ok((0..self.position_bins.len()).collect())
    }

    /// Read selected count values from Zarr into selection order.
    fn read_count_array_selection(
        &self,
        group_indices: &[usize],
        length_bin_indices: &[usize],
        position_indices: &[usize],
    ) -> Result<DenseArray3<f32>> {
        let store = Arc::new(
            FilesystemStore::new(&self.path)
                .with_context(|| format!("open midpoint Zarr store {}", self.path.display()))?,
        );
        let counts_array = Array::open(store, "/counts")?;
        let group_span = contiguous_index_span(group_indices);
        let length_span = contiguous_index_span(length_bin_indices);
        let position_span = contiguous_index_span(position_indices);
        // Fully contiguous selections can be read as one Zarr block with no
        // reordering after the read
        if let (Some(group_span), Some(length_span), Some(position_span)) =
            (group_span, length_span, position_span)
        {
            return read_contiguous_count_block(
                &counts_array,
                group_span,
                length_span,
                position_span,
            );
        }
        // Contiguous group and length axes still benefit from one larger block
        // read, followed by position-axis filtering in memory
        if let (Some(group_span), Some(length_span)) = (group_span, length_span) {
            return read_contiguous_group_length_count_block(
                &counts_array,
                group_span,
                length_span,
                position_indices,
            );
        }
        // Non-contiguous group or length selections are read profile by
        // profile so the returned array preserves requested selector order
        read_ordered_count_selection(
            &counts_array,
            group_indices,
            length_bin_indices,
            position_indices,
        )
    }
}

/// Builder for selecting group, fragment length, and position axes from a midpoint output.
///
/// The builder starts with all values selected on every axis. Set at most one
/// selector per axis; for example, use `positions()` or `position_range()`, not
/// both. Conflicting selector calls are reported by `read()`. Count values are
/// read from Zarr only when `read()` is called.
#[derive(Debug, Clone)]
pub struct MidpointsSelector<'a> {
    output: &'a MidpointsOutput,
    groups: MidpointGroupSelector,
    lengths: MidpointLengthSelector,
    positions: MidpointPositionSelector,
    selection_error: Option<String>,
}

impl<'a> MidpointsSelector<'a> {
    /// Start a selector with all groups, length bins, and positions selected.
    fn new(output: &'a MidpointsOutput) -> Self {
        Self {
            output,
            groups: MidpointGroupSelector::All,
            lengths: MidpointLengthSelector::All,
            positions: MidpointPositionSelector::All,
            selection_error: None,
        }
    }

    /// Select groups by zero-based group index.
    ///
    /// Parameters
    /// ----------
    /// - `group_indices`:
    ///     Group axis indices in count-array order. The returned selection
    ///     keeps this order and rejects duplicates.
    pub fn groups(self, group_indices: &[usize]) -> Self {
        self.set_groups(
            MidpointGroupSelector::Indices(group_indices.to_vec()),
            "groups",
        )
    }

    /// Select groups by group name.
    ///
    /// Parameters
    /// ----------
    /// - `group_names`:
    ///     Group labels from the midpoint output metadata. The returned
    ///     selection follows this order and rejects duplicates.
    pub fn groups_by_name<S: AsRef<str>>(self, group_names: &[S]) -> Self {
        self.set_groups(
            MidpointGroupSelector::Names(
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
    ///     Length-bin axis indices in count-array order. The returned selection
    ///     keeps this order and rejects duplicates.
    pub fn length_bins(self, length_bin_indices: &[usize]) -> Self {
        self.set_lengths(
            MidpointLengthSelector::Indices(length_bin_indices.to_vec()),
            "length_bins",
        )
    }

    /// Select fragment length bins overlapping a half-open length range.
    ///
    /// Parameters
    /// ----------
    /// - `range`:
    ///     Half-open fragment length interval `[start, end)` in bp. Every
    ///     length bin that intersects the interval is selected.
    pub fn length_range(self, range: Interval<u32>) -> Self {
        self.set_lengths(MidpointLengthSelector::Range(range), "length_range")
    }

    /// Select positions by zero-based position-bin index.
    ///
    /// Passing an empty slice is an error. Midpoint selections need at least
    /// one position bin.
    ///
    /// Parameters
    /// ----------
    /// - `position_indices`:
    ///     Position-bin axis indices in count-array order. The returned
    ///     selection keeps this order and rejects duplicates.
    pub fn positions(self, position_indices: &[usize]) -> Self {
        self.set_positions(
            MidpointPositionSelector::Indices(position_indices.to_vec()),
            "positions",
        )
    }

    /// Select positions overlapping a half-open interval-relative range.
    ///
    /// Parameters
    /// ----------
    /// - `range`:
    ///     Half-open interval-relative position range `[start, end)` in bp.
    ///     Every position bin that intersects the interval is selected.
    pub fn position_range(self, range: Interval<u32>) -> Self {
        self.set_positions(MidpointPositionSelector::Range(range), "position_range")
    }

    /// Set the group selector or record a group-axis selector conflict.
    fn set_groups(mut self, selector: MidpointGroupSelector, selector_name: &'static str) -> Self {
        if let Some(previous_selector_name) = self.groups.selector_name() {
            self.record_axis_conflict("group", previous_selector_name, selector_name);
        } else {
            self.groups = selector;
        }
        self
    }

    /// Set the length selector or record a fragment length-axis selector conflict.
    fn set_lengths(
        mut self,
        selector: MidpointLengthSelector,
        selector_name: &'static str,
    ) -> Self {
        if let Some(previous_selector_name) = self.lengths.selector_name() {
            self.record_axis_conflict("fragment length", previous_selector_name, selector_name);
        } else {
            self.lengths = selector;
        }
        self
    }

    /// Set the position selector or record a position-axis selector conflict.
    fn set_positions(
        mut self,
        selector: MidpointPositionSelector,
        selector_name: &'static str,
    ) -> Self {
        if let Some(previous_selector_name) = self.positions.selector_name() {
            self.record_axis_conflict("position", previous_selector_name, selector_name);
        } else {
            self.positions = selector;
        }
        self
    }

    /// Store the first selector conflict so `read()` reports it as an error.
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

    /// Read the selected count array with selected axis metadata.
    pub fn read(self) -> OutputLoaderResult<MidpointCountSelection> {
        self.ensure_no_selector_conflict()?;
        let (length_bin_indices, length_range) = match self.lengths {
            MidpointLengthSelector::All => (None, None),
            MidpointLengthSelector::Indices(indices) => (Some(indices), None),
            MidpointLengthSelector::Range(range) => (None, Some(range)),
        };
        let length_bin_indices = length_bin_indices.as_deref();

        let (position_indices, position_range) = match self.positions {
            MidpointPositionSelector::All => (None, None),
            MidpointPositionSelector::Indices(indices) => (Some(indices), None),
            MidpointPositionSelector::Range(range) => (None, Some(range)),
        };
        let position_indices = position_indices.as_deref();

        let selection = match self.groups {
            MidpointGroupSelector::All => self.output.select_counts(
                None,
                length_bin_indices,
                length_range,
                position_indices,
                position_range,
            ),
            MidpointGroupSelector::Indices(indices) => self.output.select_counts(
                Some(indices.as_slice()),
                length_bin_indices,
                length_range,
                position_indices,
                position_range,
            ),
            MidpointGroupSelector::Names(names) => self.output.select_group_counts_by_name(
                Some(names.as_slice()),
                length_bin_indices,
                length_range,
                position_indices,
                position_range,
            ),
        }?;
        Ok(selection)
    }
}

/// Group-axis selector state recorded by `MidpointsSelector`.
#[derive(Debug, Clone)]
enum MidpointGroupSelector {
    /// Select all groups.
    All,
    /// Select groups by group-axis index.
    Indices(Vec<usize>),
    /// Select groups by public group label.
    Names(Vec<String>),
}

impl MidpointGroupSelector {
    /// Return the public selector method name for conflict messages.
    fn selector_name(&self) -> Option<&'static str> {
        match self {
            Self::All => None,
            Self::Indices(_) => Some("groups"),
            Self::Names(_) => Some("groups_by_name"),
        }
    }
}

/// Fragment length-axis selector state recorded by `MidpointsSelector`.
#[derive(Debug, Clone)]
enum MidpointLengthSelector {
    /// Select all fragment length bins.
    All,
    /// Select fragment length bins by index.
    Indices(Vec<usize>),
    /// Select fragment length bins overlapping a half-open length range.
    Range(Interval<u32>),
}

impl MidpointLengthSelector {
    /// Return the public selector method name for conflict messages.
    fn selector_name(&self) -> Option<&'static str> {
        match self {
            Self::All => None,
            Self::Indices(_) => Some("length_bins"),
            Self::Range(_) => Some("length_range"),
        }
    }
}

/// Position-axis selector state recorded by `MidpointsSelector`.
#[derive(Debug, Clone)]
enum MidpointPositionSelector {
    /// Select all position bins.
    All,
    /// Select position bins by index.
    Indices(Vec<usize>),
    /// Select position bins overlapping a half-open position range.
    Range(Interval<u32>),
}

impl MidpointPositionSelector {
    /// Return the public selector method name for conflict messages.
    fn selector_name(&self) -> Option<&'static str> {
        match self {
            Self::All => None,
            Self::Indices(_) => Some("positions"),
            Self::Range(_) => Some("position_range"),
        }
    }
}

/// Metadata for one midpoint profile group.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MidpointGroupRow {
    /// Zero-based index into the first dimension of the midpoint count array.
    pub index: usize,
    /// Public group name from the grouped input intervals.
    pub name: String,
    /// Number of profile-eligible intervals retained for this group.
    pub eligible_intervals: u64,
}

/// Selected midpoint profile counts and axis metadata.
#[derive(Debug, Clone, PartialEq)]
pub struct MidpointCountSelection {
    group_indices: Vec<usize>,
    groups: Vec<MidpointGroupRow>,
    length_bins: Vec<LengthBin>,
    position_bins: Vec<MidpointPositionBin>,
    counts: DenseArray3<f32>,
}

impl MidpointCountSelection {
    /// Return selected source group indices in selection order.
    pub fn group_indices(&self) -> &[usize] {
        &self.group_indices
    }

    /// Return selected group metadata in selection order.
    pub fn groups(&self) -> &[MidpointGroupRow] {
        &self.groups
    }

    /// Return selected fragment length bins in selection order.
    pub fn length_bins(&self) -> &[LengthBin] {
        &self.length_bins
    }

    /// Return selected position bins in selection order.
    pub fn position_bins(&self) -> &[MidpointPositionBin] {
        &self.position_bins
    }

    /// Return the selected array shape as `(groups, length_bins, positions)`.
    pub fn shape(&self) -> (usize, usize, usize) {
        self.counts.shape()
    }

    /// Return the number of selected groups.
    pub fn group_count(&self) -> usize {
        self.groups.len()
    }

    /// Return the number of selected fragment length bins.
    pub fn length_bin_count(&self) -> usize {
        self.length_bins.len()
    }

    /// Return the number of selected position bins.
    pub fn position_bin_count(&self) -> usize {
        self.position_bins.len()
    }

    /// Return selected counts as a dense `(group, length_bin, position)` array.
    pub fn counts(&self) -> &DenseArray3<f32> {
        &self.counts
    }

    /// Return selected count values in row-major order.
    pub fn counts_row_major(&self) -> &[f32] {
        self.counts.values_row_major()
    }

    /// Return one selected count value, if all selection indices are in bounds.
    ///
    /// Parameters
    /// ----------
    /// - `group_index`:
    ///     Zero-based group index within the selected array.
    /// - `length_bin_index`:
    ///     Zero-based fragment length-bin index within the selected array.
    /// - `position_index`:
    ///     Zero-based position-bin index within the selected array.
    pub fn count(
        &self,
        group_index: usize,
        length_bin_index: usize,
        position_index: usize,
    ) -> Option<f32> {
        self.counts
            .get(group_index, length_bin_index, position_index)
            .copied()
    }

    /// Return one selected profile over positions.
    ///
    /// Parameters
    /// ----------
    /// - `group_index`:
    ///     Zero-based group index within the selected array.
    /// - `length_bin_index`:
    ///     Zero-based fragment length-bin index within the selected array.
    pub fn profile(&self, group_index: usize, length_bin_index: usize) -> Option<&[f32]> {
        self.counts
            .values_along_third_axis(group_index, length_bin_index)
    }
}

/// Parser for one `cfdna midpoints` Zarr store.
///
/// The parser validates root metadata and reads axis metadata arrays during
/// loading. It does not read the `/counts` array into memory, because
/// `MidpointsOutput` keeps the path and reads count blocks only through
/// selection methods.
struct MidpointsParser {
    path: PathBuf,
}

impl MidpointsParser {
    /// Store the input path until `load()` opens it.
    fn new(path: &Path) -> Self {
        Self {
            path: path.to_path_buf(),
        }
    }

    /// Validate metadata axes and build a midpoint output handle.
    fn load(&self) -> Result<MidpointsOutput> {
        validate_zarr_store_path(&self.path)?;
        let root_attributes = read_zarr_root_attributes(&self.path)
            .with_context(|| format!("read midpoint Zarr metadata {}", self.path.display()))?;
        validate_root_metadata(&root_attributes)?;
        let store = Arc::new(
            FilesystemStore::new(&self.path)
                .with_context(|| format!("open midpoint Zarr store {}", self.path.display()))?,
        );

        let counts_array = Array::open(store.clone(), "/counts")?;
        let count_shape = count_shape(&counts_array)?;
        let groups = read_groups(&self.path, store.clone(), count_shape.0)?;
        let group_name_indices = build_group_name_indices(&groups)?;
        let length_bins = read_indexed_axis(
            store.clone(),
            "/length_bin",
            "/length_start_bp",
            "/length_end_bp",
            count_shape.1,
            "length_bin",
            IndexedAxisKind::FragmentLength,
        )?;
        let position_bins = read_indexed_axis(
            store,
            "/position",
            "/position_bin_start_bp",
            "/position_bin_end_bp",
            count_shape.2,
            "position",
            IndexedAxisKind::Position,
        )?;

        Ok(MidpointsOutput {
            path: self.path.clone(),
            groups,
            group_name_indices,
            length_bins,
            position_bins,
            count_shape,
        })
    }
}

/// Validate that a path points to a readable midpoint Zarr store directory.
fn validate_zarr_store_path(path: &Path) -> Result<()> {
    ensure!(
        path.exists(),
        "midpoint Zarr store does not exist: {}",
        path.display()
    );
    ensure!(
        path.is_dir(),
        "midpoint Zarr store is not a directory: {}",
        path.display()
    );
    ensure!(
        path.file_name()
            .and_then(|name| name.to_str())
            .is_some_and(|name| name.ends_with(".zarr")),
        "midpoint Zarr store path must end with .zarr: {}",
        path.display()
    );
    ensure!(
        path.join("zarr.json").is_file(),
        "midpoint Zarr store is missing root zarr.json: {}",
        path.display()
    );
    Ok(())
}

/// Validate root attributes that identify the midpoint output schema.
fn validate_root_metadata(attributes: &Value) -> Result<()> {
    ensure!(
        string_attr(attributes, "cfdnalab_schema")? == "midpoint_profiles",
        "midpoint Zarr schema mismatch"
    );
    ensure!(
        u64_attr(attributes, "cfdnalab_schema_version")? == MIDPOINT_SCHEMA_VERSION,
        "midpoint Zarr schema version mismatch: expected {}",
        MIDPOINT_SCHEMA_VERSION
    );
    ensure!(
        string_attr(attributes, "primary_array")? == "counts",
        "midpoint Zarr primary_array must be counts"
    );
    ensure!(
        string_attr(attributes, "count_units")? == "weighted_midpoint_count",
        "midpoint Zarr count_units must be weighted_midpoint_count"
    );
    Ok(())
}

/// Read and validate the `(group, length_bin, position)` count-array shape.
fn count_shape(counts_array: &Array<FilesystemStore>) -> Result<(usize, usize, usize)> {
    let shape = counts_array
        .shape()
        .iter()
        .map(|dimension| usize::try_from(*dimension).context("Zarr dimension exceeds usize"))
        .collect::<Result<Vec<_>>>()?;
    ensure!(
        shape.len() == 3,
        "midpoint counts must be rank 3, found rank {}",
        shape.len()
    );
    ensure!(
        shape.iter().all(|dimension| *dimension > 0),
        "midpoint counts shape must be non-empty, found {:?}",
        shape
    );
    Ok((shape[0], shape[1], shape[2]))
}

/// Read group-axis metadata and validate it against the count-array shape.
fn read_groups(
    root_path: &Path,
    store: Arc<FilesystemStore>,
    expected_group_count: usize,
) -> Result<Vec<MidpointGroupRow>> {
    // The group coordinate axis defines the first count-array dimension
    let group = read_zarr_array1::<i32>(store.clone(), "/group")?;
    validate_zero_based_axis(&group, "group")?;
    ensure_same_len(&group, expected_group_count, "group")?;
    // Group names and eligible interval counts have one metadata value per
    // group coordinate
    let group_names = read_zarr_labels(root_path, "group", "group_name", expected_group_count)?;
    let eligible_intervals = read_zarr_array1::<i32>(store, "/eligible_intervals")?;
    ensure_same_len(
        &eligible_intervals,
        expected_group_count,
        "eligible_intervals",
    )?;

    (0..expected_group_count)
        .map(|group_index| {
            Ok(MidpointGroupRow {
                index: group_index,
                name: group_names[group_index].clone(),
                eligible_intervals: u64_from_i32(
                    eligible_intervals[group_index],
                    "eligible_intervals",
                )?,
            })
        })
        .collect()
}

/// Validation policy for one indexed interval axis.
#[derive(Debug, Clone, Copy)]
enum IndexedAxisKind {
    /// Fragment length bins share command-wide length bounds.
    FragmentLength,
    /// Position bins only need valid, contiguous intervals.
    Position,
}

/// Read one indexed interval axis from index, start, and end arrays.
fn read_indexed_axis(
    store: Arc<FilesystemStore>,
    index_array_path: &str,
    start_array_path: &str,
    end_array_path: &str,
    expected_len: usize,
    axis_name: &str,
    axis_kind: IndexedAxisKind,
) -> Result<Vec<IndexedInterval<u32, usize>>> {
    // Validate the coordinate axis first so starts and ends can be indexed by
    // zero-based axis position
    let indices = read_zarr_array1::<i32>(store.clone(), index_array_path)?;
    validate_zero_based_axis(&indices, axis_name)?;
    ensure_same_len(&indices, expected_len, axis_name)?;
    // Starts and ends carry the biological interval represented by each axis
    // coordinate
    let starts = read_zarr_array1::<i32>(store.clone(), start_array_path)?;
    let ends = read_zarr_array1::<i32>(store, end_array_path)?;
    ensure_same_len(&starts, expected_len, start_array_path)?;
    ensure_same_len(&ends, expected_len, end_array_path)?;

    let mut intervals = Vec::with_capacity(expected_len);
    let mut previous_end = None;
    for axis_index in 0..expected_len {
        let start = u32_from_i32(starts[axis_index], start_array_path)?;
        let end = u32_from_i32(ends[axis_index], end_array_path)?;
        validate_indexed_axis_interval(axis_name, axis_index, start, end, axis_kind)?;
        if let Some(previous_end_bp) = previous_end {
            ensure!(
                start == previous_end_bp,
                "{axis_name} intervals must be contiguous and sorted: index {axis_index} starts at {start}, previous end was {previous_end_bp}"
            );
        }
        let interval = IndexedInterval::new(start, end, axis_index).map_err(|error| {
            anyhow::anyhow!("{axis_name} index {axis_index} has invalid interval: {error}")
        })?;
        previous_end = Some(end);
        intervals.push(interval);
    }
    Ok(intervals)
}

/// Validate one midpoint axis interval before it becomes public metadata.
fn validate_indexed_axis_interval(
    axis_name: &str,
    axis_index: usize,
    start: u32,
    end: u32,
    axis_kind: IndexedAxisKind,
) -> Result<()> {
    ensure!(
        end > start,
        "{axis_name} index {axis_index} has invalid interval: end {end} <= start {start}"
    );
    if matches!(axis_kind, IndexedAxisKind::FragmentLength) {
        ensure!(
            start >= MIN_ACGT_BASES_FOR_GC_FRACTION,
            "{axis_name} index {axis_index} starts below minimum supported fragment length {}",
            MIN_ACGT_BASES_FOR_GC_FRACTION
        );
        ensure!(
            end <= MAX_SUPPORTED_FRAGMENT_LENGTH + 1,
            "{axis_name} index {axis_index} ends above maximum supported fragment length edge {}",
            MAX_SUPPORTED_FRAGMENT_LENGTH + 1
        );
    }
    Ok(())
}

/// Build a group-name lookup while rejecting duplicate group names.
fn build_group_name_indices(groups: &[MidpointGroupRow]) -> Result<FxHashMap<String, usize>> {
    let mut group_name_indices =
        FxHashMap::with_capacity_and_hasher(groups.len(), Default::default());
    for group in groups {
        ensure!(
            group_name_indices
                .insert(group.name.clone(), group.index)
                .is_none(),
            "midpoint group name is not unique: {}",
            group.name
        );
    }
    Ok(group_name_indices)
}

/// Copy selected axis metadata values in requested selector order.
fn selected_axis_values<T: Copy>(values: &[T], indices: &[usize], label: &str) -> Result<Vec<T>> {
    indices
        .iter()
        .map(|&index| {
            values
                .get(index)
                .copied()
                .with_context(|| format!("{label} index {index} is outside 0..{}", values.len()))
        })
        .collect()
}

/// Read one fully contiguous `(group, length_bin, position)` count block.
fn read_contiguous_count_block(
    counts_array: &Array<FilesystemStore>,
    group_span: (usize, usize),
    length_span: (usize, usize),
    position_span: (usize, usize),
) -> Result<DenseArray3<f32>> {
    let subset = count_subset(group_span, length_span, position_span)?;
    let values = counts_array
        .retrieve_array_subset::<Vec<f32>>(&subset)
        .context("read midpoint count block")?;
    dense_midpoint_counts_from_row_major(
        values,
        group_span.1 - group_span.0,
        length_span.1 - length_span.0,
        position_span.1 - position_span.0,
    )
}

/// Read contiguous group and length axes, then filter requested positions.
fn read_contiguous_group_length_count_block(
    counts_array: &Array<FilesystemStore>,
    group_span: (usize, usize),
    length_span: (usize, usize),
    position_indices: &[usize],
) -> Result<DenseArray3<f32>> {
    // Read the smallest position span that covers all requested positions, then
    // copy only requested positions into the final row-major output order
    let position_span = bounding_index_span(position_indices)?;
    let subset = count_subset(group_span, length_span, position_span)?;
    let block_values = counts_array
        .retrieve_array_subset::<Vec<f32>>(&subset)
        .context("read midpoint count block")?;
    let group_count = group_span.1 - group_span.0;
    let length_bin_count = length_span.1 - length_span.0;
    let position_span_len = position_span.1 - position_span.0;
    let mut values = Vec::with_capacity(
        group_count
            .saturating_mul(length_bin_count)
            .saturating_mul(position_indices.len()),
    );

    for group_offset in 0..group_count {
        for length_bin_offset in 0..length_bin_count {
            let profile_start =
                (group_offset * length_bin_count + length_bin_offset) * position_span_len;
            for &position_index in position_indices {
                values.push(block_values[profile_start + position_index - position_span.0]);
            }
        }
    }

    dense_midpoint_counts_from_row_major(
        values,
        group_count,
        length_bin_count,
        position_indices.len(),
    )
}

/// Read non-contiguous group or length selections in requested output order.
fn read_ordered_count_selection(
    counts_array: &Array<FilesystemStore>,
    group_indices: &[usize],
    length_bin_indices: &[usize],
    position_indices: &[usize],
) -> Result<DenseArray3<f32>> {
    // Each selected group and length-bin pair corresponds to one midpoint
    // profile across the position axis
    let mut values = Vec::with_capacity(
        group_indices
            .len()
            .saturating_mul(length_bin_indices.len())
            .saturating_mul(position_indices.len()),
    );
    let position_span = bounding_index_span(position_indices)?;
    for &group_index in group_indices {
        for &length_bin_index in length_bin_indices {
            let profile_subset = count_subset(
                (group_index, group_index + 1),
                (length_bin_index, length_bin_index + 1),
                position_span,
            )?;
            let profile_values = counts_array
                .retrieve_array_subset::<Vec<f32>>(&profile_subset)
                .with_context(|| {
                    format!(
                        "read midpoint counts for group {group_index}, length bin {length_bin_index}"
                    )
                })?;
            for &position_index in position_indices {
                let offset = position_index - position_span.0;
                values.push(profile_values[offset]);
            }
        }
    }
    dense_midpoint_counts_from_row_major(
        values,
        group_indices.len(),
        length_bin_indices.len(),
        position_indices.len(),
    )
}

/// Build a midpoint count array after rejecting non-finite selected values.
fn dense_midpoint_counts_from_row_major(
    values: Vec<f32>,
    group_count: usize,
    length_bin_count: usize,
    position_count: usize,
) -> Result<DenseArray3<f32>> {
    ensure_finite_midpoint_counts(&values)?;
    DenseArray3::from_row_major(values, group_count, length_bin_count, position_count)
}

/// Validate midpoint count values while allowing finite negative smoothed values.
fn ensure_finite_midpoint_counts(values: &[f32]) -> Result<()> {
    for (value_index, &value) in values.iter().enumerate() {
        ensure!(
            value.is_finite(),
            "midpoint counts contain non-finite value at selected row-major index {value_index}: {value}"
        );
    }
    Ok(())
}

/// Build a Zarr subset from half-open spans on all three count axes.
fn count_subset(
    group_span: (usize, usize),
    length_span: (usize, usize),
    position_span: (usize, usize),
) -> Result<ArraySubset> {
    Ok(ArraySubset::new_with_ranges(&[
        to_u64_range(group_span, "group")?,
        to_u64_range(length_span, "length bin")?,
        to_u64_range(position_span, "position")?,
    ]))
}

/// Convert a usize half-open span to the u64 range type expected by zarrs.
fn to_u64_range(span: (usize, usize), label: &str) -> Result<Range<u64>> {
    ensure!(
        span.0 <= span.1,
        "{label} span start {} is greater than end {}",
        span.0,
        span.1
    );
    Ok(
        u64::try_from(span.0).context("Zarr subset start exceeds u64")?
            ..u64::try_from(span.1).context("Zarr subset end exceeds u64")?,
    )
}

/// Return the smallest half-open span that contains all selected indices.
fn bounding_index_span(indices: &[usize]) -> Result<(usize, usize)> {
    let first_index = indices
        .first()
        .copied()
        .context("cannot build a midpoint count selection with an empty position axis")?;
    let mut min_index = first_index;
    let mut max_index = first_index;
    for &index in &indices[1..] {
        min_index = min_index.min(index);
        max_index = max_index.max(index);
    }
    Ok((min_index, max_index + 1))
}

/// Read one rank-1 Zarr array into a vector.
fn read_zarr_array1<T>(store: Arc<FilesystemStore>, array_path: &str) -> Result<Vec<T>>
where
    T: ElementOwned,
{
    let array = Array::open(store, array_path)?;
    let shape = array
        .shape()
        .iter()
        .map(|dimension| usize::try_from(*dimension).context("Zarr dimension exceeds usize"))
        .collect::<Result<Vec<_>>>()?;
    ensure!(
        shape.len() == 1,
        "Zarr array {array_path} must be rank 1, found rank {}",
        shape.len()
    );
    Ok(array.retrieve_array_subset(&array.subset_all())?)
}

/// Read string labels from a Zarr array's `labels` attribute.
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

/// Read the `attributes` object from one Zarr array metadata file.
fn read_zarr_array_attributes(root_path: &Path, array_path: &str) -> Result<Value> {
    let metadata_path = zarr_metadata_path(root_path, array_path);
    let metadata: Value = serde_json::from_str(
        &fs::read_to_string(&metadata_path)
            .with_context(|| format!("read Zarr metadata {}", metadata_path.display()))?,
    )
    .with_context(|| format!("parse Zarr metadata {}", metadata_path.display()))?;
    Ok(metadata.get("attributes").cloned().with_context(|| {
        format!(
            "Zarr metadata {} is missing attributes",
            metadata_path.display()
        )
    })?)
}

/// Return the metadata path for one Zarr array inside a store.
fn zarr_metadata_path(root_path: &Path, array_path: &str) -> PathBuf {
    let mut path = root_path.to_path_buf();
    for component in array_path.trim_start_matches('/').split('/') {
        path.push(component);
    }
    path.join("zarr.json")
}

/// Require an axis coordinate array to equal `0..len`.
fn validate_zero_based_axis(values: &[i32], axis_name: &str) -> Result<()> {
    for (index, &value) in values.iter().enumerate() {
        ensure!(
            value == i32::try_from(index).context("axis index exceeds i32")?,
            "{axis_name} must be a zero-based coordinate axis"
        );
    }
    Ok(())
}

/// Require an array to have the expected axis length.
fn ensure_same_len<T>(values: &[T], expected_len: usize, array_name: &str) -> Result<()> {
    ensure!(
        values.len() == expected_len,
        "{array_name} has {} entries, expected {expected_len}",
        values.len()
    );
    Ok(())
}

/// Read a required string attribute from Zarr metadata.
fn string_attr<'a>(attributes: &'a Value, name: &str) -> Result<&'a str> {
    attributes
        .get(name)
        .and_then(Value::as_str)
        .with_context(|| format!("midpoint Zarr metadata is missing string attribute {name}"))
}

/// Read a required unsigned integer attribute from Zarr metadata.
fn u64_attr(attributes: &Value, name: &str) -> Result<u64> {
    attributes
        .get(name)
        .and_then(Value::as_u64)
        .with_context(|| format!("midpoint Zarr metadata is missing integer attribute {name}"))
}

/// Convert a non-negative i32 Zarr value to u32.
fn u32_from_i32(value: i32, field_name: &str) -> Result<u32> {
    ensure!(value >= 0, "{field_name} must be non-negative, got {value}");
    Ok(u32::try_from(value).expect("non-negative i32 always fits u32"))
}

/// Convert a non-negative i32 Zarr value to u64.
fn u64_from_i32(value: i32, field_name: &str) -> Result<u64> {
    ensure!(value >= 0, "{field_name} must be non-negative, got {value}");
    Ok(u64::try_from(value).expect("non-negative i32 always fits u64"))
}
