/// Axis bounds used for rendering one positional track.
///
/// The axis is inclusive on both ends because the selections operate over
/// discrete base indices rather than floating-point coordinates.
#[derive(Debug, Clone)]
pub(crate) struct AxisBounds {
    pub(crate) start: i32,
    pub(crate) end: i32,
}

impl AxisBounds {
    /// Build a new inclusive axis range.
    ///
    /// Parameters
    /// ----------
    /// - `start`:
    ///   Inclusive axis start
    /// - `end`:
    ///   Inclusive axis end
    ///
    /// Returns
    /// -------
    /// - `AxisBounds`:
    ///   Axis object that can be reused across renderers
    pub(crate) fn new(start: i32, end: i32) -> Self {
        Self { start, end }
    }
}

/// One logical coordinate track for positional visualization.
///
/// Each track has a name, an axis definition, and the selected indices that
/// should be highlighted on that axis.
#[derive(Debug, Clone)]
pub(crate) struct Track {
    pub(crate) name: String,
    #[cfg_attr(not(feature = "cmd_visualize_positions"), allow(dead_code))]
    pub(crate) axis: AxisBounds,
    pub(crate) selected_indices: Vec<i32>,
}

impl Track {
    /// Check whether the track highlights any positions.
    ///
    /// Parameters
    /// ----------
    /// - `self`:
    ///   Track to inspect
    ///
    /// Returns
    /// -------
    /// - `bool`:
    ///   `true` when no indices are selected on the track
    #[cfg(any(feature = "cmd_visualize_positions"))]
    pub(crate) fn is_empty(&self) -> bool {
        self.selected_indices.is_empty()
    }
}
