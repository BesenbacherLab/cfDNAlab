/// Axis bounds used for rendering one positional track.
///
/// The axis is inclusive on both ends because the selections operate over
/// discrete base indices rather than floating-point coordinates.
#[derive(Debug, Clone)]
pub struct AxisBounds {
    pub start: i32,
    pub end: i32,
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
    pub fn new(start: i32, end: i32) -> Self {
        Self { start, end }
    }

    /// Measure the span of the axis.
    ///
    /// Parameters
    /// ----------
    /// - `self`:
    ///   Axis to measure
    ///
    /// Returns
    /// -------
    /// - `i32`:
    ///   End minus start in axis units
    pub fn length(&self) -> i32 {
        self.end - self.start
    }
}

/// One logical coordinate track for positional visualization.
///
/// Each track has a name, an axis definition, and the selected indices that
/// should be highlighted on that axis.
#[derive(Debug, Clone)]
pub struct Track {
    pub name: String,
    pub axis: AxisBounds,
    pub selected_indices: Vec<i32>,
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
    pub fn is_empty(&self) -> bool {
        self.selected_indices.is_empty()
    }
}
