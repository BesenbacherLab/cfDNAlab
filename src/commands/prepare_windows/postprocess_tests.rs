use super::build_coverage_index;
use crate::commands::prepare_windows::{
    config::CoordinateSet, labels::LabelTuple, prepare_windows::Window,
};
use std::sync::Arc;

fn win(chrom: &str, start: u32, end: u32) -> Window {
    Window::from_bounds(
        Arc::<str>::from(chrom.to_string()),
        start,
        end,
        start,
        end,
        vec![LabelTuple::new("group".to_string())],
        "group".to_string(),
        None,
    )
    .expect("test window should be valid")
}

#[test]
fn build_coverage_index_combines_same_position_deltas_without_losing_coverage() {
    // Human verification status: unverified
    // Arrange
    // The boundaries at 15 and 20 contain both starts and ends. In particular, the net delta at
    // 15 is zero, so the segment split must be preserved without changing the running depth.
    let windows = vec![
        win("chr1", 10, 20),
        win("chr1", 10, 15),
        win("chr1", 15, 20),
        win("chr1", 20, 30),
    ];

    // Act
    let (boundaries, coverage_by_segment, coverage_prefix) =
        build_coverage_index(&windows, CoordinateSet::Resized);

    // Assert
    assert_eq!(boundaries, vec![10, 15, 20, 30]);
    assert_eq!(coverage_by_segment, vec![2, 2, 1]);
    assert_eq!(coverage_prefix, vec![0, 10, 20, 30]);
}
