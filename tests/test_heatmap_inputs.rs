#![cfg(feature = "plotters")]

use anyhow::Result;
use cfdnalab::shared::plotters::heatmap::{HeatmapFormat, HeatmapUpsample, write_heatmap};
use ndarray::Array2;
use tempfile::TempDir;

#[test]
fn renders_constant_heatmap_when_values_are_flat() -> Result<()> {
    // Human verification status: unverified
    // Arrange: create a flat matrix and matching edges
    let values = Array2::from_elem((2, 3), 0.0);
    let x_edges: Vec<f64> = (0..=3).map(|idx| idx as f64).collect();
    let y_edges: Vec<f64> = (0..=2).map(|idx| idx as f64).collect();
    let temp_dir = TempDir::new()?;
    let out_path = temp_dir.path().join("flat.png");

    // Act: attempt to render with collapsed limits (min == max)
    write_heatmap(
        &out_path,
        "Flat heatmap",
        "x",
        "y",
        &values,
        Some(&x_edges),
        Some(&y_edges),
        Some(0.0),
        None,
        None,
        None,
        None,
        false,
        1,
        HeatmapUpsample::Nearest,
        160,
        120,
        HeatmapFormat::Png,
    )?;

    // Assert: file is written and non-empty
    let metadata = std::fs::metadata(&out_path)?;
    assert!(metadata.is_file());
    assert!(metadata.len() > 0);
    Ok(())
}
