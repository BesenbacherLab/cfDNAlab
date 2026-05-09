use super::{WpsTileTempOutput, collect_wps_aggregate_tile_outputs_by_chromosome};
use anyhow::Result;
use std::path::PathBuf;

#[test]
fn wps_reducer_inputs_are_derived_from_returned_typed_tile_outputs() -> Result<()> {
    // Arrange
    // WPS reuses the fcoverage aggregate reducers, but WPS orchestration must still pass the typed
    // tile outputs returned by tile processing. Positional outputs and unrelated aggregate files
    // must not be rediscovered from the temp directory by prefix.
    let tile_outputs = vec![
        WpsTileTempOutput::Positional {
            chromosome: "chr1".to_string(),
            tile_index: 0,
            path: PathBuf::from("tile0.positional"),
        },
        WpsTileTempOutput::AggregatesBySize {
            chromosome: "chr2".to_string(),
            tile_index: 1,
            partials_path: PathBuf::from("chr2.tile1.size_partials"),
            cross_index_path: Some(PathBuf::from("chr2.tile1.size_cross")),
        },
        WpsTileTempOutput::AggregatesBySize {
            chromosome: "chr1".to_string(),
            tile_index: 0,
            partials_path: PathBuf::from("chr1.tile0.size_partials"),
            cross_index_path: None,
        },
        WpsTileTempOutput::AggregatesByBed {
            chromosome: "chr1".to_string(),
            tile_index: 9,
            partials_path: PathBuf::from("wrong_mode.bed_partials"),
            cross_index_path: Some(PathBuf::from("wrong_mode.bed_cross")),
        },
    ];

    // Act
    let by_chromosome = collect_wps_aggregate_tile_outputs_by_chromosome(
        &tile_outputs,
        WpsTileTempOutput::is_size_aggregate,
    )?;

    // Assert
    assert_eq!(by_chromosome.len(), 2);
    assert_eq!(by_chromosome["chr1"].len(), 1);
    assert_eq!(by_chromosome["chr1"][0].tile_index, 0);
    assert_eq!(
        by_chromosome["chr1"][0].partials_path,
        PathBuf::from("chr1.tile0.size_partials")
    );
    assert_eq!(by_chromosome["chr1"][0].cross_index_path, None);
    assert_eq!(by_chromosome["chr2"].len(), 1);
    assert_eq!(by_chromosome["chr2"][0].tile_index, 1);
    assert_eq!(
        by_chromosome["chr2"][0].partials_path,
        PathBuf::from("chr2.tile1.size_partials")
    );
    assert_eq!(
        by_chromosome["chr2"][0].cross_index_path,
        Some(PathBuf::from("chr2.tile1.size_cross"))
    );

    Ok(())
}
