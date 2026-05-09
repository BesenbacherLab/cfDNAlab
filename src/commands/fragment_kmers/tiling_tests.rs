use crate::{
    commands::cli_common::WindowSpec,
    shared::{
        interval::IndexedInterval,
        tiled_run::{Tile, TileWindowSpan},
        window_fetch::{BedFetchPolicy, fetch_span_for_tile},
    },
};
use anyhow::Result;
use std::path::PathBuf;

#[test]
fn fetch_span_for_tile_core_overlap_keeps_fragment_halo_for_bed_windows_at_chromosome_start()
-> Result<()> {
    // Arrange:
    // Use a chromosome-start tile whose core is [0, 100) and whose fetch band was built with a
    // 60 bp fragment halo, truncated on the left chromosome edge:
    //   core  = [0, 100)
    //   fetch = [0, 160)
    //
    // Keep one BED window [40, 50) inside that tile. `fragment-kmers` counts selected k-mer
    // positions from fragments, so a fragment may start left of the BED window but still
    // contribute positions inside it. With a max-fragment-length halo of 60, the correct
    // narrowed fetch span is therefore:
    //   [40 - 60, 50 + 60) clamped to tile fetch band = [0, 110)
    let tile = Tile::from_coords("chr1".to_string(), 0, 0, 0, 100, 0, 160)?;
    let windows = vec![IndexedInterval::new(40, 50, 0)?];
    let tile_window_span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 1,
    };
    let window_spec = WindowSpec::Bed(PathBuf::from("windows.bed"));
    // Act
    let fetch_span = fetch_span_for_tile(
        &tile,
        Some(&tile_window_span),
        Some(&windows),
        &window_spec,
        200,
        60,
        BedFetchPolicy::CoreOverlap,
    )?
    .expect("the tile overlaps the BED window");

    // Assert
    assert_eq!(fetch_span.as_tuple(), (0, 110));

    Ok(())
}

#[test]
fn fetch_span_for_tile_core_overlap_returns_none_for_halo_only_bed_windows() -> Result<()> {
    // Arrange:
    // - `fragment-kmers` uses core-overlap BED semantics
    // - the cached span may still cover halo-only candidates from another model, but this helper
    //   must ignore them and skip the tile when no BED window actually overlaps the core
    let tile = Tile::from_coords("chr1".to_string(), 0, 0, 10, 20, 6, 24)?;
    let windows = vec![IndexedInterval::new(22, 23, 0)?];
    let tile_window_span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 1,
    };
    let window_spec = WindowSpec::Bed(PathBuf::from("windows.bed"));
    let fetch_span = fetch_span_for_tile(
        &tile,
        Some(&tile_window_span),
        Some(&windows),
        &window_spec,
        200,
        4,
        BedFetchPolicy::CoreOverlap,
    )?;

    assert!(fetch_span.is_none());
    Ok(())
}

#[test]
fn fetch_span_for_tile_core_overlap_ignores_halo_only_candidates_when_a_core_window_is_present()
-> Result<()> {
    // Arrange:
    // - BED window [10,11) overlaps the core and should define the narrowed fetch.
    // - BED window [22,23) is halo-only for this core-overlap command and must not widen the
    //   fetch span even if the cached index span includes it.
    // - With fragment halo 4, [10,11) widens to [6,15) inside tile fetch [6,24).
    let tile = Tile::from_coords("chr1".to_string(), 0, 0, 10, 20, 6, 24)?;
    let windows = vec![
        IndexedInterval::new(10, 11, 0)?,
        IndexedInterval::new(22, 23, 1)?,
    ];
    let tile_window_span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 2,
    };
    let window_spec = WindowSpec::Bed(PathBuf::from("windows.bed"));
    let fetch_span = fetch_span_for_tile(
        &tile,
        Some(&tile_window_span),
        Some(&windows),
        &window_spec,
        200,
        4,
        BedFetchPolicy::CoreOverlap,
    )?
    .expect("core-overlap window should keep a fetch span");

    assert_eq!(fetch_span.as_tuple(), (6, 15));
    Ok(())
}
