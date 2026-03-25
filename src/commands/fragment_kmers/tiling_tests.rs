use super::determine_fetch_span;
use crate::{
    commands::{cli_common::WindowSpec, fragment_kmers::windows::WindowContext},
    shared::{
        interval::IndexedInterval,
        tiled_run::{Tile, TileWindowSpan},
    },
};
use anyhow::Result;
use std::path::PathBuf;

#[test]
fn determine_fetch_span_keeps_fragment_halo_for_bed_windows_at_chromosome_start() -> Result<()> {
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
    let window_ctx = WindowContext {
        spec: &window_spec,
        windows: Some(&windows),
        chr_idx_offset: 0,
    };

    // Act
    let fetch_span = determine_fetch_span(&tile, &window_ctx, Some(&tile_window_span), 200, 60)?
        .expect("the tile overlaps the BED window");

    // Assert
    assert_eq!(fetch_span.as_tuple(), (0, 110));

    Ok(())
}
