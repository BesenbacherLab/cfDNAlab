use super::*;
use crate::{
    commands::cli_common::WindowSpec,
    shared::{
        interval::{IndexedInterval, Interval},
        tiled_run::{Tile, TileWindowSpan},
    },
};
use std::path::PathBuf;

fn make_template() -> GCCounts {
    GCCounts::new(1, 1, 0, (0, 0)).expect("failed to build template")
}

#[test]
fn prepares_fixed_size_streaming_buffers_for_last_partial_window() -> Result<()> {
    // Human verification status: unverified
    // Chromosome length leaves a final partial 100 kb window
    // [114300000,114364328). The tile core lies in that last window, so there is no
    // valid window after it.
    let template = make_template();
    let tile = Tile::from_coords(
        "chr1".to_string(),
        0,
        11,
        114_300_000,
        114_364_328,
        114_299_000,
        114_364_328,
    )
    .expect("test tile should be valid");

    let prepared = prepare_fixed_size_streaming_buffers(
        100_000,
        114_364_328,
        tile.core.try_to_u64()?,
        &template,
    );

    assert!(
        prepared.is_ok(),
        "last partial window should not require an out-of-bounds next buffer: {prepared:?}"
    );
    let (current, next) = prepared?;
    assert_eq!(current.start(), 114_300_000);
    assert_eq!(current.end(), 114_364_328);
    assert!(
        next.is_none(),
        "the last partial window should not invent a following window"
    );
    Ok(())
}

#[test]
fn advances_fixed_size_streaming_buffers_into_last_partial_window() -> Result<()> {
    // Human verification status: unverified
    // Current window is the second-to-last 100 kb bin, next is the last partial bin,
    // and advancing once more must not try to build [114400000,114364328).
    let template = make_template();
    let chrom_len = 114_364_328_u64;
    let window_bp = 100_000_u64;
    let core_interval = Interval::new(110_000_000_u64, chrom_len)?;

    let current = window_state_from_idx(1142, window_bp, chrom_len, core_interval, &template)?;
    let next = window_state_from_idx(1143, window_bp, chrom_len, core_interval, &template)?;

    let advanced = advance_fixed_size_streaming_buffers(
        current,
        next,
        window_bp,
        chrom_len,
        core_interval,
        &template,
    );

    assert!(
        advanced.is_ok(),
        "advancing into the last partial window should not construct an invalid interval: {advanced:?}"
    );
    let (current, next) = advanced?;
    assert_eq!(current.start(), 114_300_000);
    assert_eq!(current.end(), 114_364_328);
    assert!(
        next.is_none(),
        "advancing past the last real window should leave no next window"
    );
    Ok(())
}

#[test]
fn rejects_fixed_size_window_index_past_chromosome_end() -> Result<()> {
    // Human verification status: unverified
    let err = fixed_size_window_interval(1144, 100_000, 114_364_328)
        .expect_err("out-of-range fixed window index should fail");

    assert!(
        format!("{err}").contains("beyond chromosome length"),
        "unexpected error message: {err}"
    );
    Ok(())
}

#[test]
fn prepare_tile_windows_bed_keeps_a_right_halo_only_window_for_fragment_owned_models() -> Result<()> {
    // Human verification status: unverified
    // Manual derivation:
    // - `gc_bias` owns fragments by aligned start in tile core [10,20).
    // - A fragment starting at 19 with aligned length 4 can overlap BED window [22,23).
    // - BED preparation must therefore keep that right-halo-only window, matching the fragment
    //   ownership model already used by the fixed-size path.
    let template = make_template();
    let tile = Tile::from_coords("chr1".to_string(), 0, 0, 10, 20, 6, 24)
        .expect("test tile should be valid");
    let windows = vec![IndexedInterval::new(22, 23, 0)?];
    let span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 1,
    };
    let window_opt = WindowSpec::Bed(PathBuf::from("windows.bed"));

    let prepared =
        prepare_tile_windows(&window_opt, Some(&windows), &tile, Some(&span), 200, &template)?;

    assert!(!prepared.skip_tile);
    assert!(prepared.streaming_buffers.is_none());
    assert_eq!(prepared.windows.len(), 1);
    assert_eq!(prepared.windows[0].idx, 0);
    assert_eq!(prepared.windows[0].interval, Interval::new(22, 23)?);
    assert!(!prepared.windows[0].contained);
    Ok(())
}

#[test]
fn prepare_tile_windows_bed_keeps_core_and_right_halo_windows_together_for_fragment_owned_models()
-> Result<()> {
    // Human verification status: unverified
    // Manual derivation:
    // - BED window [10,11) overlaps the core and is fully contained.
    // - BED window [22,23) lies outside the core but is still reachable by a tile-owned fragment.
    // - `gc_bias` BED preparation should therefore build two window states, one contained and one
    //   boundary/crossing.
    let template = make_template();
    let tile = Tile::from_coords("chr1".to_string(), 0, 0, 10, 20, 6, 24)
        .expect("test tile should be valid");
    let windows = vec![IndexedInterval::new(10, 11, 0)?, IndexedInterval::new(22, 23, 1)?];
    let span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 2,
    };
    let window_opt = WindowSpec::Bed(PathBuf::from("windows.bed"));

    let prepared =
        prepare_tile_windows(&window_opt, Some(&windows), &tile, Some(&span), 200, &template)?;

    assert!(!prepared.skip_tile);
    assert!(prepared.streaming_buffers.is_none());
    assert_eq!(prepared.windows.len(), 2);
    assert_eq!(prepared.windows[0].idx, 0);
    assert_eq!(prepared.windows[0].interval, Interval::new(10, 11)?);
    assert!(prepared.windows[0].contained);
    assert_eq!(prepared.windows[1].idx, 1);
    assert_eq!(prepared.windows[1].interval, Interval::new(22, 23)?);
    assert!(!prepared.windows[1].contained);
    Ok(())
}
