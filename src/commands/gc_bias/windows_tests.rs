use super::*;
use crate::{
    commands::cli_common::WindowSpec,
    commands::gc_bias::counting::build_gc_prefixes,
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
    let err = fixed_size_window_interval(1144, 100_000, 114_364_328)
        .expect_err("out-of-range fixed window index should fail");

    assert!(
        format!("{err}").contains("beyond chromosome length"),
        "unexpected error message: {err}"
    );
    Ok(())
}

#[test]
fn prepare_tile_windows_bed_keeps_a_right_halo_only_window_for_fragment_owned_models() -> Result<()>
{
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

    let prepared = prepare_tile_windows(
        &window_opt,
        Some(&windows),
        &tile,
        Some(&span),
        200,
        &template,
    )?;

    assert!(!prepared.skip_tile);
    assert!(prepared.streaming_buffers.is_none());
    assert_eq!(prepared.windows.len(), 1);
    assert_eq!(prepared.windows[0].idx, 0);
    assert_eq!(prepared.windows[0].interval, Interval::new(22, 23)?);
    assert!(!prepared.windows[0].contained);
    Ok(())
}

#[test]
fn set_window_acgt_uses_sequence_interval_as_prefix_origin() -> Result<()> {
    // The prefix arrays are built from reference slice [900,961), not from chromosome start.
    // The observed interval [930,941) should therefore be counted as local [30,41).
    //
    // The repeated ACGT reference starts at A at coordinate 900. In [930,941), all 11 bases are
    // A/C/G/T, so the ACGT support numerator and denominator should both be 11.
    let seq = b"ACGT".repeat(16);
    let prefixes = build_gc_prefixes(&seq);
    let sequence_interval = Interval::new(900_u64, 961_u64)?;
    let observed_interval = Interval::new(930_u64, 941_u64)?;
    let mut window = WindowState::new(0, observed_interval, true, &make_template())?;

    // Act
    set_window_acgt_in_observed_interval(
        &mut window,
        &prefixes,
        observed_interval,
        sequence_interval,
    )?;

    // Assert
    assert_eq!(window.counts.num_acgt_out_of, (11, 11));
    Ok(())
}

#[test]
fn set_window_acgt_clips_observed_interval_to_loaded_sequence() -> Result<()> {
    // Manual derivation:
    // - Prefixes cover reference slice [900,910), all A/C/G/T.
    // - The observed interval [895,905) overlaps only [900,905) of that loaded slice.
    // - The support denominator is therefore the clipped length 5, not the original 10.
    let prefixes = build_gc_prefixes(&vec![b'C'; 10]);
    let sequence_interval = Interval::new(900_u64, 910_u64)?;
    let observed_interval = Interval::new(895_u64, 905_u64)?;
    let mut window = WindowState::new(0, observed_interval, true, &make_template())?;

    set_window_acgt_in_observed_interval(
        &mut window,
        &prefixes,
        observed_interval,
        sequence_interval,
    )?;

    assert_eq!(window.counts.num_acgt_out_of, (5, 5));
    Ok(())
}

#[test]
fn set_window_acgt_errors_when_observed_interval_misses_loaded_sequence() -> Result<()> {
    // Manual derivation:
    // - Prefixes cover reference slice [900,910).
    // - Observed interval [880,890) has no overlap with that loaded slice.
    // - This is a caller error for support bookkeeping and should produce a clear error instead
    //   of silently recording zero support.
    let prefixes = build_gc_prefixes(&vec![b'C'; 10]);
    let sequence_interval = Interval::new(900_u64, 910_u64)?;
    let observed_interval = Interval::new(880_u64, 890_u64)?;
    let mut window = WindowState::new(0, observed_interval, true, &make_template())?;

    let error = set_window_acgt_in_observed_interval(
        &mut window,
        &prefixes,
        observed_interval,
        sequence_interval,
    )
    .expect_err("missing loaded-sequence overlap should error");

    assert!(
        error.to_string().contains("does not overlap sequence"),
        "unexpected error: {error}"
    );
    Ok(())
}

#[test]
fn prepare_tile_windows_bed_keeps_core_and_right_halo_windows_together_for_fragment_owned_models()
-> Result<()> {
    // Manual derivation:
    // - BED window [10,11) overlaps the core and is fully contained.
    // - BED window [22,23) lies outside the core but is still reachable by a tile-owned fragment.
    // - `gc_bias` BED preparation should therefore build two window states, one contained and one
    //   boundary/crossing.
    let template = make_template();
    let tile = Tile::from_coords("chr1".to_string(), 0, 0, 10, 20, 6, 24)
        .expect("test tile should be valid");
    let windows = vec![
        IndexedInterval::new(10, 11, 0)?,
        IndexedInterval::new(22, 23, 1)?,
    ];
    let span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 2,
    };
    let window_opt = WindowSpec::Bed(PathBuf::from("windows.bed"));

    let prepared = prepare_tile_windows(
        &window_opt,
        Some(&windows),
        &tile,
        Some(&span),
        200,
        &template,
    )?;

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

#[test]
fn prepare_tile_windows_bed_errors_when_windows_are_missing() -> Result<()> {
    // BED mode requires a chromosome window slice. The helper should fail loudly instead of
    // silently treating missing windows as an empty BED file.
    let template = make_template();
    let tile = Tile::from_coords("chr1".to_string(), 0, 0, 10, 20, 6, 24)
        .expect("test tile should be valid");
    let span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 1,
    };
    let window_opt = WindowSpec::Bed(PathBuf::from("windows.bed"));

    let err = prepare_tile_windows(&window_opt, None, &tile, Some(&span), 200, &template)
        .expect_err("BED mode without windows should fail");

    assert!(
        format!("{err}").contains("no windows provided"),
        "unexpected error message: {err}"
    );
    Ok(())
}

#[test]
fn prepare_tile_windows_bed_skips_tile_for_empty_cached_span() -> Result<()> {
    // An explicitly empty cached candidate span means BED precomputation already proved that the
    // tile has no relevant windows, so the helper should return skip=true without building states.
    let template = make_template();
    let tile = Tile::from_coords("chr1".to_string(), 0, 0, 10, 20, 6, 24)
        .expect("test tile should be valid");
    let windows = vec![IndexedInterval::new(10, 11, 0)?];
    let span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 0,
    };
    let window_opt = WindowSpec::Bed(PathBuf::from("windows.bed"));

    let prepared = prepare_tile_windows(
        &window_opt,
        Some(&windows),
        &tile,
        Some(&span),
        200,
        &template,
    )?;

    assert!(prepared.skip_tile);
    assert!(prepared.windows.is_empty());
    assert!(prepared.streaming_buffers.is_none());
    Ok(())
}

#[test]
fn prepare_tile_windows_bed_skips_tile_for_empty_window_slice_even_with_nonempty_span() -> Result<()>
{
    // This is a helper-level defensive case: the caller supplies a non-empty cached span but the
    // chromosome window slice itself is empty. The helper clamps nothing and reports skip=true.
    let template = make_template();
    let tile = Tile::from_coords("chr1".to_string(), 0, 0, 10, 20, 6, 24)
        .expect("test tile should be valid");
    let windows: Vec<IndexedInterval<u64>> = Vec::new();
    let span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 2,
    };
    let window_opt = WindowSpec::Bed(PathBuf::from("windows.bed"));

    let prepared = prepare_tile_windows(
        &window_opt,
        Some(&windows),
        &tile,
        Some(&span),
        200,
        &template,
    )?;

    assert!(prepared.skip_tile);
    assert!(prepared.windows.is_empty());
    assert!(prepared.streaming_buffers.is_none());
    Ok(())
}

#[test]
fn prepare_tile_windows_bed_clamps_candidate_span_to_available_window_slice() -> Result<()> {
    // Manual derivation:
    // - The cached span says [1,4), but the chromosome slice has only two windows.
    // - The helper clamps that to the available suffix [1,2), so only the second window is built.
    let template = make_template();
    let tile = Tile::from_coords("chr1".to_string(), 0, 0, 10, 20, 6, 24)
        .expect("test tile should be valid");
    let windows = vec![
        IndexedInterval::new(10, 11, 0)?,
        IndexedInterval::new(22, 23, 1)?,
    ];
    let span = TileWindowSpan {
        first_idx: 1,
        last_idx_exclusive: 4,
    };
    let window_opt = WindowSpec::Bed(PathBuf::from("windows.bed"));

    let prepared = prepare_tile_windows(
        &window_opt,
        Some(&windows),
        &tile,
        Some(&span),
        200,
        &template,
    )?;

    assert!(!prepared.skip_tile);
    assert!(prepared.streaming_buffers.is_none());
    assert_eq!(prepared.windows.len(), 1);
    assert_eq!(prepared.windows[0].idx, 1);
    assert_eq!(prepared.windows[0].interval, Interval::new(22, 23)?);
    assert!(!prepared.windows[0].contained);
    Ok(())
}

#[test]
fn prepare_tile_windows_bed_returns_empty_windows_when_clamped_span_starts_past_slice_end()
-> Result<()> {
    // Manual derivation:
    // - The cached span says [5,8), but the chromosome slice has length 2.
    // - Clamping both bounds to len=2 yields [2,2), so no windows are built.
    // - This is a defensive helper behavior, distinct from skip=true on a truly empty span.
    let template = make_template();
    let tile = Tile::from_coords("chr1".to_string(), 0, 0, 10, 20, 6, 24)
        .expect("test tile should be valid");
    let windows = vec![
        IndexedInterval::new(10, 11, 0)?,
        IndexedInterval::new(22, 23, 1)?,
    ];
    let span = TileWindowSpan {
        first_idx: 5,
        last_idx_exclusive: 8,
    };
    let window_opt = WindowSpec::Bed(PathBuf::from("windows.bed"));

    let prepared = prepare_tile_windows(
        &window_opt,
        Some(&windows),
        &tile,
        Some(&span),
        200,
        &template,
    )?;

    assert!(!prepared.skip_tile);
    assert!(prepared.streaming_buffers.is_none());
    assert!(prepared.windows.is_empty());
    Ok(())
}

#[test]
fn prepare_tile_windows_bed_marks_exact_core_boundary_windows_as_contained() -> Result<()> {
    // Manual derivation:
    // - Containment uses `start >= core_start && end <= core_end`.
    // - Therefore windows that start exactly at the core start or end exactly at the core end are
    //   still contained.
    let template = make_template();
    let tile = Tile::from_coords("chr1".to_string(), 0, 0, 10, 20, 6, 24)
        .expect("test tile should be valid");
    let windows = vec![
        IndexedInterval::new(10, 12, 0)?,
        IndexedInterval::new(18, 20, 1)?,
    ];
    let span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 2,
    };
    let window_opt = WindowSpec::Bed(PathBuf::from("windows.bed"));

    let prepared = prepare_tile_windows(
        &window_opt,
        Some(&windows),
        &tile,
        Some(&span),
        200,
        &template,
    )?;

    assert!(!prepared.skip_tile);
    assert_eq!(prepared.windows.len(), 2);
    assert!(prepared.windows[0].contained);
    assert!(prepared.windows[1].contained);
    Ok(())
}
