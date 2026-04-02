use super::*;
use crate::shared::{
    interval::{IndexedInterval, Interval},
    tiled_run::{Tile, TileWindowSpan, clamp_fetch_to_window_span},
};

fn indexed_windows(entries: &[(u64, u64, u64)]) -> Vec<IndexedInterval<u64>> {
    entries
        .iter()
        .map(|&(start, end, original_index)| {
            IndexedInterval::new(start, end, original_index)
                .expect("test windows should be valid non-empty intervals")
        })
        .collect()
}

fn make_tile(core_start: u32, core_end: u32, fetch_start: u32, fetch_end: u32, index: u32) -> Tile {
    Tile::from_coords(
        "chr1".to_string(),
        0,
        index,
        core_start,
        core_end,
        fetch_start,
        fetch_end,
    )
    .expect("test tile should be valid")
}

#[test]
fn aligned_window_extent_for_core_overlap_bed_keeps_one_core_window() -> Result<()> {
    let tile = make_tile(10, 20, 6, 24, 0);
    let windows = indexed_windows(&[(10, 11, 0)]);
    let candidate_span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 1,
    };

    let window_extent =
        aligned_window_extent_for_core_overlap_bed(&windows, &tile, Some(&candidate_span))?
            .expect("one core-overlap window should produce an aligned window extent");

    assert_eq!(window_extent, Interval::new(10, 11)?);
    Ok(())
}

#[test]
fn aligned_window_extent_for_core_overlap_bed_ignores_halo_only_candidates_in_mixed_cached_span()
-> Result<()> {
    let tile = make_tile(10, 20, 6, 24, 0);
    let windows = indexed_windows(&[(8, 9, 0), (10, 11, 1), (22, 23, 2)]);
    let candidate_span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 3,
    };

    let window_extent =
        aligned_window_extent_for_core_overlap_bed(&windows, &tile, Some(&candidate_span))?
            .expect("mixed cached span should still keep the core-overlap window extent");

    assert_eq!(window_extent, Interval::new(10, 11)?);
    Ok(())
}

#[test]
fn aligned_window_extent_for_core_overlap_bed_returns_none_when_cached_span_has_no_true_core_window()
-> Result<()> {
    let tile = make_tile(10, 20, 6, 24, 0);
    let windows = indexed_windows(&[(8, 9, 0), (22, 23, 1)]);
    let candidate_span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 2,
    };

    let window_extent =
        aligned_window_extent_for_core_overlap_bed(&windows, &tile, Some(&candidate_span))?;

    assert!(window_extent.is_none());
    Ok(())
}

#[test]
fn aligned_window_extent_for_core_overlap_bed_widens_monotonically_when_more_core_windows_are_added()
-> Result<()> {
    let tile = make_tile(10, 20, 6, 24, 0);
    let windows = indexed_windows(&[(10, 11, 0), (18, 19, 1)]);

    let narrow = aligned_window_extent_for_core_overlap_bed(
        &windows,
        &tile,
        Some(&TileWindowSpan {
            first_idx: 0,
            last_idx_exclusive: 1,
        }),
    )?
    .expect("single core window should produce a window extent");
    let wide = aligned_window_extent_for_core_overlap_bed(
        &windows,
        &tile,
        Some(&TileWindowSpan {
            first_idx: 0,
            last_idx_exclusive: 2,
        }),
    )?
    .expect("two core windows should produce a window extent");

    assert_eq!(narrow, Interval::new(10, 11)?);
    assert_eq!(wide, Interval::new(10, 19)?);
    assert!(wide.start() <= narrow.start());
    assert!(wide.end() >= narrow.end());
    Ok(())
}

#[test]
fn aligned_window_extent_for_bed_candidates_keeps_one_right_halo_only_window() -> Result<()> {
    let windows = indexed_windows(&[(22, 23, 0)]);
    let candidate_span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 1,
    };

    let window_extent = aligned_window_extent_for_bed_candidates(&windows, Some(&candidate_span))?
        .expect("one right halo-only candidate should produce an aligned window extent");

    assert_eq!(window_extent, Interval::new(22, 23)?);
    Ok(())
}

#[test]
fn aligned_window_extent_for_bed_candidates_uses_the_full_candidate_extent_in_mixed_case()
-> Result<()> {
    let windows = indexed_windows(&[(10, 11, 0), (22, 23, 1)]);
    let candidate_span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 2,
    };

    let window_extent = aligned_window_extent_for_bed_candidates(&windows, Some(&candidate_span))?
        .expect("mixed fragment-reach candidates should produce an aligned window extent");

    assert_eq!(window_extent, Interval::new(10, 23)?);
    Ok(())
}

#[test]
fn candidate_window_extent_and_clamp_keep_reads_starting_reach_bases_earlier() -> Result<()> {
    // Manual derivation:
    // - Tile core is [10,20) and the tile fetch band is [6,24).
    // - Candidate BED window [22,23) is kept because aligned fragment reach is 4 bp.
    // - The window-extent helper should return [22,23).
    // - Clamping with that same aligned reach must then widen left to 18 and right to 24.
    // - If the later caller used a smaller reach, the owned read starting at 19 would be lost.
    let tile = make_tile(10, 20, 6, 24, 0);
    let windows = indexed_windows(&[(22, 23, 0)]);
    let candidate_span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 1,
    };

    let window_extent = aligned_window_extent_for_bed_candidates(&windows, Some(&candidate_span))?
        .expect("fragment-reach candidate should produce an aligned window extent");
    let narrowed_fetch = clamp_fetch_to_window_span(&tile, 200, window_extent, 4)?
        .expect("fragment-reach window extent should produce a clamped fetch interval");

    assert_eq!(narrowed_fetch, Interval::new(18, 24)?);
    Ok(())
}

#[test]
fn candidate_window_extent_and_clamp_widen_monotonically_when_more_candidates_are_added()
-> Result<()> {
    let tile = make_tile(10, 20, 6, 30, 0);
    let windows = indexed_windows(&[(10, 11, 0), (22, 23, 1), (25, 26, 2)]);

    let narrow_window_extent = aligned_window_extent_for_bed_candidates(
        &windows,
        Some(&TileWindowSpan {
            first_idx: 0,
            last_idx_exclusive: 2,
        }),
    )?
    .expect("narrow candidate set should produce a window extent");
    let wide_window_extent = aligned_window_extent_for_bed_candidates(
        &windows,
        Some(&TileWindowSpan {
            first_idx: 0,
            last_idx_exclusive: 3,
        }),
    )?
    .expect("wide candidate set should produce a window extent");

    let narrow_fetch = clamp_fetch_to_window_span(&tile, 200, narrow_window_extent, 4)?
        .expect("narrow window extent should clamp to a fetch interval");
    let wide_fetch = clamp_fetch_to_window_span(&tile, 200, wide_window_extent, 4)?
        .expect("wide window extent should clamp to a fetch interval");

    assert!(wide_fetch.start() <= narrow_fetch.start());
    assert!(wide_fetch.end() >= narrow_fetch.end());
    Ok(())
}

#[test]
fn raw_candidate_window_extent_and_clamp_uses_explicit_halo_when_it_exceeds_tile_fetch_halo()
-> Result<()> {
    // Human verification status: Verified
    // Manual derivation:
    // - Tile core [10,20), fetch band [7,26).
    // - Raw left-only candidate window [8,9) sits 2 bp left of the core.
    // - Tile-carried fetch halos are 3 bp on the left and 6 bp on the right, both smaller than
    //   the explicit clamp halo 14.
    // - `clamp_fetch_to_window_span(...)` therefore uses the explicit halo on both sides.
    // - That yields [8-14, 9+14) = [0,23), then clamps to tile.fetch -> [7,23).
    let tile = make_tile(10, 20, 7, 26, 0);
    let windows = indexed_windows(&[(8, 9, 0)]);
    let candidate_span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 1,
    };

    let window_extent = aligned_window_extent_for_bed_candidates(&windows, Some(&candidate_span))?
        .expect("raw left-only candidate should produce an aligned window extent");

    // The single window
    assert_eq!(window_extent, Interval::new(8, 9)?);

    let narrowed_fetch = clamp_fetch_to_window_span(&tile, 200, window_extent, 14)?
        .expect("raw left-only candidate should produce a narrowed fetch interval");

    assert_eq!(narrowed_fetch, Interval::new(7, 23)?);
    assert!(narrowed_fetch.end() < tile.fetch_end() as u64);
    Ok(())
}

#[test]
fn raw_candidate_window_extent_and_clamp_keeps_larger_tile_fetch_halo_when_it_exceeds_explicit_halo()
-> Result<()> {
    // Human verification status: Verified
    // Manual derivation:
    // - Tile core [10,20), fetch band [6,40).
    // - Raw left-only candidate window [8,9) sits 2 bp left of the core.
    // - The explicit clamp halo is 14, but the tile already carries a larger right fetch halo:
    //   20 bp on the right versus only 4 bp on the left.
    // - `clamp_fetch_to_window_span(...)` uses the larger of tile halo and explicit halo on each
    //   side, so it uses 14 on the left and 20 on the right.
    // - That yields [8-14, 9+20) = [0,29), then clamps to tile.fetch -> [6,29).
    let tile = make_tile(10, 20, 6, 40, 0);
    let windows = indexed_windows(&[(8, 9, 0)]);
    let candidate_span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 1,
    };

    let window_extent = aligned_window_extent_for_bed_candidates(&windows, Some(&candidate_span))?
        .expect("raw left-only candidate should produce an aligned window extent");

    assert_eq!(window_extent, Interval::new(8, 9)?);

    let narrowed_fetch = clamp_fetch_to_window_span(&tile, 200, window_extent, 14)?
        .expect("raw left-only candidate should produce a narrowed fetch interval");

    assert_eq!(narrowed_fetch, Interval::new(6, 29)?);
    assert!(narrowed_fetch.end() < tile.fetch_end() as u64);
    Ok(())
}

#[test]
fn raw_candidate_window_extent_and_clamp_narrows_far_right_raw_only_fetch() -> Result<()> {
    // Manual derivation:
    // - Tile core [10,20), fetch band [6,40).
    // - Raw far-right candidate window [32,33) is outside the core but still reachable.
    // - Using the same raw reach halo 14 gives [32-14, 33+14) = [18,47), clamped to [18,40).
    let tile = make_tile(10, 20, 6, 40, 0);
    let windows = indexed_windows(&[(32, 33, 0)]);
    let candidate_span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 1,
    };

    let window_extent = aligned_window_extent_for_bed_candidates(&windows, Some(&candidate_span))?
        .expect("raw far-right candidate should produce an aligned window extent");
    let narrowed_fetch = clamp_fetch_to_window_span(&tile, 200, window_extent, 14)?
        .expect("raw far-right candidate should produce a narrowed fetch interval");

    assert_eq!(narrowed_fetch, Interval::new(18, 40)?);
    assert!(narrowed_fetch.start() > tile.fetch_start() as u64);
    Ok(())
}

#[test]
fn raw_candidate_window_extent_and_clamp_matches_full_tile_when_candidates_span_most_of_it()
-> Result<()> {
    let tile = make_tile(10, 20, 6, 40, 0);
    let windows = indexed_windows(&[(10, 11, 0), (32, 33, 1)]);
    let candidate_span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 2,
    };

    let window_extent = aligned_window_extent_for_bed_candidates(&windows, Some(&candidate_span))?
        .expect("mixed raw candidates should produce an aligned window extent");
    let narrowed_fetch = clamp_fetch_to_window_span(&tile, 200, window_extent, 14)?
        .expect("mixed raw candidates should produce a narrowed fetch interval");
    let full_fetch =
        full_tile_fetch_span(&tile, 200)?.expect("full-tile helper should return tile.fetch");

    assert_eq!(narrowed_fetch, Interval::new(6, 40)?);
    assert_eq!(narrowed_fetch, full_fetch);
    Ok(())
}
