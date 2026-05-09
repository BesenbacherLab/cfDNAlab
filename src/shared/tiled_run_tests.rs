use super::*;

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
fn temp_dir_guard_removes_non_empty_directory_on_drop() -> anyhow::Result<()> {
    let base_dir = tempfile::TempDir::new()?;
    let guarded_path = {
        let guard = TempDirGuard::new(base_dir.path(), "guard_drop")?;
        let guarded_path = guard.path().to_path_buf();
        std::fs::write(guarded_path.join("tile.tmp"), b"temporary tile payload")?;
        assert!(guarded_path.exists());
        guarded_path
    };

    assert!(!guarded_path.exists());
    Ok(())
}

#[test]
fn temp_dir_guard_remove_is_idempotent() -> anyhow::Result<()> {
    let base_dir = tempfile::TempDir::new()?;
    let mut guard = TempDirGuard::new(base_dir.path(), "guard_remove")?;
    let guarded_path = guard.path().to_path_buf();

    guard.remove()?;
    guard.remove()?;

    assert!(!guarded_path.exists());
    Ok(())
}

#[test]
fn precompute_tile_window_spans_keeps_left_halo_only_windows_for_raw_end_reach() {
    // Raw clipping can move the counted left endpoint left of the aligned fragment start.
    // For a tile core [10,20) and max-soft-clip reach of 2 bp, a BED window [8,9) must stay
    // visible even though it does not overlap the core itself.
    let tiles = vec![make_tile(10, 20, 6, 24, 0)];
    let windows = indexed_windows(&[(8, 9, 0), (10, 11, 1)]);

    let spans = precompute_tile_window_spans(&tiles, |_| windows.as_slice(), 2, 4);
    let span = spans[0].expect("left raw halo window should keep a candidate span");

    assert_eq!(span.first_idx, 0);
    assert_eq!(span.last_idx_exclusive, 2);
}

#[test]
fn precompute_tile_window_spans_excludes_windows_too_far_left_for_raw_end_reach() {
    // Manual derivation:
    // - Tile core is [10,20).
    // - Raw left reach is 2 bp, so the left candidate bound is 8.
    // - Window [7,8) ends exactly at that bound and is therefore too far left.
    // - Window [8,9) still overlaps the reachable raw-left halo and must stay.
    let tiles = vec![make_tile(10, 20, 6, 24, 0)];
    let windows = indexed_windows(&[(7, 8, 0), (8, 9, 1), (10, 11, 2)]);

    let spans = precompute_tile_window_spans(&tiles, |_| windows.as_slice(), 2, 4);
    let span = spans[0].expect("reachable windows should keep a candidate span");

    assert_eq!(span.first_idx, 1);
    assert_eq!(span.last_idx_exclusive, 3);
}

#[test]
fn precompute_tile_window_spans_keeps_far_right_halo_windows_for_raw_end_reach() {
    // If a fragment starts inside core [10,20), aligned length is at most 4 bp, and the right end
    // may extend 10 bp farther in raw mode, then windows starting before 20 + 4 + 10 = 34 must
    // stay visible. A BED window [32,33) is therefore relevant even though it sits far outside
    // the tile core and beyond the aligned-length-only right halo.
    let tiles = vec![make_tile(10, 20, 6, 24, 0)];
    let windows = indexed_windows(&[(10, 11, 0), (32, 33, 1)]);

    let spans = precompute_tile_window_spans(&tiles, |_| windows.as_slice(), 0, 14);
    let span = spans[0].expect("far-right raw halo window should keep a candidate span");

    assert_eq!(span.first_idx, 0);
    assert_eq!(span.last_idx_exclusive, 2);
}

#[test]
fn precompute_tile_window_spans_excludes_windows_too_far_right_for_raw_end_reach() {
    // Manual derivation:
    // - Tile core is [10,20).
    // - Right raw reach is 14 bp, so candidate windows may start before 34.
    // - Window [33,34) is still reachable and must stay.
    // - Window [34,35) starts exactly at the exclusive bound and must be dropped.
    let tiles = vec![make_tile(10, 20, 6, 24, 0)];
    let windows = indexed_windows(&[(10, 11, 0), (33, 34, 1), (34, 35, 2)]);

    let spans = precompute_tile_window_spans(&tiles, |_| windows.as_slice(), 0, 14);
    let span = spans[0].expect("reachable windows should keep a candidate span");

    assert_eq!(span.first_idx, 0);
    assert_eq!(span.last_idx_exclusive, 2);
}

#[test]
fn precompute_tile_window_spans_has_no_left_reach_for_aligned_fragment_models() {
    // Manual derivation:
    // - For aligned fragment-owned models like `lengths`, fragments start inside the core and have
    //   no left reach beyond that aligned start.
    // - With tile core [10,20) and max fragment length 4, a left halo of 0 means [9,10) must be
    //   dropped while the right reachable window [22,23) stays.
    let tiles = vec![make_tile(10, 20, 6, 24, 0)];
    let windows = indexed_windows(&[(9, 10, 0), (10, 11, 1), (22, 23, 2)]);

    let spans = precompute_tile_window_spans(&tiles, |_| windows.as_slice(), 0, 4);
    let span = spans[0].expect("aligned fragment reach should keep the core and right windows");

    assert_eq!(span.first_idx, 1);
    assert_eq!(span.last_idx_exclusive, 3);
}

#[test]
fn precompute_tile_window_spans_preserves_or_extends_span_in_example_when_right_reach_grows() {
    // Human verification status: Verified
    // Manual derivation:
    // - Tile core is [10,20). The fetch band is wider, [4,30), because fetch halo is an aligned
    //   BAM-reading concern carried by the tile.
    // - Candidate selection here does not use that fetch halo. It uses only the tile core plus
    //   the explicit fragment-reach halo passed to `precompute_tile_window_spans(...)`.
    // - With right reach 4, the candidate region is [10,24), so [10,11) and [23,24) stay while
    //   [25,26) and [28,29) are still too far right.
    // - With right reach 8, the candidate region is [10,28), so [25,26) becomes newly reachable
    //   while [28,29) still starts exactly at the exclusive bound and must stay out.
    // - This is one concrete monotonicity example. A wider right reach should preserve the old
    //   candidates and may add more on the right, but this fixture does not prove that invariant
    //   for all possible inputs.
    let tiles = vec![make_tile(10, 20, 4, 30, 0)];
    let windows = indexed_windows(&[(10, 11, 0), (23, 24, 1), (25, 26, 2), (28, 29, 3)]);

    let narrow = precompute_tile_window_spans(&tiles, |_| windows.as_slice(), 0, 4);
    let wide = precompute_tile_window_spans(&tiles, |_| windows.as_slice(), 0, 8);
    let narrow_span = narrow[0].expect("narrow reach should keep a candidate span");
    let wide_span = wide[0].expect("wide reach should keep a candidate span");

    // General asserts (conceptually important though technically redundant)
    assert!(
        wide_span.first_idx <= narrow_span.first_idx,
        "wider right reach must not move the first candidate index to the right"
    );
    assert!(
        wide_span.last_idx_exclusive >= narrow_span.last_idx_exclusive,
        "wider right reach must not shrink the candidate span"
    );

    // Specific asserts
    assert_eq!(narrow_span.first_idx, 0);
    assert_eq!(narrow_span.last_idx_exclusive, 2);
    assert_eq!(wide_span.first_idx, 0);
    assert_eq!(wide_span.last_idx_exclusive, 3);
}

#[test]
fn precompute_tile_window_spans_keeps_boundary_crossing_windows_for_both_neighboring_tiles() {
    // Manual derivation:
    // - Two neighboring tile cores are [10,20) and [20,30).
    // - Window [12,13) lies fully inside the first tile only.
    // - Window [18,22) crosses the shared boundary and should stay visible to both tiles.
    // - Window [22,23) lies fully inside the second tile only.
    // - With zero left/right reach, the candidate spans should therefore be:
    //   tile 0 -> windows 0 and 1
    //   tile 1 -> windows 1 and 2
    let tiles = vec![make_tile(10, 20, 10, 20, 0), make_tile(20, 30, 20, 30, 1)];
    let windows = indexed_windows(&[(12, 13, 0), (18, 22, 1), (22, 23, 2)]);

    let spans = precompute_tile_window_spans(&tiles, |_| windows.as_slice(), 0, 0);
    let first_span = spans[0].expect("first tile should keep a candidate span");
    let second_span = spans[1].expect("second tile should keep a candidate span");

    assert_eq!(first_span.first_idx, 0);
    assert_eq!(first_span.last_idx_exclusive, 2);
    assert_eq!(second_span.first_idx, 1);
    assert_eq!(second_span.last_idx_exclusive, 3);
}

/* candidate_window_span for core-overlap */

#[test]
fn candidate_window_span_for_tile_core_overlap_keeps_fully_internal_window() {
    // Human verification status: Verified
    let tile = make_tile(10, 20, 6, 24, 0);
    let windows = indexed_windows(&[(12, 13, 0)]);

    let span = candidate_window_span_for_tile_core_overlap(&windows, &tile)
        .expect("internal core-overlap window should produce a span");

    assert_eq!(span.first_idx, 0);
    assert_eq!(span.last_idx_exclusive, 1);
}

#[test]
fn candidate_window_span_for_tile_core_overlap_keeps_left_boundary_crossing_window() {
    // Human verification status: Verified
    let tile = make_tile(10, 20, 6, 24, 0);
    let windows = indexed_windows(&[(8, 12, 0)]);

    let span = candidate_window_span_for_tile_core_overlap(&windows, &tile)
        .expect("left boundary-crossing window should produce a span");

    assert_eq!(span.first_idx, 0);
    assert_eq!(span.last_idx_exclusive, 1);
}

#[test]
fn candidate_window_span_for_tile_core_overlap_keeps_right_boundary_crossing_window() {
    // Human verification status: Verified
    let tile = make_tile(10, 20, 6, 24, 0);
    let windows = indexed_windows(&[(18, 22, 0)]);

    let span = candidate_window_span_for_tile_core_overlap(&windows, &tile)
        .expect("right boundary-crossing window should produce a span");

    assert_eq!(span.first_idx, 0);
    assert_eq!(span.last_idx_exclusive, 1);
}

#[test]
fn candidate_window_span_for_tile_core_overlap_drops_window_ending_at_core_start() {
    // Human verification status: Verified
    let tile = make_tile(10, 20, 6, 24, 0);
    let windows = indexed_windows(&[(9, 10, 0)]);

    let span = candidate_window_span_for_tile_core_overlap(&windows, &tile);

    assert!(span.is_none());
}

#[test]
fn candidate_window_span_for_tile_core_overlap_drops_window_starting_at_core_end() {
    // Human verification status: Verified
    let tile = make_tile(10, 20, 6, 24, 0);
    let windows = indexed_windows(&[(20, 21, 0)]);

    let span = candidate_window_span_for_tile_core_overlap(&windows, &tile);

    assert!(span.is_none());
}

#[test]
fn candidate_window_span_for_tile_core_overlap_keeps_only_true_core_windows_in_mixed_case() {
    // Human verification status: Verified
    let tile = make_tile(10, 20, 6, 24, 0);
    let windows = indexed_windows(&[(8, 9, 0), (10, 11, 1), (18, 22, 2), (22, 23, 3)]);

    let span = candidate_window_span_for_tile_core_overlap(&windows, &tile)
        .expect("mixed core-overlap case should keep the true core windows");

    assert_eq!(span.first_idx, 1);
    assert_eq!(span.last_idx_exclusive, 3);
}

/* candidate_window_span for fragment reach

Note that it's the left and right reach values that determine these tests, not the fetch halos.
*/

#[test]
fn candidate_window_span_for_tile_fragment_reach_has_no_left_reach_in_aligned_models() {
    // Human verification status: Verified
    let tile = make_tile(10, 20, 6, 24, 0);
    let windows = indexed_windows(&[(9, 10, 0), (10, 11, 1), (22, 23, 2)]);

    // The relevant window span reach "halo" are the ones specified here (0 and 4)
    let span = candidate_window_span_for_tile_fragment_reach(&windows, &tile, 0, 4)
        .expect("aligned fragment reach should keep a span");

    assert_eq!(span.first_idx, 1);
    assert_eq!(span.last_idx_exclusive, 3);
}

#[test]
fn candidate_window_span_for_tile_fragment_reach_keeps_one_right_halo_only_window() {
    // Human verification status: Verified
    let tile = make_tile(10, 20, 6, 24, 0);
    let windows = indexed_windows(&[(22, 23, 0)]);

    let span = candidate_window_span_for_tile_fragment_reach(&windows, &tile, 0, 4)
        .expect("right halo-only fragment-reach window should keep a span");

    assert_eq!(span.first_idx, 0);
    assert_eq!(span.last_idx_exclusive, 1);
}

#[test]
fn candidate_window_span_for_tile_fragment_reach_drops_window_at_exclusive_right_bound() {
    // Human verification status: Verified
    let tile = make_tile(10, 20, 6, 28, 0);
    let windows = indexed_windows(&[(24, 25, 0)]);

    let span = candidate_window_span_for_tile_fragment_reach(&windows, &tile, 0, 4);

    assert!(span.is_none());
}

#[test]
fn candidate_window_span_for_tile_fragment_reach_keeps_core_and_right_halo_windows_together() {
    // Human verification status: Verified
    let tile = make_tile(10, 20, 6, 21, 0);
    let windows = indexed_windows(&[(10, 11, 0), (22, 23, 1)]);

    let span = candidate_window_span_for_tile_fragment_reach(&windows, &tile, 0, 4)
        .expect("core and right halo windows should share one fragment-reach span");

    assert_eq!(span.first_idx, 0);
    assert_eq!(span.last_idx_exclusive, 2);
}

#[test]
fn candidate_window_span_for_tile_fragment_reach_drops_left_non_overlap_in_mixed_case() {
    // Human verification status: Verified
    let tile = make_tile(10, 20, 6, 24, 0);
    let windows = indexed_windows(&[(8, 9, 0), (10, 11, 1), (22, 23, 2)]);

    let span = candidate_window_span_for_tile_fragment_reach(&windows, &tile, 0, 4)
        .expect("mixed fragment-reach case should keep a span");

    assert_eq!(span.first_idx, 1);
    assert_eq!(span.last_idx_exclusive, 3);
}

#[test]
fn candidate_window_span_for_tile_fragment_reach_selects_downstream_window_also_selected_in_next_tile()
 {
    // Human verification status: Verified
    // Manual derivation:
    // - Tile 0 core is [10,20) with right reach 4, so its candidate region is [10,24).
    // - Window [22,23) is therefore reachable from tile 0, while [25,28) and [38,42) are too far
    //   right for that tile.
    // - Tile 1 core is [20,30) with the same right reach, so its candidate region is [20,34).
    // - Windows [22,23) and [25,28) are both relevant for tile 1, while [38,42) starts beyond
    //   the exclusive reachable bound.
    // - Candidate selection is about relevance, not exclusive window ownership, so [22,23) should
    //   appear in both spans.
    let tiles = vec![make_tile(10, 20, 6, 24, 0), make_tile(20, 30, 16, 34, 1)];
    let windows = indexed_windows(&[(22, 23, 0), (25, 28, 1), (38, 42, 2)]);

    let first_span = candidate_window_span_for_tile_fragment_reach(&windows, &tiles[0], 0, 4)
        .expect("window in downstream tile within fragment reach should be included ");
    let second_span = candidate_window_span_for_tile_fragment_reach(&windows, &tiles[1], 0, 4)
        .expect("on-tile windows should be found");

    assert_eq!(first_span.first_idx, 0);
    assert_eq!(first_span.last_idx_exclusive, 1);
    assert_eq!(second_span.first_idx, 0);
    assert_eq!(second_span.last_idx_exclusive, 2);
}

#[test]
fn candidate_window_span_for_tile_raw_fragment_reach_keeps_left_halo_only_window() {
    // Human verification status: Verified
    let tile = make_tile(10, 20, 10, 24, 0);
    let windows = indexed_windows(&[(8, 9, 0)]);

    let span = candidate_window_span_for_tile_fragment_reach(&windows, &tile, 2, 4)
        .expect("raw left halo-only window should keep a span");

    assert_eq!(span.first_idx, 0);
    assert_eq!(span.last_idx_exclusive, 1);
}

#[test]
fn candidate_window_span_for_tile_raw_fragment_reach_drops_window_at_left_exclusive_bound() {
    // Human verification status: Verified
    let tile = make_tile(10, 20, 6, 24, 0);
    let windows = indexed_windows(&[(7, 8, 0)]);

    let span = candidate_window_span_for_tile_fragment_reach(&windows, &tile, 2, 4);

    assert!(span.is_none());
}

#[test]
fn candidate_window_span_for_tile_raw_fragment_reach_keeps_left_core_and_far_right_windows_together()
 {
    // Human verification status: Verified
    let tile = make_tile(10, 20, 6, 24, 0);
    let windows = indexed_windows(&[(8, 9, 0), (10, 11, 1), (32, 33, 2)]);

    let span = candidate_window_span_for_tile_fragment_reach(&windows, &tile, 2, 14)
        .expect("raw fragment reach should keep left/core/right mixed candidates");

    assert_eq!(span.first_idx, 0);
    assert_eq!(span.last_idx_exclusive, 3);
}

#[test]
fn candidate_window_span_for_tile_raw_fragment_reach_uses_asymmetric_left_and_right_bounds() {
    // Human verification status: Verified
    // Manual derivation:
    // - Tile core is [10,20).
    // - Raw left reach is 2 bp, so [7,8) is too far left and must be dropped.
    // - Raw right reach is 14 bp, so [32,33) is still reachable and must stay.
    // - A symmetric implementation would either keep both or drop both.
    let tile = make_tile(10, 20, 6, 24, 0);
    let windows = indexed_windows(&[(7, 8, 0), (10, 11, 1), (32, 33, 2)]);

    let span = candidate_window_span_for_tile_fragment_reach(&windows, &tile, 2, 14)
        .expect("asymmetric raw fragment reach should keep the core and far-right windows");

    assert_eq!(span.first_idx, 1);
    assert_eq!(span.last_idx_exclusive, 3);
}
