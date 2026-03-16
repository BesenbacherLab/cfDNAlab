mod tests {
    use cfdnalab::commands::fcoverage::tiling::adapt_fetch_to_extreme_windows;
    use cfdnalab::shared::bam::Contigs;
    use cfdnalab::shared::interval::IndexedInterval;
    use cfdnalab::shared::tiled_run::{
        Tile, TileMode, TileWindowSpan, build_tiles, clamp_fetch_to_window_span,
        precompute_tile_window_spans, tile_window_min_max,
    };
    use fxhash::FxHashMap;
    use std::path::PathBuf;

    fn indexed_windows(entries: &[(u64, u64, u64)]) -> Vec<IndexedInterval<u64>> {
        entries
            .iter()
            .map(|&(start, end, original_index)| {
                IndexedInterval::new(start, end, original_index)
                    .expect("test windows should be valid non-empty intervals")
            })
            .collect()
    }

    fn make_tile(
        core_start: u32,
        core_end: u32,
        fetch_start: u32,
        fetch_end: u32,
        index: u32,
    ) -> Tile {
        Tile::new(
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
    fn parse_tile_index_basic() {
        use cfdnalab::shared::tiled_run::parse_tile_index;
        assert_eq!(parse_tile_index("coverage.pos.chr1.12.tsv"), Some(12));
        assert_eq!(
            parse_tile_index("coverage.pos.chr10.000123.bedgraph.zst"),
            Some(123)
        );
        assert_eq!(
            parse_tile_index("coverage.part.chrX.7.part.tsv.zst"),
            Some(7)
        );
        assert_eq!(parse_tile_index("weird.noindex.zst"), None);
    }

    #[test]
    fn clamp_fetch_respects_halo_and_chrom() {
        let tile = make_tile(100, 150, 80, 170, 0);
        // Windows span 90..160; halos are 20 left, 20 right; chrom len 155
        let res = clamp_fetch_to_window_span(&tile, 155, 90, 160, 0).unwrap();
        // Left: 90 - 20 = 70, clamped to fetch_start 80 => 80
        // Right: 160 + 20 = 180, clamped to fetch_end 170 and chrom_len 155 => 155
        assert_eq!(res, (80, 155));
    }

    #[test]
    fn clamp_fetch_returns_none_on_empty_span() {
        let tile = make_tile(100, 150, 80, 170, 0);
        // min_ws >= max_we should return None
        assert!(clamp_fetch_to_window_span(&tile, 200, 120, 120, 0).is_none());
    }

    #[test]
    fn clamp_fetch_clamps_to_fetch_start_when_windows_left_of_tile() {
        let tile = make_tile(100, 150, 90, 200, 0);
        // Windows far left. Even after adding halos the span ends before fetch_start,
        // so there is nothing to fetch and the range is discarded
        let res = clamp_fetch_to_window_span(&tile, 300, 10, 20, 0);
        assert!(res.is_none());
    }

    #[test]
    fn clamp_fetch_clamps_to_chrom_when_windows_right_of_chrom() {
        let tile = make_tile(50, 70, 40, 120, 0);
        // Windows extend beyond chrom_len=100
        let res = clamp_fetch_to_window_span(&tile, 100, 80, 150, 0).unwrap();
        // Left bound follows the nearest window minus halo (80-10) = 70, not the original fetch_start,
        // because there is no window support between 40 and 70; right clamp hits chrom_len=100
        assert_eq!(res, (70, 100));
    }

    #[test]
    fn clamp_fetch_returns_none_when_windows_right_of_tile() {
        let tile = make_tile(100, 150, 90, 200, 0);
        // Windows sit to the right of the tile; even after halo expansion the span begins at 210-10=200,
        // matching fetch_end, so start >= end and the fetch range is discarded
        let res = clamp_fetch_to_window_span(&tile, 230, 210, 220, 0);
        assert!(res.is_none());
    }

    #[test]
    fn tile_window_min_max_returns_extremes() {
        let tile = make_tile(50, 150, 40, 160, 0);
        let windows = indexed_windows(&[(0, 40, 0), (40, 60, 1), (120, 200, 2), (300, 400, 3)]);
        let span = TileWindowSpan {
            first_idx: 1,
            last_idx_exclusive: 3,
        };
        let (min_ws, max_we) = tile_window_min_max(&windows, &tile, Some(&span)).unwrap();
        assert_eq!((min_ws, max_we), (40, 200));
    }

    #[test]
    fn precompute_tile_window_spans_filters_left_windows() {
        let tiles = vec![
            make_tile(100, 150, 90, 170, 0),
            make_tile(150, 200, 140, 220, 1),
        ];
        let windows = indexed_windows(&[(50, 80, 0), (80, 120, 1), (140, 160, 2), (200, 240, 3)]);
        let spans = precompute_tile_window_spans(&tiles, |_| windows.as_slice(), 0, 0);
        let span0 = spans[0].unwrap();
        let span1 = spans[1].unwrap();
        // For tile0, window 0 ends before core start (100) and should be dropped; expect windows 1..3
        assert_eq!(span0.first_idx, 1);
        assert_eq!(span0.last_idx_exclusive, 3);
        // For tile1, window 1 ends before core start (150) and should be dropped; expect windows 2..3
        // because window 3 starts at the core end and does not overlap the half-open core
        // Halos are zero here, so filtering and overlap are based on the core interval only
        assert_eq!(span1.first_idx, 2);
        assert_eq!(span1.last_idx_exclusive, 3);
    }

    #[test]
    fn build_tiles_respects_chrom_end_and_halo() {
        let mut contigs = FxHashMap::default();
        contigs.insert("chr1".to_string(), (0, 95u32));
        let contigs = Contigs { contigs };
        let (tiles, aligned) =
            build_tiles(&vec!["chr1".to_string()], &contigs, 40, 10, None).unwrap();
        assert!(!aligned);
        // Expect 3 tiles: cores [0,40), [40,80), [80,95); halos extend but clamp to chrom len
        assert_eq!(tiles.len(), 3);
        assert_eq!(tiles[2].core_start(), 80);
        assert_eq!(tiles[2].core_end(), 95);
        // TODO: Fetch ends should of course not extend outside the chromosome end, we must fix that if its the case
        assert!(tiles[2].fetch_end() == 95);
        assert!(tiles[2].fetch_start() == 70);
    }

    #[test]
    fn build_tiles_clamps_left_halo_and_zero_halo_matches_core() {
        let mut contigs = FxHashMap::default();
        contigs.insert("chr1".to_string(), (0, 50u32));
        contigs.insert("chr2".to_string(), (1, 30u32));
        let contigs = Contigs { contigs };

        let (tiles, aligned) = build_tiles(
            &vec!["chr1".to_string(), "chr2".to_string()],
            &contigs,
            30,
            20,
            None,
        )
        .unwrap();
        assert!(!aligned);
        // First tile halo would go negative (0-halo), fetch_start is clamped at chromosome start
        assert_eq!(tiles[0].core_start(), 0);
        assert_eq!(tiles[0].fetch_start(), 0);
        assert_eq!(tiles[0].fetch_end(), 50);

        // Zero halo keeps fetch identical to the core on the second chromosome
        let (tiles_zero_halo, aligned_zero) =
            build_tiles(&vec!["chr2".to_string()], &contigs, 15, 0, None).unwrap();
        assert!(!aligned_zero);
        for t in tiles_zero_halo {
            assert_eq!(t.fetch_start(), t.core_start());
            assert_eq!(t.fetch_end(), t.core_end());
        }
    }

    #[test]
    fn precompute_tile_window_spans_expands_with_halos() {
        let tiles = vec![
            make_tile(100, 150, 80, 190, 0),
            make_tile(150, 200, 130, 230, 1),
        ];
        // Window 1 starts inside the core, window 2 starts inside the right halo of tile 0,
        // window 3 starts inside the right halo of tile 1
        let windows = indexed_windows(&[
            (60, 82, 0),
            (110, 120, 1),
            (160, 170, 2),
            (210, 220, 3),
            (240, 250, 4),
        ]);
        // Left halo: 10, Right halo: 15
        let spans = precompute_tile_window_spans(&tiles, |_| windows.as_slice(), 10, 15);
        let span0 = spans[0].unwrap();
        let span1 = spans[1].unwrap();
        // Left halo removes window 0 whose end is before 90, right halo keeps window 2
        assert_eq!(span0.first_idx, 1);
        assert_eq!(span0.last_idx_exclusive, 3);
        // Next tile drops windows ending before 140 and keeps the right-halo window 3
        assert_eq!(span1.first_idx, 2);
        assert_eq!(span1.last_idx_exclusive, 4);
    }

    #[test]
    fn tile_window_min_max_returns_none_for_empty_span() {
        let tile = make_tile(10, 20, 0, 30, 0);
        let windows: Vec<IndexedInterval<u64>> = Vec::new();
        let span = TileWindowSpan {
            first_idx: 0,
            last_idx_exclusive: 0,
        };
        assert!(tile_window_min_max(&windows, &tile, Some(&span)).is_none());
    }

    #[test]
    fn clamp_fetch_returns_none_when_clamping_collapses_span() {
        let tile = make_tile(10, 20, 0, 100, 0);
        // Window span sits right of the chromosome, clamping pulls end left of start
        let res = clamp_fetch_to_window_span(&tile, 120, 150, 160, 0);
        assert!(res.is_none());
    }

    #[test]
    fn clamp_fetch_uses_explicit_halo_when_tile_has_no_inferred_halo() {
        let tile = make_tile(0, 200, 0, 200, 0);

        // This tile has no inferred right halo because the core already reaches the chromosome end.
        // The old logic therefore narrowed a window span of [0, 40) to exactly [0, 40), which is
        // too short to reconstruct fragments extending beyond the window. The explicit halo must
        // widen the fetch to [0, 60) while still staying inside the original tile fetch band.
        let narrowed_fetch = clamp_fetch_to_window_span(&tile, 200, 0, 40, 20).unwrap();

        assert_eq!(narrowed_fetch, (0, 60));
    }

    #[test]
    fn adapt_fetch_keeps_fragment_context_for_bed_aggregate_tiles() {
        let tile = make_tile(0, 200, 0, 200, 0);
        let windows = indexed_windows(&[(0, 40, 0)]);
        let mode = TileMode::AggregatesByBed {
            windows: windows.as_slice(),
            masked: false,
            partials_out: PathBuf::from("partials.tsv.zst"),
            cross_idx_out: PathBuf::from("cross.tsv.zst"),
        };

        // The overlapping BED window itself ends at 40, but fragment reconstruction for this mode
        // still needs an extra halo. With halo_bp=20, the narrowed fetch must keep [0, 60).
        let narrowed_fetch = adapt_fetch_to_extreme_windows(&tile, None, &mode, 200, 20).unwrap();

        assert_eq!(narrowed_fetch, (0, 60));
    }
}
