mod tests_blacklist_windowing_2 {
    use super::super::*;
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
    fn prepare_tile_windows_bed_keeps_a_right_halo_only_window_for_fragment_owned_models()
    -> Result<()> {
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
        let prefixes = build_gc_prefixes(&[b'C'; 10]);
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
        let prefixes = build_gc_prefixes(&[b'C'; 10]);
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
    fn prepare_tile_windows_bed_skips_tile_for_empty_window_slice_even_with_nonempty_span()
    -> Result<()> {
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
}

mod tests_blacklist_windowing {
    use crate::commands::cli_common::WindowSpec;
    use crate::commands::gc_bias::counting::{GCCounts, build_gc_prefixes};
    use crate::commands::gc_bias::windows::{
        WindowState, compute_window_stats, set_window_acgt_in_observed_interval,
    };
    use crate::shared::bam::Contigs;
    use crate::shared::bed::Windows;
    use crate::shared::interval::{IndexedInterval, Interval};
    use fxhash::FxHashMap;

    fn build_contigs(entries: &[(&str, u32)]) -> Contigs {
        let mut contigs = FxHashMap::default();
        for (name, length) in entries {
            contigs.insert((*name).to_string(), (0, *length));
        }
        Contigs { contigs }
    }

    fn build_windows(entries: Vec<(u64, u64, u64)>) -> Windows {
        let windows = entries
            .into_iter()
            .map(|(start, end, idx)| {
                IndexedInterval::new(start, end, idx).expect("test windows should be valid")
            })
            .collect();
        Windows::from_sorted(windows)
    }

    fn make_window_state(start: u64, end: u64) -> WindowState {
        let template = GCCounts::new(1, 1, 0, (0, 0)).expect("failed to create template");
        WindowState::new(
            0,
            Interval::new(start, end).expect("test interval should be valid"),
            true,
            &template,
        )
        .expect("failed to build window state")
    }

    mod test_compute_window_stats {
        use super::*;

        #[test]
        fn returns_average_span_and_total_for_bed_windows() {
            // Three windows across two chromosomes: lengths 100, 60, 50
            let window_spec = WindowSpec::Bed("dummy.bed".into());
            let mut map = FxHashMap::default();
            map.insert(
                "chr1".to_string(),
                build_windows(vec![(0, 100, 0), (200, 260, 1)]),
            );
            map.insert("chr2".to_string(), build_windows(vec![(0, 50, 2)]));
            let contigs = build_contigs(&[("chr1", 300), ("chr2", 100)]);
            let chromosomes = vec!["chr1".to_string(), "chr2".to_string()];

            let stats =
                compute_window_stats(&window_spec, Some(&map), &contigs, &chromosomes).unwrap();

            assert!(
                (stats.avg_span - 70.0).abs() < f64::EPSILON,
                "expected average span of 70.0, got {}",
                stats.avg_span
            );
            assert_eq!(
                stats.total_windows, 3,
                "expected three BED windows across both chromosomes"
            );
        }

        #[test]
        fn returns_genome_span_and_count_for_global_windows() {
            // Global mode should sum chromosome lengths regardless of window data
            let window_spec = WindowSpec::Global;
            let contigs = build_contigs(&[("chr1", 1000), ("chr2", 1500)]);
            let chromosomes = vec!["chr1".to_string(), "chr2".to_string()];

            let stats = compute_window_stats(&window_spec, None, &contigs, &chromosomes).unwrap();

            assert!(
                (stats.avg_span - 2500.0).abs() < f64::EPSILON,
                "expected combined span of 2500.0, got {}",
                stats.avg_span
            );
            assert_eq!(
                stats.total_windows, 1,
                "global mode should report one window"
            );
        }

        #[test]
        fn averages_and_counts_fixed_size_windows() {
            // Fixed windows create ceil(len/size) windows per contig: chr1 -> 4, chr2 -> 2
            let window_spec = WindowSpec::Size(300);
            let contigs = build_contigs(&[("chr1", 1000), ("chr2", 500)]);
            let chromosomes = vec!["chr1".to_string(), "chr2".to_string()];

            let stats = compute_window_stats(&window_spec, None, &contigs, &chromosomes).unwrap();

            assert!(
                (stats.avg_span - 250.0).abs() < f64::EPSILON,
                "expected average window span of 250.0, got {}",
                stats.avg_span
            );
            assert_eq!(
                stats.total_windows, 6,
                "expected four windows on chr1 and two on chr2"
            );
        }
    }

    mod tests_set_window_acgt_in_observed_interval {
        use super::*;

        #[test]
        fn counts_only_the_supplied_observed_subinterval() {
            // The window spans 0-6, but we only observe the middle 2-5 segment.
            // Sequence A C N G T A
            //               ^^^^^
            // observed 2-5 = N G T -> 2 ACGT bases over length 3
            let seq = b"ACNGTA";
            let prefixes = build_gc_prefixes(seq);
            let mut window = make_window_state(0, 6);
            let observed_interval =
                Interval::new(2, 5).expect("test observed interval should be valid");

            set_window_acgt_in_observed_interval(
                &mut window,
                &prefixes,
                observed_interval,
                Interval::new(0, seq.len() as u64).expect("test sequence interval should be valid"),
            )
            .unwrap();

            assert_eq!(
                window.counts.num_acgt_out_of,
                (2, 3),
                "expected support to be counted from the observed subinterval, not the full window"
            );
        }

        #[test]
        fn counts_acgt_when_window_overlaps_sequence() {
            // Sequence AACGTN has four ACGT bases inside window 1-5 (ACGT)
            let seq = b"AACGTN";
            let prefixes = build_gc_prefixes(seq);
            let mut window = make_window_state(1, 5); // covers "ACGT"
            let observed_interval = window.interval;

            set_window_acgt_in_observed_interval(
                &mut window,
                &prefixes,
                observed_interval,
                Interval::new(0, seq.len() as u64).expect("test sequence interval should be valid"),
            )
            .unwrap();

            assert_eq!(
                window.counts.num_acgt_out_of,
                (4, 4),
                "expected 4 ACGT bases over length 4"
            );
        }

        #[test]
        fn counts_acgt_detect_n() {
            // Sequence AANGTN has three ACGT bases inside window 1-5 (ANGT)
            let seq = b"AANGTN";
            let prefixes = build_gc_prefixes(seq);
            let mut window = make_window_state(1, 5); // covers "ANGT"
            let observed_interval = window.interval;

            set_window_acgt_in_observed_interval(
                &mut window,
                &prefixes,
                observed_interval,
                Interval::new(0, seq.len() as u64).expect("test sequence interval should be valid"),
            )
            .unwrap();

            assert_eq!(
                window.counts.num_acgt_out_of,
                (3, 4),
                "expected 3 ACGT bases over length 4"
            );
        }

        #[test]
        fn errors_when_window_has_no_overlap() {
            // Window 10-12 sits entirely outside the available sequence 0-4 so should error
            let seq = b"ACGT";
            let prefixes = build_gc_prefixes(seq);
            let mut window = make_window_state(10, 12);
            let observed_interval = window.interval;

            let err = set_window_acgt_in_observed_interval(
                &mut window,
                &prefixes,
                observed_interval,
                Interval::new(0, seq.len() as u64).expect("test sequence interval should be valid"),
            )
            .expect_err("expected an overlap error");

            let msg = format!("{err}");
            assert!(
                msg.contains("does not overlap sequence"),
                "unexpected error message: {msg}"
            );
        }

        #[test]
        fn errors_when_window_exceeds_prefix_bounds() {
            // Window 0-6 extends past the computed prefix bounds of a 0-6 sequence slice, expect bounds error
            let seq = b"ACGT";
            let prefixes = build_gc_prefixes(seq);
            let mut window = make_window_state(0, 6);
            let observed_interval = window.interval;

            let err = set_window_acgt_in_observed_interval(
                &mut window,
                &prefixes,
                observed_interval,
                Interval::new(0, 6).expect("test sequence interval should be valid"),
            )
            .expect_err("expected a prefix bounds error");

            let msg = format!("{err}");
            let chain: Vec<String> = err.chain().map(|cause| cause.to_string()).collect();
            assert!(
                msg.contains("counting ACGT support for observed interval [0, 6)"),
                "expected high-level window context, got: {msg}"
            );
            assert!(
                chain
                    .iter()
                    .any(|cause| cause.contains("ACGT interval [0, 6) out of bounds")),
                "expected low-level prefix-bounds cause in error chain, got: {chain:?}"
            );
        }
    }

    mod tests_prepare_tile_windows {
        use crate::{
            commands::{
                cli_common::WindowSpec,
                gc_bias::{counting::GCCounts, windows::prepare_tile_windows},
            },
            shared::{
                interval::IndexedInterval,
                tiled_run::{Tile, precompute_tile_window_spans},
            },
        };
        use anyhow::Result;
        use std::path::PathBuf;

        fn make_template() -> GCCounts {
            GCCounts::new(1, 1, 0, (0, 0)).expect("failed to build template")
        }

        fn make_tile() -> Tile {
            Tile::from_coords("chr1".to_string(), 0, 0, 100, 190, 80, 210)
                .expect("test tile should be valid")
        }

        fn indexed_windows(entries: &[(u64, u64, u64)]) -> Vec<IndexedInterval<u64>> {
            entries
                .iter()
                .map(|&(start, end, original_index)| {
                    IndexedInterval::new(start, end, original_index)
                        .expect("test windows should be valid non-empty intervals")
                })
                .collect()
        }

        #[test]
        fn builds_bed_windows_for_tile_reach() -> Result<()> {
            let template = make_template();
            let tile = make_tile();
            // Manual derivation:
            // - Tile core is [100,190).
            // - `gc_bias` candidate spans use left_halo=0 and right_halo=max_fragment_length.
            // - With max_fragment_length=20, the reachable right bound is 210.
            // - Therefore windows starting before 210 stay, and windows starting at/after 210 drop.
            // - [90,140), [120,180), and [200,240) stay; [210,250) and [220,260) drops.
            let windows = indexed_windows(&[
                (90, 140, 0),
                (120, 180, 1),
                (200, 240, 2),
                (210, 250, 3),
                (220, 260, 4),
            ]);
            let tiles = vec![tile.clone()];
            let spans = precompute_tile_window_spans(&tiles, |_| windows.as_slice(), 0, 20);
            let span = spans[0]
                .as_ref()
                .expect("fragment-reach precompute should keep the first three windows");
            assert_eq!(span.first_idx, 0);
            assert_eq!(span.last_idx_exclusive, 3);

            let prepared = prepare_tile_windows(
                &WindowSpec::Bed(PathBuf::from("dummy.bed")),
                Some(&windows),
                &tile,
                Some(span),
                500,
                &template,
            )?;

            assert!(!prepared.skip_tile);
            assert!(prepared.streaming_buffers.is_none());
            assert_eq!(prepared.windows.len(), 3);
            assert_eq!(prepared.windows[0].idx, 0);
            assert_eq!(prepared.windows[1].idx, 1);
            assert_eq!(prepared.windows[2].idx, 2);
            assert!(!prepared.windows[0].contained);
            assert!(prepared.windows[1].contained);
            assert!(!prepared.windows[2].contained);
            Ok(())
        }

        #[test]
        fn skips_tile_when_no_bed_windows_available() -> Result<()> {
            let template = make_template();
            let tile = make_tile();
            // Empty BED slice should return skip=true so caller can bail out early
            let windows: Vec<IndexedInterval<u64>> = Vec::new();

            let prepared = prepare_tile_windows(
                &WindowSpec::Bed(PathBuf::from("empty.bed")),
                Some(&windows),
                &tile,
                None,
                500,
                &template,
            )?;

            assert!(prepared.skip_tile);
            assert!(prepared.windows.is_empty());
            assert!(prepared.streaming_buffers.is_none());
            Ok(())
        }

        #[test]
        fn prepares_streaming_buffers_for_fixed_windows() -> Result<()> {
            let template = make_template();
            let tile = Tile::from_coords("chr1".to_string(), 0, 0, 250, 450, 230, 470)
                .expect("test tile should be valid");

            // Fixed-size windows use rolling buffers (current and next) instead of per-window allocation
            let prepared =
                prepare_tile_windows(&WindowSpec::Size(200), None, &tile, None, 1000, &template)?;

            let windows = prepared.windows;
            assert!(windows.is_empty());
            assert!(!prepared.skip_tile);

            let (window_bp, current, next) = prepared
                .streaming_buffers
                .expect("expected streaming buffers for fixed-size windows");
            assert_eq!(window_bp, 200);

            assert_eq!(current.idx, 1);
            assert_eq!(current.start(), 200);
            assert_eq!(current.end(), 400);
            assert!(!current.contained);

            let next = next.expect("expected next fixed-size window buffer");
            assert_eq!(next.idx, 2);
            assert_eq!(next.start(), 400);
            assert_eq!(next.end(), 600);
            assert!(!next.contained);
            Ok(())
        }

        #[test]
        fn builds_global_window_for_tile_core() -> Result<()> {
            let template = make_template();
            let tile = make_tile();

            let prepared =
                prepare_tile_windows(&WindowSpec::Global, None, &tile, None, 500, &template)?;

            assert!(!prepared.skip_tile);
            assert!(prepared.streaming_buffers.is_none());
            assert_eq!(prepared.windows.len(), 1);
            let window = &prepared.windows[0];
            assert_eq!(window.idx, 0);
            assert_eq!(window.start(), tile.core_start() as u64);
            assert_eq!(window.end(), tile.core_end() as u64);
            assert!(window.contained);
            Ok(())
        }
    }

    mod tests_gc_bias_window_logic {
        use anyhow::Result;
        use ndarray::{Array1, Array2};
        use tempfile::tempdir;

        use crate::commands::{
            cli_common::{ChromosomeArgs, IOCArgs},
            gc_bias::{
                config::GCConfig, counting::GCCounts, cross_tile_parts::stream_crossing_files,
                gc_bias::process_window,
            },
        };

        fn make_config(tmp: &tempfile::TempDir) -> GCConfig {
            let ioc = IOCArgs {
                bam: tmp.path().join("dummy.bam"),
                output_dir: tmp.path().join("out"),
                n_threads: 1,
            };
            let mut cfg = GCConfig::new(
                ioc,
                tmp.path().join("ref.2bit"),
                tmp.path().join("ref_gc_package.zarr"),
                ChromosomeArgs::default(),
            );
            cfg.set_min_window_acgt_pct(0);
            cfg
        }

        #[test]
        fn scales_window_by_mean_and_acgt_coverage() -> Result<()> {
            // Arrange: One length row (effective length 10 -> 11 GC bins). Only two bins set (2 and 4),
            // so mean = (2+4) / 11 = 0.54545...
            // Scale factor = (1/mean) * (num_acgt/avg_span) = (1/0.54545) * (40/100) = 0.73333...
            let tmp = tempdir()?;
            let cfg = make_config(&tmp);

            let mut counts = GCCounts::new(10, 10, 0, (40, 50))?;
            counts.set(10, 0, 2.0);
            counts.set(10, 1, 4.0);

            let scaled = process_window(counts, &cfg, 100.0)?.expect("window should be retained");

            // Assert
            let c0 = scaled.get(10, 0).unwrap();
            let c1 = scaled.get(10, 1).unwrap();
            assert!((c0 - 1.4666667).abs() < 1e-6);
            assert!((c1 - 2.9333334).abs() < 1e-6);
            Ok(())
        }

        #[test]
        fn drops_window_when_acgt_fraction_below_threshold() -> Result<()> {
            // Arrange: Only 25% of the positions are ACGT, below the 50% threshold
            let tmp = tempdir()?;
            let mut cfg = make_config(&tmp);
            cfg.set_min_window_acgt_pct(50);

            let mut counts = GCCounts::new(10, 10, 0, (5, 20))?;
            counts.set(10, 0, 5.0);

            // Act
            let result = process_window(counts, &cfg, 100.0)?;

            // Assert
            assert!(
                result.is_none(),
                "window with low ACGT should be filtered out"
            );
            Ok(())
        }

        // TODO: Validate this
        #[test]
        fn merges_crossing_files_and_scales_once_per_window() -> Result<()> {
            // Arrange: two crossing chunks for the same window idx=3, counts 2 and 3, acgt 20 and 30.
            // Merged counts=5, num_acgt=50 -> mean=5/11=0.45454..., scale=(1/0.45454)*(50/20)=5.5, final count=27.5.
            let tmp = tempdir()?;
            let cfg = make_config(&tmp);
            let template = GCCounts::new(10, 10, 0, (0, 0))?;
            let counts_len = template.borrow_raw_counts().len();

            let file1 = tmp.path().join("cross.1.npz");
            let mut npz1 = ndarray_npy::NpzWriter::new(std::fs::File::create(&file1)?);
            npz1.add_array("idx", &Array1::from(vec![3u64]))?;
            npz1.add_array("acgt0", &Array1::from(vec![20u64]))?;
            npz1.add_array("acgt1", &Array1::from(vec![20u64]))?;
            let mut counts_arr1 = Array2::zeros((1, counts_len));
            counts_arr1[[0, 0]] = 2.0;
            npz1.add_array("counts", &counts_arr1)?;
            npz1.finish()?;

            let file2 = tmp.path().join("cross.2.npz");
            let mut npz2 = ndarray_npy::NpzWriter::new(std::fs::File::create(&file2)?);
            npz2.add_array("idx", &Array1::from(vec![3u64]))?;
            npz2.add_array("acgt0", &Array1::from(vec![30u64]))?;
            npz2.add_array("acgt1", &Array1::from(vec![30u64]))?;
            let mut counts_arr2 = Array2::zeros((1, counts_len));
            counts_arr2[[0, 0]] = 3.0;
            npz2.add_array("counts", &counts_arr2)?;
            npz2.finish()?;

            let (merged, weight) =
                stream_crossing_files(vec![file1, file2], &template, &cfg, 20.0)?;

            // Assert
            assert_eq!(weight, 1, "one window should contribute once");
            let v = merged.get(10, 0).unwrap();
            assert!((v - 27.5).abs() < 1e-6);
            Ok(())
        }
    }
}
