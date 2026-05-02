#![cfg(feature = "cmd_gc_bias")]

mod tests_stream_helpers_and_finalizer {
    use anyhow::Result;
    use tempfile::tempdir;

    use cfdnalab::commands::{
        cli_common::{ChromosomeArgs, IOCArgs},
        gc_bias::{
            config::GCConfig,
            counting::{GCCounts, build_gc_prefixes},
            gc_bias::finalize_window_buffer,
            windows::{WindowState, fixed_size_window_interval, overlap_length},
        },
    };
    use cfdnalab::shared::interval::Interval;

    fn make_config(tmp: &tempfile::TempDir) -> GCConfig {
        let ioc = IOCArgs {
            bam: tmp.path().join("dummy.bam"),
            output_dir: tmp.path().join("out"),
            n_threads: 1,
        };
        let mut cfg = GCConfig::new(
            ioc,
            tmp.path().join("ref.2bit"),
            tmp.path().join("ref_gc_package.npz"),
            ChromosomeArgs::default(),
        );
        cfg.set_min_window_acgt_pct(0);
        cfg
    }

    #[test]
    fn window_bounds_caps_at_chrom_len() {
        // Human verification status: unverified
        let interval = fixed_size_window_interval(8, 100, 850).expect("interval should be valid");
        assert_eq!(interval.start(), 800);
        assert_eq!(interval.end(), 850);
    }

    #[test]
    fn overlap_length_returns_expected_span() {
        // Human verification status: unverified
        assert_eq!(
            overlap_length(
                Interval::new(0, 10).expect("test interval should be valid"),
                Interval::new(5, 15).expect("test interval should be valid")
            ),
            5
        );
        assert_eq!(
            overlap_length(
                Interval::new(0, 10).expect("test interval should be valid"),
                Interval::new(10, 20).expect("test interval should be valid")
            ),
            0
        );
    }

    #[test]
    fn finalize_window_buffer_scales_and_merges() -> Result<()> {
        // Human verification status: unverified
        let tmp = tempdir()?;
        let cfg = make_config(&tmp);

        let template = GCCounts::new(1, 1, 0, (0, 0))?;
        let mut out = WindowState::new(0, Interval::new(0, 1)?, true, &template)?;
        let mut crossing_parts = Vec::new();

        let seq = vec![b'A', b'A'];
        let prefixes = build_gc_prefixes(&seq);

        let mut buf = WindowState::new(0, Interval::new(0, 2)?, true, &template)?;
        buf.counts.set(1, 0, 2.0);
        buf.has_counts = true;

        finalize_window_buffer(
            &mut buf,
            &prefixes,
            Interval::new(0, 2)?,
            Interval::new(0, 2)?,
            true,
            true,
            &cfg,
            2.0,
            &mut out,
            &mut crossing_parts,
            0,
        )?;

        assert_eq!(out.weight, 1);
        let val = out.counts.get(1, 0).unwrap();
        assert!((val - 2.0).abs() < 1e-6);
        assert!(crossing_parts.is_empty());
        Ok(())
    }

    #[test]
    fn finalize_window_buffer_ignores_empty_placeholder_window_outside_loaded_sequence()
    -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Model the fixed-size streaming "next" placeholder that can exist beyond the loaded
        // sequence span for the current tile:
        // - tile core interval:       [0, 95)
        // - loaded sequence interval: [0, 99)
        //   (the fetch band extends 4 bp beyond the core)
        // - placeholder window:       [100, 200)
        //
        // The buffer is empty (`has_counts = false`), so finalization should no-op instead of
        // trying to compute ACGT support for a window that does not overlap the loaded sequence.
        let tmp = tempdir()?;
        let cfg = make_config(&tmp);

        let template = GCCounts::new(1, 1, 0, (0, 0))?;
        let mut out = WindowState::new(0, Interval::new(0, 1)?, true, &template)?;
        let mut crossing_parts = Vec::new();

        let seq = vec![b'A'; 99];
        let prefixes = build_gc_prefixes(&seq);

        let mut buf = WindowState::new(1, Interval::new(100, 200)?, false, &template)?;
        buf.has_counts = false;

        // Act
        finalize_window_buffer(
            &mut buf,
            &prefixes,
            Interval::new(0, 99)?,
            Interval::new(0, 95)?,
            false,
            true,
            &cfg,
            100.0,
            &mut out,
            &mut crossing_parts,
            0,
        )?;

        // Assert
        assert_eq!(out.weight, 0);
        assert_eq!(out.counts.sum(), 0.0);
        assert!(crossing_parts.is_empty());
        assert_eq!(buf.counts.sum(), 0.0);
        assert!(!buf.has_counts);

        Ok(())
    }

    #[test]
    fn finalize_window_buffer_spills_crossing_support_without_double_counting_fetch_halo()
    -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Model one fixed window [0,100) that crosses a tile boundary at 60, with a 10 bp fetch
        // halo on each tile:
        // - left  tile core  [0,60),  loaded sequence [0,70)
        // - right tile core [60,100), loaded sequence [50,100)
        //
        // The two loaded sequence intervals overlap on [50,70), so computing support from the
        // fetch spans would double count 20 bp when the crossing parts are merged:
        //   wrong support = 70 + 50 = 120
        //
        // The correct support must be counted from the tile-owned segments only:
        //   left  contribution = [0,60)   -> 60 bp
        //   right contribution = [60,100) -> 40 bp
        //   merged support      = 100 bp
        //
        // The right tile is empty (`has_counts = false`), so the merged window still depends on
        // support contributed by an otherwise empty neighboring tile.
        let tmp = tempdir()?;
        let cfg = make_config(&tmp);

        let template = GCCounts::new(1, 1, 0, (0, 0))?;
        let mut out = WindowState::new(999, Interval::new(0, 1)?, true, &template)?;
        let mut crossing_parts = Vec::new();

        let mut left_buf = WindowState::new(0, Interval::new(0, 100)?, false, &template)?;
        left_buf.counts.set(1, 0, 1.0);
        left_buf.has_counts = true;

        let mut right_buf = WindowState::new(0, Interval::new(0, 100)?, false, &template)?;
        right_buf.has_counts = false;

        let left_seq = vec![b'A'; 70];
        let right_seq = vec![b'A'; 50];
        let left_prefixes = build_gc_prefixes(&left_seq);
        let right_prefixes = build_gc_prefixes(&right_seq);

        // Act
        finalize_window_buffer(
            &mut left_buf,
            &left_prefixes,
            Interval::new(0, 70)?,
            Interval::new(0, 60)?,
            false,
            true,
            &cfg,
            100.0,
            &mut out,
            &mut crossing_parts,
            0,
        )?;
        finalize_window_buffer(
            &mut right_buf,
            &right_prefixes,
            Interval::new(50, 100)?,
            Interval::new(60, 100)?,
            false,
            true,
            &cfg,
            100.0,
            &mut out,
            &mut crossing_parts,
            0,
        )?;

        // Assert
        assert_eq!(crossing_parts.len(), 2);
        assert_eq!(crossing_parts[0].idx, 0);
        assert_eq!(crossing_parts[1].idx, 0);
        assert_eq!(crossing_parts[0].counts.num_acgt_out_of, (60, 60));
        assert_eq!(crossing_parts[1].counts.num_acgt_out_of, (40, 40));

        let merged_acgt: (u64, u64) = crossing_parts.iter().fold((0, 0), |acc, part| {
            (
                acc.0 + part.counts.num_acgt_out_of.0,
                acc.1 + part.counts.num_acgt_out_of.1,
            )
        });
        assert_eq!(
            merged_acgt,
            (100, 100),
            "crossing support must sum to the true window span, not double count overlapping fetch halos"
        );

        Ok(())
    }

    #[test]
    fn finalize_window_buffer_spills_counted_next_window_with_zero_owned_support() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Model the fixed-size streaming "next" window [100,200) while the current tile core is
        // still [0,95) and the fetched sequence extends only to 99. Fragments starting in the
        // current tile can still contribute counts to the next window, but this tile owns none of
        // that window's support span.
        let tmp = tempdir()?;
        let cfg = make_config(&tmp);

        let template = GCCounts::new(1, 1, 0, (0, 0))?;
        let mut out = WindowState::new(0, Interval::new(0, 1)?, true, &template)?;
        let mut crossing_parts = Vec::new();

        let seq = vec![b'A'; 99];
        let prefixes = build_gc_prefixes(&seq);

        let mut buf = WindowState::new(1, Interval::new(100, 200)?, false, &template)?;
        buf.counts.set(1, 0, 2.0);
        buf.has_counts = true;

        // Act
        finalize_window_buffer(
            &mut buf,
            &prefixes,
            Interval::new(0, 99)?,
            Interval::new(0, 95)?,
            false,
            true,
            &cfg,
            100.0,
            &mut out,
            &mut crossing_parts,
            0,
        )?;

        // Assert
        assert_eq!(out.weight, 0);
        assert_eq!(crossing_parts.len(), 1);
        assert_eq!(crossing_parts[0].idx, 1);
        assert_eq!(crossing_parts[0].counts.get(1, 0).unwrap(), 2.0);
        assert_eq!(crossing_parts[0].counts.num_acgt_out_of, (0, 0));
        assert_eq!(buf.counts.sum(), 0.0);
        assert!(!buf.has_counts);

        Ok(())
    }
}
mod tests_streaming_parts {

    use anyhow::Result;
    use tempfile::tempdir;

    use cfdnalab::commands::{
        cli_common::{ChromosomeArgs, IOCArgs},
        gc_bias::{
            config::GCConfig,
            counting::GCCounts,
            cross_tile_parts::{CrossingPart, stream_crossing_files, write_crossing_parts},
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
            tmp.path().join("ref_gc_package.npz"),
            ChromosomeArgs::default(),
        );
        cfg.set_min_window_acgt_pct(0);
        cfg
    }

    fn make_template() -> GCCounts {
        GCCounts::new(1, 1, 0, (0, 0)).expect("failed to create template")
    }

    fn counts_from(values: [f64; 2], acgt: (u64, u64)) -> GCCounts {
        GCCounts::from_parts(values.to_vec(), 1, 1, 0, acgt).expect("failed to build counts")
    }

    #[test]
    fn finalizes_windows_missing_from_next_tile() -> Result<()> {
        // Human verification status: unverified
        let tmp = tempdir()?;
        let cfg = make_config(&tmp);
        let template = make_template();
        let temp_dir = tmp.path().to_path_buf();

        // Tile 0 contributes two windows; idx0 vanishes in tile1 so it should flush immediately
        let file_a = write_crossing_parts(
            &temp_dir,
            0,
            &template,
            &[
                CrossingPart {
                    idx: 0,
                    counts: counts_from([2.0, 0.0], (4, 4)),
                },
                CrossingPart {
                    idx: 1,
                    counts: counts_from([1.0, 0.0], (2, 2)),
                },
            ],
        )?
        .expect("first file missing");

        // Tile 1 only continues idx1, signaling idx0 is complete
        let file_b = write_crossing_parts(
            &temp_dir,
            1,
            &template,
            &[CrossingPart {
                idx: 1,
                counts: counts_from([3.0, 0.0], (4, 4)),
            }],
        )?
        .expect("second file missing");

        let (sum, weight) = stream_crossing_files(vec![file_a, file_b], &template, &cfg, 4.0)?;

        assert_eq!(weight, 2, "expected two finalized windows");
        let gc0 = sum.get(1, 0).unwrap();
        // Scaling derivation:
        // - idx0 only appears in tile0: counts [2,0], num_acgt=4 -> mean=1, scale=(1/1)*(4/4)=1 so contribution=2.0
        // - idx1 merges tile0 [1,0] (acgt=2) + tile1 [3,0] (acgt=4) -> counts [4,0], num_acgt=6
        //   mean=2, scale=(1/2)*(6/4)=0.75 so contribution=4*0.75=3.0
        // Total scaled contribution: 2.0 + 3.0 = 5.0
        assert!(
            (gc0 - 5.0).abs() < 1e-9,
            "expected combined scaled count of 5.0, got {gc0}"
        );
        Ok(())
    }

    #[test]
    fn merges_parts_seen_in_every_file() -> Result<()> {
        // Human verification status: unverified
        let tmp = tempdir()?;
        let cfg = make_config(&tmp);
        let template = make_template();
        let temp_dir = tmp.path().to_path_buf();

        // Both tiles carry idx0 so parts must be merged before scaling
        let file_a = write_crossing_parts(
            &temp_dir,
            0,
            &template,
            &[CrossingPart {
                idx: 0,
                counts: counts_from([2.0, 0.0], (4, 4)),
            }],
        )?
        .expect("first file missing");

        let file_b = write_crossing_parts(
            &temp_dir,
            1,
            &template,
            &[CrossingPart {
                idx: 0,
                counts: counts_from([1.0, 1.0], (4, 4)),
            }],
        )?
        .expect("second file missing");

        let (sum, weight) = stream_crossing_files(vec![file_a, file_b], &template, &cfg, 4.0)?;

        assert_eq!(weight, 1, "expected one finalized window");
        let gc0 = sum.get(1, 0).unwrap();
        let gc1 = sum.get(1, 1).unwrap();
        // Counts add to [3,1] and scale factor is 1.0, so scaled counts remain [3,1]
        assert!(
            (gc0 - 3.0).abs() < 1e-9 && (gc1 - 1.0).abs() < 1e-9,
            "expected merged scaled counts of [3.0, 1.0], got [{gc0}, {gc1}]"
        );
        Ok(())
    }

    #[test]
    fn scales_window_when_mean_exceeds_one() -> Result<()> {
        // Human verification status: unverified
        let tmp = tempdir()?;
        let cfg = make_config(&tmp);
        let template = make_template();
        let temp_dir = tmp.path().to_path_buf();

        // Single window across one tile with counts summing to 4 and mean=2
        let file = write_crossing_parts(
            &temp_dir,
            0,
            &template,
            &[CrossingPart {
                idx: 0,
                counts: counts_from([2.0, 2.0], (4, 4)),
            }],
        )?
        .expect("missing crossing file");

        let (sum, weight) = stream_crossing_files(vec![file], &template, &cfg, 4.0)?;

        assert_eq!(weight, 1, "expected one finalized window");
        let gc0 = sum.get(1, 0).unwrap();
        let gc1 = sum.get(1, 1).unwrap();
        // Scaling derivation: counts [2,2], num_acgt=4, mean=2, avg_window_span=4
        // scale = (1/mean)*(num_acgt/avg_window_span) = (1/2)*(4/4) = 0.5
        // scaled counts should be [1,1]
        assert!(
            (gc0 - 1.0).abs() < 1e-9 && (gc1 - 1.0).abs() < 1e-9,
            "expected scaled counts of [1.0, 1.0], got [{gc0}, {gc1}]"
        );
        Ok(())
    }

    #[test]
    fn merges_zero_support_counted_part_with_later_owned_support() -> Result<()> {
        // Human verification status: unverified
        let tmp = tempdir()?;
        let cfg = make_config(&tmp);
        let template = make_template();
        let temp_dir = tmp.path().to_path_buf();

        // Tile 0 is the synthetic "next window" case: counts are present, but this tile owns no
        // support for the window, so it spills with num_acgt_out_of = (0,0).
        let file_a = write_crossing_parts(
            &temp_dir,
            0,
            &template,
            &[CrossingPart {
                idx: 7,
                counts: counts_from([2.0, 0.0], (0, 0)),
            }],
        )?
        .expect("first file missing");

        // Tile 1 owns the support span for the same window and contributes it through the crossing
        // reducer path. After merging:
        // - counts = [2,0]
        // - num_acgt = 100
        // - mean count = (2+0)/2 = 1
        // - scale = (1/1) * (100/100) = 1
        // - final scaled gc0 contribution = 2
        let file_b = write_crossing_parts(
            &temp_dir,
            1,
            &template,
            &[CrossingPart {
                idx: 7,
                counts: counts_from([0.0, 0.0], (100, 100)),
            }],
        )?
        .expect("second file missing");

        let (sum, weight) = stream_crossing_files(vec![file_a, file_b], &template, &cfg, 100.0)?;

        assert_eq!(weight, 1, "merged crossing window should survive reduction");
        let gc0 = sum.get(1, 0).unwrap();
        assert!(
            (gc0 - 2.0).abs() < 1e-9,
            "expected zero-support counts to reunite with later support and scale to 2.0, got {gc0}"
        );
        Ok(())
    }
}
