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
            tmp.path().join("ref_gc"),
            ChromosomeArgs::default(),
        );
        cfg.set_min_window_acgt_pct(0);
        cfg
    }

    #[test]
    fn window_bounds_caps_at_chrom_len() {
        let interval = fixed_size_window_interval(8, 100, 850).expect("interval should be valid");
        assert_eq!(interval.start(), 800);
        assert_eq!(interval.end(), 850);
    }

    #[test]
    fn overlap_length_returns_expected_span() {
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
            0,
            2,
            0,
            2,
            true,
            true,
            &cfg,
            2.0,
            &mut out,
            &mut crossing_parts,
        )?;

        assert_eq!(out.weight, 1);
        let val = out.counts.get(1, 0).unwrap();
        assert!((val - 2.0).abs() < 1e-6);
        assert!(crossing_parts.is_empty());
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
            tmp.path().join("ref_gc"),
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
}
