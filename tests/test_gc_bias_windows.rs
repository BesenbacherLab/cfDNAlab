use cfdnalab::commands::cli_common::WindowSpec;
use cfdnalab::commands::gc_bias::counting::{GCCounts, build_gc_prefixes};
use cfdnalab::commands::gc_bias::windows::{
    WindowState, compute_window_acgt, compute_window_stats,
};
use cfdnalab::shared::bam::Contigs;
use cfdnalab::shared::bed::Windows;
use fxhash::FxHashMap;

fn build_contigs(entries: &[(&str, u32)]) -> Contigs {
    let mut contigs = FxHashMap::default();
    for (name, length) in entries {
        contigs.insert((*name).to_string(), (0, *length));
    }
    Contigs { contigs }
}

fn build_windows(entries: Vec<(u64, u64, u64)>) -> Windows {
    Windows::from_sorted(entries)
}

fn make_window_state(start: u64, end: u64) -> WindowState {
    let template = GCCounts::new(1, 1, 0, (0, 0)).expect("failed to create template");
    WindowState::new(0, start, end, true, &template).expect("failed to build window state")
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

        let stats = compute_window_stats(&window_spec, Some(&map), &contigs, &chromosomes).unwrap();

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

mod tests_compute_window_acgt {
    use super::*;

    #[test]
    fn counts_acgt_when_window_overlaps_sequence() {
        // Sequence AACGTN has four ACGT bases inside window 1-5 (ACGT)
        let seq = b"AACGTN";
        let prefixes = build_gc_prefixes(seq);
        let mut window = make_window_state(1, 5); // covers "ACGT"

        compute_window_acgt(&mut window, &prefixes, 0, seq.len() as u64).unwrap();

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

        compute_window_acgt(&mut window, &prefixes, 0, seq.len() as u64).unwrap();

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

        let err = compute_window_acgt(&mut window, &prefixes, 0, seq.len() as u64)
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

        let err = compute_window_acgt(&mut window, &prefixes, 0, 6)
            .expect_err("expected a prefix bounds error");

        let msg = format!("{err}");
        assert!(
            msg.contains("exceeds prefix length"),
            "unexpected error message: {msg}"
        );
    }
}

mod tests_prepare_tile_windows {
    use anyhow::Result;
    use cfdnalab::{
        commands::{
            cli_common::WindowSpec,
            gc_bias::{counting::GCCounts, windows::prepare_tile_windows},
        },
        shared::tiled_run::{Tile, TileWindowSpan},
    };
    use std::path::PathBuf;

    fn make_template() -> GCCounts {
        GCCounts::new(1, 1, 0, (0, 0)).expect("failed to build template")
    }

    fn make_tile() -> Tile {
        Tile {
            chr: "chr1".to_string(),
            tid: 0,
            index: 0,
            core_start: 100,
            core_end: 190,
            fetch_start: 80,
            fetch_end: 210,
        }
    }

    #[test]
    fn builds_bed_windows_for_tile_core() -> Result<()> {
        let template = make_template();
        let tile = make_tile();
        // Span covers three windows, and the last ends after the core and must be filtered out
        let windows: Vec<(u64, u64, u64)> = vec![(90, 140, 0), (120, 180, 1), (200, 240, 2)];
        let span = TileWindowSpan {
            first_idx: 0,
            last_idx_exclusive: windows.len(),
        };

        let prepared = prepare_tile_windows(
            &WindowSpec::Bed(PathBuf::from("dummy.bed")),
            Some(&windows),
            &tile,
            Some(&span),
            500,
            &template,
        )?;

        assert!(!prepared.skip_tile);
        assert!(prepared.streaming_buffers.is_none());
        assert_eq!(prepared.windows.len(), 2);
        assert_eq!(prepared.windows[0].idx, 0);
        assert_eq!(prepared.windows[1].idx, 1);
        assert!(!prepared.windows[0].contained);
        assert!(prepared.windows[1].contained);
        Ok(())
    }

    #[test]
    fn skips_tile_when_no_bed_windows_available() -> Result<()> {
        let template = make_template();
        let tile = make_tile();
        // Empty BED slice should return skip=true so caller can bail out early
        let windows: Vec<(u64, u64, u64)> = Vec::new();

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
        let tile = Tile {
            core_start: 250,
            core_end: 450,
            ..make_tile()
        };

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
        assert_eq!(current.start, 200);
        assert_eq!(current.end, 400);
        assert!(!current.contained);

        assert_eq!(next.idx, 2);
        assert_eq!(next.start, 400);
        assert_eq!(next.end, 600);
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
        assert_eq!(window.start, tile.core_start as u64);
        assert_eq!(window.end, tile.core_end as u64);
        assert!(window.contained);
        Ok(())
    }
}

mod tests_gc_bias_window_logic {
    use anyhow::Result;
    use ndarray::{Array1, Array2};
    use tempfile::tempdir;

    use cfdnalab::commands::{
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
            tmp.path().join("ref_gc"),
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

        let scaled = process_window(counts, &cfg, Some(100.0))?.expect("window should be retained");

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
        let result = process_window(counts, &cfg, Some(100.0))?;

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

        let (merged, weight) = stream_crossing_files(vec![file1, file2], &template, &cfg, 20.0)?;

        // Assert
        assert_eq!(weight, 1, "one window should contribute once");
        let v = merged.get(10, 0).unwrap();
        assert!((v - 27.5).abs() < 1e-6);
        Ok(())
    }
}
