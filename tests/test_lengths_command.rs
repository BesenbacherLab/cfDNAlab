#![cfg(feature = "cmd_lengths")]

mod fixtures;

mod tests_lengths_command {

    use super::*;

    use anyhow::Result;
    use cfdnalab::commands::cli_common::{
        AssignToWindowArgs, ChromosomeArgs, IOCArgs, WindowAssigner, WindowsArgs,
    };
    use cfdnalab::commands::lengths::config::LengthsConfig;
    use cfdnalab::commands::lengths::lengths::run;
    use cfdnalab::shared::indel_mode::IndelMode;
    use fixtures::simple_inward_bam;
    use ndarray::Array2;
    use ndarray_npy::read_npy;
    use tempfile::TempDir;

    fn base_chromosomes(chrs: &[&str]) -> ChromosomeArgs {
        ChromosomeArgs {
            chromosomes: Some(chrs.iter().map(|c| c.to_string()).collect()),
            chromosomes_file: None,
        }
    }

    #[test]
    fn counts_reference_lengths_global_window() -> Result<()> {
        let bam = simple_inward_bam()?;
        let out_dir = TempDir::new()?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 2,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(WindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 200;
        }

        run(&cfg)?;

        let prefix = cfg.output_prefix.trim();
        let npy_path = out_dir.path().join(format!("{prefix}.length_counts.npy"));
        assert!(npy_path.exists());
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.shape(), &[1, 191]);
        let len60_idx = 60 - 10; // min_fragment_length
        assert!((arr[(0, len60_idx)] - 1.0).abs() < 1e-6);
        assert_eq!(arr[(0, len60_idx - 1)], 0.0);

        Ok(())
    }

    #[test]
    fn counts_reference_lengths_size_single_window_misaligned_tiles() -> Result<()> {
        let bam = simple_inward_bam()?;
        let out_dir = TempDir::new()?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 2,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(WindowsArgs {
            by_size: Some(500),
            by_bed: None,
        });
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_tile_size(50);
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 200;
        }

        run(&cfg)?;

        let prefix = cfg.output_prefix.trim();
        let npy_path = out_dir.path().join(format!("{prefix}.length_counts.npy"));
        assert!(npy_path.exists());
        let arr: Array2<f64> = read_npy(&npy_path)?;
        // Chromosome length 200, window size 500 -> one window
        assert_eq!(arr.shape(), &[1, 191]);
        let len60_idx = 60 - 10;
        assert!((arr[(0, len60_idx)] - 1.0).abs() < 1e-6);
        assert!((arr.sum() - 1.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn counts_reference_lengths_bed_single_window() -> Result<()> {
        let bam = simple_inward_bam()?;
        let out_dir = TempDir::new()?;

        let bed_path = out_dir.path().join("windows.bed");
        fixtures::write_bed(&bed_path, &[("chr1", 0, 200, "w0")])?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 2,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(WindowsArgs {
            by_size: None,
            by_bed: Some(bed_path.clone()),
        });
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 200;
        }

        run(&cfg)?;

        let prefix = cfg.output_prefix.trim();
        let npy_path = out_dir.path().join(format!("{prefix}.length_counts.npy"));
        assert!(npy_path.exists());
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.shape(), &[1, 191]);
        let len60_idx = 60 - 10;
        assert!((arr[(0, len60_idx)] - 1.0).abs() < 1e-6);
        assert!((arr.sum() - 1.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn counts_apply_scaling_factors() -> Result<()> {
        let bam = simple_inward_bam()?;
        let out_dir = TempDir::new()?;

        let scaling_path = out_dir.path().join("scaling.tsv");
        fixtures::write_scaling_factors(&scaling_path, &[("chr1", 0, 200, 2.0)])?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 2,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(WindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_scaling_factors(Some(scaling_path.clone()));
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 200;
        }

        run(&cfg)?;

        let prefix = cfg.output_prefix.trim();
        let npy_path = out_dir.path().join(format!("{prefix}.length_counts.npy"));
        assert!(npy_path.exists());
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.shape(), &[1, 191]);
        let len60_idx = 60 - 10;
        assert!((arr[(0, len60_idx)] - 2.0).abs() < 1e-6);
        assert!((arr.sum() - 2.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn counts_are_zero_when_blacklisted() -> Result<()> {
        let bam = simple_inward_bam()?;
        let out_dir = TempDir::new()?;

        let blacklist_path = out_dir.path().join("blacklist.bed");
        fixtures::write_bed(&blacklist_path, &[("chr1", 0, 200, "bl")])?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 2,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(WindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.blacklist = Some(vec![blacklist_path.clone()]);
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 200;
        }

        run(&cfg)?;

        let prefix = cfg.output_prefix.trim();
        let npy_path = out_dir.path().join(format!("{prefix}.length_counts.npy"));
        assert!(npy_path.exists());
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.shape(), &[1, 191]);
        assert!((arr.sum()).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn assignment_modes_produce_distinct_counts() -> Result<()> {
        let bam = simple_inward_bam()?;
        let window_bp = 40u64;
        let len_idx = 60 - 10;

        let run_with_mode = |assign_by: WindowAssigner| -> Result<Array2<f64>> {
            let out_dir = TempDir::new()?;
            let mut cfg = LengthsConfig::new(
                IOCArgs {
                    bam: bam.bam.clone(),
                    output_dir: out_dir.path().to_path_buf(),
                    n_threads: 2,
                },
                base_chromosomes(&["chr1"]),
            );
            cfg.set_indel_mode(IndelMode::Ignore);
            cfg.set_windows(WindowsArgs {
                by_size: Some(window_bp),
                by_bed: None,
            });
            cfg.set_window_assignment(AssignToWindowArgs { assign_by });
            cfg.set_min_mapq(0);
            cfg.set_require_proper_pair(false);
            {
                let frag = cfg.fragment_lengths_mut();
                frag.min_fragment_length = 10;
                frag.max_fragment_length = 200;
            }

            run(&cfg)?;

            let prefix = cfg.output_prefix.trim();
            let npy_path = out_dir.path().join(format!("{prefix}.length_counts.npy"));
            let arr: Array2<f64> = read_npy(&npy_path)?;
            Ok(arr)
        };

        let arr_all = run_with_mode(WindowAssigner::All)?;
        let arr_mid = run_with_mode(WindowAssigner::Midpoint)?;
        let arr_prop = run_with_mode(WindowAssigner::Proportion(0.3))?;

        // Single fragment spans [20, 80) so it touches window 0 (0-40) and window 1 (40-80)
        // ALL requires full containment in one window, so nothing counts
        assert_eq!(arr_all.shape(), &[5, 191]);
        assert_eq!(arr_mid.shape(), &[5, 191]);
        assert_eq!(arr_prop.shape(), &[5, 191]);

        // ALL drops the fragment because it crosses the window boundary
        assert!((arr_all.sum()).abs() < 1e-6);
        // MIDPOINT picks a random center at 49 or 50, both live in window 1
        assert!((arr_mid[(1, len_idx)] - 1.0).abs() < 1e-6);
        assert!((arr_mid.sum() - 1.0).abs() < 1e-6);
        // PROPORTION=0.3 counts windows with at least 30% overlap:
        // window 0 overlap is 20/60 ≈ 0.33, window 1 overlap is 40/60 ≈ 0.67 
        // (overlap is fragment positions inside the window, not window-bases covered)
        assert!((arr_prop[(0, len_idx)] - 1.0).abs() < 1e-6);
        assert!((arr_prop[(1, len_idx)] - 1.0).abs() < 1e-6);
        assert!((arr_prop.sum() - 2.0).abs() < 1e-6);
        Ok(())
    }
}

mod tests_lengths_tiling_reducer {

    #![cfg(feature = "cmd_lengths")]

    use anyhow::Result;
    use cfdnalab::commands::lengths::counting::LengthCounts;
    use cfdnalab::commands::lengths::tiling::{
        reduce_partials_for_chr, write_cross_npy, write_partials_npz,
    };
    use ndarray::{Array1, Array2, ShapeBuilder};
    use ndarray_npy::NpzWriter;
    use std::fs::File;
    use tempfile::TempDir;

    fn template_counts() -> LengthCounts {
        let lc = LengthCounts::new(10, 10);
        lc
    }

    fn counts_with_value(val: f64) -> LengthCounts {
        let mut lc = template_counts();
        lc.counts[0] = val;
        lc
    }

    #[test]
    fn reducer_accepts_contained_only() -> Result<()> {
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let template = template_counts();

        let counts = vec![counts_with_value(2.0)];
        let contained = vec![true];
        write_partials_npz(dir, "partials", "chr1", 0, &[0], &contained, &counts)?;
        // No cross file because window is contained

        let reduced = reduce_partials_for_chr("chr1", dir, "partials", "cross", 1, &template)?;
        assert_eq!(reduced.len(), 1);
        assert!((reduced[0].counts[0] - 2.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn reducer_counts_multiple_crossing_tiles() -> Result<()> {
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let template = template_counts();

        let counts = vec![counts_with_value(1.0)];
        let contained = vec![false];
        // Two tiles, both crossing the same window
        write_partials_npz(dir, "partials", "chr1", 0, &[0], &contained, &counts)?;
        write_partials_npz(dir, "partials", "chr1", 1, &[0], &contained, &counts)?;
        write_cross_npy(dir, "cross", "chr1", 0, &[0])?;
        write_cross_npy(dir, "cross", "chr1", 1, &[0])?;

        let reduced = reduce_partials_for_chr("chr1", dir, "partials", "cross", 1, &template)?;
        assert_eq!(reduced.len(), 1);
        assert!((reduced[0].counts[0] - 2.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn reducer_combines_contained_and_cross() -> Result<()> {
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let template = template_counts();

        let contained_counts = vec![counts_with_value(1.0)];
        let crossing_counts = vec![counts_with_value(3.0)];
        write_partials_npz(dir, "partials", "chr1", 0, &[0], &[true], &contained_counts)?;
        write_partials_npz(dir, "partials", "chr1", 1, &[0], &[false], &crossing_counts)?;
        write_cross_npy(dir, "cross", "chr1", 1, &[0])?;

        let reduced = reduce_partials_for_chr("chr1", dir, "partials", "cross", 1, &template)?;
        assert_eq!(reduced.len(), 1);
        // Expect 1 contained contribution and 1 crossing contribution => sum counts
        assert!((reduced[0].counts[0] - 4.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn reducer_errors_when_contribution_missing() {
        let tmp = TempDir::new().unwrap();
        let dir = tmp.path();
        let template = template_counts();

        // No partials written -> zero contributions
        let err = reduce_partials_for_chr("chr1", dir, "partials", "cross", 1, &template)
            .expect_err("should fail when contributions are missing");
        assert!(err.to_string().contains("expected 1"));
    }

    #[test]
    fn reducer_errors_on_mismatched_counts() {
        let tmp = TempDir::new().unwrap();
        let dir = tmp.path();
        let template = template_counts();

        // Cross file claims one contribution, but no partial exists
        write_cross_npy(dir, "cross", "chr1", 0, &[0]).unwrap();

        let err = reduce_partials_for_chr("chr1", dir, "partials", "cross", 1, &template)
            .expect_err("should fail when expected contributions not met");
        assert!(err.to_string().contains("expected 1"));
    }

    #[test]
    fn reducer_errors_on_counts_width_mismatch() {
        let tmp = TempDir::new().unwrap();
        let dir = tmp.path();
        let template = template_counts(); // width 1

        let idxs = Array1::from(vec![0u64]);
        let contained = Array1::from(vec![1u8]);
        let counts = Array2::from_shape_vec((1, 2), vec![1.0, 0.5]).unwrap();
        let path = dir.join("partials.chr1.0.npz");
        let file = File::create(&path).unwrap();
        let mut npz = NpzWriter::new(file);
        npz.add_array("window_idx_chr", &idxs).unwrap();
        npz.add_array("contained", &contained).unwrap();
        npz.add_array("counts", &counts).unwrap();
        npz.finish().unwrap();

        let err = reduce_partials_for_chr("chr1", dir, "partials", "cross", 1, &template)
            .expect_err("should fail on counts width mismatch");
        assert!(err.to_string().contains("counts width mismatch"));
    }

    #[test]
    fn reducer_errors_on_non_contiguous_counts_rows() {
        let tmp = TempDir::new().unwrap();
        let dir = tmp.path();
        let template = LengthCounts::new(10, 11); // two-length template

        // Two rows, both targeting window 0, but stored in Fortran order so each
        // row slice is non-contiguous and should be rejected.
        let idxs = Array1::from(vec![0u64, 0u64]);
        let contained = Array1::from(vec![1u8, 1u8]);
        let counts = Array2::from_shape_vec((2, 2).f(), vec![1.0, 0.5, 2.0, 1.5]).unwrap();
        let path = dir.join("partials.chr1.0.npz");
        let file = File::create(&path).unwrap();
        let mut npz = NpzWriter::new(file);
        npz.add_array("window_idx_chr", &idxs).unwrap();
        npz.add_array("contained", &contained).unwrap();
        npz.add_array("counts", &counts).unwrap();
        npz.finish().unwrap();

        let err = reduce_partials_for_chr("chr1", dir, "partials", "cross", 1, &template)
            .expect_err("should fail on non-contiguous counts rows");
        assert!(err.to_string().contains("counts row not contiguous"));
    }

    #[test]
    fn reducer_ignores_files_from_other_chromosomes() -> Result<()> {
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let template = template_counts();

        let counts = vec![counts_with_value(1.0)];
        let contained = vec![true];
        write_partials_npz(dir, "partials", "chr1", 0, &[0], &contained, &counts)?;

        // Stray files for chr2 that should be filtered out by the chr-specific prefix match
        write_partials_npz(dir, "partials", "chr2", 0, &[0], &contained, &counts)?;
        write_cross_npy(dir, "cross", "chr2", 0, &[0])?;

        let reduced = reduce_partials_for_chr("chr1", dir, "partials", "cross", 1, &template)?;
        assert_eq!(reduced.len(), 1);
        assert!((reduced[0].counts[0] - 1.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn write_partials_rejects_mismatched_contained() {
        let tmp = TempDir::new().unwrap();
        let dir = tmp.path();
        let template = template_counts();
        let counts = vec![template];
        let contained = vec![true, false]; // Wrong length

        let err = write_partials_npz(dir, "partials", "chr1", 0, &[0], &contained, &counts)
            .expect_err("should error on contained/idx length mismatch");
        assert!(err.to_string().contains("contained flags length mismatch"));
    }

    #[test]
    fn reducer_errors_on_out_of_bounds_partial_idx() {
        let tmp = TempDir::new().unwrap();
        let dir = tmp.path();
        let template = template_counts();

        let counts = vec![counts_with_value(1.0)];
        let contained = vec![false];
        // Write a partial with idx outside n_windows=1
        write_partials_npz(dir, "partials", "chr1", 0, &[2], &contained, &counts).unwrap();

        let err = reduce_partials_for_chr("chr1", dir, "partials", "cross", 1, &template)
            .expect_err("should fail on out-of-bounds idx");
        assert!(err.to_string().contains("out of bounds"));
    }

    #[test]
    fn reducer_errors_on_out_of_bounds_cross_idx() {
        let tmp = TempDir::new().unwrap();
        let dir = tmp.path();
        let template = template_counts();

        write_cross_npy(dir, "cross", "chr1", 0, &[3]).unwrap();
        let err = reduce_partials_for_chr("chr1", dir, "partials", "cross", 1, &template)
            .expect_err("should fail on cross index out of bounds");
        assert!(err.to_string().contains("Cross index"));
    }

    #[test]
    fn reducer_separates_windows() -> Result<()> {
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let template = template_counts();

        let counts0 = vec![counts_with_value(1.0)];
        let counts1 = vec![counts_with_value(2.0)];
        let contained = vec![true];

        // Window 0 contained in tile 0, window 1 contained in tile 1
        write_partials_npz(dir, "partials", "chr1", 0, &[0], &contained, &counts0)?;
        write_partials_npz(dir, "partials", "chr1", 1, &[1], &contained, &counts1)?;

        let reduced = reduce_partials_for_chr("chr1", dir, "partials", "cross", 2, &template)?;
        assert_eq!(reduced.len(), 2);
        assert!((reduced[0].counts[0] - 1.0).abs() < 1e-6);
        assert!((reduced[1].counts[0] - 2.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn write_partials_skips_empty() -> Result<()> {
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let template = template_counts();
        let res = write_partials_npz(dir, "partials", "chr1", 0, &[], &[], &[])?;
        assert!(res.is_none());
        // Ensure reducer still errors because nothing was written
        let err = reduce_partials_for_chr("chr1", dir, "partials", "cross", 1, &template)
            .expect_err("should fail when nothing written");
        assert!(err.to_string().contains("expected 1"));
        Ok(())
    }

    #[test]
    fn write_cross_skips_empty() -> Result<()> {
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let res = write_cross_npy(dir, "cross", "chr1", 0, &[])?;
        assert!(res.is_none());
        Ok(())
    }
}

mod tests_lengths_tiling_helpers {

    use cfdnalab::commands::cli_common::WindowSpec;
    use cfdnalab::commands::lengths::tiling::fetch_span_for_tile;
    use cfdnalab::shared::bam::Contigs;
    use cfdnalab::shared::tiled_run::{Tile, TileWindowSpan, build_tiles};
    use fxhash::FxHashMap;
    use std::path::PathBuf;

    #[test]
    fn fetch_span_size_mode_clamps_to_halo_and_chrom() {
        // Tile: core 50-150, fetch 30-200 (halo 20 left, 50 right), chrom len 180
        let tile = Tile {
            chr: "chr1".to_string(),
            tid: 0,
            index: 0,
            core_start: 50,
            core_end: 150,
            fetch_start: 30,
            fetch_end: 200,
        };
        let span = fetch_span_for_tile(&tile, None, None, &WindowSpec::Size(100), 180)
            .expect("span expected");
        // Window span touching core: 0..200, after halo clamp -> 30..180
        assert_eq!(span, (30, 180));
    }

    #[test]
    fn build_tiles_aligns_to_bin_when_divisible() {
        let mut contigs = FxHashMap::default();
        contigs.insert("chr1".to_string(), (0, 100u32));
        let contigs = Contigs { contigs };
        let (tiles, aligned) =
            build_tiles(&vec!["chr1".to_string()], &contigs, 30, 0, Some(10)).unwrap();
        assert!(aligned);
        // Cores should start on multiples of 10
        for t in &tiles {
            assert_eq!((t.core_start as u64) % 10, 0);
        }
        // Expect four tiles: 0-30,30-60,60-90,90-100
        assert_eq!(tiles.len(), 4);
        assert_eq!(tiles[0].core_end, 30);
        assert_eq!(tiles[3].core_start, 90);
        assert_eq!(tiles[3].core_end, 100);
    }

    #[test]
    fn build_tiles_not_aligned_when_too_few_bins() {
        let mut contigs = FxHashMap::default();
        contigs.insert("chr1".to_string(), (0, 50u32));
        let contigs = Contigs { contigs };
        // With tile_bp=15 and align_bp=10, only one full 10bp bin fits,
        // and build_tiles requires at least 10 bins (k >= 10) before rounding down.
        // So alignment should be disabled and tiles keep the original 15bp size.
        let (_tiles, aligned) =
            build_tiles(&vec!["chr1".to_string()], &contigs, 15, 0, Some(10)).unwrap();
        assert!(!aligned);
    }

    #[test]
    fn fetch_span_for_tile_global_clamps_to_chrom() {
        let tile = Tile {
            chr: "chr1".to_string(),
            tid: 0,
            index: 0,
            core_start: 0,
            core_end: 50,
            fetch_start: 0,
            fetch_end: 200,
        };
        let span = fetch_span_for_tile(&tile, None, None, &WindowSpec::Global, 120).expect("span");
        assert_eq!(span, (0, 120));
    }

    #[test]
    fn fetch_span_for_tile_bed_with_overlap() {
        let tile = Tile {
            chr: "chr1".to_string(),
            tid: 0,
            index: 0,
            core_start: 100,
            core_end: 160,
            fetch_start: 80,
            fetch_end: 200,
        };
        let windows: Vec<(u64, u64, u64)> = vec![(90, 110, 0), (150, 170, 1), (250, 300, 2)];
        let span = TileWindowSpan {
            first_idx: 0,
            last_idx_exclusive: 2,
        };
        let res = fetch_span_for_tile(
            &tile,
            Some(&span),
            Some(&windows),
            &WindowSpec::Bed(PathBuf::from("dummy")),
            500,
        )
        .expect("span");
        // min_ws=90, max_we=170, halos: left 20, right 40 -> widened to 70..210, clamped to fetch
        assert_eq!(res, (80, 200));
    }

    #[test]
    fn fetch_span_bed_none_when_no_overlap() {
        let tile = Tile {
            chr: "chr1".to_string(),
            tid: 0,
            index: 0,
            core_start: 100,
            core_end: 150,
            fetch_start: 80,
            fetch_end: 170,
        };
        // No windows overlap tile
        let windows: [(u64, u64, u64); 0] = [];
        let span = TileWindowSpan {
            first_idx: 0,
            last_idx_exclusive: 0,
        };
        let res = fetch_span_for_tile(
            &tile,
            Some(&span),
            Some(&windows),
            &WindowSpec::Bed(PathBuf::from("dummy")),
            200,
        );
        assert!(res.is_none());
    }

    #[test]
    fn fetch_span_size_mode_none_when_tile_right_of_chromosome() {
        let tile = Tile {
            chr: "chr1".to_string(),
            tid: 0,
            index: 0,
            core_start: 250,
            core_end: 260,
            fetch_start: 230,
            fetch_end: 270,
        };
        let res = fetch_span_for_tile(&tile, None, None, &WindowSpec::Size(50), 200);
        assert!(res.is_none());
    }

    #[test]
    fn fetch_span_size_mode_none_for_empty_core() {
        let tile = Tile {
            chr: "chr1".to_string(),
            tid: 0,
            index: 0,
            core_start: 100,
            core_end: 100,
            fetch_start: 80,
            fetch_end: 120,
        };
        let res = fetch_span_for_tile(&tile, None, None, &WindowSpec::Size(50), 150);
        assert!(res.is_none());
    }
}
