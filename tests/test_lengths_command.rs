#![cfg(feature = "cmd_lengths")]

mod fixtures;

mod tests_lengths_command {

    use super::*;

    use anyhow::Result;
    use cfdnalab::commands::cli_common::{
        AssignToWindowArgs, ChromosomeArgs, IOCArgs, WindowAssigner, WindowsArgs,
    };
    #[cfg(feature = "cmd_coverage_weights")]
    use cfdnalab::commands::coverage_weights::{
        config::CoverageWeightsConfig, coverage_weights::run as run_coverage_weights,
    };
    use cfdnalab::commands::gc_bias::{
        GC_CORRECTION_SCHEMA_VERSION, correct::MarginalizeLengthsWeightingScheme,
        package::GCCorrectionPackage,
    };
    use cfdnalab::commands::lengths::config::LengthsConfig;
    use cfdnalab::commands::lengths::lengths::run;
    use cfdnalab::shared::blacklist::strategy::BlacklistStrategy;
    use cfdnalab::shared::indel_mode::IndelMode;
    use cfdnalab::shared::io::dot_join;
    use fixtures::{
        BamFixture, FragmentSpec, ReadSpec, bam_from_specs, build_real_neutral_gc_package,
        build_real_non_neutral_gc_package, simple_inward_bam, simple_reference_twobit, write_bed,
        write_scaling_factors,
    };
    use ndarray::Array2;
    use ndarray::array;
    use ndarray_npy::read_npy;
    use tempfile::TempDir;

    fn base_chromosomes(chrs: &[&str]) -> ChromosomeArgs {
        ChromosomeArgs {
            chromosomes: Some(chrs.iter().map(|c| c.to_string()).collect()),
            chromosomes_file: None,
        }
    }

    fn fragment_on_tid(tid: usize, start: i64, fragment_len: i64, read_len: i64) -> FragmentSpec {
        let mut fragment = fixtures::paired_fragment(start, fragment_len, read_len);
        fragment.forward.tid = tid;
        fragment.reverse.tid = tid;
        fragment.forward.mate_tid = Some(tid);
        fragment.reverse.mate_tid = Some(tid);
        fragment
    }

    fn single_read_fragment_bam(name: &str, fragment_start: i64, fragment_len: u32) -> Result<BamFixture> {
        bam_from_specs(
            vec![("chr1".to_string(), 200)],
            Vec::new(),
            vec![ReadSpec {
                tid: 0,
                pos: fragment_start,
                cigar: vec![('M', fragment_len)],
                seq: vec![b'A'; fragment_len as usize],
                qual: 40,
                is_reverse: false,
                mapq: 60,
                flags: 0,
                mate_tid: None,
                mate_pos: None,
                insert_size: 0,
            }],
            name,
        )
    }

    fn three_chrom_length_fixture(name: &str) -> Result<BamFixture> {
        bam_from_specs(
            vec![
                ("chr1".to_string(), 200),
                ("chr2".to_string(), 200),
                ("chr3".to_string(), 200),
            ],
            vec![
                fragment_on_tid(0, 20, 60, 20),
                fragment_on_tid(1, 30, 80, 20),
                fragment_on_tid(2, 40, 100, 20),
            ],
            Vec::new(),
            name,
        )
    }

    #[cfg(feature = "cmd_coverage_weights")]
    fn make_simple_coverage_weights_config(
        out_dir: &std::path::Path,
        bam: &std::path::Path,
    ) -> CoverageWeightsConfig {
        let mut cfg = CoverageWeightsConfig::new(
            IOCArgs {
                bam: bam.to_path_buf(),
                output_dir: out_dir.to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_bin_size(20);
        cfg.set_stride(20);
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_output_prefix("coverage".to_string());
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 200;
        }
        cfg
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
        let npy_path = out_dir
            .path()
            .join(dot_join(&[prefix, "length_counts.npy"]));
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
        let npy_path = out_dir
            .path()
            .join(dot_join(&[prefix, "length_counts.npy"]));
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
    fn counts_reference_lengths_size_aligned_tiles_reduce_cross_tile_bins() -> Result<()> {
        let bam = bam_from_specs(
            vec![("chr1".to_string(), 200)],
            vec![fixtures::paired_fragment(95, 40, 20)],
            Vec::new(),
            "lengths_size_aligned_cross_tile",
        )?;
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
            by_size: Some(10),
            by_bed: None,
        });
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_tile_size(100);
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 50;
        }

        run(&cfg)?;

        let prefix = cfg.output_prefix.trim();
        let npy_path = out_dir
            .path()
            .join(dot_join(&[prefix, "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;

        // The fragment spans [95, 135), so in count-overlap mode it contributes:
        //  5/40 to 90-100
        // 10/40 to 100-110
        // 10/40 to 110-120
        // 10/40 to 120-130
        //  5/40 to 130-140
        // The three middle bins sit fully inside tile 1, but tile 0 still reaches them.
        let len40_idx = 40 - 10;
        assert_eq!(arr.shape(), &[20, 41]);
        assert!((arr[(9, len40_idx)] - 0.125).abs() < 1e-6);
        assert!((arr[(10, len40_idx)] - 0.25).abs() < 1e-6);
        assert!((arr[(11, len40_idx)] - 0.25).abs() < 1e-6);
        assert!((arr[(12, len40_idx)] - 0.25).abs() < 1e-6);
        assert!((arr[(13, len40_idx)] - 0.125).abs() < 1e-6);
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
        let npy_path = out_dir
            .path()
            .join(dot_join(&[prefix, "length_counts.npy"]));
        assert!(npy_path.exists());
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.shape(), &[1, 191]);
        let len60_idx = 60 - 10;
        assert!((arr[(0, len60_idx)] - 1.0).abs() < 1e-6);
        assert!((arr.sum() - 1.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn global_by_size_and_bed_full_chromosome_windows_match_exactly() -> Result<()> {
        // Arrange:
        // `simple_inward_bam()` contains one fragment spanning [20, 80), length 60, on a single
        // 200 bp chromosome.
        //
        // Compare three logically equivalent window specifications:
        // - global mode
        // - one by-size window [0, 200)
        // - one BED window [0, 200)
        //
        // In all three cases there is only one logical window covering the full chromosome, so the
        // fragment is fully assigned to that window. The expected output matrix is therefore:
        // - shape [1, 191] for lengths 10..=200
        // - one count in length bin 60
        // - zero elsewhere
        //
        // The stronger contract is exact equality of the full arrays, not just equality of the
        // occupied bin.
        let bam = simple_inward_bam()?;
        let global_out = TempDir::new()?;
        let by_size_out = TempDir::new()?;
        let by_bed_out = TempDir::new()?;
        let bed_path = by_bed_out.path().join("whole_chr_window.bed");
        fixtures::write_bed(&bed_path, &[("chr1", 0, 200, "whole_chr")])?;

        let make_cfg = |out_dir: &std::path::Path, windows: WindowsArgs| {
            let mut cfg = LengthsConfig::new(
                IOCArgs {
                    bam: bam.bam.clone(),
                    output_dir: out_dir.to_path_buf(),
                    n_threads: 1,
                },
                base_chromosomes(&["chr1"]),
            );
            cfg.set_indel_mode(IndelMode::Ignore);
            cfg.set_windows(windows);
            cfg.set_window_assignment(AssignToWindowArgs::default());
            cfg.set_min_mapq(0);
            cfg.set_require_proper_pair(false);
            {
                let frag = cfg.fragment_lengths_mut();
                frag.min_fragment_length = 10;
                frag.max_fragment_length = 200;
            }
            cfg
        };

        let global_cfg = make_cfg(global_out.path(), WindowsArgs::default());
        let by_size_cfg = make_cfg(
            by_size_out.path(),
            WindowsArgs {
                by_size: Some(200),
                by_bed: None,
            },
        );
        let by_bed_cfg = make_cfg(
            by_bed_out.path(),
            WindowsArgs {
                by_size: None,
                by_bed: Some(bed_path),
            },
        );

        // Act
        run(&global_cfg)?;
        run(&by_size_cfg)?;
        run(&by_bed_cfg)?;

        // Assert
        let read_counts = |dir: &TempDir| -> Result<Array2<f64>> {
            let npy_path = dir.path().join(dot_join(&["", "length_counts.npy"]));
            read_npy(&npy_path).map_err(Into::into)
        };

        let global_arr = read_counts(&global_out)?;
        let by_size_arr = read_counts(&by_size_out)?;
        let by_bed_arr = read_counts(&by_bed_out)?;

        let len60_idx = 60 - 10;
        assert_eq!(global_arr.shape(), &[1, 191]);
        assert_eq!(global_arr, by_size_arr);
        assert_eq!(global_arr, by_bed_arr);
        assert!((global_arr[(0, len60_idx)] - 1.0).abs() < 1e-12);
        assert!((global_arr.sum() - 1.0).abs() < 1e-12);

        Ok(())
    }

    #[test]
    fn default_min_mapq_matches_explicit_thirty_and_differs_from_explicit_zero() -> Result<()> {
        // Arrange:
        // Use three fragments with distinct lengths and MAPQ:
        // - [20, 80): length 60, MAPQ 60
        // - [100, 170): length 70, MAPQ 0
        // - [180, 260): length 80, MAPQ 30
        //
        // In global mode the output is one row of length counts, so the expected
        // distributions are:
        // - default `min_mapq = 30`: lengths 60 and 80 count once each
        // - explicit `min_mapq = 30`: identical to default
        // - explicit `min_mapq = 0`: lengths 60, 70, and 80 count once each
        let fragment_with_mapq =
            |start: i64, fragment_len: i64, mapq: u8| -> fixtures::FragmentSpec {
                let mut fragment = fixtures::paired_fragment(start, fragment_len, 20);
                fragment.forward.mapq = mapq;
                fragment.reverse.mapq = mapq;
                fragment
            };
        let bam = bam_from_specs(
            vec![("chr1".to_string(), 300)],
            vec![
                fragment_with_mapq(20, 60, 60),
                fragment_with_mapq(100, 70, 0),
                fragment_with_mapq(180, 80, 30),
            ],
            Vec::new(),
            "lengths_default_min_mapq",
        )?;
        let out_default = TempDir::new()?;
        let out_thirty = TempDir::new()?;
        let out_zero = TempDir::new()?;

        let make_cfg = |out_dir: &std::path::Path, prefix: &str| {
            let mut cfg = LengthsConfig::new(
                IOCArgs {
                    bam: bam.bam.clone(),
                    output_dir: out_dir.to_path_buf(),
                    n_threads: 1,
                },
                base_chromosomes(&["chr1"]),
            );
            cfg.output_prefix = prefix.to_string();
            cfg.set_indel_mode(IndelMode::Ignore);
            cfg.set_windows(WindowsArgs::default());
            cfg.set_window_assignment(AssignToWindowArgs::default());
            cfg.set_require_proper_pair(false);
            {
                let frag = cfg.fragment_lengths_mut();
                frag.min_fragment_length = 50;
                frag.max_fragment_length = 120;
            }
            cfg
        };

        let default_cfg = make_cfg(out_default.path(), "default");
        let mut explicit_thirty_cfg = make_cfg(out_thirty.path(), "explicit_thirty");
        explicit_thirty_cfg.set_min_mapq(30);
        let mut explicit_zero_cfg = make_cfg(out_zero.path(), "explicit_zero");
        explicit_zero_cfg.set_min_mapq(0);

        // Act
        run(&default_cfg)?;
        run(&explicit_thirty_cfg)?;
        run(&explicit_zero_cfg)?;

        // Assert
        let read_counts = |dir: &TempDir, prefix: &str| -> Result<Array2<f64>> {
            let npy_path = dir.path().join(dot_join(&[prefix, "length_counts.npy"]));
            read_npy(&npy_path).map_err(Into::into)
        };

        let default_arr = read_counts(&out_default, "default")?;
        let explicit_thirty_arr = read_counts(&out_thirty, "explicit_thirty")?;
        let explicit_zero_arr = read_counts(&out_zero, "explicit_zero")?;

        let len60_idx = 60 - 50;
        let len70_idx = 70 - 50;
        let len80_idx = 80 - 50;

        assert_eq!(default_arr.shape(), &[1, 71]);
        assert_eq!(default_arr, explicit_thirty_arr);
        assert!((default_arr[(0, len60_idx)] - 1.0).abs() < 1e-6);
        assert_eq!(default_arr[(0, len70_idx)], 0.0);
        assert!((default_arr[(0, len80_idx)] - 1.0).abs() < 1e-6);
        assert!((default_arr.sum() - 2.0).abs() < 1e-6);

        assert!((explicit_zero_arr[(0, len60_idx)] - 1.0).abs() < 1e-6);
        assert!((explicit_zero_arr[(0, len70_idx)] - 1.0).abs() < 1e-6);
        assert!((explicit_zero_arr[(0, len80_idx)] - 1.0).abs() < 1e-6);
        assert!((explicit_zero_arr.sum() - 3.0).abs() < 1e-6);

        Ok(())
    }

    #[test]
    fn counts_reference_lengths_global_window_across_three_chromosomes() -> Result<()> {
        let bam = three_chrom_length_fixture("lengths_three_chr_global")?;
        let out_dir = TempDir::new()?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 2,
            },
            base_chromosomes(&["chr1", "chr2", "chr3"]),
        );
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(WindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 120;
        }

        run(&cfg)?;

        let prefix = cfg.output_prefix.trim();
        let npy_path = out_dir
            .path()
            .join(dot_join(&[prefix, "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;

        // Global mode collapses the selected chromosomes into one combined length distribution.
        // This fixture contributes one fragment of lengths 60, 80, and 100 across the three
        // chromosomes, so the single output row should contain all three counts.
        assert_eq!(arr.shape(), &[1, 111]);
        assert!((arr[(0, 60 - 10)] - 1.0).abs() < 1e-6);
        assert!((arr[(0, 80 - 10)] - 1.0).abs() < 1e-6);
        assert!((arr[(0, 100 - 10)] - 1.0).abs() < 1e-6);
        assert!((arr.sum() - 3.0).abs() < 1e-6);

        Ok(())
    }

    #[test]
    fn unpaired_single_read_matches_paired_fragment_length_count_for_same_span() -> Result<()> {
        // Arrange:
        // Compare two representations of the same physical fragment span [20, 80):
        // - paired-end fixture `simple_inward_bam()`
        // - one unpaired read with aligned span [20, 80)
        //
        // In both commands:
        // - paired mode defines fragment span as [forward.pos, reverse.end)
        // - unpaired `reads_are_fragments` mode defines fragment span as [read.pos, read.end)
        //
        // So both inputs represent one fragment of length 60 and must produce the same global
        // length distribution: one count in length bin 60, zero elsewhere.
        let paired_bam = simple_inward_bam()?;
        let unpaired_bam = single_read_fragment_bam("lengths_unpaired_single_read", 20, 60)?;
        let paired_out = TempDir::new()?;
        let unpaired_out = TempDir::new()?;

        let make_cfg = |bam_path: &std::path::Path, out_dir: &std::path::Path, unpaired: bool| {
            let mut cfg = LengthsConfig::new(
                IOCArgs {
                    bam: bam_path.to_path_buf(),
                    output_dir: out_dir.to_path_buf(),
                    n_threads: 1,
                },
                base_chromosomes(&["chr1"]),
            );
            cfg.set_indel_mode(IndelMode::Ignore);
            cfg.set_windows(WindowsArgs::default());
            cfg.set_window_assignment(AssignToWindowArgs::default());
            cfg.set_min_mapq(0);
            cfg.set_require_proper_pair(false);
            cfg.set_unpaired(cfdnalab::commands::cli_common::UnpairedArgs {
                reads_are_fragments: unpaired,
            });
            {
                let frag = cfg.fragment_lengths_mut();
                frag.min_fragment_length = 10;
                frag.max_fragment_length = 200;
            }
            cfg
        };

        let paired_cfg = make_cfg(&paired_bam.bam, paired_out.path(), false);
        let unpaired_cfg = make_cfg(&unpaired_bam.bam, unpaired_out.path(), true);

        // Act
        run(&paired_cfg)?;
        run(&unpaired_cfg)?;

        // Assert
        let paired_arr: Array2<f64> =
            read_npy(paired_out.path().join(dot_join(&["", "length_counts.npy"])))?;
        let unpaired_arr: Array2<f64> =
            read_npy(unpaired_out.path().join(dot_join(&["", "length_counts.npy"])))?;

        let len60_idx = 60 - 10;
        assert_eq!(paired_arr, unpaired_arr);
        assert_eq!(paired_arr.shape(), &[1, 191]);
        assert!((paired_arr[(0, len60_idx)] - 1.0).abs() < 1e-12);
        assert!((paired_arr.sum() - 1.0).abs() < 1e-12);

        Ok(())
    }

    #[test]
    fn counts_reference_lengths_size_single_window_across_three_chromosomes() -> Result<()> {
        let bam = three_chrom_length_fixture("lengths_three_chr_size")?;
        let out_dir = TempDir::new()?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 2,
            },
            base_chromosomes(&["chr1", "chr2", "chr3"]),
        );
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(WindowsArgs {
            by_size: Some(200),
            by_bed: None,
        });
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_tile_size(30);
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 120;
        }

        run(&cfg)?;

        let prefix = cfg.output_prefix.trim();
        let npy_path = out_dir
            .path()
            .join(dot_join(&[prefix, "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;

        assert_eq!(arr.shape(), &[3, 111]);
        assert!((arr[(0, 60 - 10)] - 1.0).abs() < 1e-6);
        assert!((arr[(1, 80 - 10)] - 1.0).abs() < 1e-6);
        assert!((arr[(2, 100 - 10)] - 1.0).abs() < 1e-6);
        assert!((arr.sum() - 3.0).abs() < 1e-6);

        Ok(())
    }

    #[test]
    fn counts_reference_lengths_bed_single_window_across_three_chromosomes() -> Result<()> {
        let bam = three_chrom_length_fixture("lengths_three_chr_bed")?;
        let out_dir = TempDir::new()?;
        let bed_path = out_dir.path().join("windows_three_chr.bed");
        fixtures::write_bed(
            &bed_path,
            &[
                ("chr1", 0, 200, "chr1_window"),
                ("chr2", 0, 200, "chr2_window"),
                ("chr3", 0, 200, "chr3_window"),
            ],
        )?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 2,
            },
            base_chromosomes(&["chr1", "chr2", "chr3"]),
        );
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(WindowsArgs {
            by_size: None,
            by_bed: Some(bed_path),
        });
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 120;
        }

        run(&cfg)?;

        let prefix = cfg.output_prefix.trim();
        let npy_path = out_dir
            .path()
            .join(dot_join(&[prefix, "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;

        assert_eq!(arr.shape(), &[3, 111]);
        assert!((arr[(0, 60 - 10)] - 1.0).abs() < 1e-6);
        assert!((arr[(1, 80 - 10)] - 1.0).abs() < 1e-6);
        assert!((arr[(2, 100 - 10)] - 1.0).abs() < 1e-6);
        assert!((arr.sum() - 3.0).abs() < 1e-6);

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
        let npy_path = out_dir
            .path()
            .join(dot_join(&[prefix, "length_counts.npy"]));
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
        let npy_path = out_dir
            .path()
            .join(dot_join(&[prefix, "length_counts.npy"]));
        assert!(npy_path.exists());
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.shape(), &[1, 191]);
        assert!((arr.sum()).abs() < 1e-6);
        Ok(())
    }

    fn build_gc_package(path: &std::path::Path, end_offset: u64) -> Result<()> {
        // Two length bins: [10,60) and [60,200]; two GC bins: [0,50) and [50,101]
        let correction_matrix = array![[1.0_f64, 1.0_f64], [2.0_f64, 10.0_f64]];
        let package = GCCorrectionPackage {
            version: GC_CORRECTION_SCHEMA_VERSION,
            end_offset,
            length_edges: vec![10, 60, 200],
            gc_edges: vec![0, 50, 101],
            length_bin_frequencies: array![1.0_f64, 3.0_f64],
            correction_matrix,
        };
        package.write_npz(path)?;
        Ok(())
    }

    #[test]
    fn applies_gc_correction_weighting_modes() -> Result<()> {
        let bam = simple_inward_bam()?;
        let ref_twobit = simple_reference_twobit()?;
        let gc_dir = TempDir::new()?;
        let gc_path = gc_dir.path().join("gc_pkg.npz");
        build_gc_package(&gc_path, 0)?;

        let expected = |scheme: MarginalizeLengthsWeightingScheme| -> f64 {
            match scheme {
                MarginalizeLengthsWeightingScheme::Equal => 5.5, // mean of rows at GC bin 50
                MarginalizeLengthsWeightingScheme::Coverage => 7.75, // weighted by [1,3]
                MarginalizeLengthsWeightingScheme::MaxCoverage => 10.0, // most frequent row
            }
        };

        let run_with_scheme = |scheme: MarginalizeLengthsWeightingScheme| -> Result<f64> {
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
            cfg.set_gc_length_weighting(scheme);
            cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
                gc_file: Some(gc_path.clone()),
                drop_invalid_gc: false,
            });
            cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
            {
                let frag = cfg.fragment_lengths_mut();
                frag.min_fragment_length = 10;
                frag.max_fragment_length = 200;
            }

            run(&cfg)?;

            let prefix = cfg.output_prefix.trim();
            let npy_path = out_dir
                .path()
                .join(dot_join(&[prefix, "length_counts.npy"]));
            let arr: Array2<f64> = read_npy(&npy_path)?;
            let len60_idx = 60 - 10;
            Ok(arr[(0, len60_idx)])
        };

        for scheme in [
            MarginalizeLengthsWeightingScheme::Equal,
            MarginalizeLengthsWeightingScheme::Coverage,
            MarginalizeLengthsWeightingScheme::MaxCoverage,
        ] {
            let observed = run_with_scheme(scheme)?;
            let exp = expected(scheme);
            assert!(
                (observed - exp).abs() < 1e-6,
                "scheme {:?}: expected {}, got {}",
                scheme,
                exp,
                observed
            );
        }

        Ok(())
    }

    #[test]
    fn real_ref_gc_bias_then_gc_bias_package_is_neutral_in_single_bin_case_for_lengths() -> Result<()>
    {
        let bam = simple_inward_bam()?;
        let ref_twobit = simple_reference_twobit()?;
        let out_dir = TempDir::new()?;
        let gc_path =
            build_real_neutral_gc_package(&bam.bam, &ref_twobit.path, out_dir.path(), 60)?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(WindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
            gc_file: Some(gc_path),
            drop_invalid_gc: false,
        });
        cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
        cfg.set_gc_length_weighting(MarginalizeLengthsWeightingScheme::Equal);
        {
            let fragment_lengths = cfg.fragment_lengths_mut();
            fragment_lengths.min_fragment_length = 60;
            fragment_lengths.max_fragment_length = 60;
        }

        // Manual expectations:
        // - `ref-gc-bias` is run for exactly one fragment length: 60 bp.
        // - On the simple reference ("ACGT" repeated), every 60 bp fragment has GC%=50,
        //   so the reference package has exactly one populated GC-by-length cell.
        // - `gc-bias` on `simple_inward_bam` also places all cfDNA mass in that same single cell.
        // - A 1x1 normalized cfDNA count divided by a 1x1 normalized reference count gives 1.0,
        //   so the produced correction package is neutral.
        // - `lengths` therefore receives GC weight 1.0 for the only fragment and must match the
        //   uncorrected global result: one count at fragment length 60 and zero elsewhere.
        run(&cfg)?;

        let npy_path = out_dir
            .path()
            .join(dot_join(&[cfg.output_prefix.trim(), "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.dim(), (1, 91));

        let len60_idx = 60 - 10;
        assert!((arr[(0, len60_idx)] - 1.0).abs() < 1e-12);
        for idx in 0..arr.ncols() {
            if idx == len60_idx {
                continue;
            }
            assert!(
                arr[(0, idx)].abs() < 1e-12,
                "expected only length 60 to be occupied, but column {idx} had {}",
                arr[(0, idx)]
            );
        }

        Ok(())
    }

    #[test]
    fn real_ref_gc_bias_then_gc_bias_package_changes_lengths_in_expected_direction() -> Result<()> {
        // Arrange:
        // Use the same real non-neutral producer setup as the corresponding `gc-bias` test:
        // - A/C split reference with pure-start BED windows on the reference side
        // - one A-only sample fragment and nine C-only sample fragments, all length 10
        //
        // The resulting real correction package is hand-derived as:
        // - weight 5.0   for the GC%=0 bin
        // - weight 5/9   for the GC%=100 bin
        //
        // `lengths` counts fragment mass by length after applying GC weight. All ten fragments
        // have the same length 10, so the only occupied output cell must be:
        //   1 * 5.0 + 9 * (5/9) = 10.0
        let reference = fixtures::twobit_from_sequences(
            "lengths_real_non_neutral_reference",
            vec![(
                "chr1".to_string(),
                format!("{}{}", "A".repeat(100), "C".repeat(100)),
            )],
        )?;
        let starts = [10_i64, 110, 120, 130, 140, 150, 160, 170, 180, 190];
        let fragments = starts
            .into_iter()
            .map(|start| fixtures::paired_fragment(start, 10, 5))
            .collect();
        let bam = bam_from_specs(
            vec![("chr1".to_string(), 200)],
            fragments,
            Vec::new(),
            "lengths_real_non_neutral_bam",
        )?;
        let out_dir = TempDir::new()?;
        let gc_path = build_real_non_neutral_gc_package(
            &bam.bam,
            &reference.path,
            out_dir.path(),
            10,
            "chr1\t0\t91\nchr1\t100\t191\n",
            // Chromosome length 200 and fragment length 10 give:
            //   200 - 10 + 1 = 191 valid starts.
            191,
        )?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(WindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
            gc_file: Some(gc_path),
            drop_invalid_gc: false,
        });
        cfg.set_ref_2bit(Some(reference.path.clone()));
        cfg.set_gc_length_weighting(MarginalizeLengthsWeightingScheme::Equal);
        {
            let fragment_lengths = cfg.fragment_lengths_mut();
            fragment_lengths.min_fragment_length = 10;
            fragment_lengths.max_fragment_length = 10;
        }

        // Act
        run(&cfg)?;

        // Assert
        let npy_path = out_dir
            .path()
            .join(dot_join(&[cfg.output_prefix.trim(), "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.dim(), (1, 1));
        assert!(
            (arr[(0, 0)] - 10.0).abs() < 1e-12,
            "expected total weighted length-10 count 10.0, got {}",
            arr[(0, 0)]
        );

        Ok(())
    }

    #[test]
    fn gc_file_rejects_package_when_fragment_length_range_is_outside_supported_range() -> Result<()> {
        // Arrange:
        // `simple_inward_bam()` contains one fragment of length 60.
        // Give `lengths` a GC package that only covers fragment lengths 10..=59:
        //   length_edges = [10, 59]
        // Restrict the command to [60, 60] so the shared GC loader must reject the package before
        // any reference lookup or length counting begins.
        let bam = simple_inward_bam()?;
        let ref_twobit = simple_reference_twobit()?;
        let out_dir = TempDir::new()?;
        let gc_path = out_dir.path().join("gc_pkg_short.npz");
        let package = GCCorrectionPackage {
            version: GC_CORRECTION_SCHEMA_VERSION,
            end_offset: 0,
            length_edges: vec![10, 59],
            gc_edges: vec![0, 101],
            length_bin_frequencies: array![1.0_f64],
            correction_matrix: array![[1.0_f64]],
        };
        package.write_npz(&gc_path)?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(WindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
            gc_file: Some(gc_path),
            drop_invalid_gc: false,
        });
        cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 60;
            frag.max_fragment_length = 60;
        }

        // Act
        let err = run(&cfg).expect_err("out-of-range GC package should fail");

        // Assert
        let msg = err.to_string();
        assert!(
            msg.contains("fragment length range [60-60] is outside the range covered by the correction package [10-59]"),
            "unexpected error message: {msg}"
        );

        Ok(())
    }

    #[test]
    fn gc_file_rejects_package_with_schema_version_mismatch() -> Result<()> {
        // Arrange:
        // Build a minimal GC correction package with an intentionally incompatible schema version.
        // `lengths` should fail while loading the package, before any reference lookup or length
        // counting begins.
        let bam = simple_inward_bam()?;
        let ref_twobit = simple_reference_twobit()?;
        let out_dir = TempDir::new()?;
        let gc_path = out_dir.path().join("gc_pkg_bad_version.npz");
        let package = GCCorrectionPackage {
            version: GC_CORRECTION_SCHEMA_VERSION + 1,
            end_offset: 0,
            length_edges: vec![10, 200],
            gc_edges: vec![0, 101],
            length_bin_frequencies: array![1.0_f64],
            correction_matrix: array![[1.0_f64]],
        };
        package.write_npz(&gc_path)?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(WindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
            gc_file: Some(gc_path),
            drop_invalid_gc: false,
        });
        cfg.set_ref_2bit(Some(ref_twobit.path.clone()));

        // Act
        let err = run(&cfg).expect_err("schema version mismatch should fail");

        // Assert
        let msg = err.to_string();
        assert!(
            msg.contains("GC correction package schema version mismatch"),
            "unexpected error message: {msg}"
        );

        Ok(())
    }

    #[test]
    fn gc_requires_ref_2bit_errors() -> Result<()> {
        let bam = simple_inward_bam()?;
        let gc_dir = TempDir::new()?;
        let gc_path = gc_dir.path().join("gc_pkg.npz");
        build_gc_package(&gc_path, 0)?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: gc_dir.path().to_path_buf(),
                n_threads: 2,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(WindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
            gc_file: Some(gc_path.clone()),
            drop_invalid_gc: false,
        });
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 200;
        }
        // Intentionally omit ref_2bit

        let err = run(&cfg).expect_err("missing ref_2bit should error");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("--ref-2bit"),
            "unexpected error message: {msg}"
        );
        Ok(())
    }

    #[test]
    fn gc_drop_invalid_reports_end_offset_validation_error() -> Result<()> {
        let bam = simple_inward_bam()?;
        let ref_twobit = simple_reference_twobit()?;
        let gc_dir = TempDir::new()?;
        let gc_path = gc_dir.path().join("gc_pkg.npz");
        // Choose large end_offset so offset_start >= offset_end, causing GC weight failure
        build_gc_package(&gc_path, 40)?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: gc_dir.path().to_path_buf(),
                n_threads: 2,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(WindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
            gc_file: Some(gc_path.clone()),
            drop_invalid_gc: true,
        });
        cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 200;
        }

        let err = run(&cfg).expect_err("should fail validation when end-offset too large");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("must exceed twice the end-offset"),
            "unexpected error: {msg}"
        );
        Ok(())
    }

    fn indel_bam_fixture() -> Result<BamFixture> {
        // Reuse edge-case fixture with one clean fragment, one insertion, one deletion.
        fixtures::fragment_kmers_edge_bam()
    }

    #[test]
    fn indel_adjust_counts_adjusted_length_and_skip_drops() -> Result<()> {
        let bam = indel_bam_fixture()?;
        let out_dir = TempDir::new()?;

        let mut base_cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 2,
            },
            base_chromosomes(&["chr1"]),
        );
        base_cfg.set_windows(WindowsArgs::default());
        base_cfg.set_window_assignment(AssignToWindowArgs::default());
        base_cfg.set_min_mapq(0);
        base_cfg.set_require_proper_pair(false);
        {
            let frag = base_cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 100;
        }

        // Adjust mode: expect all fragments counted with indel-aware lengths
        let mut adjust_cfg = base_cfg.clone();
        adjust_cfg.set_indel_mode(IndelMode::Adjust);
        run(&adjust_cfg)?;
        let npy_path = out_dir.path().join(dot_join(&[
            adjust_cfg.output_prefix.trim(),
            "length_counts.npy",
        ]));
        let arr: Array2<f64> = read_npy(&npy_path)?;
        // Expected adjusted lengths from fixture:
        //   frag0 (no indel): len 24
        //   frag1 (insertion): len 17
        //   frag2 (deletion): len 10
        let l24 = 24 - 10;
        let l17 = 17 - 10;
        let l10 = 10 - 10;
        assert!((arr[(0, l24)] - 1.0).abs() < 1e-6);
        assert!((arr[(0, l17)] - 1.0).abs() < 1e-6);
        assert!((arr[(0, l10)] - 1.0).abs() < 1e-6);
        assert!((arr.sum() - 3.0).abs() < 1e-6);

        // Skip mode: fragment carries indels, so nothing counted
        let mut skip_cfg = base_cfg.clone();
        skip_cfg.set_indel_mode(IndelMode::Skip);
        run(&skip_cfg)?;
        let skip_path = out_dir.path().join(dot_join(&[
            skip_cfg.output_prefix.trim(),
            "length_counts.npy",
        ]));
        let skip_arr: Array2<f64> = read_npy(&skip_path)?;
        // Only the indel-free fragment remains
        let l24 = 24 - 10;
        assert!((skip_arr[(0, l24)] - 1.0).abs() < 1e-6);
        assert!((skip_arr.sum() - 1.0).abs() < 1e-6);

        Ok(())
    }

    #[test]
    fn indel_adjust_bins_by_adjusted_length_but_scales_over_reference_span() -> Result<()> {
        let bam = indel_bam_fixture()?;
        let out_dir = TempDir::new()?;
        let scaling_path = out_dir.path().join("indel_adjust_scaling.tsv");
        write_scaling_factors(
            &scaling_path,
            &[("chr1", 0, 10, 1.0_f32), ("chr1", 10, 40, 3.0_f32)],
        )?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 2,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_indel_mode(IndelMode::Adjust);
        cfg.set_windows(WindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_scaling_factors(Some(scaling_path));
        {
            let frag = cfg.fragment_lengths_mut();
            // Keep only the insertion-bearing fragment from the edge fixture.
            frag.min_fragment_length = 17;
            frag.max_fragment_length = 17;
        }

        // Manual expectations:
        // - In adjust mode, the insertion fragment is counted in adjusted-length bin 17.
        // - Its reference span is still [5, 21), which is 16 bp long.
        // - Scaling is intentionally computed over that full reference span:
        //   [5, 10): 5 bp at factor 1
        //   [10, 21): 11 bp at factor 3
        //   average scaling = (5*1 + 11*3) / 16 = 38 / 16 = 2.375
        run(&cfg)?;

        let npy_path = out_dir
            .path()
            .join(dot_join(&[cfg.output_prefix.trim(), "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.shape(), &[1, 1]);
        assert!((arr[(0, 0)] - 2.375).abs() < 1e-6);
        assert!((arr.sum() - 2.375).abs() < 1e-6);

        Ok(())
    }

    #[test]
    fn indel_adjust_blacklist_uses_full_reference_span_not_only_adjusted_length() -> Result<()> {
        let bam = indel_bam_fixture()?;
        let out_dir = TempDir::new()?;
        let blacklist_path = out_dir.path().join("indel_adjust_blacklist.bed");
        fixtures::write_bed(&blacklist_path, &[("chr1", 19, 20, "deleted_base")])?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 2,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_indel_mode(IndelMode::Adjust);
        cfg.set_windows(WindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.blacklist = Some(vec![blacklist_path]);
        {
            let frag = cfg.fragment_lengths_mut();
            // Keep only the deletion-bearing fragment from the edge fixture.
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 10;
        }

        // Manual expectations:
        // - In adjust mode, the deletion fragment is counted in adjusted-length bin 10.
        // - Its full reference span is [16, 27), even though one deleted base means the
        //   adjusted length is shorter than the reference interval.
        // - Blacklisting [19, 20) therefore still overlaps the fragment and must exclude it.
        run(&cfg)?;

        let npy_path = out_dir
            .path()
            .join(dot_join(&[cfg.output_prefix.trim(), "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.shape(), &[1, 1]);
        assert!((arr.sum()).abs() < 1e-6);

        Ok(())
    }

    #[test]
    fn indel_adjust_blacklist_proportion_uses_reference_span_denominator() -> Result<()> {
        let bam = indel_bam_fixture()?;
        let out_dir = TempDir::new()?;
        let blacklist_path = out_dir.path().join("indel_adjust_blacklist_proportion.bed");
        fixtures::write_bed(&blacklist_path, &[("chr1", 19, 20, "deleted_base")])?;

        let mut base_cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 2,
            },
            base_chromosomes(&["chr1"]),
        );
        base_cfg.set_windows(WindowsArgs::default());
        base_cfg.set_window_assignment(AssignToWindowArgs::default());
        base_cfg.set_min_mapq(0);
        base_cfg.set_require_proper_pair(false);
        base_cfg.blacklist = Some(vec![blacklist_path]);
        base_cfg.blacklist_strategy = BlacklistStrategy::Proportion(0.095);

        // Manual expectations:
        // - The deletion fragment has adjusted length 10 but reference span [16, 27), length 11.
        // - The blacklist contributes exactly 1 bp of overlap: [19, 20).
        // - With proportion=0.095:
        //   expected behavior uses the full reference span:
        //   1 / 11 = 0.090909... < 0.095                        -> survives
        //   buggy adjusted-length denominator would instead use:
        //   1 / 10 = 0.1 >= 0.095                              -> blacklisted
        // - With proportion=0.09:
        //   expected behavior still uses the full reference span:
        //   1 / 11 = 0.090909... >= 0.09                        -> blacklisted
        // - Therefore both Ignore and Adjust must:
        //   survive at 0.095, and
        //   be blacklisted at 0.09.
        // - If Adjust used the adjusted length as denominator, it would already be dropped at 0.095.

        let run_sum = |indel_mode: IndelMode, threshold: f64| -> Result<f64> {
            let mut cfg = base_cfg.clone();
            cfg.set_indel_mode(indel_mode);
            cfg.blacklist_strategy = BlacklistStrategy::Proportion(threshold);
            {
                let frag = cfg.fragment_lengths_mut();
                frag.min_fragment_length = 10;
                frag.max_fragment_length = 11;
            }

            run(&cfg)?;
            let npy_path = out_dir
                .path()
                .join(dot_join(&[cfg.output_prefix.trim(), "length_counts.npy"]));
            let arr: Array2<f64> = read_npy(&npy_path)?;
            Ok(arr.sum())
        };

        let ignore_survives = run_sum(IndelMode::Ignore, 0.095)?;
        let adjust_survives = run_sum(IndelMode::Adjust, 0.095)?;
        let ignore_blacklisted = run_sum(IndelMode::Ignore, 0.09)?;
        let adjust_blacklisted = run_sum(IndelMode::Adjust, 0.09)?;

        assert!((ignore_survives - 1.0).abs() < 1e-6);
        assert!((adjust_survives - 1.0).abs() < 1e-6);
        assert!((ignore_blacklisted).abs() < 1e-6);
        assert!((adjust_blacklisted).abs() < 1e-6);

        Ok(())
    }

    #[test]
    fn scaling_overlapping_bins_error() -> Result<()> {
        let bam = simple_inward_bam()?;
        let out_dir = TempDir::new()?;

        let scaling_path = out_dir.path().join("scaling.tsv");
        // Overlapping bins: 0-150 and 100-200 on chr1 (chr len 200)
        write_scaling_factors(
            &scaling_path,
            &[("chr1", 0, 150, 1.0_f32), ("chr1", 100, 200, 1.0_f32)],
        )?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 2,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_scaling_factors(Some(scaling_path.clone()));
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_indel_mode(IndelMode::Ignore);
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 200;
        }

        let err = run(&cfg).expect_err("overlapping scaling bins should fail");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("not contiguous"),
            "unexpected error message: {msg}"
        );
        Ok(())
    }

    #[test]
    fn custom_output_prefix_is_used() -> Result<()> {
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
        cfg.output_prefix = "custom_lengths".to_string();
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_windows(WindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 200;
        }

        run(&cfg)?;

        let npy_path = out_dir
            .path()
            .join(dot_join(&[cfg.output_prefix.trim(), "length_counts.npy"]));
        assert!(npy_path.exists(), "expected {}", npy_path.display());
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.shape(), &[1, 191]);
        Ok(())
    }

    fn multi_chrom_simple_bam() -> Result<BamFixture> {
        // Different contig lengths and fragment lengths to catch duplicated chr handling
        let chroms = vec![("chr1".to_string(), 140u32), ("chr2".to_string(), 200u32)];

        let paired_chr1 = FragmentSpec {
            forward: ReadSpec {
                tid: 0,
                pos: 20,
                cigar: vec![('M', 20)],
                seq: vec![b'A'; 20],
                qual: 35,
                is_reverse: false,
                mapq: 60,
                flags: 0x1 | 0x2 | 0x40 | 0x20,
                mate_tid: Some(0),
                mate_pos: Some(60),
                insert_size: 60,
            },
            reverse: ReadSpec {
                tid: 0,
                pos: 60,
                cigar: vec![('M', 20)],
                seq: vec![b'T'; 20],
                qual: 35,
                is_reverse: true,
                mapq: 60,
                flags: 0x1 | 0x2 | 0x80,
                mate_tid: Some(0),
                mate_pos: Some(20),
                insert_size: -60,
            },
        };

        let paired_chr2 = FragmentSpec {
            forward: ReadSpec {
                tid: 1,
                pos: 40,
                cigar: vec![('M', 20)],
                seq: vec![b'C'; 20],
                qual: 35,
                is_reverse: false,
                mapq: 60,
                flags: 0x1 | 0x2 | 0x40 | 0x20,
                mate_tid: Some(1),
                mate_pos: Some(100),
                insert_size: 80,
            },
            reverse: ReadSpec {
                tid: 1,
                pos: 100,
                cigar: vec![('M', 20)],
                seq: vec![b'G'; 20],
                qual: 35,
                is_reverse: true,
                mapq: 60,
                flags: 0x1 | 0x2 | 0x80,
                mate_tid: Some(1),
                mate_pos: Some(40),
                insert_size: -80,
            },
        };

        let fragments = vec![paired_chr1, paired_chr2];
        bam_from_specs(chroms, fragments, Vec::new(), "multi_chrom_simple")
    }

    #[test]
    fn multi_chrom_size_counts_mass_conserved() -> Result<()> {
        let bam = multi_chrom_simple_bam()?;
        let out_dir = TempDir::new()?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 2,
            },
            ChromosomeArgs {
                chromosomes: Some(vec!["chr1".to_string(), "chr2".to_string()]),
                chromosomes_file: None,
            },
        );
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_tile_size(30); // force multiple tiles per chromosome
        cfg.set_indel_mode(IndelMode::Ignore);
        // Use a large size bin so each chromosome produces exactly one window, avoiding global collapse
        cfg.set_windows(WindowsArgs {
            by_size: Some(200),
            by_bed: None,
        });
        cfg.set_window_assignment(AssignToWindowArgs::default());
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 100;
        }

        run(&cfg)?;

        let npy_path = out_dir
            .path()
            .join(dot_join(&[cfg.output_prefix.trim(), "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;
        // Two chromosomes -> two windows (one per chr because by_size is large)
        assert_eq!(arr.shape(), &[2, 91]);
        let len60_idx = 60 - 10; // chr1 fragment length
        let len80_idx = 80 - 10; // chr2 fragment length
        assert!((arr[(0, len60_idx)] - 1.0).abs() < 1e-6);
        assert!((arr[(1, len80_idx)] - 1.0).abs() < 1e-6);
        assert!((arr.sum() - 2.0).abs() < 1e-6);
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
            let npy_path = out_dir
                .path()
                .join(dot_join(&[prefix, "length_counts.npy"]));
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

    #[test]
    fn scaling_tsv_must_cover_requested_chromosome_end_in_lengths() -> Result<()> {
        // Arrange:
        // `simple_inward_bam()` uses chr1 length 200.
        // A valid scaling TSV must cover every requested chromosome contiguously from 0 up to
        // the exact chromosome length. Here we intentionally stop at 100:
        //   [0,100) factor 2.0
        //
        // The fragment itself lies inside that span, but the artifact is still malformed for the
        // requested chromosome and should be rejected before any counting starts.
        let bam = simple_inward_bam()?;
        let out_dir = TempDir::new()?;
        let scaling_path = out_dir.path().join("truncated_scaling.tsv");
        write_scaling_factors(&scaling_path, &[("chr1", 0, 100, 2.0_f32)])?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(WindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_scaling_factors(Some(scaling_path));
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 200;
        }

        // Act
        let err = run(&cfg).expect_err("truncated scaling TSV should fail");

        // Assert:
        // `lengths` wraps the shared loader with `load scaling factors`, so inspect the full
        // error chain rather than only the top-level context string.
        let msg = format!("{err:#}");
        assert!(
            msg.contains("scaling TSV: bins on 'chr1' must end at chrom_len=200 (got end=100)"),
            "unexpected error message: {msg}"
        );

        Ok(())
    }

    #[cfg(feature = "cmd_coverage_weights")]
    #[test]
    fn coverage_weights_tsv_in_count_overlap_mode_uses_overlap_span_scaling() -> Result<()> {
        // Arrange:
        // Producer BAM:
        // - `simple_inward_bam()` has one fragment [20, 80).
        // - Run `coverage-weights` with `bin_size = stride = 20`, so the overlap kernel is the
        //   identity and the written scaling factors are:
        //     [0,20):  0
        //     [20,40): 1
        //     [40,60): 1
        //     [60,80): 1
        //     [80,200): 0
        //
        // Consumer BAM:
        // - One odd-length fragment [20, 81), length 61.
        //
        // Count window:
        // - BED window [45, 56), length 11.
        // - The fragment fully covers this window, so in `count-overlap` mode:
        //     overlap fraction = 11 / 61.
        //
        // Crucial scaling derivation for `lengths` in `count-overlap` mode:
        // - This mode averages scaling only over the fragment/window overlap span, not the full
        //   fragment.
        // - The overlap span is exactly [45, 56), which lies entirely inside scaling bin [40, 60)
        //   with factor 1.
        // - The scaling multiplier is therefore exactly 1.0.
        // - Final contribution to length 61 is:
        //     (overlap fraction) * (overlap-span scaling) = (11 / 61) * 1 = 11 / 61.
        //
        // If the implementation incorrectly averaged scaling over the full fragment, this test
        // would instead observe:
        //     (11 / 61) * (60 / 61),
        // because the fragment spends 60 bp in factor-1 bins and 1 bp in a factor-0 bin.
        let producer_bam = simple_inward_bam()?;
        let consumer_bam = bam_from_specs(
            vec![("chr1".to_string(), 200)],
            vec![fragment_on_tid(0, 20, 61, 20)],
            Vec::new(),
            "lengths_scaling_consumer",
        )?;
        let out_dir = TempDir::new()?;
        let weights_out_dir = out_dir.path().join("coverage_weights");
        std::fs::create_dir_all(&weights_out_dir)?;
        let scaling_cfg = make_simple_coverage_weights_config(&weights_out_dir, &producer_bam.bam);
        let bed_path = out_dir.path().join("windows.bed");
        write_bed(&bed_path, &[("chr1", 45, 56, "w0")])?;

        // Act
        run_coverage_weights(&scaling_cfg)?;
        let scaling_path = weights_out_dir.join("coverage.scaling_factors.tsv");

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: consumer_bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(WindowsArgs {
            by_size: None,
            by_bed: Some(bed_path),
        });
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::CountOverlap,
        });
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_scaling_factors(Some(scaling_path));
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 61;
            frag.max_fragment_length = 61;
        }

        run(&cfg)?;

        // Assert
        let npy_path = out_dir
            .path()
            .join(dot_join(&[cfg.output_prefix.trim(), "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.shape(), &[1, 1]);

        let expected = 11.0_f64 / 61.0_f64;
        assert!(
            (arr[(0, 0)] - expected).abs() <= 1e-9,
            "expected weighted count {expected}, got {}",
            arr[(0, 0)]
        );
        assert!(
            (arr.sum() - expected).abs() <= 1e-9,
            "expected total mass {expected}, got {}",
            arr.sum()
        );
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
    use cfdnalab::shared::interval::IndexedInterval;
    use cfdnalab::shared::tiled_run::{Tile, TileWindowSpan, build_tiles};
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

    #[test]
    fn fetch_span_size_mode_clamps_to_halo_and_chrom() {
        // Tile: core 50-150, fetch 30-200 (halo 20 left, 50 right), chrom len 180
        let tile = Tile::from_coords("chr1".to_string(), 0, 0, 50, 150, 30, 200)
            .expect("test tile should be valid");
        let span = fetch_span_for_tile(&tile, None, None, &WindowSpec::Size(100), 180)
            .expect("span expected")
            .expect("fetch span expected");
        // Window span touching core: 0..200, after halo clamp -> 30..180
        assert_eq!(span.start(), 30);
        assert_eq!(span.end(), 180);
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
            assert_eq!((t.core_start() as u64) % 10, 0);
        }
        // Expect four tiles: 0-30,30-60,60-90,90-100
        assert_eq!(tiles.len(), 4);
        assert_eq!(tiles[0].core_end(), 30);
        assert_eq!(tiles[3].core_start(), 90);
        assert_eq!(tiles[3].core_end(), 100);
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
        let tile = Tile::from_coords("chr1".to_string(), 0, 0, 0, 50, 0, 200)
            .expect("test tile should be valid");
        let span = fetch_span_for_tile(&tile, None, None, &WindowSpec::Global, 120)
            .expect("span")
            .expect("fetch span expected");
        assert_eq!(span.start(), 0);
        assert_eq!(span.end(), 120);
    }

    #[test]
    fn fetch_span_for_tile_bed_with_overlap() {
        let tile = Tile::from_coords("chr1".to_string(), 0, 0, 100, 160, 80, 200)
            .expect("test tile should be valid");
        let windows = indexed_windows(&[(90, 110, 0), (150, 170, 1), (250, 300, 2)]);
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
        .expect("span")
        .expect("fetch span expected");
        // min_ws=90, max_we=170, halos: left 20, right 40 -> widened to 70..210, clamped to fetch
        assert_eq!(res.start(), 80);
        assert_eq!(res.end(), 200);
    }

    #[test]
    fn fetch_span_bed_none_when_no_overlap() {
        let tile = Tile::from_coords("chr1".to_string(), 0, 0, 100, 150, 80, 170)
            .expect("test tile should be valid");
        // No windows overlap tile
        let windows: [IndexedInterval<u64>; 0] = [];
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
        )
        .expect("fetch span computation should succeed");
        assert!(res.is_none());
    }

    #[test]
    fn fetch_span_size_mode_none_when_tile_right_of_chromosome() {
        let tile = Tile::from_coords("chr1".to_string(), 0, 0, 250, 260, 230, 270)
            .expect("test tile should be valid");
        let res = fetch_span_for_tile(&tile, None, None, &WindowSpec::Size(50), 200)
            .expect("fetch span computation should succeed");
        assert!(res.is_none());
    }

    #[test]
    fn tile_constructor_rejects_empty_core() {
        let err = Tile::from_coords("chr1".to_string(), 0, 0, 100, 100, 80, 120).unwrap_err();
        assert!(format!("{err}").contains("interval end (100) must be greater than start (100)"));
    }
}
