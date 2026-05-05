#![cfg(feature = "cmd_lengths")]

mod fixtures;

mod tests_lengths_command {

    use super::*;

    use anyhow::Result;
    use cfdnalab::commands::cli_common::{
        AssignToWindowArgs, ChromosomeArgs, DistributionWindowsArgs, IOCArgs, UnpairedArgs,
        WindowAssigner,
    };
    #[cfg(feature = "cmd_coverage_weights")]
    use cfdnalab::commands::coverage_weights::{
        config::CoverageWeightsConfig, coverage_weights::run as run_coverage_weights,
    };
    use cfdnalab::commands::gc_bias::{
        correct::{GCLengthRange, MarginalizeLengthsWeightingScheme},
        package::GCCorrectionPackage,
    };
    use cfdnalab::commands::lengths::config::LengthsConfig;
    use cfdnalab::commands::lengths::lengths::run;
    use cfdnalab::shared::blacklist::strategy::BlacklistStrategy;
    use cfdnalab::shared::constants::GC_CORRECTION_SCHEMA_VERSION;
    use cfdnalab::shared::io::dot_join;
    use cfdnalab::shared::{
        clip_mode::ClipMode,
        indel_mode::IndelMode,
        reference::{ContigFootprintEntry, twobit_contig_footprint},
    };
    use fixtures::{
        BamFixture, FragmentSpec, ReadSpec, bam_from_specs, build_real_neutral_gc_package,
        build_real_neutral_gc_package_for_range, build_real_non_neutral_gc_package,
        late_origin_gc_reference_sequence, simple_inward_bam, simple_reference_twobit,
        twobit_from_sequences, write_bed, write_scaling_factors, write_two_bin_gc_package,
    };
    use ndarray::Array2;
    use ndarray::array;
    use ndarray_npy::read_npy;
    use serde_json::Value;
    use tempfile::TempDir;

    fn parse_group_index_rows(text: &str) -> Vec<(u64, String, f64)> {
        let mut lines = text.lines();
        let header = lines.next().expect("group index TSV must have a header");
        assert_eq!(header, "group_idx\tgroup_name\tblacklisted_fraction");

        lines
            .map(|line| {
                let mut fields = line.split('\t');
                let group_idx = fields
                    .next()
                    .expect("group index row must contain group_idx")
                    .parse::<u64>()
                    .expect("group_idx must parse as u64");
                let group_name = fields
                    .next()
                    .expect("group index row must contain group_name")
                    .to_string();
                let blacklisted_fraction = fields
                    .next()
                    .expect("group index row must contain blacklisted_fraction")
                    .parse::<f64>()
                    .expect("blacklisted_fraction must parse as f64");
                assert!(
                    fields.next().is_none(),
                    "group index row must contain exactly three columns"
                );
                (group_idx, group_name, blacklisted_fraction)
            })
            .collect()
    }

    fn parse_group_index_tsv(text: &str) -> Vec<(u64, String)> {
        let mut lines = text.lines();
        let header = lines.next().expect("group index TSV must have a header");
        let expected_column_count = match header {
            "group_idx\tgroup_name" => 2,
            "group_idx\tgroup_name\tblacklisted_fraction" => 3,
            _ => panic!("unexpected group index TSV header: {header}"),
        };

        lines
            .map(|line| {
                let fields: Vec<&str> = line.split('\t').collect();
                assert_eq!(
                    fields.len(),
                    expected_column_count,
                    "group index row must match the header column count"
                );
                let group_idx = fields[0]
                    .parse::<u64>()
                    .expect("group_idx must parse as u64");
                let group_name = fields[1].to_string();
                (group_idx, group_name)
            })
            .collect()
    }

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

    fn single_read_fragment_bam(
        name: &str,
        fragment_start: i64,
        fragment_len: u32,
    ) -> Result<BamFixture> {
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

    fn single_read_fragment_bam_with_cigar(
        name: &str,
        fragment_start: i64,
        cigar: Vec<(char, u32)>,
        seq: Vec<u8>,
    ) -> Result<BamFixture> {
        bam_from_specs(
            vec![("chr1".to_string(), 200)],
            Vec::new(),
            vec![ReadSpec {
                tid: 0,
                pos: fragment_start,
                cigar,
                seq,
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
    fn default_length_bins_preserve_per_bp_distribution_shape() -> Result<()> {
        // Arrange:
        // The default `--length-bins 30:1001:1` creates one column for each integer
        // length from 30 through 1000. The simple fixture has one length-60 fragment.
        let bam = simple_inward_bam()?;
        let out_dir = TempDir::new()?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);

        // Act
        run(&cfg)?;

        // Assert
        let counts_path = out_dir.path().join("length_counts.npy");
        let arr: Array2<f64> = read_npy(&counts_path)?;
        assert_eq!(arr.shape(), &[1, 971]);
        assert_eq!(arr[(0, 60 - 30)], 1.0);
        assert_eq!(arr.sum(), 1.0);

        let settings_text =
            std::fs::read_to_string(out_dir.path().join("fragment_length_settings.json"))?;
        let settings: Value = serde_json::from_str(&settings_text)?;
        assert_eq!(settings["length_axis"]["column_intervals"], "half_open");
        assert_eq!(settings["length_axis"]["min_fragment_length"], 30);
        assert_eq!(settings["length_axis"]["max_fragment_length"], 1000);
        assert_eq!(settings["length_axis"]["n_bins"], 971);
        assert_eq!(settings["length_axis"]["single_bp_bins"], true);
        assert_eq!(
            settings["length_axis"]["bin_definition"]["kind"],
            "stepped_range"
        );
        assert_eq!(settings["length_axis"]["bin_definition"]["start"], 30);
        assert_eq!(settings["length_axis"]["bin_definition"]["end"], 1001);
        assert_eq!(settings["length_axis"]["bin_definition"]["step"], 1);
        assert!(settings["length_axis"].get("edges").is_none());
        assert_eq!(settings["aggregation_level"], "global");
        assert!(settings.get("row_semantics").is_none());
        assert_eq!(settings["window_mode"], "global");
        assert_eq!(settings["indel_mode"], "ignore");
        assert_eq!(settings["clip_mode"], "aligned");
        assert_eq!(settings["max_soft_clips"], 256);
        assert_eq!(settings["max_deletion_bases"], 100);
        assert_eq!(settings["assign_by"], "count-overlap");
        assert_eq!(settings["gc_length_weighting"], "equal");
        assert_eq!(settings["gc_length_range"], "requested");
        assert_eq!(settings["gc_length_trim_rare"], 0.0);
        assert_eq!(settings["gc_correction_used"], false);
        assert_eq!(settings["scaling_factors_used"], false);
        assert!(settings.get("min_mapq").is_none());
        assert!(settings.get("require_proper_pair").is_none());
        Ok(())
    }

    #[test]
    fn wider_length_bins_collapse_multiple_lengths_into_one_column() -> Result<()> {
        // Arrange:
        // Two fragments have different exact lengths but both fall in [30,40).
        let bam = bam_from_specs(
            vec![("chr1".to_string(), 200)],
            vec![
                fixtures::paired_fragment(20, 35, 10),
                fixtures::paired_fragment(80, 39, 10),
            ],
            Vec::new(),
            "lengths_wider_bins_collapse",
        )?;
        let out_dir = TempDir::new()?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_length_bins(vec![30, 40, 50]);

        // Act
        run(&cfg)?;

        // Assert
        let arr: Array2<f64> = read_npy(out_dir.path().join("length_counts.npy"))?;
        assert_eq!(arr.shape(), &[1, 2]);
        assert_eq!(arr[(0, 0)], 2.0);
        assert_eq!(arr[(0, 1)], 0.0);
        assert_eq!(arr.sum(), 2.0);
        Ok(())
    }

    #[test]
    fn length_bins_filter_at_final_exclusive_edge() -> Result<()> {
        // Arrange:
        // [10,20) includes length 19 and excludes length 20.
        let bam = bam_from_specs(
            vec![("chr1".to_string(), 200)],
            vec![
                fixtures::paired_fragment(20, 19, 10),
                fixtures::paired_fragment(80, 20, 10),
            ],
            Vec::new(),
            "lengths_final_exclusive_edge",
        )?;
        let out_dir = TempDir::new()?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_length_bins(vec![10, 20]);

        // Act
        run(&cfg)?;

        // Assert
        let arr: Array2<f64> = read_npy(out_dir.path().join("length_counts.npy"))?;
        assert_eq!(arr.shape(), &[1, 1]);
        assert_eq!(arr[(0, 0)], 1.0);
        assert_eq!(arr.sum(), 1.0);
        Ok(())
    }

    #[test]
    fn counts_reference_lengths_global_window() -> Result<()> {
        // Human verification status: unverified
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
        cfg.set_windows(DistributionWindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 200);

        run(&cfg)?;

        let prefix = cfg.output_prefix.trim();
        let npy_path = out_dir
            .path()
            .join(dot_join(&[prefix, "length_counts.npy"]));
        assert!(npy_path.exists());
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.shape(), &[1, 191]);
        let len60_idx = 60 - 10; // min_fragment_length
        let row = arr.row(0);
        for (length_index, &value) in row.iter().enumerate() {
            let expected_value = if length_index == len60_idx { 1.0 } else { 0.0 };
            assert!(
                (value - expected_value).abs() < 1e-6,
                "expected value {expected_value} at length index {length_index}, got {value}"
            );
        }
        assert!((row.sum() - 1.0).abs() < 1e-6);

        Ok(())
    }

    #[test]
    fn counts_reference_lengths_size_single_window_misaligned_tiles() -> Result<()> {
        // Human verification status: unverified
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
        cfg.set_windows(DistributionWindowsArgs {
            by_size: Some(500),
            by_bed: None,
            by_grouped_bed: None,
        });
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_tile_size(50);
        cfg.set_per_bp_length_bins(10, 200);

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
        let row = arr.row(0);
        for (length_index, &value) in row.iter().enumerate() {
            let expected_value = if length_index == len60_idx { 1.0 } else { 0.0 };
            assert!(
                (value - expected_value).abs() < 1e-6,
                "expected value {expected_value} at length index {length_index}, got {value}"
            );
        }
        assert!((row.sum() - 1.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn counts_reference_lengths_size_aligned_tiles_reduce_cross_tile_bins() -> Result<()> {
        // Human verification status: unverified
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
        cfg.set_windows(DistributionWindowsArgs {
            by_size: Some(10),
            by_bed: None,
            by_grouped_bed: None,
        });
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_tile_size(100);
        cfg.set_per_bp_length_bins(10, 50);

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
        // Tile 0's fetch halo still sees the fragment start at 95, so tile 0 processes the
        // fragment and assigns overlap mass into windows [100,110), [110,120), and [120,130)
        // even though those windows begin in tile 1. The test verifies that the tile-1 pass
        // does not double-count those middle windows.
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
        // Human verification status: unverified
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
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(bed_path.clone()),
            by_grouped_bed: None,
        });
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 200);

        run(&cfg)?;

        let prefix = cfg.output_prefix.trim();
        let npy_path = out_dir
            .path()
            .join(dot_join(&[prefix, "length_counts.npy"]));
        assert!(npy_path.exists());
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.shape(), &[1, 191]);
        let len60_idx = 60 - 10;
        let row = arr.row(0);
        for (length_index, &value) in row.iter().enumerate() {
            let expected_value = if length_index == len60_idx { 1.0 } else { 0.0 };
            assert!(
                (value - expected_value).abs() < 1e-6,
                "expected value {expected_value} at length index {length_index}, got {value}"
            );
        }
        assert!((row.sum() - 1.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn bed_windowing_counts_a_right_halo_only_window_reached_by_an_owned_fragment() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // - one unpaired fragment span [19,29), length 10
        // - tile size 10 gives the owning tile core [10,20)
        // - BED window [28,29) sits entirely in the right halo, not in the core
        // - with assign-by=any, that halo window should still receive the fragment once
        let bam = single_read_fragment_bam("lengths_right_halo_only_bed", 19, 10)?;
        let out_dir = TempDir::new()?;
        let bed_path = out_dir.path().join("right_halo_only.bed");
        write_bed(&bed_path, &[("chr1", 28, 29, "halo_only")])?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_unpaired(cfdnalab::commands::cli_common::UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(bed_path),
            by_grouped_bed: None,
        });
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::Any,
        });
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_tile_size(10);
        cfg.set_per_bp_length_bins(10, 10);

        // Act
        run(&cfg)?;

        // Assert
        let prefix = cfg.output_prefix.trim();
        let npy_path = out_dir
            .path()
            .join(dot_join(&[prefix, "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.shape(), &[1, 1]);
        assert!((arr[(0, 0)] - 1.0).abs() < 1e-6);
        assert!((arr.sum() - 1.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn bed_windowing_does_not_count_a_window_starting_at_fragment_end() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // - one unpaired fragment span [19,29), length 10
        // - tile size 10 again gives the owning tile core [10,20)
        // - BED window [29,30) still sits inside the reachable right-side candidate envelope
        // - but the fragment ends exactly at 29, so the half-open overlap is zero
        let bam = single_read_fragment_bam("lengths_right_boundary_open", 19, 10)?;
        let out_dir = TempDir::new()?;
        let bed_path = out_dir.path().join("right_boundary_open.bed");
        write_bed(&bed_path, &[("chr1", 29, 30, "touches_end_only")])?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_unpaired(cfdnalab::commands::cli_common::UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(bed_path),
            by_grouped_bed: None,
        });
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::Any,
        });
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_tile_size(10);
        cfg.set_per_bp_length_bins(10, 10);

        // Act
        run(&cfg)?;

        // Assert
        let prefix = cfg.output_prefix.trim();
        let npy_path = out_dir
            .path()
            .join(dot_join(&[prefix, "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.shape(), &[1, 1]);
        assert_eq!(arr[(0, 0)], 0.0);
        assert_eq!(arr.sum(), 0.0);
        Ok(())
    }

    #[test]
    fn bed_windowing_must_not_shrink_fetch_to_unrelated_core_windows() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // - one unpaired fragment span [19,29), length 10
        // - tile size 10 again gives the owning tile core [10,20)
        // - BED window [10,11) overlaps the core but not the fragment
        // - BED window [28,29) is the real target window in the right halo
        // The correct result is two output rows with counts [0.0, 1.0].
        let bam = single_read_fragment_bam("lengths_right_halo_with_core_window", 19, 10)?;
        let out_dir = TempDir::new()?;
        let bed_path = out_dir.path().join("mixed_windows.bed");
        write_bed(
            &bed_path,
            &[
                ("chr1", 10, 11, "core_only"),
                ("chr1", 28, 29, "halo_target"),
            ],
        )?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_unpaired(cfdnalab::commands::cli_common::UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(bed_path),
            by_grouped_bed: None,
        });
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::Any,
        });
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_tile_size(10);
        cfg.set_per_bp_length_bins(10, 10);

        // Act
        run(&cfg)?;

        // Assert
        let prefix = cfg.output_prefix.trim();
        let npy_path = out_dir
            .path()
            .join(dot_join(&[prefix, "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.shape(), &[2, 1]);
        assert_eq!(arr[(0, 0)], 0.0);
        assert!((arr[(1, 0)] - 1.0).abs() < 1e-6);
        assert!((arr.sum() - 1.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn bed_windowing_right_halo_only_count_is_tile_size_invariant() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // - one unpaired fragment span [19,29), length 10
        // - BED window [28,29) is the only counted row
        // - with tile_size=10 the window lies outside the owning tile core; with tile_size=1000
        //   it lies in the same tile
        // The final BED output must be identical across those decompositions.
        let bam = single_read_fragment_bam("lengths_right_halo_tile_invariance", 19, 10)?;
        let tile_sizes = [10_u32, 1_000_u32];
        let mut outputs = Vec::new();

        for tile_size in tile_sizes {
            let out_dir = TempDir::new()?;
            let bed_path = out_dir.path().join(format!("right_halo_{tile_size}.bed"));
            write_bed(&bed_path, &[("chr1", 28, 29, "halo_only")])?;

            let mut cfg = LengthsConfig::new(
                IOCArgs {
                    bam: bam.bam.clone(),
                    output_dir: out_dir.path().to_path_buf(),
                    n_threads: 1,
                },
                base_chromosomes(&["chr1"]),
            );
            cfg.set_unpaired(cfdnalab::commands::cli_common::UnpairedArgs {
                reads_are_fragments: true,
            });
            cfg.set_indel_mode(IndelMode::Ignore);
            cfg.set_windows(DistributionWindowsArgs {
                by_size: None,
                by_bed: Some(bed_path),
                by_grouped_bed: None,
            });
            cfg.set_window_assignment(AssignToWindowArgs {
                assign_by: WindowAssigner::Any,
            });
            cfg.set_min_mapq(0);
            cfg.set_require_proper_pair(false);
            cfg.set_tile_size(tile_size);
            cfg.set_per_bp_length_bins(10, 10);

            run(&cfg)?;

            let prefix = cfg.output_prefix.trim();
            let npy_path = out_dir
                .path()
                .join(dot_join(&[prefix, "length_counts.npy"]));
            outputs.push(read_npy::<_, Array2<f64>>(&npy_path)?);
        }

        assert_eq!(outputs[0], outputs[1]);
        assert_eq!(outputs[0].shape(), &[1, 1]);
        assert_eq!(outputs[0][[0, 0]], 1.0);
        Ok(())
    }

    #[test]
    fn bed_windowing_mixed_core_and_right_halo_rows_are_tile_size_invariant() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // - one unpaired fragment span [19,29), length 10
        // - BED row [10,11) is upstream and must stay zero
        // - BED row [28,29) is the true counted row and must equal 1
        // - with tile_size=10 these rows fall in different tiles; with tile_size=1000 they do not
        let bam = single_read_fragment_bam("lengths_mixed_bed_tile_invariance", 19, 10)?;
        let tile_sizes = [10_u32, 1_000_u32];
        let mut outputs = Vec::new();

        for tile_size in tile_sizes {
            let out_dir = TempDir::new()?;
            let bed_path = out_dir.path().join(format!("mixed_{tile_size}.bed"));
            write_bed(
                &bed_path,
                &[
                    ("chr1", 10, 11, "core_only"),
                    ("chr1", 28, 29, "halo_target"),
                ],
            )?;

            let mut cfg = LengthsConfig::new(
                IOCArgs {
                    bam: bam.bam.clone(),
                    output_dir: out_dir.path().to_path_buf(),
                    n_threads: 1,
                },
                base_chromosomes(&["chr1"]),
            );
            cfg.set_unpaired(cfdnalab::commands::cli_common::UnpairedArgs {
                reads_are_fragments: true,
            });
            cfg.set_indel_mode(IndelMode::Ignore);
            cfg.set_windows(DistributionWindowsArgs {
                by_size: None,
                by_bed: Some(bed_path),
                by_grouped_bed: None,
            });
            cfg.set_window_assignment(AssignToWindowArgs {
                assign_by: WindowAssigner::Any,
            });
            cfg.set_min_mapq(0);
            cfg.set_require_proper_pair(false);
            cfg.set_tile_size(tile_size);
            cfg.set_per_bp_length_bins(10, 10);

            run(&cfg)?;

            let prefix = cfg.output_prefix.trim();
            let npy_path = out_dir
                .path()
                .join(dot_join(&[prefix, "length_counts.npy"]));
            outputs.push(read_npy::<_, Array2<f64>>(&npy_path)?);
        }

        assert_eq!(outputs[0], outputs[1]);
        assert_eq!(outputs[0].shape(), &[2, 1]);
        assert_eq!(outputs[0][[0, 0]], 0.0);
        assert_eq!(outputs[0][[1, 0]], 1.0);
        Ok(())
    }

    #[test]
    fn bed_windowing_right_boundary_zero_is_tile_size_invariant() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // - one unpaired fragment span [19,29), length 10
        // - BED window [29,30) touches fragment end only and must stay zero
        // - tile_size changes the decomposition but not the final half-open overlap result
        let bam = single_read_fragment_bam("lengths_right_boundary_tile_invariance", 19, 10)?;
        let tile_sizes = [10_u32, 1_000_u32];
        let mut outputs = Vec::new();

        for tile_size in tile_sizes {
            let out_dir = TempDir::new()?;
            let bed_path = out_dir.path().join(format!("boundary_{tile_size}.bed"));
            write_bed(&bed_path, &[("chr1", 29, 30, "touches_end_only")])?;

            let mut cfg = LengthsConfig::new(
                IOCArgs {
                    bam: bam.bam.clone(),
                    output_dir: out_dir.path().to_path_buf(),
                    n_threads: 1,
                },
                base_chromosomes(&["chr1"]),
            );
            cfg.set_unpaired(cfdnalab::commands::cli_common::UnpairedArgs {
                reads_are_fragments: true,
            });
            cfg.set_indel_mode(IndelMode::Ignore);
            cfg.set_windows(DistributionWindowsArgs {
                by_size: None,
                by_bed: Some(bed_path),
                by_grouped_bed: None,
            });
            cfg.set_window_assignment(AssignToWindowArgs {
                assign_by: WindowAssigner::Any,
            });
            cfg.set_min_mapq(0);
            cfg.set_require_proper_pair(false);
            cfg.set_tile_size(tile_size);
            cfg.set_per_bp_length_bins(10, 10);

            run(&cfg)?;

            let prefix = cfg.output_prefix.trim();
            let npy_path = out_dir
                .path()
                .join(dot_join(&[prefix, "length_counts.npy"]));
            outputs.push(read_npy::<_, Array2<f64>>(&npy_path)?);
        }

        assert_eq!(outputs[0], outputs[1]);
        assert_eq!(outputs[0].shape(), &[1, 1]);
        assert_eq!(outputs[0][[0, 0]], 0.0);
        Ok(())
    }

    #[test]
    fn global_by_size_and_bed_full_chromosome_windows_match_exactly() -> Result<()> {
        // Human verification status: unverified
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

        let make_cfg = |out_dir: &std::path::Path, windows: DistributionWindowsArgs| {
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
            cfg.set_per_bp_length_bins(10, 200);
            cfg
        };

        let global_cfg = make_cfg(global_out.path(), DistributionWindowsArgs::default());
        let by_size_cfg = make_cfg(
            by_size_out.path(),
            DistributionWindowsArgs {
                by_size: Some(200),
                by_bed: None,
                by_grouped_bed: None,
            },
        );
        let by_bed_cfg = make_cfg(
            by_bed_out.path(),
            DistributionWindowsArgs {
                by_size: None,
                by_bed: Some(bed_path),
                by_grouped_bed: None,
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
    fn lengths_default_min_mapq_matches_explicit_thirty_and_differs_from_explicit_zero()
    -> Result<()> {
        // Human verification status: unverified
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
            cfg.set_windows(DistributionWindowsArgs::default());
            cfg.set_window_assignment(AssignToWindowArgs::default());
            cfg.set_require_proper_pair(false);
            cfg.set_per_bp_length_bins(50, 120);
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
        // Human verification status: unverified
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
        cfg.set_windows(DistributionWindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 120);

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
        let row = arr.row(0);
        let occupied_indices = [60 - 10, 80 - 10, 100 - 10];
        for (length_index, &value) in row.iter().enumerate() {
            let expected_value = if occupied_indices.contains(&length_index) {
                1.0
            } else {
                0.0
            };
            assert!(
                (value - expected_value).abs() < 1e-6,
                "expected value {expected_value} at length index {length_index}, got {value}"
            );
        }
        assert!((row.sum() - 3.0).abs() < 1e-6);

        Ok(())
    }

    #[test]
    fn unpaired_single_read_matches_paired_fragment_length_count_for_same_span() -> Result<()> {
        // Human verification status: unverified
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
            cfg.set_windows(DistributionWindowsArgs::default());
            cfg.set_window_assignment(AssignToWindowArgs::default());
            cfg.set_min_mapq(0);
            cfg.set_require_proper_pair(false);
            cfg.set_unpaired(cfdnalab::commands::cli_common::UnpairedArgs {
                reads_are_fragments: unpaired,
            });
            cfg.set_per_bp_length_bins(10, 200);
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
        let unpaired_arr: Array2<f64> = read_npy(
            unpaired_out
                .path()
                .join(dot_join(&["", "length_counts.npy"])),
        )?;

        let len60_idx = 60 - 10;
        assert_eq!(paired_arr, unpaired_arr);
        assert_eq!(paired_arr.shape(), &[1, 191]);
        assert!((paired_arr[(0, len60_idx)] - 1.0).abs() < 1e-12);
        assert!((paired_arr.sum() - 1.0).abs() < 1e-12);

        Ok(())
    }

    #[test]
    fn counts_reference_lengths_size_single_window_across_three_chromosomes() -> Result<()> {
        // Human verification status: unverified
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
        cfg.set_windows(DistributionWindowsArgs {
            by_size: Some(200),
            by_bed: None,
            by_grouped_bed: None,
        });
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_tile_size(30);
        cfg.set_per_bp_length_bins(10, 120);

        run(&cfg)?;

        let prefix = cfg.output_prefix.trim();
        let npy_path = out_dir
            .path()
            .join(dot_join(&[prefix, "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;

        assert_eq!(arr.shape(), &[3, 111]);
        let expected_rows = [60 - 10, 80 - 10, 100 - 10];
        for (row_index, row) in arr.outer_iter().enumerate() {
            for (length_index, &value) in row.iter().enumerate() {
                let expected_value = if length_index == expected_rows[row_index] {
                    1.0
                } else {
                    0.0
                };
                assert!(
                    (value - expected_value).abs() < 1e-6,
                    "row {row_index}: expected value {expected_value} at length index {length_index}, got {value}"
                );
            }
            assert!((row.sum() - 1.0).abs() < 1e-6);
        }
        assert!((arr.sum() - 3.0).abs() < 1e-6);

        Ok(())
    }

    #[test]
    fn counts_reference_lengths_bed_single_window_across_three_chromosomes() -> Result<()> {
        // Human verification status: unverified
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
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(bed_path),
            by_grouped_bed: None,
        });
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 120);

        run(&cfg)?;

        let prefix = cfg.output_prefix.trim();
        let npy_path = out_dir
            .path()
            .join(dot_join(&[prefix, "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;

        assert_eq!(arr.shape(), &[3, 111]);
        let expected_rows = [60 - 10, 80 - 10, 100 - 10];
        for (row_index, row) in arr.outer_iter().enumerate() {
            for (length_index, &value) in row.iter().enumerate() {
                let expected_value = if length_index == expected_rows[row_index] {
                    1.0
                } else {
                    0.0
                };
                assert!(
                    (value - expected_value).abs() < 1e-6,
                    "row {row_index}: expected value {expected_value} at length index {length_index}, got {value}"
                );
            }
            assert!((row.sum() - 1.0).abs() < 1e-6);
        }
        assert!((arr.sum() - 3.0).abs() < 1e-6);

        Ok(())
    }

    #[test]
    fn by_size_and_bed_equivalent_full_chromosome_windows_match_across_three_chromosomes()
    -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // - three chromosomes of length 200
        // - one fragment per chromosome with lengths 60, 80, and 100
        // - by-size 200 creates one full-chromosome row per chromosome
        // - BED windows [0,200) for each chromosome describe the exact same row partition
        // The two modes must therefore produce identical row-wise length matrices.
        let bam = three_chrom_length_fixture("lengths_three_chr_bed_vs_size")?;
        let by_size_out = TempDir::new()?;
        let bed_out = TempDir::new()?;
        let bed_path = bed_out.path().join("windows_three_chr_full.bed");
        fixtures::write_bed(
            &bed_path,
            &[
                ("chr1", 0, 200, "chr1_window"),
                ("chr2", 0, 200, "chr2_window"),
                ("chr3", 0, 200, "chr3_window"),
            ],
        )?;

        let make_cfg = |output_dir: &std::path::Path, windows: DistributionWindowsArgs| {
            let mut cfg = LengthsConfig::new(
                IOCArgs {
                    bam: bam.bam.clone(),
                    output_dir: output_dir.to_path_buf(),
                    n_threads: 1,
                },
                base_chromosomes(&["chr1", "chr2", "chr3"]),
            );
            cfg.set_indel_mode(IndelMode::Ignore);
            cfg.set_windows(windows);
            cfg.set_window_assignment(AssignToWindowArgs::default());
            cfg.set_min_mapq(0);
            cfg.set_require_proper_pair(false);
            cfg.set_tile_size(50);
            cfg.set_per_bp_length_bins(10, 120);
            cfg
        };

        let by_size_cfg = make_cfg(
            by_size_out.path(),
            DistributionWindowsArgs {
                by_size: Some(200),
                by_bed: None,
                by_grouped_bed: None,
            },
        );
        let bed_cfg = make_cfg(
            bed_out.path(),
            DistributionWindowsArgs {
                by_size: None,
                by_bed: Some(bed_path),
                by_grouped_bed: None,
            },
        );

        // Act
        run(&by_size_cfg)?;
        run(&bed_cfg)?;

        // Assert
        let read_counts = |dir: &std::path::Path| -> Result<Array2<f64>> {
            Ok(read_npy(dir.join(dot_join(&["", "length_counts.npy"])))?)
        };
        let by_size_arr = read_counts(by_size_out.path())?;
        let bed_arr = read_counts(bed_out.path())?;

        assert_eq!(by_size_arr.shape(), &[3, 111]);
        assert_eq!(bed_arr.shape(), &[3, 111]);
        assert_eq!(by_size_arr, bed_arr);
        assert!((by_size_arr.sum() - 3.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn counts_apply_scaling_factors() -> Result<()> {
        // Human verification status: unverified
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
        cfg.set_windows(DistributionWindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_scaling_factors(Some(scaling_path.clone()));
        cfg.set_per_bp_length_bins(10, 200);

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
        // Human verification status: unverified
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
        cfg.set_windows(DistributionWindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.blacklist = Some(vec![blacklist_path.clone()]);
        cfg.set_per_bp_length_bins(10, 200);

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

    fn build_gc_package(
        path: &std::path::Path,
        end_offset: u64,
        reference_contig_footprint: Vec<ContigFootprintEntry>,
    ) -> Result<()> {
        // Two length bins: [10,60) and [60,200]; two GC bins: [0,50) and [50,101]
        let correction_matrix = array![[1.0_f64, 1.0_f64], [2.0_f64, 10.0_f64]];
        let package = GCCorrectionPackage {
            version: GC_CORRECTION_SCHEMA_VERSION,
            end_offset,
            length_edges: vec![10, 60, 200],
            gc_edges: vec![0, 50, 101],
            length_bin_frequencies: array![1.0_f64, 3.0_f64],
            reference_contig_footprint,
            correction_matrix,
        };
        package.write_npz(path)?;
        Ok(())
    }

    #[test]
    fn applies_gc_correction_weighting_modes() -> Result<()> {
        // Human verification status: unverified
        let bam = simple_inward_bam()?;
        let ref_twobit = simple_reference_twobit()?;
        let gc_dir = TempDir::new()?;
        let gc_path = gc_dir.path().join("gc_pkg.npz");
        build_gc_package(&gc_path, 0, twobit_contig_footprint(&ref_twobit.path)?)?;

        let expected = |scheme: MarginalizeLengthsWeightingScheme| -> f64 {
            match scheme {
                MarginalizeLengthsWeightingScheme::Equal => 5.5, // mean of rows at GC bin 50
                MarginalizeLengthsWeightingScheme::Frequency => 7.75, // weighted by [1,3]
                MarginalizeLengthsWeightingScheme::MaxFrequency => 10.0, // most frequent row
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
            cfg.set_windows(DistributionWindowsArgs::default());
            cfg.set_window_assignment(AssignToWindowArgs::default());
            cfg.set_min_mapq(0);
            cfg.set_require_proper_pair(false);
            cfg.set_gc_length_weighting(scheme);
            cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
                gc_file: Some(gc_path.clone()),
                neutralize_invalid_gc: false,
            });
            cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
            cfg.set_per_bp_length_bins(10, 200);

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
            MarginalizeLengthsWeightingScheme::Frequency,
            MarginalizeLengthsWeightingScheme::MaxFrequency,
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
    fn gc_length_trim_rare_removes_low_frequency_rows_before_equal_weighting() -> Result<()> {
        // Package rows:
        // - [10,60): GC bin 50 correction = 1, frequency = 1
        // - [60,200]: GC bin 50 correction = 10, frequency = 3
        //
        // With `--gc-length-trim-rare 0.25`, the trim budget is 25% of total
        // frequency: (1 + 3) * 0.25 = 1. The rare row with frequency 1 is
        // removed exactly. Equal weighting over the retained rows therefore
        // uses only correction 10.
        //
        // The simple fixture contains one 60 bp fragment with GC%=50, so the
        // output cell for length 60 should be 1 fragment * correction 10 = 10.
        let bam = simple_inward_bam()?;
        let ref_twobit = simple_reference_twobit()?;
        let gc_dir = TempDir::new()?;
        let gc_path = gc_dir.path().join("gc_pkg_trim_rare.npz");
        build_gc_package(&gc_path, 0, twobit_contig_footprint(&ref_twobit.path)?)?;
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
        cfg.set_windows(DistributionWindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_gc_length_weighting(MarginalizeLengthsWeightingScheme::Equal);
        cfg.set_gc_length_trim_rare(0.25);
        cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
            gc_file: Some(gc_path),
            neutralize_invalid_gc: false,
        });
        cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
        cfg.set_per_bp_length_bins(10, 200);

        run(&cfg)?;

        let prefix = cfg.output_prefix.trim();
        let npy_path = out_dir
            .path()
            .join(dot_join(&[prefix, "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;
        let len60_idx = 60 - 10;
        assert!(
            (arr[(0, len60_idx)] - 10.0).abs() < 1e-6,
            "expected rare-row-trimmed GC correction to give length-60 count 10, got {}",
            arr[(0, len60_idx)]
        );

        Ok(())
    }

    #[test]
    fn gc_length_range_controls_which_package_rows_are_marginalized() -> Result<()> {
        // Human verification status: unverified
        // Package rows:
        // - [10,60): GC bin 50 correction = 1
        // - [60,200]: GC bin 50 correction = 10
        //
        // The simple fixture has one 60 bp fragment with GC%=50.
        // With requested range [60,60], default `requested` selection uses only the second row.
        // With `package`, equal weighting uses both rows: (1 + 10) / 2 = 5.5.
        let bam = simple_inward_bam()?;
        let ref_twobit = simple_reference_twobit()?;
        let gc_dir = TempDir::new()?;
        let gc_path = gc_dir.path().join("gc_pkg_range.npz");
        build_gc_package(&gc_path, 0, twobit_contig_footprint(&ref_twobit.path)?)?;

        let run_with_range = |gc_length_range: GCLengthRange| -> Result<f64> {
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
            cfg.set_windows(DistributionWindowsArgs::default());
            cfg.set_window_assignment(AssignToWindowArgs::default());
            cfg.set_min_mapq(0);
            cfg.set_require_proper_pair(false);
            cfg.set_gc_length_range(gc_length_range);
            cfg.set_gc_length_weighting(MarginalizeLengthsWeightingScheme::Equal);
            cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
                gc_file: Some(gc_path.clone()),
                neutralize_invalid_gc: false,
            });
            cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
            cfg.set_per_bp_length_bins(60, 60);

            run(&cfg)?;

            let prefix = cfg.output_prefix.trim();
            let npy_path = out_dir
                .path()
                .join(dot_join(&[prefix, "length_counts.npy"]));
            let arr: Array2<f64> = read_npy(&npy_path)?;
            Ok(arr[(0, 0)])
        };

        let requested_value = run_with_range(GCLengthRange::Requested)?;
        let package_value = run_with_range(GCLengthRange::Package)?;
        assert!(
            (requested_value - 10.0).abs() < 1e-6,
            "requested range should use only the 60 bp package row, got {requested_value}"
        );
        assert!(
            (package_value - 5.5).abs() < 1e-6,
            "package range should average both package rows, got {package_value}"
        );

        Ok(())
    }

    #[test]
    fn real_ref_gc_bias_then_gc_bias_package_is_neutral_in_single_bin_case_for_lengths()
    -> Result<()> {
        // Human verification status: unverified
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
        cfg.set_windows(DistributionWindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
            gc_file: Some(gc_path),
            neutralize_invalid_gc: false,
        });
        cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
        cfg.set_gc_length_weighting(MarginalizeLengthsWeightingScheme::Equal);
        cfg.set_per_bp_length_bins(60, 60);

        // Manual expectations:
        // - `ref-gc-bias` is run for exactly one fragment length: 60 bp.
        // - On the simple reference ("ACGT" repeated), every 60 bp fragment has GC%=50,
        //   so the reference package has exactly one populated GC-by-length cell.
        // - `gc-bias` on `simple_inward_bam` also places all cfDNA mass in that same single cell.
        // - A 1x1 normalized cfDNA count divided by a 1x1 normalized reference count gives 1.0,
        //   so the produced correction package is neutral.
        // - `lengths` therefore receives GC weight 1.0 for the only fragment.
        // - This test also constrains the counted length range to exactly 60 bp, so the output has
        //   one row and one column, with the single cell equal to 1.0.
        run(&cfg)?;

        let npy_path = out_dir
            .path()
            .join(dot_join(&[cfg.output_prefix.trim(), "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.dim(), (1, 1));
        assert!((arr[(0, 0)] - 1.0).abs() < 1e-12);
        assert!((arr.sum() - 1.0).abs() < 1e-12);

        Ok(())
    }

    #[test]
    fn real_ref_gc_bias_then_gc_bias_package_changes_lengths_in_expected_direction() -> Result<()> {
        // Human verification status: unverified
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
            // Chromosome length 200 and fragment length 10 give 191 valid starts in total. Under
            // the `ref-gc-bias` fit rule these BED rows still contribute balanced pure-A and
            // pure-C support, so the expected downstream weights remain 5.0 and 5/9.
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
        cfg.set_windows(DistributionWindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
            gc_file: Some(gc_path),
            neutralize_invalid_gc: false,
        });
        cfg.set_ref_2bit(Some(reference.path.clone()));
        cfg.set_gc_length_weighting(MarginalizeLengthsWeightingScheme::Equal);
        cfg.set_per_bp_length_bins(10, 10);

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

    #[cfg(feature = "cmd_coverage_weights")]
    #[test]
    fn gc_file_and_scaling_tsv_weights_multiply_in_lengths() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Producer BAM:
        // - `simple_inward_bam()` contains one fragment [20, 80) on chr1.
        // - Run `coverage-weights` with a neutral real GC package over its full configured
        //   fragment length range so the written scaling TSV is GC-compatible with the consumer
        //   command, but the numerical scaling profile stays unchanged.
        // - With `bin_size = stride = 20`, the written scaling profile is therefore still the
        //   identity over the covered stride bins:
        //     [20,40): 1
        //     [40,60): 1
        //     [60,80): 1
        //     everything else: 0
        //
        // Consumer BAM:
        // - One 61 bp fragment [20, 81), so the only occupied length bin is 61.
        //
        // Scaling derivation in global `lengths` mode:
        // - Global mode has exactly one count-window, so the count weight is 1.0.
        // - For non-`count-overlap` assignment, `lengths` averages scaling over the full fragment:
        //     [20,40): 20 bp at factor 1
        //     [40,60): 20 bp at factor 1
        //     [60,80): 20 bp at factor 1
        //     [80,81):  1 bp at factor 0
        // - Average scaling over the fragment is therefore:
        //     (20 + 20 + 20 + 0) / 61 = 60 / 61.
        //
        // GC derivation:
        // - Use the smallest valid GC package for the only supported fragment length 61:
        //     length_edges = [61, 62]
        //     gc_edges     = [0, 101]
        //     correction_matrix = [[3.0]]
        // - Every accepted fragment therefore gets GC weight 3.0.
        //
        // Final contract:
        // - `lengths` multiplies scaling and GC correction for the fragment before incrementing the
        //   length bin.
        // - The only occupied cell must therefore be:
        //     3.0 * (60 / 61) = 180 / 61.
        let producer_bam = simple_inward_bam()?;
        let consumer_bam = bam_from_specs(
            vec![("chr1".to_string(), 200)],
            vec![fragment_on_tid(0, 20, 61, 20)],
            Vec::new(),
            "lengths_gc_and_scaling_consumer",
        )?;
        let ref_twobit = simple_reference_twobit()?;
        let out_dir = TempDir::new()?;
        let weights_out_dir = out_dir.path().join("coverage_weights");
        std::fs::create_dir_all(&weights_out_dir)?;
        let mut scaling_cfg =
            make_simple_coverage_weights_config(&weights_out_dir, &producer_bam.bam);
        let weights_gc_path = build_real_neutral_gc_package_for_range(
            &producer_bam.bam,
            &ref_twobit.path,
            out_dir.path(),
            10,
            200,
        )?;
        let gc_path = out_dir.path().join("constant_gc_pkg.npz");
        let package = GCCorrectionPackage {
            version: GC_CORRECTION_SCHEMA_VERSION,
            end_offset: 0,
            length_edges: vec![61, 62],
            gc_edges: vec![0, 101],
            length_bin_frequencies: array![1.0_f64],
            reference_contig_footprint: twobit_contig_footprint(&ref_twobit.path)?,
            correction_matrix: array![[3.0_f64]],
        };
        package.write_npz(&gc_path)?;
        scaling_cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgs {
            gc_file: Some(weights_gc_path),
            gc_tag: None,
            neutralize_invalid_gc: false,
        });
        scaling_cfg.set_ref_2bit(Some(ref_twobit.path.clone()));

        // Act
        run_coverage_weights(&scaling_cfg)?;
        let scaling_path = weights_out_dir.join("coverage.coverage.scaling_factors.tsv");

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: consumer_bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(DistributionWindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_scaling_factors(Some(scaling_path));
        cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
            gc_file: Some(gc_path),
            neutralize_invalid_gc: false,
        });
        cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
        cfg.set_per_bp_length_bins(61, 61);

        run(&cfg)?;

        // Assert
        let npy_path = out_dir
            .path()
            .join(dot_join(&[cfg.output_prefix.trim(), "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.shape(), &[1, 1]);

        let expected = 180.0_f64 / 61.0_f64;
        assert!(
            (arr[(0, 0)] - expected).abs() <= 1e-9,
            "expected weighted length count {expected}, got {}",
            arr[(0, 0)]
        );
        assert!(
            (arr.sum() - expected).abs() <= 1e-9,
            "expected total mass {expected}, got {}",
            arr.sum()
        );

        Ok(())
    }

    #[test]
    fn gc_file_rejects_package_when_fragment_length_range_is_outside_supported_range() -> Result<()>
    {
        // Human verification status: unverified
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
            reference_contig_footprint: twobit_contig_footprint(&ref_twobit.path)?,
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
        cfg.set_windows(DistributionWindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
            gc_file: Some(gc_path),
            neutralize_invalid_gc: false,
        });
        cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
        cfg.set_per_bp_length_bins(60, 60);

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
        // Human verification status: unverified
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
            reference_contig_footprint: twobit_contig_footprint(&ref_twobit.path)?,
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
        cfg.set_windows(DistributionWindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
            gc_file: Some(gc_path),
            neutralize_invalid_gc: false,
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
        // Human verification status: unverified
        let bam = simple_inward_bam()?;
        let gc_dir = TempDir::new()?;
        let gc_path = gc_dir.path().join("gc_pkg.npz");
        let ref_twobit = simple_reference_twobit()?;
        build_gc_package(&gc_path, 0, twobit_contig_footprint(&ref_twobit.path)?)?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: gc_dir.path().to_path_buf(),
                n_threads: 2,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(DistributionWindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
            gc_file: Some(gc_path.clone()),
            neutralize_invalid_gc: false,
        });
        cfg.set_per_bp_length_bins(10, 200);
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
        // Human verification status: unverified
        let bam = simple_inward_bam()?;
        let ref_twobit = simple_reference_twobit()?;
        let gc_dir = TempDir::new()?;
        let gc_path = gc_dir.path().join("gc_pkg.npz");
        // Choose large end_offset so offset_start >= offset_end, causing GC weight failure
        build_gc_package(&gc_path, 40, twobit_contig_footprint(&ref_twobit.path)?)?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: gc_dir.path().to_path_buf(),
                n_threads: 2,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(DistributionWindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
            gc_file: Some(gc_path.clone()),
            neutralize_invalid_gc: true,
        });
        cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
        cfg.set_per_bp_length_bins(10, 200);

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
        // Human verification status: unverified
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
        base_cfg.set_windows(DistributionWindowsArgs::default());
        base_cfg.set_window_assignment(AssignToWindowArgs::default());
        base_cfg.set_min_mapq(0);
        base_cfg.set_require_proper_pair(false);
        base_cfg.set_per_bp_length_bins(10, 100);

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
    fn max_deletion_bases_filters_indel_adjusted_fragments_before_counting() -> Result<()> {
        // Human verification status: unverified
        // Reuse the indel fixture with one clean fragment, one insertion-bearing fragment,
        // and one deletion-bearing fragment.
        //
        // Mental derivation:
        // - In adjust mode, the deletion-bearing fragment has one deleted reference base.
        // - max_deletion_bases=0 therefore drops only that fragment.
        // - The clean fragment remains in adjusted-length bin 24.
        // - The insertion-bearing fragment remains in adjusted-length bin 17.
        let bam = indel_bam_fixture()?;
        let out_dir = TempDir::new()?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 2,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_indel_mode(IndelMode::Adjust);
        cfg.max_deletion_bases = 0;
        cfg.set_windows(DistributionWindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 100);

        run(&cfg)?;

        let npy_path = out_dir
            .path()
            .join(dot_join(&[cfg.output_prefix.trim(), "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;
        let clean_length_bin = 24 - 10;
        let insertion_length_bin = 17 - 10;
        let deletion_length_bin = 10 - 10;
        assert!((arr[(0, clean_length_bin)] - 1.0).abs() < 1e-6);
        assert!((arr[(0, insertion_length_bin)] - 1.0).abs() < 1e-6);
        assert!((arr[(0, deletion_length_bin)]).abs() < 1e-6);
        assert!((arr.sum() - 2.0).abs() < 1e-6);

        Ok(())
    }

    #[test]
    fn clip_adjust_counts_adjusted_length_and_clip_skip_drops() -> Result<()> {
        // Human verification status: unverified
        // One unpaired read-as-fragment with cigar 2S10M2S at pos 10.
        //
        // Mental derivation:
        // - aligned fragment span is [10,20), so aligned mode counts length 10
        // - clip-adjust mode uses 10 + 2 + 2 = 14
        // - clip-skip mode rejects the fragment because it has soft clipping
        let bam = single_read_fragment_bam_with_cigar(
            "lengths_clip_modes_counting",
            10,
            vec![('S', 2), ('M', 10), ('S', 2)],
            b"TTAAAAAAAAAAAA".to_vec(),
        )?;

        let build_cfg = |out_dir: &std::path::Path, clip_mode: ClipMode| {
            let mut cfg = LengthsConfig::new(
                IOCArgs {
                    bam: bam.bam.clone(),
                    output_dir: out_dir.to_path_buf(),
                    n_threads: 1,
                },
                base_chromosomes(&["chr1"]),
            );
            cfg.set_indel_mode(IndelMode::Ignore);
            cfg.clip_mode = clip_mode;
            cfg.set_unpaired(UnpairedArgs {
                reads_are_fragments: true,
            });
            cfg.set_windows(DistributionWindowsArgs::default());
            cfg.set_window_assignment(AssignToWindowArgs::default());
            cfg.set_min_mapq(0);
            cfg.set_require_proper_pair(false);
            cfg.set_per_bp_length_bins(10, 14);
            cfg
        };

        let aligned_out = TempDir::new()?;
        let adjust_out = TempDir::new()?;
        let skip_out = TempDir::new()?;

        let aligned_cfg = build_cfg(aligned_out.path(), ClipMode::Aligned);
        let adjust_cfg = build_cfg(adjust_out.path(), ClipMode::Adjust);
        let skip_cfg = build_cfg(skip_out.path(), ClipMode::Skip);

        run(&aligned_cfg)?;
        run(&adjust_cfg)?;
        run(&skip_cfg)?;

        let aligned_path = aligned_out.path().join(dot_join(&[
            aligned_cfg.output_prefix.trim(),
            "length_counts.npy",
        ]));
        let adjust_path = adjust_out.path().join(dot_join(&[
            adjust_cfg.output_prefix.trim(),
            "length_counts.npy",
        ]));
        let skip_path = skip_out.path().join(dot_join(&[
            skip_cfg.output_prefix.trim(),
            "length_counts.npy",
        ]));

        let aligned_arr: Array2<f64> = read_npy(&aligned_path)?;
        let adjust_arr: Array2<f64> = read_npy(&adjust_path)?;
        let skip_arr: Array2<f64> = read_npy(&skip_path)?;

        assert_eq!(aligned_arr.shape(), &[1, 5]);
        assert_eq!(adjust_arr.shape(), &[1, 5]);
        assert_eq!(skip_arr.shape(), &[1, 5]);
        assert!((aligned_arr[(0, 0)] - 1.0).abs() < 1e-6);
        assert!((aligned_arr.sum() - 1.0).abs() < 1e-6);
        assert!((adjust_arr[(0, 4)] - 1.0).abs() < 1e-6);
        assert!((adjust_arr.sum() - 1.0).abs() < 1e-6);
        assert!((skip_arr.sum()).abs() < 1e-6);

        Ok(())
    }

    #[test]
    fn clip_adjust_count_overlap_uses_the_adjusted_assignment_interval() -> Result<()> {
        // Human verification status: unverified
        // One unpaired 2S10M2S fragment at pos 10.
        //
        // Aligned mode:
        // - aligned interval [10,20) sits fully inside bin [10,20)
        // - length bin is 10
        //
        // Adjust mode:
        // - assignment interval expands to [8,22), length 14
        // - overlap with 10 bp bins is:
        //   [0,10):  2 / 14
        //   [10,20): 10 / 14
        //   [20,30): 2 / 14
        let bam = single_read_fragment_bam_with_cigar(
            "lengths_clip_adjust_overlap",
            10,
            vec![('S', 2), ('M', 10), ('S', 2)],
            b"TTAAAAAAAAAAAA".to_vec(),
        )?;

        let build_cfg = |out_dir: &std::path::Path, clip_mode: ClipMode| {
            let mut cfg = LengthsConfig::new(
                IOCArgs {
                    bam: bam.bam.clone(),
                    output_dir: out_dir.to_path_buf(),
                    n_threads: 1,
                },
                base_chromosomes(&["chr1"]),
            );
            cfg.set_indel_mode(IndelMode::Ignore);
            cfg.clip_mode = clip_mode;
            cfg.set_unpaired(UnpairedArgs {
                reads_are_fragments: true,
            });
            cfg.set_windows(DistributionWindowsArgs {
                by_size: Some(10),
                by_bed: None,
                by_grouped_bed: None,
            });
            cfg.set_window_assignment(AssignToWindowArgs {
                assign_by: WindowAssigner::CountOverlap,
            });
            cfg.set_min_mapq(0);
            cfg.set_require_proper_pair(false);
            cfg.set_tile_size(10);
            cfg.set_per_bp_length_bins(10, 14);
            cfg
        };

        let aligned_out = TempDir::new()?;
        let adjust_out = TempDir::new()?;

        let aligned_cfg = build_cfg(aligned_out.path(), ClipMode::Aligned);
        let adjust_cfg = build_cfg(adjust_out.path(), ClipMode::Adjust);

        run(&aligned_cfg)?;
        run(&adjust_cfg)?;

        let aligned_path = aligned_out.path().join(dot_join(&[
            aligned_cfg.output_prefix.trim(),
            "length_counts.npy",
        ]));
        let adjust_path = adjust_out.path().join(dot_join(&[
            adjust_cfg.output_prefix.trim(),
            "length_counts.npy",
        ]));

        let aligned_arr: Array2<f64> = read_npy(&aligned_path)?;
        let adjust_arr: Array2<f64> = read_npy(&adjust_path)?;

        assert!((aligned_arr[(1, 0)] - 1.0).abs() < 1e-6);
        assert!((aligned_arr.row(0).sum()).abs() < 1e-6);
        assert!((aligned_arr.row(2).sum()).abs() < 1e-6);
        assert!((aligned_arr.sum() - 1.0).abs() < 1e-6);

        assert!((adjust_arr[(0, 4)] - (2.0 / 14.0)).abs() < 1e-6);
        assert!((adjust_arr[(1, 4)] - (10.0 / 14.0)).abs() < 1e-6);
        assert!((adjust_arr[(2, 4)] - (2.0 / 14.0)).abs() < 1e-6);
        assert!((adjust_arr.sum() - 1.0).abs() < 1e-6);

        Ok(())
    }

    #[test]
    fn clip_adjust_bins_by_adjusted_length_but_scales_over_aligned_span() -> Result<()> {
        // Human verification status: unverified
        // One unpaired 2S10M2S fragment at pos 10.
        //
        // Mental derivation:
        // - clip-adjust mode bins this fragment at adjusted length 14
        // - scaling must still use the aligned span [10,20)
        // - with factors [0,10):1, [10,20):3, [20,200):1, the aligned-span average is exactly 3
        // - a buggy assignment-interval average over [8,22) would instead be 34 / 14
        let bam = single_read_fragment_bam_with_cigar(
            "lengths_clip_adjust_scaling",
            10,
            vec![('S', 2), ('M', 10), ('S', 2)],
            b"TTAAAAAAAAAAAA".to_vec(),
        )?;
        let out_dir = TempDir::new()?;
        let scaling_path = out_dir.path().join("clip_adjust_scaling.tsv");
        write_scaling_factors(
            &scaling_path,
            &[
                ("chr1", 0, 10, 1.0_f32),
                ("chr1", 10, 20, 3.0_f32),
                ("chr1", 20, 200, 1.0_f32),
            ],
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
        cfg.clip_mode = ClipMode::Adjust;
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.set_windows(DistributionWindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_scaling_factors(Some(scaling_path));
        cfg.set_per_bp_length_bins(14, 14);

        run(&cfg)?;

        let npy_path = out_dir
            .path()
            .join(dot_join(&[cfg.output_prefix.trim(), "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.shape(), &[1, 1]);
        assert!((arr[(0, 0)] - 3.0).abs() < 1e-6);
        assert!((arr.sum() - 3.0).abs() < 1e-6);

        Ok(())
    }

    #[test]
    fn clip_adjust_count_overlap_scaling_uses_the_nearest_aligned_base_for_clipped_only_windows()
    -> Result<()> {
        // Human verification status: unverified
        // One unpaired 2S10M2S fragment at pos 10, counted against three 10 bp windows.
        //
        // Assignment in clip-adjust mode uses [8,22), so the count fractions are:
        // - [0,10):  2 / 14
        // - [10,20): 10 / 14
        // - [20,30): 2 / 14
        //
        // Scaling should remain reference-based:
        // - the middle window samples its aligned overlap [10,20) => weight 3
        // - the clipped-only left/right windows have no aligned overlap, so they should use the
        //   nearest aligned base instead of the flanking bins
        // - with flanking bins set to 11 and 13, all three windows should still use weight 3
        let bam = single_read_fragment_bam_with_cigar(
            "lengths_clip_adjust_scaling_nearest_base",
            10,
            vec![('S', 2), ('M', 10), ('S', 2)],
            b"TTAAAAAAAAAAAA".to_vec(),
        )?;
        let out_dir = TempDir::new()?;
        let bed_path = out_dir.path().join("clip_adjust_scaling_windows.bed");
        let scaling_path = out_dir.path().join("clip_adjust_scaling.tsv");
        write_bed(
            &bed_path,
            &[
                ("chr1", 0, 10, "left"),
                ("chr1", 10, 20, "middle"),
                ("chr1", 20, 30, "right"),
            ],
        )?;
        write_scaling_factors(
            &scaling_path,
            &[
                ("chr1", 0, 10, 11.0_f32),
                ("chr1", 10, 20, 3.0_f32),
                ("chr1", 20, 200, 13.0_f32),
            ],
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
        cfg.clip_mode = ClipMode::Adjust;
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(bed_path),
            by_grouped_bed: None,
        });
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_scaling_factors(Some(scaling_path));
        cfg.set_per_bp_length_bins(14, 14);

        run(&cfg)?;

        let npy_path = out_dir
            .path()
            .join(dot_join(&[cfg.output_prefix.trim(), "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.shape(), &[3, 1]);
        assert!((arr[(0, 0)] - (2.0 / 14.0) * 3.0).abs() < 1e-6);
        assert!((arr[(1, 0)] - (10.0 / 14.0) * 3.0).abs() < 1e-6);
        assert!((arr[(2, 0)] - (2.0 / 14.0) * 3.0).abs() < 1e-6);
        assert!((arr.sum() - 3.0).abs() < 1e-6);

        Ok(())
    }

    #[test]
    fn max_soft_clips_filters_lengths_fragments_before_counting_clip_adjusted_lengths() -> Result<()>
    {
        // Human verification status: unverified
        // One unpaired 2S10M fragment at pos 10.
        //
        // Mental derivation:
        // - adjusted length is 12
        // - max_soft_clips=2 keeps it because the left end equals the threshold
        // - max_soft_clips=1 drops it before counting
        let bam = single_read_fragment_bam_with_cigar(
            "lengths_max_soft_clips",
            10,
            vec![('S', 2), ('M', 10)],
            b"TTAAAAAAAAAT".to_vec(),
        )?;

        let build_cfg = |out_dir: &std::path::Path, max_soft_clips: u16| {
            let mut cfg = LengthsConfig::new(
                IOCArgs {
                    bam: bam.bam.clone(),
                    output_dir: out_dir.to_path_buf(),
                    n_threads: 1,
                },
                base_chromosomes(&["chr1"]),
            );
            cfg.set_indel_mode(IndelMode::Ignore);
            cfg.clip_mode = ClipMode::Adjust;
            cfg.max_soft_clips = max_soft_clips;
            cfg.set_unpaired(UnpairedArgs {
                reads_are_fragments: true,
            });
            cfg.set_windows(DistributionWindowsArgs::default());
            cfg.set_window_assignment(AssignToWindowArgs::default());
            cfg.set_min_mapq(0);
            cfg.set_require_proper_pair(false);
            cfg.set_per_bp_length_bins(12, 12);
            cfg
        };

        let keep_out = TempDir::new()?;
        let drop_out = TempDir::new()?;

        let keep_cfg = build_cfg(keep_out.path(), 2);
        let drop_cfg = build_cfg(drop_out.path(), 1);

        run(&keep_cfg)?;
        run(&drop_cfg)?;

        let keep_path = keep_out.path().join(dot_join(&[
            keep_cfg.output_prefix.trim(),
            "length_counts.npy",
        ]));
        let drop_path = drop_out.path().join(dot_join(&[
            drop_cfg.output_prefix.trim(),
            "length_counts.npy",
        ]));

        let keep_arr: Array2<f64> = read_npy(&keep_path)?;
        let drop_arr: Array2<f64> = read_npy(&drop_path)?;

        assert_eq!(keep_arr.shape(), &[1, 1]);
        assert_eq!(drop_arr.shape(), &[1, 1]);
        assert!((keep_arr[(0, 0)] - 1.0).abs() < 1e-6);
        assert!((keep_arr.sum() - 1.0).abs() < 1e-6);
        assert!((drop_arr.sum()).abs() < 1e-6);

        Ok(())
    }

    #[test]
    fn indel_adjust_bins_by_adjusted_length_but_scales_over_reference_span() -> Result<()> {
        // Human verification status: unverified
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
        cfg.set_windows(DistributionWindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_scaling_factors(Some(scaling_path));
        // Keep only the insertion-bearing fragment from the edge fixture.
        cfg.set_per_bp_length_bins(17, 17);

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
        // Human verification status: unverified
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
        cfg.set_windows(DistributionWindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.blacklist = Some(vec![blacklist_path]);
        // Keep only the deletion-bearing fragment from the edge fixture.
        cfg.set_per_bp_length_bins(10, 10);

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
        // Human verification status: unverified
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
        base_cfg.set_windows(DistributionWindowsArgs::default());
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
            cfg.set_per_bp_length_bins(10, 11);

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
        // Human verification status: unverified
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
        cfg.set_per_bp_length_bins(10, 200);

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
        // Human verification status: unverified
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
        cfg.set_windows(DistributionWindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_per_bp_length_bins(10, 200);

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
        // Human verification status: unverified
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
        cfg.set_windows(DistributionWindowsArgs {
            by_size: Some(200),
            by_bed: None,
            by_grouped_bed: None,
        });
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_per_bp_length_bins(10, 100);

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
        // Human verification status: unverified
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
            cfg.set_windows(DistributionWindowsArgs {
                by_size: Some(window_bp),
                by_bed: None,
                by_grouped_bed: None,
            });
            cfg.set_window_assignment(AssignToWindowArgs { assign_by });
            cfg.set_min_mapq(0);
            cfg.set_require_proper_pair(false);
            cfg.set_per_bp_length_bins(10, 200);

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
    fn midpoint_assignment_on_even_length_boundary_counts_exactly_one_adjacent_window() -> Result<()>
    {
        // Arrange:
        // One even-length fragment spans [40,50), so midpoint assignment randomizes between:
        //   44 and 45
        //
        // Put the window boundary exactly between those two central bases:
        //   window 0 -> [0,45)
        //   window 1 -> [45,90)
        //
        // With `assign_by=midpoint`, exactly one of those windows must receive the fragment. This
        // locks in the released contract near the midpoint seam without pretending the random tie
        // chooses one side deterministically.
        let bam = bam_from_specs(
            vec![("chr1".to_string(), 100)],
            vec![fixtures::paired_fragment(40, 10, 5)],
            Vec::new(),
            "lengths_even_midpoint_boundary",
        )?;
        let out_dir = TempDir::new()?;
        let bed_path = out_dir.path().join("windows.bed");
        write_bed(
            &bed_path,
            &[("chr1", 0, 45, "window0"), ("chr1", 45, 90, "window1")],
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
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(bed_path),
            by_grouped_bed: None,
        });
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::Midpoint,
        });
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 10);

        // Act
        run(&cfg)?;

        // Assert
        let npy_path = out_dir
            .path()
            .join(dot_join(&[cfg.output_prefix.trim(), "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.shape(), &[2, 1]);
        assert!((arr.sum() - 1.0).abs() < 1e-6);
        let first_window = arr[(0, 0)];
        let second_window = arr[(1, 0)];
        let is_valid_one_hot = (first_window - 1.0).abs() < 1e-6 && second_window.abs() < 1e-6
            || first_window.abs() < 1e-6 && (second_window - 1.0).abs() < 1e-6;
        assert!(
            is_valid_one_hot,
            "midpoint tie at the window edge must count exactly one adjacent window, got [{first_window}, {second_window}]"
        );

        Ok(())
    }

    #[test]
    fn scaling_tsv_must_cover_requested_chromosome_end_in_lengths() -> Result<()> {
        // Human verification status: unverified
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
        cfg.set_windows(DistributionWindowsArgs::default());
        cfg.set_window_assignment(AssignToWindowArgs::default());
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_scaling_factors(Some(scaling_path));
        cfg.set_per_bp_length_bins(10, 200);

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
        // Human verification status: unverified
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
        let scaling_path = weights_out_dir.join("coverage.coverage.scaling_factors.tsv");

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: consumer_bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(bed_path),
            by_grouped_bed: None,
        });
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::CountOverlap,
        });
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_scaling_factors(Some(scaling_path));
        cfg.set_per_bp_length_bins(61, 61);

        run(&cfg)?;

        // Assert
        let npy_path = out_dir
            .path()
            .join(dot_join(&[cfg.output_prefix.trim(), "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&npy_path)?;
        assert_eq!(arr.shape(), &[1, 1]);

        // `count-overlap` keeps overlap fractions as `f64`, so the stable contract here is the
        // direct `f64` arithmetic for 11 / 61.
        let expected = 11.0_f64 / 61.0_f64;
        assert!(
            (arr[(0, 0)] - expected).abs() <= 1e-12,
            "expected weighted count {expected}, got {}",
            arr[(0, 0)]
        );
        assert!(
            (arr.sum() - expected).abs() <= 1e-12,
            "expected total mass {expected}, got {}",
            arr.sum()
        );
        Ok(())
    }

    #[test]
    fn bed_windowed_runs_write_prefixed_bins_tsv_with_exact_blacklisted_fractions() -> Result<()> {
        // Arrange:
        // - `simple_inward_bam()` gives one 60 bp fragment on chr1 spanning [20,80).
        // - The two BED windows are [10,20) and [20,30).
        // - The blacklist interval [15,20) overlaps only the first window for 5 of its 10 bases.
        // - Therefore the persisted window metadata must be:
        //     chr1  10  20  0.5
        //     chr1  20  30  0
        let bam = simple_inward_bam()?;
        let out_dir = TempDir::new()?;
        let windows_bed = out_dir.path().join("windows.bed");
        let blacklist_bed = out_dir.path().join("blacklist.bed");
        write_bed(
            &windows_bed,
            &[("chr1", 10, 20, "left"), ("chr1", 20, 30, "right")],
        )?;
        write_bed(&blacklist_bed, &[("chr1", 15, 20, "masked")])?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.output_prefix = "sampleA".to_string();
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(windows_bed),
            by_grouped_bed: None,
        });
        cfg.blacklist = Some(vec![blacklist_bed]);
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(60, 60);

        // Act
        run(&cfg)?;
        let bins_tsv =
            std::fs::read_to_string(out_dir.path().join(dot_join(&["sampleA", "bins.tsv"])))?;

        // Assert
        assert_eq!(
            bins_tsv,
            concat!(
                "chrom\tstart\tend\tblacklisted_fraction\n",
                "chr1\t10\t20\t0.5\n",
                "chr1\t20\t30\t0\n"
            )
        );
        Ok(())
    }

    #[test]
    fn grouped_bed_count_overlap_uses_first_group_occurrence_keeps_zero_rows_and_writes_group_metadata()
    -> Result<()> {
        // Arrange:
        // - one unpaired fragment spans [10,20) with exact fragment length 10
        // - grouped intervals use mixed widths on purpose:
        //     [10,15) beta  -> overlap 5/10 = 0.5
        //     [40,50) beta  -> contributes no counts, but 5/10 of the group's bp are blacklisted
        //     [15,20) alpha -> overlap 5/10 = 0.5
        //     [30,35) gamma -> zero-count group, fully blacklisted in metadata
        // - first group occurrence order in the BED scan is beta, alpha, gamma
        // - grouped `count-overlap` therefore yields:
        //     beta  -> 0.5 in the lone length-10 column
        //     alpha -> 0.5
        //     gamma -> 0.0
        // - group-level blacklist fractions are width-weighted across loaded intervals:
        //     beta  -> 5 / (5 + 10) = 1/3
        //     alpha -> 0
        //     gamma -> 1
        let bam = single_read_fragment_bam("lengths_grouped_count_overlap", 10, 10)?;
        let out_dir = TempDir::new()?;
        let grouped_bed = out_dir.path().join("grouped_windows.bed");
        let blacklist_bed = out_dir.path().join("blacklist.bed");
        write_bed(
            &grouped_bed,
            &[
                ("chr1", 10, 15, "beta"),
                ("chr1", 40, 50, "beta"),
                ("chr1", 15, 20, "alpha"),
                ("chr1", 30, 35, "gamma"),
            ],
        )?;
        write_bed(
            &blacklist_bed,
            &[
                ("chr1", 45, 50, "masked_beta"),
                ("chr1", 30, 35, "masked_gamma"),
            ],
        )?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.output_prefix = "sampleA".to_string();
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        });
        cfg.blacklist = Some(vec![blacklist_bed]);
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::CountOverlap,
        });
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 10);

        // Act
        run(&cfg)?;
        let counts_path = out_dir
            .path()
            .join(dot_join(&["sampleA", "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&counts_path)?;
        let group_index = std::fs::read_to_string(
            out_dir
                .path()
                .join(dot_join(&["sampleA", "group_index.tsv"])),
        )?;
        let settings_text = std::fs::read_to_string(
            out_dir
                .path()
                .join(dot_join(&["sampleA", "fragment_length_settings.json"])),
        )?;
        let settings: Value = serde_json::from_str(&settings_text)?;

        // Assert
        assert_eq!(arr.shape(), &[3, 1]);
        assert!((arr[(0, 0)] - 0.5).abs() < 1e-12);
        assert!((arr[(1, 0)] - 0.5).abs() < 1e-12);
        assert_eq!(arr[(2, 0)], 0.0);
        assert!((arr.sum() - 1.0).abs() < 1e-12);

        assert_eq!(
            parse_group_index_rows(&group_index),
            vec![
                (0, "beta".to_string(), 1.0 / 3.0),
                (1, "alpha".to_string(), 0.0),
                (2, "gamma".to_string(), 1.0),
            ]
        );
        assert_eq!(settings["aggregation_level"], "groups");
        assert_eq!(settings["window_mode"], "by-grouped-bed");
        assert!(!out_dir.path().join("sampleA.bins.tsv").exists());
        assert!(!out_dir.path().join("sampleA.grouped_windows.tsv").exists());
        Ok(())
    }

    #[test]
    fn grouped_bed_count_overlap_sums_same_group_window_weights_above_one() -> Result<()> {
        // Arrange:
        // - one unpaired fragment spans [10,20) with exact fragment length 10
        // - grouped intervals:
        //     [10,15) beta  -> overlap 5/10 = 0.5
        //     [12,18) beta  -> overlap 6/10 = 0.6
        //     [15,20) alpha -> overlap 5/10 = 0.5
        //     [30,35) gamma -> zero row
        // - grouped `count-overlap` must therefore yield:
        //     beta  -> 1.1 in the lone length-10 column
        //     alpha -> 0.5
        //     gamma -> 0.0
        // This fails if grouped mode unions same-group windows or normalizes their combined
        // overlap mass back down to one fragment.
        let bam = single_read_fragment_bam("lengths_grouped_count_overlap_above_one", 10, 10)?;
        let out_dir = TempDir::new()?;
        let grouped_bed = out_dir.path().join("grouped_windows_overlap_mass.bed");
        write_bed(
            &grouped_bed,
            &[
                ("chr1", 10, 15, "beta"),
                ("chr1", 12, 18, "beta"),
                ("chr1", 15, 20, "alpha"),
                ("chr1", 30, 35, "gamma"),
            ],
        )?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        });
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::CountOverlap,
        });
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 10);

        // Act
        run(&cfg)?;
        let counts_path = out_dir.path().join("length_counts.npy");
        let arr: Array2<f64> = read_npy(&counts_path)?;
        let group_index = std::fs::read_to_string(out_dir.path().join("group_index.tsv"))?;

        // Assert
        assert_eq!(
            parse_group_index_tsv(&group_index),
            vec![
                (0, "beta".to_string()),
                (1, "alpha".to_string()),
                (2, "gamma".to_string()),
            ]
        );
        assert_eq!(arr.shape(), &[3, 1]);
        assert!((arr[(0, 0)] - 1.1).abs() < 1e-12);
        assert!((arr[(1, 0)] - 0.5).abs() < 1e-12);
        assert_eq!(arr[(2, 0)], 0.0);
        assert!((arr.sum() - 1.6).abs() < 1e-12);
        Ok(())
    }

    #[test]
    fn grouped_bed_group_index_omits_blacklisted_fraction_without_blacklist() -> Result<()> {
        // Arrange:
        // - grouped output without any blacklist input should describe only the group index mapping
        // - the metadata file should therefore have exactly two columns:
        //     group_idx, group_name
        // - it should not include a synthetic zero-filled blacklist column
        let bam = single_read_fragment_bam("lengths_grouped_group_index_no_blacklist", 10, 10)?;
        let out_dir = TempDir::new()?;
        let grouped_bed = out_dir.path().join("grouped_windows_no_blacklist.bed");
        write_bed(
            &grouped_bed,
            &[("chr1", 10, 20, "beta"), ("chr1", 30, 40, "gamma")],
        )?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        });
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 10);

        // Act
        run(&cfg)?;
        let group_index = std::fs::read_to_string(out_dir.path().join("group_index.tsv"))?;

        // Assert
        let rows: Vec<&str> = group_index.lines().collect();
        assert_eq!(rows[0], "group_idx\tgroup_name");
        assert_eq!(rows[1], "0\tbeta");
        assert_eq!(rows[2], "1\tgamma");
        assert_eq!(rows.len(), 3);
        assert_eq!(
            parse_group_index_tsv(&group_index),
            vec![(0, "beta".to_string()), (1, "gamma".to_string())]
        );
        Ok(())
    }

    #[test]
    fn grouped_bed_any_counts_same_group_intervals_separately() -> Result<()> {
        // Arrange:
        // - one unpaired fragment spans [10,20) with length 10
        // - under `assign-by=any`, every overlapping interval gets the full fragment count
        // - grouped intervals:
        //     [10,15) beta
        //     [12,18) beta
        //     [15,20) alpha
        //     [30,35) gamma
        // - expected grouped counts:
        //     beta  -> 2.0
        //     alpha -> 1.0
        //     gamma -> 0.0
        // This fails if grouped mode unions same-group intervals or divides by the number of
        // windows in a group.
        let bam = single_read_fragment_bam("lengths_grouped_any", 10, 10)?;
        let out_dir = TempDir::new()?;
        let grouped_bed = out_dir.path().join("grouped_windows_any.bed");
        write_bed(
            &grouped_bed,
            &[
                ("chr1", 10, 15, "beta"),
                ("chr1", 12, 18, "beta"),
                ("chr1", 15, 20, "alpha"),
                ("chr1", 30, 35, "gamma"),
            ],
        )?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        });
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::Any,
        });
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 10);

        // Act
        run(&cfg)?;
        let counts_path = out_dir.path().join("length_counts.npy");
        let arr: Array2<f64> = read_npy(&counts_path)?;
        let group_index = std::fs::read_to_string(out_dir.path().join("group_index.tsv"))?;

        // Assert
        assert_eq!(
            parse_group_index_tsv(&group_index),
            vec![
                (0, "beta".to_string()),
                (1, "alpha".to_string()),
                (2, "gamma".to_string()),
            ]
        );
        assert_eq!(arr.shape(), &[3, 1]);
        assert_eq!(arr[(0, 0)], 2.0);
        assert_eq!(arr[(1, 0)], 1.0);
        assert_eq!(arr[(2, 0)], 0.0);
        assert_eq!(arr.sum(), 3.0);
        Ok(())
    }

    #[test]
    fn grouped_bed_any_aggregates_shared_groups_across_chromosomes_and_uses_prefixed_sidecars()
    -> Result<()> {
        // Arrange:
        // The three-chrom fixture contributes exactly one fragment per chromosome:
        // - chr1 -> length 60
        // - chr2 -> length 80
        // - chr3 -> length 100
        //
        // Grouped windows intentionally place `beta` first, reuse it on chr3, and keep `gamma`
        // empty, so row order and exact grouped counts must be:
        // - row 0 / beta  -> lengths 80 and 100
        // - row 1 / alpha -> length 60
        // - row 2 / gamma -> all zeros
        let bam = three_chrom_length_fixture("lengths_grouped_three_chr")?;
        let out_dir = TempDir::new()?;
        let grouped_bed = out_dir.path().join("grouped_three_chr.bed");
        write_bed(
            &grouped_bed,
            &[
                ("chr2", 0, 200, "beta"),
                ("chr1", 0, 200, "alpha"),
                ("chr3", 0, 200, "beta"),
                ("chr1", 150, 160, "gamma"),
            ],
        )?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1", "chr2", "chr3"]),
        );
        cfg.output_prefix = "sampleA".to_string();
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        });
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::Any,
        });
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 120);

        // Act
        run(&cfg)?;
        let counts_path = out_dir
            .path()
            .join(dot_join(&["sampleA", "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&counts_path)?;
        let group_index = std::fs::read_to_string(
            out_dir
                .path()
                .join(dot_join(&["sampleA", "group_index.tsv"])),
        )?;
        // Assert
        assert_eq!(
            parse_group_index_tsv(&group_index),
            vec![
                (0, "beta".to_string()),
                (1, "alpha".to_string()),
                (2, "gamma".to_string()),
            ]
        );
        assert_eq!(arr.shape(), &[3, 111]);
        assert_eq!(arr[(0, 80 - 10)], 1.0);
        assert_eq!(arr[(0, 100 - 10)], 1.0);
        assert_eq!(arr.row(0).sum(), 2.0);
        assert_eq!(arr[(1, 60 - 10)], 1.0);
        assert_eq!(arr.row(1).sum(), 1.0);
        assert_eq!(arr.row(2).sum(), 0.0);
        assert_eq!(arr.sum(), 3.0);

        assert!(out_dir.path().join("sampleA.group_index.tsv").exists());
        assert!(!out_dir.path().join("sampleA.grouped_windows.tsv").exists());
        assert!(!out_dir.path().join("group_index.tsv").exists());
        assert!(!out_dir.path().join("grouped_windows.tsv").exists());
        assert!(!out_dir.path().join("sampleA.bins.tsv").exists());
        Ok(())
    }

    #[test]
    fn grouped_bed_scaling_factors_aggregate_shared_groups_across_chromosomes() -> Result<()> {
        // Arrange:
        // The three-chrom fixture contributes exactly one fragment per chromosome:
        // - chr1 -> length 60
        // - chr2 -> length 80
        // - chr3 -> length 100
        //
        // Grouped windows intentionally reuse `beta` across chr2 and chr3 and place `alpha`
        // second in the BED. Per-chromosome scaling factors are:
        // - chr1 -> 1.5
        // - chr2 -> 2.0
        // - chr3 -> 3.0
        //
        // So grouped output must be:
        // - row 0 / beta  -> length 80 = 2.0, length 100 = 3.0
        // - row 1 / alpha -> length 60 = 1.5
        // - row 2 / gamma -> all zeros
        let bam = three_chrom_length_fixture("lengths_grouped_three_chr_scaling")?;
        let out_dir = TempDir::new()?;
        let grouped_bed = out_dir.path().join("grouped_three_chr_scaling.bed");
        let scaling_path = out_dir.path().join("grouped_three_chr_scaling.tsv");
        write_bed(
            &grouped_bed,
            &[
                ("chr2", 0, 200, "beta"),
                ("chr1", 0, 200, "alpha"),
                ("chr3", 0, 200, "beta"),
                ("chr1", 150, 160, "gamma"),
            ],
        )?;
        write_scaling_factors(
            &scaling_path,
            &[
                ("chr1", 0, 200, 1.5),
                ("chr2", 0, 200, 2.0),
                ("chr3", 0, 200, 3.0),
            ],
        )?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1", "chr2", "chr3"]),
        );
        cfg.output_prefix = "sampleA".to_string();
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_tile_size(50);
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        });
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::Any,
        });
        cfg.set_scaling_factors(Some(scaling_path));
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 120);

        // Act
        run(&cfg)?;
        let counts_path = out_dir
            .path()
            .join(dot_join(&["sampleA", "length_counts.npy"]));
        let arr: Array2<f64> = read_npy(&counts_path)?;
        let group_index = std::fs::read_to_string(
            out_dir
                .path()
                .join(dot_join(&["sampleA", "group_index.tsv"])),
        )?;

        // Assert
        assert_eq!(
            parse_group_index_tsv(&group_index),
            vec![
                (0, "beta".to_string()),
                (1, "alpha".to_string()),
                (2, "gamma".to_string()),
            ]
        );
        assert_eq!(arr.shape(), &[3, 111]);
        assert_eq!(arr[(0, 80 - 10)], 2.0);
        assert_eq!(arr[(0, 100 - 10)], 3.0);
        assert_eq!(arr.row(0).sum(), 5.0);
        assert_eq!(arr[(1, 60 - 10)], 1.5);
        assert_eq!(arr.row(1).sum(), 1.5);
        assert_eq!(arr.row(2).sum(), 0.0);
        assert_eq!(arr.sum(), 6.5);
        Ok(())
    }

    #[test]
    fn grouped_bed_any_aggregates_shared_groups_across_tiles_on_same_chromosome() -> Result<()> {
        // Arrange:
        // - two unpaired fragments both have length 10 and start in different tiles when
        //   `tile_size=50`:
        //     [10,20) in tile [0,50)
        //     [60,70) in tile [50,100)
        // - grouped intervals place the same group on both tiles:
        //     [10,20) beta
        //     [60,70) beta
        //     [120,130) gamma
        // - under `assign-by=any`, each fragment contributes exactly 1.0 to its overlapping
        //   window, so grouped output must aggregate to:
        //     beta  -> 2.0 in the lone length-10 column
        //     gamma -> 0.0
        let bam = bam_from_specs(
            vec![("chr1".to_string(), 200)],
            Vec::new(),
            vec![
                ReadSpec {
                    tid: 0,
                    pos: 10,
                    cigar: vec![('M', 10)],
                    seq: vec![b'A'; 10],
                    qual: 40,
                    is_reverse: false,
                    mapq: 60,
                    flags: 0,
                    mate_tid: None,
                    mate_pos: None,
                    insert_size: 0,
                },
                ReadSpec {
                    tid: 0,
                    pos: 60,
                    cigar: vec![('M', 10)],
                    seq: vec![b'A'; 10],
                    qual: 40,
                    is_reverse: false,
                    mapq: 60,
                    flags: 0,
                    mate_tid: None,
                    mate_pos: None,
                    insert_size: 0,
                },
            ],
            "lengths_grouped_same_chr_tiles",
        )?;
        let out_dir = TempDir::new()?;
        let grouped_bed = out_dir.path().join("grouped_same_chr_tiles.bed");
        write_bed(
            &grouped_bed,
            &[
                ("chr1", 10, 20, "beta"),
                ("chr1", 60, 70, "beta"),
                ("chr1", 120, 130, "gamma"),
            ],
        )?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_tile_size(50);
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        });
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::Any,
        });
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 10);

        // Act
        run(&cfg)?;
        let counts_path = out_dir.path().join("length_counts.npy");
        let arr: Array2<f64> = read_npy(&counts_path)?;
        let group_index = std::fs::read_to_string(out_dir.path().join("group_index.tsv"))?;

        // Assert
        assert_eq!(
            parse_group_index_tsv(&group_index),
            vec![(0, "beta".to_string()), (1, "gamma".to_string())]
        );
        assert_eq!(arr.shape(), &[2, 1]);
        assert_eq!(arr[(0, 0)], 2.0);
        assert_eq!(arr[(1, 0)], 0.0);
        assert_eq!(arr.sum(), 2.0);
        Ok(())
    }

    #[test]
    fn grouped_bed_scaling_factors_weight_each_grouped_count() -> Result<()> {
        // Arrange:
        // - one unpaired fragment spans [10,20) with length 10
        // - grouped windows keep that fragment in `beta` and one explicit zero row in `gamma`
        // - a chromosome-wide scaling factor of 2.0 should double the grouped count mass
        let bam = single_read_fragment_bam("lengths_grouped_scaling", 10, 10)?;
        let out_dir = TempDir::new()?;
        let grouped_bed = out_dir.path().join("grouped_scaling.bed");
        let scaling_path = out_dir.path().join("grouped_scaling.tsv");
        write_bed(
            &grouped_bed,
            &[("chr1", 10, 20, "beta"), ("chr1", 30, 40, "gamma")],
        )?;
        write_scaling_factors(&scaling_path, &[("chr1", 0, 200, 2.0)])?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        });
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::Any,
        });
        cfg.set_scaling_factors(Some(scaling_path));
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 10);

        // Act
        run(&cfg)?;
        let counts_path = out_dir.path().join("length_counts.npy");
        let arr: Array2<f64> = read_npy(&counts_path)?;
        let group_index = std::fs::read_to_string(out_dir.path().join("group_index.tsv"))?;

        // Assert
        assert_eq!(
            parse_group_index_tsv(&group_index),
            vec![(0, "beta".to_string()), (1, "gamma".to_string())]
        );
        assert_eq!(arr.shape(), &[2, 1]);
        assert_eq!(arr[(0, 0)], 2.0);
        assert_eq!(arr[(1, 0)], 0.0);
        assert_eq!(arr.sum(), 2.0);
        Ok(())
    }

    #[test]
    fn grouped_bed_blacklist_filtering_drops_matching_fragments_before_grouping() -> Result<()> {
        // Arrange:
        // - one unpaired fragment spans [10,20) with length 10
        // - grouped windows keep `beta` as the only would-be hit and preserve `gamma` as a zero row
        // - blacklisting [15,16) with `blacklist-strategy=any` must drop the fragment entirely
        //   before any grouped counting happens
        // - metadata still reflects the grouped BED geometry:
        //     beta  -> 1 / 10 = 0.1 blacklisted
        //     gamma -> 0
        let bam = single_read_fragment_bam("lengths_grouped_blacklist", 10, 10)?;
        let out_dir = TempDir::new()?;
        let grouped_bed = out_dir.path().join("grouped_blacklist.bed");
        let blacklist_bed = out_dir.path().join("grouped_blacklist_mask.bed");
        write_bed(
            &grouped_bed,
            &[("chr1", 10, 20, "beta"), ("chr1", 30, 40, "gamma")],
        )?;
        write_bed(&blacklist_bed, &[("chr1", 15, 16, "masked_beta")])?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        });
        cfg.blacklist = Some(vec![blacklist_bed]);
        cfg.blacklist_strategy = BlacklistStrategy::Any;
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::Any,
        });
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 10);

        // Act
        run(&cfg)?;
        let counts_path = out_dir.path().join("length_counts.npy");
        let arr: Array2<f64> = read_npy(&counts_path)?;
        let group_index = std::fs::read_to_string(out_dir.path().join("group_index.tsv"))?;

        // Assert
        assert_eq!(
            parse_group_index_rows(&group_index),
            vec![(0, "beta".to_string(), 0.1), (1, "gamma".to_string(), 0.0),]
        );
        assert_eq!(arr.shape(), &[2, 1]);
        assert_eq!(arr[(0, 0)], 0.0);
        assert_eq!(arr[(1, 0)], 0.0);
        assert_eq!(arr.sum(), 0.0);
        Ok(())
    }

    #[test]
    fn grouped_bed_gc_correction_weights_each_grouped_count() -> Result<()> {
        // Arrange:
        // - one unpaired fragment spans [10,20) with length 10 on the simple ACGT reference
        // - that fragment has GC%=50
        // - `lengths` uses the length-agnostic GC corrector over the full package range, which
        //   averages the two length rows below with the default `equal` weighting:
        //     GC bin [0,51): (3.0 + 1.0) / 2 = 2.0
        // - grouped windows keep `beta` as the counted row and preserve `gamma` as an explicit
        //   zero row
        let bam = single_read_fragment_bam("lengths_grouped_gc", 10, 10)?;
        let reference = simple_reference_twobit()?;
        let out_dir = TempDir::new()?;
        let grouped_bed = out_dir.path().join("grouped_gc.bed");
        let gc_path = out_dir.path().join("grouped_gc_package.npz");
        let package = GCCorrectionPackage {
            version: GC_CORRECTION_SCHEMA_VERSION,
            end_offset: 0,
            length_edges: vec![10, 11, 20],
            gc_edges: vec![0, 51, 100],
            length_bin_frequencies: array![1.0_f64, 1.0_f64],
            reference_contig_footprint: twobit_contig_footprint(&reference.path)?,
            correction_matrix: array![[3.0_f64, 1.0_f64], [1.0_f64, 1.0_f64]],
        };
        package.write_npz(&gc_path)?;
        write_bed(
            &grouped_bed,
            &[("chr1", 10, 20, "beta"), ("chr1", 30, 40, "gamma")],
        )?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        });
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::Any,
        });
        cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
            gc_file: Some(gc_path),
            neutralize_invalid_gc: false,
        });
        cfg.set_gc_length_range(GCLengthRange::Package);
        cfg.set_ref_2bit(Some(reference.path.clone()));
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 10);

        // Act
        run(&cfg)?;
        let counts_path = out_dir.path().join("length_counts.npy");
        let arr: Array2<f64> = read_npy(&counts_path)?;
        let group_index = std::fs::read_to_string(out_dir.path().join("group_index.tsv"))?;

        // Assert
        assert_eq!(
            parse_group_index_tsv(&group_index),
            vec![(0, "beta".to_string()), (1, "gamma".to_string())]
        );
        assert_eq!(arr.shape(), &[2, 1]);
        assert_eq!(arr[(0, 0)], 2.0);
        assert_eq!(arr[(1, 0)], 0.0);
        assert_eq!(arr.sum(), 2.0);
        Ok(())
    }

    #[test]
    fn gc_file_late_tile_window_uses_reference_coordinates_after_fetch_narrowing() -> Result<()> {
        // Arrange:
        // - The counted BED window [930,941) is far from the tile fetch origin 0.
        // - The reference is shorter than the BAM chromosome, but long enough for the narrowed
        //   window-derived fetch span. Reading the full tile reference would overrun the reference.
        // - The one 61 bp unpaired fragment overlaps the window.
        // - Its reference interval [900,961) is all C, so it lands in the high-GC correction bin
        //   with weight 7.0. Using prefix-local origin 0 would see A-only sequence instead.
        let bam = bam_from_specs(
            vec![("chr1".to_string(), 1_500)],
            Vec::new(),
            vec![ReadSpec {
                tid: 0,
                pos: 900,
                cigar: vec![('M', 61)],
                seq: vec![b'A'; 61],
                qual: 40,
                is_reverse: false,
                mapq: 60,
                flags: 0,
                mate_tid: None,
                mate_pos: None,
                insert_size: 0,
            }],
            "lengths_late_tile_gc_origin",
        )?;
        let reference = twobit_from_sequences(
            "lengths_late_tile_gc_origin_ref",
            vec![("chr1".to_string(), late_origin_gc_reference_sequence())],
        )?;
        let out_dir = TempDir::new()?;
        let bed_path = out_dir.path().join("late_window.bed");
        let gc_path = out_dir.path().join("two_bin_gc_package.npz");
        write_bed(&bed_path, &[("chr1", 930, 941, "late")])?;
        write_two_bin_gc_package(
            &gc_path,
            61,
            2.0,
            7.0,
            twobit_contig_footprint(&reference.path)?,
        )?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(bed_path),
            by_grouped_bed: None,
        });
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::Any,
        });
        cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
            gc_file: Some(gc_path),
            neutralize_invalid_gc: false,
        });
        cfg.set_ref_2bit(Some(reference.path.clone()));
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(61, 61);

        // Act
        run(&cfg)?;
        let counts_path = out_dir.path().join("length_counts.npy");
        let arr: Array2<f64> = read_npy(&counts_path)?;

        // Assert
        assert_eq!(arr.shape(), &[1, 1]);
        assert_eq!(arr[(0, 0)], 7.0);
        assert_eq!(arr.sum(), 7.0);
        Ok(())
    }

    #[test]
    fn gc_file_chromosome_end_window_keeps_clamped_fetch_halo() -> Result<()> {
        // Manual derivation:
        // - The BAM and reference chromosome both end at 1022.
        // - The one unpaired fragment spans [961,1022), touching chromosome end and overlapping
        //   the BED window [1010,1022).
        // - Fetch narrowing must clamp the right halo to chrom_len=1022 while still retaining the
        //   full 61 bp fragment needed for GC correction.
        // - The reference interval [961,1022) is all A in `late_origin_gc_reference_sequence()`,
        //   so the fragment lands in the low-GC bin with weight 2.0.
        let bam = bam_from_specs(
            vec![("chr1".to_string(), 1_022)],
            Vec::new(),
            vec![ReadSpec {
                tid: 0,
                pos: 961,
                cigar: vec![('M', 61)],
                seq: vec![b'A'; 61],
                qual: 40,
                is_reverse: false,
                mapq: 60,
                flags: 0,
                mate_tid: None,
                mate_pos: None,
                insert_size: 0,
            }],
            "lengths_chrom_end_gc_fetch_halo",
        )?;
        let reference = twobit_from_sequences(
            "lengths_chrom_end_gc_fetch_halo_ref",
            vec![("chr1".to_string(), late_origin_gc_reference_sequence())],
        )?;
        let out_dir = TempDir::new()?;
        let bed_path = out_dir.path().join("right_end_window.bed");
        let gc_path = out_dir.path().join("two_bin_gc_package.npz");
        write_bed(&bed_path, &[("chr1", 1010, 1022, "right_end")])?;
        write_two_bin_gc_package(
            &gc_path,
            61,
            2.0,
            7.0,
            twobit_contig_footprint(&reference.path)?,
        )?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(bed_path),
            by_grouped_bed: None,
        });
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::Any,
        });
        cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
            gc_file: Some(gc_path),
            neutralize_invalid_gc: false,
        });
        cfg.set_ref_2bit(Some(reference.path.clone()));
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(61, 61);

        run(&cfg)?;
        let arr: Array2<f64> = read_npy(out_dir.path().join("length_counts.npy"))?;

        assert_eq!(arr.shape(), &[1, 1]);
        assert_eq!(arr[(0, 0)], 2.0);
        assert_eq!(arr.sum(), 2.0);
        Ok(())
    }

    #[test]
    fn grouped_bed_assign_when_all_counts_each_containing_window_in_group() -> Result<()> {
        // Arrange:
        // - one unpaired fragment spans [10,20) with length 10
        // - under `assign-by=all`, a window counts only if it fully contains the fragment
        // - grouped intervals:
        //     [10,20) beta  -> contains the fragment
        //     [0,20)  beta  -> also contains the fragment
        //     [10,19) alpha -> does not fully contain the fragment
        //     [30,40) gamma -> zero row
        // - expected grouped counts:
        //     beta  -> 2.0
        //     alpha -> 0.0
        //     gamma -> 0.0
        let bam = single_read_fragment_bam("lengths_grouped_all", 10, 10)?;
        let out_dir = TempDir::new()?;
        let grouped_bed = out_dir.path().join("grouped_all.bed");
        write_bed(
            &grouped_bed,
            &[
                ("chr1", 10, 20, "beta"),
                ("chr1", 0, 20, "beta"),
                ("chr1", 10, 19, "alpha"),
                ("chr1", 30, 40, "gamma"),
            ],
        )?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        });
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::All,
        });
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 10);

        // Act
        run(&cfg)?;
        let counts_path = out_dir.path().join("length_counts.npy");
        let arr: Array2<f64> = read_npy(&counts_path)?;
        let group_index = std::fs::read_to_string(out_dir.path().join("group_index.tsv"))?;

        // Assert
        assert_eq!(
            parse_group_index_tsv(&group_index),
            vec![
                (0, "beta".to_string()),
                (1, "alpha".to_string()),
                (2, "gamma".to_string()),
            ]
        );
        assert_eq!(arr.shape(), &[3, 1]);
        assert_eq!(arr[(0, 0)], 2.0);
        assert_eq!(arr[(1, 0)], 0.0);
        assert_eq!(arr[(2, 0)], 0.0);
        assert_eq!(arr.sum(), 2.0);
        Ok(())
    }

    #[test]
    fn grouped_bed_keeps_zero_group_rows_when_only_filtered_fragments_would_match() -> Result<()> {
        // Arrange:
        // The three-chrom fixture has one fragment per chromosome with lengths 60, 80, and 100.
        // Restricting the allowed range to length 100 should therefore keep only the chr3 group,
        // while chr1 and chr2 still remain as explicit zero rows in grouped mode.
        let bam = three_chrom_length_fixture("lengths_grouped_filtered_zero_rows")?;
        let out_dir = TempDir::new()?;
        let grouped_bed = out_dir.path().join("grouped_filtered_zero_rows.bed");
        write_bed(
            &grouped_bed,
            &[
                ("chr1", 0, 200, "alpha"),
                ("chr2", 0, 200, "beta"),
                ("chr3", 0, 200, "gamma"),
            ],
        )?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1", "chr2", "chr3"]),
        );
        cfg.set_indel_mode(IndelMode::Ignore);
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        });
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::Any,
        });
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(100, 100);

        // Act
        run(&cfg)?;
        let counts_path = out_dir.path().join("length_counts.npy");
        let arr: Array2<f64> = read_npy(&counts_path)?;
        let group_index = std::fs::read_to_string(out_dir.path().join("group_index.tsv"))?;

        // Assert
        assert_eq!(
            parse_group_index_tsv(&group_index),
            vec![
                (0, "alpha".to_string()),
                (1, "beta".to_string()),
                (2, "gamma".to_string()),
            ]
        );
        assert_eq!(arr.shape(), &[3, 1]);
        assert_eq!(arr[(0, 0)], 0.0);
        assert_eq!(arr[(1, 0)], 0.0);
        assert_eq!(arr[(2, 0)], 1.0);
        assert_eq!(arr.sum(), 1.0);
        Ok(())
    }

    #[test]
    fn grouped_bed_assign_when_midpoint_counts_exactly_one_adjacent_group_at_boundary() -> Result<()>
    {
        // Arrange:
        // - unpaired fragment [40,50) has even midpoint 44 or 45
        // - grouped windows split that seam into adjacent groups:
        //     [44,45) alpha
        //     [45,46) beta
        // - endpoint-only control groups must stay zero in midpoint mode
        // - exactly one of alpha/beta must receive the fragment
        let bam = single_read_fragment_bam("lengths_grouped_midpoint_boundary", 40, 10)?;
        let out_dir = TempDir::new()?;
        let grouped_bed = out_dir.path().join("grouped_midpoint_boundary.bed");
        write_bed(
            &grouped_bed,
            &[
                ("chr1", 40, 41, "left_endpoint"),
                ("chr1", 44, 45, "alpha"),
                ("chr1", 45, 46, "beta"),
                ("chr1", 49, 50, "right_endpoint"),
            ],
        )?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        });
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::Midpoint,
        });
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 10);

        // Act
        run(&cfg)?;
        let counts_path = out_dir.path().join("length_counts.npy");
        let arr: Array2<f64> = read_npy(&counts_path)?;

        // Assert
        assert_eq!(arr.shape(), &[4, 1]);
        assert_eq!(arr[(0, 0)], 0.0);
        assert_eq!(arr[(3, 0)], 0.0);
        let alpha_hot = (arr[(1, 0)] - 1.0).abs() < 1e-12 && arr[(2, 0)].abs() < 1e-12;
        let beta_hot = arr[(1, 0)].abs() < 1e-12 && (arr[(2, 0)] - 1.0).abs() < 1e-12;
        assert!(
            alpha_hot || beta_hot,
            "midpoint seam must count exactly one adjacent group, got [{}, {}]",
            arr[(1, 0)],
            arr[(2, 0)]
        );
        assert_eq!(arr.sum(), 1.0);
        Ok(())
    }

    #[test]
    fn grouped_bed_assign_when_midpoint_counts_only_midpoint_group_not_endpoint_groups()
    -> Result<()> {
        // Arrange:
        // - unpaired fragment [10,20) has midpoint 14 or 15, both inside [14,16)
        // - endpoint-only groups are added as negative controls
        // - midpoint mode must therefore count only the midpoint group
        let bam = single_read_fragment_bam("lengths_grouped_midpoint_only_mid", 10, 10)?;
        let out_dir = TempDir::new()?;
        let grouped_bed = out_dir.path().join("grouped_midpoint_only_mid.bed");
        write_bed(
            &grouped_bed,
            &[
                ("chr1", 10, 11, "left_endpoint"),
                ("chr1", 14, 16, "mid"),
                ("chr1", 19, 20, "right_endpoint"),
                ("chr1", 30, 31, "gamma"),
            ],
        )?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        });
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::Midpoint,
        });
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 10);

        // Act
        run(&cfg)?;
        let counts_path = out_dir.path().join("length_counts.npy");
        let arr: Array2<f64> = read_npy(&counts_path)?;
        let group_index = std::fs::read_to_string(out_dir.path().join("group_index.tsv"))?;

        // Assert
        assert_eq!(
            parse_group_index_tsv(&group_index),
            vec![
                (0, "left_endpoint".to_string()),
                (1, "mid".to_string()),
                (2, "right_endpoint".to_string()),
                (3, "gamma".to_string()),
            ]
        );
        assert_eq!(arr.shape(), &[4, 1]);
        assert_eq!(arr[(0, 0)], 0.0);
        assert_eq!(arr[(1, 0)], 1.0);
        assert_eq!(arr[(2, 0)], 0.0);
        assert_eq!(arr[(3, 0)], 0.0);
        assert_eq!(arr.sum(), 1.0);
        Ok(())
    }

    #[test]
    fn grouped_bed_assign_when_proportion_counts_only_groups_meeting_threshold() -> Result<()> {
        // Arrange:
        // - one unpaired fragment spans [10,20) with length 10
        // - with `proportion=0.5`, windows with at least 5 fragment bp qualify
        // - grouped intervals:
        //     [10,15) beta  -> 5/10, counts
        //     [15,20) beta  -> 5/10, counts
        //     [10,14) alpha -> 4/10, rejected
        //     [30,40) gamma -> zero row
        // - expected grouped counts:
        //     beta  -> 2.0
        //     alpha -> 0.0
        //     gamma -> 0.0
        let bam = single_read_fragment_bam("lengths_grouped_proportion", 10, 10)?;
        let out_dir = TempDir::new()?;
        let grouped_bed = out_dir.path().join("grouped_proportion.bed");
        write_bed(
            &grouped_bed,
            &[
                ("chr1", 10, 15, "beta"),
                ("chr1", 15, 20, "beta"),
                ("chr1", 10, 14, "alpha"),
                ("chr1", 30, 40, "gamma"),
            ],
        )?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        });
        cfg.set_window_assignment(AssignToWindowArgs {
            assign_by: WindowAssigner::Proportion(0.5),
        });
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 10);

        // Act
        run(&cfg)?;
        let counts_path = out_dir.path().join("length_counts.npy");
        let arr: Array2<f64> = read_npy(&counts_path)?;
        let group_index = std::fs::read_to_string(out_dir.path().join("group_index.tsv"))?;

        // Assert
        assert_eq!(
            parse_group_index_tsv(&group_index),
            vec![
                (0, "beta".to_string()),
                (1, "alpha".to_string()),
                (2, "gamma".to_string()),
            ]
        );
        assert_eq!(arr.shape(), &[3, 1]);
        assert_eq!(arr[(0, 0)], 2.0);
        assert_eq!(arr[(1, 0)], 0.0);
        assert_eq!(arr[(2, 0)], 0.0);
        assert_eq!(arr.sum(), 2.0);
        Ok(())
    }

    #[test]
    fn grouped_bed_errors_when_group_name_column_is_missing() -> Result<()> {
        // Arrange:
        // - grouped BED mode requires a fourth column naming the group
        // - a three-column BED row is therefore invalid and should fail loudly
        let bam = single_read_fragment_bam("lengths_grouped_missing_group_name", 10, 10)?;
        let out_dir = TempDir::new()?;
        let grouped_bed = out_dir.path().join("grouped_missing_name.bed");
        std::fs::write(&grouped_bed, "chr1\t10\t20\n")?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        });
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 10);

        // Act
        let err = run(&cfg).expect_err("grouped BED without a group column should fail");

        // Assert
        assert!(err.to_string().contains("missing group name"));
        Ok(())
    }

    #[test]
    fn grouped_bed_errors_when_no_windows_survive_selected_chromosomes() -> Result<()> {
        // Arrange:
        // - grouped BED contains a valid group, but only on chr2
        // - the run is restricted to chr1, so grouped mode has no usable groups at all
        let bam = single_read_fragment_bam("lengths_grouped_empty_after_filtering", 10, 10)?;
        let out_dir = TempDir::new()?;
        let grouped_bed = out_dir.path().join("grouped_wrong_chr.bed");
        write_bed(&grouped_bed, &[("chr2", 10, 20, "beta")])?;

        let mut cfg = LengthsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.path().to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        });
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.set_per_bp_length_bins(10, 10);

        // Act
        let err =
            run(&cfg).expect_err("grouped BED with no selected-chromosome windows should fail");

        // Assert
        assert!(err.to_string().contains(
            "grouped BED file did not contain any valid windows on the selected chromosomes"
        ));
        Ok(())
    }
}

mod tests_lengths_tiling_reducer {

    #![cfg(feature = "cmd_lengths")]

    use anyhow::Result;
    use cfdnalab::commands::lengths::counting::{LengthAxis, LengthCounts};
    use cfdnalab::commands::lengths::tiling::{
        reduce_partials_for_chr, write_cross_npy, write_partials_npz,
    };
    use ndarray::{Array1, Array2, ShapeBuilder};
    use ndarray_npy::NpzWriter;
    use std::sync::Arc;
    use std::{fs::File, path::PathBuf};
    use tempfile::TempDir;

    fn exact_axis(min_length: u32, max_length: u32) -> Arc<LengthAxis> {
        let edges: Vec<u32> = (min_length..=max_length + 1).collect();
        Arc::new(LengthAxis::new(edges).expect("test length axis should be valid"))
    }

    fn template_counts() -> LengthCounts {
        LengthCounts::new(exact_axis(10, 10))
    }

    fn counts_with_value(val: f64) -> LengthCounts {
        let mut lc = template_counts();
        lc.counts[0] = val;
        lc
    }

    fn expect_written_path(path: Option<PathBuf>, label: &str) -> PathBuf {
        path.unwrap_or_else(|| panic!("{label} should have been written"))
    }

    #[test]
    fn reducer_accepts_contained_only() -> Result<()> {
        // Human verification status: unverified
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let template = template_counts();

        let counts = vec![counts_with_value(2.0)];
        let contained = vec![true];
        let partial_path = expect_written_path(
            write_partials_npz(dir, "partials", "chr1", 0, &[0], &contained, &counts)?,
            "contained partial",
        );
        // No cross file because window is contained

        let reduced = reduce_partials_for_chr("chr1", &[partial_path], &[], 1, &template)?;
        assert_eq!(reduced.len(), 1);
        assert!((reduced[0].counts[0] - 2.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn partial_writer_rejects_count_rows_without_matching_window_indices() {
        // Human verification status: unverified
        let tmp = TempDir::new().unwrap();
        let dir = tmp.path();

        let counts = vec![counts_with_value(2.0)];
        let contained = vec![true, false];

        let err = write_partials_npz(dir, "partials", "chr1", 0, &[0, 1], &contained, &counts)
            .expect_err("partial writer should reject mismatched row metadata");

        assert!(err.to_string().contains("counts length mismatch"));
    }

    #[test]
    fn reducer_counts_multiple_crossing_tiles() -> Result<()> {
        // Human verification status: unverified
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let template = template_counts();

        let counts = vec![counts_with_value(1.0)];
        let contained = vec![false];
        // Two tiles, both crossing the same window
        let partial_paths = vec![
            expect_written_path(
                write_partials_npz(dir, "partials", "chr1", 0, &[0], &contained, &counts)?,
                "first crossing partial",
            ),
            expect_written_path(
                write_partials_npz(dir, "partials", "chr1", 1, &[0], &contained, &counts)?,
                "second crossing partial",
            ),
        ];
        let cross_paths = vec![
            expect_written_path(
                write_cross_npy(dir, "cross", "chr1", 0, &[0])?,
                "first cross",
            ),
            expect_written_path(
                write_cross_npy(dir, "cross", "chr1", 1, &[0])?,
                "second cross",
            ),
        ];

        let reduced = reduce_partials_for_chr(
            "chr1",
            partial_paths.as_slice(),
            cross_paths.as_slice(),
            1,
            &template,
        )?;
        assert_eq!(reduced.len(), 1);
        assert!((reduced[0].counts[0] - 2.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn reducer_combines_contained_and_cross() -> Result<()> {
        // Human verification status: unverified
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let template = template_counts();

        let contained_counts = vec![counts_with_value(1.0)];
        let crossing_counts = vec![counts_with_value(3.0)];
        let partial_paths = vec![
            expect_written_path(
                write_partials_npz(dir, "partials", "chr1", 0, &[0], &[true], &contained_counts)?,
                "contained partial",
            ),
            expect_written_path(
                write_partials_npz(dir, "partials", "chr1", 1, &[0], &[false], &crossing_counts)?,
                "crossing partial",
            ),
        ];
        let cross_paths = vec![expect_written_path(
            write_cross_npy(dir, "cross", "chr1", 1, &[0])?,
            "cross index",
        )];

        let reduced = reduce_partials_for_chr(
            "chr1",
            partial_paths.as_slice(),
            cross_paths.as_slice(),
            1,
            &template,
        )?;
        assert_eq!(reduced.len(), 1);
        // Expect 1 contained contribution and 1 crossing contribution => sum counts
        assert!((reduced[0].counts[0] - 4.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn reducer_errors_when_contribution_missing() {
        // Human verification status: unverified
        let template = template_counts();

        // No partials written -> zero contributions
        let err = reduce_partials_for_chr("chr1", &[], &[], 1, &template)
            .expect_err("should fail when contributions are missing");
        assert!(err.to_string().contains("expected 1"));
    }

    #[test]
    fn reducer_errors_on_mismatched_counts() {
        // Human verification status: unverified
        let tmp = TempDir::new().unwrap();
        let dir = tmp.path();
        let template = template_counts();

        // Cross file claims one contribution, but no partial exists
        let cross_path = expect_written_path(
            write_cross_npy(dir, "cross", "chr1", 0, &[0]).unwrap(),
            "cross",
        );

        let err = reduce_partials_for_chr("chr1", &[], &[cross_path], 1, &template)
            .expect_err("should fail when expected contributions not met");
        assert!(err.to_string().contains("expected 1"));
    }

    #[test]
    fn reducer_errors_on_counts_width_mismatch() {
        // Human verification status: unverified
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

        let err = reduce_partials_for_chr("chr1", &[path], &[], 1, &template)
            .expect_err("should fail on counts width mismatch");
        assert!(err.to_string().contains("counts width mismatch"));
    }

    #[test]
    fn reducer_errors_on_non_contiguous_counts_rows() {
        // Human verification status: unverified
        let tmp = TempDir::new().unwrap();
        let dir = tmp.path();
        let template = LengthCounts::new(exact_axis(10, 11)); // two-length template

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

        let err = reduce_partials_for_chr("chr1", &[path], &[], 1, &template)
            .expect_err("should fail on non-contiguous counts rows");
        assert!(err.to_string().contains("counts row not contiguous"));
    }

    #[test]
    fn reducer_ignores_files_from_other_chromosomes() -> Result<()> {
        // Human verification status: unverified
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let template = template_counts();

        let counts = vec![counts_with_value(1.0)];
        let contained = vec![true];
        let partial_path = expect_written_path(
            write_partials_npz(dir, "partials", "chr1", 0, &[0], &contained, &counts)?,
            "chr1 partial",
        );

        // Stray files for another chromosome are not passed to the reducer because the command
        // now records exact output paths from each tile.
        write_partials_npz(dir, "partials", "chr2", 0, &[0], &contained, &counts)?;
        write_cross_npy(dir, "cross", "chr2", 0, &[0])?;

        let reduced = reduce_partials_for_chr("chr1", &[partial_path], &[], 1, &template)?;
        assert_eq!(reduced.len(), 1);
        assert!((reduced[0].counts[0] - 1.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn reducer_uses_explicit_paths_for_overlapping_chromosome_names() -> Result<()> {
        // Human verification status: verified
        // Manual expectations:
        // - `chr1` contributes one contained count with value 1.
        // - `chr1.extra` contributes one contained count with value 5.
        // - Old dotted-substring discovery for `chr1` could also match `partials.chr1.extra.0.npz`.
        // - Passing explicit paths means the `chr1.extra` file is ignored unless the caller provides
        //   it, so the reduced `chr1` count remains 1.
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let template = template_counts();

        let contained = vec![true];
        let chr1_counts = vec![counts_with_value(1.0)];
        let chr1_extra_counts = vec![counts_with_value(5.0)];
        let chr1_partial_path = expect_written_path(
            write_partials_npz(dir, "partials", "chr1", 0, &[0], &contained, &chr1_counts)?,
            "chr1 partial",
        );
        let _chr1_extra_partial_path = expect_written_path(
            write_partials_npz(
                dir,
                "partials",
                "chr1.extra",
                0,
                &[0],
                &contained,
                &chr1_extra_counts,
            )?,
            "chr1.extra partial",
        );

        let reduced = reduce_partials_for_chr("chr1", &[chr1_partial_path], &[], 1, &template)?;

        assert_eq!(reduced.len(), 1);
        assert!((reduced[0].counts[0] - 1.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn write_partials_rejects_mismatched_contained() {
        // Human verification status: unverified
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
        // Human verification status: unverified
        let tmp = TempDir::new().unwrap();
        let dir = tmp.path();
        let template = template_counts();

        let counts = vec![counts_with_value(1.0)];
        let contained = vec![false];
        // Write a partial with idx outside n_windows=1
        let partial_path = expect_written_path(
            write_partials_npz(dir, "partials", "chr1", 0, &[2], &contained, &counts).unwrap(),
            "out-of-bounds partial",
        );

        let err = reduce_partials_for_chr("chr1", &[partial_path], &[], 1, &template)
            .expect_err("should fail on out-of-bounds idx");
        assert!(err.to_string().contains("out of bounds"));
    }

    #[test]
    fn reducer_errors_on_out_of_bounds_cross_idx() {
        // Human verification status: unverified
        let tmp = TempDir::new().unwrap();
        let dir = tmp.path();
        let template = template_counts();

        let cross_path = expect_written_path(
            write_cross_npy(dir, "cross", "chr1", 0, &[3]).unwrap(),
            "out-of-bounds cross",
        );
        let err = reduce_partials_for_chr("chr1", &[], &[cross_path], 1, &template)
            .expect_err("should fail on cross index out of bounds");
        assert!(err.to_string().contains("Cross index"));
    }

    #[test]
    fn reducer_separates_windows() -> Result<()> {
        // Human verification status: unverified
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let template = template_counts();

        let counts0 = vec![counts_with_value(1.0)];
        let counts1 = vec![counts_with_value(2.0)];
        let contained = vec![true];

        // Window 0 contained in tile 0, window 1 contained in tile 1
        let partial_paths = vec![
            expect_written_path(
                write_partials_npz(dir, "partials", "chr1", 0, &[0], &contained, &counts0)?,
                "first window partial",
            ),
            expect_written_path(
                write_partials_npz(dir, "partials", "chr1", 1, &[1], &contained, &counts1)?,
                "second window partial",
            ),
        ];

        let reduced = reduce_partials_for_chr("chr1", partial_paths.as_slice(), &[], 2, &template)?;
        assert_eq!(reduced.len(), 2);
        assert!((reduced[0].counts[0] - 1.0).abs() < 1e-6);
        assert!((reduced[1].counts[0] - 2.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn write_partials_skips_empty() -> Result<()> {
        // Human verification status: unverified
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let template = template_counts();
        let res = write_partials_npz(dir, "partials", "chr1", 0, &[], &[], &[])?;
        assert!(res.is_none());
        // Ensure reducer still errors because nothing was written
        let err = reduce_partials_for_chr("chr1", &[], &[], 1, &template)
            .expect_err("should fail when nothing written");
        assert!(err.to_string().contains("expected 1"));
        Ok(())
    }

    #[test]
    fn write_cross_skips_empty() -> Result<()> {
        // Human verification status: unverified
        let tmp = TempDir::new()?;
        let dir = tmp.path();
        let res = write_cross_npy(dir, "cross", "chr1", 0, &[])?;
        assert!(res.is_none());
        Ok(())
    }
}

mod tests_lengths_tiling_helpers {

    use cfdnalab::commands::cli_common::WindowSpec;
    use cfdnalab::shared::bam::Contigs;
    use cfdnalab::shared::interval::IndexedInterval;
    use cfdnalab::shared::tiled_run::{Tile, TileWindowSpan, build_tiles};
    use cfdnalab::shared::window_fetch::{BedFetchPolicy, fetch_span_for_tile};
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
        // Human verification status: unverified
        // Tile: core 50-150, fetch 30-200 (halo 20 left, 50 right), chrom len 180
        let tile = Tile::from_coords("chr1".to_string(), 0, 0, 50, 150, 30, 200)
            .expect("test tile should be valid");
        let span = fetch_span_for_tile(
            &tile,
            None,
            None,
            &WindowSpec::Size(100),
            180,
            0,
            BedFetchPolicy::CandidateWindowExtent,
        )
        .expect("span expected")
        .expect("fetch span expected");
        // Window span touching core: 0..200, after halo clamp -> 30..180
        assert_eq!(span.start(), 30);
        assert_eq!(span.end(), 180);
    }

    #[test]
    fn build_tiles_aligns_to_bin_when_divisible() {
        // Human verification status: unverified
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
        // Human verification status: unverified
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
        // Human verification status: unverified
        let tile = Tile::from_coords("chr1".to_string(), 0, 0, 0, 50, 0, 200)
            .expect("test tile should be valid");
        let span = fetch_span_for_tile(
            &tile,
            None,
            None,
            &WindowSpec::Global,
            120,
            0,
            BedFetchPolicy::CandidateWindowExtent,
        )
        .expect("span")
        .expect("fetch span expected");
        assert_eq!(span.start(), 0);
        assert_eq!(span.end(), 120);
    }

    #[test]
    fn fetch_span_for_tile_bed_with_overlap() {
        // Human verification status: unverified
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
            0,
            BedFetchPolicy::CoreOverlap,
        )
        .expect("span")
        .expect("fetch span expected");
        // min_ws=90, max_we=170, halos: left 20, right 40 -> widened to 70..210, clamped to fetch
        assert_eq!(res.start(), 80);
        assert_eq!(res.end(), 200);
    }

    #[test]
    fn fetch_span_bed_none_when_no_overlap() {
        // Human verification status: unverified
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
            0,
            BedFetchPolicy::CoreOverlap,
        )
        .expect("fetch span computation should succeed");
        assert!(res.is_none());
    }

    #[test]
    fn fetch_span_size_mode_none_when_tile_right_of_chromosome() {
        // Human verification status: unverified
        let tile = Tile::from_coords("chr1".to_string(), 0, 0, 250, 260, 230, 270)
            .expect("test tile should be valid");
        let res = fetch_span_for_tile(
            &tile,
            None,
            None,
            &WindowSpec::Size(50),
            200,
            0,
            BedFetchPolicy::CandidateWindowExtent,
        )
        .expect("fetch span computation should succeed");
        assert!(res.is_none());
    }

    #[test]
    fn tile_constructor_rejects_empty_core() {
        // Human verification status: unverified
        let err = Tile::from_coords("chr1".to_string(), 0, 0, 100, 100, 80, 120).unwrap_err();
        assert!(format!("{err}").contains("interval end (100) must be greater than start (100)"));
    }
}
