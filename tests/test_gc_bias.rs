#![cfg(all(feature = "cmd_gc_bias", feature = "cmd_ref_gc_bias"))]

mod fixtures;

mod tests_gc_bias {
    use crate::fixtures;
    use anyhow::Result;
    use fxhash::FxHashMap;
    use ndarray::array;
    use ndarray_npy::{NpzWriter, read_npy};
    use tempfile::{TempDir, tempdir};

    use cfdnalab::commands::{
        cli_common::{ChromosomeArgs, GCWindowsArgs, IOCArgs, LoggingArgs, Ref2BitRequiredArgs},
        gc_bias::{
            GC_CORRECTION_SCHEMA_VERSION,
            binning::{BinnedAxis, bins_from_edges, compute_bin_edges},
            config::GCConfig,
            correct::{GCCorrector, LengthAgnosticGCCorrector, MarginalizeLengthsWeightingScheme},
            counting::gc_percent_widths,
            gc_bias::{interpolate_masked_corrections, run as run_gc_bias},
            load_reference_bias::load_reference_gc_data,
            outliers::{
                OutlierAction, OutlierRule, OutlierScope, OutlierStats, apply_outliers_to_matrix,
                interpolated_quantile, outlier_bounds,
            },
            package::GCCorrectionPackage,
            support_masking::build_extreme_bins_support_mask,
        },
        ref_gc_bias::{config::RefGCBiasConfig, ref_gc_bias::run as run_ref_gc_bias},
    };

    const GC_COMMAND_F64_TOL: f64 = 1e-6;

    fn assert_gc_command_close(actual: f64, expected: f64, context: &str) {
        // The outlier helpers estimate quantiles/bounds in `f32` and only then write the matrix
        // back as `f64`, so the stable contract here is "matches the hand-derived value within the
        // command's float precision", not bit-exact `f64` arithmetic on the ideal decimal.
        assert!(
            (actual - expected).abs() <= GC_COMMAND_F64_TOL,
            "{context}: expected {expected}, got {actual}"
        );
    }

    #[test]
    fn masks_extreme_gc_bins_per_side_in_square_matrix() {
        // Human verification status: unverified
        // Arrange: 6x6 matrix with two extreme GC bins on each side.
        let expected = array![
            [false, false, true, true, false, false],
            [false, false, true, true, false, false],
            [false, false, true, true, false, false],
            [false, false, true, true, false, false],
            [false, false, true, true, false, false],
            [false, false, true, true, false, false],
        ];

        // Act: build the support mask after binning.
        let mask = build_extreme_bins_support_mask((6, 6), 2, 0);

        // Assert: the central two GC bins remain supported across all lengths.
        assert_eq!(mask, expected);
    }

    #[test]
    fn masks_shortest_length_bins_in_matrix() {
        // Human verification status: unverified
        // Arrange: 5x4 matrix with one shortest length bin masked.
        let expected = array![
            [false, false, false, false],
            [true, true, true, true],
            [true, true, true, true],
            [true, true, true, true],
            [true, true, true, true],
        ];

        // Act: build the support mask after binning.
        let mask = build_extreme_bins_support_mask((5, 4), 0, 1);

        // Assert: the central three length bins remain supported across all GC bins.
        assert_eq!(mask, expected);
    }

    #[test]
    fn interpolates_masked_short_length_row() -> Result<()> {
        // Human verification status: unverified
        // Arrange: first length row is masked; other rows are supported.
        let mut matrix = array![
            [0.0_f64, 0.0_f64],
            [2.0_f64, 2.0_f64],
            [4.0_f64, 4.0_f64],
            [6.0_f64, 6.0_f64],
        ];
        let mask = build_extreme_bins_support_mask((4, 2), 0, 1);

        // Act: interpolate masked bins.
        interpolate_masked_corrections(&mut matrix, &mask)?;

        // Assert:
        // - the masked first row is filled from the nearest supported row
        // - the supported rows remain unchanged
        let expected = array![
            [2.0_f64, 2.0_f64],
            [2.0_f64, 2.0_f64],
            [4.0_f64, 4.0_f64],
            [6.0_f64, 6.0_f64],
        ];
        assert_eq!(matrix, expected);
        Ok(())
    }

    #[test]
    fn round_trips_bins_to_edges_and_back() {
        // Human verification status: unverified
        // Arrange: build a simple BinnedAxis where bins group indices as [0-1], [2-4], and [5-7].
        let mut index_to_bin = FxHashMap::default();
        let mut bin_to_indices = FxHashMap::default();
        let bins: [Vec<usize>; 3] = [vec![0, 1], vec![2, 3, 4], vec![5, 6, 7]];
        for (bin_idx, indices) in bins.iter().enumerate() {
            bin_to_indices.insert(bin_idx, indices.clone());
            for &idx in indices {
                index_to_bin.insert(idx, bin_idx);
            }
        }
        let axis = BinnedAxis {
            index_to_bin,
            bin_to_indices,
            num_bins: 3,
        };

        // Act: compute edges then reconstruct the bins.
        let edges = compute_bin_edges(&axis, 0, 7).expect("edges should be computed");
        let reconstructed_axis = bins_from_edges(edges.as_slice()).expect("rebuild should work");

        // Assert: the derived edges match the expected bin boundaries, and the reconstructed
        // axis matches the original bin layout.
        assert_eq!(edges, vec![0, 2, 5, 7]);
        assert_eq!(reconstructed_axis.num_bins, axis.num_bins);
        assert_eq!(reconstructed_axis.bin_to_indices, axis.bin_to_indices);
        assert_eq!(reconstructed_axis.index_to_bin, axis.index_to_bin);
    }

    #[test]
    fn provides_expected_weights_after_roundtrip() -> Result<()> {
        // Human verification status: unverified
        // Arrange: a package whose edges start at non-zero values so offset logic is exercised.
        let length_edges = vec![30, 34, 40];
        let gc_edges = vec![10, 60, 90];
        let correction_matrix = array![[1.5_f64, 2.0_f64], [0.5_f64, 0.75_f64]];
        let package = GCCorrectionPackage {
            version: GC_CORRECTION_SCHEMA_VERSION,
            end_offset: 3,
            length_edges: length_edges.clone(),
            gc_edges: gc_edges.clone(),
            correction_matrix,
            length_bin_frequencies: array![1.0_f64, 1.0_f64],
        };
        let tmp_dir = tempdir()?;
        let pkg_path = tmp_dir.path().join("gc_package.npz");
        package.write_npz(&pkg_path)?;

        // Act: load the package and build a corrector.
        let loaded = GCCorrectionPackage::from_file(&pkg_path)?;
        let corrector = GCCorrector::from_package(&loaded)?;

        // Assert: fragments landing in each bin retrieve the expected weights.
        let weight_len31_gc20 = corrector.get_correction_weight(31, 20)?;
        assert!(
            (weight_len31_gc20 - 1.5).abs() < f64::EPSILON,
            "length 31 / GC 20 should map to 1.5"
        );

        let weight_len32_gc70 = corrector.get_correction_weight(32, 70)?;
        assert!(
            (weight_len32_gc70 - 2.0).abs() < f64::EPSILON,
            "length 32 / GC 70 should map to 2.0"
        );

        let weight_len38_gc55 = corrector.get_correction_weight(38, 55)?;
        assert!(
            (weight_len38_gc55 - 0.5).abs() < f64::EPSILON,
            "length 38 / GC 55 should map to 0.5"
        );

        let weight_len39_gc80 = corrector.get_correction_weight(39, 80)?;
        assert!(
            (weight_len39_gc80 - 0.75).abs() < f64::EPSILON,
            "length 39 / GC 80 should map to 0.75"
        );

        Ok(())
    }

    #[test]
    fn gc_correction_package_rejects_missing_path_before_opening() -> Result<()> {
        // Human verification status: unverified
        // Arrange: point the loader at a `.npz` path that does not exist.
        let tmp_dir = tempdir()?;
        let missing_path = tmp_dir.path().join("missing_gc_package.npz");

        // Act: try to load the missing package.
        let err = GCCorrectionPackage::from_file(&missing_path)
            .expect_err("missing GC correction package should fail");

        // Assert: the user gets the shared "existing .npz file" contract directly instead of a
        // lower-level IO or NPZ parsing error.
        let msg = err.to_string();
        assert!(
            msg.contains("must point to an existing .npz file"),
            "unexpected error message: {msg}"
        );
        assert!(
            msg.contains("missing_gc_package.npz"),
            "unexpected error message: {msg}"
        );

        Ok(())
    }

    #[test]
    fn gc_correction_package_rejects_non_npz_extension_before_parsing() -> Result<()> {
        // Human verification status: unverified
        // Arrange: create a regular file with the wrong extension.
        let tmp_dir = tempdir()?;
        let wrong_extension_path = tmp_dir.path().join("gc_package.txt");
        std::fs::write(&wrong_extension_path, b"not an npz archive")?;

        // Act: try to load the non-`.npz` file.
        let err = GCCorrectionPackage::from_file(&wrong_extension_path)
            .expect_err("wrong extension should fail");

        // Assert: extension validation runs before the NPZ reader.
        let msg = err.to_string();
        assert!(
            msg.contains("must point to a .npz file"),
            "unexpected error message: {msg}"
        );
        assert!(
            msg.contains("gc_package.txt"),
            "unexpected error message: {msg}"
        );

        Ok(())
    }

    fn make_length_agnostic_package() -> GCCorrectionPackage {
        let correction_matrix = array![[1.0_f64, 2.0_f64], [3.0_f64, 5.0_f64]];
        GCCorrectionPackage {
            version: GC_CORRECTION_SCHEMA_VERSION,
            end_offset: 0,
            length_edges: vec![20, 30, 40],
            gc_edges: vec![0, 50, 100],
            correction_matrix,
            length_bin_frequencies: array![0.2_f64, 0.8_f64],
        }
    }

    #[test]
    fn length_agnostic_equal_weighting_means_rows() -> Result<()> {
        // Human verification status: unverified
        let package = make_length_agnostic_package();
        let corrector = GCCorrector::from_package(&package)?;
        let agnostic = LengthAgnosticGCCorrector::from_gc_corrector(
            &corrector,
            &MarginalizeLengthsWeightingScheme::Equal,
        )?;

        for gc_pct in 0..50 {
            assert!(
                (agnostic.get_correction_weight(gc_pct)? - 2.0).abs() < 1e-12,
                "equal weighting should map GC% {gc_pct} into the first GC bin"
            );
        }
        for gc_pct in 50..=100 {
            assert!(
                (agnostic.get_correction_weight(gc_pct)? - 3.5).abs() < 1e-12,
                "equal weighting should map GC% {gc_pct} into the second GC bin"
            );
        }
        Ok(())
    }

    #[test]
    fn length_agnostic_coverage_weighting_uses_frequencies() -> Result<()> {
        // Human verification status: unverified
        let package = make_length_agnostic_package();
        let corrector = GCCorrector::from_package(&package)?;
        let agnostic = LengthAgnosticGCCorrector::from_gc_corrector(
            &corrector,
            &MarginalizeLengthsWeightingScheme::Coverage,
        )?;

        // Weighted average with frequencies [0.2, 0.8]
        for gc_pct in 0..50 {
            assert!(
                (agnostic.get_correction_weight(gc_pct)? - 2.6).abs() < 1e-12,
                "coverage weighting should map GC% {gc_pct} into the first GC bin"
            );
        }
        for gc_pct in 50..=100 {
            assert!(
                (agnostic.get_correction_weight(gc_pct)? - 4.4).abs() < 1e-12,
                "coverage weighting should map GC% {gc_pct} into the second GC bin"
            );
        }
        Ok(())
    }

    #[test]
    fn length_agnostic_max_coverage_picks_most_frequent_row() -> Result<()> {
        // Human verification status: unverified
        let package = make_length_agnostic_package();
        let corrector = GCCorrector::from_package(&package)?;
        let agnostic = LengthAgnosticGCCorrector::from_gc_corrector(
            &corrector,
            &MarginalizeLengthsWeightingScheme::MaxCoverage,
        )?;

        // Row with highest frequency is [3.0, 5.0]
        for gc_pct in 0..50 {
            assert!(
                (agnostic.get_correction_weight(gc_pct)? - 3.0).abs() < 1e-12,
                "max-coverage weighting should map GC% {gc_pct} into the first GC bin"
            );
        }
        for gc_pct in 50..=100 {
            assert!(
                (agnostic.get_correction_weight(gc_pct)? - 5.0).abs() < 1e-12,
                "max-coverage weighting should map GC% {gc_pct} into the second GC bin"
            );
        }
        Ok(())
    }

    fn write_reference_package_for_single_length(
        reference_path: &std::path::Path,
        out_dir: &TempDir,
        fragment_length: u32,
        end_offset: u8,
    ) -> Result<()> {
        let cfg = RefGCBiasConfig {
            ref_genome: Ref2BitRequiredArgs {
                ref_2bit: reference_path.to_path_buf(),
            },
            output_dir: out_dir.path().to_path_buf(),
            output_prefix: String::new(),
            n_threads: 1,
            // These tests use small synthetic references, for example:
            // - `simple_reference_twobit()` is 256 bp
            // - the smallest custom GC-bias fixture here is 200 bp
            //
            // `ref-gc-bias` samples fragment starts from the set of valid start positions,
            // so `n_positions` must stay below that count:
            //   valid_starts = chrom_len - fragment_length + 1
            //
            // Using 100 keeps the helper valid for all current fixtures while still exercising
            // the full producer -> consumer path. The exact number is not part of the behavior
            // under test in this file.
            n_positions: 100,
            seed: Some(7),
            windows: Default::default(),
            chromosomes: ChromosomeArgs {
                chromosomes: Some(vec!["chr1".to_string()]),
                chromosomes_file: None,
            },
            blacklist: None,
            fragment_lengths: cfdnalab::commands::cli_common::FragmentLengthArgs {
                min_fragment_length: fragment_length,
                max_fragment_length: fragment_length,
            },
            end_offset,
            skip_interpolation: true,
            smoothing_sigma: 0.55,
            smoothing_radius: 2,
            skip_smoothing: true,
            tile_size: 1_000_000,
            logging: LoggingArgs::default(),
        };
        run_ref_gc_bias(&cfg)
    }

    fn make_gc_bias_cfg(
        bam_path: &std::path::Path,
        reference_path: &std::path::Path,
        ref_gc_dir: &std::path::Path,
        output_dir: &std::path::Path,
    ) -> GCConfig {
        let ioc = IOCArgs {
            bam: bam_path.to_path_buf(),
            output_dir: output_dir.to_path_buf(),
            n_threads: 1,
        };
        let mut cfg = GCConfig::new(
            ioc,
            reference_path.to_path_buf(),
            ref_gc_dir.join("ref_gc_package.npz"),
            ChromosomeArgs {
                chromosomes: Some(vec!["chr1".to_string()]),
                chromosomes_file: None,
            },
        );
        cfg.set_min_mapq(0);
        cfg.set_tile_size(1_000_000);
        cfg.set_min_window_acgt_pct(0);
        cfg.set_save_intermediates(true);
        cfg
    }

    #[test]
    fn gc_bias_default_min_mapq_matches_explicit_thirty_and_differs_from_explicit_zero()
    -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Use the repeated 256 bp ACGT reference from `simple_reference_twobit()`.
        // For fragment length 60, every 60 bp fragment contains exactly 30 GC bases because:
        // - the reference repeats a 4 bp cycle with 2 GC bases per cycle
        // - 60 = 15 * 4, so each fragment spans exactly 15 full cycles
        // - therefore GC fraction is 30 / 60 = 50%, i.e. GC% bin 50
        //
        // We then place three 60 bp fragments with different MAPQ:
        // - fragment A: MAPQ 60
        // - fragment B: MAPQ 0
        // - fragment C: MAPQ 30
        //
        // Run in global mode with `save_intermediates=true`, because global mode keeps one
        // combined raw count table instead of per-window mean scaling. The saved
        // `avg_cfdna_counts` matrix should therefore be raw surviving counts:
        // - default `min_mapq = 30`: 2 counts at GC% 50
        // - explicit `min_mapq = 30`: same as default
        // - explicit `min_mapq = 0`: 3 counts at GC% 50
        let reference = fixtures::simple_reference_twobit()?;
        let fragment_with_mapq = |start: i64, mapq: u8| {
            let mut fragment = fixtures::paired_fragment(start, 60, 20);
            fragment.forward.mapq = mapq;
            fragment.reverse.mapq = mapq;
            fragment
        };
        let bam = fixtures::bam_from_specs(
            vec![("chr1".to_string(), 256)],
            vec![
                fragment_with_mapq(20, 60),
                fragment_with_mapq(80, 0),
                fragment_with_mapq(140, 30),
            ],
            Vec::new(),
            "gc_bias_default_min_mapq",
        )?;
        let ref_gc_dir = TempDir::new()?;
        write_reference_package_for_single_length(&reference.path, &ref_gc_dir, 60, 0)?;

        let out_default = TempDir::new()?;
        let out_thirty = TempDir::new()?;
        let out_zero = TempDir::new()?;

        let make_cfg = |output_dir: &std::path::Path| {
            let ioc = IOCArgs {
                bam: bam.bam.clone(),
                output_dir: output_dir.to_path_buf(),
                n_threads: 1,
            };
            let mut cfg = GCConfig::new(
                ioc,
                reference.path.clone(),
                ref_gc_dir.path().join("ref_gc_package.npz"),
                ChromosomeArgs {
                    chromosomes: Some(vec!["chr1".to_string()]),
                    chromosomes_file: None,
                },
            );
            cfg.set_windows(GCWindowsArgs {
                by_size: None,
                by_bed: None,
                global: true,
            });
            cfg.set_tile_size(1_000_000);
            cfg.set_min_window_acgt_pct(0);
            cfg.set_num_extreme_gc_bins(0);
            cfg.set_num_short_length_bins(0);
            cfg.outlier_method = cfdnalab::commands::gc_bias::config::OutlierMethodArg::None;
            cfg.set_save_intermediates(true);
            cfg
        };

        let default_cfg = make_cfg(out_default.path());
        let mut explicit_thirty_cfg = make_cfg(out_thirty.path());
        explicit_thirty_cfg.set_min_mapq(30);
        let mut explicit_zero_cfg = make_cfg(out_zero.path());
        explicit_zero_cfg.set_min_mapq(0);

        // Act
        run_gc_bias(&default_cfg)?;
        run_gc_bias(&explicit_thirty_cfg)?;
        run_gc_bias(&explicit_zero_cfg)?;

        // Assert
        let read_avg_counts = |dir: &TempDir| -> Result<ndarray::Array2<f64>> {
            read_npy(dir.path().join("gc_bias.avg_cfdna_counts.0.npy")).map_err(Into::into)
        };

        let default_counts = read_avg_counts(&out_default)?;
        let explicit_thirty_counts = read_avg_counts(&out_thirty)?;
        let explicit_zero_counts = read_avg_counts(&out_zero)?;

        assert_eq!(default_counts.dim(), (1, 101));
        assert_eq!(default_counts, explicit_thirty_counts);
        assert!((default_counts[(0, 50)] - 2.0).abs() < 1e-12);
        assert!((default_counts.sum() - 2.0).abs() < 1e-12);

        assert_eq!(explicit_zero_counts.dim(), (1, 101));
        assert!((explicit_zero_counts[(0, 50)] - 3.0).abs() < 1e-12);
        assert!((explicit_zero_counts.sum() - 3.0).abs() < 1e-12);

        Ok(())
    }

    #[test]
    fn default_windows_match_explicit_by_size_and_differ_from_global() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Build a two-window genome with the command's default GC window size.
        // - chr1[0,100000) is all A, so any 10 bp fragment there has GC=0 -> GC%=0.
        // - chr1[100000,200000) is all C, so any 10 bp fragment there has GC=10 -> GC%=100.
        //
        // Sample fragments:
        // - one fragment in the left window
        // - nine fragments in the right window
        //
        // This makes the default `by-size 100000` path scientifically different from `--global`:
        // - windowed mode scales each window by its own mean count before averaging windows
        // - global mode keeps one combined raw count table
        //
        // For length 10 with end_offset 0 there are 11 reachable GC-count cells (0..=10), so:
        //
        // Windowed/default path:
        // - left window raw counts: one fragment at gc=0
        //   mean_count = 1 / 11
        //   scale      = 1 / (1/11) * (100000 / 100000) = 11
        //   scaled row = 11 at GC%=0
        // - right window raw counts: nine fragments at gc=10
        //   mean_count = 9 / 11
        //   scale      = 1 / (9/11) * (100000 / 100000) = 11/9
        //   scaled row = 11 at GC%=100
        // - average across two windows -> 5.5 at GC%=0 and 5.5 at GC%=100
        //
        // Global path:
        // - `gc-bias` does not apply per-window scaling in global mode.
        // - There is one combined raw count table with:
        //     1 count at GC%=0
        //     9 counts at GC%=100
        // - So the saved `avg_cfdna_counts` array is exactly those raw counts.
        //
        // The width correction is a no-op here because length 10 maps exactly to GC% bins
        // {0,10,20,...,100}, each with width 1.
        let reference = fixtures::twobit_from_sequences(
            "gc_bias_two_window_reference",
            vec![(
                "chr1".to_string(),
                format!("{}{}", "A".repeat(100_000), "C".repeat(100_000)),
            )],
        )?;
        let starts = [
            100_i64, 100_100, 100_120, 100_140, 100_160, 100_180, 100_200, 100_220, 100_240,
            100_260,
        ];
        let fragments = starts
            .into_iter()
            .map(|start| fixtures::paired_fragment(start, 10, 5))
            .collect();
        let bam = fixtures::bam_from_specs(
            vec![("chr1".to_string(), 200_000)],
            fragments,
            Vec::new(),
            "gc_bias_two_window_bam",
        )?;
        let ref_gc_dir = TempDir::new()?;
        write_reference_package_for_single_length(&reference.path, &ref_gc_dir, 10, 0)?;

        let default_out = TempDir::new()?;
        let explicit_out = TempDir::new()?;
        let global_out = TempDir::new()?;

        let default_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            default_out.path(),
        );

        let mut explicit_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            explicit_out.path(),
        );
        explicit_cfg.set_windows(GCWindowsArgs {
            by_size: Some(100_000),
            by_bed: None,
            global: false,
        });

        let mut global_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            global_out.path(),
        );
        global_cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });

        // Act
        run_gc_bias(&default_cfg)?;
        run_gc_bias(&explicit_cfg)?;
        run_gc_bias(&global_cfg)?;

        let default_counts: ndarray::Array2<f64> =
            read_npy(default_out.path().join("gc_bias.avg_cfdna_counts.0.npy"))?;
        let explicit_counts: ndarray::Array2<f64> =
            read_npy(explicit_out.path().join("gc_bias.avg_cfdna_counts.0.npy"))?;
        let global_counts: ndarray::Array2<f64> =
            read_npy(global_out.path().join("gc_bias.avg_cfdna_counts.0.npy"))?;

        // Assert:
        // The implicit default must be exactly the same behavior as explicit `--by-size 100000`.
        assert_eq!(default_counts, explicit_counts);
        assert_eq!(default_counts.dim(), (1, 101));
        assert_eq!(global_counts.dim(), (1, 101));

        for (gc_pct, &value) in default_counts.row(0).iter().enumerate() {
            match gc_pct {
                0 | 100 => assert!(
                    (value - 5.5).abs() < 1e-12,
                    "default by-size expected 5.5 at GC% {gc_pct}, got {value}"
                ),
                _ => assert!(
                    value.abs() < 1e-12,
                    "default by-size expected 0 outside GC% 0/100, got bin {gc_pct}={value}"
                ),
            }
        }

        for (gc_pct, &value) in global_counts.row(0).iter().enumerate() {
            match gc_pct {
                0 => assert!(
                    (value - 1.0).abs() < 1e-12,
                    "global expected 1.0 at GC% 0, got {value}"
                ),
                100 => assert!(
                    (value - 9.0).abs() < 1e-12,
                    "global expected 9.0 at GC% 100, got {value}"
                ),
                _ => assert!(
                    value.abs() < 1e-12,
                    "global expected 0 outside GC% 0/100, got bin {gc_pct}={value}"
                ),
            }
        }

        assert_ne!(
            default_counts, global_counts,
            "default gc-bias windowing should not collapse to the global calculation here"
        );

        Ok(())
    }

    #[test]
    fn errors_when_blacklist_removes_all_usable_gc_support() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // The blacklist masks the entire chromosome to N for gc-bias counting.
        // The command then cannot compute GC for any fragment, so no window contributes counts.
        // The scientifically correct outcome is a hard error rather than silently writing a
        // degenerate correction package.
        let bam = fixtures::simple_inward_bam()?;
        let reference = fixtures::simple_reference_twobit()?;
        let ref_gc_dir = TempDir::new()?;
        write_reference_package_for_single_length(&reference.path, &ref_gc_dir, 60, 0)?;

        let baseline_out_dir = TempDir::new()?;
        let baseline_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            baseline_out_dir.path(),
        );
        run_gc_bias(&baseline_cfg)?;
        let baseline_package =
            GCCorrectionPackage::from_file(baseline_out_dir.path().join("gc_bias_correction.npz"))?;
        assert_eq!(baseline_package.correction_matrix.dim(), (1, 1));
        assert!((baseline_package.correction_matrix[(0, 0)] - 1.0).abs() < 1e-12);

        let out_dir = TempDir::new()?;
        let blacklist_path = out_dir.path().join("blacklist.bed");
        fixtures::write_bed(&blacklist_path, &[("chr1", 0, 256, "all_masked")])?;

        let mut cfg =
            make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());
        cfg.set_blacklist(Some(vec![blacklist_path]));

        // Act
        let err = run_gc_bias(&cfg).expect_err("fully masked input should fail");

        // Assert:
        // All fragments lose their usable ACGT support after masking, so `scaled_weight` stays 0
        // and the command must fail at the explicit guardrail in `run()`.
        let msg = format!("{err:#}");
        assert!(
            msg.contains("No usable GC bias windows produced counts"),
            "unexpected error message: {msg}"
        );

        Ok(())
    }

    #[test]
    fn correction_package_propagates_reference_end_offset_for_single_length() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Use a reference package and cfDNA run that both allow exactly one fragment length (60 bp)
        // and trim 2 bp from each end when computing GC.
        //
        // The correction package should therefore:
        // - preserve the schema version
        // - preserve `end_offset = 2`
        // - keep exactly one length bin, whose inclusive edge encoding is [60, 60]
        // - assign all length-bin frequency mass to that single bin
        let bam = fixtures::simple_inward_bam()?;
        let reference = fixtures::simple_reference_twobit()?;
        let ref_gc_dir = TempDir::new()?;
        write_reference_package_for_single_length(&reference.path, &ref_gc_dir, 60, 2)?;

        let out_dir = TempDir::new()?;
        let mut cfg =
            make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());
        cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });

        // Act
        run_gc_bias(&cfg)?;
        let package =
            GCCorrectionPackage::from_file(out_dir.path().join("gc_bias_correction.npz"))?;

        // Assert:
        // The reference is `ACGT` repeated and the only cfDNA fragment spans [20,80). With
        // `end_offset = 2`, every counted 56 bp span still contains exactly 28 GC bases:
        //   GC% = 28 / 56 = 50
        // The reference-side counts are built from the same repeated sequence and same effective
        // length, so all mass also lands at GC% 50. With the default 1% GC-bin mass threshold,
        // greedy GC binning therefore collapses the whole GC axis into one bin [0,100], and the
        // correction factor in that only cell must be exactly neutral (1.0).
        assert_eq!(package.version, GC_CORRECTION_SCHEMA_VERSION);
        assert_eq!(package.end_offset, 2);
        assert_eq!(package.length_edges, vec![60, 60]);
        assert_eq!(package.length_bin_frequencies.len(), 1);
        assert!((package.length_bin_frequencies[0] - 1.0).abs() < 1e-12);
        assert_eq!(package.correction_matrix.nrows(), 1);
        assert_eq!(package.correction_matrix.ncols(), 1);
        assert_eq!(package.gc_edges, vec![0, 100]);
        assert!((package.correction_matrix[(0, 0)] - 1.0).abs() < 1e-12);

        Ok(())
    }

    #[test]
    fn overlapping_and_touching_bed_windows_does_not_match_explicitly_merged_gc_bias_run()
    -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Build one chromosome with a sharp A->C transition:
        // - chr1[0,100000)   = all A
        // - chr1[100000,250000) = all C
        //
        // Compare two BED descriptions:
        // - split BED:
        //     [0,100000)     -- A-only
        //     [50000,150000) -- overlaps the first window and touches the third at 150000
        //     [150000,250000) -- C-only
        // - merged BED:
        //     [0,250000)
        //
        // Sample fragments:
        // - one A-only fragment at 60000
        // - nine C-only fragments in the third window at distinct starts
        //
        // With exact BED-window semantics:
        // - split window 1 -> scaled row 11 at GC%=0
        // - split window 2 -> scaled row 11 at GC%=0
        // - split window 3 -> scaled row 11 at GC%=100
        // - average across the three windows:
        //     GC% 0   -> 22/3
        //     GC% 100 -> 11/3
        //
        // With merged-window semantics:
        // - one window with raw row [1, 0, ..., 0, 9]
        // - for length 10 there are 11 reachable GC-count states, so the mean is 10/11
        // - scaled row becomes:
        //     GC% 0   -> 1.1
        //     GC% 100 -> 9.9
        //
        // Therefore the split and merged runs must differ. If `gc-bias` ever starts merging BED
        // windows again, the split run will collapse to the merged result and this test will fail.
        let reference = fixtures::twobit_from_sequences(
            "gc_bias_overlapping_touching_vs_merged_reference",
            vec![(
                "chr1".to_string(),
                format!("{}{}", "A".repeat(100_000), "C".repeat(150_000)),
            )],
        )?;
        let bam = fixtures::bam_from_specs(
            vec![("chr1".to_string(), 250_000)],
            {
                let mut fragments = vec![fixtures::paired_fragment(60_000, 10, 5)];
                for fragment_idx in 0..9 {
                    fragments.push(fixtures::paired_fragment(
                        160_000 + (fragment_idx as i64 * 20),
                        10,
                        5,
                    ));
                }
                fragments
            },
            Vec::new(),
            "gc_bias_overlapping_touching_vs_merged_bam",
        )?;
        let ref_gc_dir = TempDir::new()?;
        write_reference_package_for_single_length(&reference.path, &ref_gc_dir, 10, 0)?;

        let split_out = TempDir::new()?;
        let merged_out = TempDir::new()?;
        let split_bed = split_out.path().join("split_windows.bed");
        let merged_bed = merged_out.path().join("merged_windows.bed");
        std::fs::write(
            &split_bed,
            "chr1\t0\t100000\nchr1\t50000\t150000\nchr1\t150000\t250000\n",
        )?;
        std::fs::write(&merged_bed, "chr1\t0\t250000\n")?;

        let mut split_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            split_out.path(),
        );
        split_cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: Some(split_bed),
            global: false,
        });
        split_cfg.set_min_gc_bin_mass(1.0);
        split_cfg.set_min_length_bin_mass(0.0);
        split_cfg.set_min_length_bin_width(1);
        // Disable package-stage support masking so the split-vs-merged difference survives all
        // the way into the written correction matrix.
        split_cfg.set_num_extreme_gc_bins(0);
        split_cfg.set_num_short_length_bins(0);

        let mut merged_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            merged_out.path(),
        );
        merged_cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: Some(merged_bed),
            global: false,
        });
        merged_cfg.set_min_gc_bin_mass(1.0);
        merged_cfg.set_min_length_bin_mass(0.0);
        merged_cfg.set_min_length_bin_width(1);
        merged_cfg.set_num_extreme_gc_bins(0);
        merged_cfg.set_num_short_length_bins(0);

        // Act
        run_gc_bias(&split_cfg)?;
        run_gc_bias(&merged_cfg)?;

        // Assert
        let split_avg: ndarray::Array2<f64> =
            read_npy(split_out.path().join("gc_bias.avg_cfdna_counts.0.npy"))?;
        let merged_avg: ndarray::Array2<f64> =
            read_npy(merged_out.path().join("gc_bias.avg_cfdna_counts.0.npy"))?;
        assert_eq!(split_avg.dim(), (1, 101));
        assert_eq!(merged_avg.dim(), (1, 101));

        for (gc_pct, &value) in split_avg.row(0).iter().enumerate() {
            match gc_pct {
                0 => assert!(
                    (value - (22.0 / 3.0)).abs() < 1e-12,
                    "expected 22/3 at GC% 0 in split-window mode, got {value}"
                ),
                100 => assert!(
                    (value - (11.0 / 3.0)).abs() < 1e-12,
                    "expected 11/3 at GC% 100 in split-window mode, got {value}"
                ),
                _ => assert!(
                    value.abs() < 1e-12,
                    "expected 0 outside GC% 0/100 in split-window mode, got bin {gc_pct}={value}"
                ),
            }
        }
        for (gc_pct, &value) in merged_avg.row(0).iter().enumerate() {
            match gc_pct {
                0 => assert!(
                    (value - 1.1).abs() < 1e-12,
                    "expected 1.1 at GC% 0 in merged-window mode, got {value}"
                ),
                100 => assert!(
                    (value - 9.9).abs() < 1e-12,
                    "expected 9.9 at GC% 100 in merged-window mode, got {value}"
                ),
                _ => assert!(
                    value.abs() < 1e-12,
                    "expected 0 outside GC% 0/100 in merged-window mode, got bin {gc_pct}={value}"
                ),
            }
        }

        let split_pkg =
            GCCorrectionPackage::from_file(split_out.path().join("gc_bias_correction.npz"))?;
        let merged_pkg =
            GCCorrectionPackage::from_file(merged_out.path().join("gc_bias_correction.npz"))?;
        assert_eq!(split_pkg.version, merged_pkg.version);
        assert_eq!(split_pkg.end_offset, merged_pkg.end_offset);
        assert_eq!(split_pkg.length_edges, merged_pkg.length_edges);
        assert_eq!(split_pkg.gc_edges, merged_pkg.gc_edges);
        assert_eq!(
            split_pkg.length_bin_frequencies,
            merged_pkg.length_bin_frequencies
        );
        assert_ne!(split_pkg.correction_matrix, merged_pkg.correction_matrix);

        Ok(())
    }

    #[test]
    fn by_size_gc_bias_is_invariant_to_aligned_vs_misaligned_tile_sizes() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Build a 2 Mb two-state genome in 100 kb windows:
        // - windows with even index are all A  -> GC%=0 for 10 bp fragments
        // - windows with odd  index are all C  -> GC%=100 for 10 bp fragments
        //
        // Place one 10 bp fragment in each 100 kb window, 20 windows total.
        //
        // This means the logical GC-bias result is completely determined by the windowing, not by
        // the tile cuts:
        // - 10 windows contribute one fragment at GC%=0
        // - 10 windows contribute one fragment at GC%=100
        // - every per-window scaled row is therefore identical within its GC state
        //
        // We then run the same explicit `--by-size 100000` command with:
        // - tile_size = 1,000,000  -> aligned to the 100 kb window grid
        // - tile_size =   950,000  -> misaligned (ratio < 10, so no alignment correction)
        //
        // The chromosome is long enough to span multiple tiles in both cases, so both the aligned
        // and the general cross-tile reducer paths are real. Since tiling is an execution detail,
        // the final averaged cfDNA counts and correction package must be identical.
        let window_bp = 100_000usize;
        let num_windows = 20usize;
        let mut sequence = String::with_capacity(window_bp * num_windows);
        let mut fragments = Vec::with_capacity(num_windows);
        for window_idx in 0..num_windows {
            let base = if window_idx % 2 == 0 { 'A' } else { 'C' };
            sequence.push_str(&base.to_string().repeat(window_bp));
            let start = (window_idx * window_bp + 100) as i64;
            fragments.push(fixtures::paired_fragment(start, 10, 5));
        }
        let chrom_len = (window_bp * num_windows) as u32;
        let reference = fixtures::twobit_from_sequences(
            "gc_bias_tile_invariance_reference",
            vec![("chr1".to_string(), sequence)],
        )?;
        let bam = fixtures::bam_from_specs(
            vec![("chr1".to_string(), chrom_len)],
            fragments,
            Vec::new(),
            "gc_bias_tile_invariance_bam",
        )?;
        let ref_gc_dir = TempDir::new()?;
        write_reference_package_for_single_length(&reference.path, &ref_gc_dir, 10, 0)?;

        let aligned_out = TempDir::new()?;
        let misaligned_out = TempDir::new()?;

        let mut aligned_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            aligned_out.path(),
        );
        aligned_cfg.set_windows(GCWindowsArgs {
            by_size: Some(100_000),
            by_bed: None,
            global: false,
        });
        aligned_cfg.set_tile_size(1_000_000);

        let mut misaligned_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            misaligned_out.path(),
        );
        misaligned_cfg.set_windows(GCWindowsArgs {
            by_size: Some(100_000),
            by_bed: None,
            global: false,
        });
        misaligned_cfg.set_tile_size(950_000);

        // Act
        run_gc_bias(&aligned_cfg)?;
        run_gc_bias(&misaligned_cfg)?;

        // Assert
        let aligned_avg: ndarray::Array2<f64> =
            read_npy(aligned_out.path().join("gc_bias.avg_cfdna_counts.0.npy"))?;
        let misaligned_avg: ndarray::Array2<f64> =
            read_npy(misaligned_out.path().join("gc_bias.avg_cfdna_counts.0.npy"))?;
        assert_eq!(aligned_avg, misaligned_avg);
        assert_eq!(aligned_avg.dim(), (1, 101));
        for (gc_pct, &value) in aligned_avg.row(0).iter().enumerate() {
            match gc_pct {
                0 | 100 => assert!(
                    (value - 5.5).abs() < 1e-12,
                    "expected 5.5 at GC% {gc_pct}, got {value}"
                ),
                _ => assert!(
                    value.abs() < 1e-12,
                    "expected 0 outside GC% 0/100, got bin {gc_pct}={value}"
                ),
            }
        }

        let aligned_pkg =
            GCCorrectionPackage::from_file(aligned_out.path().join("gc_bias_correction.npz"))?;
        let misaligned_pkg =
            GCCorrectionPackage::from_file(misaligned_out.path().join("gc_bias_correction.npz"))?;
        assert_eq!(aligned_pkg.version, misaligned_pkg.version);
        assert_eq!(aligned_pkg.end_offset, misaligned_pkg.end_offset);
        assert_eq!(aligned_pkg.length_edges, misaligned_pkg.length_edges);
        assert_eq!(aligned_pkg.gc_edges, misaligned_pkg.gc_edges);
        assert_eq!(
            aligned_pkg.length_bin_frequencies,
            misaligned_pkg.length_bin_frequencies
        );
        assert_eq!(
            aligned_pkg.correction_matrix,
            misaligned_pkg.correction_matrix
        );
        assert_eq!(aligned_pkg.version, GC_CORRECTION_SCHEMA_VERSION);
        assert_eq!(aligned_pkg.end_offset, 0);
        assert_eq!(aligned_pkg.length_edges, vec![10, 10]);
        assert_eq!(aligned_pkg.gc_edges, vec![0, 1, 100]);
        assert_eq!(aligned_pkg.length_bin_frequencies.len(), 1);
        assert!((aligned_pkg.length_bin_frequencies[0] - 1.0).abs() < 1e-12);
        assert_eq!(aligned_pkg.correction_matrix.dim(), (1, 2));
        assert!((aligned_pkg.correction_matrix[(0, 0)] - 1.0).abs() < 1e-12);
        assert!((aligned_pkg.correction_matrix[(0, 1)] - 1.0).abs() < 1e-12);

        Ok(())
    }

    #[test]
    fn multi_chromosome_by_size_gc_bias_accumulates_windows_across_chromosomes() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Build two 100 bp chromosomes and count one 10 bp fragment in one fixed-size window on
        // each chromosome:
        // - chr1 = all A, so the fragment contributes only to GC% 0
        // - chr2 = all C, so the fragment contributes only to GC% 100
        //
        // Sample-side window scaling in `gc-bias` is per counted window:
        // - a 10 bp fragment has 11 reachable GC-count states (0..10), so one observed count in a
        //   pure window is scaled by 11 after mean normalization
        // - each chromosome contributes exactly one counted 100 bp window
        // - the saved `avg_cfdna_counts` is therefore the mean of two scaled windows:
        //     GC% 0   -> 11 / 2 = 5.5
        //     GC% 100 -> 11 / 2 = 5.5
        //     all other GC bins -> 0
        //
        // This is the smallest real multi-chromosome `gc-bias::run()` fixture that proves:
        // - per-chromosome fixed-size windows are all visited
        // - counts accumulate across chromosomes before the final averaging
        // - the default single-length package stage neutralizes to the exact 1x2 row [1, 1]
        //
        // That package derivation is deterministic here:
        // - greedy GC binning on the saved sample counts closes bins at [0] and [1..100]
        // - there is only one length bin, so the default `num_short_length_bins = 1` masks the
        //   entire row
        // - on a 2-bin GC axis, the default `num_extreme_gc_bins = 1` also masks both columns
        // - the package pipeline therefore falls back to the neutral multiplicative row [1, 1]
        let reference = fixtures::twobit_from_sequences(
            "gc_bias_multi_chr_reference",
            vec![
                ("chr1".to_string(), "A".repeat(100)),
                ("chr2".to_string(), "C".repeat(100)),
            ],
        )?;
        let bam = fixtures::bam_from_specs(
            vec![("chr1".to_string(), 100), ("chr2".to_string(), 100)],
            vec![fixtures::paired_fragment(10, 10, 5), {
                let mut fragment = fixtures::paired_fragment(10, 10, 5);
                fragment.forward.tid = 1;
                fragment.reverse.tid = 1;
                fragment.forward.mate_tid = Some(1);
                fragment.reverse.mate_tid = Some(1);
                fragment
            }],
            Vec::new(),
            "gc_bias_multi_chr_bam",
        )?;
        let ref_gc_dir = TempDir::new()?;
        let ref_cfg = RefGCBiasConfig {
            ref_genome: Ref2BitRequiredArgs {
                ref_2bit: reference.path.clone(),
            },
            output_dir: ref_gc_dir.path().to_path_buf(),
            output_prefix: String::new(),
            n_threads: 1,
            // Two chromosomes of length 100 with fragment length 10 give:
            //   (100 - 10 + 1) * 2 = 91 * 2 = 182 valid starts.
            // Sample all of them so the reference package is deterministic and balanced.
            n_positions: 182,
            seed: Some(7),
            windows: Default::default(),
            chromosomes: ChromosomeArgs {
                chromosomes: Some(vec!["chr1".to_string(), "chr2".to_string()]),
                chromosomes_file: None,
            },
            blacklist: None,
            fragment_lengths: cfdnalab::commands::cli_common::FragmentLengthArgs {
                min_fragment_length: 10,
                max_fragment_length: 10,
            },
            end_offset: 0,
            skip_interpolation: true,
            smoothing_sigma: 0.55,
            smoothing_radius: 2,
            skip_smoothing: true,
            tile_size: 1_000_000,
            logging: LoggingArgs::default(),
        };
        run_ref_gc_bias(&ref_cfg)?;

        let out_dir = TempDir::new()?;
        let ioc = IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.path().to_path_buf(),
            n_threads: 1,
        };
        let mut cfg = GCConfig::new(
            ioc,
            reference.path.clone(),
            ref_gc_dir.path().join("ref_gc_package.npz"),
            ChromosomeArgs {
                chromosomes: Some(vec!["chr1".to_string(), "chr2".to_string()]),
                chromosomes_file: None,
            },
        );
        cfg.set_windows(GCWindowsArgs {
            by_size: Some(100),
            by_bed: None,
            global: false,
        });
        cfg.set_min_mapq(0);
        cfg.set_tile_size(1_000_000);
        cfg.set_min_window_acgt_pct(0);
        cfg.set_save_intermediates(true);

        // Act
        run_gc_bias(&cfg)?;

        // Assert
        let avg_counts: ndarray::Array2<f64> =
            read_npy(out_dir.path().join("gc_bias.avg_cfdna_counts.0.npy"))?;
        assert_eq!(avg_counts.dim(), (1, 101));
        for (gc_pct, &value) in avg_counts.row(0).iter().enumerate() {
            match gc_pct {
                0 | 100 => assert!(
                    (value - 5.5).abs() < 1e-12,
                    "expected 5.5 at GC% {gc_pct}, got {value}"
                ),
                _ => assert!(
                    value.abs() < 1e-12,
                    "expected no mass outside GC% 0/100, got bin {gc_pct}={value}"
                ),
            }
        }

        let package =
            GCCorrectionPackage::from_file(out_dir.path().join("gc_bias_correction.npz"))?;
        assert_eq!(package.version, GC_CORRECTION_SCHEMA_VERSION);
        assert_eq!(package.end_offset, 0);
        assert_eq!(package.length_edges, vec![10, 10]);
        assert_eq!(package.gc_edges, vec![0, 1, 100]);
        assert_eq!(package.length_bin_frequencies.len(), 1);
        assert!((package.length_bin_frequencies[0] - 1.0).abs() < 1e-12);
        assert_eq!(package.correction_matrix.dim(), (1, 2));
        assert!((package.correction_matrix[(0, 0)] - 1.0).abs() < 1e-12);
        assert!((package.correction_matrix[(0, 1)] - 1.0).abs() < 1e-12);

        Ok(())
    }

    #[test]
    fn empty_middle_tile_matches_single_tile_gc_bias_run() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Use three logical 100 bp windows on one chromosome:
        // - [0,100)   = all A, contains one 10 bp fragment -> GC% 0
        // - [100,200) = all A, contains no fragments       -> empty counted window
        // - [200,300) = all C, contains one 10 bp fragment -> GC% 100
        //
        // For every counted pure window, sample-side window scaling gives 11 at its observed GC
        // bin, exactly as in the two-window tests. The middle window has zero counts, so
        // `gc-bias` drops it before the global averaging step. The saved `avg_cfdna_counts` must
        // therefore be:
        //   GC% 0   -> 11 / 2 = 5.5
        //   GC% 100 -> 11 / 2 = 5.5
        //   all other bins -> 0
        //
        // The point of this test is the execution layout:
        // - `tile_size = 95` forces a four-tile run where tile 1 has no fragments at all
        // - `tile_size = 1_000` keeps the whole chromosome in one tile
        // Both runs must produce the same saved counts and the same default single-length package.
        //
        // As in the multi-chromosome test above, the exact package is still hand-derivable:
        // greedy GC binning creates [0] and [1..100], and the default single-length masking then
        // neutralizes the only row to multiplicative weights [1, 1].
        let reference = fixtures::twobit_from_sequences(
            "gc_bias_empty_tile_reference",
            vec![(
                "chr1".to_string(),
                format!("{}{}{}", "A".repeat(100), "A".repeat(100), "C".repeat(100)),
            )],
        )?;
        let bam = fixtures::bam_from_specs(
            vec![("chr1".to_string(), 300)],
            vec![
                fixtures::paired_fragment(10, 10, 5),
                fixtures::paired_fragment(210, 10, 5),
            ],
            Vec::new(),
            "gc_bias_empty_tile_bam",
        )?;
        let ref_gc_dir = TempDir::new()?;
        write_reference_package_for_single_length(&reference.path, &ref_gc_dir, 10, 0)?;

        let multi_tile_out = TempDir::new()?;
        let single_tile_out = TempDir::new()?;

        let mut multi_tile_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            multi_tile_out.path(),
        );
        multi_tile_cfg.set_windows(GCWindowsArgs {
            by_size: Some(100),
            by_bed: None,
            global: false,
        });
        multi_tile_cfg.set_tile_size(95);

        let mut single_tile_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            single_tile_out.path(),
        );
        single_tile_cfg.set_windows(GCWindowsArgs {
            by_size: Some(100),
            by_bed: None,
            global: false,
        });
        single_tile_cfg.set_tile_size(1_000);

        // Act
        run_gc_bias(&multi_tile_cfg)?;
        run_gc_bias(&single_tile_cfg)?;

        // Assert
        let multi_tile_avg: ndarray::Array2<f64> =
            read_npy(multi_tile_out.path().join("gc_bias.avg_cfdna_counts.0.npy"))?;
        let single_tile_avg: ndarray::Array2<f64> = read_npy(
            single_tile_out
                .path()
                .join("gc_bias.avg_cfdna_counts.0.npy"),
        )?;
        assert_eq!(multi_tile_avg, single_tile_avg);
        assert_eq!(multi_tile_avg.dim(), (1, 101));
        for (gc_pct, &value) in multi_tile_avg.row(0).iter().enumerate() {
            match gc_pct {
                0 | 100 => assert!(
                    (value - 5.5).abs() < 1e-12,
                    "expected 5.5 at GC% {gc_pct}, got {value}"
                ),
                _ => assert!(
                    value.abs() < 1e-12,
                    "expected no mass outside GC% 0/100, got bin {gc_pct}={value}"
                ),
            }
        }

        let multi_tile_pkg =
            GCCorrectionPackage::from_file(multi_tile_out.path().join("gc_bias_correction.npz"))?;
        let single_tile_pkg =
            GCCorrectionPackage::from_file(single_tile_out.path().join("gc_bias_correction.npz"))?;
        assert_eq!(multi_tile_pkg.version, single_tile_pkg.version);
        assert_eq!(multi_tile_pkg.end_offset, single_tile_pkg.end_offset);
        assert_eq!(multi_tile_pkg.length_edges, single_tile_pkg.length_edges);
        assert_eq!(multi_tile_pkg.gc_edges, single_tile_pkg.gc_edges);
        assert_eq!(
            multi_tile_pkg.length_bin_frequencies,
            single_tile_pkg.length_bin_frequencies
        );
        assert_eq!(
            multi_tile_pkg.correction_matrix,
            single_tile_pkg.correction_matrix
        );
        assert_eq!(multi_tile_pkg.version, GC_CORRECTION_SCHEMA_VERSION);
        assert_eq!(multi_tile_pkg.end_offset, 0);
        assert_eq!(multi_tile_pkg.length_edges, vec![10, 10]);
        assert_eq!(multi_tile_pkg.gc_edges, vec![0, 1, 100]);
        assert_eq!(multi_tile_pkg.length_bin_frequencies.len(), 1);
        assert!((multi_tile_pkg.length_bin_frequencies[0] - 1.0).abs() < 1e-12);
        assert_eq!(multi_tile_pkg.correction_matrix.dim(), (1, 2));
        assert!((multi_tile_pkg.correction_matrix[(0, 0)] - 1.0).abs() < 1e-12);
        assert!((multi_tile_pkg.correction_matrix[(0, 1)] - 1.0).abs() < 1e-12);

        Ok(())
    }

    #[test]
    fn touching_bed_windows_match_by_size_counts_and_default_single_length_packages() -> Result<()>
    {
        // Human verification status: unverified
        // Arrange:
        // Use the same two-window A/C genome as the default-window test, but compare two
        // different *window representations* of the same logical partition:
        //   by-size 100000
        //   by-bed  [0,100000), [100000,200000)
        //
        // The chromosome is split at the same 100 kb boundary in both runs, and the fragment
        // placement is identical:
        // - one 10 bp fragment in the left A-only window -> GC%=0
        // - nine 10 bp fragments in the right C-only window -> GC%=100
        //
        // With window scaling enabled, the expected per-window scaled rows are exactly the same
        // as in `default_windows_match_explicit_by_size_and_differ_from_global`:
        // - left window  -> scaled count 11 at GC%=0
        // - right window -> scaled count 11 at GC%=100
        // - averaged over the two windows -> 5.5 at GC%=0 and 5.5 at GC%=100
        //
        // We deliberately use `tile_size = 95_000` so neither 100 kb window is tile-contained.
        // That makes the fixed-size streaming path and the explicit BED-window path both exercise
        // real cross-tile counting/reduction instead of a degenerate one-tile case.
        let reference = fixtures::twobit_from_sequences(
            "gc_bias_by_size_vs_bed_reference",
            vec![(
                "chr1".to_string(),
                format!("{}{}", "A".repeat(100_000), "C".repeat(100_000)),
            )],
        )?;
        let starts = [
            100_i64, 100_100, 100_120, 100_140, 100_160, 100_180, 100_200, 100_220, 100_240,
            100_260,
        ];
        let fragments = starts
            .into_iter()
            .map(|start| fixtures::paired_fragment(start, 10, 5))
            .collect();
        let bam = fixtures::bam_from_specs(
            vec![("chr1".to_string(), 200_000)],
            fragments,
            Vec::new(),
            "gc_bias_by_size_vs_bed_bam",
        )?;
        let ref_gc_dir = TempDir::new()?;
        write_reference_package_for_single_length(&reference.path, &ref_gc_dir, 10, 0)?;

        let by_size_out = TempDir::new()?;
        let by_bed_out = TempDir::new()?;
        let bed_path = by_bed_out.path().join("windows.bed");
        std::fs::write(&bed_path, "chr1\t0\t100000\nchr1\t100000\t200000\n")?;

        let mut by_size_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            by_size_out.path(),
        );
        by_size_cfg.set_windows(GCWindowsArgs {
            by_size: Some(100_000),
            by_bed: None,
            global: false,
        });
        by_size_cfg.set_tile_size(95_000);

        let mut by_bed_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            by_bed_out.path(),
        );
        by_bed_cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: Some(bed_path),
            global: false,
        });
        by_bed_cfg.set_tile_size(95_000);

        // Act
        run_gc_bias(&by_size_cfg)?;
        run_gc_bias(&by_bed_cfg)?;

        // Assert:
        // `gc-bias` does not merge touching BED windows on load. These two runs therefore encode
        // the same logical partition and must produce the same averaged counts and the same
        // default single-length neutral package.
        let by_size_avg: ndarray::Array2<f64> =
            read_npy(by_size_out.path().join("gc_bias.avg_cfdna_counts.0.npy"))?;
        let by_bed_avg: ndarray::Array2<f64> =
            read_npy(by_bed_out.path().join("gc_bias.avg_cfdna_counts.0.npy"))?;
        assert_eq!(by_size_avg, by_bed_avg);
        assert_eq!(by_size_avg.dim(), (1, 101));
        for (gc_pct, &value) in by_size_avg.row(0).iter().enumerate() {
            match gc_pct {
                0 | 100 => assert!(
                    (value - 5.5).abs() < 1e-12,
                    "expected 5.5 at GC% {gc_pct}, got {value}"
                ),
                _ => assert!(
                    value.abs() < 1e-12,
                    "expected 0 outside GC% 0/100, got bin {gc_pct}={value}"
                ),
            }
        }

        let by_size_pkg =
            GCCorrectionPackage::from_file(by_size_out.path().join("gc_bias_correction.npz"))?;
        let by_bed_pkg =
            GCCorrectionPackage::from_file(by_bed_out.path().join("gc_bias_correction.npz"))?;
        assert_eq!(by_size_pkg.version, by_bed_pkg.version);
        assert_eq!(by_size_pkg.end_offset, by_bed_pkg.end_offset);
        assert_eq!(by_size_pkg.length_edges, by_bed_pkg.length_edges);
        assert_eq!(by_size_pkg.gc_edges, by_bed_pkg.gc_edges);
        assert_eq!(
            by_size_pkg.length_bin_frequencies,
            by_bed_pkg.length_bin_frequencies
        );
        assert_eq!(by_size_pkg.correction_matrix, by_bed_pkg.correction_matrix);
        assert_eq!(by_size_pkg.version, GC_CORRECTION_SCHEMA_VERSION);
        assert_eq!(by_size_pkg.end_offset, 0);
        assert_eq!(by_size_pkg.length_edges, vec![10, 10]);
        assert_eq!(by_size_pkg.gc_edges, vec![0, 1, 100]);
        assert_eq!(by_size_pkg.length_bin_frequencies.len(), 1);
        assert!((by_size_pkg.length_bin_frequencies[0] - 1.0).abs() < 1e-12);
        assert_eq!(by_size_pkg.correction_matrix.dim(), (1, 2));
        assert!((by_size_pkg.correction_matrix[(0, 0)] - 1.0).abs() < 1e-12);
        assert!((by_size_pkg.correction_matrix[(0, 1)] - 1.0).abs() < 1e-12);

        Ok(())
    }

    #[test]
    fn real_ref_gc_bias_then_gc_bias_package_is_non_neutral_in_two_bin_case() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Build a reference and cfDNA sample where the real producer->consumer workflow must
        // create a non-neutral correction package with exactly two GC bins.
        //
        // Reference genome:
        // - chr1[0,100)   = all A
        // - chr1[100,200) = all C
        //
        // Reference-GC producer setup:
        // - fragment length is fixed at 10
        // - valid starts are 0..=190, so there are exactly 191 valid start positions
        // - we set `n_positions = 191`, so every valid start is sampled
        // - BED windows keep only the pure-A and pure-C start ranges that also leave enough room
        //   for the full 10 bp fragment inside the same BED interval:
        //     [0,91)    -> starts 0..=81   -> GC%=0
        //     [100,191) -> starts 100..=181 -> GC%=100
        // - starts 82..=99 and 182..=190 are inside a BED row but do not fit before its right
        //   edge, so `ref-gc-bias` excludes them
        //
        // Therefore the written reference counts are exactly balanced:
        // - 82 counts at GC%=0
        // - 82 counts at GC%=100
        //
        // Sample BAM:
        // - one 10 bp fragment in the A-only region  -> GC%=0
        // - nine 10 bp fragments in the C-only region -> GC%=100
        //
        // `gc-bias` is run in global mode with:
        // - no extreme-bin masking
        // - no short-length masking
        // - no outlier handling
        // - default `min_gc_bin_mass = 1.0`
        //
        // Global sample counts are therefore:
        // - raw avg counts: 1 at GC%=0, 9 at GC%=100
        // - reference-supported mean = (1 + 9) / 2 = 5
        // - normalized avg counts = 0.2 at GC%=0, 1.8 at GC%=100
        //
        // GC binning derivation with `min_gc_bin_mass = 1%`:
        // - total normalized mass = 2.0
        // - minimum mass per bin = 0.02
        // - index 0 alone already exceeds 0.02, so the first GC bin is [0]
        // - indices 1..100 then accumulate until index 100 contributes 1.8, so the second GC
        //   bin is [1..100]
        // - therefore the package edges must be [0, 1, 100]
        //
        // Reference-side binned counts are balanced, so after per-length normalization they are:
        // - [1.0, 1.0]
        //
        // cfDNA-side binned counts are:
        // - [0.2, 1.8]
        //
        // Raw correction row:
        // - [0.2 / 1.0, 1.8 / 1.0] = [0.2, 1.8]
        //
        // Its mean is already 1.0, so the final re-centering step changes nothing.
        // The package stores multiplicative correction factors, so the row is inverted:
        // - [1 / 0.2, 1 / 1.8] = [5.0, 5/9]
        let reference = fixtures::twobit_from_sequences(
            "gc_bias_real_non_neutral_reference",
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
        let bam = fixtures::bam_from_specs(
            vec![("chr1".to_string(), 200)],
            fragments,
            Vec::new(),
            "gc_bias_real_non_neutral_bam",
        )?;

        let ref_gc_dir = TempDir::new()?;
        let bed_path = ref_gc_dir.path().join("pure_windows.bed");
        std::fs::write(&bed_path, "chr1\t0\t91\nchr1\t100\t191\n")?;
        let ref_cfg = RefGCBiasConfig {
            ref_genome: Ref2BitRequiredArgs {
                ref_2bit: reference.path.clone(),
            },
            output_dir: ref_gc_dir.path().to_path_buf(),
            output_prefix: String::new(),
            n_threads: 1,
            n_positions: 191,
            seed: Some(23),
            windows: cfdnalab::commands::ref_gc_bias::config::RefGCWindowsArgs {
                by_bed: Some(bed_path),
            },
            chromosomes: ChromosomeArgs {
                chromosomes: Some(vec!["chr1".to_string()]),
                chromosomes_file: None,
            },
            blacklist: None,
            fragment_lengths: cfdnalab::commands::cli_common::FragmentLengthArgs {
                min_fragment_length: 10,
                max_fragment_length: 10,
            },
            end_offset: 0,
            skip_interpolation: true,
            smoothing_sigma: 0.55,
            smoothing_radius: 2,
            skip_smoothing: true,
            tile_size: 1_000_000,
            logging: LoggingArgs::default(),
        };

        let out_dir = TempDir::new()?;
        let mut cfg =
            make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());
        cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });
        cfg.set_num_extreme_gc_bins(0);
        cfg.set_num_short_length_bins(0);
        cfg.set_min_gc_bin_mass(1.0);
        cfg.set_min_length_bin_mass(0.0);
        cfg.set_min_length_bin_width(1);
        cfg.outlier_method = cfdnalab::commands::gc_bias::config::OutlierMethodArg::None;

        // Act
        run_ref_gc_bias(&ref_cfg)?;
        run_gc_bias(&cfg)?;

        // Assert
        let avg_counts: ndarray::Array2<f64> =
            read_npy(out_dir.path().join("gc_bias.avg_cfdna_counts.0.npy"))?;
        assert_eq!(avg_counts.dim(), (1, 101));
        for (gc_pct, &value) in avg_counts.row(0).iter().enumerate() {
            match gc_pct {
                0 => assert!(
                    (value - 1.0).abs() < 1e-12,
                    "expected raw global count 1.0 at GC% 0, got {value}"
                ),
                100 => assert!(
                    (value - 9.0).abs() < 1e-12,
                    "expected raw global count 9.0 at GC% 100, got {value}"
                ),
                _ => assert!(
                    value.abs() < 1e-12,
                    "expected no mass outside GC% 0/100, got bin {gc_pct}={value}"
                ),
            }
        }

        let package =
            GCCorrectionPackage::from_file(out_dir.path().join("gc_bias_correction.npz"))?;
        assert_eq!(package.length_edges, vec![10, 10]);
        assert_eq!(package.gc_edges, vec![0, 1, 100]);
        assert_eq!(package.correction_matrix.dim(), (1, 2));
        assert!((package.correction_matrix[(0, 0)] - 5.0).abs() < 1e-12);
        assert!((package.correction_matrix[(0, 1)] - (5.0 / 9.0)).abs() < 1e-12);
        assert_eq!(package.length_bin_frequencies.len(), 1);
        assert!((package.length_bin_frequencies[0] - 1.0).abs() < 1e-12);

        Ok(())
    }

    #[test]
    fn real_ref_gc_bias_smoothing_metadata_drives_gc_bias_smoothed_avg_counts() -> Result<()> {
        // Arrange:
        // Reuse the three isolated 10 bp windows from the ref-gc-bias smoothing test:
        // - one A-only fragment start  -> GC%=0
        // - one mixed 5C/5A start      -> GC%=50
        // - one C-only fragment start  -> GC%=100
        // - plus two unused trailing A bases so the chromosome length is 52 bp instead of 50,
        //   avoiding the upstream `.2bit` partial-byte tail bug while keeping `[40,50)` intact
        //
        // The sample BAM contains exactly those same three 10 bp fragments, so the pre-smoothing
        // cfDNA counts row is also:
        //   gc_count 0  -> 1
        //   gc_count 5  -> 1
        //   gc_count 10 -> 1
        //
        // The real reference package is produced with:
        //   sigma = sqrt(1 / (2 ln 2)), radius = 1
        // so the 1D smoothing kernel is exactly [1/4, 1/2, 1/4].
        //
        // `gc-bias` must respect the written metadata and smooth its own cfDNA counts in the same
        // way before converting to GC percentages. Therefore the saved `avg_cfdna_counts` matrix
        // must be:
        //   GC% 0   -> 1/2
        //   GC% 10  -> 1/4
        //   GC% 40  -> 1/4
        //   GC% 50  -> 1/2
        //   GC% 60  -> 1/4
        //   GC% 90  -> 1/4
        //   GC% 100 -> 1/2
        // and zero elsewhere.
        let sigma = (1.0_f64 / (2.0 * std::f64::consts::LN_2)).sqrt();
        let reference = fixtures::twobit_from_sequences(
            "gc_bias_real_smoothed_reference",
            vec![(
                "chr1".to_string(),
                format!(
                    "{}{}{}{}{}{}",
                    "A".repeat(10),
                    "T".repeat(10),
                    "C".repeat(5) + &"A".repeat(5),
                    "T".repeat(10),
                    "C".repeat(10),
                    "A".repeat(2)
                ),
            )],
        )?;
        let bam = fixtures::bam_from_specs(
            vec![("chr1".to_string(), 52)],
            vec![
                fixtures::paired_fragment(0, 10, 5),
                fixtures::paired_fragment(20, 10, 5),
                fixtures::paired_fragment(40, 10, 5),
            ],
            Vec::new(),
            "gc_bias_real_smoothed_bam",
        )?;

        let ref_gc_dir = TempDir::new()?;
        let bed_path = ref_gc_dir.path().join("windows.bed");
        std::fs::write(&bed_path, "chr1\t0\t10\nchr1\t20\t30\nchr1\t40\t50\n")?;
        let ref_cfg = RefGCBiasConfig {
            ref_genome: Ref2BitRequiredArgs {
                ref_2bit: reference.path.clone(),
            },
            output_dir: ref_gc_dir.path().to_path_buf(),
            output_prefix: String::new(),
            n_threads: 1,
            n_positions: 41,
            seed: Some(11),
            windows: cfdnalab::commands::ref_gc_bias::config::RefGCWindowsArgs {
                by_bed: Some(bed_path),
            },
            chromosomes: ChromosomeArgs {
                chromosomes: Some(vec!["chr1".to_string()]),
                chromosomes_file: None,
            },
            blacklist: None,
            fragment_lengths: cfdnalab::commands::cli_common::FragmentLengthArgs {
                min_fragment_length: 10,
                max_fragment_length: 10,
            },
            end_offset: 0,
            skip_interpolation: true,
            smoothing_sigma: sigma,
            smoothing_radius: 1,
            skip_smoothing: false,
            tile_size: 1_000_000,
            logging: LoggingArgs::default(),
        };

        let out_dir = TempDir::new()?;
        let mut cfg =
            make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());
        cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });
        cfg.set_num_extreme_gc_bins(0);
        cfg.set_num_short_length_bins(0);
        cfg.outlier_method = cfdnalab::commands::gc_bias::config::OutlierMethodArg::None;

        // Act
        run_ref_gc_bias(&ref_cfg)?;
        run_gc_bias(&cfg)?;

        // Assert
        let avg_counts: ndarray::Array2<f64> =
            read_npy(out_dir.path().join("gc_bias.avg_cfdna_counts.0.npy"))?;
        let expected_non_zero = [
            (0_usize, 0.5_f64),
            (10, 0.25),
            (40, 0.25),
            (50, 0.5),
            (60, 0.25),
            (90, 0.25),
            (100, 0.5),
        ];
        assert_eq!(avg_counts.dim(), (1, 101));
        for (gc_pct, &value) in avg_counts.row(0).iter().enumerate() {
            let expected_value = expected_non_zero
                .iter()
                .find_map(|(expected_gc_pct, expected_value)| {
                    (*expected_gc_pct == gc_pct).then_some(*expected_value)
                })
                .unwrap_or(0.0);
            assert_gc_command_close(
                value,
                expected_value,
                &format!("smoothed avg counts at GC% {gc_pct}"),
            );
        }

        Ok(())
    }

    #[test]
    fn real_ref_gc_bias_interpolation_metadata_drives_gc_bias_interpolated_counts() -> Result<()> {
        // Arrange:
        // Use the same three 10 bp starts as above, but produce the reference package with
        // interpolation enabled and smoothing disabled.
        // The chromosome again has two unused trailing A bases so the critical `[40,50)` interval
        // is preserved without relying on the upstream `.2bit` partial-byte tail behavior.
        //
        // The cfDNA sample again has raw global counts:
        //   GC% 0   -> 1
        //   GC% 50  -> 1
        //   GC% 100 -> 1
        //
        // The real reference package writes `skip_interpolation = false` and an empirical support
        // mask that is true only at GC% 0, 50, and 100. `gc-bias` first mean-scales over those
        // supported cells, so `normalized_avg_cfdna_counts` still has:
        //   1 at GC% 0, 50, 100
        //   0 elsewhere
        //
        // Interpolation then sees three equal anchors and must fill every unsupported GC% bin with
        // the constant value 1.0. The saved `interpolated_cfdna_counts` matrix should therefore be
        // exactly all ones.
        let reference = fixtures::twobit_from_sequences(
            "gc_bias_real_interpolated_reference",
            vec![(
                "chr1".to_string(),
                format!(
                    "{}{}{}{}{}{}",
                    "A".repeat(10),
                    "T".repeat(10),
                    "C".repeat(5) + &"A".repeat(5),
                    "T".repeat(10),
                    "C".repeat(10),
                    "A".repeat(2)
                ),
            )],
        )?;
        let bam = fixtures::bam_from_specs(
            vec![("chr1".to_string(), 52)],
            vec![
                fixtures::paired_fragment(0, 10, 5),
                fixtures::paired_fragment(20, 10, 5),
                fixtures::paired_fragment(40, 10, 5),
            ],
            Vec::new(),
            "gc_bias_real_interpolated_bam",
        )?;

        let ref_gc_dir = TempDir::new()?;
        let bed_path = ref_gc_dir.path().join("windows.bed");
        std::fs::write(&bed_path, "chr1\t0\t10\nchr1\t20\t30\nchr1\t40\t50\n")?;
        let ref_cfg = RefGCBiasConfig {
            ref_genome: Ref2BitRequiredArgs {
                ref_2bit: reference.path.clone(),
            },
            output_dir: ref_gc_dir.path().to_path_buf(),
            output_prefix: String::new(),
            n_threads: 1,
            n_positions: 41,
            seed: Some(11),
            windows: cfdnalab::commands::ref_gc_bias::config::RefGCWindowsArgs {
                by_bed: Some(bed_path),
            },
            chromosomes: ChromosomeArgs {
                chromosomes: Some(vec!["chr1".to_string()]),
                chromosomes_file: None,
            },
            blacklist: None,
            fragment_lengths: cfdnalab::commands::cli_common::FragmentLengthArgs {
                min_fragment_length: 10,
                max_fragment_length: 10,
            },
            end_offset: 0,
            skip_interpolation: false,
            smoothing_sigma: 0.55,
            smoothing_radius: 2,
            skip_smoothing: true,
            tile_size: 1_000_000,
            logging: LoggingArgs::default(),
        };

        let out_dir = TempDir::new()?;
        let mut cfg =
            make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());
        cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });
        cfg.set_num_extreme_gc_bins(0);
        cfg.set_num_short_length_bins(0);
        cfg.outlier_method = cfdnalab::commands::gc_bias::config::OutlierMethodArg::None;

        // Act
        run_ref_gc_bias(&ref_cfg)?;
        run_gc_bias(&cfg)?;

        // Assert
        let normalized_counts: ndarray::Array2<f64> = read_npy(
            out_dir
                .path()
                .join("gc_bias.normalized_avg_cfdna_counts.1.npy"),
        )?;
        assert_eq!(normalized_counts.dim(), (1, 101));
        for (gc_pct, &value) in normalized_counts.row(0).iter().enumerate() {
            let expected_value = if matches!(gc_pct, 0 | 50 | 100) {
                1.0
            } else {
                0.0
            };
            assert_gc_command_close(
                value,
                expected_value,
                &format!("normalized avg counts at GC% {gc_pct}"),
            );
        }

        let interpolated_counts: ndarray::Array2<f64> = read_npy(
            out_dir
                .path()
                .join("gc_bias.interpolated_cfdna_counts.2.npy"),
        )?;
        assert_eq!(interpolated_counts.dim(), (1, 101));
        for (gc_pct, &value) in interpolated_counts.row(0).iter().enumerate() {
            assert_gc_command_close(value, 1.0, &format!("interpolated counts at GC% {gc_pct}"));
        }

        Ok(())
    }

    #[test]
    fn save_intermediates_writes_expected_sequence_and_mean_scaled_average_counts() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Use a single global window and a reference package that already disables smoothing and
        // interpolation. In that configuration `gc-bias` should save exactly six intermediate
        // arrays:
        //   0 avg_cfdna_counts
        //   1 normalized_avg_cfdna_counts
        //   2 binned_ref_counts
        //   3 binned_cfdna_counts
        //   4 normalized_binned_cfdna_counts
        //   5 normalized_binned_ref_counts
        //
        // The strongest low-level coherence check in this branch is the first normalization step:
        // `normalized_avg_cfdna_counts` must equal `avg_cfdna_counts / supported_mean`, where the
        // mean is taken only over the reference outlier-support mask.
        let bam = fixtures::simple_inward_bam()?;
        let reference = fixtures::simple_reference_twobit()?;
        let ref_gc_dir = TempDir::new()?;
        write_reference_package_for_single_length(&reference.path, &ref_gc_dir, 60, 0)?;

        let out_dir = TempDir::new()?;
        let mut cfg =
            make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());
        cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });
        cfg.set_save_intermediates(true);

        // Act
        run_gc_bias(&cfg)?;

        // Assert:
        // No interpolation/smoothing intermediates should exist for this reference package, so the
        // numbering must stay dense across exactly six saved arrays.
        let mut intermediate_files: Vec<String> = std::fs::read_dir(out_dir.path())?
            .filter_map(|entry| entry.ok())
            .filter_map(|entry| {
                let name = entry.file_name().into_string().ok()?;
                if name.starts_with("gc_bias.") && name.ends_with(".npy") {
                    Some(name)
                } else {
                    None
                }
            })
            .collect();
        intermediate_files.sort();
        assert_eq!(
            intermediate_files,
            vec![
                "gc_bias.avg_cfdna_counts.0.npy".to_string(),
                "gc_bias.binned_cfdna_counts.3.npy".to_string(),
                "gc_bias.binned_ref_counts.2.npy".to_string(),
                "gc_bias.normalized_avg_cfdna_counts.1.npy".to_string(),
                "gc_bias.normalized_binned_cfdna_counts.4.npy".to_string(),
                "gc_bias.normalized_binned_ref_counts.5.npy".to_string(),
            ]
        );

        let avg_counts: ndarray::Array2<f64> =
            read_npy(out_dir.path().join("gc_bias.avg_cfdna_counts.0.npy"))?;
        let normalized_avg: ndarray::Array2<f64> = read_npy(
            out_dir
                .path()
                .join("gc_bias.normalized_avg_cfdna_counts.1.npy"),
        )?;
        let reference_data = load_reference_gc_data(&ref_gc_dir.path().join("ref_gc_package.npz"))?;

        // The support mask defines exactly which cells contribute to the mean-scaling denominator.
        let mut supported_sum = 0.0_f64;
        let mut supported_count = 0usize;
        for (value, supported) in avg_counts
            .iter()
            .zip(reference_data.outliers_support_mask.iter())
        {
            if *supported {
                supported_sum += *value;
                supported_count += 1;
            }
        }
        assert!(
            supported_count > 0,
            "fixture must have supported reference bins"
        );
        let supported_mean = supported_sum / supported_count as f64;
        assert!(
            supported_mean > 0.0,
            "supported mean must be positive for mean scaling"
        );

        for ((row_idx, col_idx), avg_value) in avg_counts.indexed_iter() {
            let expected = *avg_value / supported_mean;
            let actual = normalized_avg[(row_idx, col_idx)];
            assert!(
                (actual - expected).abs() < 1e-12,
                "normalized avg mismatch at ({row_idx}, {col_idx}): expected {expected}, got {actual}"
            );
        }

        Ok(())
    }

    #[test]
    fn min_window_acgt_pct_excludes_mostly_blacklisted_window_but_keeps_clean_window() -> Result<()>
    {
        // Human verification status: unverified
        // Arrange:
        // Two explicit 100 bp windows on a 200 bp chromosome:
        // - left  window [0,100)   is all A
        // - right window [100,200) is all C
        //
        // Blacklist [0,85), leaving only 15 usable ACGT bases in the left window.
        // We place one valid 10 bp fragment fully inside the surviving 15 bp tail:
        //   [85,95) -> GC%=0
        // and one valid 10 bp fragment in the clean right window:
        //   [110,120) -> GC%=100
        //
        // This lets us separate fragment-level validity from window-level validity:
        // - the left fragment is valid because its own 10 bp span is entirely unmasked
        // - but the left window has only 15 / 100 = 15% usable ACGT, so it should be rejected
        //   when `min_window_acgt_pct = 20`
        //
        // Window scaling math for length 10 (11 reachable GC-count cells):
        //
        // Threshold = 0:
        // - left window:
        //     raw count at GC%=0 = 1
        //     mean_count         = 1 / 11
        //     usable ACGT        = 15
        //     scale              = (1 / (1/11)) * (15 / 100) = 11 * 0.15 = 1.65
        //     scaled row         = 1.65 at GC%=0
        // - right window:
        //     raw count at GC%=100 = 1
        //     mean_count           = 1 / 11
        //     usable ACGT          = 100
        //     scale                = (1 / (1/11)) * (100 / 100) = 11
        //     scaled row           = 11 at GC%=100
        // - average over two kept windows:
        //     GC%=0   -> 1.65 / 2 = 0.825
        //     GC%=100 -> 11   / 2 = 5.5
        //
        // Threshold = 20:
        // - left window is rejected because 15% < 20%
        // - right window remains, so the final average is just:
        //     GC%=100 -> 11
        //
        // As elsewhere for length 10, width correction is a no-op because GC%=0 and 100 are exact
        // reachable percentage bins with width 1.
        let reference = fixtures::twobit_from_sequences(
            "gc_bias_min_window_acgt_reference",
            vec![(
                "chr1".to_string(),
                format!("{}{}", "A".repeat(100), "C".repeat(100)),
            )],
        )?;
        let bam = fixtures::bam_from_specs(
            vec![("chr1".to_string(), 200)],
            vec![
                fixtures::paired_fragment(85, 10, 5),
                fixtures::paired_fragment(110, 10, 5),
            ],
            Vec::new(),
            "gc_bias_min_window_acgt_bam",
        )?;
        let ref_gc_dir = TempDir::new()?;
        write_reference_package_for_single_length(&reference.path, &ref_gc_dir, 10, 0)?;

        let threshold0_out = TempDir::new()?;
        let threshold20_out = TempDir::new()?;
        let blacklist_path = threshold20_out.path().join("blacklist.bed");
        std::fs::write(&blacklist_path, "chr1\t0\t85\n")?;

        let mut threshold0_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            threshold0_out.path(),
        );
        threshold0_cfg.set_windows(GCWindowsArgs {
            by_size: Some(100),
            by_bed: None,
            global: false,
        });
        threshold0_cfg.set_tile_size(95);
        threshold0_cfg.set_blacklist(Some(vec![blacklist_path.clone()]));
        threshold0_cfg.set_min_window_acgt_pct(0);

        let mut threshold20_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            threshold20_out.path(),
        );
        threshold20_cfg.set_windows(GCWindowsArgs {
            by_size: Some(100),
            by_bed: None,
            global: false,
        });
        threshold20_cfg.set_tile_size(95);
        threshold20_cfg.set_blacklist(Some(vec![blacklist_path]));
        threshold20_cfg.set_min_window_acgt_pct(20);

        // Act
        run_gc_bias(&threshold0_cfg)?;
        run_gc_bias(&threshold20_cfg)?;

        // Assert
        let threshold0_counts: ndarray::Array2<f64> =
            read_npy(threshold0_out.path().join("gc_bias.avg_cfdna_counts.0.npy"))?;
        let threshold20_counts: ndarray::Array2<f64> = read_npy(
            threshold20_out
                .path()
                .join("gc_bias.avg_cfdna_counts.0.npy"),
        )?;
        assert_eq!(threshold0_counts.dim(), (1, 101));
        assert_eq!(threshold20_counts.dim(), (1, 101));

        for (gc_pct, &value) in threshold0_counts.row(0).iter().enumerate() {
            match gc_pct {
                0 => assert!(
                    (value - 0.825).abs() < 1e-12,
                    "threshold 0 expected 0.825 at GC% 0, got {value}"
                ),
                100 => assert!(
                    (value - 5.5).abs() < 1e-12,
                    "threshold 0 expected 5.5 at GC% 100, got {value}"
                ),
                _ => assert!(
                    value.abs() < 1e-12,
                    "threshold 0 expected 0 outside GC% 0/100, got bin {gc_pct}={value}"
                ),
            }
        }

        for (gc_pct, &value) in threshold20_counts.row(0).iter().enumerate() {
            match gc_pct {
                100 => assert!(
                    (value - 11.0).abs() < 1e-12,
                    "threshold 20 expected 11.0 at GC% 100, got {value}"
                ),
                _ => assert!(
                    value.abs() < 1e-12,
                    "threshold 20 expected 0 outside GC% 100, got bin {gc_pct}={value}"
                ),
            }
        }

        Ok(())
    }

    #[test]
    fn multiple_blacklist_files_with_touching_intervals_match_single_merged_gc_bias_run()
    -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // `gc-bias` uses the shared blacklist loader with:
        // - `min_size = 1`
        // - `halo_bp = 0`
        //
        // Therefore, touching intervals from separate files:
        //   [40, 60) and [60, 80)
        // must be merged into the same effective masked span as:
        //   [40, 80)
        //
        // We compare two otherwise identical global runs. To keep the fixture realistic, place the
        // touching blacklist away from the only fragment in `simple_inward_bam()`:
        // - fragment span is [20,80)
        // - touching masked span is [120,160)
        //
        // So the blacklist is a real loaded artifact, but it does not invalidate the only
        // fragment's GC context. Because the effective masked genomic coordinates are logically
        // identical, both the saved `avg_cfdna_counts` matrix and the final correction package
        // must match exactly.
        let bam = fixtures::simple_inward_bam()?;
        let reference = fixtures::simple_reference_twobit()?;
        let ref_gc_dir = TempDir::new()?;
        write_reference_package_for_single_length(&reference.path, &ref_gc_dir, 60, 0)?;

        let split_out_dir = TempDir::new()?;
        let merged_out_dir = TempDir::new()?;
        let blacklist_dir = TempDir::new()?;
        let split_a = blacklist_dir.path().join("blacklist_a.bed");
        let split_b = blacklist_dir.path().join("blacklist_b.bed");
        let merged = blacklist_dir.path().join("blacklist_merged.bed");
        std::fs::write(&split_a, "chr1\t120\t140\n")?;
        std::fs::write(&split_b, "chr1\t140\t160\n")?;
        std::fs::write(&merged, "chr1\t120\t160\n")?;

        let make_cfg =
            |output_dir: &std::path::Path, blacklist: Vec<std::path::PathBuf>| -> GCConfig {
                let mut cfg =
                    make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), output_dir);
                cfg.set_windows(GCWindowsArgs {
                    by_size: None,
                    by_bed: None,
                    global: true,
                });
                cfg.set_blacklist(Some(blacklist));
                cfg
            };

        // Act
        run_gc_bias(&make_cfg(
            split_out_dir.path(),
            vec![split_a.clone(), split_b.clone()],
        ))?;
        run_gc_bias(&make_cfg(merged_out_dir.path(), vec![merged.clone()]))?;

        // Assert
        let split_avg: ndarray::Array2<f64> =
            read_npy(split_out_dir.path().join("gc_bias.avg_cfdna_counts.0.npy"))?;
        let merged_avg: ndarray::Array2<f64> =
            read_npy(merged_out_dir.path().join("gc_bias.avg_cfdna_counts.0.npy"))?;
        assert_eq!(split_avg, merged_avg);
        assert_eq!(split_avg.dim(), (1, 101));
        for (gc_pct, &value) in split_avg.row(0).iter().enumerate() {
            match gc_pct {
                50 => assert!(
                    (value - 1.0).abs() < 1e-12,
                    "expected exactly one surviving GC%=50 fragment, got {value}"
                ),
                _ => assert!(
                    value.abs() < 1e-12,
                    "expected no mass outside GC% 50, got bin {gc_pct}={value}"
                ),
            }
        }

        let split_package =
            GCCorrectionPackage::from_file(split_out_dir.path().join("gc_bias_correction.npz"))?;
        let merged_package =
            GCCorrectionPackage::from_file(merged_out_dir.path().join("gc_bias_correction.npz"))?;
        assert_eq!(split_package.version, merged_package.version);
        assert_eq!(split_package.end_offset, merged_package.end_offset);
        assert_eq!(split_package.length_edges, merged_package.length_edges);
        assert_eq!(split_package.gc_edges, merged_package.gc_edges);
        assert_eq!(
            split_package.length_bin_frequencies,
            merged_package.length_bin_frequencies
        );
        assert_eq!(
            split_package.correction_matrix,
            merged_package.correction_matrix
        );
        assert_eq!(split_package.version, GC_CORRECTION_SCHEMA_VERSION);
        assert_eq!(split_package.end_offset, 0);
        assert_eq!(split_package.length_edges, vec![60, 60]);
        assert_eq!(split_package.gc_edges, vec![0, 100]);
        assert_eq!(split_package.length_bin_frequencies.len(), 1);
        assert!((split_package.length_bin_frequencies[0] - 1.0).abs() < 1e-12);
        assert_eq!(split_package.correction_matrix.dim(), (1, 1));
        assert!((split_package.correction_matrix[(0, 0)] - 1.0).abs() < 1e-12);

        Ok(())
    }

    #[test]
    fn gc_bias_run_rejects_reference_package_with_non_scalar_metadata_array() -> Result<()> {
        // Human verification status: unverified
        let bam = fixtures::simple_inward_bam()?;
        let reference = fixtures::simple_reference_twobit()?;
        let ref_gc_dir = TempDir::new()?;
        write_reference_gc_package_fixture(
            ref_gc_dir.path(),
            &[GC_CORRECTION_SCHEMA_VERSION],
            &[false],
            &[2],
            &[0.55],
            &[true, false],
        )?;
        let out_dir = TempDir::new()?;
        let cfg = make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());

        // Manual expectations:
        // - The reference package is malformed before the command ever starts sample counting:
        //   `skip_smoothing` is written as a length-2 array instead of a scalar metadata field.
        // - `gc-bias` loads the reference package at the start of `run()`, so the correct
        //   behavior is an immediate loader failure with the scalar-shape guardrail message.
        let err = run_gc_bias(&cfg).expect_err("non-scalar reference metadata should fail");
        let msg = err.to_string();
        assert!(
            msg.contains("skip_smoothing should be length 1"),
            "unexpected error message: {msg}"
        );

        Ok(())
    }

    #[test]
    fn gc_bias_run_rejects_reference_package_with_schema_version_mismatch() -> Result<()> {
        // Human verification status: unverified
        let bam = fixtures::simple_inward_bam()?;
        let reference = fixtures::simple_reference_twobit()?;
        let ref_gc_dir = TempDir::new()?;
        write_reference_gc_package_fixture(
            ref_gc_dir.path(),
            &[GC_CORRECTION_SCHEMA_VERSION + 1],
            &[false],
            &[2],
            &[0.55],
            &[true],
        )?;
        let out_dir = TempDir::new()?;
        let cfg = make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());

        // Manual expectations:
        // - The package schema version is intentionally incompatible.
        // - `gc-bias` must fail while loading the reference-GC artifact, before producing any
        //   sample-side intermediates or correction output.
        let err = run_gc_bias(&cfg).expect_err("schema version mismatch should fail");
        let msg = err.to_string();
        assert!(
            msg.contains("Reference GC package schema version mismatch"),
            "unexpected error message: {msg}"
        );

        Ok(())
    }

    #[test]
    fn apply_outliers_per_length_winsorizes_rows() {
        // Human verification status: unverified
        let mut matrix = array![[1.0_f64, 2.0_f64, 100.0_f64], [1.0_f64, 5.0_f64, 6.0_f64]];
        let mask = array![[true, true, true], [true, true, true]];

        let stats = apply_outliers_to_matrix(
            &mut matrix,
            Some(&mask),
            OutlierScope::PerLength,
            OutlierRule::Quantile {
                lower: 0.0,
                upper: 0.5,
            },
            OutlierAction::Winsorize,
        );

        assert_eq!(matrix[[0, 0]], 1.0);
        assert_eq!(matrix[[0, 1]], 2.0);
        assert_eq!(matrix[[0, 2]], 2.0); // Clamped
        assert_eq!(matrix[[1, 0]], 1.0);
        assert_eq!(matrix[[1, 1]], 5.0);
        assert_eq!(matrix[[1, 2]], 5.0); // Clamped
        assert_eq!(
            stats,
            OutlierStats {
                total_examined: 6,
                total_outliers_handled: 2,
                unsupported_examined: 0,
                unsupported_outliers_handled: 0,
                hard_clamped: 0
            }
        );
    }

    #[test]
    fn quantile_outliers_symmetry_clamps_extremes() {
        // Human verification status: unverified
        let mut matrix = array![[1.0_f64, 1.0_f64, 100.0_f64]];

        apply_outliers_to_matrix(
            &mut matrix,
            None,
            OutlierScope::Global,
            OutlierRule::Quantile {
                lower: 0.25,
                upper: 0.75,
            },
            OutlierAction::Winsorize,
        );

        assert!((matrix[[0, 0]] - 1.0).abs() < 1e-6);
        assert!((matrix[[0, 1]] - 1.0).abs() < 1e-6);
        assert!((matrix[[0, 2]] - 50.5).abs() < 1e-6);
    }

    #[test]
    fn masked_cells_are_clamped_but_not_counted() {
        // Human verification status: unverified
        let mut matrix = array![[1.0_f64, 2.0_f64, 100.0_f64]];
        let mask = array![[true, true, false]];

        let stats = apply_outliers_to_matrix(
            &mut matrix,
            Some(&mask),
            OutlierScope::Global,
            OutlierRule::TukeyIqr { k: 1.0 },
            OutlierAction::Winsorize,
        );

        assert!((matrix[[0, 0]] - 1.0).abs() < 1e-6);
        assert!((matrix[[0, 1]] - 2.0).abs() < 1e-6);
        assert!((matrix[[0, 2]] - 2.25).abs() < 1e-6); // Unsupported cell still clamped
        assert_eq!(
            stats,
            OutlierStats {
                total_examined: 2,
                total_outliers_handled: 0,
                unsupported_examined: 1,
                unsupported_outliers_handled: 1,
                hard_clamped: 0
            }
        );
    }

    #[test]
    fn interpolated_quantile_weights_neighbors_by_offset() {
        // Human verification status: unverified
        // Arrange
        let values = vec![0.0_f32, 10.0_f32, 20.0_f32, 30.0_f32, 40.0_f32];

        // Act
        let p_0 = interpolated_quantile(&values, 0.0);
        let p_05 = interpolated_quantile(&values, 0.5);
        let p_06 = interpolated_quantile(&values, 0.6);
        let p_08 = interpolated_quantile(&values, 0.8);
        let p_1 = interpolated_quantile(&values, 1.0);

        // Assert
        assert!((p_0 - 0.0).abs() < 1e-6);
        assert!((p_05 - 20.0).abs() < 1e-6);
        assert!((p_06 - 24.0).abs() < 1e-6); // 40% from 20 to 30
        assert!((p_08 - 32.0).abs() < 1e-6); // 20% from 30 to 40
        assert!((p_1 - 40.0).abs() < 1e-6);
    }

    #[test]
    fn quantile_bounds_interpolate_between_indices() {
        // Human verification status: unverified
        // Arrange: Percentiles fall between indices, so bounds should blend neighbors
        let values = vec![0.0_f32, 10.0_f32, 20.0_f32, 30.0_f32, 40.0_f32];

        // Act: compute bounds for percentiles that require interpolation.
        let bounds = outlier_bounds(
            &values,
            OutlierRule::Quantile {
                lower: 0.6,
                upper: 0.8,
            },
        )
        .expect("quantile bounds should exist");

        // Assert: 0.6 is 40% from element 2 (20) to 3 (30); 0.8 is 20% from 3 (30) to 4 (40)
        assert!((bounds.0 - 24.0).abs() < 1e-6);
        assert!((bounds.1 - 32.0).abs() < 1e-6);
    }

    #[test]
    fn iqr_outliers_per_length_clamps_high_values() {
        // Human verification status: unverified
        let mut matrix = array![[1.0_f64, 2.0_f64, 8.0_f64]];

        apply_outliers_to_matrix(
            &mut matrix,
            None,
            OutlierScope::PerLength,
            OutlierRule::TukeyIqr { k: 0.5 },
            OutlierAction::Winsorize,
        );

        assert!((matrix[[0, 0]] - 1.0).abs() < 1e-6);
        assert!((matrix[[0, 1]] - 2.0).abs() < 1e-6);
        assert!((matrix[[0, 2]] - 6.75).abs() < 1e-6);
    }

    #[test]
    fn stddev_outliers_global_clamps_tail() {
        // Human verification status: unverified
        let mut matrix = array![[1.0_f64, 1.0_f64, 10.0_f64]];

        apply_outliers_to_matrix(
            &mut matrix,
            None,
            OutlierScope::Global,
            OutlierRule::StdDev { k: 1.0 },
            OutlierAction::Winsorize,
        );

        assert!((matrix[[0, 2]] - 8.2426405).abs() < 1e-5);
    }

    #[test]
    fn mad_outliers_symmetrically_clamp() {
        // Human verification status: unverified
        let mut matrix = array![[1.0_f64, 2.0_f64, 3.0_f64, 9.0_f64]];

        apply_outliers_to_matrix(
            &mut matrix,
            None,
            OutlierScope::Global,
            OutlierRule::Mad { k: 1.0 },
            OutlierAction::Winsorize,
        );

        assert!((matrix[[0, 0]] - 1.0174).abs() < 1e-4);
        assert!((matrix[[0, 1]] - 2.0).abs() < 1e-6);
        assert!((matrix[[0, 2]] - 3.0).abs() < 1e-6);
        assert!((matrix[[0, 3]] - 3.9826).abs() < 1e-4);
    }

    #[test]
    fn per_length_scope_differs_from_global() {
        // Human verification status: unverified
        let mut matrix = array![[1.0_f64, 100.0_f64], [1.0_f64, 1.0_f64]];

        apply_outliers_to_matrix(
            &mut matrix,
            None,
            OutlierScope::PerLength,
            OutlierRule::Quantile {
                lower: 0.25,
                upper: 0.75,
            },
            OutlierAction::Winsorize,
        );

        assert!((matrix[[0, 0]] - 25.75).abs() < 1e-6);
        assert!((matrix[[0, 1]] - 75.25).abs() < 1e-6);
        assert!((matrix[[1, 0]] - 1.0).abs() < 1e-6);
        assert!((matrix[[1, 1]] - 1.0).abs() < 1e-6);
    }

    fn write_reference_gc_package_fixture(
        out_dir: &std::path::Path,
        version: &[u32],
        skip_interpolation: &[bool],
        smoothing_radius: &[u32],
        smoothing_sigma: &[f64],
        skip_smoothing: &[bool],
    ) -> Result<()> {
        let package_path = out_dir.join("ref_gc_package.npz");
        let counts = array![[1.0_f64, 2.0_f64], [3.0_f64, 4.0_f64]];
        let support_unobservables = array![[true, false], [true, true]];
        let support_outliers = array![[true, true], [false, true]];
        let gc_percent_widths = array![[10_u16, 20_u16], [30_u16, 40_u16]];

        let file = std::fs::File::create(&package_path)?;
        let mut npz = NpzWriter::new(file);
        npz.add_array("counts", &counts)?;
        npz.add_array("support_mask_unobservables", &support_unobservables)?;
        npz.add_array("support_mask_outliers", &support_outliers)?;
        npz.add_array("gc_percent_widths", &gc_percent_widths)?;
        npz.add_array("version", &ndarray::Array1::from(version.to_vec()))?;
        npz.add_array("length_range", &ndarray::Array1::from(vec![30_u32, 40_u32]))?;
        npz.add_array("end_offset", &ndarray::Array1::from(vec![10_u32]))?;
        npz.add_array(
            "skip_interpolation",
            &ndarray::Array1::from(skip_interpolation.to_vec()),
        )?;
        npz.add_array(
            "smoothing_radius",
            &ndarray::Array1::from(smoothing_radius.to_vec()),
        )?;
        npz.add_array(
            "smoothing_sigma",
            &ndarray::Array1::from(smoothing_sigma.to_vec()),
        )?;
        npz.add_array(
            "skip_smoothing",
            &ndarray::Array1::from(skip_smoothing.to_vec()),
        )?;
        npz.finish()?;
        Ok(())
    }

    fn write_reference_gc_package_with_shape_mismatch(out_dir: &std::path::Path) -> Result<()> {
        let package_path = out_dir.join("ref_gc_package.npz");
        let counts = array![[1.0_f64, 2.0_f64], [3.0_f64, 4.0_f64]];
        let support_unobservables = array![[true, false]];
        let support_outliers = array![[true, true]];
        let gc_percent_widths = array![[10_u16, 20_u16], [30_u16, 40_u16]];

        let file = std::fs::File::create(&package_path)?;
        let mut npz = NpzWriter::new(file);
        npz.add_array("counts", &counts)?;
        npz.add_array("support_mask_unobservables", &support_unobservables)?;
        npz.add_array("support_mask_outliers", &support_outliers)?;
        npz.add_array("gc_percent_widths", &gc_percent_widths)?;
        npz.add_array(
            "version",
            &ndarray::Array1::from(vec![GC_CORRECTION_SCHEMA_VERSION]),
        )?;
        npz.add_array("length_range", &ndarray::Array1::from(vec![30_u32, 40_u32]))?;
        npz.add_array("end_offset", &ndarray::Array1::from(vec![10_u32]))?;
        npz.add_array("skip_interpolation", &ndarray::Array1::from(vec![false]))?;
        npz.add_array("smoothing_radius", &ndarray::Array1::from(vec![2_u32]))?;
        npz.add_array("smoothing_sigma", &ndarray::Array1::from(vec![0.55_f64]))?;
        npz.add_array("skip_smoothing", &ndarray::Array1::from(vec![true]))?;
        npz.finish()?;
        Ok(())
    }

    fn write_two_bin_reference_gc_package(
        out_dir: &std::path::Path,
        length_range: (u32, u32),
    ) -> Result<()> {
        let package_path = out_dir.join("ref_gc_package.npz");
        let n_lengths = (length_range.1 - length_range.0 + 1) as usize;
        let mut counts = ndarray::Array2::<f64>::zeros((n_lengths, 101));
        let mut support_outliers = ndarray::Array2::<bool>::from_elem((n_lengths, 101), false);
        let gc_percent_widths =
            gc_percent_widths(length_range.0 as usize, length_range.1 as usize, 0);
        let support_unobservables = gc_percent_widths.mapv(|width| width > 0);

        for row_idx in 0..n_lengths {
            counts[(row_idx, 0)] = 1.0;
            counts[(row_idx, 100)] = 1.0;
            support_outliers[(row_idx, 0)] = true;
            support_outliers[(row_idx, 100)] = true;
        }

        let file = std::fs::File::create(&package_path)?;
        let mut npz = NpzWriter::new(file);
        npz.add_array("counts", &counts)?;
        npz.add_array("support_mask_unobservables", &support_unobservables)?;
        npz.add_array("support_mask_outliers", &support_outliers)?;
        npz.add_array("gc_percent_widths", &gc_percent_widths)?;
        npz.add_array(
            "version",
            &ndarray::Array1::from(vec![GC_CORRECTION_SCHEMA_VERSION]),
        )?;
        npz.add_array(
            "length_range",
            &ndarray::Array1::from(vec![length_range.0, length_range.1]),
        )?;
        npz.add_array("end_offset", &ndarray::Array1::from(vec![0_u32]))?;
        npz.add_array("skip_interpolation", &ndarray::Array1::from(vec![true]))?;
        npz.add_array("smoothing_radius", &ndarray::Array1::from(vec![2_u32]))?;
        npz.add_array("smoothing_sigma", &ndarray::Array1::from(vec![0.55_f64]))?;
        npz.add_array("skip_smoothing", &ndarray::Array1::from(vec![true]))?;
        npz.finish()?;
        Ok(())
    }

    fn write_balanced_two_length_reference_gc_package(out_dir: &std::path::Path) -> Result<()> {
        let package_path = out_dir.join("ref_gc_package.npz");

        // Hand-built but still realistic reference package for run-level outlier tests.
        //
        // The package covers two fragment lengths, 10 and 11. We intentionally place reference
        // mass only at GC% 0 and GC% 100 for both rows, because the paired sample fixture also
        // lives only in those two extreme classes.
        //
        // But the metadata must still look like a real package:
        // - `support_mask_unobservables` follows the true theoretical GC%-reachability geometry
        //   for lengths 10 and 11
        // - `gc_percent_widths` uses the real rounding widths for those lengths
        // - `support_mask_outliers` marks only the empirically populated bins (0 and 100)
        //
        // That keeps the reference-side normalization denominator restricted to the two populated
        // bins while still preserving realistic theoretical metadata.
        //
        // With that support mask, once `gc-bias` bins the GC axis into `[0]` and `[1..100]`, the
        // reference-side binned rows are perfectly balanced:
        //   length 10 -> [1, 1]
        //   length 11 -> [1, 1]
        //
        // That keeps the run-level derivation simple: after per-length normalization the raw
        // correction matrix is driven entirely by the sample BAM's within-row imbalance.
        let mut counts = ndarray::Array2::<f64>::zeros((2, 101));
        counts[(0, 0)] = 1.0;
        counts[(0, 100)] = 1.0;
        counts[(1, 0)] = 1.0;
        counts[(1, 100)] = 1.0;

        let gc_percent_widths = gc_percent_widths(10, 11, 0);
        let support_unobservables = gc_percent_widths.mapv(|width| width > 0);
        let mut support_outliers = ndarray::Array2::<bool>::from_elem((2, 101), false);
        support_outliers[(0, 0)] = true;
        support_outliers[(0, 100)] = true;
        support_outliers[(1, 0)] = true;
        support_outliers[(1, 100)] = true;

        let file = std::fs::File::create(&package_path)?;
        let mut npz = NpzWriter::new(file);
        npz.add_array("counts", &counts)?;
        npz.add_array("support_mask_unobservables", &support_unobservables)?;
        npz.add_array("support_mask_outliers", &support_outliers)?;
        npz.add_array("gc_percent_widths", &gc_percent_widths)?;
        npz.add_array(
            "version",
            &ndarray::Array1::from(vec![GC_CORRECTION_SCHEMA_VERSION]),
        )?;
        npz.add_array("length_range", &ndarray::Array1::from(vec![10_u32, 11_u32]))?;
        npz.add_array("end_offset", &ndarray::Array1::from(vec![0_u32]))?;
        npz.add_array("skip_interpolation", &ndarray::Array1::from(vec![true]))?;
        npz.add_array("smoothing_radius", &ndarray::Array1::from(vec![2_u32]))?;
        npz.add_array("smoothing_sigma", &ndarray::Array1::from(vec![0.55_f64]))?;
        npz.add_array("skip_smoothing", &ndarray::Array1::from(vec![true]))?;
        npz.finish()?;
        Ok(())
    }

    fn write_three_bin_reference_gc_package(out_dir: &std::path::Path) -> Result<()> {
        let package_path = out_dir.join("ref_gc_package.npz");

        // One fragment length, with reference mass only at GC% 0, 50, and 100.
        //
        // As above, keep the metadata realistic:
        // - theoretical support and width correction follow the true length-10 rounding geometry
        // - empirical outlier support is restricted to the three populated GC bins
        //
        // That keeps the reference-side normalization easy to reason about for GC binning tests:
        // every empirically supported GC point starts with the same mass 1.0.
        let mut counts = ndarray::Array2::<f64>::zeros((1, 101));
        counts[(0, 0)] = 1.0;
        counts[(0, 50)] = 1.0;
        counts[(0, 100)] = 1.0;

        let gc_percent_widths = gc_percent_widths(10, 10, 0);
        let support_unobservables = gc_percent_widths.mapv(|width| width > 0);
        let mut support_outliers = ndarray::Array2::<bool>::from_elem((1, 101), false);
        support_outliers[(0, 0)] = true;
        support_outliers[(0, 50)] = true;
        support_outliers[(0, 100)] = true;

        let file = std::fs::File::create(&package_path)?;
        let mut npz = NpzWriter::new(file);
        npz.add_array("counts", &counts)?;
        npz.add_array("support_mask_unobservables", &support_unobservables)?;
        npz.add_array("support_mask_outliers", &support_outliers)?;
        npz.add_array("gc_percent_widths", &gc_percent_widths)?;
        npz.add_array(
            "version",
            &ndarray::Array1::from(vec![GC_CORRECTION_SCHEMA_VERSION]),
        )?;
        npz.add_array("length_range", &ndarray::Array1::from(vec![10_u32, 10_u32]))?;
        npz.add_array("end_offset", &ndarray::Array1::from(vec![0_u32]))?;
        npz.add_array("skip_interpolation", &ndarray::Array1::from(vec![true]))?;
        npz.add_array("smoothing_radius", &ndarray::Array1::from(vec![2_u32]))?;
        npz.add_array("smoothing_sigma", &ndarray::Array1::from(vec![0.55_f64]))?;
        npz.add_array("skip_smoothing", &ndarray::Array1::from(vec![true]))?;
        npz.finish()?;
        Ok(())
    }

    fn make_two_length_outlier_fixture() -> Result<(fixtures::TwoBitFixture, fixtures::BamFixture)>
    {
        let reference = fixtures::twobit_from_sequences(
            "gc_bias_two_length_outlier_reference",
            vec![(
                "chr1".to_string(),
                format!("{}{}", "A".repeat(100), "C".repeat(100)),
            )],
        )?;

        // Length-10 row:
        // - one pure-A fragment  -> GC%=0
        // - nine pure-C fragments -> GC%=100
        //
        // Length-11 row:
        // - five pure-A fragments  -> GC%=0
        // - five pure-C fragments -> GC%=100
        //
        // Later, after global mean-scaling and GC binning into `[0]` and `[1..100]`, the sample
        // binned rows are proportional to:
        //   length 10 -> [1, 9] -> normalized to [0.2, 1.8]
        //   length 11 -> [5, 5] -> normalized to [1.0, 1.0]
        let mut fragments = Vec::new();
        fragments.push(fixtures::paired_fragment(10, 10, 5));
        for start in [110_i64, 120, 130, 140, 150, 160, 170, 180, 190] {
            fragments.push(fixtures::paired_fragment(start, 10, 5));
        }
        for start in [20_i64, 30, 40, 50, 60] {
            fragments.push(fixtures::paired_fragment(start, 11, 5));
        }
        for start in [120_i64, 130, 140, 150, 160] {
            fragments.push(fixtures::paired_fragment(start, 11, 5));
        }

        // Several length-10 and length-11 fragments deliberately share the same left start. Use
        // strict BAM identity here so those stacked molecules do not collapse onto one qname.
        let bam = fixtures::bam_from_specs_strict_identity(
            vec![("chr1".to_string(), 200)],
            fragments,
            Vec::new(),
            "gc_bias_two_length_outlier_bam",
        )?;
        Ok((reference, bam))
    }

    fn make_two_length_low_mass_tail_fixture()
    -> Result<(fixtures::TwoBitFixture, fixtures::BamFixture)> {
        let reference = fixtures::twobit_from_sequences(
            "gc_bias_two_length_low_mass_tail_reference",
            vec![(
                "chr1".to_string(),
                format!("{}{}", "A".repeat(100), "C".repeat(100)),
            )],
        )?;

        // Length-10 row:
        // - one pure-A fragment  -> GC%=0
        // - nine pure-C fragments -> GC%=100
        //
        // Length-11 row:
        // - one pure-A fragment -> GC%=0
        // - one pure-C fragment -> GC%=100
        //
        // So before any length binning the per-row normalized correction rows are:
        //   length 10 -> [0.2, 1.8]
        //   length 11 -> [1.0, 1.0]
        //
        // But the total row masses are deliberately unequal:
        //   length 10 -> 10
        //   length 11 -> 2
        //
        // That makes length 11 a clean "low-mass tail" row for testing greedy length binning by
        // percentage mass.
        let mut fragments = Vec::new();
        fragments.push(fixtures::paired_fragment(10, 10, 5));
        for start in [110_i64, 120, 130, 140, 150, 160, 170, 180, 190] {
            fragments.push(fixtures::paired_fragment(start, 10, 5));
        }
        fragments.push(fixtures::paired_fragment(20, 11, 5));
        fragments.push(fixtures::paired_fragment(120, 11, 5));

        // The low-mass tail fixture reuses start 120 across two distinct fragment lengths, so
        // each synthetic fragment needs its own qname.
        let bam = fixtures::bam_from_specs_strict_identity(
            vec![("chr1".to_string(), 200)],
            fragments,
            Vec::new(),
            "gc_bias_two_length_low_mass_tail_bam",
        )?;
        Ok((reference, bam))
    }

    fn make_three_gc_bin_fixture() -> Result<(fixtures::TwoBitFixture, fixtures::BamFixture)> {
        let reference = fixtures::twobit_from_sequences(
            "gc_bias_three_gc_bin_reference",
            vec![(
                "chr1".to_string(),
                format!(
                    "{}{}{}",
                    "A".repeat(40),
                    "CCCCCAAAAA".repeat(4),
                    "C".repeat(80)
                ),
            )],
        )?;

        // One fragment length only: 10 bp.
        //
        // Chosen starts land in three exact GC classes:
        // - start 10  -> AAAAAAAAAA   -> GC%=0
        // - start 40  -> CCCCCAAAAA   -> GC%=50
        // - start 100 -> CCCCCCCCCC   -> GC%=100
        //
        // Counts are deliberately imbalanced:
        //   GC%=0   -> 1 fragment
        //   GC%=50  -> 5 fragments
        //   GC%=100 -> 9 fragments
        let mut fragments = Vec::new();
        fragments.push(fixtures::paired_fragment(10, 10, 5));
        for _ in 0..5 {
            fragments.push(fixtures::paired_fragment(40, 10, 5));
        }
        for _ in 0..9 {
            fragments.push(fixtures::paired_fragment(100, 10, 5));
        }

        // This fixture intentionally stacks five fragments at GC%=50 and nine at GC%=100. Use
        // strict identity so repeated starts still represent repeated molecules.
        let bam = fixtures::bam_from_specs_strict_identity(
            vec![("chr1".to_string(), 160)],
            fragments,
            Vec::new(),
            "gc_bias_three_gc_bin_bam",
        )?;
        Ok((reference, bam))
    }

    #[test]
    fn loads_versioned_reference_gc_package() -> Result<()> {
        // Human verification status: unverified
        // Arrange: write a minimal reference package with the current schema version and scalar
        // metadata fields.
        let tmp = tempdir()?;
        write_reference_gc_package_fixture(
            tmp.path(),
            &[GC_CORRECTION_SCHEMA_VERSION],
            &[false],
            &[2],
            &[0.55],
            &[true],
        )?;

        // Act
        let loaded = load_reference_gc_data(&tmp.path().join("ref_gc_package.npz"))?;

        // Assert: arrays and scalar metadata survive round-trip exactly.
        assert_eq!(
            loaded.counts,
            array![[1.0_f64, 2.0_f64], [3.0_f64, 4.0_f64]]
        );
        assert_eq!(
            loaded.unobservables_support_mask,
            array![[true, false], [true, true]]
        );
        assert_eq!(
            loaded.outliers_support_mask,
            array![[true, true], [false, true]]
        );
        assert_eq!(
            loaded.gc_percent_widths,
            array![[10_u16, 20_u16], [30_u16, 40_u16]]
        );
        assert_eq!(loaded.metadata.min_fragment_length, 30);
        assert_eq!(loaded.metadata.max_fragment_length, 40);
        assert_eq!(loaded.metadata.end_offset, 10);
        assert!(!loaded.metadata.skip_interpolation);
        assert_eq!(loaded.metadata.smoothing_radius, 2);
        assert!((loaded.metadata.smoothing_sigma - 0.55).abs() < 1e-12);
        assert!(loaded.metadata.skip_smoothing);
        Ok(())
    }

    #[test]
    fn rejects_reference_gc_package_with_non_scalar_metadata_array() -> Result<()> {
        // Human verification status: unverified
        // Arrange: `skip_smoothing` is written with two values. This should fail cleanly instead of
        // indexing `[0]` and panicking.
        let tmp = tempdir()?;
        write_reference_gc_package_fixture(
            tmp.path(),
            &[GC_CORRECTION_SCHEMA_VERSION],
            &[false],
            &[2],
            &[0.55],
            &[true, false],
        )?;

        // Act
        let error = load_reference_gc_data(&tmp.path().join("ref_gc_package.npz"))
            .expect_err("expected scalar-length error");

        // Assert
        assert!(
            error
                .to_string()
                .contains("skip_smoothing should be length 1")
        );
        Ok(())
    }

    #[test]
    fn rejects_reference_gc_package_with_schema_version_mismatch() -> Result<()> {
        // Human verification status: unverified
        // Arrange: same package shape, but an incompatible version number.
        let tmp = tempdir()?;
        write_reference_gc_package_fixture(
            tmp.path(),
            &[GC_CORRECTION_SCHEMA_VERSION + 1],
            &[false],
            &[2],
            &[0.55],
            &[true],
        )?;

        // Act
        let error = load_reference_gc_data(&tmp.path().join("ref_gc_package.npz"))
            .expect_err("expected schema version mismatch");

        // Assert
        assert!(
            error
                .to_string()
                .contains("Reference GC package schema version mismatch")
        );
        Ok(())
    }

    #[test]
    fn gc_bias_run_rejects_reference_package_with_incompatible_support_mask_shape() -> Result<()> {
        // Human verification status: unverified
        let bam = fixtures::simple_inward_bam()?;
        let reference = fixtures::simple_reference_twobit()?;
        let ref_gc_dir = TempDir::new()?;
        write_reference_gc_package_with_shape_mismatch(ref_gc_dir.path())?;
        let out_dir = TempDir::new()?;
        let cfg = make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());

        // Manual expectations:
        // - The reference package is syntactically present, but the support masks are the wrong
        //   shape for the count matrix:
        //     counts                     = (2, 2)
        //     support_mask_unobservables = (1, 2)
        //     support_mask_outliers      = (1, 2)
        // - `gc-bias` must reject this immediately while loading the artifact rather than trying
        //   to continue with inconsistent masking semantics.
        let err = run_gc_bias(&cfg).expect_err("shape-mismatched reference package should fail");
        let msg = err.to_string();
        assert!(
            msg.contains("Reference counts") && msg.contains("incompatible shapes"),
            "unexpected error message: {msg}"
        );

        Ok(())
    }

    #[test]
    fn quantile_outlier_method_changes_real_command_correction_matrix_in_expected_way() -> Result<()>
    {
        // Human verification status: unverified
        // Arrange:
        // Use the synthetic two-length reference package and BAM fixture defined above.
        //
        // Reference package after GC binning:
        // - GC bins are `[0]` and `[1..100]`
        // - both length rows are balanced, so normalized reference rows are:
        //     length 10 -> [1.0, 1.0]
        //     length 11 -> [1.0, 1.0]
        //
        // Sample BAM after the same GC binning:
        // - length 10 raw row is [1, 9], so normalized row is [0.2, 1.8]
        // - length 11 raw row is [5, 5], so normalized row is [1.0, 1.0]
        //
        // Therefore, before outlier handling, the raw correction matrix is:
        //   [[0.2, 1.8],
        //    [1.0, 1.0]]
        //
        // `--outlier-method none`:
        // - no winsorization
        // - no hard clamp, since all values already lie inside [0.1, 10]
        // - inversion gives:
        //     [[5.0, 5/9],
        //      [1.0, 1.0]]
        //
        // `--outlier-method quantile --outlier-scope per-length --outlier-quantiles 0.25,0.75`:
        // - length-10 row sorted values are [0.2, 1.8]
        // - with the command's linear interpolation:
        //     Q25 = 0.2 + 0.25 * (1.8 - 0.2) = 0.6
        //     Q75 = 0.2 + 0.75 * (1.8 - 0.2) = 1.4
        // - winsorized length-10 row becomes [0.6, 1.4]
        // - length-11 row stays [1.0, 1.0]
        // - inversion gives:
        //     [[5/3, 5/7],
        //      [1.0, 1.0]]
        let (reference, bam) = make_two_length_outlier_fixture()?;
        let ref_gc_dir = TempDir::new()?;
        write_balanced_two_length_reference_gc_package(ref_gc_dir.path())?;

        let out_none = TempDir::new()?;
        let out_quantile = TempDir::new()?;

        let mut none_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            out_none.path(),
        );
        none_cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });
        none_cfg.set_min_length_bin_mass(0.0);
        none_cfg.set_min_length_bin_width(1);
        none_cfg.set_min_gc_bin_mass(1.0);
        none_cfg.set_num_extreme_gc_bins(0);
        none_cfg.set_num_short_length_bins(0);
        none_cfg.outlier_method = cfdnalab::commands::gc_bias::config::OutlierMethodArg::None;

        let mut quantile_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            out_quantile.path(),
        );
        quantile_cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });
        quantile_cfg.set_min_length_bin_mass(0.0);
        quantile_cfg.set_min_length_bin_width(1);
        quantile_cfg.set_min_gc_bin_mass(1.0);
        quantile_cfg.set_num_extreme_gc_bins(0);
        quantile_cfg.set_num_short_length_bins(0);
        quantile_cfg.outlier_method =
            cfdnalab::commands::gc_bias::config::OutlierMethodArg::Quantile;
        quantile_cfg.outlier_scope =
            cfdnalab::commands::gc_bias::config::OutlierScopeArg::PerLength;
        quantile_cfg.outlier_quantiles = vec![0.25, 0.75];

        // Act
        run_gc_bias(&none_cfg)?;
        run_gc_bias(&quantile_cfg)?;

        // Assert
        let package_none =
            GCCorrectionPackage::from_file(out_none.path().join("gc_bias_correction.npz"))?;
        let package_quantile =
            GCCorrectionPackage::from_file(out_quantile.path().join("gc_bias_correction.npz"))?;

        assert_eq!(package_none.correction_matrix.dim(), (2, 2));
        assert_eq!(package_quantile.correction_matrix.dim(), (2, 2));
        assert_eq!(package_none.gc_edges, vec![0, 1, 100]);
        assert_eq!(package_quantile.gc_edges, vec![0, 1, 100]);
        assert_eq!(package_none.length_bin_frequencies.len(), 2);
        assert_eq!(package_quantile.length_bin_frequencies.len(), 2);
        assert!((package_none.length_bin_frequencies[0] - 0.5).abs() < 1e-12);
        assert!((package_none.length_bin_frequencies[1] - 0.5).abs() < 1e-12);
        assert!((package_quantile.length_bin_frequencies[0] - 0.5).abs() < 1e-12);
        assert!((package_quantile.length_bin_frequencies[1] - 0.5).abs() < 1e-12);

        assert!((package_none.correction_matrix[(0, 0)] - 5.0).abs() < 1e-12);
        assert!((package_none.correction_matrix[(0, 1)] - (5.0 / 9.0)).abs() < 1e-12);
        assert!((package_none.correction_matrix[(1, 0)] - 1.0).abs() < 1e-12);
        assert!((package_none.correction_matrix[(1, 1)] - 1.0).abs() < 1e-12);

        assert_gc_command_close(
            package_quantile.correction_matrix[(0, 0)],
            5.0 / 3.0,
            "quantile row0 col0",
        );
        assert_gc_command_close(
            package_quantile.correction_matrix[(0, 1)],
            5.0 / 7.0,
            "quantile row0 col1",
        );
        assert_gc_command_close(
            package_quantile.correction_matrix[(1, 0)],
            1.0,
            "quantile row1 col0",
        );
        assert_gc_command_close(
            package_quantile.correction_matrix[(1, 1)],
            1.0,
            "quantile row1 col1",
        );

        Ok(())
    }

    #[test]
    fn quantile_outlier_scope_global_differs_from_per_length_in_real_command() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Reuse the same raw correction matrix derivation as the previous test:
        //   [[0.2, 1.8],
        //    [1.0, 1.0]]
        //
        // With `quantile` and explicit `Q25/Q75 = 0.25/0.75`:
        //
        // Per-length scope:
        // - length 10 row [0.2, 1.8] -> [0.6, 1.4]
        // - length 11 row [1.0, 1.0] -> [1.0, 1.0]
        // - final weights:
        //     [[5/3, 5/7],
        //      [1.0, 1.0]]
        //
        // Global scope:
        // - full sorted matrix values are [0.2, 1.0, 1.0, 1.8]
        // - linear interpolation gives:
        //     Q25 = 0.2 + 0.75 * (1.0 - 0.2) = 0.8
        //     Q75 = 1.0 + 0.25 * (1.8 - 1.0) = 1.2
        // - winsorized matrix becomes:
        //     [[0.8, 1.2],
        //      [1.0, 1.0]]
        // - final weights are therefore:
        //     [[1.25, 5/6],
        //      [1.0, 1.0]]
        let (reference, bam) = make_two_length_outlier_fixture()?;
        let ref_gc_dir = TempDir::new()?;
        write_balanced_two_length_reference_gc_package(ref_gc_dir.path())?;

        let out_per_length = TempDir::new()?;
        let out_global = TempDir::new()?;

        let make_cfg = |output_dir: &std::path::Path| {
            let mut cfg =
                make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), output_dir);
            cfg.set_windows(GCWindowsArgs {
                by_size: None,
                by_bed: None,
                global: true,
            });
            cfg.set_min_length_bin_mass(0.0);
            cfg.set_min_length_bin_width(1);
            cfg.set_min_gc_bin_mass(1.0);
            cfg.set_num_extreme_gc_bins(0);
            cfg.set_num_short_length_bins(0);
            cfg.outlier_method = cfdnalab::commands::gc_bias::config::OutlierMethodArg::Quantile;
            cfg.outlier_quantiles = vec![0.25, 0.75];
            cfg
        };

        let mut per_length_cfg = make_cfg(out_per_length.path());
        per_length_cfg.outlier_scope =
            cfdnalab::commands::gc_bias::config::OutlierScopeArg::PerLength;

        let mut global_cfg = make_cfg(out_global.path());
        global_cfg.outlier_scope = cfdnalab::commands::gc_bias::config::OutlierScopeArg::Global;

        // Act
        run_gc_bias(&per_length_cfg)?;
        run_gc_bias(&global_cfg)?;

        // Assert
        let package_per_length =
            GCCorrectionPackage::from_file(out_per_length.path().join("gc_bias_correction.npz"))?;
        let package_global =
            GCCorrectionPackage::from_file(out_global.path().join("gc_bias_correction.npz"))?;

        assert_eq!(package_per_length.correction_matrix.dim(), (2, 2));
        assert_eq!(package_global.correction_matrix.dim(), (2, 2));

        assert_gc_command_close(
            package_per_length.correction_matrix[(0, 0)],
            5.0 / 3.0,
            "per-length quantile row0 col0",
        );
        assert_gc_command_close(
            package_per_length.correction_matrix[(0, 1)],
            5.0 / 7.0,
            "per-length quantile row0 col1",
        );
        assert_gc_command_close(
            package_per_length.correction_matrix[(1, 0)],
            1.0,
            "per-length quantile row1 col0",
        );
        assert_gc_command_close(
            package_per_length.correction_matrix[(1, 1)],
            1.0,
            "per-length quantile row1 col1",
        );

        assert_gc_command_close(
            package_global.correction_matrix[(0, 0)],
            1.25,
            "global quantile row0 col0",
        );
        assert_gc_command_close(
            package_global.correction_matrix[(0, 1)],
            5.0 / 6.0,
            "global quantile row0 col1",
        );
        assert_gc_command_close(
            package_global.correction_matrix[(1, 0)],
            1.0,
            "global quantile row1 col0",
        );
        assert_gc_command_close(
            package_global.correction_matrix[(1, 1)],
            1.0,
            "global quantile row1 col1",
        );

        Ok(())
    }

    #[test]
    fn iqr_outlier_method_changes_real_command_correction_matrix_in_expected_way() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Reuse the same raw correction matrix as the quantile tests:
        //   [[0.2, 1.8],
        //    [1.0, 1.0]]
        //
        // We now use `--outlier-method iqr --outlier-scope per-length --outlier-k 0.25`.
        //
        // For the skewed length-10 row [0.2, 1.8]:
        // - Q1 = 0.6
        // - Q3 = 1.4
        // - IQR = 0.8
        // - Tukey bounds with k=0.25 are:
        //     lower = 0.6 - 0.25 * 0.8 = 0.4
        //     upper = 1.4 + 0.25 * 0.8 = 1.6
        // - the row is winsorized to [0.4, 1.6]
        //
        // The balanced length-11 row [1.0, 1.0] has IQR=0 and therefore stays [1.0, 1.0].
        //
        // After the command's final inversion step, the package must store:
        //   [[1 / 0.4, 1 / 1.6],
        //    [1 / 1.0, 1 / 1.0]]
        // =
        //   [[2.5, 0.625],
        //    [1.0, 1.0]]
        let (reference, bam) = make_two_length_outlier_fixture()?;
        let ref_gc_dir = TempDir::new()?;
        write_balanced_two_length_reference_gc_package(ref_gc_dir.path())?;

        let out_dir = TempDir::new()?;
        let mut cfg =
            make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());
        cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });
        cfg.set_min_length_bin_mass(0.0);
        cfg.set_min_length_bin_width(1);
        cfg.set_min_gc_bin_mass(1.0);
        cfg.set_num_extreme_gc_bins(0);
        cfg.set_num_short_length_bins(0);
        cfg.outlier_method = cfdnalab::commands::gc_bias::config::OutlierMethodArg::Iqr;
        cfg.outlier_scope = cfdnalab::commands::gc_bias::config::OutlierScopeArg::PerLength;
        cfg.outlier_k = 0.25;

        // Act
        run_gc_bias(&cfg)?;

        // Assert
        let package =
            GCCorrectionPackage::from_file(out_dir.path().join("gc_bias_correction.npz"))?;
        assert_eq!(package.correction_matrix.dim(), (2, 2));
        assert_eq!(package.gc_edges, vec![0, 1, 100]);

        assert_gc_command_close(package.correction_matrix[(0, 0)], 2.5, "iqr row0 col0");
        assert_gc_command_close(package.correction_matrix[(0, 1)], 0.625, "iqr row0 col1");
        assert_gc_command_close(package.correction_matrix[(1, 0)], 1.0, "iqr row1 col0");
        assert_gc_command_close(package.correction_matrix[(1, 1)], 1.0, "iqr row1 col1");

        Ok(())
    }

    #[test]
    fn stddev_outlier_method_changes_real_command_correction_matrix_in_expected_way() -> Result<()>
    {
        // Human verification status: unverified
        // Arrange:
        // Reuse the same raw correction matrix as the other run-level outlier tests:
        //   [[0.2, 1.8],
        //    [1.0, 1.0]]
        //
        // We now use `--outlier-method stddev --outlier-scope per-length --outlier-k 0.6`.
        //
        // For the skewed length-10 row [0.2, 1.8]:
        // - mean = (0.2 + 1.8) / 2 = 1.0
        // - sd   = sqrt((0.8^2 + 0.8^2) / 2) = 0.8
        // - bounds are:
        //     lower = 1.0 - 0.6 * 0.8 = 0.52 = 13/25
        //     upper = 1.0 + 0.6 * 0.8 = 1.48 = 37/25
        // - winsorized row becomes [0.52, 1.48]
        //
        // The balanced length-11 row [1.0, 1.0] has sd=0 and therefore stays [1.0, 1.0].
        //
        // After the command's final inversion step, the package must store:
        //   [[1 / (13/25), 1 / (37/25)],
        //    [1 / 1.0,     1 / 1.0]]
        // =
        //   [[25/13, 25/37],
        //    [1.0,   1.0]]
        let (reference, bam) = make_two_length_outlier_fixture()?;
        let ref_gc_dir = TempDir::new()?;
        write_balanced_two_length_reference_gc_package(ref_gc_dir.path())?;

        let out_dir = TempDir::new()?;
        let mut cfg =
            make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());
        cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });
        cfg.set_min_length_bin_mass(0.0);
        cfg.set_min_length_bin_width(1);
        cfg.set_min_gc_bin_mass(1.0);
        cfg.set_num_extreme_gc_bins(0);
        cfg.set_num_short_length_bins(0);
        cfg.outlier_method = cfdnalab::commands::gc_bias::config::OutlierMethodArg::Stddev;
        cfg.outlier_scope = cfdnalab::commands::gc_bias::config::OutlierScopeArg::PerLength;
        cfg.outlier_k = 0.6;

        // Act
        run_gc_bias(&cfg)?;

        // Assert
        let package =
            GCCorrectionPackage::from_file(out_dir.path().join("gc_bias_correction.npz"))?;
        assert_eq!(package.correction_matrix.dim(), (2, 2));
        assert_eq!(package.gc_edges, vec![0, 1, 100]);

        assert_gc_command_close(
            package.correction_matrix[(0, 0)],
            25.0 / 13.0,
            "stddev row0 col0",
        );
        assert_gc_command_close(
            package.correction_matrix[(0, 1)],
            25.0 / 37.0,
            "stddev row0 col1",
        );
        assert_gc_command_close(package.correction_matrix[(1, 0)], 1.0, "stddev row1 col0");
        assert_gc_command_close(package.correction_matrix[(1, 1)], 1.0, "stddev row1 col1");

        Ok(())
    }

    #[test]
    fn mad_outlier_method_changes_real_command_correction_matrix_in_expected_way() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Reuse the same raw correction matrix as the other run-level outlier tests:
        //   [[0.2, 1.8],
        //    [1.0, 1.0]]
        //
        // We now use `--outlier-method mad --outlier-scope per-length --outlier-k 0.5`.
        //
        // For the skewed length-10 row [0.2, 1.8]:
        // - median = 1.0
        // - absolute deviations are [0.8, 0.8]
        // - the implementation scales MAD by 1.4826, so:
        //     scaled_mad = 1.4826 * 0.8 = 1.18608
        // - bounds are:
        //     lower = 1.0 - 0.5 * 1.18608 = 0.40696
        //     upper = 1.0 + 0.5 * 1.18608 = 1.59304
        // - winsorized row becomes [0.40696, 1.59304]
        //
        // The balanced length-11 row [1.0, 1.0] has zero MAD and therefore stays [1.0, 1.0].
        //
        // After the command's final inversion step, the package must store:
        //   [[1 / 0.40696, 1 / 1.59304],
        //    [1 / 1.0,     1 / 1.0]]
        // =
        //   [[2.457244200903038, 0.627730626600084],
        //    [1.0,               1.0]]
        let (reference, bam) = make_two_length_outlier_fixture()?;
        let ref_gc_dir = TempDir::new()?;
        write_balanced_two_length_reference_gc_package(ref_gc_dir.path())?;

        let out_dir = TempDir::new()?;
        let mut cfg =
            make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());
        cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });
        cfg.set_min_length_bin_mass(0.0);
        cfg.set_min_length_bin_width(1);
        cfg.set_min_gc_bin_mass(1.0);
        cfg.set_num_extreme_gc_bins(0);
        cfg.set_num_short_length_bins(0);
        cfg.outlier_method = cfdnalab::commands::gc_bias::config::OutlierMethodArg::Mad;
        cfg.outlier_scope = cfdnalab::commands::gc_bias::config::OutlierScopeArg::PerLength;
        cfg.outlier_k = 0.5;

        // Act
        run_gc_bias(&cfg)?;

        // Assert
        let package =
            GCCorrectionPackage::from_file(out_dir.path().join("gc_bias_correction.npz"))?;
        assert_eq!(package.correction_matrix.dim(), (2, 2));
        assert_eq!(package.gc_edges, vec![0, 1, 100]);

        assert_gc_command_close(
            package.correction_matrix[(0, 0)],
            2.457244200903038,
            "mad row0 col0",
        );
        assert_gc_command_close(
            package.correction_matrix[(0, 1)],
            0.627730626600084,
            "mad row0 col1",
        );
        assert_gc_command_close(package.correction_matrix[(1, 0)], 1.0, "mad row1 col0");
        assert_gc_command_close(package.correction_matrix[(1, 1)], 1.0, "mad row1 col1");

        Ok(())
    }

    #[test]
    fn hard_clamp_changes_real_command_correction_matrix_in_expected_way() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Use a hand-built reference package that supports only GC% 0 and GC% 100 for length 10.
        // That keeps the mean-scaling denominator restricted to the two relevant cells.
        //
        // Sample BAM:
        // - one fragment in the A-only region  -> GC%=0
        // - 999 fragments in the C-only region -> GC%=100
        //
        // Before GC binning, global raw counts on the supported cells are [1, 999], so:
        // - supported mean = (1 + 999) / 2 = 500
        // - normalized supported counts = [0.002, 1.998]
        //
        // With `min_gc_bin_mass = 0.09%`, the total normalized mass is exactly 2.0, so:
        // - min bin mass = 0.0018
        // - the first GC bin `[0]` already exceeds threshold on its own
        // - the second bin is `[1..100]`
        //
        // The reference package is balanced across those same two bins, so after reference-side
        // normalization the raw correction row is:
        //   [0.002, 1.998]
        //
        // Outlier handling is disabled, so only the hard safety clamp applies:
        // - low cell 0.002 is clamped up to 0.1
        // - high cell 1.998 is unchanged
        // giving:
        //   [0.1, 1.998]
        //
        // The command then re-normalizes the row to mean 1.0 before inversion:
        // - mean = (0.1 + 1.998) / 2 = 1.049
        // - normalized row = [0.1 / 1.049, 1.998 / 1.049]
        // - final multiplicative weights after inversion are:
        //     [1.049 / 0.1, 1.049 / 1.998]
        //   = [10.49, 0.525025025025025]
        //
        // This is the important contract: the hard clamp happens *before* the final re-centering,
        // so the written package can end up slightly outside the nominal clamp range afterwards.
        let reference = fixtures::twobit_from_sequences(
            "gc_bias_hard_clamp_reference",
            vec![(
                "chr1".to_string(),
                format!("{}{}", "A".repeat(100), "C".repeat(100)),
            )],
        )?;

        let mut fragments = Vec::new();
        fragments.push(fixtures::paired_fragment(10, 10, 5));
        for _ in 0..999 {
            fragments.push(fixtures::paired_fragment(120, 10, 5));
        }
        // The hard-clamp setup stacks 999 fragments at the same C-only start, so it must opt
        // into unique qnames to keep those molecules distinct in the paired-end parser.
        let bam = fixtures::bam_from_specs_strict_identity(
            vec![("chr1".to_string(), 200)],
            fragments,
            Vec::new(),
            "gc_bias_hard_clamp_bam",
        )?;

        let ref_gc_dir = TempDir::new()?;
        write_two_bin_reference_gc_package(ref_gc_dir.path(), (10, 10))?;

        let out_dir = TempDir::new()?;
        let mut cfg =
            make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());
        cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });
        cfg.set_min_length_bin_mass(0.0);
        cfg.set_min_length_bin_width(1);
        // Stay slightly below the exact 0.002 boundary so the GC bin split is not sensitive to
        // floating-point equality at the threshold.
        cfg.set_min_gc_bin_mass(0.09);
        cfg.set_num_extreme_gc_bins(0);
        cfg.set_num_short_length_bins(0);
        cfg.outlier_method = cfdnalab::commands::gc_bias::config::OutlierMethodArg::None;

        // Act
        run_gc_bias(&cfg)?;

        // Assert
        let package =
            GCCorrectionPackage::from_file(out_dir.path().join("gc_bias_correction.npz"))?;
        assert_eq!(package.correction_matrix.dim(), (1, 2));
        assert_eq!(package.gc_edges, vec![0, 1, 100]);
        assert!((package.correction_matrix[(0, 0)] - 10.49).abs() < 1e-12);
        assert!((package.correction_matrix[(0, 1)] - 0.525025025025025).abs() < 1e-12);

        Ok(())
    }

    #[test]
    fn min_length_bin_width_merges_two_lengths_into_one_binned_correction_row() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Reuse the same two-length fixture as the run-level outlier tests. Before any length
        // binning, the normalized cfDNA rows are exactly:
        //   length 10 -> [0.2, 1.8]
        //   length 11 -> [1.0, 1.0]
        //
        // and the balanced handcrafted reference package gives:
        //   reference rows -> [1.0, 1.0] for both lengths
        //
        // With:
        // - `min_length_bin_mass = 0.0`
        // - `min_length_bin_width = 2`
        // the greedy length binning must merge the two adjacent lengths into one bin, because a
        // bin cannot close until it has width >= 2.
        //
        // The merged binned rows are then the simple mean of the two source rows:
        //   merged cfDNA row = ([0.2, 1.8] + [1.0, 1.0]) / 2 = [0.6, 1.4]
        //   merged ref row   = ([1.0, 1.0] + [1.0, 1.0]) / 2 = [1.0, 1.0]
        //
        // So the raw merged correction row is [0.6, 1.4]. Its mean is already 1.0, and with
        // outlier handling disabled no other transform changes it before inversion.
        //
        // The final written correction row must therefore be:
        //   [1 / 0.6, 1 / 1.4] = [5/3, 5/7]
        let (reference, bam) = make_two_length_outlier_fixture()?;
        let ref_gc_dir = TempDir::new()?;
        write_balanced_two_length_reference_gc_package(ref_gc_dir.path())?;

        let out_dir = TempDir::new()?;
        let mut cfg =
            make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());
        cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });
        cfg.set_min_length_bin_mass(0.0);
        cfg.set_min_length_bin_width(2);
        cfg.set_min_gc_bin_mass(1.0);
        cfg.set_num_extreme_gc_bins(0);
        cfg.set_num_short_length_bins(0);
        cfg.outlier_method = cfdnalab::commands::gc_bias::config::OutlierMethodArg::None;

        // Act
        run_gc_bias(&cfg)?;

        // Assert
        let package =
            GCCorrectionPackage::from_file(out_dir.path().join("gc_bias_correction.npz"))?;
        assert_eq!(package.correction_matrix.dim(), (1, 2));
        assert_eq!(package.length_edges, vec![10, 11]);
        assert_eq!(package.gc_edges, vec![0, 1, 100]);
        assert_eq!(package.length_bin_frequencies.len(), 1);
        assert!((package.length_bin_frequencies[0] - 1.0).abs() < 1e-12);
        assert!((package.correction_matrix[(0, 0)] - (5.0 / 3.0)).abs() < 1e-12);
        assert!((package.correction_matrix[(0, 1)] - (5.0 / 7.0)).abs() < 1e-12);

        Ok(())
    }

    #[test]
    fn num_short_length_bins_neutralizes_the_shortest_length_row_in_real_command() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Reuse the same two-length fixture as the other run-level `gc-bias` tests.
        //
        // Baseline, with no short-length masking:
        //   length 10 -> raw correction [0.2, 1.8] -> final weights [5, 5/9]
        //   length 11 -> raw correction [1.0, 1.0] -> final weights [1, 1]
        //
        // Now set `num_short_length_bins = 1`. After length binning there are exactly two
        // length rows, so the shortest row is the entire length-10 row.
        //
        // The support mask then becomes:
        //   row 0 (length 10): unsupported everywhere
        //   row 1 (length 11): supported everywhere
        //
        // The pipeline contract is:
        // 1. Unsupported entries in the normalized cfDNA and reference matrices are set to 1.0.
        // 2. The raw correction row for the masked shortest length therefore becomes [1, 1].
        // 3. Re-centering and inversion keep [1, 1] unchanged.
        //
        // So the shortest row must change from the informative baseline `[5, 5/9]` to the
        // neutral row `[1, 1]`, while the longer row stays `[1, 1]`.
        let (reference, bam) = make_two_length_outlier_fixture()?;
        let ref_gc_dir = TempDir::new()?;
        write_balanced_two_length_reference_gc_package(ref_gc_dir.path())?;

        let baseline_out = TempDir::new()?;
        let masked_out = TempDir::new()?;

        let mut baseline_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            baseline_out.path(),
        );
        baseline_cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });
        baseline_cfg.set_min_length_bin_mass(0.0);
        baseline_cfg.set_min_length_bin_width(1);
        baseline_cfg.set_min_gc_bin_mass(1.0);
        baseline_cfg.set_num_extreme_gc_bins(0);
        baseline_cfg.set_num_short_length_bins(0);
        baseline_cfg.outlier_method = cfdnalab::commands::gc_bias::config::OutlierMethodArg::None;

        let mut masked_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            masked_out.path(),
        );
        masked_cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });
        masked_cfg.set_min_length_bin_mass(0.0);
        masked_cfg.set_min_length_bin_width(1);
        masked_cfg.set_min_gc_bin_mass(1.0);
        masked_cfg.set_num_extreme_gc_bins(0);
        masked_cfg.set_num_short_length_bins(1);
        masked_cfg.outlier_method = cfdnalab::commands::gc_bias::config::OutlierMethodArg::None;

        // Act
        run_gc_bias(&baseline_cfg)?;
        run_gc_bias(&masked_cfg)?;

        // Assert
        let baseline_package =
            GCCorrectionPackage::from_file(baseline_out.path().join("gc_bias_correction.npz"))?;
        let masked_package =
            GCCorrectionPackage::from_file(masked_out.path().join("gc_bias_correction.npz"))?;

        assert_eq!(baseline_package.correction_matrix.dim(), (2, 2));
        assert_eq!(masked_package.correction_matrix.dim(), (2, 2));
        assert_eq!(baseline_package.gc_edges, vec![0, 1, 100]);
        assert_eq!(masked_package.gc_edges, vec![0, 1, 100]);

        assert!((baseline_package.correction_matrix[(0, 0)] - 5.0).abs() < 1e-12);
        assert!((baseline_package.correction_matrix[(0, 1)] - (5.0 / 9.0)).abs() < 1e-12);
        assert!((baseline_package.correction_matrix[(1, 0)] - 1.0).abs() < 1e-12);
        assert!((baseline_package.correction_matrix[(1, 1)] - 1.0).abs() < 1e-12);

        assert!((masked_package.correction_matrix[(0, 0)] - 1.0).abs() < 1e-12);
        assert!((masked_package.correction_matrix[(0, 1)] - 1.0).abs() < 1e-12);
        assert!((masked_package.correction_matrix[(1, 0)] - 1.0).abs() < 1e-12);
        assert!((masked_package.correction_matrix[(1, 1)] - 1.0).abs() < 1e-12);

        Ok(())
    }

    #[test]
    fn num_extreme_gc_bins_neutralizes_a_two_bin_gc_axis_in_real_command() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Reuse the same two-length fixture again. With `min_gc_bin_mass = 1.0`, the sparse
        // counts collapse to exactly two GC bins:
        //   left  bin = GC% 0
        //   right bin = GC% 1..100
        //
        // Baseline, with `num_extreme_gc_bins = 0`:
        //   length 10 -> [5, 5/9]
        //   length 11 -> [1, 1]
        //
        // Now set `num_extreme_gc_bins = 1`. On a 2-bin GC axis, "one extreme bin from each side"
        // masks both columns:
        //   leftmost column  -> masked
        //   rightmost column -> masked
        //
        // So every correction cell is unsupported. The pipeline then sets all masked normalized
        // cfDNA/reference counts to 1.0 before division, yielding a raw correction matrix of all
        // ones. Re-centering and inversion keep it at all ones.
        //
        // This is an important boundary contract: on a 2-bin GC axis, one extreme bin per side
        // completely neutralizes the matrix rather than leaving any informative GC correction.
        let (reference, bam) = make_two_length_outlier_fixture()?;
        let ref_gc_dir = TempDir::new()?;
        write_balanced_two_length_reference_gc_package(ref_gc_dir.path())?;

        let baseline_out = TempDir::new()?;
        let masked_out = TempDir::new()?;

        let mut baseline_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            baseline_out.path(),
        );
        baseline_cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });
        baseline_cfg.set_min_length_bin_mass(0.0);
        baseline_cfg.set_min_length_bin_width(1);
        baseline_cfg.set_min_gc_bin_mass(1.0);
        baseline_cfg.set_num_extreme_gc_bins(0);
        baseline_cfg.set_num_short_length_bins(0);
        baseline_cfg.outlier_method = cfdnalab::commands::gc_bias::config::OutlierMethodArg::None;

        let mut masked_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            masked_out.path(),
        );
        masked_cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });
        masked_cfg.set_min_length_bin_mass(0.0);
        masked_cfg.set_min_length_bin_width(1);
        masked_cfg.set_min_gc_bin_mass(1.0);
        masked_cfg.set_num_extreme_gc_bins(1);
        masked_cfg.set_num_short_length_bins(0);
        masked_cfg.outlier_method = cfdnalab::commands::gc_bias::config::OutlierMethodArg::None;

        // Act
        run_gc_bias(&baseline_cfg)?;
        run_gc_bias(&masked_cfg)?;

        // Assert
        let baseline_package =
            GCCorrectionPackage::from_file(baseline_out.path().join("gc_bias_correction.npz"))?;
        let masked_package =
            GCCorrectionPackage::from_file(masked_out.path().join("gc_bias_correction.npz"))?;

        assert_eq!(baseline_package.correction_matrix.dim(), (2, 2));
        assert_eq!(masked_package.correction_matrix.dim(), (2, 2));
        assert_eq!(baseline_package.gc_edges, vec![0, 1, 100]);
        assert_eq!(masked_package.gc_edges, vec![0, 1, 100]);

        assert!((baseline_package.correction_matrix[(0, 0)] - 5.0).abs() < 1e-12);
        assert!((baseline_package.correction_matrix[(0, 1)] - (5.0 / 9.0)).abs() < 1e-12);
        assert!((baseline_package.correction_matrix[(1, 0)] - 1.0).abs() < 1e-12);
        assert!((baseline_package.correction_matrix[(1, 1)] - 1.0).abs() < 1e-12);

        for value in masked_package.correction_matrix.iter() {
            assert!((*value - 1.0).abs() < 1e-12);
        }

        Ok(())
    }

    #[test]
    fn min_length_bin_mass_merges_a_sparse_tail_length_into_the_previous_bin() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Use a two-length fixture where the shorter length has mass 10 and the longer tail length
        // has mass 2.
        //
        // Baseline, with `min_length_bin_mass = 0` and width 1:
        //   length 10 -> [0.2, 1.8] -> final weights [5, 5/9]
        //   length 11 -> [1.0, 1.0] -> final weights [1, 1]
        //
        // Now set `min_length_bin_mass = 20%`.
        // Total binned mass is 12, so the minimum bin mass is:
        //   12 * 0.20 = 2.4
        //
        // Greedy binning over the length axis then behaves as follows:
        // - row for length 10 has mass 10, so it can close a bin by itself
        // - row for length 11 has mass 2, so it cannot form its own bin
        // - the final underweight tail row is therefore appended to the previous bin
        //
        // The important implementation detail is that greedy *length* merging happens before the
        // later per-row mean-scaling step in `gc-bias`.
        //
        // In this fixture the globally normalized GC rows are:
        //   length 10 -> [1/3, 3]
        //   length 11 -> [1/3, 1/3]
        //
        // So after greedy length merging, the single merged row is the arithmetic mean of those
        // pre-row-normalized rows:
        //   ([1/3, 3] + [1/3, 1/3]) / 2 = [1/3, 5/3]
        //
        // The reference side is still balanced after the same merge:
        //   [1, 1]
        //
        // Per-row mean scaling then keeps the raw correction row at:
        //   [1/3, 5/3]
        //
        // With no outlier handling, the final multiplicative row is therefore:
        //   [1 / (1/3), 1 / (5/3)] = [3, 3/5]
        let (reference, bam) = make_two_length_low_mass_tail_fixture()?;
        let ref_gc_dir = TempDir::new()?;
        write_balanced_two_length_reference_gc_package(ref_gc_dir.path())?;

        let baseline_out = TempDir::new()?;
        let merged_out = TempDir::new()?;

        let mut baseline_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            baseline_out.path(),
        );
        baseline_cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });
        baseline_cfg.set_min_length_bin_mass(0.0);
        baseline_cfg.set_min_length_bin_width(1);
        baseline_cfg.set_min_gc_bin_mass(1.0);
        baseline_cfg.set_num_extreme_gc_bins(0);
        baseline_cfg.set_num_short_length_bins(0);
        baseline_cfg.outlier_method = cfdnalab::commands::gc_bias::config::OutlierMethodArg::None;

        let mut merged_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            merged_out.path(),
        );
        merged_cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });
        merged_cfg.set_min_length_bin_mass(20.0);
        merged_cfg.set_min_length_bin_width(1);
        merged_cfg.set_min_gc_bin_mass(1.0);
        merged_cfg.set_num_extreme_gc_bins(0);
        merged_cfg.set_num_short_length_bins(0);
        merged_cfg.outlier_method = cfdnalab::commands::gc_bias::config::OutlierMethodArg::None;

        // Act
        run_gc_bias(&baseline_cfg)?;
        run_gc_bias(&merged_cfg)?;

        // Assert
        let baseline_package =
            GCCorrectionPackage::from_file(baseline_out.path().join("gc_bias_correction.npz"))?;
        let merged_package =
            GCCorrectionPackage::from_file(merged_out.path().join("gc_bias_correction.npz"))?;

        assert_eq!(baseline_package.correction_matrix.dim(), (2, 2));
        assert_eq!(merged_package.correction_matrix.dim(), (1, 2));
        assert_eq!(merged_package.length_edges, vec![10, 11]);
        assert_eq!(merged_package.gc_edges, vec![0, 1, 100]);
        assert_eq!(merged_package.length_bin_frequencies.len(), 1);
        assert!((merged_package.length_bin_frequencies[0] - 1.0).abs() < 1e-12);

        assert!((baseline_package.correction_matrix[(0, 0)] - 5.0).abs() < 1e-12);
        assert!((baseline_package.correction_matrix[(0, 1)] - (5.0 / 9.0)).abs() < 1e-12);
        assert!((baseline_package.correction_matrix[(1, 0)] - 1.0).abs() < 1e-12);
        assert!((baseline_package.correction_matrix[(1, 1)] - 1.0).abs() < 1e-12);

        assert!((merged_package.correction_matrix[(0, 0)] - 3.0).abs() < 1e-12);
        assert!((merged_package.correction_matrix[(0, 1)] - (3.0 / 5.0)).abs() < 1e-12);

        Ok(())
    }

    #[test]
    fn min_gc_bin_mass_greedily_merges_sparse_gc_tail_bins_in_real_command() -> Result<()> {
        // Human verification status: unverified
        // Arrange:
        // Use a one-length fixture with exact GC-class masses:
        //   GC%=0   -> 1
        //   GC%=50  -> 5
        //   GC%=100 -> 9
        //
        // The handcrafted reference package marks those same three GC points as equally likely.
        //
        // Baseline, with `min_gc_bin_mass = 1%`:
        // - total sample mass = 15
        // - min bin mass = 0.15
        // - greedy GC binning therefore closes bins at:
        //     [0], [1..50], [51..100]
        // - sample binned row is [1, 5, 9]
        // - reference binned row is [1, 1, 1]
        // - normalized correction row is already [0.2, 1.0, 1.8]
        // - final weights are therefore [5, 1, 5/9]
        //
        // Now set `min_gc_bin_mass = 25%`:
        // - min bin mass = 15 * 0.25 = 3.75
        // - the first sparse GC point (mass 1 at GC%=0) cannot close a bin by itself
        // - it gets merged with the next non-zero point at GC%=50
        // - the resulting GC bins are:
        //     [0..50], [51..100]
        //
        // So the binned rows become:
        //   sample    -> [1 + 5, 9] = [6, 9]
        //   reference -> [1 + 1, 1] = [2, 1]
        //
        // Per-row normalization gives:
        //   sample    -> [0.8, 1.2]
        //   reference -> [4/3, 2/3]
        //
        // Raw correction is therefore:
        //   [0.8 / (4/3), 1.2 / (2/3)] = [0.6, 1.8]
        //
        // Re-centering that row to mean 1.0 gives [0.5, 1.5], and inversion yields:
        //   [2, 2/3]
        let (reference, bam) = make_three_gc_bin_fixture()?;
        let ref_gc_dir = TempDir::new()?;
        write_three_bin_reference_gc_package(ref_gc_dir.path())?;

        let baseline_out = TempDir::new()?;
        let merged_out = TempDir::new()?;

        let mut baseline_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            baseline_out.path(),
        );
        baseline_cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });
        baseline_cfg.set_min_length_bin_mass(0.0);
        baseline_cfg.set_min_length_bin_width(1);
        baseline_cfg.set_min_gc_bin_mass(1.0);
        baseline_cfg.set_num_extreme_gc_bins(0);
        baseline_cfg.set_num_short_length_bins(0);
        baseline_cfg.outlier_method = cfdnalab::commands::gc_bias::config::OutlierMethodArg::None;

        let mut merged_cfg = make_gc_bias_cfg(
            &bam.bam,
            &reference.path,
            ref_gc_dir.path(),
            merged_out.path(),
        );
        merged_cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: None,
            global: true,
        });
        merged_cfg.set_min_length_bin_mass(0.0);
        merged_cfg.set_min_length_bin_width(1);
        merged_cfg.set_min_gc_bin_mass(25.0);
        merged_cfg.set_num_extreme_gc_bins(0);
        merged_cfg.set_num_short_length_bins(0);
        merged_cfg.outlier_method = cfdnalab::commands::gc_bias::config::OutlierMethodArg::None;

        // Act
        run_gc_bias(&baseline_cfg)?;
        run_gc_bias(&merged_cfg)?;

        // Assert
        let baseline_package =
            GCCorrectionPackage::from_file(baseline_out.path().join("gc_bias_correction.npz"))?;
        let merged_package =
            GCCorrectionPackage::from_file(merged_out.path().join("gc_bias_correction.npz"))?;

        assert_eq!(baseline_package.correction_matrix.dim(), (1, 3));
        assert_eq!(baseline_package.gc_edges, vec![0, 1, 51, 100]);
        assert!((baseline_package.correction_matrix[(0, 0)] - 5.0).abs() < 1e-12);
        assert!((baseline_package.correction_matrix[(0, 1)] - 1.0).abs() < 1e-12);
        assert!((baseline_package.correction_matrix[(0, 2)] - (5.0 / 9.0)).abs() < 1e-12);

        assert_eq!(merged_package.correction_matrix.dim(), (1, 2));
        assert_eq!(merged_package.gc_edges, vec![0, 51, 100]);
        assert!((merged_package.correction_matrix[(0, 0)] - 2.0).abs() < 1e-12);
        assert!((merged_package.correction_matrix[(0, 1)] - (2.0 / 3.0)).abs() < 1e-12);

        Ok(())
    }
}

mod tests_counts_end_offset {
    use cfdnalab::commands::gc_bias::counting::GCCounts;

    #[test]
    fn should_use_effective_length_when_binning_to_gc_percent_with_end_offset() {
        // Human verification status: unverified
        // Arrange: one 30bp fragment with 20 GC bases after trimming 5bp from each end
        let mut counts = GCCounts::new(30, 30, 5, (0, 0)).expect("counts init");
        counts.incr(30, 20);

        // Act
        let grid = counts.to_gc_percent_grid(0, 100).expect("gc percent grid");

        // Assert: value lands in the 100% bin, not in the 67% bin (which used full length)
        assert_eq!(grid[(0, 100)], 1.0);
        assert_eq!(grid[(0, 67)], 0.0);
    }

    #[test]
    fn should_not_smooth_into_gc_counts_beyond_effective_length() {
        // Human verification status: unverified
        // Arrange: length=6, end_offset=2 -> effective length is 2bp, so gc>2 is unreachable.
        let mut counts = GCCounts::new(6, 6, 2, (0, 0)).expect("counts init");
        counts.set(6, 2, 10.0);

        // Act: smooth only the reachable portion of the row.
        counts.smooth_length_rows_in_place(1.0, 1);

        // Assert: unreachable GC counts are absent and storage matches the effective length.
        assert!(counts.get(6, 3).is_none());
        assert_eq!(counts.borrow_raw_counts().len(), 3);
    }
}

mod tests_gc_percent_grid {

    use cfdnalab::commands::gc_bias::counting::GCCounts;

    #[test]
    fn should_place_gc_counts_in_matching_percent_bins() {
        // Human verification status: unverified
        // Arrange: one length row with distinct weights per GC count.
        let mut counts = GCCounts::new(10, 10, 0, (0, 0)).expect("counts init");
        for gc in 0..=10 {
            counts.set(10, gc, (gc + 1) as f64); // unique weight per bin
        }

        // Act
        let grid = counts.to_gc_percent_grid(0, 100).expect("gc percent grid");
        let row = grid.row(0);

        // Assert: each GC count lands in its integer percent bin.
        for gc in 0..=10 {
            let pct_bin = (gc * 10) as usize;
            assert!(
                (row[pct_bin] - (gc + 1) as f64).abs() < 1e-12,
                "gc {} expected at pct {}, got {}",
                gc,
                pct_bin,
                row[pct_bin]
            );
        }
    }

    #[test]
    fn should_round_half_up_for_fractional_percentages() {
        // Human verification status: unverified
        // Arrange: length=3 has fractional percentages for gc=1 and gc=2.
        let mut counts = GCCounts::new(3, 3, 0, (0, 0)).expect("counts init");
        counts.set(3, 1, 2.0); // 33.3...% -> 33 via half-up
        counts.set(3, 2, 3.0); // 66.6...% -> 67 via half-up

        // Act
        let grid = counts.to_gc_percent_grid(0, 100).expect("gc percent grid");
        let row = grid.row(0);

        // Assert: derive the half-up bins explicitly
        // calculate_gc_bin does round_half_up(100 * gc / effective_length) via (100 * gc + len/2) / len
        // Effective length is 3 (no end trimming)
        // gc=1 -> (100 * 1 + 3/2) / 3 = (100 + 1) / 3 = 33
        // gc=2 -> (100 * 2 + 3/2) / 3 = (200 + 1) / 3 = 67
        // Mass must land only in those bins
        for (idx, &val) in row.iter().enumerate() {
            match idx {
                33 => assert!(
                    (val - 2.0).abs() < 1e-12,
                    "bin 33 expected 2.0, got {}",
                    val
                ),
                67 => assert!(
                    (val - 3.0).abs() < 1e-12,
                    "bin 67 expected 3.0, got {}",
                    val
                ),
                _ => assert!(val.abs() < 1e-12, "bin {} expected 0, got {}", idx, val),
            }
        }
    }

    #[test]
    fn should_propagate_acgt_totals_and_length_metadata() {
        // Human verification status: unverified
        // Arrange
        let mut counts = GCCounts::new(5, 6, 1, (8, 12)).expect("counts init");
        counts.set(5, 2, 1.0);
        counts.set(6, 3, 2.0);

        // Act
        let grid = counts.to_gc_percent_grid(0, 100).expect("gc percent grid");

        // Assert: shapes match the two length bins and 101 GC bins
        assert_eq!(grid.nrows(), 2);
        assert_eq!(grid.ncols(), 101);

        let row_len5 = grid.row(0);
        let row_len6 = grid.row(1);
        // Derivation with end offsets
        // End offset is 1 so effective length = length - 2
        // calculate_gc_bin uses (100 * gc + eff_len/2) / eff_len
        // len5 -> eff3: gc=2 gives (100 * 2 + 3/2) / 3 = (200 + 1) / 3 = 67
        // len6 -> eff4: gc=3 gives (100 * 3 + 4/2) / 4 = (300 + 2) / 4 = 75
        assert!((row_len5[67] - 1.0).abs() < 1e-12);
        assert!((row_len6[75] - 2.0).abs() < 1e-12);
    }
}

mod tests_length_bounds {

    use cfdnalab::commands::gc_bias::counting::GCCounts;

    #[test]
    fn reports_offsets_based_on_effective_length() {
        // Human verification status: unverified
        // length_min=3, length_max=5, end_offset=1 -> effective lengths: 1,2,3
        let counts = GCCounts::new(3, 5, 1, (0, 0)).expect("init counts");

        let bounds_len3 = counts.length_bounds(3).expect("len3 bounds");
        let bounds_len4 = counts.length_bounds(4).expect("len4 bounds");
        let bounds_len5 = counts.length_bounds(5).expect("len5 bounds");

        assert_eq!(bounds_len3, (0, 2)); // size 2 for effective len 1 (gc 0..1)
        assert_eq!(bounds_len4, (2, 5)); // size 3 for effective len 2 (gc 0..2)
        assert_eq!(bounds_len5, (5, 9)); // size 4 for effective len 3 (gc 0..3)

        // Verify the slice lengths match the effective length + 1
        assert_eq!(bounds_len3.1 - bounds_len3.0, 2);
        assert_eq!(bounds_len4.1 - bounds_len4.0, 3);
        assert_eq!(bounds_len5.1 - bounds_len5.0, 4);
    }

    #[test]
    fn row_bounds_errors_outside_length_range() {
        // Human verification status: unverified
        let counts = GCCounts::new(10, 12, 0, (0, 0)).expect("init counts");
        assert!(counts.length_bounds(9).is_err());
        assert!(counts.length_bounds(13).is_err());
    }
}

mod tests_helpers {

    #[cfg(test)]
    mod tests {
        use cfdnalab::commands::gc_bias::gc_bias::mean_scale_per_length_array;
        use ndarray::array;

        #[test]
        fn leaves_zero_rows_untouched_in_mean_scaling() {
            // Human verification status: unverified
            // Arrange: first length row has no mass; second has values that should be mean-scaled.
            let counts = array![[0.0, 0.0], [2.0, 4.0]];
            let mask = array![[true, true], [true, true]];

            // Act
            let scaled = mean_scale_per_length_array(&counts, 0.0, Some(&mask));

            // Assert: empty row stays zero; non-empty row divides by its mean (3.0).
            assert!(
                scaled.row(0).iter().all(|&v| v == 0.0),
                "zero row should remain zero after scaling"
            );
            assert!((scaled[(1, 0)] - 2.0 / 3.0).abs() < 1e-12);
            assert!((scaled[(1, 1)] - 4.0 / 3.0).abs() < 1e-12);
        }
    }
}
