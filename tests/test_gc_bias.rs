#![cfg(all(feature = "cmd_gc_bias", feature = "cmd_ref_gc_bias"))]

mod fixtures;

mod tests_gc_bias {
    use crate::fixtures;
    use anyhow::Result;
    use ndarray::array;
    use ndarray_npy::read_npy;
    use tempfile::{TempDir, tempdir};

    use cfdnalab::RunOptions;
    use cfdnalab::constants::GC_CORRECTION_SCHEMA_VERSION;
    use cfdnalab::gc_bias::{
        GCCorrectionPackage, GCCorrector, GCLengthRange, LengthAgnosticGCCorrector,
        MarginalizeLengthsWeightingScheme, load_gc_corrector, load_length_agnostic_gc_corrector,
    };
    use cfdnalab::reference::twobit_contig_footprint;
    use cfdnalab::run_like_cli::{
        common::{
            ChromosomeArgs, FragmentLengthArgs, GCWindowsArgs, IOCArgs, LoggingArgs,
            Ref2BitRequiredArgs,
        },
        gc_bias::{GCConfig, OutlierMethodArg, run_gc_bias as run_gc_bias_command},
        ref_gc_bias::{
            RefGCBiasConfig, RefGCWindowsArgs, run_ref_gc_bias as run_ref_gc_bias_command,
        },
    };

    const GC_COMMAND_F64_TOL: f64 = 1e-6;

    fn run_gc_bias(config: &GCConfig) -> Result<()> {
        run_gc_bias_command(config, RunOptions::new_quiet()).map(|_| ())
    }

    fn run_ref_gc_bias(config: &RefGCBiasConfig) -> Result<()> {
        run_ref_gc_bias_command(config, RunOptions::new_quiet()).map(|_| ())
    }

    fn assert_gc_command_close(actual: f64, expected: f64, context: &str) {
        // The outlier helpers estimate quantiles/bounds in `f32` and only then write the matrix
        // back as `f64`, so the stable contract here is "matches the hand-derived value within the
        // command's float precision", not bit-exact `f64` arithmetic on the ideal decimal.
        assert!(
            (actual - expected).abs() <= GC_COMMAND_F64_TOL,
            "{context}: expected {expected}, got {actual}"
        );
    }

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn gc_bias_late_tile_bed_window_uses_nonzero_sequence_interval_origin() -> Result<()> {
        // Manual derivation:
        // - The chromosome is 1022 bp and tile_size is 500, so the tile owning [900,961) has
        //   core [500,1000) and fetch origin 439 after the 61 bp halo.
        // - The observed fragment [900,961) lies in the all-C part of
        //   `late_origin_gc_reference_sequence()`, so its GC percentage is 100.
        // - In BED-window mode, `gc-bias` scales a single supported length-61 count by the number
        //   of reachable GC-count cells: 0..=61 gives 62 cells.
        // - Therefore the saved observed average count matrix has exactly 62.0 at GC% 100 and
        //   zero everywhere else. If the tile's non-zero sequence interval is ignored, the
        //   fragment cannot be mapped to the loaded prefix coordinates correctly.
        let reference = fixtures::twobit_from_sequences(
            "gc_bias_late_tile_origin_reference",
            vec![(
                "chr1".to_string(),
                fixtures::late_origin_gc_reference_sequence(),
            )],
        )?;
        let bam = fixtures::bam_from_specs(
            vec![("chr1".to_string(), 1_022)],
            vec![fixtures::paired_fragment(900, 61, 20)],
            Vec::new(),
            "gc_bias_late_tile_origin_bam",
        )?;
        let ref_gc_dir = TempDir::new()?;
        let bed_path = ref_gc_dir.path().join("late_window.bed");
        std::fs::write(&bed_path, "chr1\t900\t961\n")?;
        let ref_cfg = RefGCBiasConfig {
            ref_genome: Ref2BitRequiredArgs {
                ref_2bit: reference.path.clone(),
            },
            output_dir: ref_gc_dir.path().to_path_buf(),
            output_prefix: String::new(),
            n_threads: 1,
            n_positions: 962,
            seed: Some(23),
            windows: RefGCWindowsArgs {
                by_bed: Some(bed_path.clone()),
            },
            chromosomes: ChromosomeArgs {
                chromosomes: Some(vec!["chr1".to_string()]),
                chromosomes_file: None,
            },
            blacklist: None,
            fragment_lengths: FragmentLengthArgs {
                min_fragment_length: 61,
                max_fragment_length: 61,
            },
            end_offset: 0,
            skip_interpolation: true,
            smoothing_sigma: 0.55,
            smoothing_radius: 2,
            skip_smoothing: true,
            tile_size: 500,
            logging: LoggingArgs::default(),
        };

        let out_dir = TempDir::new()?;
        let mut cfg =
            make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());
        cfg.set_tile_size(500);
        cfg.set_windows(GCWindowsArgs {
            by_size: None,
            by_bed: Some(bed_path),
            global: false,
        });
        cfg.set_num_extreme_gc_bins(0);
        cfg.set_num_short_length_bins(0);
        cfg.set_min_gc_bin_mass(1.0);
        cfg.set_min_length_bin_mass(0.0);
        cfg.set_min_length_bin_width(1);
        cfg.outlier_method = OutlierMethodArg::None;

        run_ref_gc_bias(&ref_cfg)?;
        run_gc_bias(&cfg)?;

        let avg_counts: ndarray::Array2<f64> =
            read_npy(out_dir.path().join("gc_bias.avg_cfdna_counts.0.npy"))?;
        assert_eq!(avg_counts.dim(), (1, 101));
        for (gc_pct, &value) in avg_counts.row(0).iter().enumerate() {
            match gc_pct {
                100 => assert!(
                    (value - 62.0).abs() < 1e-12,
                    "expected scaled observed count 62.0 at GC% 100, got {value}"
                ),
                _ => assert!(
                    value.abs() < 1e-12,
                    "expected no observed mass outside GC% 100, got bin {gc_pct}={value}"
                ),
            }
        }
        Ok(())
    }

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn provides_expected_weights_after_roundtrip() -> Result<()> {
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
            reference_contig_footprint: Vec::new(),
        };
        let tmp_dir = tempdir()?;
        let pkg_path = tmp_dir.path().join("gc_package.zarr");
        package.write_zarr(&pkg_path)?;

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

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn gc_correction_loaders_reject_reference_footprint_mismatch() -> Result<()> {
        // Human verification status: verified
        // Manual expectations:
        // - The correction package carries the footprint from `reference_a`.
        // - The loaders are asked to apply it with `reference_b`, whose chr1 length differs.
        // - Both loaders should fail before returning a usable correction matrix.
        let reference_a = fixtures::simple_reference_twobit()?;
        let reference_b = fixtures::twobit_from_sequences(
            "gc_correction_loader_reference_mismatch",
            vec![("chr1".to_string(), "ACGT".repeat(80))],
        )?;
        let package = GCCorrectionPackage {
            version: GC_CORRECTION_SCHEMA_VERSION,
            end_offset: 0,
            length_edges: vec![30, 31],
            gc_edges: vec![0, 101],
            correction_matrix: array![[1.0_f64]],
            length_bin_frequencies: array![1.0_f64],
            reference_contig_footprint: twobit_contig_footprint(&reference_a.path)?,
        };
        let tmp_dir = tempdir()?;
        let package_path = tmp_dir.path().join("gc_package.zarr");
        package.write_zarr(&package_path)?;

        let standard_error =
            load_gc_corrector(Some(&package_path), Some(&reference_b.path), 30, 30)
                .expect_err("mismatched package should fail standard GC correction loading");
        let length_agnostic_error = load_length_agnostic_gc_corrector(
            Some(&package_path),
            Some(&reference_b.path),
            &MarginalizeLengthsWeightingScheme::Equal,
            GCLengthRange::Package,
            0.0,
            30,
            30,
        )
        .expect_err("mismatched package should fail length-agnostic GC correction loading");

        for error in [standard_error, length_agnostic_error] {
            let message = error.to_string();
            assert!(
                message.contains(
                    "GC correction package was built against a different reference contig than --ref-2bit."
                ),
                "unexpected error message: {message}"
            );
        }

        Ok(())
    }

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn gc_correction_package_rejects_missing_path_before_opening() -> Result<()> {
        // Arrange: point the loader at a `.zarr` path that does not exist.
        let tmp_dir = tempdir()?;
        let missing_path = tmp_dir.path().join("missing_gc_package.zarr");

        // Act: try to load the missing package.
        let err = GCCorrectionPackage::from_file(&missing_path)
            .expect_err("missing GC correction package should fail");

        // Assert: the user gets the shared "existing .zarr directory" contract directly instead
        // of a lower-level storage error.
        let msg = err.to_string();
        assert!(
            msg.contains("must point to an existing .zarr directory"),
            "unexpected error message: {msg}"
        );
        assert!(
            msg.contains("missing_gc_package.zarr"),
            "unexpected error message: {msg}"
        );

        Ok(())
    }

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn gc_correction_package_rejects_regular_file_before_parsing() -> Result<()> {
        // Arrange: create a regular file with the wrong extension.
        let tmp_dir = tempdir()?;
        let wrong_extension_path = tmp_dir.path().join("gc_package.txt");
        std::fs::write(&wrong_extension_path, b"not a zarr store")?;

        // Act: try to load the non-directory path.
        let err = GCCorrectionPackage::from_file(&wrong_extension_path)
            .expect_err("wrong extension should fail");

        // Assert: directory validation runs before opening the Zarr store.
        let msg = err.to_string();
        assert!(
            msg.contains("must point to an existing .zarr directory"),
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
            reference_contig_footprint: Vec::new(),
        }
    }

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn length_agnostic_equal_weighting_means_rows() -> Result<()> {
        let package = make_length_agnostic_package();
        let corrector = GCCorrector::from_package(&package)?;
        let agnostic = LengthAgnosticGCCorrector::from_gc_corrector(
            &corrector,
            &MarginalizeLengthsWeightingScheme::Equal,
            GCLengthRange::Package,
            0.0,
            20,
            40,
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

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn length_agnostic_frequency_weighting_uses_frequencies() -> Result<()> {
        let package = make_length_agnostic_package();
        let corrector = GCCorrector::from_package(&package)?;
        let agnostic = LengthAgnosticGCCorrector::from_gc_corrector(
            &corrector,
            &MarginalizeLengthsWeightingScheme::Frequency,
            GCLengthRange::Package,
            0.0,
            20,
            40,
        )?;

        // Weighted average with frequencies [0.2, 0.8]
        for gc_pct in 0..50 {
            assert!(
                (agnostic.get_correction_weight(gc_pct)? - 2.6).abs() < 1e-12,
                "frequency weighting should map GC% {gc_pct} into the first GC bin"
            );
        }
        for gc_pct in 50..=100 {
            assert!(
                (agnostic.get_correction_weight(gc_pct)? - 4.4).abs() < 1e-12,
                "frequency weighting should map GC% {gc_pct} into the second GC bin"
            );
        }
        Ok(())
    }

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn length_agnostic_max_frequency_picks_most_frequent_row() -> Result<()> {
        let package = make_length_agnostic_package();
        let corrector = GCCorrector::from_package(&package)?;
        let agnostic = LengthAgnosticGCCorrector::from_gc_corrector(
            &corrector,
            &MarginalizeLengthsWeightingScheme::MaxFrequency,
            GCLengthRange::Package,
            0.0,
            20,
            40,
        )?;

        // Row with highest frequency is [3.0, 5.0]
        for gc_pct in 0..50 {
            assert!(
                (agnostic.get_correction_weight(gc_pct)? - 3.0).abs() < 1e-12,
                "max-frequency weighting should map GC% {gc_pct} into the first GC bin"
            );
        }
        for gc_pct in 50..=100 {
            assert!(
                (agnostic.get_correction_weight(gc_pct)? - 5.0).abs() < 1e-12,
                "max-frequency weighting should map GC% {gc_pct} into the second GC bin"
            );
        }
        Ok(())
    }

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn length_agnostic_requested_range_uses_only_overlapping_length_rows() -> Result<()> {
        // The package has two length rows:
        // - [20,30): correction [1,2]
        // - [30,40]: correction [3,5]
        //
        // Requested range [30,30] overlaps only the second row, so equal weighting should
        // collapse to [3,5] rather than the full-package mean [2,3.5].
        let package = make_length_agnostic_package();
        let corrector = GCCorrector::from_package(&package)?;
        let agnostic = LengthAgnosticGCCorrector::from_gc_corrector(
            &corrector,
            &MarginalizeLengthsWeightingScheme::Equal,
            GCLengthRange::Requested,
            0.0,
            30,
            30,
        )?;

        for gc_pct in 0..50 {
            assert!(
                (agnostic.get_correction_weight(gc_pct)? - 3.0).abs() < 1e-12,
                "requested range should map GC% {gc_pct} into the selected first GC bin"
            );
        }
        for gc_pct in 50..=100 {
            assert!(
                (agnostic.get_correction_weight(gc_pct)? - 5.0).abs() < 1e-12,
                "requested range should map GC% {gc_pct} into the selected second GC bin"
            );
        }
        Ok(())
    }

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn length_agnostic_package_range_keeps_full_package_rows() -> Result<()> {
        // Even with requested range [30,30], package range selection should keep both
        // package rows, preserving the full-package equal-weighted mean [2,3.5].
        let package = make_length_agnostic_package();
        let corrector = GCCorrector::from_package(&package)?;
        let agnostic = LengthAgnosticGCCorrector::from_gc_corrector(
            &corrector,
            &MarginalizeLengthsWeightingScheme::Equal,
            GCLengthRange::Package,
            0.0,
            30,
            30,
        )?;

        for gc_pct in 0..50 {
            assert!(
                (agnostic.get_correction_weight(gc_pct)? - 2.0).abs() < 1e-12,
                "package range should map GC% {gc_pct} into the full-package first GC bin"
            );
        }
        for gc_pct in 50..=100 {
            assert!(
                (agnostic.get_correction_weight(gc_pct)? - 3.5).abs() < 1e-12,
                "package range should map GC% {gc_pct} into the full-package second GC bin"
            );
        }
        Ok(())
    }

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn length_agnostic_requested_range_collapses_multiple_selected_rows() -> Result<()> {
        // Manual derivation:
        // - Length bins are [20,30), [30,40), and [40,50].
        // - Requested range [35,45] overlaps only the second and third rows.
        // - The first row has large corrections and 0.5 frequency, so including it by mistake
        //   would change every expected value.
        // - Selected-row frequencies are [0.125, 0.375], so frequency weighting must
        //   renormalize over the selected rows before averaging:
        //     (2 * 0.125 + 6 * 0.375) / 0.5 = 5
        //     (20 * 0.125 + 60 * 0.375) / 0.5 = 50
        //     (200 * 0.125 + 600 * 0.375) / 0.5 = 500
        let package = GCCorrectionPackage {
            version: GC_CORRECTION_SCHEMA_VERSION,
            end_offset: 0,
            length_edges: vec![20, 30, 40, 50],
            gc_edges: vec![0, 25, 50, 100],
            correction_matrix: array![
                [1000.0_f64, 10000.0_f64, 100000.0_f64],
                [2.0_f64, 20.0_f64, 200.0_f64],
                [6.0_f64, 60.0_f64, 600.0_f64],
            ],
            length_bin_frequencies: array![0.5_f64, 0.125_f64, 0.375_f64],
            reference_contig_footprint: Vec::new(),
        };
        let corrector = GCCorrector::from_package(&package)?;

        let expected_by_scheme = [
            (
                MarginalizeLengthsWeightingScheme::Equal,
                [4.0_f64, 40.0_f64, 400.0_f64],
            ),
            (
                MarginalizeLengthsWeightingScheme::Frequency,
                [5.0_f64, 50.0_f64, 500.0_f64],
            ),
            (
                MarginalizeLengthsWeightingScheme::MaxFrequency,
                [6.0_f64, 60.0_f64, 600.0_f64],
            ),
        ];

        for (scheme, expected) in expected_by_scheme {
            let agnostic = LengthAgnosticGCCorrector::from_gc_corrector(
                &corrector,
                &scheme,
                GCLengthRange::Requested,
                0.0,
                35,
                45,
            )?;

            let checked_gc_bins = [
                (0_usize, expected[0]),
                (24_usize, expected[0]),
                (25_usize, expected[1]),
                (49_usize, expected[1]),
                (50_usize, expected[2]),
                (100_usize, expected[2]),
            ];
            for (gc_pct, expected_weight) in checked_gc_bins {
                let observed = agnostic.get_correction_weight(gc_pct)?;
                assert!(
                    (observed - expected_weight).abs() < 1e-12,
                    "scheme {:?}, GC% {}: expected {}, got {}",
                    scheme,
                    gc_pct,
                    expected_weight,
                    observed
                );
            }
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
            fragment_lengths: FragmentLengthArgs {
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
            ref_gc_dir.join("ref_gc_package.zarr"),
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

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn gc_bias_default_min_mapq_matches_explicit_thirty_and_differs_from_explicit_zero()
    -> Result<()> {
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
                ref_gc_dir.path().join("ref_gc_package.zarr"),
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
            cfg.outlier_method = OutlierMethodArg::None;
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

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn default_windows_match_explicit_by_size_and_differ_from_global() -> Result<()> {
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

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn errors_when_blacklist_removes_all_usable_gc_support() -> Result<()> {
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
        let baseline_package = GCCorrectionPackage::from_file(
            baseline_out_dir.path().join("gc_bias_correction.zarr"),
        )?;
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

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn correction_package_propagates_reference_end_offset_for_single_length() -> Result<()> {
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
            GCCorrectionPackage::from_file(out_dir.path().join("gc_bias_correction.zarr"))?;

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

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn overlapping_and_touching_bed_windows_does_not_match_explicitly_merged_gc_bias_run()
    -> Result<()> {
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
            GCCorrectionPackage::from_file(split_out.path().join("gc_bias_correction.zarr"))?;
        let merged_pkg =
            GCCorrectionPackage::from_file(merged_out.path().join("gc_bias_correction.zarr"))?;
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

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn by_size_gc_bias_is_invariant_to_aligned_vs_misaligned_tile_sizes() -> Result<()> {
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
            GCCorrectionPackage::from_file(aligned_out.path().join("gc_bias_correction.zarr"))?;
        let misaligned_pkg =
            GCCorrectionPackage::from_file(misaligned_out.path().join("gc_bias_correction.zarr"))?;
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

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn multi_chromosome_by_size_gc_bias_accumulates_windows_across_chromosomes() -> Result<()> {
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
            fragment_lengths: FragmentLengthArgs {
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
            ref_gc_dir.path().join("ref_gc_package.zarr"),
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
            GCCorrectionPackage::from_file(out_dir.path().join("gc_bias_correction.zarr"))?;
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

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn empty_middle_tile_matches_single_tile_gc_bias_run() -> Result<()> {
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
            GCCorrectionPackage::from_file(multi_tile_out.path().join("gc_bias_correction.zarr"))?;
        let single_tile_pkg =
            GCCorrectionPackage::from_file(single_tile_out.path().join("gc_bias_correction.zarr"))?;
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

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn touching_bed_windows_match_by_size_counts_and_default_single_length_packages() -> Result<()>
    {
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
            GCCorrectionPackage::from_file(by_size_out.path().join("gc_bias_correction.zarr"))?;
        let by_bed_pkg =
            GCCorrectionPackage::from_file(by_bed_out.path().join("gc_bias_correction.zarr"))?;
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

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn real_ref_gc_bias_then_gc_bias_package_is_non_neutral_in_two_bin_case() -> Result<()> {
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
            windows: RefGCWindowsArgs {
                by_bed: Some(bed_path),
            },
            chromosomes: ChromosomeArgs {
                chromosomes: Some(vec!["chr1".to_string()]),
                chromosomes_file: None,
            },
            blacklist: None,
            fragment_lengths: FragmentLengthArgs {
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
        cfg.outlier_method = OutlierMethodArg::None;

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
            GCCorrectionPackage::from_file(out_dir.path().join("gc_bias_correction.zarr"))?;
        assert_eq!(package.length_edges, vec![10, 10]);
        assert_eq!(package.gc_edges, vec![0, 1, 100]);
        assert_eq!(package.correction_matrix.dim(), (1, 2));
        assert!((package.correction_matrix[(0, 0)] - 5.0).abs() < 1e-12);
        assert!((package.correction_matrix[(0, 1)] - (5.0 / 9.0)).abs() < 1e-12);
        assert_eq!(package.length_bin_frequencies.len(), 1);
        assert!((package.length_bin_frequencies[0] - 1.0).abs() < 1e-12);

        Ok(())
    }

    // KEEP-IN-TESTS: public API or command artifact behavior.
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
            windows: RefGCWindowsArgs {
                by_bed: Some(bed_path),
            },
            chromosomes: ChromosomeArgs {
                chromosomes: Some(vec!["chr1".to_string()]),
                chromosomes_file: None,
            },
            blacklist: None,
            fragment_lengths: FragmentLengthArgs {
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
        cfg.outlier_method = OutlierMethodArg::None;

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

    // KEEP-IN-TESTS: public API or command artifact behavior.
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
            windows: RefGCWindowsArgs {
                by_bed: Some(bed_path),
            },
            chromosomes: ChromosomeArgs {
                chromosomes: Some(vec!["chr1".to_string()]),
                chromosomes_file: None,
            },
            blacklist: None,
            fragment_lengths: FragmentLengthArgs {
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
        cfg.outlier_method = OutlierMethodArg::None;

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

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn save_intermediates_uses_output_prefixes_in_shared_output_directory() -> Result<()> {
        // Arrange:
        // Use one output directory for two runs with different prefixes. With this reference
        // package, each run writes exactly six intermediate arrays:
        //   0 avg_cfdna_counts
        //   1 normalized_avg_cfdna_counts
        //   2 binned_ref_counts
        //   3 binned_cfdna_counts
        //   4 normalized_binned_cfdna_counts
        //   5 normalized_binned_ref_counts
        //
        // The prefixes should make those twelve paths distinct in the shared directory.
        let bam = fixtures::simple_inward_bam()?;
        let reference = fixtures::simple_reference_twobit()?;
        let ref_gc_dir = TempDir::new()?;
        write_reference_package_for_single_length(&reference.path, &ref_gc_dir, 60, 0)?;

        let out_dir = TempDir::new()?;
        for prefix in ["sampleA", "sampleB"] {
            let mut cfg =
                make_gc_bias_cfg(&bam.bam, &reference.path, ref_gc_dir.path(), out_dir.path());
            cfg.set_windows(GCWindowsArgs {
                by_size: None,
                by_bed: None,
                global: true,
            });
            cfg.set_save_intermediates(true);
            cfg.set_output_prefix(prefix.to_string());

            run_gc_bias(&cfg)?;
        }

        // Assert:
        let mut intermediate_files: Vec<String> = std::fs::read_dir(out_dir.path())?
            .filter_map(|entry| entry.ok())
            .filter_map(|entry| {
                let name = entry.file_name().into_string().ok()?;
                if name.ends_with(".npy") {
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
                "sampleA.gc_bias.avg_cfdna_counts.0.npy".to_string(),
                "sampleA.gc_bias.binned_cfdna_counts.3.npy".to_string(),
                "sampleA.gc_bias.binned_ref_counts.2.npy".to_string(),
                "sampleA.gc_bias.normalized_avg_cfdna_counts.1.npy".to_string(),
                "sampleA.gc_bias.normalized_binned_cfdna_counts.4.npy".to_string(),
                "sampleA.gc_bias.normalized_binned_ref_counts.5.npy".to_string(),
                "sampleB.gc_bias.avg_cfdna_counts.0.npy".to_string(),
                "sampleB.gc_bias.binned_cfdna_counts.3.npy".to_string(),
                "sampleB.gc_bias.binned_ref_counts.2.npy".to_string(),
                "sampleB.gc_bias.normalized_avg_cfdna_counts.1.npy".to_string(),
                "sampleB.gc_bias.normalized_binned_cfdna_counts.4.npy".to_string(),
                "sampleB.gc_bias.normalized_binned_ref_counts.5.npy".to_string(),
            ]
        );

        let sample_a_avg: ndarray::Array2<f64> = read_npy(
            out_dir
                .path()
                .join("sampleA.gc_bias.avg_cfdna_counts.0.npy"),
        )?;
        let sample_b_avg: ndarray::Array2<f64> = read_npy(
            out_dir
                .path()
                .join("sampleB.gc_bias.avg_cfdna_counts.0.npy"),
        )?;
        assert_eq!(sample_a_avg, sample_b_avg);

        Ok(())
    }

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn min_window_acgt_pct_excludes_mostly_blacklisted_window_but_keeps_clean_window() -> Result<()>
    {
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

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn multiple_blacklist_files_with_touching_intervals_match_single_merged_gc_bias_run()
    -> Result<()> {
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
            GCCorrectionPackage::from_file(split_out_dir.path().join("gc_bias_correction.zarr"))?;
        let merged_package =
            GCCorrectionPackage::from_file(merged_out_dir.path().join("gc_bias_correction.zarr"))?;
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

    // KEEP-IN-TESTS: public API or command artifact behavior.
    #[test]
    fn gc_bias_run_rejects_reference_gc_package_from_different_reference_footprint() -> Result<()> {
        // Human verification status: verified
        // Manual expectations:
        // - The reference GC package is built from `reference_a`, whose 2bit footprint is stored
        //   in the package.
        // - The run uses `reference_b`, which has the same selected chromosome name but a different
        //   contig footprint.
        // - `gc-bias` should reject the mismatch before creating the output directory or writing a
        //   downstream correction package.
        let reference_a = fixtures::simple_reference_twobit()?;
        let reference_b = fixtures::twobit_from_sequences(
            "gc_bias_reference_footprint_mismatch",
            vec![("chr1".to_string(), "ACGT".repeat(80))],
        )?;
        let bam = fixtures::simple_inward_bam()?;
        let ref_gc_dir = TempDir::new()?;
        write_reference_package_for_single_length(&reference_a.path, &ref_gc_dir, 60, 0)?;
        let output_parent = TempDir::new()?;
        let output_dir = output_parent.path().join("not_created");
        let cfg = make_gc_bias_cfg(&bam.bam, &reference_b.path, ref_gc_dir.path(), &output_dir);

        let err = run_gc_bias(&cfg).expect_err("reference-footprint mismatch should fail");
        let msg = err.to_string();

        assert!(
            msg.contains(
                "Reference GC package was built against a different reference contig footprint"
            ),
            "unexpected error message: {msg}"
        );
        assert!(
            !output_dir.exists(),
            "reference-footprint validation should fail before creating output_dir"
        );

        Ok(())
    }
}
