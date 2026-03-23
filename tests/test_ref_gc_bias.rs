use anyhow::Result;
use ndarray::array;
use ndarray_npy::NpzReader;
use std::path::Path;
use tempfile::TempDir;

use cfdnalab::{
    commands::cli_common::ChromosomeArgs,
    commands::gc_bias::counting::{
        GCCounts, build_gc_prefixes, count_reference_gc_and_length_by_window,
    },
    commands::gc_bias::support_masking::create_support_mask_threshold_per_mb,
    commands::ref_gc_bias::{config::RefGCBiasConfig, ref_gc_bias::run},
    shared::{
        blacklist::apply_blacklist_mask_to_seq,
        interval::{IndexedInterval, Interval},
    },
};
mod fixtures;

fn load_ref_gc_package_arrays(
    package_path: &Path,
) -> Result<(
    ndarray::Array2<f64>,
    ndarray::Array2<bool>,
    ndarray::Array2<bool>,
    ndarray::Array2<u16>,
)> {
    let file = std::fs::File::open(package_path)?;
    let mut npz = NpzReader::new(file)?;
    let counts: ndarray::Array2<f64> = npz.by_name("counts")?;
    let support_unobservables: ndarray::Array2<bool> =
        npz.by_name("support_mask_unobservables")?;
    let support_outliers: ndarray::Array2<bool> = npz.by_name("support_mask_outliers")?;
    let gc_percent_widths: ndarray::Array2<u16> = npz.by_name("gc_percent_widths")?;
    Ok((
        counts,
        support_unobservables,
        support_outliers,
        gc_percent_widths,
    ))
}

#[test]
fn counts_gc_for_each_window_with_end_offset() -> Result<()> {
    // Arrange: Two windows of equal size with start positions seeded at each window start
    // End offset trims one base on each side so GC is counted on the inner span.
    let seq = b"ACGTACGTACGT".to_vec();
    let prefixes = build_gc_prefixes(&seq);
    let windows = IndexedInterval::from_tuples(&[(0, 6, 0), (6, 12, 1)])?;
    let starts = vec![0usize, 6usize];
    let mut counts_by_bin = vec![
        GCCounts::new(4, 6, 1, (0, 0))?,
        GCCounts::new(4, 6, 1, (0, 0))?,
    ];

    // Act
    count_reference_gc_and_length_by_window(
        &mut counts_by_bin,
        &prefixes,
        (4, 7),
        windows.as_slice(),
        starts.as_slice(),
        seq.len() as u64,
        1.0,
        1,
        1,
    );

    // Assert: Window 0 sees GC=2 for lengths 4, 5, and 6; window 1 sees an increasing GC profile
    let window0 = &counts_by_bin[0];
    assert_eq!(window0.get(4, 2).unwrap(), 1.0);
    assert_eq!(window0.get(5, 2).unwrap(), 1.0);
    assert_eq!(window0.get(6, 2).unwrap(), 1.0);

    let window1 = &counts_by_bin[1];
    assert_eq!(window1.get(4, 0).unwrap(), 1.0);
    assert_eq!(window1.get(5, 1).unwrap(), 1.0);
    assert_eq!(window1.get(6, 2).unwrap(), 1.0);
    Ok(())
}

#[test]
fn skips_counts_after_blacklist_removes_acgt_support() -> Result<()> {
    // Arrange: Blacklist the middle of the fragment so only half the bases remain ACGT
    let mut seq = b"ACGT".to_vec();
    let blacklist_intervals = Interval::from_tuples(&[(1, 3)])?;
    apply_blacklist_mask_to_seq(&mut seq, &blacklist_intervals, 0);
    let prefixes = build_gc_prefixes(&seq);
    let windows = IndexedInterval::from_tuples(&[(0, 4, 0)])?;
    let starts = vec![0usize];
    let mut counts_by_bin = vec![GCCounts::new(4, 4, 0, (0, 0))?];

    // Act
    count_reference_gc_and_length_by_window(
        &mut counts_by_bin,
        &prefixes,
        (4, 5),
        windows.as_slice(),
        starts.as_slice(),
        seq.len() as u64,
        1.0,
        1,
        0,
    );

    // Assert: Masking drops the ACGT fraction below 1.0 so no counts are emitted
    assert_eq!(counts_by_bin[0].sum(), 0.0);
    Ok(())
}

#[test]
fn support_threshold_per_mb_steps_at_hundred_million_positions() -> Result<()> {
    // Arrange:
    // Use exactly 1 Mb of valid ACGT positions so the absolute threshold equals `threshold_per_mb`.
    // The row values are chosen so each threshold step changes exactly one additional support bit.
    let counts = array![[0.5_f64, 1.5_f64, 2.5_f64, 3.5_f64]];
    let num_acgt_positions = 1_000_000_u64;

    let scenarios = [
        (
            99_999_999_usize,
            1_usize,
            vec![false, true, true, true],
        ),
        (
            100_000_000_usize,
            2_usize,
            vec![false, false, true, true],
        ),
        (
            200_000_000_usize,
            3_usize,
            vec![false, false, false, true],
        ),
    ];

    for (n_positions, expected_threshold_per_mb, expected_mask_row) in scenarios {
        // Act:
        // The command computes:
        //   threshold_per_mb = 1 + n_positions / 100_000_000
        // with integer division.
        let threshold_per_mb = 1 + n_positions / 100_000_000;
        let mask = create_support_mask_threshold_per_mb(
            &[counts.clone()],
            num_acgt_positions,
            threshold_per_mb as f64,
        )
        .expect("support mask should be created");

        // Assert:
        // Because num_acgt_positions is exactly 1 Mb, the usable-count threshold is:
        //   threshold = 1.0 * threshold_per_mb
        // So the support mask should flip at the exact crossover boundaries.
        assert_eq!(threshold_per_mb, expected_threshold_per_mb);
        assert_eq!(mask.row(0).to_vec(), expected_mask_row);
    }

    Ok(())
}

#[test]
fn ref_gc_bias_run_writes_expected_package_metadata_and_shapes() -> Result<()> {
    let reference = fixtures::simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let cfg = RefGCBiasConfig {
        ref_genome: cfdnalab::commands::cli_common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: out_dir.path().to_path_buf(),
        n_threads: 1,
        n_positions: 100,
        seed: Some(7),
        windows: Default::default(),
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: None,
        fragment_lengths: cfdnalab::commands::cli_common::FragmentLengthArgs {
            min_fragment_length: 10,
            max_fragment_length: 12,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 1_000_000,
    };
    cfg.check_smoothing_settings()?;

    // Manual expectations:
    // - We request fragment lengths 10, 11, and 12, so the output must have 3 rows.
    // - The package stores GC as integer percentages 0..=100, so the counts and both support masks
    //   must have 101 columns.
    // - We explicitly disable smoothing and interpolation, so the written scalar metadata must
    //   reflect those choices exactly.
    // - This is a command-level test of the written artifact contract, not of the exact sampled counts.
    run(&cfg)?;

    let package_path = out_dir.path().join("ref_gc_package.npz");
    let file = std::fs::File::open(&package_path)?;
    let mut npz = NpzReader::new(file)?;

    let counts: ndarray::Array2<f64> = npz.by_name("counts")?;
    let support_unobservables: ndarray::Array2<bool> =
        npz.by_name("support_mask_unobservables")?;
    let support_outliers: ndarray::Array2<bool> = npz.by_name("support_mask_outliers")?;
    let gc_percent_widths: ndarray::Array2<u16> = npz.by_name("gc_percent_widths")?;
    let version: ndarray::Array1<u32> = npz.by_name("version")?;
    let length_range: ndarray::Array1<u32> = npz.by_name("length_range")?;
    let end_offset: ndarray::Array1<u32> = npz.by_name("end_offset")?;
    let skip_interpolation: ndarray::Array1<bool> = npz.by_name("skip_interpolation")?;
    let smoothing_radius: ndarray::Array1<u32> = npz.by_name("smoothing_radius")?;
    let smoothing_sigma: ndarray::Array1<f64> = npz.by_name("smoothing_sigma")?;
    let skip_smoothing: ndarray::Array1<bool> = npz.by_name("skip_smoothing")?;

    assert_eq!(counts.dim(), (3, 101));
    assert_eq!(support_unobservables.dim(), (3, 101));
    assert_eq!(support_outliers.dim(), (3, 101));
    assert_eq!(gc_percent_widths.dim(), (3, 101));

    assert_eq!(version.to_vec(), vec![cfdnalab::commands::gc_bias::GC_CORRECTION_SCHEMA_VERSION]);
    assert_eq!(length_range.to_vec(), vec![10, 12]);
    assert_eq!(end_offset.to_vec(), vec![0]);
    assert_eq!(skip_interpolation.to_vec(), vec![true]);
    assert_eq!(smoothing_radius.to_vec(), vec![2]);
    assert_eq!(smoothing_sigma.to_vec(), vec![0.55]);
    assert_eq!(skip_smoothing.to_vec(), vec![true]);

    assert!(
        counts.iter().all(|value| value.is_finite() && *value >= 0.0),
        "reference GC counts should be finite and non-negative"
    );

    Ok(())
}

#[test]
fn overlapping_and_touching_bed_windows_match_explicitly_merged_ref_gc_bias_run() -> Result<()> {
    let reference = fixtures::simple_reference_twobit()?;
    let split_out_dir = TempDir::new()?;
    let merged_out_dir = TempDir::new()?;
    let bed_split = split_out_dir.path().join("split_windows.bed");
    let bed_merged = merged_out_dir.path().join("merged_windows.bed");

    std::fs::write(&bed_split, "chr1\t0\t40\nchr1\t20\t60\nchr1\t60\t80\n")?;
    std::fs::write(&bed_merged, "chr1\t0\t80\n")?;

    let make_cfg = |output_dir: &Path, bed_path: &Path| RefGCBiasConfig {
        ref_genome: cfdnalab::commands::cli_common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: output_dir.to_path_buf(),
        n_threads: 1,
        n_positions: 100,
        seed: Some(11),
        windows: cfdnalab::commands::ref_gc_bias::config::RefGCWindowsArgs {
            by_bed: Some(bed_path.to_path_buf()),
        },
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: None,
        fragment_lengths: cfdnalab::commands::cli_common::FragmentLengthArgs {
            min_fragment_length: 10,
            max_fragment_length: 12,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 1_000_000,
    };

    // Manual expectations:
    // - `ref-gc-bias` flattens overlapping and touching BED windows to unique positions before
    //   counting.
    // - The split BED:
    //     [0,40), [20,60), [60,80)
    //   has the same unique covered positions as the explicitly merged BED:
    //     [0,80)
    // - Sampling is independent of windowing in this command, and the seed is fixed.
    // - Therefore both runs must produce exactly the same reference package arrays:
    //   counts, both support masks, and GC-percent widths.
    run(&make_cfg(split_out_dir.path(), &bed_split))?;
    run(&make_cfg(merged_out_dir.path(), &bed_merged))?;

    let split_package = split_out_dir.path().join("ref_gc_package.npz");
    let merged_package = merged_out_dir.path().join("ref_gc_package.npz");
    let (
        split_counts,
        split_support_unobservables,
        split_support_outliers,
        split_gc_percent_widths,
    ) = load_ref_gc_package_arrays(&split_package)?;
    let (
        merged_counts,
        merged_support_unobservables,
        merged_support_outliers,
        merged_gc_percent_widths,
    ) = load_ref_gc_package_arrays(&merged_package)?;

    assert_eq!(split_counts, merged_counts);
    assert_eq!(split_support_unobservables, merged_support_unobservables);
    assert_eq!(split_support_outliers, merged_support_outliers);
    assert_eq!(split_gc_percent_widths, merged_gc_percent_widths);

    Ok(())
}

#[test]
fn full_chromosome_bed_window_matches_global_ref_gc_bias_run() -> Result<()> {
    // Arrange:
    // `ref-gc-bias` samples candidate fragment starts independently of windowing, then counts only
    // starts that land inside the configured windows. A single BED window spanning the whole
    // chromosome therefore describes exactly the same logical region as global mode.
    //
    // With a fixed seed, identical length range, and identical blacklist settings, the full
    // written reference package must match exactly: counts, both support masks, and GC-percent
    // widths.
    let reference = fixtures::simple_reference_twobit()?;
    let global_out_dir = TempDir::new()?;
    let bed_out_dir = TempDir::new()?;
    let bed_path = bed_out_dir.path().join("whole_chr.bed");
    std::fs::write(&bed_path, "chr1\t0\t256\n")?;

    let make_global_cfg = |output_dir: &Path| RefGCBiasConfig {
        ref_genome: cfdnalab::commands::cli_common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: output_dir.to_path_buf(),
        n_threads: 1,
        n_positions: 100,
        seed: Some(13),
        windows: Default::default(),
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: None,
        fragment_lengths: cfdnalab::commands::cli_common::FragmentLengthArgs {
            min_fragment_length: 10,
            max_fragment_length: 12,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 80,
    };
    let make_bed_cfg = |output_dir: &Path| RefGCBiasConfig {
        windows: cfdnalab::commands::ref_gc_bias::config::RefGCWindowsArgs {
            by_bed: Some(bed_path.clone()),
        },
        ..make_global_cfg(output_dir)
    };

    // Act
    run(&make_global_cfg(global_out_dir.path()))?;
    run(&make_bed_cfg(bed_out_dir.path()))?;

    // Assert
    let global_package = global_out_dir.path().join("ref_gc_package.npz");
    let bed_package = bed_out_dir.path().join("ref_gc_package.npz");
    let (
        global_counts,
        global_support_unobservables,
        global_support_outliers,
        global_gc_percent_widths,
    ) = load_ref_gc_package_arrays(&global_package)?;
    let (
        bed_counts,
        bed_support_unobservables,
        bed_support_outliers,
        bed_gc_percent_widths,
    ) = load_ref_gc_package_arrays(&bed_package)?;

    assert_eq!(global_counts, bed_counts);
    assert_eq!(global_support_unobservables, bed_support_unobservables);
    assert_eq!(global_support_outliers, bed_support_outliers);
    assert_eq!(global_gc_percent_widths, bed_gc_percent_widths);

    Ok(())
}

#[test]
fn fixed_seed_ref_gc_bias_is_invariant_to_tile_size() -> Result<()> {
    // Arrange:
    // With a fixed seed, tile size is only an execution detail. Changing it changes how sampled
    // starts are partitioned and merged across workers, but not the final logical reference-GC
    // package.
    //
    // We compare two multi-tile runs on the 256 bp fixture chromosome:
    // - tile_size = 128 -> 2 tiles
    // - tile_size =  80 -> 4 tiles
    //
    // Both runs keep everything else identical, including the seed, so the written package arrays
    // must match exactly.
    let reference = fixtures::simple_reference_twobit()?;
    let large_tile_out = TempDir::new()?;
    let small_tile_out = TempDir::new()?;

    let make_cfg = |output_dir: &Path, tile_size: u32| RefGCBiasConfig {
        ref_genome: cfdnalab::commands::cli_common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: output_dir.to_path_buf(),
        n_threads: 1,
        n_positions: 100,
        seed: Some(17),
        windows: Default::default(),
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: None,
        fragment_lengths: cfdnalab::commands::cli_common::FragmentLengthArgs {
            min_fragment_length: 10,
            max_fragment_length: 12,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size,
    };

    // Act
    run(&make_cfg(large_tile_out.path(), 128))?;
    run(&make_cfg(small_tile_out.path(), 80))?;

    // Assert
    let large_package = large_tile_out.path().join("ref_gc_package.npz");
    let small_package = small_tile_out.path().join("ref_gc_package.npz");
    let (
        large_counts,
        large_support_unobservables,
        large_support_outliers,
        large_gc_percent_widths,
    ) = load_ref_gc_package_arrays(&large_package)?;
    let (
        small_counts,
        small_support_unobservables,
        small_support_outliers,
        small_gc_percent_widths,
    ) = load_ref_gc_package_arrays(&small_package)?;

    assert_eq!(large_counts, small_counts);
    assert_eq!(large_support_unobservables, small_support_unobservables);
    assert_eq!(large_support_outliers, small_support_outliers);
    assert_eq!(large_gc_percent_widths, small_gc_percent_widths);

    Ok(())
}

#[test]
fn rejects_n_positions_when_sampling_density_would_exceed_one() -> Result<()> {
    let reference = fixtures::simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let cfg = RefGCBiasConfig {
        ref_genome: cfdnalab::commands::cli_common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: out_dir.path().to_path_buf(),
        n_threads: 1,
        // `simple_reference_twobit()` uses chr1 length 256.
        // With max fragment length 60, the number of valid starts is:
        //   256 - 60 + 1 = 197.
        // Asking for 198 positions therefore gives:
        //   sampling_density = 198 / 197 > 1.0
        // which the command must reject before counting.
        n_positions: 198,
        seed: Some(7),
        windows: Default::default(),
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: None,
        fragment_lengths: cfdnalab::commands::cli_common::FragmentLengthArgs {
            min_fragment_length: 60,
            max_fragment_length: 60,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 1_000_000,
    };

    let err = run(&cfg).expect_err("sampling density above 1.0 should fail");
    let err_text = err.to_string();
    assert!(
        err_text.contains("Sampling density") && err_text.contains("exceeds 1.0"),
        "unexpected error message: {err_text}"
    );

    Ok(())
}

#[test]
fn fixed_seed_ref_gc_bias_is_invariant_to_thread_count() -> Result<()> {
    let reference = fixtures::simple_reference_twobit()?;
    let single_thread_out = TempDir::new()?;
    let two_thread_out = TempDir::new()?;

    let make_cfg = |output_dir: &Path, n_threads: usize| RefGCBiasConfig {
        ref_genome: cfdnalab::commands::cli_common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: output_dir.to_path_buf(),
        n_threads,
        n_positions: 100,
        seed: Some(7),
        windows: Default::default(),
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: None,
        fragment_lengths: cfdnalab::commands::cli_common::FragmentLengthArgs {
            min_fragment_length: 10,
            max_fragment_length: 12,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 80,
    };

    // Manual expectations:
    // - The seed is fixed, and `ref-gc-bias` derives deterministic per-tile seeds from it before
    //   parallel work starts.
    // - Changing only `n_threads` must therefore change runtime only, not the sampled starts or
    //   the merged reference package.
    // - We use `tile_size = 80` on the 256 bp fixture chromosome so the run spans multiple tiles,
    //   making the parallel merge path real rather than vacuous.
    run(&make_cfg(single_thread_out.path(), 1))?;
    run(&make_cfg(two_thread_out.path(), 2))?;

    let single_package = single_thread_out.path().join("ref_gc_package.npz");
    let two_package = two_thread_out.path().join("ref_gc_package.npz");
    let (
        single_counts,
        single_support_unobservables,
        single_support_outliers,
        single_gc_percent_widths,
    ) = load_ref_gc_package_arrays(&single_package)?;
    let (
        two_counts,
        two_support_unobservables,
        two_support_outliers,
        two_gc_percent_widths,
    ) = load_ref_gc_package_arrays(&two_package)?;

    assert_eq!(single_counts, two_counts);
    assert_eq!(single_support_unobservables, two_support_unobservables);
    assert_eq!(single_support_outliers, two_support_outliers);
    assert_eq!(single_gc_percent_widths, two_gc_percent_widths);

    Ok(())
}
