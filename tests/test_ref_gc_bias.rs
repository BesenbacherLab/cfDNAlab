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
fn ref_gc_bias_run_counts_expected_two_bin_reference_distribution() -> Result<()> {
    // Arrange:
    // Build a 200 bp chromosome with two pure regions:
    // - chr1[0,100)   = all A
    // - chr1[100,200) = all C
    //
    // We fix fragment length at 10 and sample every valid start:
    //   valid starts = 200 - 10 + 1 = 191
    // so `n_positions = 191` means every possible start is counted exactly once.
    //
    // Then we restrict counting to BED windows:
    //   [0, 91) and [100, 191)
    // which keeps only the pure-A starts 0..=90 and pure-C starts 100..=190.
    // The boundary-crossing starts 91..=99 are still sampled, but excluded by the windows.
    //
    // Therefore the output counts are exactly:
    // - 91 counts at GC% 0
    // - 91 counts at GC% 100
    // - 0 everywhere else
    //
    // For length 10 and end_offset 0, theoretical GC% support exists only at multiples of 10:
    //   {0, 10, 20, ..., 100}
    // so the theoretical support mask must be true exactly at those 11 bins.
    //
    // The outlier-support threshold is:
    //   threshold_per_mb = 1 + 191 / 100_000_000 = 1
    // and the total covered ACGT bases are:
    //   91 + 91 = 182
    // so the absolute threshold is only 182 / 1_000_000 = 0.000182.
    // That means only the two non-zero bins pass the empirical support threshold.
    let reference = fixtures::twobit_from_sequences(
        "ref_gc_bias_two_bin_reference",
        vec![(
            "chr1".to_string(),
            format!("{}{}", "A".repeat(100), "C".repeat(100)),
        )],
    )?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("pure_windows.bed");
    std::fs::write(&bed_path, "chr1\t0\t91\nchr1\t100\t191\n")?;

    let cfg = RefGCBiasConfig {
        ref_genome: cfdnalab::commands::cli_common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: out_dir.path().to_path_buf(),
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
    };

    // Act
    run(&cfg)?;

    // Assert
    let package_path = out_dir.path().join("ref_gc_package.npz");
    let (counts, support_unobservables, support_outliers, gc_percent_widths) =
        load_ref_gc_package_arrays(&package_path)?;

    assert_eq!(counts.dim(), (1, 101));
    assert_eq!(support_unobservables.dim(), (1, 101));
    assert_eq!(support_outliers.dim(), (1, 101));
    assert_eq!(gc_percent_widths.dim(), (1, 101));

    for (gc_pct, &value) in counts.row(0).iter().enumerate() {
        match gc_pct {
            0 | 100 => assert!(
                (value - 91.0).abs() < 1e-12,
                "expected 91 counts at GC% {gc_pct}, got {value}"
            ),
            _ => assert!(
                value.abs() < 1e-12,
                "expected 0 counts outside GC% 0/100, got bin {gc_pct}={value}"
            ),
        }
    }

    for gc_pct in 0..=100 {
        let theoretical = gc_pct % 10 == 0;
        assert_eq!(
            support_unobservables[(0, gc_pct)],
            theoretical,
            "unexpected theoretical support at GC% {gc_pct}"
        );
    }

    for gc_pct in 0..=100 {
        let empirical = matches!(gc_pct, 0 | 100);
        assert_eq!(
            support_outliers[(0, gc_pct)],
            empirical,
            "unexpected empirical support at GC% {gc_pct}"
        );
    }

    // With length 10, each reachable GC% corresponds to exactly one GC-count state, so the width
    // correction array should be 1 on the reachable multiples of 10 and 0 elsewhere.
    for gc_pct in 0..=100 {
        let expected_width = if gc_pct % 10 == 0 { 1 } else { 0 };
        assert_eq!(
            gc_percent_widths[(0, gc_pct)],
            expected_width,
            "unexpected GC-percent width at GC% {gc_pct}"
        );
    }

    Ok(())
}

#[test]
fn ref_gc_bias_run_blacklist_removes_exactly_the_overlapping_start_positions() -> Result<()> {
    // Arrange:
    // Start from the same exact-count setup as the previous test:
    // - chr1[0,100)   = all A
    // - chr1[100,200) = all C
    // - fragment length fixed at 10
    // - every valid start sampled exactly once (`n_positions = 191`)
    // - BED windows keep only the pure-A starts 0..=90 and pure-C starts 100..=190
    //
    // Then apply a blacklist over [45,55). `ref-gc-bias` masks blacklisted sequence to `N` before
    // counting, and with the command's default `min_acgt_fraction = 1.0`, any fragment touching
    // that mask is dropped entirely.
    //
    // A 10 bp fragment starting at `s` spans [s, s+10), so it overlaps [45,55) exactly when:
    //   s < 55 and s + 10 > 45
    // -> s <= 54 and s >= 36
    // -> s in 36..=54
    //
    // Those 19 dropped starts all belong to the pure-A BED window [0,91). The pure-C starts
    // 100..=190 never touch the blacklist.
    //
    // Therefore the exact output counts must be:
    // - GC% 0   -> 91 - 19 = 72
    // - GC% 100 -> 91
    // - 0 everywhere else
    //
    // The empirical support threshold is still tiny:
    //   total covered ACGT bases = 72 + 91 = 163
    //   threshold = 163 / 1_000_000 = 0.000163
    // so both non-zero bins remain empirically supported and all zero bins remain unsupported.
    let reference = fixtures::twobit_from_sequences(
        "ref_gc_bias_blacklist_exact_counts_reference",
        vec![(
            "chr1".to_string(),
            format!("{}{}", "A".repeat(100), "C".repeat(100)),
        )],
    )?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("pure_windows.bed");
    let blacklist_path = out_dir.path().join("blacklist.bed");
    std::fs::write(&bed_path, "chr1\t0\t91\nchr1\t100\t191\n")?;
    std::fs::write(&blacklist_path, "chr1\t45\t55\n")?;

    let cfg = RefGCBiasConfig {
        ref_genome: cfdnalab::commands::cli_common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: out_dir.path().to_path_buf(),
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
        blacklist: Some(vec![blacklist_path]),
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
    };

    // Act
    run(&cfg)?;

    // Assert
    let package_path = out_dir.path().join("ref_gc_package.npz");
    let (counts, support_unobservables, support_outliers, gc_percent_widths) =
        load_ref_gc_package_arrays(&package_path)?;

    assert_eq!(counts.dim(), (1, 101));
    assert_eq!(support_unobservables.dim(), (1, 101));
    assert_eq!(support_outliers.dim(), (1, 101));
    assert_eq!(gc_percent_widths.dim(), (1, 101));

    for (gc_pct, &value) in counts.row(0).iter().enumerate() {
        match gc_pct {
            0 => assert!(
                (value - 72.0).abs() < 1e-12,
                "expected 72 counts at GC% 0 after blacklist masking, got {value}"
            ),
            100 => assert!(
                (value - 91.0).abs() < 1e-12,
                "expected 91 counts at GC% 100 after blacklist masking, got {value}"
            ),
            _ => assert!(
                value.abs() < 1e-12,
                "expected 0 counts outside GC% 0/100, got bin {gc_pct}={value}"
            ),
        }
    }

    for gc_pct in 0..=100 {
        let theoretical = gc_pct % 10 == 0;
        assert_eq!(
            support_unobservables[(0, gc_pct)],
            theoretical,
            "unexpected theoretical support at GC% {gc_pct}"
        );
    }

    for gc_pct in 0..=100 {
        let empirical = matches!(gc_pct, 0 | 100);
        assert_eq!(
            support_outliers[(0, gc_pct)],
            empirical,
            "unexpected empirical support at GC% {gc_pct}"
        );
    }

    for gc_pct in 0..=100 {
        let expected_width = if gc_pct % 10 == 0 { 1 } else { 0 };
        assert_eq!(
            gc_percent_widths[(0, gc_pct)],
            expected_width,
            "unexpected GC-percent width at GC% {gc_pct}"
        );
    }

    Ok(())
}

#[test]
fn ref_gc_bias_run_end_offset_counts_expected_trimmed_two_bin_distribution() -> Result<()> {
    // Arrange:
    // Build the same 200 bp pure-reference setup:
    // - chr1[0,100)   = all A
    // - chr1[100,200) = all C
    //
    // But now use:
    // - fragment length = 12
    // - end_offset = 1
    //
    // GC is therefore counted on the contracted inner span [start+1, start+11), which is 10 bp
    // long. So although the fragment length is 12, the GC-percentage geometry is the same as an
    // effective length of 10:
    // - reachable GC% bins are exactly the multiples of 10
    // - each reachable GC% still has width 1
    //
    // We sample every valid start:
    //   valid starts = 200 - 12 + 1 = 189
    //
    // To keep only pure trimmed spans:
    // - pure-A requires start + 11 <= 100  -> start <= 89
    // - pure-C requires start + 1  >= 100  -> start >= 99
    //
    // So the BED windows:
    //   [0, 90)   -> starts 0..=89   -> 90 pure-A fragments
    //   [99, 189) -> starts 99..=188 -> 90 pure-C fragments
    //
    // Starts 90..=98 are still sampled, but excluded by the BED windows because their trimmed
    // spans would cross the A/C boundary.
    //
    // Therefore the exact output counts must be:
    // - 90 counts at GC% 0
    // - 90 counts at GC% 100
    // - 0 everywhere else
    //
    // The empirical support threshold is again tiny:
    //   threshold_per_mb = 1 + 189 / 100_000_000 = 1
    //   total covered ACGT bases = 90 + 90 = 180
    //   threshold = 180 / 1_000_000 = 0.00018
    // so the two non-zero bins remain empirically supported.
    let reference = fixtures::twobit_from_sequences(
        "ref_gc_bias_end_offset_two_bin_reference",
        vec![(
            "chr1".to_string(),
            format!("{}{}", "A".repeat(100), "C".repeat(100)),
        )],
    )?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("pure_trimmed_windows.bed");
    std::fs::write(&bed_path, "chr1\t0\t90\nchr1\t99\t189\n")?;

    let cfg = RefGCBiasConfig {
        ref_genome: cfdnalab::commands::cli_common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: out_dir.path().to_path_buf(),
        n_threads: 1,
        n_positions: 189,
        seed: Some(37),
        windows: cfdnalab::commands::ref_gc_bias::config::RefGCWindowsArgs {
            by_bed: Some(bed_path),
        },
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: None,
        fragment_lengths: cfdnalab::commands::cli_common::FragmentLengthArgs {
            min_fragment_length: 12,
            max_fragment_length: 12,
        },
        end_offset: 1,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 1_000_000,
    };

    // Act
    run(&cfg)?;

    // Assert
    let package_path = out_dir.path().join("ref_gc_package.npz");
    let (counts, support_unobservables, support_outliers, gc_percent_widths) =
        load_ref_gc_package_arrays(&package_path)?;

    assert_eq!(counts.dim(), (1, 101));
    assert_eq!(support_unobservables.dim(), (1, 101));
    assert_eq!(support_outliers.dim(), (1, 101));
    assert_eq!(gc_percent_widths.dim(), (1, 101));

    for (gc_pct, &value) in counts.row(0).iter().enumerate() {
        match gc_pct {
            0 | 100 => assert!(
                (value - 90.0).abs() < 1e-12,
                "expected 90 counts at GC% {gc_pct}, got {value}"
            ),
            _ => assert!(
                value.abs() < 1e-12,
                "expected 0 counts outside GC% 0/100, got bin {gc_pct}={value}"
            ),
        }
    }

    // Theoretical support depends on effective GC-counting length 12 - 2*1 = 10, not on the raw
    // fragment length. So reachable GC% bins are still exactly the multiples of 10.
    for gc_pct in 0..=100 {
        let theoretical = gc_pct % 10 == 0;
        assert_eq!(
            support_unobservables[(0, gc_pct)],
            theoretical,
            "unexpected theoretical support at GC% {gc_pct}"
        );
    }

    for gc_pct in 0..=100 {
        let empirical = matches!(gc_pct, 0 | 100);
        assert_eq!(
            support_outliers[(0, gc_pct)],
            empirical,
            "unexpected empirical support at GC% {gc_pct}"
        );
    }

    for gc_pct in 0..=100 {
        let expected_width = if gc_pct % 10 == 0 { 1 } else { 0 };
        assert_eq!(
            gc_percent_widths[(0, gc_pct)],
            expected_width,
            "unexpected GC-percent width at GC% {gc_pct}"
        );
    }

    Ok(())
}

#[test]
fn ref_gc_bias_run_blacklist_with_end_offset_drops_only_trimmed_overlaps() -> Result<()> {
    // Arrange:
    // Reuse the trimmed-span setup from the previous test:
    // - chr1[0,100)   = all A
    // - chr1[100,200) = all C
    // - fragment length = 12
    // - end_offset = 1
    // - every valid start sampled exactly once (`n_positions = 189`)
    // - BED windows keep only pure trimmed spans:
    //     [0,90)   -> starts 0..=89
    //     [99,189) -> starts 99..=188
    //
    // Now apply blacklist [45,55). In `ref-gc-bias`, blacklist masking is evaluated on the
    // effective GC-counting interval [start+1, start+11), not on the raw 12 bp fragment span.
    //
    // A trimmed 10 bp span overlaps [45,55) exactly when:
    //   start + 1  < 55
    //   start + 11 > 45
    // -> start <= 53 and start >= 35
    // -> start in 35..=53
    //
    // So exactly 19 pure-A starts are dropped. The pure-C starts 99..=188 remain untouched.
    //
    // This is intentionally different from raw-span overlap:
    // - raw 12 bp fragment overlap would drop starts 34..=54 (21 starts)
    // - the command should *not* do that here
    //
    // Therefore the exact output counts must be:
    // - GC% 0   -> 90 - 19 = 71
    // - GC% 100 -> 90
    // - 0 everywhere else
    let reference = fixtures::twobit_from_sequences(
        "ref_gc_bias_end_offset_blacklist_reference",
        vec![(
            "chr1".to_string(),
            format!("{}{}", "A".repeat(100), "C".repeat(100)),
        )],
    )?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("pure_trimmed_windows.bed");
    let blacklist_path = out_dir.path().join("blacklist.bed");
    std::fs::write(&bed_path, "chr1\t0\t90\nchr1\t99\t189\n")?;
    std::fs::write(&blacklist_path, "chr1\t45\t55\n")?;

    let cfg = RefGCBiasConfig {
        ref_genome: cfdnalab::commands::cli_common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: out_dir.path().to_path_buf(),
        n_threads: 1,
        n_positions: 189,
        seed: Some(37),
        windows: cfdnalab::commands::ref_gc_bias::config::RefGCWindowsArgs {
            by_bed: Some(bed_path),
        },
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: Some(vec![blacklist_path]),
        fragment_lengths: cfdnalab::commands::cli_common::FragmentLengthArgs {
            min_fragment_length: 12,
            max_fragment_length: 12,
        },
        end_offset: 1,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 1_000_000,
    };

    // Act
    run(&cfg)?;

    // Assert
    let package_path = out_dir.path().join("ref_gc_package.npz");
    let (counts, support_unobservables, support_outliers, gc_percent_widths) =
        load_ref_gc_package_arrays(&package_path)?;

    assert_eq!(counts.dim(), (1, 101));
    assert_eq!(support_unobservables.dim(), (1, 101));
    assert_eq!(support_outliers.dim(), (1, 101));
    assert_eq!(gc_percent_widths.dim(), (1, 101));

    for (gc_pct, &value) in counts.row(0).iter().enumerate() {
        match gc_pct {
            0 => assert!(
                (value - 71.0).abs() < 1e-12,
                "expected 71 counts at GC% 0 after trimmed-span blacklist masking, got {value}"
            ),
            100 => assert!(
                (value - 90.0).abs() < 1e-12,
                "expected 90 counts at GC% 100 after trimmed-span blacklist masking, got {value}"
            ),
            _ => assert!(
                value.abs() < 1e-12,
                "expected 0 counts outside GC% 0/100, got bin {gc_pct}={value}"
            ),
        }
    }

    for gc_pct in 0..=100 {
        let theoretical = gc_pct % 10 == 0;
        assert_eq!(
            support_unobservables[(0, gc_pct)],
            theoretical,
            "unexpected theoretical support at GC% {gc_pct}"
        );
    }

    for gc_pct in 0..=100 {
        let empirical = matches!(gc_pct, 0 | 100);
        assert_eq!(
            support_outliers[(0, gc_pct)],
            empirical,
            "unexpected empirical support at GC% {gc_pct}"
        );
    }

    for gc_pct in 0..=100 {
        let expected_width = if gc_pct % 10 == 0 { 1 } else { 0 };
        assert_eq!(
            gc_percent_widths[(0, gc_pct)],
            expected_width,
            "unexpected GC-percent width at GC% {gc_pct}"
        );
    }

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
fn overlapping_and_touching_bed_windows_with_blacklist_match_explicitly_merged_ref_gc_bias_run()
-> Result<()> {
    // Arrange:
    // `ref-gc-bias` flattens overlapping and touching BED intervals to unique positions before
    // counting. That flattening must stay correct even when a blacklist removes part of the
    // included span afterwards.
    //
    // We compare two logically equivalent descriptions of the same included region:
    // - split BED:  [0, 40), [20, 60), [60, 80)
    // - merged BED: [0, 80)
    //
    // Then we apply the same blacklist [30, 50) to both runs. Because flattening happens before
    // counting and the blacklist is identical, the final written reference package must match
    // exactly.
    let reference = fixtures::simple_reference_twobit()?;
    let split_out_dir = TempDir::new()?;
    let merged_out_dir = TempDir::new()?;
    let bed_dir = TempDir::new()?;
    let blacklist_dir = TempDir::new()?;
    let bed_split = bed_dir.path().join("split_windows.bed");
    let bed_merged = bed_dir.path().join("merged_windows.bed");
    let blacklist_path = blacklist_dir.path().join("blacklist.bed");

    std::fs::write(&bed_split, "chr1\t0\t40\nchr1\t20\t60\nchr1\t60\t80\n")?;
    std::fs::write(&bed_merged, "chr1\t0\t80\n")?;
    std::fs::write(&blacklist_path, "chr1\t30\t50\n")?;

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
        blacklist: Some(vec![blacklist_path.clone()]),
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

    // Act
    run(&make_cfg(split_out_dir.path(), &bed_split))?;
    run(&make_cfg(merged_out_dir.path(), &bed_merged))?;

    // Assert
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
fn full_chromosome_bed_window_with_blacklist_matches_global_ref_gc_bias_run() -> Result<()> {
    // Arrange:
    // `ref-gc-bias` samples start positions first, then counts only the sampled starts that land
    // inside the configured logical region after blacklist masking. A single BED window spanning
    // the full chromosome therefore still describes the same logical region as global mode, even
    // when some bases are later excluded by the blacklist.
    //
    // We keep all execution parameters identical and compare:
    // - global mode + blacklist
    // - one BED window [0, 256) + the same blacklist
    //
    // Because both runs see the same chromosome span, the same fixed seed, the same fragment
    // lengths, and the same blacklisted bases, the written package must match exactly:
    // counts, both support masks, and GC-percent widths.
    let reference = fixtures::simple_reference_twobit()?;
    let global_out_dir = TempDir::new()?;
    let bed_out_dir = TempDir::new()?;
    let bed_dir = TempDir::new()?;
    let blacklist_dir = TempDir::new()?;
    let bed_path = bed_dir.path().join("whole_chr.bed");
    let blacklist_path = blacklist_dir.path().join("blacklist.bed");
    std::fs::write(&bed_path, "chr1\t0\t256\n")?;
    std::fs::write(&blacklist_path, "chr1\t40\t80\n")?;

    let make_global_cfg = |output_dir: &Path| RefGCBiasConfig {
        ref_genome: cfdnalab::commands::cli_common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: output_dir.to_path_buf(),
        n_threads: 1,
        n_positions: 100,
        seed: Some(19),
        windows: Default::default(),
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: Some(vec![blacklist_path.clone()]),
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
fn multiple_blacklist_files_with_touching_intervals_match_single_merged_ref_gc_bias_run()
-> Result<()> {
    // Arrange:
    // `ref-gc-bias` uses the shared blacklist loader with:
    // - `min_size = 1`
    // - `halo_bp = 0`
    //
    // Under those settings, touching intervals from separate files must be merged before any
    // sequence masking happens. So the two-file blacklist:
    //   file A: [40, 60)
    //   file B: [60, 80)
    // must be scientifically identical to the one-file blacklist:
    //   [40, 80)
    //
    // With a fixed seed and otherwise identical config, the full written reference package must
    // match exactly: counts, both support masks, and GC-percent widths.
    let reference = fixtures::simple_reference_twobit()?;
    let split_out_dir = TempDir::new()?;
    let merged_out_dir = TempDir::new()?;
    let blacklist_dir = TempDir::new()?;
    let split_a = blacklist_dir.path().join("blacklist_a.bed");
    let split_b = blacklist_dir.path().join("blacklist_b.bed");
    let merged = blacklist_dir.path().join("blacklist_merged.bed");

    std::fs::write(&split_a, "chr1\t40\t60\n")?;
    std::fs::write(&split_b, "chr1\t60\t80\n")?;
    std::fs::write(&merged, "chr1\t40\t80\n")?;

    let make_cfg = |output_dir: &Path, blacklist: Vec<std::path::PathBuf>| RefGCBiasConfig {
        ref_genome: cfdnalab::commands::cli_common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: output_dir.to_path_buf(),
        n_threads: 1,
        n_positions: 100,
        seed: Some(31),
        windows: Default::default(),
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: Some(blacklist),
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

    // Act
    run(&make_cfg(
        split_out_dir.path(),
        vec![split_a.clone(), split_b.clone()],
    ))?;
    run(&make_cfg(merged_out_dir.path(), vec![merged.clone()]))?;

    // Assert
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
fn fixed_seed_ref_gc_bias_with_blacklist_and_bed_is_invariant_to_tile_size() -> Result<()> {
    // Arrange:
    // Tile size is an execution detail, even when `ref-gc-bias` must do both:
    // - BED-window selection
    // - blacklist masking inside those windows
    //
    // This is a stronger invariance test than the plain global case because the command now has to
    // merge tile-local work while respecting both inclusion windows and excluded bases.
    //
    // We compare two fixed-seed runs over the same sparse BED selection:
    // - tile_size = 128 -> 2 tiles on the 256 bp chromosome
    // - tile_size =  80 -> 4 tiles
    //
    // The BED windows and blacklist are chosen so:
    // - counting is not equivalent to global mode
    // - one included interval is partially blacklisted
    //
    // With the same seed and all other settings identical, the final reference package arrays must
    // still match exactly.
    let reference = fixtures::simple_reference_twobit()?;
    let large_tile_out = TempDir::new()?;
    let small_tile_out = TempDir::new()?;
    let bed_dir = TempDir::new()?;
    let blacklist_dir = TempDir::new()?;
    let bed_path = bed_dir.path().join("selected_windows.bed");
    let blacklist_path = blacklist_dir.path().join("blacklist.bed");
    std::fs::write(&bed_path, "chr1\t0\t80\nchr1\t120\t220\n")?;
    std::fs::write(&blacklist_path, "chr1\t140\t160\n")?;

    let make_cfg = |output_dir: &Path, tile_size: u32| RefGCBiasConfig {
        ref_genome: cfdnalab::commands::cli_common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: output_dir.to_path_buf(),
        n_threads: 1,
        n_positions: 100,
        seed: Some(23),
        windows: cfdnalab::commands::ref_gc_bias::config::RefGCWindowsArgs {
            by_bed: Some(bed_path.clone()),
        },
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: Some(vec![blacklist_path.clone()]),
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

#[test]
fn fixed_seed_ref_gc_bias_with_blacklist_and_bed_is_invariant_to_thread_count() -> Result<()> {
    // Arrange:
    // With a fixed seed, changing only thread count must not change the logical sampled starts or
    // the merged reference package, even when the command simultaneously has to:
    // - limit counting to BED-selected regions
    // - remove blacklisted bases inside those regions
    //
    // We force a multi-tile execution on the 256 bp fixture chromosome by using tile_size = 80,
    // then compare `n_threads = 1` and `n_threads = 2` with all other inputs identical.
    let reference = fixtures::simple_reference_twobit()?;
    let single_thread_out = TempDir::new()?;
    let two_thread_out = TempDir::new()?;
    let bed_dir = TempDir::new()?;
    let blacklist_dir = TempDir::new()?;
    let bed_path = bed_dir.path().join("selected_windows.bed");
    let blacklist_path = blacklist_dir.path().join("blacklist.bed");
    std::fs::write(&bed_path, "chr1\t0\t80\nchr1\t120\t220\n")?;
    std::fs::write(&blacklist_path, "chr1\t140\t160\n")?;

    let make_cfg = |output_dir: &Path, n_threads: usize| RefGCBiasConfig {
        ref_genome: cfdnalab::commands::cli_common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: output_dir.to_path_buf(),
        n_threads,
        n_positions: 100,
        seed: Some(29),
        windows: cfdnalab::commands::ref_gc_bias::config::RefGCWindowsArgs {
            by_bed: Some(bed_path.clone()),
        },
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: Some(vec![blacklist_path.clone()]),
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

    // Act
    run(&make_cfg(single_thread_out.path(), 1))?;
    run(&make_cfg(two_thread_out.path(), 2))?;

    // Assert
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
