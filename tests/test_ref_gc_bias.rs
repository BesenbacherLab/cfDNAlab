#![cfg(feature = "cmd_ref_gc_bias")]

use anyhow::Result;
use cfdnalab::RunOptions;
use ndarray::array;
use std::path::Path;
use tempfile::TempDir;

use cfdnalab::{
    commands::cli_common::{ChromosomeArgs, LoggingArgs},
    commands::gc_bias::counting::{
        GCCounts, build_gc_prefixes, count_reference_gc_and_length_by_window,
        get_gc_integer_percentage_for_window,
    },
    commands::gc_bias::load_reference_bias::load_reference_gc_data,
    commands::gc_bias::support_masking::create_support_mask_threshold_per_mb,
    run_like_cli::ref_gc_bias::{RefGCBiasConfig, run_ref_gc_bias},
    shared::{
        blacklist::apply_blacklist_mask_to_seq,
        interval::{IndexedInterval, Interval},
        reference::twobit_contig_footprint,
    },
};
mod fixtures;

fn run(cfg: &RefGCBiasConfig) -> Result<()> {
    run_ref_gc_bias(cfg, RunOptions::new_quiet()).map(|_| ())
}

fn load_ref_gc_package_arrays(
    package_path: &Path,
) -> Result<(
    ndarray::Array2<f64>,
    ndarray::Array2<bool>,
    ndarray::Array2<bool>,
    ndarray::Array2<u16>,
)> {
    let loaded = load_reference_gc_data(package_path)?;
    Ok((
        loaded.counts,
        loaded.unobservables_support_mask,
        loaded.outliers_support_mask,
        loaded.gc_percent_widths,
    ))
}

#[test]
fn gc_prefix_helpers_return_prefix_differences_for_checked_intervals() -> Result<()> {
    // Sequence A C N G T A:
    // - [1,5) = C N G T -> GC=2, ACGT=3
    // - [2,4) = N G     -> GC=1, ACGT=1
    let seq = b"ACNGTA".to_vec();
    let prefixes = build_gc_prefixes(&seq);
    let full_interval = Interval::new(1usize, 5usize)?;
    let inner_interval = Interval::new(2usize, 4usize)?;

    assert_eq!(prefixes.gc_count(full_interval)?, 2);
    assert_eq!(prefixes.acgt_count(full_interval)?, 3);

    assert_eq!(prefixes.gc_count(inner_interval)?, 1);
    assert_eq!(prefixes.acgt_count(inner_interval)?, 1);
    Ok(())
}

#[test]
fn gc_prefix_helpers_error_when_interval_exceeds_prefix_bounds() -> Result<()> {
    // Sequence A C G T has prefix length 5, so the largest valid half-open interval end is 4.
    // Asking for [1,5) should therefore fail before any subtraction is attempted.
    let seq = b"ACGT".to_vec();
    let prefixes = build_gc_prefixes(&seq);
    let invalid_interval = Interval::new(1usize, 5usize)?;

    let gc_err = prefixes
        .gc_count(invalid_interval)
        .expect_err("expected GC prefix bounds error");
    let acgt_err = prefixes
        .acgt_count(invalid_interval)
        .expect_err("expected ACGT prefix bounds error");

    assert!(
        gc_err
            .to_string()
            .contains("GC interval [1, 5) out of bounds"),
        "unexpected GC error: {gc_err}"
    );
    assert!(
        acgt_err
            .to_string()
            .contains("ACGT interval [1, 5) out of bounds"),
        "unexpected ACGT error: {acgt_err}"
    );
    Ok(())
}

#[test]
fn gc_integer_percentage_window_returns_some_none_and_error() -> Result<()> {
    // Sequence A C N G T:
    // - [0,5) = A C N G T -> GC=2, ACGT=4, so GC% = round(200/4) = 50
    // - [1,4) = C N G     -> GC=2, ACGT=2, so this stays valid when min_acgt_count=2
    // - [2,4) = N G       -> GC=1, ACGT=1, so min_acgt_count=2 should return Ok(None)
    // - [3,6) extends past the prefix arrays for a 5 bp sequence, so it should error
    let seq = b"ACNGT".to_vec();
    let prefixes = build_gc_prefixes(&seq);

    let full_window = Interval::new(0usize, 5usize)?;
    let low_support_window = Interval::new(2usize, 4usize)?;
    let invalid_window = Interval::new(3usize, 6usize)?;

    assert_eq!(
        get_gc_integer_percentage_for_window(&prefixes, full_window, 0.0, 1)?,
        Some(50)
    );
    assert_eq!(
        get_gc_integer_percentage_for_window(&prefixes, low_support_window, 0.0, 2)?,
        None
    );

    let err = get_gc_integer_percentage_for_window(&prefixes, invalid_window, 0.0, 1)
        .expect_err("expected out-of-bounds GC window error");
    assert!(
        err.to_string().contains("GC interval [3, 6) out of bounds"),
        "unexpected GC window error: {err}"
    );
    Ok(())
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
    )?;

    // Assert:
    // Window 0 counts fragments starting at 0 inside `ACGTAC`:
    // - len4 trimmed to `CG`   -> gc=2
    // - len5 trimmed to `CGT`  -> gc=2
    // - len6 trimmed to `CGTA` -> gc=2
    //
    // Window 1 counts fragments starting at 6 inside `GTACGT`:
    // - len4 trimmed to `TA`   -> gc=0
    // - len5 trimmed to `TAC`  -> gc=1
    // - len6 trimmed to `TACG` -> gc=2
    //
    // No other `(length, gc)` cells should receive any counts.
    let expected_window0 = &[(4_usize, 2_usize), (5, 2), (6, 2)];
    let expected_window1 = &[(4_usize, 0_usize), (5, 1), (6, 2)];

    for (window_index, (window_counts, expected_non_zero_cells)) in counts_by_bin
        .iter()
        .zip([expected_window0, expected_window1].iter())
        .enumerate()
    {
        assert_eq!(
            window_counts.sum(),
            expected_non_zero_cells.len() as f64,
            "window {window_index} should contain exactly one count for each tested length"
        );

        for length in 4..=6 {
            let effective_length = length - 2;
            for gc_count in 0..=effective_length {
                let expected_value = if expected_non_zero_cells.contains(&(length, gc_count)) {
                    1.0
                } else {
                    0.0
                };
                assert_eq!(
                    window_counts.get(length, gc_count).unwrap(),
                    expected_value,
                    "window {window_index} expected count {expected_value} at length {length}, gc {gc_count}"
                );
            }
        }
    }
    Ok(())
}

#[test]
fn reference_gc_counts_use_tile_local_prefix_coordinates_after_late_sequence_load() -> Result<()> {
    // This mirrors ref-gc-bias after loading a late reference slice, e.g. absolute [900,964).
    // The count helper receives tile-local coordinates, so absolute window [930,941) is [30,41)
    // and absolute start 930 is local start 30.
    //
    // In the ACGT repeat, local [30,41) is:
    //   G T A C G T A C G T A
    // This has 5 GC bases and 11 ACGT bases, so length 11 / GC count 5 gets one count.
    let seq = b"ACGT".repeat(16);
    let prefixes = build_gc_prefixes(&seq);
    let windows = IndexedInterval::from_tuples(&[(30, 41, 0)])?;
    let starts = vec![30usize];
    let mut counts_by_bin = vec![GCCounts::new(11, 11, 0, (0, 0))?];

    // Act
    count_reference_gc_and_length_by_window(
        &mut counts_by_bin,
        &prefixes,
        (11, 12),
        windows.as_slice(),
        starts.as_slice(),
        seq.len() as u64,
        1.0,
        1,
        0,
    )?;

    // Assert
    assert_eq!(counts_by_bin[0].sum(), 1.0);
    assert_eq!(
        counts_by_bin[0]
            .get(11, 5)
            .expect("length 11 / GC count 5 should be in range"),
        1.0
    );
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
    )?;

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
        (99_999_999_usize, 1_usize, vec![false, true, true, true]),
        (100_000_000_usize, 2_usize, vec![false, false, true, true]),
        (200_000_000_usize, 3_usize, vec![false, false, false, true]),
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

// MOVE-MODULE-LOCAL: direct GC prefix/counting/support helper tests above this point.
// REWRITE-PUBLIC-TEST: ref-gc-bias artifact behavior currently uses a private package loader.
#[test]
fn ref_gc_bias_run_writes_expected_prefixed_package_metadata_and_shapes() -> Result<()> {
    let reference = fixtures::simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let output_prefix = "unit_ref_gc";

    let cfg = RefGCBiasConfig {
        ref_genome: cfdnalab::run_like_cli::common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: out_dir.path().to_path_buf(),
        output_prefix: output_prefix.to_string(),
        n_threads: 1,
        n_positions: 100,
        seed: Some(7),
        windows: Default::default(),
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: None,
        fragment_lengths: cfdnalab::run_like_cli::common::FragmentLengthArgs {
            min_fragment_length: 10,
            max_fragment_length: 12,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 1_000_000,
        logging: LoggingArgs::default(),
    };
    cfg.check_smoothing_settings()?;

    // Manual expectations:
    // - We request fragment lengths 10, 11, and 12, so the output must have 3 rows.
    // - The package stores GC as integer percentages 0..=100, so the counts and both support masks
    //   must have 101 columns.
    // - We explicitly disable smoothing and interpolation, so the written scalar metadata must
    //   reflect those choices exactly.
    // - With `output_prefix`, the command should write only the prefixed package path.
    // - This is a command-level test of the written artifact contract, not of the exact sampled counts.
    run(&cfg)?;

    let package_path = out_dir
        .path()
        .join(format!("{output_prefix}.ref_gc_package.zarr"));
    assert!(
        !out_dir.path().join("ref_gc_package.zarr").exists(),
        "Did not expect unprefixed package when output_prefix is set"
    );
    let loaded = load_reference_gc_data(&package_path)?;
    let counts = loaded.counts;
    let support_unobservables = loaded.unobservables_support_mask;
    let support_outliers = loaded.outliers_support_mask;
    let gc_percent_widths = loaded.gc_percent_widths;
    let metadata = loaded.metadata;

    assert_eq!(counts.dim(), (3, 101));
    assert_eq!(support_unobservables.dim(), (3, 101));
    assert_eq!(support_outliers.dim(), (3, 101));
    assert_eq!(gc_percent_widths.dim(), (3, 101));

    assert_eq!(metadata.min_fragment_length, 10);
    assert_eq!(metadata.max_fragment_length, 12);
    assert_eq!(metadata.end_offset, 0);
    assert!(metadata.skip_interpolation);
    assert_eq!(metadata.smoothing_radius, 2);
    assert_eq!(metadata.smoothing_sigma, 0.55);
    assert!(metadata.skip_smoothing);
    assert_eq!(metadata.chromosomes, vec!["chr1".to_string()]);
    assert_eq!(
        metadata.reference_contig_footprint,
        twobit_contig_footprint(&reference.path)?
    );

    let expected_theoretical_bins = [
        vec![0_usize, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
        vec![0_usize, 9, 18, 27, 36, 45, 55, 64, 73, 82, 91, 100],
        vec![0_usize, 8, 17, 25, 33, 42, 50, 58, 67, 75, 83, 92, 100],
    ];
    for (row_index, row) in counts.outer_iter().enumerate() {
        let supported_bins: &[usize] = match row_index {
            0 => &[40, 50, 60], // length 10 -> GC counts 4, 5, 6
            1 => &[45, 55],     // length 11 -> GC counts 5 or 6
            2 => &[50],         // length 12 -> always 50%
            _ => panic!("unexpected row index {row_index}"),
        };
        let theoretical_bins = &expected_theoretical_bins[row_index];

        assert!(
            row.iter().all(|value| value.is_finite() && *value >= 0.0),
            "reference GC counts should be finite and non-negative in row {row_index}"
        );
        assert_eq!(
            row.sum(),
            100.0,
            "each sampled start should contribute one count per fragment length"
        );

        for (gc_percent, &value) in row.iter().enumerate() {
            if !supported_bins.contains(&gc_percent) {
                assert_eq!(
                    value, 0.0,
                    "row {row_index} should not place counts outside the reachable GC% bins"
                );
            }
        }

        let unobservable_row = support_unobservables.row(row_index);
        let outlier_row = support_outliers.row(row_index);
        let widths_row = gc_percent_widths.row(row_index);
        for gc_percent in 0..=100 {
            let expected_theoretical_support = theoretical_bins.contains(&gc_percent);
            let expected_empirical_support = supported_bins.contains(&gc_percent);
            let expected_width = if expected_theoretical_support {
                1_u16
            } else {
                0_u16
            };

            assert_eq!(
                unobservable_row[gc_percent], expected_theoretical_support,
                "row {row_index} theoretical support mismatch at GC% {gc_percent}"
            );
            assert_eq!(
                outlier_row[gc_percent], expected_empirical_support,
                "row {row_index} empirical support mismatch at GC% {gc_percent}"
            );
            assert_eq!(
                widths_row[gc_percent], expected_width,
                "row {row_index} GC-percent width mismatch at GC% {gc_percent}"
            );
        }
    }

    Ok(())
}

// REWRITE-PUBLIC-TEST: ref-gc-bias artifact behavior currently uses a private package loader.
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
    //
    // Important command detail:
    // - `ref-gc-bias` does not interpret BED windows as "allowed start positions".
    // - It calls `count_reference_gc_and_length_by_window`, which requires the full fragment to
    //   fit inside the window:
    //       frag_len <= window_end - start_pos
    // - So for `frag_len = 10`, window [0, 91) counts starts 0..=81, not 0..=90.
    // - Likewise [100, 191) counts starts 100..=181, not 100..=190.
    //
    // Those counted fragments are still pure:
    // - starts 0..=81 stay fully inside the A block
    // - starts 100..=181 stay fully inside the C block
    // The boundary-crossing starts are sampled, but excluded because the fragment would not fit
    // inside the BED window.
    //
    // Therefore the output counts are exactly:
    // - 82 counts at GC% 0
    // - 82 counts at GC% 100
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
        ref_genome: cfdnalab::run_like_cli::common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: out_dir.path().to_path_buf(),
        output_prefix: String::new(),
        n_threads: 1,
        n_positions: 191,
        seed: Some(23),
        windows: cfdnalab::run_like_cli::ref_gc_bias::RefGCWindowsArgs {
            by_bed: Some(bed_path),
        },
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: None,
        fragment_lengths: cfdnalab::run_like_cli::common::FragmentLengthArgs {
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

    // Act
    run(&cfg)?;

    // Assert
    let package_path = out_dir.path().join("ref_gc_package.zarr");
    let (counts, support_unobservables, support_outliers, gc_percent_widths) =
        load_ref_gc_package_arrays(&package_path)?;

    assert_eq!(counts.dim(), (1, 101));
    assert_eq!(support_unobservables.dim(), (1, 101));
    assert_eq!(support_outliers.dim(), (1, 101));
    assert_eq!(gc_percent_widths.dim(), (1, 101));

    for (gc_pct, &value) in counts.row(0).iter().enumerate() {
        match gc_pct {
            0 | 100 => assert!(
                (value - 82.0).abs() < 1e-12,
                "expected 82 counts at GC% {gc_pct}, got {value}"
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

// REWRITE-PUBLIC-TEST: ref-gc-bias artifact behavior currently uses a private package loader.
#[test]
fn ref_gc_bias_run_blacklist_removes_exactly_the_overlapping_start_positions() -> Result<()> {
    // Arrange:
    // Start from the same exact-count setup as the previous test:
    // - chr1[0,100)   = all A
    // - chr1[100,200) = all C
    // - fragment length fixed at 10
    // - every valid start sampled exactly once (`n_positions = 191`)
    // - BED windows are still interpreted as full fragment containers, not start-only filters:
    //     [0, 91)   with length 10 counts starts 0..=81
    //     [100,191) with length 10 counts starts 100..=181
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
    // Those 19 dropped starts all belong to the counted pure-A starts 0..=81.
    // The counted pure-C starts 100..=181 never touch the blacklist.
    //
    // Therefore the exact output counts must be:
    // - GC% 0   -> 82 - 19 = 63
    // - GC% 100 -> 82
    // - 0 everywhere else
    //
    // The empirical support threshold is still tiny:
    //   total covered ACGT bases = 63 + 82 = 145
    //   threshold = 145 / 1_000_000 = 0.000145
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
        ref_genome: cfdnalab::run_like_cli::common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: out_dir.path().to_path_buf(),
        output_prefix: String::new(),
        n_threads: 1,
        n_positions: 191,
        seed: Some(23),
        windows: cfdnalab::run_like_cli::ref_gc_bias::RefGCWindowsArgs {
            by_bed: Some(bed_path),
        },
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: Some(vec![blacklist_path]),
        fragment_lengths: cfdnalab::run_like_cli::common::FragmentLengthArgs {
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

    // Act
    run(&cfg)?;

    // Assert
    let package_path = out_dir.path().join("ref_gc_package.zarr");
    let (counts, support_unobservables, support_outliers, gc_percent_widths) =
        load_ref_gc_package_arrays(&package_path)?;

    assert_eq!(counts.dim(), (1, 101));
    assert_eq!(support_unobservables.dim(), (1, 101));
    assert_eq!(support_outliers.dim(), (1, 101));
    assert_eq!(gc_percent_widths.dim(), (1, 101));

    for (gc_pct, &value) in counts.row(0).iter().enumerate() {
        match gc_pct {
            0 => assert!(
                (value - 63.0).abs() < 1e-12,
                "expected 63 counts at GC% 0 after blacklist masking, got {value}"
            ),
            100 => assert!(
                (value - 82.0).abs() < 1e-12,
                "expected 82 counts at GC% 100 after blacklist masking, got {value}"
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

// REWRITE-PUBLIC-TEST: ref-gc-bias artifact behavior currently uses a private package loader.
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
    // To keep only pure trimmed spans we choose BED windows:
    //   [0, 90) and [99, 189)
    //
    // Important command detail:
    // - even with `end_offset = 1`, BED windows are still checked against the full raw fragment
    //   span, not only against the trimmed GC-counting span.
    // - So with raw fragment length 12:
    //     [0, 90)   counts starts 0..=78
    //     [99, 189) counts starts 99..=177
    //
    // Those counted fragments still have pure trimmed spans:
    // - start 78 trims to [79, 89), still inside the A block
    // - start 177 trims to [178, 188), still inside the C block
    //
    // Starts closer to the boundary are sampled, but excluded because the raw 12 bp fragment would
    // not fit fully inside the BED window.
    //
    // Therefore the exact output counts must be:
    // - 79 counts at GC% 0
    // - 79 counts at GC% 100
    // - 0 everywhere else
    //
    // The empirical support threshold is again tiny:
    //   threshold_per_mb = 1 + 189 / 100_000_000 = 1
    //   total covered ACGT bases = 79 + 79 = 158
    //   threshold = 158 / 1_000_000 = 0.000158
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
        ref_genome: cfdnalab::run_like_cli::common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: out_dir.path().to_path_buf(),
        output_prefix: String::new(),
        n_threads: 1,
        n_positions: 189,
        seed: Some(37),
        windows: cfdnalab::run_like_cli::ref_gc_bias::RefGCWindowsArgs {
            by_bed: Some(bed_path),
        },
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: None,
        fragment_lengths: cfdnalab::run_like_cli::common::FragmentLengthArgs {
            min_fragment_length: 12,
            max_fragment_length: 12,
        },
        end_offset: 1,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 1_000_000,
        logging: LoggingArgs::default(),
    };

    // Act
    run(&cfg)?;

    // Assert
    let package_path = out_dir.path().join("ref_gc_package.zarr");
    let (counts, support_unobservables, support_outliers, gc_percent_widths) =
        load_ref_gc_package_arrays(&package_path)?;

    assert_eq!(counts.dim(), (1, 101));
    assert_eq!(support_unobservables.dim(), (1, 101));
    assert_eq!(support_outliers.dim(), (1, 101));
    assert_eq!(gc_percent_widths.dim(), (1, 101));

    for (gc_pct, &value) in counts.row(0).iter().enumerate() {
        match gc_pct {
            0 | 100 => assert!(
                (value - 79.0).abs() < 1e-12,
                "expected 79 counts at GC% {gc_pct}, got {value}"
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

// REWRITE-PUBLIC-TEST: ref-gc-bias artifact behavior currently uses a private package loader.
#[test]
fn ref_gc_bias_run_blacklist_with_end_offset_drops_only_trimmed_overlaps() -> Result<()> {
    // Arrange:
    // Reuse the trimmed-span setup from the previous test:
    // - chr1[0,100)   = all A
    // - chr1[100,200) = all C
    // - fragment length = 12
    // - end_offset = 1
    // - every valid start sampled exactly once (`n_positions = 189`)
    // - BED windows still require the full raw 12 bp fragment to fit:
    //     [0,90)   -> starts 0..=78
    //     [99,189) -> starts 99..=177
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
    // So exactly 19 counted pure-A starts are dropped. The counted pure-C starts 99..=177 remain
    // untouched.
    //
    // This is intentionally different from raw-span overlap:
    // - raw 12 bp fragment overlap would drop starts 34..=54 (21 starts)
    // - the command should *not* do that here
    //
    // Therefore the exact output counts must be:
    // - GC% 0   -> 79 - 19 = 60
    // - GC% 100 -> 79
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
        ref_genome: cfdnalab::run_like_cli::common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: out_dir.path().to_path_buf(),
        output_prefix: String::new(),
        n_threads: 1,
        n_positions: 189,
        seed: Some(37),
        windows: cfdnalab::run_like_cli::ref_gc_bias::RefGCWindowsArgs {
            by_bed: Some(bed_path),
        },
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: Some(vec![blacklist_path]),
        fragment_lengths: cfdnalab::run_like_cli::common::FragmentLengthArgs {
            min_fragment_length: 12,
            max_fragment_length: 12,
        },
        end_offset: 1,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 1_000_000,
        logging: LoggingArgs::default(),
    };

    // Act
    run(&cfg)?;

    // Assert
    let package_path = out_dir.path().join("ref_gc_package.zarr");
    let (counts, support_unobservables, support_outliers, gc_percent_widths) =
        load_ref_gc_package_arrays(&package_path)?;

    assert_eq!(counts.dim(), (1, 101));
    assert_eq!(support_unobservables.dim(), (1, 101));
    assert_eq!(support_outliers.dim(), (1, 101));
    assert_eq!(gc_percent_widths.dim(), (1, 101));

    for (gc_pct, &value) in counts.row(0).iter().enumerate() {
        match gc_pct {
            0 => assert!(
                (value - 60.0).abs() < 1e-12,
                "expected 60 counts at GC% 0 after trimmed-span blacklist masking, got {value}"
            ),
            100 => assert!(
                (value - 79.0).abs() < 1e-12,
                "expected 79 counts at GC% 100 after trimmed-span blacklist masking, got {value}"
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

// REWRITE-PUBLIC-TEST: ref-gc-bias artifact behavior currently uses a private package loader.
#[test]
fn ref_gc_bias_run_smoothing_enabled_spreads_three_gc_anchors_by_known_kernel() -> Result<()> {
    // Arrange:
    // Build three isolated 10 bp windows that each admit exactly one fragment start:
    // - [0,10)   = AAAAAAAAAA -> GC%=0
    // - [20,30)  = CCCCCAAAAA -> GC%=50
    // - [40,50)  = CCCCCCCCCC -> GC%=100
    // - plus two unused trailing A bases so the chromosome length is 52 bp instead of 50,
    //   avoiding the upstream `.2bit` partial-byte tail bug while keeping `[40,50)` intact
    //
    // The 10 bp filler blocks between them are all T so no extra counted starts can leak in.
    // With fragment length 10 and BED rows of width 10, each row contributes exactly one start,
    // so the raw GC-count row before smoothing is:
    //   gc_count 0  -> 1
    //   gc_count 5  -> 1
    //   gc_count 10 -> 1
    //
    // Choose `sigma = sqrt(1 / (2 ln 2))` and `radius = 1`.
    // Then the unnormalized Gaussian weights are:
    //   exp(-1 / (2 sigma^2)) = 1/2
    // so the normalized 3-tap kernel is exactly:
    //   [1/4, 1/2, 1/4]
    //
    // Smoothing that sparse row gives:
    //   gc_count 0  -> 1/2
    //   gc_count 1  -> 1/4
    //   gc_count 4  -> 1/4
    //   gc_count 5  -> 1/2
    //   gc_count 6  -> 1/4
    //   gc_count 9  -> 1/4
    //   gc_count 10 -> 1/2
    //
    // For effective length 10, GC counts map exactly to GC percentages in steps of 10, and the
    // width correction is 1 for every reachable GC% bin. So the written counts row must carry the
    // same values at GC% 0,10,40,50,60,90,100 and zero elsewhere.
    let sigma = (1.0_f64 / (2.0 * std::f64::consts::LN_2)).sqrt();
    let reference = fixtures::twobit_from_sequences(
        "ref_gc_bias_smoothing_three_anchor_reference",
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
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("windows.bed");
    std::fs::write(&bed_path, "chr1\t0\t10\nchr1\t20\t30\nchr1\t40\t50\n")?;

    let cfg = RefGCBiasConfig {
        ref_genome: cfdnalab::run_like_cli::common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: out_dir.path().to_path_buf(),
        output_prefix: String::new(),
        n_threads: 1,
        n_positions: 41,
        seed: Some(11),
        windows: cfdnalab::run_like_cli::ref_gc_bias::RefGCWindowsArgs {
            by_bed: Some(bed_path),
        },
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: None,
        fragment_lengths: cfdnalab::run_like_cli::common::FragmentLengthArgs {
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

    // Act
    run(&cfg)?;

    // Assert
    let package_path = out_dir.path().join("ref_gc_package.zarr");
    let (counts, support_unobservables, support_outliers, gc_percent_widths) =
        load_ref_gc_package_arrays(&package_path)?;

    assert_eq!(counts.dim(), (1, 101));
    assert_eq!(support_unobservables.dim(), (1, 101));
    assert_eq!(support_outliers.dim(), (1, 101));
    assert_eq!(gc_percent_widths.dim(), (1, 101));

    let expected_non_zero = [
        (0_usize, 0.5_f64),
        (10, 0.25),
        (40, 0.25),
        (50, 0.5),
        (60, 0.25),
        (90, 0.25),
        (100, 0.5),
    ];
    for (gc_pct, &value) in counts.row(0).iter().enumerate() {
        let expected_value = expected_non_zero
            .iter()
            .find_map(|(expected_gc_pct, expected_value)| {
                (*expected_gc_pct == gc_pct).then_some(*expected_value)
            })
            .unwrap_or(0.0);
        assert!(
            (value - expected_value).abs() <= 1e-6,
            "expected smoothed count {expected_value} at GC% {gc_pct}, got {value}"
        );
    }

    for gc_pct in 0..=100 {
        let expected_empirical = matches!(gc_pct, 0 | 10 | 40 | 50 | 60 | 90 | 100);
        assert_eq!(
            support_outliers[(0, gc_pct)],
            expected_empirical,
            "unexpected empirical support at GC% {gc_pct}"
        );
        let expected_width = if gc_pct % 10 == 0 { 1 } else { 0 };
        assert_eq!(
            gc_percent_widths[(0, gc_pct)],
            expected_width,
            "unexpected GC-percent width at GC% {gc_pct}"
        );
    }

    Ok(())
}

// REWRITE-PUBLIC-TEST: ref-gc-bias artifact behavior currently uses a private package loader.
#[test]
fn ref_gc_bias_run_interpolation_enabled_fills_between_equal_supported_anchors() -> Result<()> {
    // Arrange:
    // Reuse the same three isolated starts as above, but now disable smoothing and enable
    // interpolation. The pre-interpolation row is therefore exactly:
    //   GC% 0   -> 1
    //   GC% 50  -> 1
    //   GC% 100 -> 1
    // and zero elsewhere.
    // The chromosome again includes two unused trailing A bases so `[40,50)` stays pure-C without
    // depending on the upstream `.2bit` partial-byte tail behavior.
    //
    // The empirical support mask is true only at those three anchor bins. With three equal anchors,
    // the fitted quadratic is the constant function 1.0, so interpolation must fill every
    // unsupported GC% bin with 1.0 while leaving the support mask unchanged.
    let reference = fixtures::twobit_from_sequences(
        "ref_gc_bias_interpolation_three_anchor_reference",
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
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("windows.bed");
    std::fs::write(&bed_path, "chr1\t0\t10\nchr1\t20\t30\nchr1\t40\t50\n")?;

    let cfg = RefGCBiasConfig {
        ref_genome: cfdnalab::run_like_cli::common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: out_dir.path().to_path_buf(),
        output_prefix: String::new(),
        n_threads: 1,
        n_positions: 41,
        seed: Some(11),
        windows: cfdnalab::run_like_cli::ref_gc_bias::RefGCWindowsArgs {
            by_bed: Some(bed_path),
        },
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: None,
        fragment_lengths: cfdnalab::run_like_cli::common::FragmentLengthArgs {
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

    // Act
    run(&cfg)?;

    // Assert
    let package_path = out_dir.path().join("ref_gc_package.zarr");
    let (counts, support_unobservables, support_outliers, gc_percent_widths) =
        load_ref_gc_package_arrays(&package_path)?;

    assert_eq!(counts.dim(), (1, 101));
    assert_eq!(support_unobservables.dim(), (1, 101));
    assert_eq!(support_outliers.dim(), (1, 101));
    assert_eq!(gc_percent_widths.dim(), (1, 101));

    for (gc_pct, &value) in counts.row(0).iter().enumerate() {
        assert!(
            (value - 1.0).abs() <= 1e-6,
            "interpolation should fill GC% {gc_pct} to 1.0, got {value}"
        );
        let expected_empirical = matches!(gc_pct, 0 | 50 | 100);
        assert_eq!(
            support_outliers[(0, gc_pct)],
            expected_empirical,
            "interpolation must not rewrite the empirical support mask at GC% {gc_pct}"
        );
    }

    Ok(())
}

// REWRITE-PUBLIC-TEST: ref-gc-bias artifact behavior currently uses a private package loader.
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
        ref_genome: cfdnalab::run_like_cli::common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: output_dir.to_path_buf(),
        output_prefix: String::new(),
        n_threads: 1,
        n_positions: 100,
        seed: Some(11),
        windows: cfdnalab::run_like_cli::ref_gc_bias::RefGCWindowsArgs {
            by_bed: Some(bed_path.to_path_buf()),
        },
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: None,
        fragment_lengths: cfdnalab::run_like_cli::common::FragmentLengthArgs {
            min_fragment_length: 10,
            max_fragment_length: 12,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 1_000_000,
        logging: LoggingArgs::default(),
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

    let split_package = split_out_dir.path().join("ref_gc_package.zarr");
    let merged_package = merged_out_dir.path().join("ref_gc_package.zarr");
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

// REWRITE-PUBLIC-TEST: ref-gc-bias artifact behavior currently uses a private package loader.
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
        ref_genome: cfdnalab::run_like_cli::common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: output_dir.to_path_buf(),
        output_prefix: String::new(),
        n_threads: 1,
        n_positions: 100,
        seed: Some(11),
        windows: cfdnalab::run_like_cli::ref_gc_bias::RefGCWindowsArgs {
            by_bed: Some(bed_path.to_path_buf()),
        },
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: Some(vec![blacklist_path.clone()]),
        fragment_lengths: cfdnalab::run_like_cli::common::FragmentLengthArgs {
            min_fragment_length: 10,
            max_fragment_length: 12,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 80,
        logging: LoggingArgs::default(),
    };

    // Act
    run(&make_cfg(split_out_dir.path(), &bed_split))?;
    run(&make_cfg(merged_out_dir.path(), &bed_merged))?;

    // Assert
    let split_package = split_out_dir.path().join("ref_gc_package.zarr");
    let merged_package = merged_out_dir.path().join("ref_gc_package.zarr");
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

// REWRITE-PUBLIC-TEST: ref-gc-bias artifact behavior currently uses a private package loader.
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
        ref_genome: cfdnalab::run_like_cli::common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: output_dir.to_path_buf(),
        output_prefix: String::new(),
        n_threads: 1,
        n_positions: 100,
        seed: Some(13),
        windows: Default::default(),
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: None,
        fragment_lengths: cfdnalab::run_like_cli::common::FragmentLengthArgs {
            min_fragment_length: 10,
            max_fragment_length: 12,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 80,
        logging: LoggingArgs::default(),
    };
    let make_bed_cfg = |output_dir: &Path| RefGCBiasConfig {
        windows: cfdnalab::run_like_cli::ref_gc_bias::RefGCWindowsArgs {
            by_bed: Some(bed_path.clone()),
        },
        ..make_global_cfg(output_dir)
    };

    // Act
    run(&make_global_cfg(global_out_dir.path()))?;
    run(&make_bed_cfg(bed_out_dir.path()))?;

    // Assert
    let global_package = global_out_dir.path().join("ref_gc_package.zarr");
    let bed_package = bed_out_dir.path().join("ref_gc_package.zarr");
    let (
        global_counts,
        global_support_unobservables,
        global_support_outliers,
        global_gc_percent_widths,
    ) = load_ref_gc_package_arrays(&global_package)?;
    let (bed_counts, bed_support_unobservables, bed_support_outliers, bed_gc_percent_widths) =
        load_ref_gc_package_arrays(&bed_package)?;

    assert_eq!(global_counts, bed_counts);
    assert_eq!(global_support_unobservables, bed_support_unobservables);
    assert_eq!(global_support_outliers, bed_support_outliers);
    assert_eq!(global_gc_percent_widths, bed_gc_percent_widths);

    Ok(())
}

// REWRITE-PUBLIC-TEST: ref-gc-bias artifact behavior currently uses a private package loader.
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
        ref_genome: cfdnalab::run_like_cli::common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: output_dir.to_path_buf(),
        output_prefix: String::new(),
        n_threads: 1,
        n_positions: 100,
        seed: Some(19),
        windows: Default::default(),
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: Some(vec![blacklist_path.clone()]),
        fragment_lengths: cfdnalab::run_like_cli::common::FragmentLengthArgs {
            min_fragment_length: 10,
            max_fragment_length: 12,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 80,
        logging: LoggingArgs::default(),
    };
    let make_bed_cfg = |output_dir: &Path| RefGCBiasConfig {
        windows: cfdnalab::run_like_cli::ref_gc_bias::RefGCWindowsArgs {
            by_bed: Some(bed_path.clone()),
        },
        ..make_global_cfg(output_dir)
    };

    // Act
    run(&make_global_cfg(global_out_dir.path()))?;
    run(&make_bed_cfg(bed_out_dir.path()))?;

    // Assert
    let global_package = global_out_dir.path().join("ref_gc_package.zarr");
    let bed_package = bed_out_dir.path().join("ref_gc_package.zarr");
    let (
        global_counts,
        global_support_unobservables,
        global_support_outliers,
        global_gc_percent_widths,
    ) = load_ref_gc_package_arrays(&global_package)?;
    let (bed_counts, bed_support_unobservables, bed_support_outliers, bed_gc_percent_widths) =
        load_ref_gc_package_arrays(&bed_package)?;

    assert_eq!(global_counts, bed_counts);
    assert_eq!(global_support_unobservables, bed_support_unobservables);
    assert_eq!(global_support_outliers, bed_support_outliers);
    assert_eq!(global_gc_percent_widths, bed_gc_percent_widths);

    Ok(())
}

// REWRITE-PUBLIC-TEST: ref-gc-bias artifact behavior currently uses a private package loader.
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
        ref_genome: cfdnalab::run_like_cli::common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: output_dir.to_path_buf(),
        output_prefix: String::new(),
        n_threads: 1,
        n_positions: 100,
        seed: Some(31),
        windows: Default::default(),
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: Some(blacklist),
        fragment_lengths: cfdnalab::run_like_cli::common::FragmentLengthArgs {
            min_fragment_length: 10,
            max_fragment_length: 12,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 80,
        logging: LoggingArgs::default(),
    };

    // Act
    run(&make_cfg(
        split_out_dir.path(),
        vec![split_a.clone(), split_b.clone()],
    ))?;
    run(&make_cfg(merged_out_dir.path(), vec![merged.clone()]))?;

    // Assert
    let split_package = split_out_dir.path().join("ref_gc_package.zarr");
    let merged_package = merged_out_dir.path().join("ref_gc_package.zarr");
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

// REWRITE-PUBLIC-TEST: ref-gc-bias artifact behavior currently uses a private package loader.
#[test]
fn rejects_n_positions_when_sampling_density_would_exceed_one() -> Result<()> {
    let reference = fixtures::simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let cfg = RefGCBiasConfig {
        ref_genome: cfdnalab::run_like_cli::common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: out_dir.path().to_path_buf(),
        output_prefix: String::new(),
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
        fragment_lengths: cfdnalab::run_like_cli::common::FragmentLengthArgs {
            min_fragment_length: 60,
            max_fragment_length: 60,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 1_000_000,
        logging: LoggingArgs::default(),
    };

    let err = run(&cfg).expect_err("sampling density above 1.0 should fail");
    let err_text = err.to_string();
    assert!(
        err_text.contains("Sampling density") && err_text.contains("exceeds 1.0"),
        "unexpected error message: {err_text}"
    );

    Ok(())
}

// REWRITE-PUBLIC-TEST: ref-gc-bias artifact behavior currently uses a private package loader.
#[test]
fn fixed_seed_ref_gc_bias_is_invariant_to_thread_count() -> Result<()> {
    let reference = fixtures::simple_reference_twobit()?;
    let single_thread_out = TempDir::new()?;
    let two_thread_out = TempDir::new()?;

    let make_cfg = |output_dir: &Path, n_threads: usize| RefGCBiasConfig {
        ref_genome: cfdnalab::run_like_cli::common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: output_dir.to_path_buf(),
        output_prefix: String::new(),
        n_threads,
        n_positions: 100,
        seed: Some(7),
        windows: Default::default(),
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: None,
        fragment_lengths: cfdnalab::run_like_cli::common::FragmentLengthArgs {
            min_fragment_length: 10,
            max_fragment_length: 12,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 80,
        logging: LoggingArgs::default(),
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

    let single_package = single_thread_out.path().join("ref_gc_package.zarr");
    let two_package = two_thread_out.path().join("ref_gc_package.zarr");
    let (
        single_counts,
        single_support_unobservables,
        single_support_outliers,
        single_gc_percent_widths,
    ) = load_ref_gc_package_arrays(&single_package)?;
    let (two_counts, two_support_unobservables, two_support_outliers, two_gc_percent_widths) =
        load_ref_gc_package_arrays(&two_package)?;

    assert_eq!(single_counts, two_counts);
    assert_eq!(single_support_unobservables, two_support_unobservables);
    assert_eq!(single_support_outliers, two_support_outliers);
    assert_eq!(single_gc_percent_widths, two_gc_percent_widths);

    Ok(())
}

// REWRITE-PUBLIC-TEST: ref-gc-bias artifact behavior currently uses a private package loader.
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
        ref_genome: cfdnalab::run_like_cli::common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: output_dir.to_path_buf(),
        output_prefix: String::new(),
        n_threads,
        n_positions: 100,
        seed: Some(29),
        windows: cfdnalab::run_like_cli::ref_gc_bias::RefGCWindowsArgs {
            by_bed: Some(bed_path.clone()),
        },
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: Some(vec![blacklist_path.clone()]),
        fragment_lengths: cfdnalab::run_like_cli::common::FragmentLengthArgs {
            min_fragment_length: 10,
            max_fragment_length: 12,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 80,
        logging: LoggingArgs::default(),
    };

    // Act
    run(&make_cfg(single_thread_out.path(), 1))?;
    run(&make_cfg(two_thread_out.path(), 2))?;

    // Assert
    let single_package = single_thread_out.path().join("ref_gc_package.zarr");
    let two_package = two_thread_out.path().join("ref_gc_package.zarr");
    let (
        single_counts,
        single_support_unobservables,
        single_support_outliers,
        single_gc_percent_widths,
    ) = load_ref_gc_package_arrays(&single_package)?;
    let (two_counts, two_support_unobservables, two_support_outliers, two_gc_percent_widths) =
        load_ref_gc_package_arrays(&two_package)?;

    assert_eq!(single_counts, two_counts);
    assert_eq!(single_support_unobservables, two_support_unobservables);
    assert_eq!(single_support_outliers, two_support_outliers);
    assert_eq!(single_gc_percent_widths, two_gc_percent_widths);

    Ok(())
}

// REWRITE-PUBLIC-TEST: ref-gc-bias artifact behavior currently uses a private package loader.
#[test]
fn fixed_seed_ref_gc_bias_is_deterministic_for_same_tile_size() -> Result<()> {
    let reference = fixtures::simple_reference_twobit()?;
    let first_out = TempDir::new()?;
    let second_out = TempDir::new()?;

    let make_cfg = |output_dir: &Path| RefGCBiasConfig {
        ref_genome: cfdnalab::run_like_cli::common::Ref2BitRequiredArgs {
            ref_2bit: reference.path.clone(),
        },
        output_dir: output_dir.to_path_buf(),
        output_prefix: String::new(),
        n_threads: 1,
        n_positions: 100,
        seed: Some(41),
        windows: Default::default(),
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: None,
        fragment_lengths: cfdnalab::run_like_cli::common::FragmentLengthArgs {
            min_fragment_length: 10,
            max_fragment_length: 12,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 80,
        logging: LoggingArgs::default(),
    };

    // Manual expectations:
    // - `ref-gc-bias` derives per-tile seeds from the top-level seed before any work starts.
    // - Re-running the command with the same seed and the same tile layout must therefore produce
    //   exactly the same sampled starts and the same reference package arrays.
    // - We keep `tile_size = 80` so the run spans multiple tiles and exercises the real per-tile
    //   seeded sampling path rather than a single-tile degenerate case.
    run(&make_cfg(first_out.path()))?;
    run(&make_cfg(second_out.path()))?;

    let first_package = first_out.path().join("ref_gc_package.zarr");
    let second_package = second_out.path().join("ref_gc_package.zarr");
    let (
        first_counts,
        first_support_unobservables,
        first_support_outliers,
        first_gc_percent_widths,
    ) = load_ref_gc_package_arrays(&first_package)?;
    let (
        second_counts,
        second_support_unobservables,
        second_support_outliers,
        second_gc_percent_widths,
    ) = load_ref_gc_package_arrays(&second_package)?;

    assert_eq!(first_counts, second_counts);
    assert_eq!(first_support_unobservables, second_support_unobservables);
    assert_eq!(first_support_outliers, second_support_outliers);
    assert_eq!(first_gc_percent_widths, second_gc_percent_widths);

    Ok(())
}
