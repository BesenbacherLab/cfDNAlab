use super::*;
use crate::{
    commands::gc_bias::{
        binning::{BinnedAxis, bins_from_edges, compute_bin_edges},
        counting::{GCCounts, build_gc_prefixes},
        outliers::{
            OutlierAction, OutlierRule, OutlierScope, OutlierStats, apply_outliers_to_matrix,
            interpolated_quantile, outlier_bounds,
        },
        support_masking::build_extreme_bins_support_mask,
    },
    shared::interval::Interval,
};
use anyhow::Result;
use fxhash::FxHashMap;
use ndarray::array;

#[test]
fn get_fragment_gc_uses_sequence_interval_as_prefix_origin() -> Result<()> {
    // Manual derivation:
    // - Prefixes are built from the loaded reference slice [900,961), not from chromosome
    //   origin 0.
    // - The sequence slice is 61 C bases, so fragment [900,961) has 61 GC bases.
    // - A local-origin bug would either ask the 61 bp prefix for [900,961) or otherwise fail
    //   to count the loaded slice as the fragment interval.
    let prefixes = build_gc_prefixes(&vec![b'C'; 61]);
    let fragment_interval = Interval::new(900_u64, 961_u64)?;
    let sequence_interval = Interval::new(900_u64, 961_u64)?;

    let gc_count = get_fragment_gc(fragment_interval, sequence_interval, 0, &prefixes, 0.0)?;

    assert_eq!(gc_count, Some(61));
    Ok(())
}

#[test]
fn get_fragment_gc_returns_none_when_fragment_is_outside_loaded_sequence() -> Result<()> {
    // Manual derivation:
    // - Prefixes cover only [900,961).
    // - Fragment [961,1022) is a valid reference interval, but its contracted GC window is
    //   completely outside the loaded sequence, so this is a legitimate missing correction
    //   rather than an indexing error.
    let prefixes = build_gc_prefixes(&vec![b'C'; 61]);
    let fragment_interval = Interval::new(961_u64, 1022_u64)?;
    let sequence_interval = Interval::new(900_u64, 961_u64)?;

    let gc_count = get_fragment_gc(fragment_interval, sequence_interval, 0, &prefixes, 0.0)?;

    assert_eq!(gc_count, None);
    Ok(())
}

#[test]
fn masks_extreme_gc_bins_per_side_in_square_matrix() {
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
fn apply_outliers_per_length_winsorizes_rows() {
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

#[test]
fn should_use_effective_length_when_binning_to_gc_percent_with_end_offset() {
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
    // Arrange: length=6, end_offset=2 -> effective length is 2bp, so gc>2 is unreachable.
    let mut counts = GCCounts::new(6, 6, 2, (0, 0)).expect("counts init");
    counts.set(6, 2, 10.0);

    // Act: smooth only the reachable portion of the row.
    counts
        .smooth_length_rows_in_place(1.0, 1)
        .expect("smoothing should succeed for valid sigma and radius");

    // Assert: unreachable GC counts are absent and storage matches the effective length.
    assert!(counts.get(6, 3).is_none());
    assert_eq!(counts.borrow_raw_counts().len(), 3);
}

#[test]
fn should_place_gc_counts_in_matching_percent_bins() {
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
                "bin {} expected 2.0, got {}",
                idx,
                val
            ),
            67 => assert!(
                (val - 3.0).abs() < 1e-12,
                "bin {} expected 3.0, got {}",
                idx,
                val
            ),
            _ => assert!(val.abs() < 1e-12, "bin {} expected 0, got {}", idx, val),
        }
    }
}

#[test]
fn should_propagate_acgt_totals_and_length_metadata() {
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

#[test]
fn reports_offsets_based_on_effective_length() {
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
    let counts = GCCounts::new(10, 12, 0, (0, 0)).expect("init counts");
    assert!(counts.length_bounds(9).is_err());
    assert!(counts.length_bounds(13).is_err());
}

#[test]
fn leaves_zero_rows_untouched_in_mean_scaling() {
    // Arrange: first length row has no mass; second has values that should be mean-scaled.
    let counts = array![[0.0, 0.0], [2.0, 4.0]];
    let mask = array![[true, true], [true, true]];

    // Act
    let scaled = mean_scale_per_length_array(&counts, 0.0, Some(&mask));

    // Assert: empty row stays zero; non-empty row divides by its mean (3.0).
    assert!(
        scaled.row(0).iter().all(|&value| value == 0.0),
        "zero row should remain zero after scaling"
    );
    assert!((scaled[(1, 0)] - 2.0 / 3.0).abs() < 1e-12);
    assert!((scaled[(1, 1)] - 4.0 / 3.0).abs() < 1e-12);
}
