#[cfg(test)]
mod tests_apply_scaling {
    use cfdnalab::utils::coverage::scale_genome::apply_scaling_to_coverage_in_place;

    // Assert two slices are approximately equal within eps
    fn assert_slice_eq_eps(a: &[f32], b: &[f32], eps: f32) {
        assert_eq!(
            a.len(),
            b.len(),
            "Length mismatch: {} vs {}",
            a.len(),
            b.len()
        );
        for (i, (x, y)) in a.iter().zip(b.iter()).enumerate() {
            let diff = (x - y).abs();
            assert!(
                diff <= eps,
                "Mismatch at {}: got {}, expected {}, |diff|={}",
                i,
                x,
                y,
                diff
            );
        }
    }

    #[test]
    fn scaling_mixed_bins_partial_overlap() {
        // Tile core [105, 125), len 20
        // Bins: [90,110) sf=2.0, [110,120) sf=0.0, [120,140) sf=0.5
        // Expected (multiply-by-sf):
        //   indices 0..5   *= 2.0  -> 20
        //   indices 5..15  *= 0.0  -> 0
        //   indices 15..20 *= 0.5  -> 5
        let core_start = 105u32;
        let bins = vec![(90u64, 110u64, 2.0f32), (110, 120, 0.0), (120, 140, 0.5)];
        let mut cov = vec![10.0f32; 20];
        apply_scaling_to_coverage_in_place(&mut cov, core_start, &bins);

        let mut expected = vec![0.0f32; 20];
        for v in &mut expected[0..5] {
            *v = 20.0;
        }
        // 5..15 stay 0.0
        for v in &mut expected[15..20] {
            *v = 5.0;
        }
        assert_slice_eq_eps(&cov, &expected, 1e-6);
    }

    #[test]
    fn scaling_starts_exactly_at_bin_boundary() {
        // Tile core starts exactly at previous bin end
        // Bins: [0,100) sf=2.0, [100,200) sf=4.0
        // Tile: core_start=100, len=10 so all in second bin
        let core_start = 100u32;
        let bins = vec![(0u64, 100u64, 2.0f32), (100, 200, 4.0)];
        let mut cov = vec![8.0f32; 10];
        apply_scaling_to_coverage_in_place(&mut cov, core_start, &bins);

        let expected = vec![32.0f32; 10]; // 8 * 4.0
        assert_slice_eq_eps(&cov, &expected, 1e-6);
    }

    #[test]
    fn scaling_zero_bin_covers_entire_tile() {
        // Entire tile lies in a zero-scaled bin
        let core_start = 60u32;
        let bins = vec![(50u64, 150u64, 0.0f32)];
        let mut cov = (0..20).map(|i| i as f32 + 1.0).collect::<Vec<_>>();
        apply_scaling_to_coverage_in_place(&mut cov, core_start, &bins);

        let expected = vec![0.0f32; 20];
        assert_slice_eq_eps(&cov, &expected, 1e-6);
    }

    #[test]
    fn scaling_noop_on_empty_inputs() {
        // No bins or empty coverage should be a no-op
        let core_start = 100u32;

        let mut cov1 = vec![] as Vec<f32>;
        apply_scaling_to_coverage_in_place(&mut cov1, core_start, &[]);
        assert!(cov1.is_empty());

        let mut cov2 = vec![1.0f32, 2.0, 3.0];
        apply_scaling_to_coverage_in_place(&mut cov2, core_start, &[]);
        assert_slice_eq_eps(&cov2, &[1.0, 2.0, 3.0], 1e-6);
    }

    #[test]
    fn scaling_multiple_bins_exact_edges() {
        // Tile aligns exactly to bin edges
        // Bins: [1000,1010) sf=2, [1010,1020) sf=0, [1020,1030) sf=0.5
        let core_start = 1000u32;
        let bins = vec![
            (1000u64, 1010u64, 2.0f32),
            (1010, 1020, 0.0),
            (1020, 1030, 0.5),
        ];
        let mut cov = vec![6.0f32; 30]; // First 30 cover those bins
        apply_scaling_to_coverage_in_place(&mut cov, core_start, &bins);

        let mut expected = vec![0.0f32; 30];
        for v in &mut expected[0..10] {
            *v = 12.0; // 6 * 2
        }
        for v in &mut expected[10..20] {
            *v = 0.0; // zero bin
        }
        for v in &mut expected[20..30] {
            *v = 3.0; // 6 * 0.5
        }
        assert_slice_eq_eps(&cov, &expected, 1e-6);
    }

    #[test]
    fn scaling_bins_entirely_left_of_tile_noop() {
        // All bins end before the tile starts -> no effect
        let core_start = 500u32;
        let bins = vec![(100u64, 200u64, 2.0f32), (200, 250, 0.5)];
        let mut cov = vec![7.0f32; 10];
        apply_scaling_to_coverage_in_place(&mut cov, core_start, &bins);

        assert_slice_eq_eps(&cov, &vec![7.0f32; 10], 1e-6);
    }

    #[test]
    fn scaling_bins_entirely_right_of_tile_noop() {
        // All bins start after the tile ends -> no effect
        let core_start = 100u32;
        let bins = vec![(1000u64, 2000u64, 2.0f32)];
        let mut cov = vec![4.0f32; 8];
        apply_scaling_to_coverage_in_place(&mut cov, core_start, &bins);

        assert_slice_eq_eps(&cov, &vec![4.0f32; 8], 1e-6);
    }

    #[test]
    fn scaling_partial_left_overlap_only() {
        // Bin overlaps only the left edge of the tile
        // Tile: [200, 210), Bin: [195, 205) sf=3.0 -> indices 0..5 scaled by 3
        let core_start = 200u32;
        let bins = vec![(195u64, 205u64, 3.0f32)];
        let mut cov = vec![2.0f32; 10];
        apply_scaling_to_coverage_in_place(&mut cov, core_start, &bins);

        let mut expected = vec![2.0f32; 10];
        for v in &mut expected[0..5] {
            *v = 6.0;
        }
        assert_slice_eq_eps(&cov, &expected, 1e-6);
    }

    #[test]
    fn scaling_partial_right_overlap_only() {
        // Bin overlaps only the right edge of the tile
        // Tile: [300, 312), Bin: [308, 400) sf=0.25 -> indices 8..12 scaled by 0.25
        let core_start = 300u32;
        let bins = vec![(308u64, 400u64, 0.25f32)];
        let mut cov = vec![20.0f32; 12];
        apply_scaling_to_coverage_in_place(&mut cov, core_start, &bins);

        let mut expected = vec![20.0f32; 12];
        for v in &mut expected[8..12] {
            *v = 5.0;
        }
        assert_slice_eq_eps(&cov, &expected, 1e-6);
    }

    #[test]
    fn scaling_many_small_bins_inside_tile() {
        // Tile: [1000, 1010). Bins: per-base alternating 2.0 and 0.5.
        let core_start = 1000u32;
        let bins = (1000u64..1010u64)
            .enumerate()
            .map(|(i, s)| {
                let sf = if i % 2 == 0 { 2.0f32 } else { 0.5f32 };
                (s, s + 1, sf)
            })
            .collect::<Vec<_>>();
        let mut cov = vec![10.0f32; 10];
        apply_scaling_to_coverage_in_place(&mut cov, core_start, &bins);

        let expected: Vec<f32> = (0..10)
            .map(|i| if i % 2 == 0 { 20.0 } else { 5.0 })
            .collect();
        assert_slice_eq_eps(&cov, &expected, 1e-6);
    }
}
