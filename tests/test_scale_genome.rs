// Place this in the same module as `apply_scaling_in_place`

#[cfg(test)]
mod tests {
    use cfdnalab::utils::coverage::scale_genome::apply_scaling_in_place;

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
        // Expected:
        //   indices 0..5  divide by 2.0
        //   indices 5..15 set to 0.0
        //   indices 15..20 divide by 0.5
        let core_start = 105u32;
        let bins = vec![(90u64, 110u64, 2.0f32), (110, 120, 0.0), (120, 140, 0.5)];
        let mut cov = vec![10.0f32; 20];
        apply_scaling_in_place(&mut cov, core_start, &bins);

        let mut expected = vec![0.0f32; 20];
        // 0..5: 10 / 2.0 = 5
        for v in &mut expected[0..5] {
            *v = 5.0;
        }
        // 5..15: 0 by zero-scaled bin
        // 15..20: 10 / 0.5 = 20
        for v in &mut expected[15..20] {
            *v = 20.0;
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
        apply_scaling_in_place(&mut cov, core_start, &bins);

        let expected = vec![2.0f32; 10]; // 8 / 4.0
        assert_slice_eq_eps(&cov, &expected, 1e-6);
    }

    #[test]
    fn scaling_zero_bin_covers_entire_tile() {
        // Entire tile lies in a zero-scaled bin
        let core_start = 60u32;
        let bins = vec![(50u64, 150u64, 0.0f32)];
        let mut cov = (0..20).map(|i| i as f32 + 1.0).collect::<Vec<_>>();
        apply_scaling_in_place(&mut cov, core_start, &bins);

        let expected = vec![0.0f32; 20];
        assert_slice_eq_eps(&cov, &expected, 1e-6);
    }

    #[test]
    fn scaling_noop_on_empty_inputs() {
        // No bins or empty coverage should be a no-op
        let core_start = 100u32;

        let mut cov1 = vec![] as Vec<f32>;
        apply_scaling_in_place(&mut cov1, core_start, &[]);
        assert!(cov1.is_empty());

        let mut cov2 = vec![1.0f32, 2.0, 3.0];
        apply_scaling_in_place(&mut cov2, core_start, &[]);
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
        let mut cov = vec![6.0f32; 30]; // We will only touch first 30 entries across the three bins
        apply_scaling_in_place(&mut cov, core_start, &bins);

        let mut expected = vec![0.0f32; 30];
        for v in &mut expected[0..10] {
            *v = 3.0; // 6 / 2
        }
        for v in &mut expected[10..20] {
            *v = 0.0; // zero bin
        }
        for v in &mut expected[20..30] {
            *v = 12.0; // 6 / 0.5
        }
        assert_slice_eq_eps(&cov, &expected, 1e-6);
    }
}
