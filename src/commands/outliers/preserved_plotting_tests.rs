use super::*;

#[test]
fn coverage_bins_preserve_zeros_exclude_masks_and_ignore_chunk_boundaries() {
    let mut bins = CoverageBins::new(Interval::new(100, 120).unwrap(), 4).unwrap();
    bins.add_coverage(
        Interval::new(100, 105).unwrap(),
        &[0.0, 2.0, 8.0, 4.0, 10.0],
        Some(&[0, 0, 1, 0, 0]),
    )
    .unwrap();
    bins.add_coverage(
        Interval::new(105, 120).unwrap(),
        &[
            6.0, 0.0, 0.0, 40.0, 0.0, 2.0, 6.0, 9.0, 9.0, 9.0, 9.0, 0.0, 0.0, 0.0, 0.0,
        ],
        Some(&[0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0]),
    )
    .unwrap();

    // Eligible values are [0,2,4], [10,6], [0,2,6], [], [0,0,0,0]. The chunk split
    // falls inside the second bin. Masked high values must affect neither means nor peaks
    assert_eq!(
        bins.bins
            .iter()
            .map(|bin| bin.eligible_positions)
            .collect::<Vec<_>>(),
        vec![3, 2, 3, 0, 4]
    );
    assert_eq!(bins.bins[0].mean(), Some(2.0));
    assert_eq!(bins.bins[1].mean(), Some(8.0));
    assert!((bins.bins[2].mean().unwrap() - 8.0 / 3.0).abs() < 1e-12);
    assert_eq!(bins.bins[3].mean(), None);
    assert_eq!(bins.bins[4].mean(), Some(0.0));
    // Keep peaks at their observed eligible positions, including beside blacklisted gaps
    assert_eq!(
        bins.bins[..3]
            .iter()
            .map(|bin| bin.maximum_position)
            .collect::<Vec<_>>(),
        vec![103, 104, 111]
    );
    assert_eq!(
        bins.bins.iter().map(|bin| bin.maximum).collect::<Vec<_>>(),
        vec![4.0, 10.0, 6.0, 0.0, 0.0]
    );
    assert_eq!(
        bins.bins.iter().map(|bin| bin.coverage_sum).sum::<f64>(),
        30.0
    );
}

#[test]
fn bins_cover_partial_last_bins_and_bound_the_number_of_plot_points() {
    let bins = CoverageBins::new(Interval::new(101, 112).unwrap(), 4).unwrap();
    assert_eq!(
        bins.bins
            .iter()
            .map(|bin| bin.interval.as_tuple())
            .collect::<Vec<_>>(),
        vec![(101, 105), (105, 109), (109, 112)]
    );
    for length in [1, 2_000, 2_001, 4_001, u32::MAX] {
        let bins = CoverageBins::for_plot(Interval::new(0, length).unwrap()).unwrap();
        assert!(bins.bins.len() <= MAXIMUM_PLOT_BINS as usize);
        assert_eq!(
            bins.bins
                .iter()
                .map(|bin| u64::from(bin.interval.len()))
                .sum::<u64>(),
            u64::from(length)
        );
    }
}

#[test]
fn clipped_plot_chunks_keep_absolute_coverage_and_mask_offsets() {
    let mut bins = CoverageBins::new(Interval::new(102, 108).unwrap(), 2).unwrap();

    bins.add_coverage(
        Interval::new(100, 105).unwrap(),
        &[100.0, 200.0, 2.0, 4.0, 6.0],
        Some(&[0, 0, 0, 1, 0]),
    )
    .unwrap();
    bins.add_coverage(
        Interval::new(105, 110).unwrap(),
        &[8.0, 10.0, 12.0, 300.0, 400.0],
        Some(&[0, 0, 1, 0, 0]),
    )
    .unwrap();
    // Touching chunks contain no positions in the half-open plot interval
    bins.add_coverage(Interval::new(100, 102).unwrap(), &[900.0, 900.0], None)
        .unwrap();
    bins.add_coverage(Interval::new(108, 110).unwrap(), &[900.0, 900.0], None)
        .unwrap();

    // The bins contain [2], [6,8], [10]. Mask offsets refer to input chunks, not plot starts
    let expected = [
        (1, 2.0, 2.0, 102),
        (2, 14.0, 8.0, 105),
        (1, 10.0, 10.0, 106),
    ];
    for (bin, (positions, sum, maximum, peak_position)) in bins.bins.iter().zip(expected) {
        assert_eq!(bin.eligible_positions, positions);
        assert_eq!(bin.coverage_sum, sum);
        assert_eq!(bin.mean(), Some(sum / positions as f64));
        assert_eq!(bin.maximum, maximum);
        assert_eq!(bin.maximum_position, peak_position);
    }
}

#[test]
fn tied_plot_peaks_keep_the_first_eligible_position_across_chunks() {
    let mut bins = CoverageBins::new(Interval::new(100, 108).unwrap(), 4).unwrap();

    bins.add_coverage(
        Interval::new(100, 102).unwrap(),
        &[99.0, 5.0],
        Some(&[1, 0]),
    )
    .unwrap();
    bins.add_coverage(
        Interval::new(102, 108).unwrap(),
        &[5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        Some(&[0, 0, 1, 0, 0, 0]),
    )
    .unwrap();

    assert_eq!(bins.bins[0].maximum_position, 101);
    assert_eq!(bins.bins[0].maximum, 5.0);
    assert_eq!(bins.bins[0].eligible_positions, 3);
    assert_eq!(bins.bins[0].coverage_sum, 10.0);
    // An all-zero bin's peak must also point at an observed eligible base
    assert_eq!(bins.bins[1].maximum_position, 105);
    assert_eq!(bins.bins[1].mean(), Some(0.0));
    assert_eq!(bins.bins[1].eligible_positions, 3);
}
