use super::*;

#[test]
fn final_bin_includes_configured_maximum_length() -> Result<()> {
    let edges = build_length_bin_edges(100, 220, 3)?;
    assert_eq!(edges.first(), Some(&100.0));
    assert_eq!(edges.last(), Some(&220.0));
    let mut statistics = OverlappingLengthStatistics::new(edges)?;
    statistics.add_position(220.0, 2.0, 2)?;
    assert_eq!(statistics.length_bin_base_counts.last(), Some(&1));
    Ok(())
}

#[test]
fn merging_tiles_is_additive() -> Result<()> {
    let edges = build_length_bin_edges(100, 220, 3)?;
    let mut left = OverlappingLengthStatistics::new(edges.clone())?;
    let mut right = OverlappingLengthStatistics::new(edges)?;
    left.add_position(101.0, 2.0, 1)?;
    right.add_position(101.0, 4.0, 2)?;
    left.merge(right)?;

    assert_eq!(left.eligible_covered_bases, 2);
    assert_eq!(left.length_bin_base_counts[0], 2);
    assert_eq!(left.observed_signal_sums[0], 6.0);
    assert_eq!(left.raw_depth_frequencies.get(&1), Some(&1));
    assert_eq!(left.raw_depth_frequencies.get(&2), Some(&1));
    assert_eq!(
        left.length_bin_depth_statistics[0]
            .get(&1)
            .map(|statistics| statistics.position_count),
        Some(1)
    );
    assert_eq!(
        left.length_bin_depth_statistics[0]
            .get(&2)
            .map(|statistics| statistics.position_count),
        Some(1)
    );
    assert!((left.mean_average_overlapping_length()? - 101.0).abs() < 1.0e-12);
    Ok(())
}

#[test]
fn binwise_first_normalization_needs_no_second_positional_sweep() {
    let positional_signal = [2.0, 4.0, 8.0];
    let bin_division_factor = 2.5;
    let explicit_normalized_mean = positional_signal
        .iter()
        .map(|value| value / bin_division_factor)
        .sum::<f64>()
        / positional_signal.len() as f64;
    let sufficient_statistic_mean =
        positional_signal.iter().sum::<f64>() / positional_signal.len() as f64
            / bin_division_factor;
    assert!((explicit_normalized_mean - sufficient_statistic_mean).abs() < 1.0e-12);
}
