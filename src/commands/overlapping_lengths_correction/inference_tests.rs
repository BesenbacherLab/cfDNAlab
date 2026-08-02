use super::*;
use crate::commands::overlapping_lengths_correction::package::{
    MixtureFitParameters, OVERLAPPING_LENGTHS_CORRECTION_SCHEMA_VERSION,
};
use crate::shared::gc_tag::GCTagValue;

fn fragment(start: u32, end: u32) -> Result<FragmentWithSegments> {
    Ok(FragmentWithSegments {
        tid: 0,
        interval: Interval::new(start, end)?,
        segments: None,
        gc_tag: GCTagValue::default(),
    })
}

/// Build a structurally complete package with a hand-calculated lookup.
fn package(
    length_bin_edges: &[f64],
    combined_weights: &[f64],
    maximum_fragment_length: u32,
) -> Arc<OverlappingLengthsCorrectionPackage> {
    assert_eq!(length_bin_edges.len(), combined_weights.len() + 1);
    let length_bin_midpoints = length_bin_edges
        .windows(2)
        .map(|pair| (pair[0] + pair[1]) / 2.0)
        .collect::<Vec<_>>();
    let bin_count = combined_weights.len();
    Arc::new(OverlappingLengthsCorrectionPackage {
        version: OVERLAPPING_LENGTHS_CORRECTION_SCHEMA_VERSION,
        length_bin_edges: length_bin_edges.to_vec(),
        length_bin_midpoints,
        length_bin_base_counts: vec![1; bin_count],
        observed_signal_sums: vec![1.0; bin_count],
        raw_depths: vec![1, 2],
        raw_depth_frequencies: vec![1, 1],
        observed_bias: vec![1.0; bin_count],
        first_fitted_bias: vec![1.0; bin_count],
        first_corrected_bias: vec![1.0; bin_count],
        second_fitted_bias: vec![1.0; bin_count],
        target_bias: vec![1.0; bin_count],
        noise_division_factors: vec![1.0; bin_count],
        skew_division_factors: vec![1.0; bin_count],
        mean_shift_division_factors: vec![1.0; bin_count],
        combined_weights: combined_weights.to_vec(),
        initial_fit: MixtureFitParameters {
            scale_multiplier: 1.0,
            skewness: 0.0,
            mean_fragment_length: 110.0,
        },
        refit: MixtureFitParameters {
            scale_multiplier: 1.0,
            skewness: 0.0,
            mean_fragment_length: 110.0,
        },
        minimum_fragment_length: 100,
        maximum_fragment_length,
        minimum_mapq: 30,
        require_proper_pair: false,
        reads_are_fragments: false,
        ignore_gap: false,
        blacklist_used: false,
    })
}

#[test]
fn overlapping_fragments_produce_expected_positional_bins() -> Result<()> {
    // Arrange
    let package = package(&[100.0, 110.0, 121.0], &[2.0, 3.0], 120);
    let mut collector = PositionalOverlapLengthCollector::new(
        package,
        Interval::new(0, 300)?,
        Interval::new(100, 280)?,
    )?;

    // Act
    collector.observe(&fragment(100, 200)?)?;
    collector.observe(&fragment(150, 270)?)?;
    let positional_bins = collector.finish()?;

    // Assert
    // [100, 150) has one 100 bp fragment. [150, 200) averages 100 and 120 to 110 bp.
    // [200, 270) has only the 120 bp fragment. [270, 280) is uncovered.
    let mut expected = vec![0_u32; 50];
    expected.extend(vec![1_u32; 120]);
    expected.extend(vec![NO_OVERLAP_LENGTH_BIN; 10]);
    assert_eq!(positional_bins.bin_indices(), expected);
    Ok(())
}

#[test]
fn positional_bins_multiply_complete_coverage_at_each_position() -> Result<()> {
    // Arrange
    let package = package(&[100.0, 110.0, 121.0], &[2.0, 3.0], 120);
    let mut collector = PositionalOverlapLengthCollector::new(
        package,
        Interval::new(0, 300)?,
        Interval::new(100, 280)?,
    )?;
    collector.observe(&fragment(100, 200)?)?;
    collector.observe(&fragment(150, 270)?)?;
    let positional_bins = collector.finish()?;
    let mut coverage = vec![4.0_f32; 180];

    // Act
    positional_bins.apply_to_coverage(&mut coverage)?;

    // Assert
    assert_eq!(&coverage[..50], vec![8.0_f32; 50]);
    assert_eq!(&coverage[50..170], vec![12.0_f32; 120]);
    assert_eq!(&coverage[170..], vec![4.0_f32; 10]);
    Ok(())
}

#[test]
fn segmented_fragment_updates_only_its_counted_reference_segments() -> Result<()> {
    // Arrange
    let package = package(&[100.0, 130.0], &[2.0], 120);
    let mut collector = PositionalOverlapLengthCollector::new(
        package,
        Interval::new(0, 300)?,
        Interval::new(100, 220)?,
    )?;
    let mut segmented_fragment = fragment(100, 220)?;
    segmented_fragment.segments = Some(
        [Interval::new(100, 140)?, Interval::new(180, 220)?]
            .into_iter()
            .collect(),
    );

    // Act
    collector.observe(&segmented_fragment)?;
    let positional_bins = collector.finish()?;

    // Assert
    let mut expected = vec![0_u32; 40];
    expected.extend(vec![NO_OVERLAP_LENGTH_BIN; 40]);
    expected.extend(vec![0_u32; 40]);
    assert_eq!(positional_bins.bin_indices(), expected);
    Ok(())
}
