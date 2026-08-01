use super::*;

fn package() -> OverlappingLengthsCorrectionPackage {
    OverlappingLengthsCorrectionPackage {
        version: 1,
        length_bin_edges: vec![100.0, 103.0, 106.0],
        length_bin_midpoints: vec![101.5, 104.5],
        length_bin_base_counts: vec![1, 1],
        observed_signal_sums: vec![1.0, 1.0],
        raw_depths: vec![1],
        raw_depth_frequencies: vec![2],
        observed_bias: vec![1.0; 2],
        first_fitted_bias: vec![1.0; 2],
        first_corrected_bias: vec![1.0; 2],
        second_fitted_bias: vec![1.0; 2],
        target_bias: vec![1.0; 2],
        noise_division_factors: vec![1.0; 2],
        skew_division_factors: vec![1.0; 2],
        mean_shift_division_factors: vec![1.0; 2],
        combined_weights: vec![0.5, 2.0],
        initial_fit: MixtureFitParameters {
            scale_multiplier: 1.0,
            skewness: 0.0,
            mean_fragment_length: 103.0,
        },
        refit: MixtureFitParameters {
            scale_multiplier: 1.0,
            skewness: 0.0,
            mean_fragment_length: 103.0,
        },
        minimum_fragment_length: 100,
        maximum_fragment_length: 105,
        minimum_mapq: 30,
        require_proper_pair: false,
        reads_are_fragments: false,
        ignore_gap: false,
        gc_mode: "none".to_string(),
        scaling_enabled: false,
        blacklist_used: false,
    }
}

#[test]
fn lookup_clips_values_to_extreme_bins() -> Result<()> {
    let package = package();

    assert_eq!(package.weight_for_average_length(80.0)?, 0.5);
    assert_eq!(package.weight_for_average_length(220.0)?, 2.0);

    Ok(())
}

#[test]
fn zarr_round_trip_preserves_package_values() -> Result<()> {
    let package = package();
    let temporary_directory = tempfile::tempdir()?;
    let package_path = temporary_directory.path().join("model.zarr");

    package.write_zarr(&package_path)?;

    assert_eq!(
        OverlappingLengthsCorrectionPackage::from_file(&package_path)?,
        package
    );
    Ok(())
}
