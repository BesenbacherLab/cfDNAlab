use super::*;

fn region(start: u32, end: u32, observed_mass: f64, target_mass: f64) -> CalledRegion {
    CalledRegion {
        interval: Interval::new(start, end).expect("valid called region"),
        observed_mass,
        target_mass,
        maximum_coverage: 20,
    }
}

#[test]
fn merges_touching_regions_across_tile_boundaries_and_preserves_mass() {
    let calling_pass_tile_results = vec![
        CallingPassTileResult {
            chromosome: "chr1".to_string(),
            regions: vec![region(10, 20, 100.0, 20.0)],
            counters: FCoverageCounters::default(),
        },
        CallingPassTileResult {
            chromosome: "chr1".to_string(),
            regions: vec![region(20, 25, 50.0, 5.0), region(30, 31, 20.0, 2.0)],
            counters: FCoverageCounters::default(),
        },
    ];

    let (regions_by_chromosome, _) =
        merge_calling_pass_results(&["chr1".to_string()], calling_pass_tile_results)
            .expect("merge tile regions");
    let regions = &regions_by_chromosome["chr1"];

    assert_eq!(regions.len(), 2);
    assert_eq!(regions[0].interval.as_tuple(), (10, 25));
    assert_eq!(regions[0].observed_mass, 150.0);
    assert_eq!(regions[0].target_mass, 25.0);
    assert!((regions[0].keep_weight() - 1.0 / 6.0).abs() < 1e-12);
    assert_eq!(regions[1].interval.as_tuple(), (30, 31));
}

#[test]
fn zero_target_produces_zero_keep_weight() {
    let called_region = region(10, 11, 20.0, 0.0);

    assert_eq!(called_region.keep_weight(), 0.0);
}

#[test]
fn merged_calls_keep_independent_peaks_and_never_join_across_chromosomes() {
    let mut first = region(10, 40, 210.0, 30.0);
    first.maximum_coverage = 7;
    let mut next = region(40, 70, 330.0, 60.0);
    next.maximum_coverage = 11;
    let mut separate = region(100, 130, 90.0, 30.0);
    separate.maximum_coverage = 3;
    let mut other_chromosome = region(130, 160, 60.0, 30.0);
    other_chromosome.maximum_coverage = 2;
    let results = vec![
        CallingPassTileResult {
            chromosome: "chr1".into(),
            regions: vec![first],
            counters: FCoverageCounters::default(),
        },
        CallingPassTileResult {
            chromosome: "chr1".into(),
            regions: vec![next, separate],
            counters: FCoverageCounters::default(),
        },
        CallingPassTileResult {
            chromosome: "chr2".into(),
            regions: vec![other_chromosome],
            counters: FCoverageCounters::default(),
        },
    ];
    let chromosomes = ["chr1".into(), "chr2".into(), "chr_empty".into()];

    let (merged, _) = merge_calling_pass_results(&chromosomes, results).unwrap();

    assert_eq!(merged["chr1"].len(), 2);
    assert_eq!(merged["chr1"][0].interval.as_tuple(), (10, 70));
    assert_eq!(merged["chr1"][0].observed_mass, 540.0);
    assert_eq!(merged["chr1"][0].target_mass, 90.0);
    assert_eq!(merged["chr1"][0].maximum_coverage, 11);
    assert_eq!(merged["chr1"][1].interval.as_tuple(), (100, 130));
    assert_eq!(merged["chr1"][1].maximum_coverage, 3);
    assert_eq!(merged["chr2"].len(), 1);
    assert_eq!(merged["chr2"][0].interval.as_tuple(), (130, 160));
    assert_eq!(merged["chr2"][0].maximum_coverage, 2);
    assert!(merged["chr_empty"].is_empty());
}

#[test]
fn overlapping_tile_calls_are_rejected_instead_of_double_counting_their_mass() {
    let results = vec![
        CallingPassTileResult {
            chromosome: "chr1".into(),
            regions: vec![region(10, 40, 300.0, 30.0)],
            counters: FCoverageCounters::default(),
        },
        CallingPassTileResult {
            chromosome: "chr1".into(),
            regions: vec![region(39, 70, 310.0, 31.0)],
            counters: FCoverageCounters::default(),
        },
    ];

    let error = merge_calling_pass_results(&["chr1".into()], results).unwrap_err();

    assert!(error.to_string().contains("overlap unexpectedly"));
    assert!(error.to_string().contains("chr1"));
}

#[test]
fn blacklist_lookup_uses_half_open_boundaries() {
    let blacklist =
        Interval::from_tuples(&[(10_u64, 20_u64), (30, 40)]).expect("valid blacklist intervals");

    assert!(!position_is_blacklisted(&blacklist, 9));
    assert!(position_is_blacklisted(&blacklist, 10));
    assert!(position_is_blacklisted(&blacklist, 19));
    assert!(!position_is_blacklisted(&blacklist, 20));
    assert!(position_is_blacklisted(&blacklist, 30));
}

#[cfg(all(feature = "testing", feature = "plotters"))]
#[test]
fn actual_plot_chunk_loop_preserves_a_boundary_peak_masks_and_close_up_resolution() {
    use crate::{
        commands::cli_common::{ChromosomeArgs, IOCArgs},
        testing::{PairedFragmentSpec, TempBamBuilder},
    };
    let mut builder = TempBamBuilder::new()
        .name("plot_chunks")
        .contig("chr1", 2_000_000)
        .use_record_indexed_read_names()
        .paired_fragment(PairedFragmentSpec::new(0, 100, 30, 10));
    for _ in 0..20 {
        builder = builder.paired_fragment(PairedFragmentSpec::new(0, 999_990, 30, 10));
    }
    let bam = builder.build().unwrap();
    let output = tempfile::tempdir().unwrap();
    let mut config = OutliersConfig::new(
        IOCArgs {
            bam: bam.bam_path().to_path_buf(),
            output_dir: output.path().to_path_buf(),
            n_threads: 1,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".into()]),
            chromosomes_file: None,
        },
    );
    config.set_tile_size(1_000_000);
    config.fragment_lengths_mut().min_fragment_length = 30;
    config.fragment_lengths_mut().max_fragment_length = 30;
    let group = PreservedRegion {
        region_id: "region_000001".into(),
        chromosome: "chr1".into(),
        interval: Interval::new(999_990, 1_000_020).unwrap(),
        window_span: Interval::new(999_990, 1_099_990).unwrap(),
        plot_interval: Interval::new(0, 2_000_000).unwrap(),
        call_indices: vec![0],
        called_bases: 30,
        maximum_window_percentage: 10.0,
    };
    let blacklist = Interval::from_tuples(&[(999_999_u64, 1_000_001)]).unwrap();
    let (overview, close_up) =
        collect_preserved_plot_coverage(&config, &group, 0, 2_000_000, &blacklist).unwrap();
    // Exact raw coverage is 1 at [100,130), 20 at [999990,1000020), zero elsewhere
    // The two central bases are masked. A 1 Mb chunk split cuts this same peak into 10 and 20 bp
    let mut values = vec![0.0; 2_000_000];
    values[100..130].fill(1.0);
    values[999_990..1_000_020].fill(20.0);
    let mut mask = vec![0_u8; 2_000_000];
    mask[999_999..1_000_001].fill(1);
    let mut expected_overview = CoverageBins::for_plot(group.plot_interval).unwrap();
    expected_overview
        .add_coverage(group.plot_interval, &values, Some(&mask))
        .unwrap();
    let mut expected_close_up =
        CoverageBins::for_plot(Interval::new(999_890, 1_000_120).unwrap()).unwrap();
    expected_close_up
        .add_coverage(group.plot_interval, &values, Some(&mask))
        .unwrap();
    assert_eq!(overview, expected_overview);
    assert_eq!(close_up, expected_close_up);
}

/// Supply raw coverage ten on [999980,1000040) with different models on either side of 1 Mb.
///
/// Models are assigned explicitly so these tests isolate calling from statistical estimation.
#[cfg(feature = "testing")]
fn calling_pass_fixture() -> (
    crate::testing::TempBam,
    OutliersConfig,
    Vec<Option<CoreModel>>,
) {
    use crate::{
        commands::{
            cli_common::{ChromosomeArgs, IOCArgs},
            outliers::{
                config::OutlierModel,
                distribution::{CoverageParameters, TwoStageCoverageModel, diagnose_fit},
                striding::ModelSource,
            },
        },
        testing::{PairedFragmentSpec, TempBamBuilder},
    };
    let mut builder = TempBamBuilder::new()
        .name("calling_model_boundary")
        .contig("chr1", 1_000_300)
        .use_record_indexed_read_names();
    for _ in 0..10 {
        builder = builder.paired_fragment(PairedFragmentSpec::new(0, 999_980, 60, 20));
    }
    let bam = builder.build().unwrap();
    let mut config = OutliersConfig::new(
        IOCArgs {
            bam: bam.bam_path().to_path_buf(),
            output_dir: "unused".into(),
            n_threads: 1,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".into()]),
            chromosomes_file: None,
        },
    );
    config.set_tile_size(1_000_000);
    config.set_stride(1_000_000);
    config.set_min_mapq(0);
    config.fragment_lengths_mut().min_fragment_length = 60;
    config.fragment_lengths_mut().max_fragment_length = 60;
    config.validate().unwrap();
    let models = [(0, 1_000_000, 2.0, 5), (1_000_000, 1_000_300, 4.0, 8)]
        .into_iter()
        .map(|(start, end, mean, threshold)| {
            let parameters = CoverageParameters {
                component_mean: mean,
                dispersion: 0.0,
                zero_inflation: 0.5,
            };
            let model = TwoStageCoverageModel {
                family: OutlierModel::Zip,
                initial: parameters,
                initial_threshold: threshold,
                underlying_fit: parameters,
                final_threshold: threshold,
                second_fit_retained_positions: 100,
                initial_optimizer: None,
                underlying_optimizer: None,
            };
            let core = Interval::new(start, end).unwrap();
            Some(CoreModel {
                core,
                context: core,
                eligible_positions: 100,
                fragment_support: 10,
                source: ModelSource::Local,
                model,
                // Diagnostics are unused during calling but must still be valid fixture data
                diagnostics: diagnose_fit(&[(0, 50), (1, 50)].into_iter().collect(), model)
                    .unwrap(),
            })
        })
        .collect();
    (bam, config, models)
}

#[cfg(feature = "testing")]
#[test]
fn cross_core_call_sums_each_models_target_and_is_invariant_to_tile_splitting() {
    let (_bam, mut config, models) = calling_pass_fixture();
    let whole = Interval::new(0, 1_000_300).unwrap();
    let tile = Tile::new("chr1".into(), 0, 0, whole, whole).unwrap();
    let chromosomes = ["chr1".to_string()];
    for (target, expected_mass, expected_weight) in [
        (OutlierTarget::LocalMean, 200.0, 1.0 / 3.0),
        (OutlierTarget::Threshold, 420.0, 0.7),
        (OutlierTarget::Zero, 0.0, 0.0),
    ] {
        config.set_target(target);
        let unsplit = process_calling_pass_tile(&config, &tile, &[], &models).unwrap();
        let mut pieces = Vec::new();
        for (index, core) in [
            (0, Interval::new(0, 1_000_000).unwrap()),
            (1, Interval::new(1_000_000, 1_000_300).unwrap()),
        ] {
            let tile = Tile::new("chr1".into(), 0, index, core, whole).unwrap();
            pieces.push(process_calling_pass_tile(&config, &tile, &[], &models).unwrap());
        }
        let (split, _) = merge_calling_pass_results(&chromosomes, pieces).unwrap();

        // Ten fragments span 60 bp. Targets are 20*2+40*4, 20*5+40*8, or zero
        // Structural-zero inflation must not halve the local-mean target
        assert_eq!(unsplit.regions.len(), 1);
        assert_eq!(split["chr1"].len(), 1);
        for call in [&unsplit.regions[0], &split["chr1"][0]] {
            assert_eq!(call.interval.as_tuple(), (999_980, 1_000_040));
            assert_eq!(call.observed_mass, 600.0);
            assert_eq!(call.target_mass, expected_mass);
            assert_eq!(call.maximum_coverage, 10);
            assert!((call.keep_weight() - expected_weight).abs() < 1e-12);
        }
    }
}

#[cfg(feature = "testing")]
#[test]
fn masked_model_boundary_splits_calls_and_resets_their_masses() {
    let (_bam, mut config, models) = calling_pass_fixture();
    config.set_target(OutlierTarget::LocalMean);
    let whole = Interval::new(0, 1_000_300).unwrap();
    let tile = Tile::new("chr1".into(), 0, 0, whole, whole).unwrap();
    let blacklist = Interval::from_tuples(&[(999_999_u64, 1_000_001)]).unwrap();

    let result = process_calling_pass_tile(&config, &tile, &blacklist, &models).unwrap();

    // Removing the last left-core base and first right-core base leaves 19 and 39 called bp
    assert_eq!(result.regions.len(), 2);
    for (call, interval, observed, target, weight) in [
        (&result.regions[0], (999_980, 999_999), 190.0, 38.0, 0.2),
        (
            &result.regions[1],
            (1_000_001, 1_000_040),
            390.0,
            156.0,
            0.4,
        ),
    ] {
        assert_eq!(call.interval.as_tuple(), interval);
        assert_eq!(call.observed_mass, observed);
        assert_eq!(call.target_mass, target);
        assert_eq!(call.maximum_coverage, 10);
        assert!((call.keep_weight() - weight).abs() < 1e-12);
    }
}

#[cfg(feature = "testing")]
#[test]
fn fully_blacklisted_core_needs_no_model_and_preserves_neighboring_call_coordinates() {
    let (_bam, mut config, models) = calling_pass_fixture();
    config.set_target(OutlierTarget::LocalMean);
    let whole = Interval::new(0, 1_000_300).unwrap();
    let tile = Tile::new("chr1".into(), 0, 0, whole, whole).unwrap();

    for (masked_index, mask, expected_interval, observed_mass, target_mass) in [
        (0, (0_u64, 1_000_000), (1_000_000, 1_000_040), 400.0, 160.0),
        (1, (1_000_000, 1_000_300), (999_980, 1_000_000), 200.0, 40.0),
    ] {
        let mut model_slots = models.clone();
        model_slots[masked_index] = None;
        let blacklist = Interval::from_tuples(&[mask]).unwrap();

        let result = process_calling_pass_tile(&config, &tile, &blacklist, &model_slots).unwrap();

        // Ten 60 bp fragments cross the boundary. Masking leaves 40 right or 20 left bases
        assert_eq!(result.regions.len(), 1);
        let call = &result.regions[0];
        assert_eq!(call.interval.as_tuple(), expected_interval);
        assert_eq!(call.observed_mass, observed_mass);
        assert_eq!(call.target_mass, target_mass);
    }
}

#[cfg(feature = "testing")]
#[test]
fn eligible_zero_coverage_position_still_requires_a_model() {
    let (_bam, config, mut models) = calling_pass_fixture();
    models[0] = None;
    let whole = Interval::new(0, 1_000_300).unwrap();
    let tile = Tile::new("chr1".into(), 0, 0, whole, whole).unwrap();

    let error = process_calling_pass_tile(&config, &tile, &[], &models).unwrap_err();

    // Position zero has no fragments but remains eligible, so a missing model is an error
    assert!(error.to_string().contains("missing core model 0"));
    assert!(error.to_string().contains("eligible position 0"));
}

#[cfg(feature = "testing")]
#[test]
fn model_boundary_uses_the_new_threshold_inclusively_without_a_coverage_change() {
    let (_bam, mut config, mut models) = calling_pass_fixture();
    config.set_target(OutlierTarget::LocalMean);
    let whole = Interval::new(0, 1_000_300).unwrap();
    let tile = Tile::new("chr1".into(), 0, 0, whole, whole).unwrap();

    models[1].as_mut().unwrap().model.final_threshold = 10;
    let equal = process_calling_pass_tile(&config, &tile, &[], &models).unwrap();
    models[1].as_mut().unwrap().model.final_threshold = 11;
    let below = process_calling_pass_tile(&config, &tile, &[], &models).unwrap();

    assert_eq!(equal.regions.len(), 1);
    assert_eq!(equal.regions[0].interval.as_tuple(), (999_980, 1_000_040));
    assert_eq!(equal.regions[0].observed_mass, 600.0);
    assert_eq!(equal.regions[0].target_mass, 200.0);
    assert_eq!(below.regions.len(), 1);
    assert_eq!(below.regions[0].interval.as_tuple(), (999_980, 1_000_000));
    assert_eq!(below.regions[0].observed_mass, 200.0);
    assert_eq!(below.regions[0].target_mass, 40.0);
    assert_eq!(below.regions[0].maximum_coverage, 10);
}
