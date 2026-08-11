use super::*;

fn region(start: u32, end: u32, observed_mass: f64, target_mass: f64) -> CalledRegion {
    CalledRegion {
        interval: Interval::new(start, end).expect("valid called region"),
        observed_mass,
        target_mass,
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

    let (regions_by_chromosome, _) = merge_calling_pass_results(
        &["chr1".to_string()],
        calling_pass_tile_results,
    )
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
fn blacklist_lookup_uses_half_open_boundaries() {
    let blacklist = Interval::from_tuples(&[(10_u64, 20_u64), (30, 40)])
        .expect("valid blacklist intervals");

    assert!(!position_is_blacklisted(&blacklist, 9));
    assert!(position_is_blacklisted(&blacklist, 10));
    assert!(position_is_blacklisted(&blacklist, 19));
    assert!(!position_is_blacklisted(&blacklist, 20));
    assert!(position_is_blacklisted(&blacklist, 30));
}
