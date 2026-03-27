use super::{
    InsertionAnchor, partition_deletion_by_aligned_overlap, partition_insertion_by_aligned_overlap,
};
use crate::shared::interval::Interval;
use fxhash::FxHashMap;

#[test]
fn partition_deletion_helper_splits_left_overlap_and_right_parts() {
    // Human verification status: unverified
    // Fragment [100,220), aligned overlap [160,180), deletion [150,190):
    // - left non-overlap [150,160) => 10 bp
    // - overlap piece [160,180)
    // - right non-overlap [180,190) => 10 bp
    let fragment_interval = Interval::new(100_u32, 220_u32).expect("test fragment interval");
    let aligned_overlap_interval =
        Some(Interval::new(160_u32, 180_u32).expect("test overlap interval"));
    let deletion_interval = Interval::new(150_u32, 190_u32).expect("test deletion interval");

    let mut nonoverlap_bases_bp = 0_u32;
    let mut overlap_deletion_intervals = Vec::new();

    partition_deletion_by_aligned_overlap(
        deletion_interval,
        fragment_interval,
        aligned_overlap_interval,
        &mut nonoverlap_bases_bp,
        &mut overlap_deletion_intervals,
    );

    assert_eq!(nonoverlap_bases_bp, 20);
    assert_eq!(
        overlap_deletion_intervals,
        vec![Interval::new(160_u32, 180_u32).expect("expected overlap interval")]
    );
}

#[test]
fn partition_insertion_helper_splits_nonoverlap_and_keeps_overlap_max() {
    // Human verification status: unverified
    // Fragment [100,220), aligned overlap [160,180):
    // - insertion at 120 is inside the fragment but outside the overlap => non-overlap +3
    // - insertions at 170 inside the overlap keep the maximum per read anchor => 5
    let fragment_interval = Interval::new(100_u32, 220_u32).expect("test fragment interval");
    let aligned_overlap_interval =
        Some(Interval::new(160_u32, 180_u32).expect("test overlap interval"));

    let mut nonoverlap_bases_bp = 0_u32;
    let mut overlap_insertions_by_anchor: FxHashMap<u32, u32> = FxHashMap::default();

    partition_insertion_by_aligned_overlap(
        InsertionAnchor {
            reference_position: 120,
            inserted_length: 3,
        },
        fragment_interval,
        aligned_overlap_interval,
        &mut nonoverlap_bases_bp,
        &mut overlap_insertions_by_anchor,
    );
    partition_insertion_by_aligned_overlap(
        InsertionAnchor {
            reference_position: 170,
            inserted_length: 4,
        },
        fragment_interval,
        aligned_overlap_interval,
        &mut nonoverlap_bases_bp,
        &mut overlap_insertions_by_anchor,
    );
    partition_insertion_by_aligned_overlap(
        InsertionAnchor {
            reference_position: 170,
            inserted_length: 5,
        },
        fragment_interval,
        aligned_overlap_interval,
        &mut nonoverlap_bases_bp,
        &mut overlap_insertions_by_anchor,
    );

    assert_eq!(nonoverlap_bases_bp, 3);
    assert_eq!(overlap_insertions_by_anchor.get(&170), Some(&5));
    assert_eq!(overlap_insertions_by_anchor.len(), 1);
}
