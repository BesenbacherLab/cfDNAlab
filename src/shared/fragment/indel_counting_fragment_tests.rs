use super::{
    FragmentWithIndelCounts, InsertionAnchor, partition_deletion_by_aligned_overlap,
    partition_insertion_by_aligned_overlap,
};
use crate::shared::clip_mode::ClipMode;
use crate::shared::indel_mode::IndelMode;
use crate::shared::interval::Interval;
use fxhash::FxHashMap;

#[test]
fn partition_deletion_helper_splits_left_overlap_and_right_parts() {
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

#[test]
fn adjusted_len_applies_only_requested_indel_and_clip_adjustments() {
    // Fragment interval [100,200) has aligned length 100.
    //
    // Indel adjustments:
    // - insertions: 4 + 1 = 5
    // - deletions: 6 + 2 = 8
    // => indel-adjusted length = 100 + 5 - 8 = 97
    //
    // Soft clips:
    // - left 3 bp
    // - right 2 bp
    // => clip-adjusted length = 97 + 3 + 2 = 102
    let fragment = FragmentWithIndelCounts {
        tid: 0,
        interval: Interval::new(100, 200).expect("test fragment interval should be valid"),
        left_soft_clip_bp: 3,
        right_soft_clip_bp: 2,
        deletions_nonoverlap: 6,
        insertions_nonoverlap: 4,
        deletions_overlap_supported: 2,
        insertions_overlap_supported: 1,
    };

    assert_eq!(
        fragment.adjusted_len(IndelMode::Ignore, ClipMode::Aligned),
        100
    );
    assert_eq!(
        fragment.adjusted_len(IndelMode::Skip, ClipMode::Aligned),
        100
    );
    assert_eq!(
        fragment.adjusted_len(IndelMode::Adjust, ClipMode::Aligned),
        97
    );
    assert_eq!(
        fragment.adjusted_len(IndelMode::Ignore, ClipMode::Adjust),
        105
    );
    assert_eq!(
        fragment.adjusted_len(IndelMode::Adjust, ClipMode::Adjust),
        102
    );
}

#[test]
fn assignment_interval_with_clip_mode_shifts_only_in_adjust_mode() {
    // Aligned interval is [100,200).
    // With left/right soft clips 3 and 2:
    // - aligned/skip keep [100,200)
    // - adjust expands to [97,202)
    let fragment = FragmentWithIndelCounts {
        tid: 0,
        interval: Interval::new(100, 200).expect("test fragment interval should be valid"),
        left_soft_clip_bp: 3,
        right_soft_clip_bp: 2,
        deletions_nonoverlap: 0,
        insertions_nonoverlap: 0,
        deletions_overlap_supported: 0,
        insertions_overlap_supported: 0,
    };

    assert_eq!(
        fragment
            .assignment_interval_with_clip_mode(ClipMode::Aligned)
            .expect("aligned interval"),
        Interval::new(100_u64, 200_u64).expect("expected aligned interval")
    );
    assert_eq!(
        fragment
            .assignment_interval_with_clip_mode(ClipMode::Skip)
            .expect("skip interval"),
        Interval::new(100_u64, 200_u64).expect("expected aligned interval")
    );
    assert_eq!(
        fragment
            .assignment_interval_with_clip_mode(ClipMode::Adjust)
            .expect("adjust interval"),
        Interval::new(97_u64, 202_u64).expect("expected expanded interval")
    );
}

#[test]
fn soft_clip_limit_is_applied_independently_to_both_fragment_ends() {
    // The threshold is checked per relevant fragment end, not on the summed clipping.
    //
    // Case 1: left 3 bp and right 2 bp are both within a 4 bp limit => keep.
    // Case 2: left 5 bp exceeds a 4 bp limit => reject.
    // Case 3: left 4 bp and right 4 bp equal the threshold => keep.
    let within_limit = FragmentWithIndelCounts {
        tid: 0,
        interval: Interval::new(100, 200).expect("test fragment interval should be valid"),
        left_soft_clip_bp: 3,
        right_soft_clip_bp: 2,
        deletions_nonoverlap: 0,
        insertions_nonoverlap: 0,
        deletions_overlap_supported: 0,
        insertions_overlap_supported: 0,
    };
    let left_exceeds_limit = FragmentWithIndelCounts {
        tid: 0,
        interval: Interval::new(100, 200).expect("test fragment interval should be valid"),
        left_soft_clip_bp: 5,
        right_soft_clip_bp: 0,
        deletions_nonoverlap: 0,
        insertions_nonoverlap: 0,
        deletions_overlap_supported: 0,
        insertions_overlap_supported: 0,
    };
    let equals_limit_on_both_ends = FragmentWithIndelCounts {
        tid: 0,
        interval: Interval::new(100, 200).expect("test fragment interval should be valid"),
        left_soft_clip_bp: 4,
        right_soft_clip_bp: 4,
        deletions_nonoverlap: 0,
        insertions_nonoverlap: 0,
        deletions_overlap_supported: 0,
        insertions_overlap_supported: 0,
    };

    assert!(within_limit.soft_clips_within_limit(4));
    assert!(!left_exceeds_limit.soft_clips_within_limit(4));
    assert!(equals_limit_on_both_ends.soft_clips_within_limit(4));
}

#[test]
fn deletion_base_limit_uses_total_supported_deletion_bases() {
    // The deletion limit is applied to the fragment-level total used for length adjustment.
    //
    // Case 1: 6 non-overlap + 2 supported overlap bases = 8, equal to the limit => keep.
    // Case 2: limit 7 is below the same 8 deleted reference bases => reject.
    let fragment = FragmentWithIndelCounts {
        tid: 0,
        interval: Interval::new(100, 200).expect("test fragment interval should be valid"),
        left_soft_clip_bp: 0,
        right_soft_clip_bp: 0,
        deletions_nonoverlap: 6,
        insertions_nonoverlap: 0,
        deletions_overlap_supported: 2,
        insertions_overlap_supported: 0,
    };

    assert_eq!(fragment.deletion_bases(), 8);
    assert!(fragment.deletion_bases_within_limit(8));
    assert!(!fragment.deletion_bases_within_limit(7));
}
