#[cfg(test)]
mod test_segmented_fragments {
    use crate::shared::fragment::{
        minimal_fragment::oriented_pair_from_read_info, segment_fragment::*,
    };
    use crate::shared::interval::Interval;

    // Tiny helper to construct SegmentedReadInfo without pulling in BAM types
    fn sri(
        tid: i32,
        pos: u32,
        end: u32,
        is_reverse: bool,
        has_ref_gap: bool,
        max_ref_gap: u32,
        segs: &[(u32, u32)],
    ) -> SegmentedReadInfo {
        SegmentedReadInfo {
            tid,
            interval: Interval::new(pos, end).expect("test read interval should be valid"),
            is_reverse,
            has_ref_gap,
            max_ref_gap,
            ref_mapped_segments: segs.to_vec(),
            gc_tag: Default::default(),
        }
    }

    fn segment_tuples(segments: &[Interval<u32>]) -> Vec<(u32, u32)> {
        segments.iter().map(Interval::as_tuple).collect()
    }

    #[test]
    fn oriented_pair_generic_orders_fw_rev() {
        let fwd = sri(0, 10, 20, false, false, 0, &[]);
        let rev = sri(0, 40, 50, true, false, 0, &[]);

        let (a, b) = oriented_pair_from_read_info(&fwd, &rev).expect("pair");
        assert!(!a.is_reverse);
        assert!(b.is_reverse);

        let (a2, b2) = oriented_pair_from_read_info(&rev, &fwd).expect("pair");
        assert!(!a2.is_reverse);
        assert!(b2.is_reverse);
    }

    #[test]
    fn collect_no_gaps_exclude_inter_mate_gap_segments_present() {
        // forward 10..20, reverse 40..50, no ref gaps, do NOT include the inter-mate gap
        let fwd = sri(0, 10, 20, false, false, 0, &[]);
        let rev = sri(0, 40, 50, true, false, 0, &[]);

        let fws = collect_fragment_with_segments(&fwd, &rev, 1, false).expect("pair");
        assert_eq!(fws.start(), 10);
        assert_eq!(fws.end(), 50);

        // Expect two explicit segments: [10,20] and [40,50]
        let segs = fws.segments.as_ref().expect("segments present");
        assert_eq!(segment_tuples(segs), vec![(10u32, 20u32), (40u32, 50u32)]);
    }

    #[test]
    fn collect_no_gaps_include_inter_mate_gap_plain_span() {
        // forward 10..20, reverse 40..50, include inter-mate gap -> whole span
        let fwd = sri(0, 10, 20, false, false, 0, &[]);
        let rev = sri(0, 40, 50, true, false, 0, &[]);

        let fws = collect_fragment_with_segments(&fwd, &rev, 1, true).expect("pair");
        assert_eq!(fws.start(), 10);
        assert_eq!(fws.end(), 50);
        assert!(fws.segments.is_none());
    }

    #[test]
    fn collect_with_ref_gap_include_inter_mate_gap_merges_right() {
        // forward has a deletion making two ref-mapped blocks: [10..20], [25..30]
        // reverse is [40..50], include inter-mate gap (30..40)
        let fwd = sri(0, 10, 30, false, true, 5, &[(0, 10), (15, 5)]);
        let rev = sri(0, 40, 50, true, false, 0, &[]);

        let fws = collect_fragment_with_segments(&fwd, &rev, 1, true).expect("pair");
        assert_eq!(fws.start(), 10);
        assert_eq!(fws.end(), 50);

        // Expected segments after adding gap and merging: [10..20], [25..50]
        let segs = fws.segments.as_ref().expect("segments present");
        assert_eq!(segment_tuples(segs), vec![(10u32, 20u32), (25u32, 50u32)]);
    }

    #[test]
    fn collect_with_ref_gap_exclude_inter_mate_gap_three_blocks() {
        // Same as above but exclude the inter-mate gap
        let fwd = sri(0, 10, 30, false, true, 5, &[(0, 10), (15, 5)]);
        let rev = sri(0, 40, 50, true, false, 0, &[]);

        let fws = collect_fragment_with_segments(&fwd, &rev, 1, false).expect("pair");
        let segs = fws.segments.as_ref().expect("segments present");

        // Expected: [10..20], [25..30], [40..50]
        assert_eq!(
            segment_tuples(segs),
            vec![(10u32, 20u32), (25u32, 30u32), (40u32, 50u32)]
        );
    }

    #[test]
    fn collect_segments_uses_directional_fragment_span_not_interval_union() {
        // Forward 100..200 and reverse 150..180 define a fragment span of [100,180) by project
        // convention, even though the aligned-interval union would extend to 200.
        let fwd = sri(0, 100, 200, false, false, 0, &[]);
        let rev = sri(0, 150, 180, true, false, 0, &[]);

        let frag = collect_fragment_with_segments(&fwd, &rev, 1, true).expect("pair");
        assert_eq!(frag.start(), 100);
        assert_eq!(frag.end(), 180);
        assert!(frag.segments.is_none());
    }
}
