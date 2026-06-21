mod tests_indel_counting_fragment_2 {
    use super::super::{
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
            interval: Interval::new(100, 200).expect("test fragment interval should be valid"),
            left_soft_clip_bp: 3,
            right_soft_clip_bp: 2,
            deletions_nonoverlap: 0,
            insertions_nonoverlap: 0,
            deletions_overlap_supported: 0,
            insertions_overlap_supported: 0,
        };
        let left_exceeds_limit = FragmentWithIndelCounts {
            interval: Interval::new(100, 200).expect("test fragment interval should be valid"),
            left_soft_clip_bp: 5,
            right_soft_clip_bp: 0,
            deletions_nonoverlap: 0,
            insertions_nonoverlap: 0,
            deletions_overlap_supported: 0,
            insertions_overlap_supported: 0,
        };
        let equals_limit_on_both_ends = FragmentWithIndelCounts {
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
}

mod tests_indel_counting_fragment {
    #![cfg(feature = "cmd_lengths")]

    use crate::shared::clip_mode::ClipMode;
    use crate::shared::fragment::indel_counting_fragment::FragmentWithIndelCounts;
    use crate::shared::fragment_iterators::fragments_with_indel_counts_from_bam;
    use crate::shared::indel_mode::IndelMode;
    use rust_htslib::bam;
    use rust_htslib::bam::record::Cigar;

    fn make_record(qname: &[u8], pos: i64, len: u32, is_reverse: bool) -> bam::Record {
        let mut rec = bam::Record::new();
        rec.set(qname, None, b"AAAAA", b"IIIII");
        rec.set_cigar(Some(&bam::record::CigarString::from(vec![Cigar::Match(
            len,
        )])));
        rec.set_tid(0);
        rec.set_pos(pos);
        rec.set_flags(if is_reverse { 16 } else { 0 });
        rec
    }

    fn collect_no_indel_pair(
        inspect_cigar: bool,
    ) -> (Vec<FragmentWithIndelCounts>, (u64, u64, u64, u64, u64)) {
        let records = vec![
            Ok(make_record(b"r1", 10, 5, false)), // forward [10,15)
            Ok(make_record(b"r1", 20, 5, true)),  // reverse [20,25)
        ];
        let include_read = |_r: &bam::Record| true;
        let mut iter = fragments_with_indel_counts_from_bam(
            records.into_iter(),
            include_read,
            IndelMode::Ignore,
            inspect_cigar,
            |_f: &FragmentWithIndelCounts| true,
            false,
        )
        .with_local_counters();
        let fragments: Vec<_> = iter.by_ref().map(|f| f.unwrap()).collect();
        let counters = iter.counters_snapshot();
        (
            fragments,
            (
                counters.incoming_reads,
                counters.accepted_forward_reads,
                counters.accepted_reverse_reads,
                counters.produced_fragments,
                counters.yielded_fragments,
            ),
        )
    }

    fn fragment_signature(
        fragment: &FragmentWithIndelCounts,
    ) -> ((u32, u32), u32, u32, u32, u32, u32, u32, u32) {
        (
            fragment.interval.as_tuple(),
            fragment.left_soft_clip_bp,
            fragment.right_soft_clip_bp,
            fragment.deletions_nonoverlap,
            fragment.insertions_nonoverlap,
            fragment.deletions_overlap_supported,
            fragment.insertions_overlap_supported,
            fragment.adjusted_len(IndelMode::Ignore, ClipMode::Aligned),
        )
    }

    #[test]
    fn disabled_cigar_inspection_matches_enabled_when_pair_has_no_indels_or_clips() {
        // Arrange + Act: Both reads are pure 5M alignments, so CIGAR inspection should not change the
        // assembled fragment. The fragment span is [10,25), giving length 15.
        let (with_cigar_fragments, with_cigar_counters) = collect_no_indel_pair(true);
        let (without_cigar_fragments, without_cigar_counters) = collect_no_indel_pair(false);

        // Assert: Both paths produce the same hand-derived fragment and the same iterator counters.
        let expected_fragments = vec![((10, 25), 0, 0, 0, 0, 0, 0, 15)];
        let expected_counters = (2, 1, 1, 1, 1);
        assert_eq!(
            with_cigar_fragments
                .iter()
                .map(fragment_signature)
                .collect::<Vec<_>>(),
            expected_fragments
        );
        assert_eq!(
            without_cigar_fragments
                .iter()
                .map(fragment_signature)
                .collect::<Vec<_>>(),
            expected_fragments
        );
        assert_eq!(with_cigar_counters, expected_counters);
        assert_eq!(without_cigar_counters, expected_counters);
    }

    #[test]
    fn yields_unpaired_fragments_and_respects_filter() {
        // Arrange
        let records = vec![
            Ok(make_record(b"r1", 5, 5, false)),  // length 5
            Ok(make_record(b"r2", 30, 5, false)), // length 5
        ];
        let include_read = |_r: &bam::Record| true;
        let mut iter = fragments_with_indel_counts_from_bam(
            records.into_iter(),
            include_read,
            crate::shared::indel_mode::IndelMode::Ignore,
            true,
            // Filter out anything shorter than 5 (none here) and longer than 5 (none here)
            |f: &FragmentWithIndelCounts| f.adjusted_len(IndelMode::Ignore, ClipMode::Aligned) <= 5,
            true,
        )
        .with_local_counters();

        // Act
        let frags: Vec<_> = iter.by_ref().map(|f| f.unwrap()).collect();

        // Assert
        assert_eq!(frags.len(), 2);
        let lengths: Vec<u32> = frags
            .iter()
            .map(|f| f.adjusted_len(IndelMode::Ignore, ClipMode::Aligned))
            .collect();
        assert_eq!(lengths, vec![5, 5]);
        let snap = iter.counters_snapshot();
        assert_eq!(snap.incoming_reads, 2);
        assert_eq!(snap.yielded_fragments, 2);
    }

    #[test]
    fn pairs_reads_and_yields_single_fragment() {
        // Arrange
        let forward = Ok(make_record(b"r1", 10, 5, false)); // end 15
        let reverse = Ok(make_record(b"r1", 20, 5, true)); // end 25
        let records = vec![forward, reverse];
        let include_read = |_r: &bam::Record| true;
        let mut iter = fragments_with_indel_counts_from_bam(
            records.into_iter(),
            include_read,
            crate::shared::indel_mode::IndelMode::Ignore,
            true,
            |_f: &FragmentWithIndelCounts| true,
            false,
        )
        .with_local_counters();

        // Act
        let frags: Vec<_> = iter.by_ref().map(|f| f.unwrap()).collect();

        // Assert
        assert_eq!(frags.len(), 1);
        assert_eq!(
            frags[0].adjusted_len(IndelMode::Ignore, ClipMode::Aligned),
            15
        ); // end(reverse) - start(forward)
        let snap = iter.counters_snapshot();
        assert_eq!(snap.incoming_reads, 2);
        assert_eq!(snap.produced_fragments, 1);
        assert_eq!(snap.yielded_fragments, 1);
    }
}

mod tests_fragment_with_indel_counts_3 {
    use crate::shared::clip_mode::ClipMode;
    use crate::shared::fragment::indel_counting_fragment::*;
    use crate::shared::fragment::minimal_fragment::{
        is_inwards_oriented, oriented_pair_from_read_info,
    };
    use crate::shared::indel_mode::IndelMode;
    use crate::shared::interval::Interval;
    use rust_htslib::bam::Record;
    use rust_htslib::bam::record::{Cigar, CigarString};

    // Helper: build a BAM record with given tid/pos/strand/cigar
    fn make_rec(tid: i32, pos: u32, is_reverse: bool, cigar: Vec<Cigar>) -> Record {
        let mut r = Record::new();
        r.set_tid(tid);
        r.set_pos(pos as i64);
        r.set_mapq(60);
        let mut flags: u16 = 0;
        if is_reverse {
            flags |= 0x10; // reverse strand
        }
        r.set_flags(flags);
        r.set_cigar(Some(&CigarString(cigar)));
        r
    }

    // Helper: simple M-only cigar of length l
    fn m(len: u32) -> Vec<Cigar> {
        vec![Cigar::Match(len)]
    }

    // Helper: cigar with ref deletion D at offset off of length d, total M before/after
    fn m_del_m(m1: u32, d: u32, m2: u32) -> Vec<Cigar> {
        vec![Cigar::Match(m1), Cigar::Del(d), Cigar::Match(m2)]
    }

    // Helper: cigar with insertion I at offset off of length i, total M before/after
    fn m_ins_m(m1: u32, i: u32, m2: u32) -> Vec<Cigar> {
        vec![Cigar::Match(m1), Cigar::Ins(i), Cigar::Match(m2)]
    }

    // Helper: quirky duplicate insertions at the same ref pos: I S I
    fn m_ins_s_ins_m(m1: u32, i1: u32, sc: u32, i2: u32, m2: u32) -> Vec<Cigar> {
        vec![
            Cigar::Match(m1),
            Cigar::Ins(i1),
            Cigar::SoftClip(sc),
            Cigar::Ins(i2),
            Cigar::Match(m2),
        ]
    }

    fn collect_fragment_with_indel_counts_from_records_for_test(
        a: &Record,
        b: &Record,
        skip_indels: bool,
        count_indels: bool,
    ) -> Option<FragmentWithIndelCounts> {
        let a_info = IndelReadInfo::from_record(a, true).ok()?;
        let b_info = IndelReadInfo::from_record(b, true).ok()?;
        collect_fragment_with_indel_counts(&a_info, &b_info, skip_indels, count_indels)
    }

    #[test]
    fn indelreadinfo_parses_deletions_and_insertions() {
        // Forward read: start 100, cigar M50 D5 M45 => deletions at [150,155)
        let r = make_rec(0, 100, false, m_del_m(50, 5, 45));
        let info = IndelReadInfo::from_record(&r, true).expect("test indel read should be valid");
        assert_eq!(info.start(), 100);
        assert_eq!(info.end(), 100 + 50 + 5 + 45);
        assert_eq!(
            info.deletions,
            vec![Interval::new(150, 155).expect("test deletion should be valid")]
        );
        assert!(info.insertions.is_empty());

        // Reverse read: start 200, cigar M30 I4 M20 => insertion at ref pos 230 length 4
        let r2 = make_rec(0, 200, true, m_ins_m(30, 4, 20));
        let info2 = IndelReadInfo::from_record(&r2, true).expect("test indel read should be valid");
        assert_eq!(
            info2.insertions,
            vec![InsertionAnchor {
                reference_position: 230,
                inserted_length: 4,
            }]
        );
        assert!(info2.deletions.is_empty());
    }

    #[test]
    fn orientation_and_inward_check() {
        // Forward at 100..160, Reverse at 150..210 (inward: forward.pos <= reverse.pos)
        let f = IndelReadInfo::from_record(&make_rec(0, 100, false, m(60)), true)
            .expect("test indel read should be valid");
        let r = IndelReadInfo::from_record(&make_rec(0, 150, true, m(60)), true)
            .expect("test indel read should be valid");
        let (fwd, rev) = oriented_pair_from_read_info(&f, &r).unwrap();
        assert!(is_inwards_oriented(fwd, rev));
    }

    #[test]
    fn collect_no_indels_fast_path() {
        // No indels; expect zero adjustments.
        let f = make_rec(0, 100, false, m(60));
        let r = make_rec(0, 180, true, m(40));
        let frag =
            collect_fragment_with_indel_counts_from_records_for_test(&f, &r, false, true).unwrap();
        assert_eq!(frag.start(), 100);
        assert_eq!(frag.end(), 220);
        assert_eq!(frag.len_ref(), 120);
        assert_eq!(frag.deletions_nonoverlap, 0);
        assert_eq!(frag.insertions_nonoverlap, 0);
        assert_eq!(frag.deletions_overlap_supported, 0);
        assert_eq!(frag.insertions_overlap_supported, 0);
        assert_eq!(
            frag.adjusted_len(IndelMode::Adjust, ClipMode::Aligned),
            frag.len_ref()
        );
    }

    #[test]
    fn collect_skip_indels_filters_out() {
        // Has insertion; skip_indels=true => None.
        let f = make_rec(0, 100, false, m_ins_m(20, 3, 20));
        let r = make_rec(0, 140, true, m(40));
        assert!(
            collect_fragment_with_indel_counts_from_records_for_test(&f, &r, true, true).is_none()
        );
    }

    #[test]
    fn collect_count_indels_disabled_returns_zeroed() {
        // Indels present but count_indels=false => fragment with zero adjustments.
        let f = make_rec(0, 100, false, m_ins_m(20, 3, 20));
        let r = make_rec(0, 140, true, m_del_m(10, 4, 30));
        let frag =
            collect_fragment_with_indel_counts_from_records_for_test(&f, &r, false, false).unwrap();
        assert_eq!(frag.deletions_nonoverlap, 0);
        assert_eq!(frag.insertions_nonoverlap, 0);
        assert_eq!(frag.deletions_overlap_supported, 0);
        assert_eq!(frag.insertions_overlap_supported, 0);
    }

    #[test]
    fn nonoverlap_indels_counted_fully() {
        // Non-overlapping mates: forward 100..120, reverse 140..160.
        // Forward has D3 at [110,113), Reverse has I4 at ref pos 150.
        let f = make_rec(0, 100, false, m_del_m(10, 3, 7)); // 100..120
        let r = make_rec(0, 140, true, m_ins_m(10, 4, 10)); // 140..160
        let frag =
            collect_fragment_with_indel_counts_from_records_for_test(&f, &r, false, true).unwrap();
        // No aligned overlap -> both indels are non-overlap
        assert_eq!(frag.deletions_nonoverlap, 3);
        assert_eq!(frag.insertions_nonoverlap, 4);
        assert_eq!(frag.deletions_overlap_supported, 0);
        assert_eq!(frag.insertions_overlap_supported, 0);
        // Adjusted length = (end-start) + 4 - 3
        let expected = frag.len_ref() + 1;
        assert_eq!(
            frag.adjusted_len(IndelMode::Adjust, ClipMode::Aligned),
            expected
        );
    }

    #[test]
    fn overlap_deletion_counts_intersection_only() {
        // Overlapping mates: forward 100..180, reverse 160..220 -> overlap [160,180)
        // Forward deletion [170,175); Reverse deletion [172,178) -> intersection [172,175) len 3.
        let f = make_rec(0, 100, false, m_del_m(70, 5, 5)); // del at [170,175)
        let r = make_rec(0, 160, true, m_del_m(12, 6, 42)); // del at [172,178)
        let frag =
            collect_fragment_with_indel_counts_from_records_for_test(&f, &r, false, true).unwrap();
        assert_eq!(frag.deletions_nonoverlap, 0);
        assert_eq!(frag.deletions_overlap_supported, 3);
        assert_eq!(frag.insertions_nonoverlap, 0);
        assert_eq!(frag.insertions_overlap_supported, 0);
        // Adjusted length subtracts 3
        assert_eq!(
            frag.adjusted_len(IndelMode::Adjust, ClipMode::Aligned),
            frag.len_ref() - 3
        );
    }

    #[test]
    fn deletion_crossing_overlap_start_splits_into_left_nonoverlap_and_supported_overlap() {
        // Forward 100..180 and reverse 160..180 give aligned overlap [160,180).
        //
        // Forward deletion [150,170) crosses the overlap start:
        // - left non-overlap part: [150,160) => 10 bp
        // - overlap part:         [160,170)
        //
        // Reverse deletion [165,175) lies fully inside the overlap.
        // Supported overlap deletion is the intersection of [160,170) and [165,175),
        // which is [165,170) => 5 bp.
        let forward = IndelReadInfo {
            tid: 0,
            interval: Interval::new(100, 180).expect("test read interval should be valid"),
            is_reverse: false,
            left_soft_clip_bp: 0,
            right_soft_clip_bp: 0,
            deletions: vec![Interval::new(150, 170).expect("test deletion should be valid")],
            insertions: vec![],
        };
        let reverse = IndelReadInfo {
            tid: 0,
            interval: Interval::new(160, 180).expect("test read interval should be valid"),
            is_reverse: true,
            left_soft_clip_bp: 0,
            right_soft_clip_bp: 0,
            deletions: vec![Interval::new(165, 175).expect("test deletion should be valid")],
            insertions: vec![],
        };

        let frag = collect_fragment_with_indel_counts(&forward, &reverse, false, true).unwrap();

        assert_eq!(frag.deletions_nonoverlap, 10);
        assert_eq!(frag.deletions_overlap_supported, 5);
        assert_eq!(frag.insertions_nonoverlap, 0);
        assert_eq!(frag.insertions_overlap_supported, 0);
    }

    #[test]
    fn deletion_crossing_overlap_end_splits_into_right_nonoverlap_and_supported_overlap() {
        // Forward 100..180 and reverse 160..220 give aligned overlap [160,180).
        //
        // Reverse deletion [170,190) crosses the overlap end:
        // - overlap part:         [170,180)
        // - right non-overlap:    [180,190) => 10 bp
        //
        // Forward deletion [165,175) lies fully inside the overlap.
        // Supported overlap deletion is the intersection of [165,175) and [170,180),
        // which is [170,175) => 5 bp.
        let forward = IndelReadInfo {
            tid: 0,
            interval: Interval::new(100, 180).expect("test read interval should be valid"),
            is_reverse: false,
            left_soft_clip_bp: 0,
            right_soft_clip_bp: 0,
            deletions: vec![Interval::new(165, 175).expect("test deletion should be valid")],
            insertions: vec![],
        };
        let reverse = IndelReadInfo {
            tid: 0,
            interval: Interval::new(160, 220).expect("test read interval should be valid"),
            is_reverse: true,
            left_soft_clip_bp: 0,
            right_soft_clip_bp: 0,
            deletions: vec![Interval::new(170, 190).expect("test deletion should be valid")],
            insertions: vec![],
        };

        let frag = collect_fragment_with_indel_counts(&forward, &reverse, false, true).unwrap();

        assert_eq!(frag.deletions_nonoverlap, 10);
        assert_eq!(frag.deletions_overlap_supported, 5);
        assert_eq!(frag.insertions_nonoverlap, 0);
        assert_eq!(frag.insertions_overlap_supported, 0);
    }

    #[test]
    fn overlap_insertions_require_both_mates_same_ref_pos() {
        // Overlap [160,180).
        // Forward insertion at ref 165 len 5; Reverse insertion at ref 165 len 3 -> min = 3 counted.
        // Forward insertion at ref 170 len 2; Reverse none at 170 -> discarded because overlap
        // insertions only count when both mates report an insertion at the same reference
        // position. It does not fall back into the non-overlap counter.
        let f = make_rec(0, 100, false, {
            let mut v = m_ins_m(65, 5, 4); // ins at 165
            v.extend(m_ins_m(1, 2, 10)); // then ins at 170 (since 100 + 65 + [I] + 4 + 1 + [I] + 10)
            v
        });
        let r = make_rec(0, 160, true, m_ins_m(5, 3, 15)); // ins at 165
        let frag =
            collect_fragment_with_indel_counts_from_records_for_test(&f, &r, false, true).unwrap();
        assert_eq!(frag.insertions_nonoverlap, 0); // Unpaired overlap insertions are discarded rather than counted elsewhere
        assert_eq!(frag.insertions_overlap_supported, 3); // min(5,3)
    }

    #[test]
    fn duplicate_insertions_at_same_pos_per_read_take_max_then_min_across_mates() {
        // Overlap [160,200).
        // Forward has two insertions at ref 170: lengths 2 and 5 (separated by soft-clip); keep max=5.
        // Reverse has insertion at ref 170 length 3 -> min(5,3) = 3 counted.
        let f = make_rec(0, 150, false, m_ins_s_ins_m(20, 2, 4, 5, 26)); // two I at ref 170
        let r = make_rec(0, 160, true, m_ins_m(10, 3, 30)); // I at ref 170
        let frag =
            collect_fragment_with_indel_counts_from_records_for_test(&f, &r, false, true).unwrap();
        assert_eq!(frag.insertions_overlap_supported, 3);
        assert_eq!(frag.insertions_nonoverlap, 0);
    }
}
