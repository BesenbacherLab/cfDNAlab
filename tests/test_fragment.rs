#[cfg(test)]
mod test_minimal_fragment {
    use cfdnalab::shared::fragment::minimal_fragment::*;
    use rust_htslib::bam::record::{Cigar, CigarString, Record};

    // Helpers

    fn cigar(ops: &[(char, u32)]) -> CigarString {
        let mut v = Vec::with_capacity(ops.len());
        for (op, len) in ops {
            v.push(match *op {
                'M' => Cigar::Match(*len),
                '=' => Cigar::Equal(*len),
                'X' => Cigar::Diff(*len),
                'I' => Cigar::Ins(*len),
                'D' => Cigar::Del(*len),
                'N' => Cigar::RefSkip(*len),
                'S' => Cigar::SoftClip(*len),
                'H' => Cigar::HardClip(*len),
                'P' => Cigar::Pad(*len),
                _ => panic!("unsupported CIGAR op {}", op),
            });
        }
        CigarString(v)
    }

    /// Build a minimal BAM record with given fields.
    /// **NOTE**: Depending on your rust-htslib version, you may need to tweak the `set` calls.
    fn mk_rec(
        tid: i32,
        pos: i64,
        is_rev: bool,
        cig: CigarString,
        seq_bytes: &[u8],
        qname: &[u8],
    ) -> Record {
        use rust_htslib::bam::record::Record;
        let mut r = Record::new();

        // Set core fields
        r.set_tid(tid);
        r.set_pos(pos);
        let mut flags: u16 = 0x1; // Paired
        if is_rev {
            flags |= 0x10;
        } // Reverse strand
        r.set_flags(flags);
        r.set_mapq(60);

        // Set qname, cigar, seq, qual
        r.set(
            qname,
            Some(&cig.clone()),
            seq_bytes,
            &vec![30u8; seq_bytes.len()],
        );

        r
    }

    // Make forward and reverse mates with simple defaults
    fn mk_pair_basic(
        tid: i32,
        f_pos: i64,
        f_cigar: CigarString,
        f_seq: &[u8],
        r_pos: i64,
        r_cigar: CigarString,
        r_seq: &[u8],
    ) -> (Record, Record) {
        let qname = b"pair1";
        let f = mk_rec(tid, f_pos, false, f_cigar, f_seq, qname);
        let r = mk_rec(tid, r_pos, true, r_cigar, r_seq, qname);
        (f, r)
    }

    fn collect_fragment_from_records_for_test(a: &Record, b: &Record) -> Option<Fragment> {
        collect_fragment(
            &MinimalReadInfo::from_record_with_gc_tag(a, None).ok()?,
            &MinimalReadInfo::from_record_with_gc_tag(b, None).ok()?,
        )
    }

    // Tests for `Fragment` (simple)

    #[test]
    fn test_collect_fragment_basic() {
        // Human verification status: unverified
        // forward: 100..150 (50M), reverse: 120..170 (50M) => fragment 100..170 (len=70)
        let (f, r) = mk_pair_basic(
            0,
            100,
            cigar(&[('M', 50)]),
            &vec![b'A'; 50],
            120,
            cigar(&[('M', 50)]),
            &vec![b'A'; 50],
        );
        let frag = collect_fragment_from_records_for_test(&f, &r).expect("fragment");
        assert_eq!(frag.tid, 0);
        assert_eq!(frag.start(), 100);
        assert_eq!(frag.end(), 170);
        assert_eq!(frag.len(), 70);
    }

    #[test]
    fn test_collect_fragment_uses_directional_pair_span_not_interval_union() {
        // Human verification status: unverified
        // Forward 100..200 and reverse 150..180 are still an inward-oriented pair because the
        // forward start is left of the reverse start. Our fragment definition is directional:
        // [forward.start(), reverse.end()), so the fragment must end at 180 rather than taking
        // the union end 200.
        let (f, r) = mk_pair_basic(
            0,
            100,
            cigar(&[('M', 100)]),
            &vec![b'A'; 100],
            150,
            cigar(&[('M', 30)]),
            &vec![b'A'; 30],
        );
        let frag = collect_fragment_from_records_for_test(&f, &r).expect("fragment");
        assert_eq!(frag.start(), 100);
        assert_eq!(frag.end(), 180);
        assert_eq!(frag.len(), 80);
    }

    #[test]
    fn test_collect_fragment_invalid_orientation() {
        // Human verification status: unverified
        // Both on same strand -> None
        let f = mk_rec(0, 100, false, cigar(&[('M', 30)]), &vec![b'A'; 30], b"1");
        let r = mk_rec(0, 140, false, cigar(&[('M', 30)]), &vec![b'A'; 30], b"1");
        assert!(collect_fragment_from_records_for_test(&f, &r).is_none());

        // End <= start -> None
        let f = mk_rec(0, 200, false, cigar(&[('M', 10)]), &vec![b'A'; 10], b"2");
        let r2 = mk_rec(0, 150, true, cigar(&[('M', 10)]), &vec![b'A'; 10], b"2");
        assert!(collect_fragment_from_records_for_test(&f, &r2).is_none());
    }
}

#[cfg(test)]
mod test_segmented_fragments {

    use cfdnalab::shared::fragment::{
        minimal_fragment::oriented_pair_from_read_info, segment_fragment::*,
    };
    use cfdnalab::shared::interval::Interval;

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
        // Human verification status: unverified
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
        // Human verification status: unverified
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
        // Human verification status: unverified
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
        // Human verification status: unverified
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
        // Human verification status: unverified
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
        // Human verification status: unverified
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

#[cfg(all(test, feature = "cmd_bam_to_frag"))]
mod test_frag_file_fragment {
    use cfdnalab::shared::fragment::frag_file_fragment::{
        FragReadInfo, collect_fragment_with_frag_file_info,
    };
    use cfdnalab::shared::interval::Interval;

    #[test]
    fn collect_frag_file_fragment_uses_directional_pair_span_not_interval_union() {
        // Human verification status: unverified
        // Forward 100..200 and reverse 150..180 form an inward pair. The fragment written to the
        // frag file must follow the cfDNA definition [forward.start(), reverse.end()) rather than
        // the aligned-interval union.
        let forward = FragReadInfo {
            tid: 0,
            interval: Interval::new(100, 200).expect("test read interval should be valid"),
            is_reverse: false,
            mapq: 60,
            strand: '+',
            is_read_1: true,
        };
        let reverse = FragReadInfo {
            tid: 0,
            interval: Interval::new(150, 180).expect("test read interval should be valid"),
            is_reverse: true,
            mapq: 50,
            strand: '-',
            is_read_1: false,
        };

        let frag = collect_fragment_with_frag_file_info(&forward, &reverse).expect("fragment");
        assert_eq!(frag.start(), 100);
        assert_eq!(frag.end(), 180);
        assert_eq!(frag.len(), 80);
        assert_eq!(frag.min_mapq, 50);
        assert_eq!(frag.read1_strand, '+');
    }
}

#[cfg(all(test, feature = "cmd_lengths"))]
mod tests_fragment_with_indel_counts {
    use cfdnalab::shared::fragment::indel_counting_fragment::*;
    use cfdnalab::shared::fragment::minimal_fragment::{
        is_inwards_oriented, oriented_pair_from_read_info,
    };
    use cfdnalab::shared::interval::Interval;
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
        vec![Cigar::Match(len as u32)]
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
        let a_info = IndelReadInfo::try_from(a).ok()?;
        let b_info = IndelReadInfo::try_from(b).ok()?;
        collect_fragment_with_indel_counts(&a_info, &b_info, skip_indels, count_indels)
    }

    #[test]
    fn indelreadinfo_parses_deletions_and_insertions() {
        // Human verification status: unverified
        // Forward read: start 100, cigar M50 D5 M45 => deletions at [150,155)
        let r = make_rec(0, 100, false, m_del_m(50, 5, 45));
        let info = IndelReadInfo::try_from(&r).expect("test indel read should be valid");
        assert_eq!(info.start(), 100);
        assert_eq!(info.end(), 100 + 50 + 5 + 45);
        assert_eq!(
            info.deletions,
            vec![Interval::new(150, 155).expect("test deletion should be valid")]
        );
        assert!(info.insertions.is_empty());

        // Reverse read: start 200, cigar M30 I4 M20 => insertion at ref pos 230 length 4
        let r2 = make_rec(0, 200, true, m_ins_m(30, 4, 20));
        let info2 = IndelReadInfo::try_from(&r2).expect("test indel read should be valid");
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
        // Human verification status: unverified
        // Forward at 100..160, Reverse at 150..210 (inward: forward.pos <= reverse.pos)
        let f = IndelReadInfo::try_from(&make_rec(0, 100, false, m(60)))
            .expect("test indel read should be valid");
        let r = IndelReadInfo::try_from(&make_rec(0, 150, true, m(60)))
            .expect("test indel read should be valid");
        let (fwd, rev) = oriented_pair_from_read_info(&f, &r).unwrap();
        assert!(is_inwards_oriented(fwd, rev));
    }

    #[test]
    fn collect_no_indels_fast_path() {
        // Human verification status: unverified
        // No indels; expect zero adjustments.
        let f = make_rec(0, 100, false, m(60));
        let r = make_rec(0, 180, true, m(40));
        let frag =
            collect_fragment_with_indel_counts_from_records_for_test(&f, &r, false, true).unwrap();
        assert_eq!(frag.tid, 0);
        assert_eq!(frag.start(), 100);
        assert_eq!(frag.end(), 220);
        assert_eq!(frag.len_ref(), 120);
        assert_eq!(frag.deletions_nonoverlap, 0);
        assert_eq!(frag.insertions_nonoverlap, 0);
        assert_eq!(frag.deletions_overlap_supported, 0);
        assert_eq!(frag.insertions_overlap_supported, 0);
        assert_eq!(frag.len_indel_adjusted(), frag.len_ref());
    }

    #[test]
    fn collect_skip_indels_filters_out() {
        // Human verification status: unverified
        // Has insertion; skip_indels=true => None.
        let f = make_rec(0, 100, false, m_ins_m(20, 3, 20));
        let r = make_rec(0, 140, true, m(40));
        assert!(
            collect_fragment_with_indel_counts_from_records_for_test(&f, &r, true, true).is_none()
        );
    }

    #[test]
    fn collect_count_indels_disabled_returns_zeroed() {
        // Human verification status: unverified
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
        // Human verification status: unverified
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
        assert_eq!(frag.len_indel_adjusted(), expected);
    }

    #[test]
    fn overlap_deletion_counts_intersection_only() {
        // Human verification status: unverified
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
        assert_eq!(frag.len_indel_adjusted(), frag.len_ref() - 3);
    }

    #[test]
    fn deletion_crossing_overlap_start_splits_into_left_nonoverlap_and_supported_overlap() {
        // Human verification status: unverified
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
            deletions: vec![Interval::new(150, 170).expect("test deletion should be valid")],
            insertions: vec![],
        };
        let reverse = IndelReadInfo {
            tid: 0,
            interval: Interval::new(160, 180).expect("test read interval should be valid"),
            is_reverse: true,
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
        // Human verification status: unverified
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
            deletions: vec![Interval::new(165, 175).expect("test deletion should be valid")],
            insertions: vec![],
        };
        let reverse = IndelReadInfo {
            tid: 0,
            interval: Interval::new(160, 220).expect("test read interval should be valid"),
            is_reverse: true,
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
        // Human verification status: unverified
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
        // Human verification status: unverified
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

#[cfg(all(test, feature = "cmd_ends"))]
mod test_fragment_with_ends {
    use cfdnalab::commands::ends::config_structs::{ClipStrategy, KmerSource};
    use cfdnalab::shared::fragment::ends_fragment::{
        EndReadInfo, FragmentWithEnds, collect_fragment_with_ends,
        collect_fragment_with_ends_from_single_read,
    };
    use cfdnalab::shared::indel_mode::IndelMotifFilterPolicy;
    use rust_htslib::bam::Record;
    use rust_htslib::bam::record::{Aux, Cigar, CigarString};

    fn cigar(ops: &[(char, u32)]) -> CigarString {
        let mut v = Vec::with_capacity(ops.len());
        for (op, len) in ops {
            v.push(match *op {
                'M' => Cigar::Match(*len),
                '=' => Cigar::Equal(*len),
                'X' => Cigar::Diff(*len),
                'I' => Cigar::Ins(*len),
                'D' => Cigar::Del(*len),
                'N' => Cigar::RefSkip(*len),
                'S' => Cigar::SoftClip(*len),
                'H' => Cigar::HardClip(*len),
                'P' => Cigar::Pad(*len),
                _ => panic!("unsupported CIGAR op {}", op),
            });
        }
        CigarString(v)
    }

    fn make_record(
        qname: &[u8],
        tid: i32,
        pos: i64,
        is_reverse: bool,
        cig: CigarString,
        seq_bytes: &[u8],
        gc_weight: Option<f32>,
    ) -> Record {
        let mut r = Record::new();
        r.set_tid(tid);
        r.set_pos(pos);
        let mut flags: u16 = 0x1;
        if is_reverse {
            flags |= 0x10;
        }
        r.set_flags(flags);
        r.set_mapq(60);
        r.set(qname, Some(&cig), seq_bytes, &vec![30u8; seq_bytes.len()]);
        if let Some(weight) = gc_weight {
            r.push_aux(b"GC", Aux::Float(weight)).expect("set GC tag");
        }
        r
    }

    fn end_info(
        record: &Record,
        clip_strategy: ClipStrategy,
        k_within: usize,
        gc_tag: Option<&[u8]>,
    ) -> EndReadInfo {
        EndReadInfo::from_record_with_gc_tag(record, gc_tag, clip_strategy, k_within)
            .expect("end read info")
    }

    fn collect_pair(
        a: &Record,
        b: &Record,
        clip_strategy: ClipStrategy,
        source_within: KmerSource,
        indel_filter: IndelMotifFilterPolicy,
        k_within: usize,
        max_soft_clips: Option<u32>,
        gc_tag: Option<&[u8]>,
    ) -> Option<FragmentWithEnds> {
        let a_info = end_info(a, clip_strategy, k_within, gc_tag);
        let b_info = end_info(b, clip_strategy, k_within, gc_tag);
        collect_fragment_with_ends(
            &a_info,
            &b_info,
            clip_strategy,
            source_within,
            indel_filter,
            k_within,
            max_soft_clips,
        )
    }

    fn collect_single(
        record: &Record,
        clip_strategy: ClipStrategy,
        source_within: KmerSource,
        indel_filter: IndelMotifFilterPolicy,
        k_within: usize,
        max_soft_clips: Option<u32>,
        gc_tag: Option<&[u8]>,
    ) -> Option<FragmentWithEnds> {
        let read = end_info(record, clip_strategy, k_within, gc_tag);
        collect_fragment_with_ends_from_single_read(
            &read,
            clip_strategy,
            source_within,
            indel_filter,
            k_within,
            max_soft_clips,
        )
    }

    #[test]
    fn endreadinfo_reads_edge_clipping_hard_clips_and_gc_tag() {
        // Human verification status: unverified
        let record = make_record(
            b"edge",
            0,
            100,
            false,
            cigar(&[('H', 5), ('S', 3), ('M', 8), ('S', 2), ('H', 4)]),
            b"TTTACGTACGGGT",
            Some(2.5),
        );

        let info = end_info(&record, ClipStrategy::Aligned, 4, Some(b"GC"));
        assert_eq!(info.left_soft_clip_bp, 3);
        assert_eq!(info.right_soft_clip_bp, 2);
        assert!(info.has_hard_clip);
        assert_eq!(info.gc_tag.weight, Some(2.5));
    }

    #[test]
    fn aligned_paired_keeps_aligned_geometry_and_excludes_soft_clips_from_sequences() {
        // Human verification status: unverified
        let forward = make_record(
            b"pair_aligned",
            0,
            100,
            false,
            cigar(&[('S', 2), ('M', 6)]),
            b"TTACGTAC",
            None,
        );
        let reverse = make_record(
            b"pair_aligned",
            0,
            110,
            true,
            cigar(&[('M', 6), ('S', 2)]),
            b"CGTACGTT",
            None,
        );

        let fragment = collect_pair(
            &forward,
            &reverse,
            ClipStrategy::Aligned,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            None,
            None,
        )
        .expect("fragment");

        assert_eq!(fragment.start(), 100);
        assert_eq!(fragment.end(), 116);
        assert_eq!(fragment.assignment_start(), 100);
        assert_eq!(fragment.assignment_end(), 116);
        assert_eq!(
            fragment.left_end.as_ref().expect("left end").within_bases,
            b"ACGT".to_vec()
        );
        assert_eq!(
            fragment.right_end.as_ref().expect("right end").within_bases,
            b"TACG".to_vec()
        );
    }

    #[test]
    fn raw_paired_expands_assignment_interval_and_uses_soft_clipped_bases() {
        // Human verification status: unverified
        let forward = make_record(
            b"pair_raw",
            0,
            100,
            false,
            cigar(&[('S', 2), ('M', 6)]),
            b"TTACGTAC",
            None,
        );
        let reverse = make_record(
            b"pair_raw",
            0,
            110,
            true,
            cigar(&[('M', 6), ('S', 2)]),
            b"CGTACGTT",
            None,
        );

        let fragment = collect_pair(
            &forward,
            &reverse,
            ClipStrategy::Raw,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            None,
            None,
        )
        .expect("fragment");

        assert_eq!(fragment.start(), 100);
        assert_eq!(fragment.end(), 116);
        assert_eq!(fragment.assignment_start(), 98);
        assert_eq!(fragment.assignment_end(), 118);
        assert_eq!(
            fragment.left_end.as_ref().expect("left end").within_bases,
            b"TTAC".to_vec()
        );
        assert_eq!(
            fragment.right_end.as_ref().expect("right end").within_bases,
            b"CGTT".to_vec()
        );
    }

    #[test]
    fn aligned_kept_ends_store_aligned_boundary_positions() {
        // Human verification status: unverified
        let forward = make_record(
            b"pair_aligned_boundary",
            0,
            100,
            false,
            cigar(&[('S', 2), ('M', 6)]),
            b"TTACGTAC",
            None,
        );
        let reverse = make_record(
            b"pair_aligned_boundary",
            0,
            110,
            true,
            cigar(&[('M', 6), ('S', 2)]),
            b"CGTACGTT",
            None,
        );

        let fragment = collect_pair(
            &forward,
            &reverse,
            ClipStrategy::Aligned,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            None,
            None,
        )
        .expect("fragment");

        assert_eq!(
            fragment.left_end.as_ref().expect("left end").boundary_pos,
            100
        );
        assert_eq!(
            fragment.right_end.as_ref().expect("right end").boundary_pos,
            116
        );
    }

    #[test]
    fn raw_kept_ends_store_raw_shifted_boundary_positions() {
        // Human verification status: unverified
        let forward = make_record(
            b"pair_raw_boundary",
            0,
            100,
            false,
            cigar(&[('S', 2), ('M', 6)]),
            b"TTACGTAC",
            None,
        );
        let reverse = make_record(
            b"pair_raw_boundary",
            0,
            110,
            true,
            cigar(&[('M', 6), ('S', 2)]),
            b"CGTACGTT",
            None,
        );

        let fragment = collect_pair(
            &forward,
            &reverse,
            ClipStrategy::Raw,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            None,
            None,
        )
        .expect("fragment");

        assert_eq!(
            fragment.left_end.as_ref().expect("left end").boundary_pos,
            98
        );
        assert_eq!(
            fragment.right_end.as_ref().expect("right end").boundary_pos,
            118
        );
    }

    #[test]
    fn drop_paired_drops_one_clipped_end_but_keeps_other_end() {
        // Human verification status: unverified
        let forward = make_record(
            b"pair_drop_one",
            0,
            100,
            false,
            cigar(&[('S', 2), ('M', 6)]),
            b"TTACGTAC",
            None,
        );
        let reverse = make_record(
            b"pair_drop_one",
            0,
            110,
            true,
            cigar(&[('M', 6)]),
            b"AACCGG",
            None,
        );

        let fragment = collect_pair(
            &forward,
            &reverse,
            ClipStrategy::Drop,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            None,
            None,
        )
        .expect("fragment");

        assert!(fragment.left_end.is_none());
        assert!(fragment.right_end.is_some());
        assert_eq!(fragment.assignment_start(), 100);
        assert_eq!(fragment.assignment_end(), 116);
    }

    #[test]
    fn drop_paired_returns_none_when_both_ends_are_skipped() {
        // Human verification status: unverified
        let forward = make_record(
            b"pair_drop_both",
            0,
            100,
            false,
            cigar(&[('S', 2), ('M', 6)]),
            b"TTACGTAC",
            None,
        );
        let reverse = make_record(
            b"pair_drop_both",
            0,
            110,
            true,
            cigar(&[('M', 6), ('S', 2)]),
            b"CGTACGTT",
            None,
        );

        assert!(
            collect_pair(
                &forward,
                &reverse,
                ClipStrategy::Drop,
                KmerSource::Read,
                IndelMotifFilterPolicy::Auto,
                4,
                None,
                None,
            )
            .is_none()
        );
    }

    #[test]
    fn max_soft_clips_skips_end_and_drops_boundary_adjustment() {
        // Human verification status: unverified
        let forward = make_record(
            b"pair_max_clip",
            0,
            100,
            false,
            cigar(&[('S', 3), ('M', 6)]),
            b"TTTACGTAC",
            None,
        );
        let reverse = make_record(
            b"pair_max_clip",
            0,
            110,
            true,
            cigar(&[('M', 6)]),
            b"AACCGG",
            None,
        );

        let fragment = collect_pair(
            &forward,
            &reverse,
            ClipStrategy::Raw,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            Some(2),
            None,
        )
        .expect("fragment");

        assert!(fragment.left_end.is_none());
        assert_eq!(fragment.assignment_start(), 100);
        assert_eq!(fragment.assignment_end(), 116);
    }

    #[test]
    fn hard_clipped_pair_is_discarded() {
        // Human verification status: unverified
        let forward = make_record(
            b"pair_hard_clip",
            0,
            100,
            false,
            cigar(&[('H', 3), ('M', 6)]),
            b"ACGTAC",
            None,
        );
        let reverse = make_record(
            b"pair_hard_clip",
            0,
            110,
            true,
            cigar(&[('M', 6)]),
            b"AACCGG",
            None,
        );

        assert!(
            collect_pair(
                &forward,
                &reverse,
                ClipStrategy::Aligned,
                KmerSource::Read,
                IndelMotifFilterPolicy::Auto,
                4,
                None,
                None,
            )
            .is_none()
        );
    }

    #[test]
    fn left_indel_is_detected_inside_aligned_footprint() {
        // Human verification status: unverified
        let record = make_record(
            b"left_indel",
            0,
            100,
            false,
            cigar(&[('M', 2), ('I', 1), ('M', 5)]),
            b"ACGTACGT",
            None,
        );

        let info = end_info(&record, ClipStrategy::Aligned, 4, None);
        assert!(info.left_motif_has_indels);
        assert!(!info.right_motif_has_indels);
    }

    #[test]
    fn right_indel_is_detected_inside_aligned_footprint() {
        // Human verification status: unverified
        let record = make_record(
            b"right_indel",
            0,
            100,
            false,
            cigar(&[('M', 5), ('D', 1), ('M', 2)]),
            b"ACGTACG",
            None,
        );

        let info = end_info(&record, ClipStrategy::Aligned, 3, None);
        assert!(!info.left_motif_has_indels);
        assert!(info.right_motif_has_indels);
    }

    #[test]
    fn refskip_is_treated_as_indel_inside_the_motif_footprint() {
        // Human verification status: unverified
        // Left motif sees M2 then N3 before reaching 4 aligned bases, so it is indel-affected.
        // Right motif takes its first 4 aligned bases from the terminal M5 block and never reaches
        // the ref-skip.
        let record = make_record(
            b"refskip_left",
            0,
            100,
            false,
            cigar(&[('M', 2), ('N', 3), ('M', 5)]),
            b"ACGTACG",
            None,
        );

        let info = end_info(&record, ClipStrategy::Aligned, 4, None);
        assert!(info.left_motif_has_indels);
        assert!(!info.right_motif_has_indels);
    }

    #[test]
    fn pad_is_ignored_when_scanning_for_indels() {
        // Human verification status: unverified
        // Padding does not consume aligned motif bases and should not count as an indel.
        let record = make_record(
            b"pad_ignored",
            0,
            100,
            false,
            cigar(&[('M', 2), ('P', 3), ('M', 5)]),
            b"ACGTACG",
            None,
        );

        let info = end_info(&record, ClipStrategy::Aligned, 4, None);
        assert!(!info.left_motif_has_indels);
        assert!(!info.right_motif_has_indels);
    }

    #[test]
    fn indel_outside_aligned_footprint_is_ignored() {
        // Human verification status: unverified
        let record = make_record(
            b"indel_far",
            0,
            100,
            false,
            cigar(&[('M', 5), ('I', 1), ('M', 5)]),
            b"ACGTACGTACG",
            None,
        );

        let info = end_info(&record, ClipStrategy::Aligned, 4, None);
        assert!(!info.left_motif_has_indels);
    }

    #[test]
    fn raw_reduces_indel_footprint_after_soft_clipping() {
        // Human verification status: unverified
        let record = make_record(
            b"raw_reduces_indel",
            0,
            100,
            false,
            cigar(&[('S', 3), ('M', 2), ('I', 1), ('M', 5)]),
            b"TTTACGTACGT",
            None,
        );

        let aligned_info = end_info(&record, ClipStrategy::Aligned, 4, None);
        let raw_info = end_info(&record, ClipStrategy::Raw, 4, None);

        assert!(aligned_info.left_motif_has_indels);
        assert!(!raw_info.left_motif_has_indels);
    }

    #[test]
    fn auto_read_source_keeps_indel_affected_end() {
        // Human verification status: unverified
        let forward = make_record(
            b"auto_read",
            0,
            100,
            false,
            cigar(&[('M', 2), ('I', 1), ('M', 5)]),
            b"ACGTACGT",
            None,
        );
        let reverse = make_record(
            b"auto_read",
            0,
            110,
            true,
            cigar(&[('M', 6)]),
            b"AACCGG",
            None,
        );

        let fragment = collect_pair(
            &forward,
            &reverse,
            ClipStrategy::Aligned,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            None,
            None,
        )
        .expect("fragment");

        assert!(fragment.left_end.is_some());
        assert_eq!(fragment.assignment_start(), 100);
    }

    #[test]
    fn auto_reference_source_skips_indel_affected_end_and_keeps_aligned_boundary() {
        // Human verification status: unverified
        let forward = make_record(
            b"auto_ref",
            0,
            100,
            false,
            cigar(&[('M', 2), ('I', 1), ('M', 5)]),
            b"ACGTACGT",
            None,
        );
        let reverse = make_record(
            b"auto_ref",
            0,
            110,
            true,
            cigar(&[('M', 6)]),
            b"AACCGG",
            None,
        );

        let fragment = collect_pair(
            &forward,
            &reverse,
            ClipStrategy::Aligned,
            KmerSource::Reference,
            IndelMotifFilterPolicy::Auto,
            4,
            None,
            None,
        )
        .expect("fragment");

        assert!(fragment.left_end.is_none());
        assert_eq!(fragment.assignment_start(), 100);
    }

    #[test]
    fn skip_affected_end_in_raw_keeps_raw_assignment_boundary() {
        // Human verification status: unverified
        let forward = make_record(
            b"raw_indel_skip",
            0,
            100,
            false,
            cigar(&[('S', 2), ('M', 1), ('I', 1), ('M', 5)]),
            b"TTACGTACG",
            None,
        );
        let reverse = make_record(
            b"raw_indel_skip",
            0,
            110,
            true,
            cigar(&[('M', 6)]),
            b"AACCGG",
            None,
        );

        let fragment = collect_pair(
            &forward,
            &reverse,
            ClipStrategy::Raw,
            KmerSource::Read,
            IndelMotifFilterPolicy::SkipAffectedEnd,
            4,
            None,
            None,
        )
        .expect("fragment");

        assert!(fragment.left_end.is_none());
        assert_eq!(fragment.assignment_start(), 98);
    }

    #[test]
    fn auto_reference_source_in_raw_keeps_raw_assignment_boundary_when_skipping_end() {
        // Human verification status: unverified
        // Raw uses the soft-clipped bases in the motif but still skips this end in Auto mode when
        // the within source is reference and the aligned part of the motif crosses an insertion.
        let forward = make_record(
            b"raw_auto_ref",
            0,
            100,
            false,
            cigar(&[('S', 2), ('M', 1), ('I', 1), ('M', 5)]),
            b"TTACGTACG",
            None,
        );
        let reverse = make_record(
            b"raw_auto_ref",
            0,
            110,
            true,
            cigar(&[('M', 6)]),
            b"AACCGG",
            None,
        );

        let fragment = collect_pair(
            &forward,
            &reverse,
            ClipStrategy::Raw,
            KmerSource::Reference,
            IndelMotifFilterPolicy::Auto,
            4,
            None,
            None,
        )
        .expect("fragment");

        assert!(fragment.left_end.is_none());
        assert_eq!(fragment.assignment_start(), 98);
    }

    #[test]
    fn skip_affected_fragment_drops_whole_fragment() {
        // Human verification status: unverified
        let forward = make_record(
            b"skip_fragment",
            0,
            100,
            false,
            cigar(&[('M', 2), ('I', 1), ('M', 5)]),
            b"ACGTACGT",
            None,
        );
        let reverse = make_record(
            b"skip_fragment",
            0,
            110,
            true,
            cigar(&[('M', 6)]),
            b"AACCGG",
            None,
        );

        assert!(
            collect_pair(
                &forward,
                &reverse,
                ClipStrategy::Aligned,
                KmerSource::Read,
                IndelMotifFilterPolicy::SkipAffectedFragment,
                4,
                None,
                None,
            )
            .is_none()
        );
    }

    #[test]
    fn drop_clipping_outranks_indel_skip() {
        // Human verification status: unverified
        let forward = make_record(
            b"drop_beats_indel",
            0,
            100,
            false,
            cigar(&[('S', 2), ('M', 1), ('I', 1), ('M', 5)]),
            b"TTACGTACG",
            None,
        );
        let reverse = make_record(
            b"drop_beats_indel",
            0,
            110,
            true,
            cigar(&[('M', 6)]),
            b"AACCGG",
            None,
        );

        let fragment = collect_pair(
            &forward,
            &reverse,
            ClipStrategy::Drop,
            KmerSource::Read,
            IndelMotifFilterPolicy::SkipAffectedEnd,
            4,
            None,
            None,
        )
        .expect("fragment");

        assert!(fragment.left_end.is_none());
        assert_eq!(fragment.assignment_start(), 100);
    }

    #[test]
    fn max_soft_clips_outranks_indel_skip_and_drops_raw_assignment_boundary() {
        // Human verification status: unverified
        // The left end is both over the soft-clip threshold and indel-affected.
        // Max-soft-clips wins first, so assignment falls back to the aligned boundary rather than
        // keeping the raw outward shift.
        let forward = make_record(
            b"max_clip_beats_indel",
            0,
            100,
            false,
            cigar(&[('S', 2), ('M', 1), ('I', 1), ('M', 5)]),
            b"TTACGTACG",
            None,
        );
        let reverse = make_record(
            b"max_clip_beats_indel",
            0,
            110,
            true,
            cigar(&[('M', 6)]),
            b"AACCGG",
            None,
        );

        let fragment = collect_pair(
            &forward,
            &reverse,
            ClipStrategy::Raw,
            KmerSource::Reference,
            IndelMotifFilterPolicy::SkipAffectedEnd,
            4,
            Some(1),
            None,
        )
        .expect("fragment");

        assert!(fragment.left_end.is_none());
        assert_eq!(fragment.assignment_start(), 100);
    }

    #[test]
    fn unpaired_aligned_exposes_both_ends() {
        // Human verification status: unverified
        let record = make_record(
            b"single_aligned",
            0,
            100,
            false,
            cigar(&[('S', 2), ('M', 6), ('S', 2)]),
            b"TTACGTACGG",
            None,
        );

        let fragment = collect_single(
            &record,
            ClipStrategy::Aligned,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            None,
            None,
        )
        .expect("fragment");

        assert_eq!(fragment.start(), 100);
        assert_eq!(fragment.end(), 106);
        assert_eq!(fragment.assignment_start(), 100);
        assert_eq!(fragment.assignment_end(), 106);
        assert_eq!(
            fragment.left_end.as_ref().expect("left end").within_bases,
            b"ACGT".to_vec()
        );
        assert_eq!(
            fragment.right_end.as_ref().expect("right end").within_bases,
            b"GTAC".to_vec()
        );
    }

    #[test]
    fn unpaired_raw_expands_assignment_interval_on_both_sides() {
        // Human verification status: unverified
        let record = make_record(
            b"single_raw",
            0,
            100,
            false,
            cigar(&[('S', 2), ('M', 6), ('S', 2)]),
            b"TTACGTACGG",
            None,
        );

        let fragment = collect_single(
            &record,
            ClipStrategy::Raw,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            None,
            None,
        )
        .expect("fragment");

        assert_eq!(fragment.start(), 100);
        assert_eq!(fragment.end(), 106);
        assert_eq!(fragment.assignment_start(), 98);
        assert_eq!(fragment.assignment_end(), 108);
        assert_eq!(
            fragment.left_end.as_ref().expect("left end").within_bases,
            b"TTAC".to_vec()
        );
        assert_eq!(
            fragment.right_end.as_ref().expect("right end").within_bases,
            b"ACGG".to_vec()
        );
    }

    #[test]
    fn hard_clipped_single_read_is_discarded() {
        // Human verification status: unverified
        let record = make_record(
            b"single_hard_clip",
            0,
            100,
            false,
            cigar(&[('H', 2), ('M', 6)]),
            b"ACGTAC",
            None,
        );

        assert!(
            collect_single(
                &record,
                ClipStrategy::Aligned,
                KmerSource::Read,
                IndelMotifFilterPolicy::Auto,
                4,
                None,
                None,
            )
            .is_none()
        );
    }

    #[test]
    fn paired_gc_tag_combines_to_fragment_average() {
        // Human verification status: unverified
        let forward = make_record(
            b"pair_gc",
            0,
            100,
            false,
            cigar(&[('M', 6)]),
            b"ACGTAC",
            Some(2.0),
        );
        let reverse = make_record(
            b"pair_gc",
            0,
            110,
            true,
            cigar(&[('M', 6)]),
            b"AACCGG",
            Some(4.0),
        );

        let fragment = collect_pair(
            &forward,
            &reverse,
            ClipStrategy::Aligned,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            None,
            Some(b"GC"),
        )
        .expect("fragment");

        assert_eq!(fragment.gc_tag.weight, Some(3.0));
    }

    #[test]
    fn unpaired_gc_tag_is_preserved() {
        // Human verification status: unverified
        let record = make_record(
            b"single_gc",
            0,
            100,
            false,
            cigar(&[('M', 6)]),
            b"ACGTAC",
            Some(2.5),
        );

        let fragment = collect_single(
            &record,
            ClipStrategy::Aligned,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            None,
            Some(b"GC"),
        )
        .expect("fragment");

        assert_eq!(fragment.gc_tag.weight, Some(2.5));
    }

    #[test]
    fn outwards_or_same_strand_pairs_are_rejected() {
        // Human verification status: unverified
        let forward = make_record(
            b"bad_pair",
            0,
            120,
            false,
            cigar(&[('M', 6)]),
            b"ACGTAC",
            None,
        );
        let reverse = make_record(
            b"bad_pair",
            0,
            100,
            true,
            cigar(&[('M', 6)]),
            b"AACCGG",
            None,
        );
        assert!(
            collect_pair(
                &forward,
                &reverse,
                ClipStrategy::Aligned,
                KmerSource::Read,
                IndelMotifFilterPolicy::Auto,
                4,
                None,
                None,
            )
            .is_none()
        );

        let same_strand = make_record(
            b"bad_pair",
            0,
            130,
            false,
            cigar(&[('M', 6)]),
            b"TTTTTT",
            None,
        );
        assert!(
            collect_pair(
                &forward,
                &same_strand,
                ClipStrategy::Aligned,
                KmerSource::Read,
                IndelMotifFilterPolicy::Auto,
                4,
                None,
                None,
            )
            .is_none()
        );
    }

    #[test]
    fn insufficient_sequence_skips_end_and_returns_none_if_both_ends_fail() {
        // Human verification status: unverified
        let record = make_record(
            b"too_short",
            0,
            100,
            false,
            cigar(&[('M', 3)]),
            b"ACG",
            None,
        );

        assert!(
            collect_single(
                &record,
                ClipStrategy::Aligned,
                KmerSource::Read,
                IndelMotifFilterPolicy::Auto,
                4,
                None,
                None,
            )
            .is_none()
        );
    }
}

#[cfg(test)]
mod test_kmer_segments {
    use cfdnalab::shared::fragment::segment_kmer_fragment::{
        FragmentWithKmerSegments, KmerSegmentedReadInfo, collect_fragment_with_kmer_segments,
    };
    use cfdnalab::shared::indel_mode::IndelMode;
    use cfdnalab::shared::interval::Interval;
    use rust_htslib::bam::record::{Cigar, CigarString, Record};
    fn read_len(cigar: &[Cigar]) -> usize {
        cigar
            .iter()
            .map(|c| match *c {
                Cigar::Match(l)
                | Cigar::Equal(l)
                | Cigar::Diff(l)
                | Cigar::Ins(l)
                | Cigar::SoftClip(l) => l as usize,
                _ => 0,
            })
            .sum()
    }
    fn make_record(tid: i32, pos: i64, is_reverse: bool, cigar_ops: &[Cigar]) -> Record {
        let mut record = Record::new();
        record.set_tid(tid);
        record.set_pos(pos);
        record.set_mapq(60);
        let mut flags: u16 = 0;
        if is_reverse {
            flags |= 0x10;
        }
        record.set_flags(flags);
        let cigar = CigarString(cigar_ops.to_vec());
        let seq_len = read_len(cigar_ops);
        let seq = vec![b'A'; seq_len.max(1)];
        let qual = vec![30u8; seq.len()];
        record.set(b"pair", Some(&cigar), &seq, &qual);
        record
    }
    fn collect_pair(
        forward: &Record,
        reverse: &Record,
        indel_mode: IndelMode,
        include_gap: bool,
        end_offset: u32,
    ) -> Option<FragmentWithKmerSegments> {
        let capture_segments = matches!(indel_mode, IndelMode::Adjust);
        let f_info = KmerSegmentedReadInfo::from_record(forward, capture_segments, None)
            .expect("test k-mer segmented read should be valid");
        let r_info = KmerSegmentedReadInfo::from_record(reverse, capture_segments, None)
            .expect("test k-mer segmented read should be valid");
        collect_fragment_with_kmer_segments(&f_info, &r_info, indel_mode, include_gap, end_offset)
    }
    fn segments(frag: &FragmentWithKmerSegments) -> Vec<(u32, u32)> {
        frag.segments.iter().map(Interval::as_tuple).collect()
    }
    #[test]
    fn ignore_mode_without_gap_tracks_per_read_spans() {
        // Human verification status: unverified
        let forward = make_record(0, 100, false, &[Cigar::Match(30)]);
        let reverse = make_record(0, 140, true, &[Cigar::Match(30)]);
        let frag = collect_pair(&forward, &reverse, IndelMode::Ignore, false, 0).expect("fragment");
        assert_eq!(frag.start(), 100);
        assert_eq!(frag.end(), 170);
        assert_eq!(segments(&frag), vec![(100, 130), (140, 170)]);
    }
    #[test]
    fn ignore_mode_includes_gap_when_requested() {
        // Human verification status: unverified
        let forward = make_record(0, 100, false, &[Cigar::Match(30)]);
        let reverse = make_record(0, 140, true, &[Cigar::Match(30)]);
        let frag = collect_pair(&forward, &reverse, IndelMode::Ignore, true, 0).expect("fragment");
        assert_eq!(frag.start(), 100);
        assert_eq!(frag.end(), 170);
        assert_eq!(segments(&frag), vec![(100, 170)]);
    }
    #[test]
    fn skip_mode_filters_pairs_with_indels() {
        // Human verification status: unverified
        let forward = make_record(
            0,
            100,
            false,
            &[Cigar::Match(10), Cigar::Ins(2), Cigar::Match(10)],
        );
        let reverse = make_record(0, 140, true, &[Cigar::Match(20)]);
        assert!(collect_pair(&forward, &reverse, IndelMode::Skip, true, 0).is_none());
    }
    #[test]
    fn adjust_mode_preserves_touching_segments_from_insertions() {
        // Human verification status: unverified
        let forward = make_record(
            0,
            100,
            false,
            &[Cigar::Match(10), Cigar::Ins(2), Cigar::Match(10)],
        );
        let reverse = make_record(0, 150, true, &[Cigar::Match(15)]);
        let frag = collect_pair(&forward, &reverse, IndelMode::Adjust, false, 0).expect("fragment");
        assert_eq!(frag.start(), 100);
        assert_eq!(frag.end(), 165);
        assert_eq!(segments(&frag), vec![(100, 110), (110, 120), (150, 165)]);
    }
    #[test]
    fn inter_mate_gap_not_merged_when_border_has_insertion() {
        // Human verification status: unverified
        let forward = make_record(0, 100, false, &[Cigar::Match(20), Cigar::Ins(2)]);
        let reverse = make_record(0, 130, true, &[Cigar::Match(20)]);
        let frag = collect_pair(&forward, &reverse, IndelMode::Adjust, true, 0).expect("fragment");
        assert_eq!(frag.start(), 100);
        assert_eq!(frag.end(), 150);
        // Not merged at the left boundary but merged on the right
        assert_eq!(segments(&frag), vec![(100, 120), (120, 150)]);
    }
    #[test]
    fn end_offset_trims_segments_but_preserves_span() {
        // Human verification status: unverified
        let forward = make_record(0, 100, false, &[Cigar::Match(40)]);
        let reverse = make_record(0, 120, true, &[Cigar::Match(40)]);
        let frag = collect_pair(&forward, &reverse, IndelMode::Ignore, true, 5).expect("fragment");
        assert_eq!(frag.start(), 100);
        assert_eq!(frag.end(), 160);
        assert_eq!(segments(&frag), vec![(105, 155)]);
    }

    #[test]
    fn adjust_mode_handles_mixed_insertions_and_deletions() {
        // Human verification status: unverified
        // Segment derivation under `IndelMode::Adjust`:
        // 1. Insertions terminate the current reference segment without advancing the reference.
        // 2. Deletions terminate the current segment and then advance the reference by the deleted span.
        // 3. With `include_inter_mate_gap = false`, we keep only the per-read mapped segments.
        //
        // Forward read starts at 100 with `8M 2I 5M 3D 4M 1I 3M`:
        // 1. `8M` covers [100,108), then `2I` closes that segment -> (100,108)
        // 2. `5M` resumes at 108 and covers [108,113), then `3D` closes it and skips [113,116) -> (108,113)
        // 3. `4M` resumes at 116 and covers [116,120), then `1I` closes it -> (116,120)
        // 4. `3M` resumes at 120 and reaches read end -> (120,123)
        //
        // Reverse read starts at 130 with `6M 2D 4M 2I 5M 1D 3M`:
        // 5. `6M` covers [130,136), then `2D` closes it and skips [136,138) -> (130,136)
        // 6. `4M` resumes at 138 and covers [138,142), then `2I` closes it -> (138,142)
        // 7. `5M` resumes at 142 and covers [142,147), then `1D` closes it and skips [147,148) -> (142,147)
        // 8. `3M` resumes at 148 and reaches read end -> (148,151)
        //
        // No segments overlap or even touch across the mate boundary, so the final fragment keeps
        // these eight segment boundaries verbatim.
        let forward = make_record(
            0,
            100,
            false,
            &[
                Cigar::Match(8),
                Cigar::Ins(2),
                Cigar::Match(5),
                Cigar::Del(3),
                Cigar::Match(4),
                Cigar::Ins(1),
                Cigar::Match(3),
            ],
        );
        let reverse = make_record(
            0,
            130,
            true,
            &[
                Cigar::Match(6),
                Cigar::Del(2),
                Cigar::Match(4),
                Cigar::Ins(2),
                Cigar::Match(5),
                Cigar::Del(1),
                Cigar::Match(3),
            ],
        );

        let frag = collect_pair(&forward, &reverse, IndelMode::Adjust, false, 0).expect("fragment");

        assert_eq!(frag.start(), 100);
        assert_eq!(frag.end(), 151);
        assert_eq!(
            segments(&frag),
            vec![
                (100, 108),
                (108, 113),
                (116, 120),
                (120, 123),
                (130, 136),
                (138, 142),
                (142, 147),
                (148, 151),
            ]
        );
    }

    #[test]
    fn gap_kept_separate_when_both_mates_have_boundary_insertions() {
        // Human verification status: unverified
        let forward = make_record(0, 100, false, &[Cigar::Match(15), Cigar::Ins(2)]);
        let reverse = make_record(0, 140, true, &[Cigar::Ins(3), Cigar::Match(18)]);

        let frag = collect_pair(&forward, &reverse, IndelMode::Adjust, true, 0).expect("fragment");

        assert_eq!(frag.start(), 100);
        assert_eq!(frag.end(), 158);
        assert_eq!(segments(&frag), vec![(100, 115), (115, 140), (140, 158)]);
    }

    #[test]
    fn returns_none_for_non_inward_orientation() {
        // Human verification status: unverified
        let forward = make_record(0, 100, false, &[Cigar::Match(20)]);
        let reverse_same_strand = make_record(0, 140, false, &[Cigar::Match(20)]);
        assert!(collect_pair(&forward, &reverse_same_strand, IndelMode::Ignore, true, 0).is_none());

        let reverse_left_of_forward = make_record(0, 80, true, &[Cigar::Match(20)]);
        assert!(
            collect_pair(
                &forward,
                &reverse_left_of_forward,
                IndelMode::Ignore,
                true,
                0
            )
            .is_none()
        );
    }

    #[test]
    fn end_offset_removing_all_sequence_returns_none() {
        // Human verification status: unverified
        let forward = make_record(0, 100, false, &[Cigar::Match(20)]);
        let reverse = make_record(0, 120, true, &[Cigar::Match(20)]);

        assert!(collect_pair(&forward, &reverse, IndelMode::Ignore, true, 20).is_none());
    }

    #[test]
    fn overlap_consensus_indel_behaviour() {
        // Human verification status: unverified
        // Segment derivation under `IndelMode::Adjust`:
        // 1. Each insertion closes the current segment but does not advance the reference.
        // 2. Each deletion closes the current segment and removes the deleted bases from the safe
        //    reference coverage.
        // 3. After clipping each read into segments, we take the union of both mates' segments and
        //    merge only genuinely overlapping spans. Touching spans stay separate.
        //
        // Forward read at 100:
        // 1. `10M` -> (100,110)
        // 2. `1I` closes the first segment.
        // 3. `8M 7M` are contiguous on the reference, so they stay one segment -> (110,125)
        // 4. `1I` closes that segment.
        // 5. `5M` -> (125,130)
        // 6. `2D` skips [130,132)
        // 7. `3M` -> (132,135)
        // 8. `1I` closes that segment.
        // 9. `3M` -> (135,138)
        // 10. `1D` skips [138,139)
        // 11. trailing `1M` -> (139,140)
        //
        // Reverse read at 120:
        // 12. `5M` -> (120,125)
        // 13. `1I` closes the first segment.
        // 14. `5M` -> (125,130)
        // 15. `2D` skips [130,132)
        // 16. `8M 5M` are contiguous on the reference, so they stay one segment -> (132,145)
        // 17. `1I` closes that segment.
        // 18. trailing `15M` -> (145,160)
        //
        // Merging both mates' segments gives:
        // - (100,110) by itself
        // - forward (110,125) overlaps reverse (120,125) -> merged to (110,125)
        // - forward (125,130) overlaps reverse (125,130) -> merged to (125,130)
        // - reverse (132,145) covers forward (132,135), (135,138), and (139,140), so they merge
        //   into one span while still leaving the deleted gap [130,132) absent -> (132,145)
        // - (145,160) only touches the previous span at 145, so it stays separate
        //
        // Because the mates already overlap (`forward.end() = 140`, `reverse.start() = 120`),
        // `include_inter_mate_gap = true` cannot introduce any extra gap segment. Trimming by
        // `end_offset = 12` then clips the final union [100,160) to [112,148), removing the first
        // segment entirely and shortening the last one to (145,148).
        let forward = make_record(
            0,
            100,
            false,
            &[
                Cigar::Match(10), // 100-110
                Cigar::Ins(1),    // insertion left of overlap
                Cigar::Match(8),  // 110-118
                Cigar::Match(7),  // 118-125
                Cigar::Ins(1),    // agreed insertion in overlap
                Cigar::Match(5),  // 125-130
                Cigar::Del(2),    // agreed deletion in overlap (130-132)
                Cigar::Match(3),  // 132-135
                Cigar::Ins(1),    // non-agreed insertion in overlap
                Cigar::Match(3),  // 135-138
                Cigar::Del(1),    // non-agreed deletion in overlap (skip 138)
                Cigar::Match(1),  // 139-140
            ],
        );

        let reverse = make_record(
            0,
            120,
            true,
            &[
                Cigar::Match(5),  // 120-125
                Cigar::Ins(1),    // agreed insertion in overlap
                Cigar::Match(5),  // 125-130
                Cigar::Del(2),    // agreed deletion in overlap (130-132)
                Cigar::Match(8),  // 132-140
                Cigar::Match(5),  // 140-145
                Cigar::Ins(1),    // insertion beyond overlap (reverse tail)
                Cigar::Match(15), // 145-160
            ],
        );

        let frag = collect_pair(&forward, &reverse, IndelMode::Adjust, false, 0).expect("fragment");

        assert_eq!(frag.start(), 100);
        assert_eq!(frag.end(), 160);
        assert_eq!(
            segments(&frag),
            vec![(100, 110), (110, 125), (125, 130), (132, 145), (145, 160)]
        );

        let frag_asked_gap =
            collect_pair(&forward, &reverse, IndelMode::Adjust, true, 0).expect("fragment");
        assert_eq!(segments(&frag), segments(&frag_asked_gap));

        let frag_offset =
            collect_pair(&forward, &reverse, IndelMode::Adjust, false, 12).expect("fragment");

        assert_eq!(frag_offset.start(), 100);
        assert_eq!(frag_offset.end(), 160);
        assert_eq!(
            segments(&frag_offset),
            vec![(112, 125), (125, 130), (132, 145), (145, 148)]
        );
    }
}
