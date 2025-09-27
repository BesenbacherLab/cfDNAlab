#[cfg(test)]
mod test_minimal_fragment {
    use cfdnalab::utils::fragment::minimal_fragment::*;
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

    // Tests for `Fragment` (simple)

    #[test]
    fn test_collect_fragment_basic() {
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
        let frag = collect_fragment_from_records(&f, &r).expect("fragment");
        assert_eq!(frag.tid, 0);
        assert_eq!(frag.start, 100);
        assert_eq!(frag.end, 170);
        assert_eq!(frag.len(), 70);
    }

    #[test]
    fn test_collect_fragment_invalid_orientation() {
        // Both on same strand -> None
        let f = mk_rec(0, 100, false, cigar(&[('M', 30)]), &vec![b'A'; 30], b"1");
        let r = mk_rec(0, 140, false, cigar(&[('M', 30)]), &vec![b'A'; 30], b"1");
        assert!(collect_fragment_from_records(&f, &r).is_none());

        // End <= start -> None
        let f = mk_rec(0, 200, false, cigar(&[('M', 10)]), &vec![b'A'; 10], b"2");
        let r2 = mk_rec(0, 150, true, cigar(&[('M', 10)]), &vec![b'A'; 10], b"2");
        assert!(collect_fragment_from_records(&f, &r2).is_none());
    }

    // Tests for `FragmentOverlapMM` (Match/Mismatch-only)

    // #[test]
    // fn test_overlap_mm_basic_match() {
    //     // Forward: 100..160 (60M); reverse: 140..200 (60M)
    //     // Overlap window: 140..160 (20 bp), all 'A'
    //     let (f, r) = mk_pair_basic(
    //         0,
    //         100,
    //         cigar(&[('M', 60)]),
    //         &vec![b'C'; 60], // Only positions 40..60 will contribute to overlap
    //         140,
    //         cigar(&[('M', 60)]),
    //         &vec![b'A'; 60], // Positions 0..20 -> overlap
    //     );

    //     let ov = FragmentOverlapMM::from_pair(&f, &r).expect("overlap mm");
    //     assert_eq!(ov.overlap.start, 140);
    //     assert_eq!(ov.overlap.end, 160);
    //     assert_eq!(ov.ref_coords.len(), 20);
    //     assert_eq!(ov.left_bases.len(), 20);
    //     assert_eq!(ov.right_bases.len(), 20);

    //     // Coordinates should be contiguous 140..160
    //     for (i, &coord) in ov.ref_coords.iter().enumerate() {
    //         assert_eq!(coord, 140 + i as u32);
    //     }
    //     // Compare that we pulled aligned columns; we didn't set exact sequences for left overlapping
    //     // region; set them now to 'A' to ensure a match:
    //     // Forward read bases 40..60 must be 'A' for equality; rebuild with that:
    //     let (f2, r2) = mk_pair_basic(
    //         0,
    //         100,
    //         cigar(&[('M', 40), ('M', 20)]),
    //         &[vec![b'C'; 40], vec![b'A'; 20]].concat(),
    //         140,
    //         cigar(&[('M', 60)]),
    //         &vec![b'A'; 60],
    //     );
    //     let ov2 = FragmentOverlapMM::from_pair(&f2, &r2).unwrap();
    //     assert!(
    //         ov2.left_bases
    //             .iter()
    //             .zip(&ov2.right_bases)
    //             .all(|(l, r)| l == r)
    //     );
    //     assert!(ov2.left_bases.iter().all(|&b| b == b'A'));
    // }

    // #[test]
    // fn test_overlap_mm_drops_insertions() {
    //     // Forward: 100..120 on ref with an insertion at 110: 10M 2I 10M
    //     // Reverse: 100..122 on ref: 22M (covers entire region)
    //     // MM extractor should DROP the insertion bases and align by ref coords (20 positions).
    //     let (f, r) = mk_pair_basic(
    //         0,
    //         100,
    //         cigar(&[('M', 10), ('I', 2), ('M', 10)]),
    //         // Seq length = 10 + 2 + 10 = 22
    //         &[vec![b'G'; 10], vec![b'T'; 2], vec![b'G'; 10]].concat(),
    //         100,
    //         cigar(&[('M', 22)]),
    //         &vec![b'G'; 22],
    //     );
    //     let ov = FragmentOverlapMM::from_pair(&f, &r).expect("overlap mm");
    //     assert_eq!(ov.overlap.start, 100);
    //     assert_eq!(ov.overlap.end, 120);
    //     assert_eq!(ov.ref_coords.len(), 20); // Insertion dropped (no ref coord)
    //     // Check exact reference coordinates are 100..120 (no duplicates at the insertion boundary)
    //     let expected: Vec<u32> = (100..120).collect();
    //     assert_eq!(ov.ref_coords, expected);
    //     assert!(ov.left_bases.iter().all(|&b| b == b'G'));
    //     assert!(ov.right_bases.iter().all(|&b| b == b'G'));
    // }

    // #[test]
    // fn test_overlap_mm_none_when_no_overlap() {
    //     // Forward: 100..120; reverse: 130..150 => no overlap
    //     let (f, r) = mk_pair_basic(
    //         0,
    //         100,
    //         cigar(&[('M', 20)]),
    //         &vec![b'A'; 20],
    //         130,
    //         cigar(&[('M', 20)]),
    //         &vec![b'A'; 20],
    //     );
    //     assert!(FragmentOverlapMM::from_pair(&f, &r).is_none());
    // }

    // #[test]
    // fn test_overlap_mm_skips_deletion_positions() {
    //     // Forward: 100..122 on ref with a deletion at 110..112: 10M 2D 10M
    //     // Reverse: 100..122 on ref: 22M
    //     // MM extractor should DROP the deletion columns (no read base) and keep only aligned coords.
    //     let (f, r) = mk_pair_basic(
    //         0,
    //         100,
    //         cigar(&[('M', 10), ('D', 2), ('M', 10)]),
    //         // seq length = 10 + 10 = 20 (D consumes reference only)
    //         &[vec![b'G'; 10], vec![b'G'; 10]].concat(),
    //         100,
    //         cigar(&[('M', 22)]),
    //         &[vec![b'G'; 10], vec![b'A'; 2], vec![b'G'; 10]].concat(),
    //     );
    //     let ov = FragmentOverlapMM::from_pair(&f, &r).expect("overlap mm");
    //     assert_eq!(ov.overlap.start, 100);
    //     assert_eq!(ov.overlap.end, 122);

    //     // Expect reference coords 100..110 and 112..122 (skip 110,111 where the deletion is)
    //     let mut expected: Vec<u32> = (100..110).collect();
    //     expected.extend(112..122);
    //     assert_eq!(ov.ref_coords, expected);

    //     // Bases should align and match
    //     assert_eq!(ov.left_bases.len(), expected.len());
    //     assert_eq!(ov.right_bases.len(), expected.len());
    //     assert!(ov.left_bases.iter().all(|&b| b == b'G'));
    //     assert!(ov.right_bases.iter().all(|&b| b == b'G'));
    // }

    // // Tests for `FragmentWithSequences` (within-fragment sequences)

    // #[test]
    // fn test_fragment_with_sequences_includes_insertions_and_flags() {
    //     // Fragment: forward 100..120 (20M), reverse 105..130 (25M) => fragment 100..130
    //     // Forward read contains insertion at ref 110: 10M 2I 10M
    //     // Reverse is plain 25M
    //     let (f, r) = mk_pair_basic(
    //         0,
    //         100,
    //         cigar(&[('M', 10), ('I', 2), ('M', 10)]),
    //         &[vec![b'A'; 10], vec![b'G'; 2], vec![b'A'; 10]].concat(),
    //         105,
    //         cigar(&[('M', 25)]),
    //         &vec![b'C'; 25],
    //     );
    //     let frag_seq = FragmentWithSequences::from_pair(&f, &r).expect("frag seq");
    //     // Left sequence should include the insertion (2 bases)
    //     assert_eq!(frag_seq.left_seq.len(), 22);
    //     assert!(frag_seq.left_info.has_insertion);
    //     assert!(!frag_seq.left_info.has_deletion);
    //     assert!(!frag_seq.left_info.has_refskip);
    // }

    // #[test]
    // fn test_fragment_with_sequences_flags_deletion() {
    //     // Forward has a deletion inside the fragment: 10M 2D 10M
    //     // Reverse covers the same span with 22M
    //     let (f, r) = mk_pair_basic(
    //         0,
    //         100,
    //         cigar(&[('M', 10), ('D', 2), ('M', 10)]),
    //         &[vec![b'A'; 10], vec![b'A'; 10]].concat(),
    //         100,
    //         cigar(&[('M', 22)]),
    //         &vec![b'C'; 22],
    //     );
    //     let frag_seq = FragmentWithSequences::from_pair(&f, &r).expect("frag seq");

    //     // Fragment spans full reference window
    //     assert_eq!(frag_seq.frag.start, 100);
    //     assert_eq!(frag_seq.frag.end, 122);

    //     // Left sequence length excludes deleted reference columns (still 20 bases emitted)
    //     assert_eq!(frag_seq.left_seq.len(), 20);
    //     assert_eq!(frag_seq.right_seq.len(), 22);
    //     // Flags: deletion seen on left; no insertion/softclip
    //     assert!(frag_seq.left_info.has_deletion);
    //     assert!(!frag_seq.left_info.has_insertion);
    //     assert!(!frag_seq.left_info.has_softclip);
    //     // Right read has no deletion
    //     assert!(!frag_seq.right_info.has_deletion);
    // }

    // #[test]
    // fn test_fragment_with_sequences_softclips_excluded_but_flagged() {
    //     // Forward: 5S 20M at pos=100 -> mapped span 100..120; soft clip before POS
    //     // Reverse: 120..145 (25M) => fragment 100..145
    //     let (f, r) = mk_pair_basic(
    //         0,
    //         100,
    //         cigar(&[('S', 5), ('M', 20)]),
    //         // seq length 25; first 5 are soft-clipped and must NOT appear in left_seq
    //         &[vec![b'T'; 5], vec![b'A'; 20]].concat(),
    //         120,
    //         cigar(&[('M', 25)]),
    //         &vec![b'C'; 25],
    //     );
    //     let frag_seq = FragmentWithSequences::from_pair(&f, &r).expect("frag seq");
    //     // The left read slice inside the fragment should have only the 20 aligned A's.
    //     assert_eq!(frag_seq.left_seq, vec![b'A'; 20]);
    //     assert!(frag_seq.left_info.has_softclip);
    //     // Right has no soft clips
    //     assert!(!frag_seq.right_info.has_softclip);
    // }

    // #[test]
    // fn test_within_fragment_trimming_to_bounds() {
    //     // Forward: 100..115 (15M)
    //     // Reverse: 110..130 (20M)
    //     // Fragment: 100..130
    //     // Left slice = 100..115 (15 bases), right slice = 110..130 (20 bases)
    //     let (f, r) = mk_pair_basic(
    //         0,
    //         100,
    //         cigar(&[('M', 15)]),
    //         &vec![b'A'; 15],
    //         110,
    //         cigar(&[('M', 20)]),
    //         &vec![b'G'; 20],
    //     );
    //     let frag_seq = FragmentWithSequences::from_pair(&f, &r).expect("frag seq");
    //     assert_eq!(frag_seq.left_seq.len() as u32, 15);
    //     assert_eq!(frag_seq.right_seq.len() as u32, 20);
    //     assert_eq!(frag_seq.frag.start, 100);
    //     assert_eq!(frag_seq.frag.end, 130);
    // }
}

#[cfg(test)]
mod test_segmented_fragments {

    use cfdnalab::utils::fragment::{
        minimal_fragment::oriented_pair_from_read_info, segment_fragment::*,
    };

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
            pos,
            end,
            is_reverse,
            has_ref_gap,
            max_ref_gap,
            ref_mapped_segments: segs.to_vec(),
        }
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
        assert_eq!(fws.start, 10);
        assert_eq!(fws.end, 50);

        // Expect two explicit segments: [10,20] and [40,50]
        let segs = fws.segments.as_ref().expect("segments present");
        assert_eq!(&segs[..], &[(10u32, 20u32), (40u32, 50u32)][..]);
    }

    #[test]
    fn collect_no_gaps_include_inter_mate_gap_plain_span() {
        // forward 10..20, reverse 40..50, include inter-mate gap -> whole span
        let fwd = sri(0, 10, 20, false, false, 0, &[]);
        let rev = sri(0, 40, 50, true, false, 0, &[]);

        let fws = collect_fragment_with_segments(&fwd, &rev, 1, true).expect("pair");
        assert_eq!(fws.start, 10);
        assert_eq!(fws.end, 50);
        assert!(fws.segments.is_none());
    }

    #[test]
    fn collect_with_ref_gap_include_inter_mate_gap_merges_right() {
        // forward has a deletion making two ref-mapped blocks: [10..20], [25..30]
        // reverse is [40..50], include inter-mate gap (30..40)
        let fwd = sri(0, 10, 30, false, true, 5, &[(0, 10), (15, 5)]);
        let rev = sri(0, 40, 50, true, false, 0, &[]);

        let fws = collect_fragment_with_segments(&fwd, &rev, 1, true).expect("pair");
        assert_eq!(fws.start, 10);
        assert_eq!(fws.end, 50);

        // Expected segments after adding gap and merging: [10..20], [25..50]
        let segs = fws.segments.as_ref().expect("segments present");
        assert_eq!(&segs[..], &[(10u32, 20u32), (25u32, 50u32)][..]);
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
            &segs[..],
            &[(10u32, 20u32), (25u32, 30u32), (40u32, 50u32)][..]
        );
    }
}

#[cfg(test)]
mod tests_fragment_with_indel_counts {
    use cfdnalab::utils::fragment::indel_counting_fragment::*;
    use cfdnalab::utils::fragment::minimal_fragment::{
        is_inwards_oriented, oriented_pair_from_read_info,
    };
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

    #[test]
    fn indelreadinfo_parses_deletions_and_insertions() {
        // Forward read: start 100, cigar M50 D5 M45 => deletions at [150,155)
        let r = make_rec(0, 100, false, m_del_m(50, 5, 45));
        let info = IndelReadInfo::from(&r);
        assert_eq!(info.pos, 100);
        assert_eq!(info.end, 100 + 50 + 5 + 45);
        assert_eq!(info.deletions, vec![(150, 155)]);
        assert!(info.insertions.is_empty());

        // Reverse read: start 200, cigar M30 I4 M20 => insertion at ref pos 230 length 4
        let r2 = make_rec(0, 200, true, m_ins_m(30, 4, 20));
        let info2 = IndelReadInfo::from(&r2);
        assert_eq!(info2.insertions, vec![(230, 4)]);
        assert!(info2.deletions.is_empty());
    }

    #[test]
    fn orientation_and_inward_check() {
        // Forward at 100..160, Reverse at 150..210 (inward: forward.pos <= reverse.pos)
        let f = IndelReadInfo::from(&make_rec(0, 100, false, m(60)));
        let r = IndelReadInfo::from(&make_rec(0, 150, true, m(60)));
        let (fwd, rev) = oriented_pair_from_read_info(&f, &r).unwrap();
        assert!(is_inwards_oriented(fwd, rev));
    }

    #[test]
    fn collect_no_indels_fast_path() {
        // No indels; expect zero adjustments.
        let f = make_rec(0, 100, false, m(60));
        let r = make_rec(0, 180, true, m(40));
        let frag = collect_fragment_with_indel_counts_from_records(&f, &r, false, true).unwrap();
        assert_eq!(frag.tid, 0);
        assert_eq!(frag.start, 100);
        assert_eq!(frag.end, 220);
        assert_eq!(frag.len_ref(), 120);
        assert_eq!(frag.deletions_nonoverlap, 0);
        assert_eq!(frag.insertions_nonoverlap, 0);
        assert_eq!(frag.deletions_overlap_supported, 0);
        assert_eq!(frag.insertions_overlap_supported, 0);
        assert_eq!(frag.len_indel_adjusted(), frag.len_ref());
    }

    #[test]
    fn collect_skip_indels_filters_out() {
        // Has insertion; skip_indels=true => None.
        let f = make_rec(0, 100, false, m_ins_m(20, 3, 20));
        let r = make_rec(0, 140, true, m(40));
        assert!(collect_fragment_with_indel_counts_from_records(&f, &r, true, true).is_none());
    }

    #[test]
    fn collect_count_indels_disabled_returns_zeroed() {
        // Indels present but count_indels=false => fragment with zero adjustments.
        let f = make_rec(0, 100, false, m_ins_m(20, 3, 20));
        let r = make_rec(0, 140, true, m_del_m(10, 4, 30));
        let frag = collect_fragment_with_indel_counts_from_records(&f, &r, false, false).unwrap();
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
        let frag = collect_fragment_with_indel_counts_from_records(&f, &r, false, true).unwrap();
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
        // Overlapping mates: forward 100..180, reverse 160..220 -> overlap [160,180)
        // Forward deletion [170,175); Reverse deletion [172,178) -> intersection [172,175) len 3.
        let f = make_rec(0, 100, false, m_del_m(70, 5, 5)); // del at [170,175)
        let r = make_rec(0, 160, true, m_del_m(12, 6, 42)); // del at [172,178)
        let frag = collect_fragment_with_indel_counts_from_records(&f, &r, false, true).unwrap();
        assert_eq!(frag.deletions_nonoverlap, 0);
        assert_eq!(frag.deletions_overlap_supported, 3);
        assert_eq!(frag.insertions_nonoverlap, 0);
        assert_eq!(frag.insertions_overlap_supported, 0);
        // Adjusted length subtracts 3
        assert_eq!(frag.len_indel_adjusted(), frag.len_ref() - 3);
    }

    #[test]
    fn overlap_insertions_require_both_mates_same_ref_pos() {
        // Overlap [160,180).
        // Forward insertion at ref 165 len 5; Reverse insertion at ref 165 len 3 -> min = 3 counted.
        // Forward insertion at ref 170 len 2; Reverse none at 170 -> 0 counted in overlap.
        let f = make_rec(0, 100, false, {
            let mut v = m_ins_m(65, 5, 4); // ins at 165
            v.extend(m_ins_m(1, 2, 10)); // then ins at 170 (since 100 + 65 + [I] + 4 + 1 + [I] + 10)
            v
        });
        let r = make_rec(0, 160, true, m_ins_m(5, 3, 15)); // ins at 165
        let frag = collect_fragment_with_indel_counts_from_records(&f, &r, false, true).unwrap();
        assert_eq!(frag.insertions_nonoverlap, 0); // Second insert is in the overlap but not in both reads
        assert_eq!(frag.insertions_overlap_supported, 3); // min(5,3)
    }

    #[test]
    fn duplicate_insertions_at_same_pos_per_read_take_max_then_min_across_mates() {
        // Overlap [160,200).
        // Forward has two insertions at ref 170: lengths 2 and 5 (separated by soft-clip); keep max=5.
        // Reverse has insertion at ref 170 length 3 -> min(5,3) = 3 counted.
        let f = make_rec(0, 150, false, m_ins_s_ins_m(20, 2, 4, 5, 26)); // two I at ref 170
        let r = make_rec(0, 160, true, m_ins_m(10, 3, 30)); // I at ref 170
        let frag = collect_fragment_with_indel_counts_from_records(&f, &r, false, true).unwrap();
        assert_eq!(frag.insertions_overlap_supported, 3);
        assert_eq!(frag.insertions_nonoverlap, 0);
    }
}

#[cfg(test)]
mod test_kmer_segments {
    use cfdnalab::utils::fragment::segment_kmer_fragment::{
        FragmentWithKmerSegments, KmerSegmentedReadInfo, collect_fragment_with_kmer_segments,
    };
    use cfdnalab::utils::indel_mode::IndelMode;
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
        let f_info = KmerSegmentedReadInfo::from_record(forward, capture_segments);
        let r_info = KmerSegmentedReadInfo::from_record(reverse, capture_segments);
        collect_fragment_with_kmer_segments(&f_info, &r_info, indel_mode, include_gap, end_offset)
    }
    fn segments(frag: &FragmentWithKmerSegments) -> Vec<(u32, u32)> {
        frag.segments.iter().copied().collect()
    }
    #[test]
    fn ignore_mode_without_gap_tracks_per_read_spans() {
        let forward = make_record(0, 100, false, &[Cigar::Match(30)]);
        let reverse = make_record(0, 140, true, &[Cigar::Match(30)]);
        let frag = collect_pair(&forward, &reverse, IndelMode::Ignore, false, 0).expect("fragment");
        assert_eq!(frag.start, 100);
        assert_eq!(frag.end, 170);
        assert_eq!(segments(&frag), vec![(100, 130), (140, 170)]);
    }
    #[test]
    fn ignore_mode_includes_gap_when_requested() {
        let forward = make_record(0, 100, false, &[Cigar::Match(30)]);
        let reverse = make_record(0, 140, true, &[Cigar::Match(30)]);
        let frag = collect_pair(&forward, &reverse, IndelMode::Ignore, true, 0).expect("fragment");
        assert_eq!(frag.start, 100);
        assert_eq!(frag.end, 170);
        assert_eq!(segments(&frag), vec![(100, 170)]);
    }
    #[test]
    fn skip_mode_filters_pairs_with_indels() {
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
        let forward = make_record(
            0,
            100,
            false,
            &[Cigar::Match(10), Cigar::Ins(2), Cigar::Match(10)],
        );
        let reverse = make_record(0, 150, true, &[Cigar::Match(15)]);
        let frag = collect_pair(&forward, &reverse, IndelMode::Adjust, false, 0).expect("fragment");
        assert_eq!(frag.start, 100);
        assert_eq!(frag.end, 165);
        assert_eq!(segments(&frag), vec![(100, 110), (110, 120), (150, 165)]);
    }
    #[test]
    fn inter_mate_gap_not_merged_when_border_has_insertion() {
        let forward = make_record(0, 100, false, &[Cigar::Match(20), Cigar::Ins(2)]);
        let reverse = make_record(0, 130, true, &[Cigar::Match(20)]);
        let frag = collect_pair(&forward, &reverse, IndelMode::Adjust, true, 0).expect("fragment");
        assert_eq!(frag.start, 100);
        assert_eq!(frag.end, 150);
        // Not merged at the left boundary but merged on the right
        assert_eq!(segments(&frag), vec![(100, 120), (120, 150)]);
    }
    #[test]
    fn end_offset_trims_segments_but_preserves_span() {
        let forward = make_record(0, 100, false, &[Cigar::Match(40)]);
        let reverse = make_record(0, 120, true, &[Cigar::Match(40)]);
        let frag = collect_pair(&forward, &reverse, IndelMode::Ignore, true, 5).expect("fragment");
        assert_eq!(frag.start, 100);
        assert_eq!(frag.end, 160);
        assert_eq!(segments(&frag), vec![(105, 155)]);
    }

    #[test]
    fn adjust_mode_handles_mixed_insertions_and_deletions() {
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

        assert_eq!(frag.start, 100);
        assert_eq!(frag.end, 151);
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
        let forward = make_record(0, 100, false, &[Cigar::Match(15), Cigar::Ins(2)]);
        let reverse = make_record(0, 140, true, &[Cigar::Ins(3), Cigar::Match(18)]);

        let frag = collect_pair(&forward, &reverse, IndelMode::Adjust, true, 0).expect("fragment");

        assert_eq!(frag.start, 100);
        assert_eq!(frag.end, 158);
        assert_eq!(segments(&frag), vec![(100, 115), (115, 140), (140, 158)]);
    }

    #[test]
    fn returns_none_for_non_inward_orientation() {
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
        let forward = make_record(0, 100, false, &[Cigar::Match(20)]);
        let reverse = make_record(0, 120, true, &[Cigar::Match(20)]);

        assert!(collect_pair(&forward, &reverse, IndelMode::Ignore, true, 20).is_none());
    }

    #[test]
    fn overlap_consensus_indel_behaviour() {
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

        assert_eq!(frag.start, 100);
        assert_eq!(frag.end, 160);
        assert_eq!(
            segments(&frag),
            vec![(100, 110), (110, 125), (125, 130), (132, 145), (145, 160)]
        );

        let frag_asked_gap =
            collect_pair(&forward, &reverse, IndelMode::Adjust, true, 0).expect("fragment");
        assert_eq!(segments(&frag), segments(&frag_asked_gap));

        let frag_offset =
            collect_pair(&forward, &reverse, IndelMode::Adjust, false, 12).expect("fragment");

        assert_eq!(frag_offset.start, 100);
        assert_eq!(frag_offset.end, 160);
        assert_eq!(
            segments(&frag_offset),
            vec![(112, 125), (125, 130), (132, 145), (145, 148)]
        );
    }
}
