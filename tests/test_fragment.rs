#[cfg(test)]
mod tests {
    use cfdna_utils::cfdna_utils::fragment::*;
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
    /// NOTE: Depending on your rust-htslib version, you may need to tweak the `set` calls.
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

    #[test]
    fn test_overlap_mm_basic_match() {
        // Forward: 100..160 (60M); reverse: 140..200 (60M)
        // Overlap window: 140..160 (20 bp), all 'A'
        let (f, r) = mk_pair_basic(
            0,
            100,
            cigar(&[('M', 60)]),
            &vec![b'C'; 60], // Only positions 40..60 will contribute to overlap
            140,
            cigar(&[('M', 60)]),
            &vec![b'A'; 60], // Positions 0..20 -> overlap
        );

        let ov = FragmentOverlapMM::from_pair(&f, &r).expect("overlap mm");
        assert_eq!(ov.overlap.start, 140);
        assert_eq!(ov.overlap.end, 160);
        assert_eq!(ov.ref_coords.len(), 20);
        assert_eq!(ov.left_bases.len(), 20);
        assert_eq!(ov.right_bases.len(), 20);

        // Coordinates should be contiguous 140..160
        for (i, &coord) in ov.ref_coords.iter().enumerate() {
            assert_eq!(coord, 140 + i as u32);
        }
        // Compare that we pulled aligned columns; we didn't set exact sequences for left overlapping
        // region; set them now to 'A' to ensure a match:
        // Forward read bases 40..60 must be 'A' for equality; rebuild with that:
        let (f2, r2) = mk_pair_basic(
            0,
            100,
            cigar(&[('M', 40), ('M', 20)]),
            &[vec![b'C'; 40], vec![b'A'; 20]].concat(),
            140,
            cigar(&[('M', 60)]),
            &vec![b'A'; 60],
        );
        let ov2 = FragmentOverlapMM::from_pair(&f2, &r2).unwrap();
        assert!(
            ov2.left_bases
                .iter()
                .zip(&ov2.right_bases)
                .all(|(l, r)| l == r)
        );
        assert!(ov2.left_bases.iter().all(|&b| b == b'A'));
    }

    #[test]
    fn test_overlap_mm_drops_insertions() {
        // Forward: 100..120 on ref with an insertion at 110: 10M 2I 10M
        // Reverse: 100..122 on ref: 22M (covers entire region)
        // MM extractor should DROP the insertion bases and align by ref coords (20 positions).
        let (f, r) = mk_pair_basic(
            0,
            100,
            cigar(&[('M', 10), ('I', 2), ('M', 10)]),
            // Seq length = 10 + 2 + 10 = 22
            &[vec![b'G'; 10], vec![b'T'; 2], vec![b'G'; 10]].concat(),
            100,
            cigar(&[('M', 22)]),
            &vec![b'G'; 22],
        );
        let ov = FragmentOverlapMM::from_pair(&f, &r).expect("overlap mm");
        assert_eq!(ov.overlap.start, 100);
        assert_eq!(ov.overlap.end, 120);
        assert_eq!(ov.ref_coords.len(), 20); // Insertion dropped (no ref coord)
        // Check exact reference coordinates are 100..120 (no duplicates at the insertion boundary)
        let expected: Vec<u32> = (100..120).collect();
        assert_eq!(ov.ref_coords, expected);
        assert!(ov.left_bases.iter().all(|&b| b == b'G'));
        assert!(ov.right_bases.iter().all(|&b| b == b'G'));
    }

    #[test]
    fn test_overlap_mm_none_when_no_overlap() {
        // Forward: 100..120; reverse: 130..150 => no overlap
        let (f, r) = mk_pair_basic(
            0,
            100,
            cigar(&[('M', 20)]),
            &vec![b'A'; 20],
            130,
            cigar(&[('M', 20)]),
            &vec![b'A'; 20],
        );
        assert!(FragmentOverlapMM::from_pair(&f, &r).is_none());
    }

    #[test]
    fn test_overlap_mm_skips_deletion_positions() {
        // Forward: 100..122 on ref with a deletion at 110..112: 10M 2D 10M
        // Reverse: 100..122 on ref: 22M
        // MM extractor should DROP the deletion columns (no read base) and keep only aligned coords.
        let (f, r) = mk_pair_basic(
            0,
            100,
            cigar(&[('M', 10), ('D', 2), ('M', 10)]),
            // seq length = 10 + 10 = 20 (D consumes reference only)
            &[vec![b'G'; 10], vec![b'G'; 10]].concat(),
            100,
            cigar(&[('M', 22)]),
            &[vec![b'G'; 10], vec![b'A'; 2], vec![b'G'; 10]].concat(),
        );
        let ov = FragmentOverlapMM::from_pair(&f, &r).expect("overlap mm");
        assert_eq!(ov.overlap.start, 100);
        assert_eq!(ov.overlap.end, 122);

        // Expect reference coords 100..110 and 112..122 (skip 110,111 where the deletion is)
        let mut expected: Vec<u32> = (100..110).collect();
        expected.extend(112..122);
        assert_eq!(ov.ref_coords, expected);

        // Bases should align and match
        assert_eq!(ov.left_bases.len(), expected.len());
        assert_eq!(ov.right_bases.len(), expected.len());
        assert!(ov.left_bases.iter().all(|&b| b == b'G'));
        assert!(ov.right_bases.iter().all(|&b| b == b'G'));
    }

    // Tests for `FragmentWithSequences` (within-fragment sequences)

    #[test]
    fn test_fragment_with_sequences_includes_insertions_and_flags() {
        // Fragment: forward 100..120 (20M), reverse 105..130 (25M) => fragment 100..130
        // Forward read contains insertion at ref 110: 10M 2I 10M
        // Reverse is plain 25M
        let (f, r) = mk_pair_basic(
            0,
            100,
            cigar(&[('M', 10), ('I', 2), ('M', 10)]),
            &[vec![b'A'; 10], vec![b'G'; 2], vec![b'A'; 10]].concat(),
            105,
            cigar(&[('M', 25)]),
            &vec![b'C'; 25],
        );
        let frag_seq = FragmentWithSequences::from_pair(&f, &r).expect("frag seq");
        // Left sequence should include the insertion (2 bases)
        assert_eq!(frag_seq.left_seq.len(), 22);
        assert!(frag_seq.left_info.has_insertion);
        assert!(!frag_seq.left_info.has_deletion);
        assert!(!frag_seq.left_info.has_refskip);
    }

    #[test]
    fn test_fragment_with_sequences_flags_deletion() {
        // Forward has a deletion inside the fragment: 10M 2D 10M
        // Reverse covers the same span with 22M
        let (f, r) = mk_pair_basic(
            0,
            100,
            cigar(&[('M', 10), ('D', 2), ('M', 10)]),
            &[vec![b'A'; 10], vec![b'A'; 10]].concat(),
            100,
            cigar(&[('M', 22)]),
            &vec![b'C'; 22],
        );
        let frag_seq = FragmentWithSequences::from_pair(&f, &r).expect("frag seq");

        // Fragment spans full reference window
        assert_eq!(frag_seq.frag.start, 100);
        assert_eq!(frag_seq.frag.end, 122);

        // Left sequence length excludes deleted reference columns (still 20 bases emitted)
        assert_eq!(frag_seq.left_seq.len(), 20);
        assert_eq!(frag_seq.right_seq.len(), 22);
        // Flags: deletion seen on left; no insertion/softclip
        assert!(frag_seq.left_info.has_deletion);
        assert!(!frag_seq.left_info.has_insertion);
        assert!(!frag_seq.left_info.has_softclip);
        // Right read has no deletion
        assert!(!frag_seq.right_info.has_deletion);
    }

    #[test]
    fn test_fragment_with_sequences_softclips_excluded_but_flagged() {
        // Forward: 5S 20M at pos=100 -> mapped span 100..120; soft clip before POS
        // Reverse: 120..145 (25M) => fragment 100..145
        let (f, r) = mk_pair_basic(
            0,
            100,
            cigar(&[('S', 5), ('M', 20)]),
            // seq length 25; first 5 are soft-clipped and must NOT appear in left_seq
            &[vec![b'T'; 5], vec![b'A'; 20]].concat(),
            120,
            cigar(&[('M', 25)]),
            &vec![b'C'; 25],
        );
        let frag_seq = FragmentWithSequences::from_pair(&f, &r).expect("frag seq");
        // The left read slice inside the fragment should have only the 20 aligned A's.
        assert_eq!(frag_seq.left_seq, vec![b'A'; 20]);
        assert!(frag_seq.left_info.has_softclip);
        // Right has no soft clips
        assert!(!frag_seq.right_info.has_softclip);
    }

    #[test]
    fn test_within_fragment_trimming_to_bounds() {
        // Forward: 100..115 (15M)
        // Reverse: 110..130 (20M)
        // Fragment: 100..130
        // Left slice = 100..115 (15 bases), right slice = 110..130 (20 bases)
        let (f, r) = mk_pair_basic(
            0,
            100,
            cigar(&[('M', 15)]),
            &vec![b'A'; 15],
            110,
            cigar(&[('M', 20)]),
            &vec![b'G'; 20],
        );
        let frag_seq = FragmentWithSequences::from_pair(&f, &r).expect("frag seq");
        assert_eq!(frag_seq.left_seq.len() as u32, 15);
        assert_eq!(frag_seq.right_seq.len() as u32, 20);
        assert_eq!(frag_seq.frag.start, 100);
        assert_eq!(frag_seq.frag.end, 130);
    }
}
