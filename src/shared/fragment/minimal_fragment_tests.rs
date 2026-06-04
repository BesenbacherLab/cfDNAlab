#[cfg(test)]
mod test_minimal_fragment {
    use crate::shared::fragment::minimal_fragment::*;
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
