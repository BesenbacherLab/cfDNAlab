mod test_frag_file_fragment {
    use crate::shared::fragment::frag_file_fragment::{
        FragReadInfo, collect_fragment_with_frag_file_info,
    };
    use crate::shared::interval::Interval;

    #[test]
    fn collect_frag_file_fragment_uses_directional_pair_span_not_interval_union() {
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
