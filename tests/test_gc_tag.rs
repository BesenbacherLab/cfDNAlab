mod test_gc_tag_values {
    use rust_htslib::bam::record::{Aux, Record};

    use cfdnalab::shared::gc_tag::read_gc_tag_from_record;

    #[test]
    fn should_reject_extreme_or_invalid_gc_weights() {
        // Arrange: start with a sane weight
        let mut rec_ok = Record::new();
        rec_ok
            .push_aux(b"GC", &Aux::Float(2.5))
            .expect("set GC tag");
        let ok = read_gc_tag_from_record(&rec_ok, b"GC");

        // Assert: valid weight passes through
        assert_eq!(ok.weight, Some(2.5));
        assert!(!ok.had_invalid);

        // Arrange: record carrying a wildly high weight that should be treated as invalid
        let mut rec_high = Record::new();
        rec_high.push_aux(b"GC", &Aux::Float(1.1e3)).expect("set GC tag");
        let high = read_gc_tag_from_record(&rec_high, b"GC");

        // Assert: extreme values are rejected to avoid runaway coverage
        assert!(high.weight.is_none());
        assert!(high.had_invalid);

        // Arrange: negative weights are also invalid
        let mut rec_neg = Record::new();
        rec_neg
            .push_aux(b"GC", &Aux::Float(-3.0))
            .expect("set GC tag");
        let neg = read_gc_tag_from_record(&rec_neg, b"GC");

        assert!(neg.weight.is_none());
        assert!(neg.had_invalid);
    }
}
