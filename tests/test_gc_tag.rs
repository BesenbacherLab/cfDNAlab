mod test_gc_tag_values {
    use rust_htslib::bam::record::{Aux, Record};

    use cfdnalab::shared::gc_tag::read_gc_tag_from_record;

    #[test]
    fn should_reject_extreme_or_invalid_gc_weights() {
        // Arrange: start with a sane weight
        let mut rec_ok = Record::new();
        rec_ok.push_aux(b"GC", Aux::Float(2.5)).expect("set GC tag");
        let ok = read_gc_tag_from_record(&rec_ok, b"GC");

        // Assert: valid weight passes through
        assert_eq!(ok.weight, Some(2.5));
        assert!(!ok.had_invalid);
        assert!(!ok.was_out_of_range);

        // Arrange: record carrying a wildly high weight that should be treated as invalid
        let mut rec_high = Record::new();
        rec_high
            .push_aux(b"GC", Aux::Float(1.1e3))
            .expect("set GC tag");
        let high = read_gc_tag_from_record(&rec_high, b"GC");

        // Assert: extreme values are rejected to avoid runaway coverage
        assert!(high.weight.is_none());
        assert!(high.had_invalid);
        assert!(high.was_out_of_range);

        // Arrange: negative weights are also invalid
        let mut rec_neg = Record::new();
        rec_neg
            .push_aux(b"GC", Aux::Float(-3.0))
            .expect("set GC tag");
        let neg = read_gc_tag_from_record(&rec_neg, b"GC");

        assert!(neg.weight.is_none());
        assert!(neg.had_invalid);
        assert!(neg.was_out_of_range);

        // Arrange: NaN should be invalid but not counted as out-of-range
        let mut rec_nan = Record::new();
        rec_nan
            .push_aux(b"GC", Aux::Float(f32::NAN))
            .expect("set GC tag");
        let nan = read_gc_tag_from_record(&rec_nan, b"GC");

        assert!(nan.weight.is_none());
        assert!(nan.had_invalid);
        assert!(!nan.was_out_of_range);
    }

    #[test]
    fn combining_valid_weights_remains_in_range() {
        let mut rec_a = Record::new();
        rec_a.push_aux(b"GC", Aux::Float(2.0)).expect("set GC tag");
        let mut rec_b = Record::new();
        rec_b.push_aux(b"GC", Aux::Float(4.0)).expect("set GC tag");

        let a = read_gc_tag_from_record(&rec_a, b"GC");
        let b = read_gc_tag_from_record(&rec_b, b"GC");
        let combined = cfdnalab::shared::gc_tag::combine_gc_tag_values(&a, &b);

        assert_eq!(combined.weight, Some(3.0));
        assert!(!combined.had_invalid);
        assert!(!combined.was_out_of_range);
    }
}
