mod test_gc_tag_values {
    use rust_htslib::bam::record::{Aux, Record};

    use cfdnalab::shared::gc_tag::read_gc_tag_from_record;

    #[test]
    fn should_reject_extreme_or_invalid_gc_weights() {
        // Human verification status: unverified
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
        // Human verification status: unverified
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

#[cfg(feature = "cmd_ends")]
mod test_fragment_iterator_gc_tags {
    use anyhow::Result;
    use cfdnalab::{
        commands::ends::config_structs::{ClipStrategy, KmerSource},
        shared::{
            fragment::{
                ends_fragment::FragmentWithEnds, minimal_fragment::Fragment,
                segment_fragment::FragmentWithSegments,
                segment_kmer_fragment::FragmentWithKmerSegments,
            },
            fragment_iterators::{
                fragments_from_bam, fragments_with_ends_from_bam,
                fragments_with_kmer_segments_from_bam, fragments_with_segments_from_bam,
            },
            gc_tag::GcTagValue,
            indel_mode::{IndelMode, IndelMotifFilterPolicy},
        },
    };
    use rust_htslib::bam::record::{Aux, Cigar, CigarString, Record};

    fn assert_valid_gc_tag(observed: GcTagValue, expected_weight: f32) {
        assert_eq!(observed.weight, Some(expected_weight));
        assert!(!observed.had_invalid);
        assert!(!observed.was_out_of_range);
    }

    fn make_record(
        qname: &[u8],
        tid: i32,
        pos: i64,
        is_reverse: bool,
        seq_len: usize,
        gc_weight: f32,
    ) -> Record {
        let mut record = Record::new();
        record.set_tid(tid);
        record.set_pos(pos);
        record.set_flags(if is_reverse { 0x11 } else { 0x1 });
        record.set_mapq(60);

        let cigar = CigarString(vec![Cigar::Match(seq_len as u32)]);
        let seq = vec![b'A'; seq_len];
        let qual = vec![30u8; seq_len];
        record.set(qname, Some(&cigar), &seq, &qual);
        record
            .push_aux(b"GC", Aux::Float(gc_weight))
            .expect("set GC tag");

        record
    }

    fn first_fragment(iter: impl Iterator<Item = Result<Fragment>>) -> Fragment {
        iter.into_iter()
            .next()
            .expect("one fragment")
            .expect("valid fragment")
    }

    fn first_segment_fragment(
        iter: impl Iterator<Item = Result<FragmentWithSegments>>,
    ) -> FragmentWithSegments {
        iter.into_iter()
            .next()
            .expect("one fragment")
            .expect("valid fragment")
    }

    fn first_kmer_segment_fragment(
        iter: impl Iterator<Item = Result<FragmentWithKmerSegments>>,
    ) -> FragmentWithKmerSegments {
        iter.into_iter()
            .next()
            .expect("one fragment")
            .expect("valid fragment")
    }

    fn first_end_fragment(
        iter: impl Iterator<Item = Result<FragmentWithEnds>>,
    ) -> FragmentWithEnds {
        iter.into_iter()
            .next()
            .expect("one fragment")
            .expect("valid fragment")
    }

    #[test]
    fn basic_fragment_iterator_paired_uses_configured_gc_tag() {
        // Arrange: two mates with GC weights 2 and 4 should average to 3 on the fragment.
        let qname = b"pair_basic";
        let forward = make_record(qname, 0, 100, false, 50, 2.0);
        let reverse = make_record(qname, 0, 150, true, 50, 4.0);

        // Act: build fragments through the same iterator used by basic-fragment commands.
        let fragment = first_fragment(fragments_from_bam(
            vec![Ok(forward), Ok(reverse)].into_iter(),
            |_record| true,
            Some(b"GC"),
            |_fragment: &Fragment| true,
            false,
        ));

        // Assert: the configured GC tag is preserved and combined at fragment level.
        assert_valid_gc_tag(fragment.gc_tag, 3.0);
    }

    #[test]
    fn basic_fragment_iterator_unpaired_uses_configured_gc_tag() {
        // Arrange: a single read-as-fragment should keep its own GC-tag value.
        let record = make_record(b"single_basic", 0, 100, false, 50, 2.5);

        // Act
        let fragment = first_fragment(fragments_from_bam(
            vec![Ok(record)].into_iter(),
            |_record| true,
            Some(b"GC"),
            |_fragment: &Fragment| true,
            true,
        ));

        // Assert
        assert_valid_gc_tag(fragment.gc_tag, 2.5);
    }

    #[test]
    fn segment_fragment_iterator_paired_uses_configured_gc_tag() {
        // Arrange
        let qname = b"pair_segments";
        let forward = make_record(qname, 0, 100, false, 50, 2.0);
        let reverse = make_record(qname, 0, 150, true, 50, 4.0);

        // Act
        let fragment = first_segment_fragment(fragments_with_segments_from_bam(
            vec![Ok(forward), Ok(reverse)].into_iter(),
            |_record| true,
            1,
            true,
            Some(b"GC"),
            |_fragment: &FragmentWithSegments| true,
            false,
        ));

        // Assert
        assert_valid_gc_tag(fragment.gc_tag, 3.0);
    }

    #[test]
    fn segment_fragment_iterator_unpaired_uses_configured_gc_tag() {
        // Arrange
        let record = make_record(b"single_segments", 0, 100, false, 50, 2.5);

        // Act
        let fragment = first_segment_fragment(fragments_with_segments_from_bam(
            vec![Ok(record)].into_iter(),
            |_record| true,
            1,
            true,
            Some(b"GC"),
            |_fragment: &FragmentWithSegments| true,
            true,
        ));

        // Assert
        assert_valid_gc_tag(fragment.gc_tag, 2.5);
    }

    #[test]
    fn kmer_segment_iterator_paired_uses_configured_gc_tag() {
        // Arrange
        let qname = b"pair_kmers";
        let forward = make_record(qname, 0, 100, false, 50, 2.0);
        let reverse = make_record(qname, 0, 150, true, 50, 4.0);

        // Act
        let fragment = first_kmer_segment_fragment(fragments_with_kmer_segments_from_bam(
            vec![Ok(forward), Ok(reverse)].into_iter(),
            |_record| true,
            IndelMode::Ignore,
            true,
            0,
            Some(b"GC"),
            |_fragment: &FragmentWithKmerSegments| true,
            false,
        ));

        // Assert
        assert_valid_gc_tag(fragment.gc_tag, 3.0);
    }

    #[test]
    fn kmer_segment_iterator_unpaired_uses_configured_gc_tag() {
        // Arrange
        let record = make_record(b"single_kmers", 0, 100, false, 50, 2.5);

        // Act
        let fragment = first_kmer_segment_fragment(fragments_with_kmer_segments_from_bam(
            vec![Ok(record)].into_iter(),
            |_record| true,
            IndelMode::Ignore,
            true,
            0,
            Some(b"GC"),
            |_fragment: &FragmentWithKmerSegments| true,
            true,
        ));

        // Assert
        assert_valid_gc_tag(fragment.gc_tag, 2.5);
    }

    #[test]
    fn ends_iterator_paired_uses_configured_gc_tag() {
        // Arrange
        let qname = b"pair_ends";
        let forward = make_record(qname, 0, 100, false, 50, 2.0);
        let reverse = make_record(qname, 0, 150, true, 50, 4.0);

        // Act
        let fragment = first_end_fragment(fragments_with_ends_from_bam(
            vec![Ok(forward), Ok(reverse)].into_iter(),
            |_record| true,
            ClipStrategy::Aligned,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            None,
            Some(b"GC"),
            |_fragment: &FragmentWithEnds| true,
            false,
        ));

        // Assert
        assert_valid_gc_tag(fragment.gc_tag, 3.0);
    }

    #[test]
    fn ends_iterator_unpaired_uses_configured_gc_tag() {
        // Arrange
        let record = make_record(b"single_ends", 0, 100, false, 50, 2.5);

        // Act
        let fragment = first_end_fragment(fragments_with_ends_from_bam(
            vec![Ok(record)].into_iter(),
            |_record| true,
            ClipStrategy::Aligned,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            None,
            Some(b"GC"),
            |_fragment: &FragmentWithEnds| true,
            true,
        ));

        // Assert
        assert_valid_gc_tag(fragment.gc_tag, 2.5);
    }
}
