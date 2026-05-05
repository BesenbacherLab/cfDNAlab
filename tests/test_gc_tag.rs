mod test_gc_tag_values {
    #[cfg(feature = "cmd_gc_bias")]
    use cfdnalab::commands::gc_bias::{
        correct::GCCorrector, counting::build_gc_prefixes, package::GCCorrectionPackage,
    };
    use rust_htslib::bam::record::{Aux, Record};

    use cfdnalab::shared::base::ZEROISH_F32_TOLERANCE;
    #[cfg(feature = "cmd_gc_bias")]
    use cfdnalab::shared::constants::GC_CORRECTION_SCHEMA_VERSION;
    use cfdnalab::shared::gc_tag::{
        ClassifiedGCTagWeight, GCTagValue, MIN_REASONABLE_GC_WEIGHT, combine_gc_tag_values,
        read_gc_tag_from_record,
    };
    #[cfg(feature = "cmd_gc_bias")]
    use cfdnalab::shared::interval::Interval;
    #[cfg(feature = "cmd_gc_bias")]
    use ndarray::array;

    #[test]
    fn gc_tag_values_follow_supported_range_and_zero_snap_rules() {
        // Human verification status: unverified
        // Arrange: start with a sane weight
        let mut rec_ok = Record::new();
        rec_ok.push_aux(b"GC", Aux::Float(2.5)).expect("set GC tag");
        let ok = read_gc_tag_from_record(&rec_ok, b"GC");

        // Assert: valid weight passes through
        assert_eq!(ok.weight, Some(2.5));
        assert!(!ok.was_missing);
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

        // Arrange: meaningfully negative values are invalid, not zero-snapped.
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

        // Arrange: tiny positive values near zero are snapped to zero.
        let mut rec_tiny = Record::new();
        rec_tiny
            .push_aux(b"GC", Aux::Float(ZEROISH_F32_TOLERANCE))
            .expect("set GC tag");
        let tiny = read_gc_tag_from_record(&rec_tiny, b"GC");
        assert_eq!(tiny.weight, Some(0.0));
        assert!(!tiny.had_invalid);
        assert!(!tiny.was_out_of_range);

        // Arrange: the zero-snap window is symmetric around zero.
        let mut rec_tiny_negative = Record::new();
        rec_tiny_negative
            .push_aux(b"GC", Aux::Float(-ZEROISH_F32_TOLERANCE))
            .expect("set GC tag");
        let tiny_negative = read_gc_tag_from_record(&rec_tiny_negative, b"GC");
        assert_eq!(tiny_negative.weight, Some(0.0));
        assert!(!tiny_negative.had_invalid);
        assert!(!tiny_negative.was_out_of_range);

        // Arrange: positive values below the minimum supported GC weight are invalid.
        let mut rec_low = Record::new();
        rec_low
            .push_aux(b"GC", Aux::Float(MIN_REASONABLE_GC_WEIGHT / 10.0))
            .expect("set GC tag");
        let low = read_gc_tag_from_record(&rec_low, b"GC");
        assert!(low.weight.is_none());
        assert!(low.had_invalid);
        assert!(low.was_out_of_range);
    }

    #[test]
    fn gc_tag_values_just_below_minimum_supported_weight_are_invalid() {
        // Human verification status: unverified
        // Arrange: choose the nearest representable f32 below the supported lower bound.
        // This is a stronger boundary than "/10" because it proves the exact cutoff behavior.
        let just_below_min = f32::from_bits(MIN_REASONABLE_GC_WEIGHT.to_bits() - 1);
        let mut rec = Record::new();
        rec.push_aux(b"GC", Aux::Float(just_below_min))
            .expect("set GC tag");

        // Act
        let observed = read_gc_tag_from_record(&rec, b"GC");

        // Assert: values below 1e-3 remain invalid even when they are only one f32 step lower.
        assert!(observed.weight.is_none());
        assert!(observed.had_invalid);
        assert!(observed.was_out_of_range);
    }

    #[test]
    fn missing_gc_tag_is_reported_separately() {
        // Human verification status: unverified
        let rec = Record::new();
        let missing = read_gc_tag_from_record(&rec, b"GC");
        assert!(missing.weight.is_none());
        assert!(missing.was_missing);
        assert!(!missing.had_invalid);
        assert!(!missing.was_out_of_range);
    }

    #[test]
    fn combining_valid_weights_averages_before_final_range_check() {
        // Human verification status: unverified
        let mut rec_a = Record::new();
        rec_a.push_aux(b"GC", Aux::Float(2.0)).expect("set GC tag");
        let mut rec_b = Record::new();
        rec_b.push_aux(b"GC", Aux::Float(4.0)).expect("set GC tag");

        let a = read_gc_tag_from_record(&rec_a, b"GC");
        let b = read_gc_tag_from_record(&rec_b, b"GC");
        let combined = combine_gc_tag_values(&a, &b);

        assert_eq!(combined.weight, Some(3.0));
        assert!(!combined.had_invalid);
        assert!(!combined.was_out_of_range);
    }

    #[test]
    fn combining_paired_tags_reuses_single_usable_mate_and_keeps_zero_precedence() {
        // Human verification status: unverified
        let mut rec_zero = Record::new();
        rec_zero
            .push_aux(b"GC", Aux::Float(0.0))
            .expect("set GC tag");
        let mut rec_valid = Record::new();
        rec_valid
            .push_aux(b"GC", Aux::Float(4.0))
            .expect("set GC tag");

        let zero = read_gc_tag_from_record(&rec_zero, b"GC");
        let valid = read_gc_tag_from_record(&rec_valid, b"GC");
        let zero_combined = combine_gc_tag_values(&zero, &valid);
        assert_eq!(zero_combined.weight, Some(0.0));
        assert!(!zero_combined.had_invalid);

        let missing = read_gc_tag_from_record(&Record::new(), b"GC");
        let missing_combined = combine_gc_tag_values(&valid, &missing);
        assert_eq!(missing_combined.weight, Some(4.0));
        assert!(!missing_combined.was_missing);
        assert!(!missing_combined.had_invalid);
        assert!(!missing_combined.was_out_of_range);

        let mut rec_low = Record::new();
        rec_low
            .push_aux(b"GC", Aux::Float(MIN_REASONABLE_GC_WEIGHT / 10.0))
            .expect("set GC tag");
        let low = read_gc_tag_from_record(&rec_low, b"GC");
        let invalid_combined = combine_gc_tag_values(&valid, &low);
        assert!(invalid_combined.weight.is_none());
        assert!(invalid_combined.had_invalid);
        assert!(invalid_combined.was_out_of_range);
    }

    #[test]
    fn gc_tag_classify_exposes_one_explicit_state() {
        // Human verification status: unverified
        assert_eq!(
            GCTagValue {
                weight: Some(2.5),
                was_missing: false,
                had_invalid: false,
                was_out_of_range: false,
            }
            .classify()
            .expect("valid classification"),
            ClassifiedGCTagWeight::Usable(2.5)
        );
        assert_eq!(
            GCTagValue::missing()
                .classify()
                .expect("missing classification"),
            ClassifiedGCTagWeight::Missing
        );
        assert_eq!(
            GCTagValue {
                weight: None,
                was_missing: false,
                had_invalid: true,
                was_out_of_range: true,
            }
            .classify()
            .expect("invalid classification"),
            ClassifiedGCTagWeight::Invalid { out_of_range: true }
        );
    }

    #[test]
    fn gc_tag_classify_rejects_inconsistent_internal_state() {
        // Human verification status: unverified
        let err = GCTagValue {
            weight: None,
            was_missing: false,
            had_invalid: false,
            was_out_of_range: false,
        }
        .classify()
        .expect_err("inconsistent state should error");

        assert!(
            err.to_string().contains("inconsistent GC tag state"),
            "unexpected error: {err}"
        );
    }

    #[cfg(feature = "cmd_gc_bias")]
    #[test]
    fn gc_file_weights_follow_the_same_sanity_rules() {
        // Human verification status: unverified
        let prefixes = build_gc_prefixes(b"AAAAAAAAAA");
        let interval = Interval::new(0_u64, 10_u64).expect("valid interval");
        let scenarios = [
            ("negative_below_snap_window_is_unusable", -3.0_f64, None),
            (
                "tiny_negative_becomes_zero",
                -(ZEROISH_F32_TOLERANCE as f64),
                Some(0.0_f64),
            ),
            (
                "tiny_positive_becomes_zero",
                ZEROISH_F32_TOLERANCE as f64,
                Some(0.0_f64),
            ),
            (
                "too_small_positive_is_unusable",
                (MIN_REASONABLE_GC_WEIGHT / 10.0) as f64,
                None,
            ),
            ("too_large_positive_is_unusable", 1.1e3_f64, None),
        ];

        for (name, weight, expected) in scenarios {
            let package = GCCorrectionPackage {
                version: GC_CORRECTION_SCHEMA_VERSION,
                end_offset: 0,
                length_edges: vec![10, 11],
                gc_edges: vec![0, 101],
                length_bin_frequencies: array![1.0_f64],
                reference_contig_footprint: Vec::new(),
                correction_matrix: array![[weight]],
            };
            let corrector = GCCorrector::from_package(&package).expect("build corrector");

            let observed = corrector
                .correct_fragment(interval, &prefixes)
                .expect("correct fragment");
            assert_eq!(observed, expected, "unexpected sanitized weight for {name}");
        }
    }

    #[cfg(feature = "cmd_gc_bias")]
    #[test]
    fn gc_file_weights_just_below_minimum_supported_weight_are_invalid() {
        // Human verification status: unverified
        // Arrange: use the nearest representable f64 below the exact f64 threshold used by the
        // sanitizer so this checks the boundary itself, not just a clearly too-small value.
        let just_below_min = f64::from_bits((MIN_REASONABLE_GC_WEIGHT as f64).to_bits() - 1);
        let prefixes = build_gc_prefixes(b"AAAAAAAAAA");
        let interval = Interval::new(0_u64, 10_u64).expect("valid interval");
        let package = GCCorrectionPackage {
            version: GC_CORRECTION_SCHEMA_VERSION,
            end_offset: 0,
            length_edges: vec![10, 11],
            gc_edges: vec![0, 101],
            length_bin_frequencies: array![1.0_f64],
            reference_contig_footprint: Vec::new(),
            correction_matrix: array![[just_below_min]],
        };
        let corrector = GCCorrector::from_package(&package).expect("build corrector");

        // Act
        let observed = corrector
            .correct_fragment(interval, &prefixes)
            .expect("correct fragment");

        // Assert: the GC-file path rejects values that fall just below the accepted range.
        assert_eq!(observed, None);
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
            gc_tag::GCTagValue,
            indel_mode::{IndelMode, IndelMotifFilterPolicy},
        },
    };
    use rust_htslib::bam::record::{Aux, Cigar, CigarString, Record};

    fn assert_valid_gc_tag(observed: GCTagValue, expected_weight: f32) {
        assert_eq!(observed.weight, Some(expected_weight));
        assert!(!observed.was_missing);
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
            u32::MAX,
            &[],
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
            u32::MAX,
            &[],
            Some(b"GC"),
            |_fragment: &FragmentWithEnds| true,
            true,
        ));

        // Assert
        assert_valid_gc_tag(fragment.gc_tag, 2.5);
    }
}
