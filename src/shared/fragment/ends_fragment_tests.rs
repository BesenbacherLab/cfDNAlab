mod test_fragment_with_ends_1 {
    use super::super::*;
    use crate::commands::ends::config_structs::BaseQualityComparisonOp;
    use crate::shared::gc_tag::GCTagValue;
    use rust_htslib::bam::record::{Cigar, CigarString, Record};

    fn base_read_info(qualities: &[u8]) -> EndReadInfo {
        base_read_info_with_clips(qualities, 0, 0)
    }

    fn base_read_info_with_clips(
        qualities: &[u8],
        left_soft_clip_bp: u32,
        right_soft_clip_bp: u32,
    ) -> EndReadInfo {
        EndReadInfo {
            tid: 0,
            interval: Interval::new(10, 24).expect("valid interval"),
            is_reverse: false,
            left_soft_clip_bp,
            right_soft_clip_bp,
            left_motif_has_indels: false,
            right_motif_has_indels: false,
            has_hard_clip: false,
            seq: vec![b'A'; qualities.len()],
            qualities: Some(qualities.to_vec()),
            gc_tag: GCTagValue::default(),
        }
    }

    fn base_end(boundary_pos: u32) -> ResolvedFragmentEnd {
        ResolvedFragmentEnd {
            boundary_pos,
            inside_bases: vec![b'A', b'C', b'G'],
            inside_reference_validation_bp: 3,
        }
    }

    #[test]
    fn apply_base_quality_filters_distinguishes_min_mean_and_max_for_end_filters() {
        // Arrange: one end with qualities [35, 20, 40].
        //
        // Mental derivation:
        // - min = 20, so `min in end >= 30` fails and should drop the end
        // - mean = (35 + 20 + 40) / 3 = 31.666..., so `mean in end >= 30` passes
        // - max = 40, so `max in end < 30` fails and should drop the end
        let read = base_read_info(&[35, 20, 40]);
        let left_end = Some(base_end(10));

        let min_filter = [BaseQualityFilter {
            aggregation: BaseQualityAggregation::Min,
            scope: BaseQualityFilterScope::End,
            op: BaseQualityComparisonOp::Ge,
            threshold: 30.0,
        }];
        let mean_filter = [BaseQualityFilter {
            aggregation: BaseQualityAggregation::Mean,
            scope: BaseQualityFilterScope::End,
            op: BaseQualityComparisonOp::Ge,
            threshold: 30.0,
        }];
        let max_filter = [BaseQualityFilter {
            aggregation: BaseQualityAggregation::Max,
            scope: BaseQualityFilterScope::End,
            op: BaseQualityComparisonOp::Lt,
            threshold: 30.0,
        }];

        // Act / Assert
        assert!(
            apply_base_quality_filters(
                &read,
                &read,
                left_end.clone(),
                None,
                ClipStrategy::Aligned,
                3,
                &min_filter
            )
            .is_none()
        );
        let mean_result = apply_base_quality_filters(
            &read,
            &read,
            left_end.clone(),
            None,
            ClipStrategy::Aligned,
            3,
            &mean_filter,
        );
        assert!(mean_result.is_some());
        let (kept_left_end, kept_right_end) = mean_result.expect("mean filter should keep the end");
        assert!(kept_right_end.is_none());
        let kept_left_end = kept_left_end.expect("left end should stay present");
        let original_left_end = left_end.clone().expect("fixture should contain a left end");
        assert_eq!(kept_left_end.boundary_pos, original_left_end.boundary_pos);
        assert_eq!(kept_left_end.inside_bases, original_left_end.inside_bases);
        assert_eq!(
            kept_left_end.inside_reference_validation_bp,
            original_left_end.inside_reference_validation_bp
        );
        assert!(
            apply_base_quality_filters(
                &read,
                &read,
                left_end,
                None,
                ClipStrategy::Aligned,
                3,
                &max_filter
            )
            .is_none()
        );
    }

    #[test]
    fn apply_base_quality_filters_is_a_no_op_when_no_filters_are_present() {
        // Arrange
        let left_read = base_read_info(&[35, 20, 40]);
        let right_read = base_read_info(&[10, 20, 30]);
        let left_end = Some(base_end(10));
        let right_end = Some(base_end(20));

        // Act
        let result = apply_base_quality_filters(
            &left_read,
            &right_read,
            left_end.clone(),
            right_end.clone(),
            ClipStrategy::Aligned,
            3,
            &[],
        );

        // Assert
        assert!(result.is_some());
        let (kept_left_end, kept_right_end) = result.expect("no filters should preserve both ends");
        let kept_left_end = kept_left_end.expect("left end should stay present");
        let kept_right_end = kept_right_end.expect("right end should stay present");
        let original_left_end = left_end.expect("fixture should contain a left end");
        let original_right_end = right_end.expect("fixture should contain a right end");
        assert_eq!(kept_left_end.boundary_pos, original_left_end.boundary_pos);
        assert_eq!(kept_left_end.inside_bases, original_left_end.inside_bases);
        assert_eq!(
            kept_left_end.inside_reference_validation_bp,
            original_left_end.inside_reference_validation_bp
        );
        assert_eq!(kept_right_end.boundary_pos, original_right_end.boundary_pos);
        assert_eq!(kept_right_end.inside_bases, original_right_end.inside_bases);
        assert_eq!(
            kept_right_end.inside_reference_validation_bp,
            original_right_end.inside_reference_validation_bp
        );
    }

    #[test]
    fn fragment_scope_filters_run_before_end_scope_filters_remove_failing_ends() {
        // Arrange: left qualities [40], right qualities [10].
        //
        // Mental derivation:
        // - fragment mean = (40 + 10) / 2 = 25, so `mean in fragment >= 30` fails
        // - left end still passes `min in end >= 30`
        // - because fragment filters are applied before end filters, the fragment must be dropped
        let left_read = base_read_info(&[40]);
        let right_read = base_read_info(&[10]);
        let filters = [
            BaseQualityFilter {
                aggregation: BaseQualityAggregation::Min,
                scope: BaseQualityFilterScope::End,
                op: BaseQualityComparisonOp::Ge,
                threshold: 30.0,
            },
            BaseQualityFilter {
                aggregation: BaseQualityAggregation::Mean,
                scope: BaseQualityFilterScope::Fragment,
                op: BaseQualityComparisonOp::Ge,
                threshold: 30.0,
            },
        ];

        // Act / Assert
        assert!(
            apply_base_quality_filters(
                &left_read,
                &right_read,
                Some(base_end(10)),
                Some(base_end(20)),
                ClipStrategy::Aligned,
                1,
                &filters
            )
            .is_none()
        );
    }

    #[test]
    fn apply_base_quality_filters_drop_the_fragment_when_fragment_scope_passes_but_both_end_filters_fail()
     {
        // Arrange: both ends have quality 20.
        //
        // Mental derivation:
        // - fragment mean = (20 + 20) / 2 = 20, so `mean in fragment >= 20` passes
        // - `min in end >= 30` fails for both ends
        // - once both ends are dropped, the helper must return `None`
        let left_read = base_read_info(&[20]);
        let right_read = base_read_info(&[20]);
        let filters = [
            BaseQualityFilter {
                aggregation: BaseQualityAggregation::Mean,
                scope: BaseQualityFilterScope::Fragment,
                op: BaseQualityComparisonOp::Ge,
                threshold: 20.0,
            },
            BaseQualityFilter {
                aggregation: BaseQualityAggregation::Min,
                scope: BaseQualityFilterScope::End,
                op: BaseQualityComparisonOp::Ge,
                threshold: 30.0,
            },
        ];

        // Act / Assert
        assert!(
            apply_base_quality_filters(
                &left_read,
                &right_read,
                Some(base_end(10)),
                Some(base_end(20)),
                ClipStrategy::Aligned,
                1,
                &filters
            )
            .is_none()
        );
    }

    #[test]
    fn fragment_scope_filters_distinguish_min_mean_and_max_for_k_inside_gt_one() {
        // Arrange: the fragment has left qualities [40, 35, 30] and right qualities [20, 20, 20].
        //
        // Mental derivation across the raw candidate fragment:
        // - min = 20, so `min in fragment > 20` fails
        // - mean = (40 + 35 + 30 + 20 + 20 + 20) / 6 = 27.5, so `mean in fragment >= 27.5` passes
        // - max = 40, so `max in fragment < 35` fails
        let left_read = base_read_info(&[40, 35, 30]);
        let right_read = base_read_info(&[20, 20, 20]);
        let left_end = Some(base_end(10));
        let right_end = Some(base_end(20));

        let min_filter = [BaseQualityFilter {
            aggregation: BaseQualityAggregation::Min,
            scope: BaseQualityFilterScope::Fragment,
            op: BaseQualityComparisonOp::Gt,
            threshold: 20.0,
        }];
        let mean_filter = [BaseQualityFilter {
            aggregation: BaseQualityAggregation::Mean,
            scope: BaseQualityFilterScope::Fragment,
            op: BaseQualityComparisonOp::Ge,
            threshold: 27.5,
        }];
        let max_filter = [BaseQualityFilter {
            aggregation: BaseQualityAggregation::Max,
            scope: BaseQualityFilterScope::Fragment,
            op: BaseQualityComparisonOp::Lt,
            threshold: 35.0,
        }];

        // Act / Assert
        assert!(
            apply_base_quality_filters(
                &left_read,
                &right_read,
                left_end.clone(),
                right_end.clone(),
                ClipStrategy::Aligned,
                3,
                &min_filter
            )
            .is_none()
        );
        let mean_result = apply_base_quality_filters(
            &left_read,
            &right_read,
            left_end.clone(),
            right_end.clone(),
            ClipStrategy::Aligned,
            3,
            &mean_filter,
        );
        assert!(mean_result.is_some());
        let (kept_left_end, kept_right_end) =
            mean_result.expect("mean fragment filter should keep both ends");
        let kept_left_end = kept_left_end.expect("left end should stay present");
        let kept_right_end = kept_right_end.expect("right end should stay present");
        let original_left_end = left_end.clone().expect("fixture should contain a left end");
        let original_right_end = right_end
            .clone()
            .expect("fixture should contain a right end");
        assert_eq!(kept_left_end.boundary_pos, original_left_end.boundary_pos);
        assert_eq!(kept_left_end.inside_bases, original_left_end.inside_bases);
        assert_eq!(
            kept_left_end.inside_reference_validation_bp,
            original_left_end.inside_reference_validation_bp
        );
        assert_eq!(kept_right_end.boundary_pos, original_right_end.boundary_pos);
        assert_eq!(kept_right_end.inside_bases, original_right_end.inside_bases);
        assert_eq!(
            kept_right_end.inside_reference_validation_bp,
            original_right_end.inside_reference_validation_bp
        );
        assert!(
            apply_base_quality_filters(
                &left_read,
                &right_read,
                left_end,
                right_end,
                ClipStrategy::Aligned,
                3,
                &max_filter
            )
            .is_none()
        );
    }

    #[test]
    fn extract_inside_qualities_respects_clip_strategy_for_k_inside_gt_one() {
        // Arrange: 2S10M2S with per-base qualities 10, 20, ..., 140.
        //
        // Mental derivation for k_inside=3:
        // - aligned left skips the left clips -> [30, 40, 50]
        // - aligned right skips the right clips -> [100, 110, 120]
        // - raw left starts at the raw read edge -> [10, 20, 30]
        // - raw right ends at the raw read edge -> [120, 130, 140]
        let read = base_read_info_with_clips(
            &[10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140],
            2,
            2,
        );

        // Act / Assert
        assert_eq!(
            extract_inside_qualities(&read, FragmentEndSide::Left, ClipStrategy::Aligned, 3),
            Some(vec![30, 40, 50])
        );
        assert_eq!(
            extract_inside_qualities(&read, FragmentEndSide::Right, ClipStrategy::Aligned, 3),
            Some(vec![100, 110, 120])
        );
        assert_eq!(
            extract_inside_qualities(
                &read,
                FragmentEndSide::Left,
                ClipStrategy::IncludeAtAlignedBoundary,
                3
            ),
            Some(vec![10, 20, 30])
        );
        assert_eq!(
            extract_inside_qualities(
                &read,
                FragmentEndSide::Right,
                ClipStrategy::IncludeAtAlignedBoundary,
                3
            ),
            Some(vec![120, 130, 140])
        );
        assert_eq!(
            extract_inside_qualities(
                &read,
                FragmentEndSide::Left,
                ClipStrategy::IncludeAtShiftedBoundary,
                3
            ),
            Some(vec![10, 20, 30])
        );
        assert_eq!(
            extract_inside_qualities(
                &read,
                FragmentEndSide::Right,
                ClipStrategy::IncludeAtShiftedBoundary,
                3
            ),
            Some(vec![120, 130, 140])
        );
    }

    #[test]
    fn inside_slice_bounds_return_the_expected_indices_for_each_strategy_and_side() {
        // Arrange: len=14 with 2 soft-clipped bases on the left and 3 on the right.
        //
        // Mental derivation for k_inside=3:
        // - aligned/skip left starts after the left clips -> [2, 5)
        // - aligned/skip right ends before the right clips -> [8, 11)
        // - include-at-aligned-boundary/include-at-shifted-boundary left starts at the raw read edge -> [0, 3)
        // - include-at-aligned-boundary/include-at-shifted-boundary right ends at the raw read edge -> [11, 14)
        let len = 14;
        let left_soft_clip_bp = 2;
        let right_soft_clip_bp = 3;
        let k_inside = 3;

        // Act / Assert
        assert_eq!(
            inside_slice_bounds(
                len,
                left_soft_clip_bp,
                right_soft_clip_bp,
                FragmentEndSide::Left,
                ClipStrategy::Aligned,
                k_inside
            ),
            Some((2, 5))
        );
        assert_eq!(
            inside_slice_bounds(
                len,
                left_soft_clip_bp,
                right_soft_clip_bp,
                FragmentEndSide::Right,
                ClipStrategy::Aligned,
                k_inside
            ),
            Some((8, 11))
        );
        assert_eq!(
            inside_slice_bounds(
                len,
                left_soft_clip_bp,
                right_soft_clip_bp,
                FragmentEndSide::Left,
                ClipStrategy::Skip,
                k_inside
            ),
            Some((2, 5))
        );
        assert_eq!(
            inside_slice_bounds(
                len,
                left_soft_clip_bp,
                right_soft_clip_bp,
                FragmentEndSide::Right,
                ClipStrategy::Skip,
                k_inside
            ),
            Some((8, 11))
        );
        assert_eq!(
            inside_slice_bounds(
                len,
                left_soft_clip_bp,
                right_soft_clip_bp,
                FragmentEndSide::Left,
                ClipStrategy::IncludeAtAlignedBoundary,
                k_inside
            ),
            Some((0, 3))
        );
        assert_eq!(
            inside_slice_bounds(
                len,
                left_soft_clip_bp,
                right_soft_clip_bp,
                FragmentEndSide::Right,
                ClipStrategy::IncludeAtAlignedBoundary,
                k_inside
            ),
            Some((11, 14))
        );
        assert_eq!(
            inside_slice_bounds(
                len,
                left_soft_clip_bp,
                right_soft_clip_bp,
                FragmentEndSide::Left,
                ClipStrategy::IncludeAtShiftedBoundary,
                k_inside
            ),
            Some((0, 3))
        );
        assert_eq!(
            inside_slice_bounds(
                len,
                left_soft_clip_bp,
                right_soft_clip_bp,
                FragmentEndSide::Right,
                ClipStrategy::IncludeAtShiftedBoundary,
                k_inside
            ),
            Some((11, 14))
        );
    }

    #[test]
    fn inside_slice_bounds_return_zero_width_slices_when_k_inside_is_zero() {
        // Arrange: zero inside bases should produce an empty half-open slice at the side-specific
        // start/end position rather than failing.
        let len = 14;
        let left_soft_clip_bp = 2;
        let right_soft_clip_bp = 3;

        // Act / Assert
        assert_eq!(
            inside_slice_bounds(
                len,
                left_soft_clip_bp,
                right_soft_clip_bp,
                FragmentEndSide::Left,
                ClipStrategy::Aligned,
                0
            ),
            Some((2, 2))
        );
        assert_eq!(
            inside_slice_bounds(
                len,
                left_soft_clip_bp,
                right_soft_clip_bp,
                FragmentEndSide::Right,
                ClipStrategy::Aligned,
                0
            ),
            Some((11, 11))
        );
        assert_eq!(
            inside_slice_bounds(
                len,
                left_soft_clip_bp,
                right_soft_clip_bp,
                FragmentEndSide::Left,
                ClipStrategy::IncludeAtShiftedBoundary,
                0
            ),
            Some((0, 0))
        );
        assert_eq!(
            inside_slice_bounds(
                len,
                left_soft_clip_bp,
                right_soft_clip_bp,
                FragmentEndSide::Right,
                ClipStrategy::IncludeAtShiftedBoundary,
                0
            ),
            Some((14, 14))
        );
    }

    #[test]
    fn inside_slice_bounds_return_none_when_the_requested_inside_span_does_not_fit() {
        // Arrange / Assert
        //
        // Mental derivation:
        // - len=4, left clip=2, aligned left start is 2, so k_inside=3 would need [2, 5) and fail
        // - len=4, right clip=2, aligned right end is 2, so k_inside=3 would need [-1, 2) and fail
        // - len=2 with raw boundaries cannot serve k_inside=3 from either end
        assert_eq!(
            inside_slice_bounds(4, 2, 0, FragmentEndSide::Left, ClipStrategy::Aligned, 3),
            None
        );
        assert_eq!(
            inside_slice_bounds(4, 0, 2, FragmentEndSide::Right, ClipStrategy::Aligned, 3),
            None
        );
        assert_eq!(
            inside_slice_bounds(
                2,
                0,
                0,
                FragmentEndSide::Left,
                ClipStrategy::IncludeAtAlignedBoundary,
                3
            ),
            None
        );
        assert_eq!(
            inside_slice_bounds(
                2,
                0,
                0,
                FragmentEndSide::Right,
                ClipStrategy::IncludeAtShiftedBoundary,
                3
            ),
            None
        );
    }

    #[test]
    fn resolve_fragment_end_does_not_require_read_bases_for_reference_backed_inside_motifs() {
        // Arrange: the read is only 10 bp long, but the requested inside motif is 30 bp. Read-backed
        // inside motifs cannot be sliced from this read. Reference-backed inside motifs should still
        // keep the end, because the inside bases will be loaded from the reference during motif
        // encoding.
        let read = base_read_info(&[40; 10]);
        let k_inside = 30;

        // Act
        let read_backed_left = resolve_fragment_end(
            &read,
            FragmentEndSide::Left,
            ClipStrategy::Aligned,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            k_inside,
            0,
        );
        let reference_backed_left = resolve_fragment_end(
            &read,
            FragmentEndSide::Left,
            ClipStrategy::Aligned,
            KmerSource::Reference,
            IndelMotifFilterPolicy::Auto,
            k_inside,
            0,
        );
        let reference_backed_right = resolve_fragment_end(
            &read,
            FragmentEndSide::Right,
            ClipStrategy::Aligned,
            KmerSource::Reference,
            IndelMotifFilterPolicy::Auto,
            k_inside,
            0,
        );

        // Assert
        assert!(matches!(
            read_backed_left,
            ResolvedEndOutcome::SkipEndDropAssignmentBoundary
        ));

        let ResolvedEndOutcome::KeepEnd {
            assignment_boundary_pos,
            end,
        } = reference_backed_left
        else {
            panic!("reference-backed left end should be kept");
        };
        assert_eq!(assignment_boundary_pos, read.start());
        assert_eq!(end.boundary_pos, read.start());
        assert!(end.inside_bases.is_empty());
        assert_eq!(end.inside_reference_validation_bp, k_inside);

        let ResolvedEndOutcome::KeepEnd {
            assignment_boundary_pos,
            end,
        } = reference_backed_right
        else {
            panic!("reference-backed right end should be kept");
        };
        assert_eq!(assignment_boundary_pos, read.end());
        assert_eq!(end.boundary_pos, read.end());
        assert!(end.inside_bases.is_empty());
        assert_eq!(end.inside_reference_validation_bp, k_inside);
    }

    #[test]
    fn from_record_with_gc_tag_skips_loading_qualities_when_not_requested() {
        // Arrange: qualities of 255 denote missing QVs in BAM, but this should not matter when no
        // base-quality filter is active and the hot path intentionally avoids loading them.
        let mut record = Record::new();
        record.set_tid(0);
        record.set_pos(10);
        record.set(
            b"missing_quals_allowed_without_filter",
            Some(&CigarString(vec![Cigar::Match(4)])),
            b"ACGT",
            &[255, 255, 255, 255],
        );

        // Act
        let read_info =
            EndReadInfo::from_record_with_gc_tag(&record, None, ClipStrategy::Aligned, 1, false)
                .expect("qualities should stay unloaded when BQ filters are absent");

        // Assert
        assert_eq!(read_info.qualities, None);
    }

    #[test]
    fn from_record_with_gc_tag_errors_on_missing_base_qualities_when_requested() {
        // Arrange: BAM encodes missing qualities as 255 placeholders.
        let mut record = Record::new();
        record.set_tid(0);
        record.set_pos(10);
        record.set(
            b"missing_quals",
            Some(&CigarString(vec![Cigar::Match(4)])),
            b"ACGT",
            &[255, 255, 255, 255],
        );

        // Act
        let error =
            EndReadInfo::from_record_with_gc_tag(&record, None, ClipStrategy::Aligned, 1, true)
                .expect_err("missing BAM qualities should fail when BQ filters are active");

        // Assert
        assert!(error.to_string().contains("missing base qualities"));
        assert!(error.to_string().contains("--bq-filter"));
    }
}

mod test_fragment_with_ends_2 {
    use crate::commands::ends::config_structs::{ClipStrategy, KmerSource};
    use crate::shared::fragment::ends_fragment::{
        EndReadInfo, FragmentWithEnds, collect_fragment_with_ends,
        collect_fragment_with_ends_from_single_read,
    };
    use crate::shared::indel_mode::IndelMotifFilterPolicy;
    use rust_htslib::bam::Record;
    use rust_htslib::bam::record::{Aux, Cigar, CigarString};

    const NO_MAX_SOFT_CLIP_LIMIT: u32 = u32::MAX;

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
        k_inside: usize,
        gc_tag: Option<&[u8]>,
    ) -> EndReadInfo {
        EndReadInfo::from_record_with_gc_tag(record, gc_tag, clip_strategy, k_inside, false)
            .expect("end read info")
    }

    fn collect_pair(
        a: &Record,
        b: &Record,
        clip_strategy: ClipStrategy,
        source_inside: KmerSource,
        indel_filter: IndelMotifFilterPolicy,
        k_inside: usize,
        max_soft_clips: u32,
        gc_tag: Option<&[u8]>,
    ) -> Option<FragmentWithEnds> {
        let a_info = end_info(a, clip_strategy, k_inside, gc_tag);
        let b_info = end_info(b, clip_strategy, k_inside, gc_tag);
        collect_fragment_with_ends(
            &a_info,
            &b_info,
            clip_strategy,
            source_inside,
            indel_filter,
            k_inside,
            max_soft_clips,
            &[],
        )
    }

    fn collect_single(
        record: &Record,
        clip_strategy: ClipStrategy,
        source_inside: KmerSource,
        indel_filter: IndelMotifFilterPolicy,
        k_inside: usize,
        max_soft_clips: u32,
        gc_tag: Option<&[u8]>,
    ) -> Option<FragmentWithEnds> {
        let read = end_info(record, clip_strategy, k_inside, gc_tag);
        collect_fragment_with_ends_from_single_read(
            &read,
            clip_strategy,
            source_inside,
            indel_filter,
            k_inside,
            max_soft_clips,
            &[],
        )
    }

    #[test]
    fn endreadinfo_reads_edge_clipping_hard_clips_and_gc_tag() {
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
            NO_MAX_SOFT_CLIP_LIMIT,
            None,
        )
        .expect("fragment");

        assert_eq!(fragment.start(), 100);
        assert_eq!(fragment.end(), 116);
        assert_eq!(fragment.assignment_start(), 100);
        assert_eq!(fragment.assignment_end(), 116);
        assert_eq!(
            fragment.left_end.as_ref().expect("left end").inside_bases,
            b"ACGT".to_vec()
        );
        assert_eq!(
            fragment.right_end.as_ref().expect("right end").inside_bases,
            b"TACG".to_vec()
        );
    }

    #[test]
    fn raw_paired_expands_assignment_interval_and_uses_soft_clipped_bases() {
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
            ClipStrategy::IncludeAtShiftedBoundary,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            NO_MAX_SOFT_CLIP_LIMIT,
            None,
        )
        .expect("fragment");

        assert_eq!(fragment.start(), 100);
        assert_eq!(fragment.end(), 116);
        assert_eq!(fragment.assignment_start(), 98);
        assert_eq!(fragment.assignment_end(), 118);
        assert_eq!(
            fragment.left_end.as_ref().expect("left end").inside_bases,
            b"TTAC".to_vec()
        );
        assert_eq!(
            fragment.right_end.as_ref().expect("right end").inside_bases,
            b"CGTT".to_vec()
        );
    }

    #[test]
    fn aligned_kept_ends_store_aligned_boundary_positions() {
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
            NO_MAX_SOFT_CLIP_LIMIT,
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
    fn raw_kept_ends_store_include_at_shifted_boundary_positions() {
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
            ClipStrategy::IncludeAtShiftedBoundary,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            NO_MAX_SOFT_CLIP_LIMIT,
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
    fn include_at_aligned_boundary_keeps_aligned_positions_but_uses_raw_inside_bases() {
        let forward = make_record(
            b"pair_include_at_aligned_boundary",
            0,
            100,
            false,
            cigar(&[('S', 2), ('M', 6)]),
            b"TTACGTAC",
            None,
        );
        let reverse = make_record(
            b"pair_include_at_aligned_boundary",
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
            ClipStrategy::IncludeAtAlignedBoundary,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            NO_MAX_SOFT_CLIP_LIMIT,
            None,
        )
        .expect("fragment");

        assert_eq!(fragment.assignment_start(), 100);
        assert_eq!(fragment.assignment_end(), 116);
        assert_eq!(
            fragment.left_end.as_ref().expect("left end").boundary_pos,
            100
        );
        assert_eq!(
            fragment.right_end.as_ref().expect("right end").boundary_pos,
            116
        );
        assert_eq!(
            fragment.left_end.as_ref().expect("left end").inside_bases,
            b"TTAC".to_vec()
        );
        assert_eq!(
            fragment.right_end.as_ref().expect("right end").inside_bases,
            b"CGTT".to_vec()
        );
    }

    #[test]
    fn drop_paired_drops_one_clipped_end_but_keeps_other_end() {
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
            ClipStrategy::Skip,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            NO_MAX_SOFT_CLIP_LIMIT,
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
                ClipStrategy::Skip,
                KmerSource::Read,
                IndelMotifFilterPolicy::Auto,
                4,
                NO_MAX_SOFT_CLIP_LIMIT,
                None,
            )
            .is_none()
        );
    }

    #[test]
    fn max_soft_clips_skips_end_and_drops_boundary_adjustment() {
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
            ClipStrategy::IncludeAtShiftedBoundary,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            2,
            None,
        )
        .expect("fragment");

        assert!(fragment.left_end.is_none());
        assert_eq!(fragment.assignment_start(), 100);
        assert_eq!(fragment.assignment_end(), 116);
    }

    #[test]
    fn hard_clipped_pair_is_discarded() {
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
                NO_MAX_SOFT_CLIP_LIMIT,
                None,
            )
            .is_none()
        );
    }

    #[test]
    fn left_indel_is_detected_inside_aligned_footprint() {
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
        let include_at_aligned_boundary_info =
            end_info(&record, ClipStrategy::IncludeAtAlignedBoundary, 4, None);
        let include_at_shifted_boundary_info =
            end_info(&record, ClipStrategy::IncludeAtShiftedBoundary, 4, None);

        assert!(aligned_info.left_motif_has_indels);
        assert!(!include_at_aligned_boundary_info.left_motif_has_indels);
        assert!(!include_at_shifted_boundary_info.left_motif_has_indels);
    }

    #[test]
    fn auto_read_source_keeps_indel_affected_end() {
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
            NO_MAX_SOFT_CLIP_LIMIT,
            None,
        )
        .expect("fragment");

        assert!(fragment.left_end.is_some());
        assert_eq!(fragment.assignment_start(), 100);
    }

    #[test]
    fn auto_reference_source_skips_indel_affected_end_and_keeps_aligned_boundary() {
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
            NO_MAX_SOFT_CLIP_LIMIT,
            None,
        )
        .expect("fragment");

        assert!(fragment.left_end.is_none());
        assert_eq!(fragment.assignment_start(), 100);
    }

    #[test]
    fn skip_affected_end_in_raw_keeps_raw_assignment_boundary() {
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
            ClipStrategy::IncludeAtShiftedBoundary,
            KmerSource::Read,
            IndelMotifFilterPolicy::SkipAffectedEnd,
            4,
            NO_MAX_SOFT_CLIP_LIMIT,
            None,
        )
        .expect("fragment");

        assert!(fragment.left_end.is_none());
        assert_eq!(fragment.assignment_start(), 98);
    }

    #[test]
    fn auto_reference_source_in_raw_keeps_raw_assignment_boundary_when_skipping_end() {
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
            ClipStrategy::IncludeAtShiftedBoundary,
            KmerSource::Reference,
            IndelMotifFilterPolicy::Auto,
            4,
            NO_MAX_SOFT_CLIP_LIMIT,
            None,
        )
        .expect("fragment");

        assert!(fragment.left_end.is_none());
        assert_eq!(fragment.assignment_start(), 98);
    }

    #[test]
    fn skip_affected_fragment_drops_whole_fragment() {
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
                NO_MAX_SOFT_CLIP_LIMIT,
                None,
            )
            .is_none()
        );
    }

    #[test]
    fn drop_clipping_outranks_indel_skip() {
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
            ClipStrategy::Skip,
            KmerSource::Read,
            IndelMotifFilterPolicy::SkipAffectedEnd,
            4,
            NO_MAX_SOFT_CLIP_LIMIT,
            None,
        )
        .expect("fragment");

        assert!(fragment.left_end.is_none());
        assert_eq!(fragment.assignment_start(), 100);
    }

    #[test]
    fn max_soft_clips_outranks_indel_skip_and_drops_raw_assignment_boundary() {
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
            ClipStrategy::IncludeAtShiftedBoundary,
            KmerSource::Reference,
            IndelMotifFilterPolicy::SkipAffectedEnd,
            4,
            1,
            None,
        )
        .expect("fragment");

        assert!(fragment.left_end.is_none());
        assert_eq!(fragment.assignment_start(), 100);
    }

    #[test]
    fn unpaired_aligned_exposes_both_ends() {
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
            NO_MAX_SOFT_CLIP_LIMIT,
            None,
        )
        .expect("fragment");

        assert_eq!(fragment.start(), 100);
        assert_eq!(fragment.end(), 106);
        assert_eq!(fragment.assignment_start(), 100);
        assert_eq!(fragment.assignment_end(), 106);
        assert_eq!(
            fragment.left_end.as_ref().expect("left end").inside_bases,
            b"ACGT".to_vec()
        );
        assert_eq!(
            fragment.right_end.as_ref().expect("right end").inside_bases,
            b"GTAC".to_vec()
        );
    }

    #[test]
    fn unpaired_raw_expands_assignment_interval_on_both_sides() {
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
            ClipStrategy::IncludeAtShiftedBoundary,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            NO_MAX_SOFT_CLIP_LIMIT,
            None,
        )
        .expect("fragment");

        assert_eq!(fragment.start(), 100);
        assert_eq!(fragment.end(), 106);
        assert_eq!(fragment.assignment_start(), 98);
        assert_eq!(fragment.assignment_end(), 108);
        assert_eq!(
            fragment.left_end.as_ref().expect("left end").inside_bases,
            b"TTAC".to_vec()
        );
        assert_eq!(
            fragment.right_end.as_ref().expect("right end").inside_bases,
            b"ACGG".to_vec()
        );
    }

    #[test]
    fn unpaired_include_at_aligned_boundary_keeps_assignment_interval_aligned() {
        let record = make_record(
            b"single_include_at_aligned_boundary",
            0,
            100,
            false,
            cigar(&[('S', 2), ('M', 6), ('S', 2)]),
            b"TTACGTACGG",
            None,
        );

        let fragment = collect_single(
            &record,
            ClipStrategy::IncludeAtAlignedBoundary,
            KmerSource::Read,
            IndelMotifFilterPolicy::Auto,
            4,
            NO_MAX_SOFT_CLIP_LIMIT,
            None,
        )
        .expect("fragment");

        assert_eq!(fragment.start(), 100);
        assert_eq!(fragment.end(), 106);
        assert_eq!(fragment.assignment_start(), 100);
        assert_eq!(fragment.assignment_end(), 106);
        assert_eq!(
            fragment.left_end.as_ref().expect("left end").inside_bases,
            b"TTAC".to_vec()
        );
        assert_eq!(
            fragment.right_end.as_ref().expect("right end").inside_bases,
            b"ACGG".to_vec()
        );
    }

    #[test]
    fn hard_clipped_single_read_is_discarded() {
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
                NO_MAX_SOFT_CLIP_LIMIT,
                None,
            )
            .is_none()
        );
    }

    #[test]
    fn paired_gc_tag_combines_to_fragment_average() {
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
            NO_MAX_SOFT_CLIP_LIMIT,
            Some(b"GC"),
        )
        .expect("fragment");

        assert_eq!(fragment.gc_tag.weight, Some(3.0));
    }

    #[test]
    fn unpaired_gc_tag_is_preserved() {
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
            NO_MAX_SOFT_CLIP_LIMIT,
            Some(b"GC"),
        )
        .expect("fragment");

        assert_eq!(fragment.gc_tag.weight, Some(2.5));
    }

    #[test]
    fn outwards_or_same_strand_pairs_are_rejected() {
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
                NO_MAX_SOFT_CLIP_LIMIT,
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
                NO_MAX_SOFT_CLIP_LIMIT,
                None,
            )
            .is_none()
        );
    }

    #[test]
    fn insufficient_sequence_skips_end_and_returns_none_if_both_ends_fail() {
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
                NO_MAX_SOFT_CLIP_LIMIT,
                None,
            )
            .is_none()
        );
    }
}
