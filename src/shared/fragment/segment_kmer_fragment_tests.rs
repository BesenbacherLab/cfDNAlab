#[cfg(test)]
mod test_kmer_segments {
    use crate::shared::fragment::segment_kmer_fragment::{
        FragmentWithKmerSegments, KmerSegmentedReadInfo, collect_fragment_with_kmer_segments,
    };
    use crate::shared::indel_mode::IndelMode;
    use crate::shared::interval::Interval;
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
        let forward = make_record(0, 100, false, &[Cigar::Match(30)]);
        let reverse = make_record(0, 140, true, &[Cigar::Match(30)]);
        let frag = collect_pair(&forward, &reverse, IndelMode::Ignore, false, 0).expect("fragment");
        assert_eq!(frag.start(), 100);
        assert_eq!(frag.end(), 170);
        assert_eq!(segments(&frag), vec![(100, 130), (140, 170)]);
    }
    #[test]
    fn ignore_mode_includes_gap_when_requested() {
        let forward = make_record(0, 100, false, &[Cigar::Match(30)]);
        let reverse = make_record(0, 140, true, &[Cigar::Match(30)]);
        let frag = collect_pair(&forward, &reverse, IndelMode::Ignore, true, 0).expect("fragment");
        assert_eq!(frag.start(), 100);
        assert_eq!(frag.end(), 170);
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
        assert_eq!(frag.start(), 100);
        assert_eq!(frag.end(), 165);
        assert_eq!(segments(&frag), vec![(100, 110), (110, 120), (150, 165)]);
    }
    #[test]
    fn inter_mate_gap_not_merged_when_border_has_insertion() {
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
        let forward = make_record(0, 100, false, &[Cigar::Match(40)]);
        let reverse = make_record(0, 120, true, &[Cigar::Match(40)]);
        let frag = collect_pair(&forward, &reverse, IndelMode::Ignore, true, 5).expect("fragment");
        assert_eq!(frag.start(), 100);
        assert_eq!(frag.end(), 160);
        assert_eq!(segments(&frag), vec![(105, 155)]);
    }

    #[test]
    fn adjust_mode_handles_mixed_insertions_and_deletions() {
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
        let forward = make_record(0, 100, false, &[Cigar::Match(15), Cigar::Ins(2)]);
        let reverse = make_record(0, 140, true, &[Cigar::Ins(3), Cigar::Match(18)]);

        let frag = collect_pair(&forward, &reverse, IndelMode::Adjust, true, 0).expect("fragment");

        assert_eq!(frag.start(), 100);
        assert_eq!(frag.end(), 158);
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
