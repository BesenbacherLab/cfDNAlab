#![cfg(feature = "cmd_lengths")]

use cfdnalab::shared::clip_mode::ClipMode;
use cfdnalab::shared::fragment::indel_counting_fragment::FragmentWithIndelCounts;
use cfdnalab::shared::fragment_iterators::fragments_with_indel_counts_from_bam;
use cfdnalab::shared::indel_mode::IndelMode;
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

fn make_record(qname: &[u8], pos: i64, len: u32, is_reverse: bool) -> bam::Record {
    let mut rec = bam::Record::new();
    rec.set(qname, None, b"AAAAA", b"IIIII");
    rec.set_cigar(Some(&bam::record::CigarString::from(vec![Cigar::Match(
        len,
    )])));
    rec.set_tid(0);
    rec.set_pos(pos);
    rec.set_flags(if is_reverse { 16 } else { 0 });
    rec
}

fn collect_no_indel_pair(
    inspect_cigar: bool,
) -> (Vec<FragmentWithIndelCounts>, (u64, u64, u64, u64, u64, u64)) {
    let records = vec![
        Ok(make_record(b"r1", 10, 5, false)), // forward [10,15)
        Ok(make_record(b"r1", 20, 5, true)),  // reverse [20,25)
    ];
    let include_read = |_r: &bam::Record| true;
    let mut iter = fragments_with_indel_counts_from_bam(
        records.into_iter(),
        include_read,
        IndelMode::Ignore,
        inspect_cigar,
        |_f: &FragmentWithIndelCounts| true,
        false,
    )
    .with_local_counters();
    let fragments: Vec<_> = iter.by_ref().map(|f| f.unwrap()).collect();
    let counters = iter.counters_snapshot();
    (
        fragments,
        (
            counters.incoming_reads,
            counters.incoming_fragments,
            counters.accepted_forward_reads,
            counters.accepted_reverse_reads,
            counters.produced_fragments,
            counters.yielded_fragments,
        ),
    )
}

fn fragment_signature(
    fragment: &FragmentWithIndelCounts,
) -> (i32, (u32, u32), u32, u32, u32, u32, u32, u32, u32) {
    (
        fragment.tid,
        fragment.interval.as_tuple(),
        fragment.left_soft_clip_bp,
        fragment.right_soft_clip_bp,
        fragment.deletions_nonoverlap,
        fragment.insertions_nonoverlap,
        fragment.deletions_overlap_supported,
        fragment.insertions_overlap_supported,
        fragment.adjusted_len(IndelMode::Ignore, ClipMode::Aligned),
    )
}

#[test]
fn disabled_cigar_inspection_matches_enabled_when_pair_has_no_indels_or_clips() {
    // Arrange + Act: Both reads are pure 5M alignments, so CIGAR inspection should not change the
    // assembled fragment. The fragment span is [10,25), giving length 15.
    let (with_cigar_fragments, with_cigar_counters) = collect_no_indel_pair(true);
    let (without_cigar_fragments, without_cigar_counters) = collect_no_indel_pair(false);

    // Assert: Both paths produce the same hand-derived fragment and the same iterator counters.
    let expected_fragments = vec![(0, (10, 25), 0, 0, 0, 0, 0, 0, 15)];
    let expected_counters = (2, 0, 1, 1, 1, 1);
    assert_eq!(
        with_cigar_fragments
            .iter()
            .map(fragment_signature)
            .collect::<Vec<_>>(),
        expected_fragments
    );
    assert_eq!(
        without_cigar_fragments
            .iter()
            .map(fragment_signature)
            .collect::<Vec<_>>(),
        expected_fragments
    );
    assert_eq!(with_cigar_counters, expected_counters);
    assert_eq!(without_cigar_counters, expected_counters);
}

#[test]
fn yields_unpaired_fragments_and_respects_filter() {
    // Arrange
    let records = vec![
        Ok(make_record(b"r1", 5, 5, false)),  // length 5
        Ok(make_record(b"r2", 30, 5, false)), // length 5
    ];
    let include_read = |_r: &bam::Record| true;
    let mut iter = fragments_with_indel_counts_from_bam(
        records.into_iter(),
        include_read,
        cfdnalab::shared::indel_mode::IndelMode::Ignore,
        true,
        // Filter out anything shorter than 5 (none here) and longer than 5 (none here)
        |f: &FragmentWithIndelCounts| f.adjusted_len(IndelMode::Ignore, ClipMode::Aligned) <= 5,
        true,
    )
    .with_local_counters();

    // Act
    let frags: Vec<_> = iter.by_ref().map(|f| f.unwrap()).collect();

    // Assert
    assert_eq!(frags.len(), 2);
    let lengths: Vec<u32> = frags
        .iter()
        .map(|f| f.adjusted_len(IndelMode::Ignore, ClipMode::Aligned))
        .collect();
    assert_eq!(lengths, vec![5, 5]);
    let snap = iter.counters_snapshot();
    assert_eq!(snap.incoming_reads, 2);
    assert_eq!(snap.yielded_fragments, 2);
}

#[test]
fn pairs_reads_and_yields_single_fragment() {
    // Arrange
    let forward = Ok(make_record(b"r1", 10, 5, false)); // end 15
    let reverse = Ok(make_record(b"r1", 20, 5, true)); // end 25
    let records = vec![forward, reverse];
    let include_read = |_r: &bam::Record| true;
    let mut iter = fragments_with_indel_counts_from_bam(
        records.into_iter(),
        include_read,
        cfdnalab::shared::indel_mode::IndelMode::Ignore,
        true,
        |_f: &FragmentWithIndelCounts| true,
        false,
    )
    .with_local_counters();

    // Act
    let frags: Vec<_> = iter.by_ref().map(|f| f.unwrap()).collect();

    // Assert
    assert_eq!(frags.len(), 1);
    assert_eq!(
        frags[0].adjusted_len(IndelMode::Ignore, ClipMode::Aligned),
        15
    ); // end(reverse) - start(forward)
    let snap = iter.counters_snapshot();
    assert_eq!(snap.incoming_reads, 2);
    assert_eq!(snap.produced_fragments, 1);
    assert_eq!(snap.yielded_fragments, 1);
}
