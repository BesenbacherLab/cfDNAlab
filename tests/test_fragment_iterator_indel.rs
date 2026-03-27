use cfdnalab::shared::fragment::indel_counting_fragment::FragmentWithIndelCounts;
use cfdnalab::shared::fragment_iterator::fragments_with_indel_counts_from_bam;
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

#[test]
fn yields_unpaired_fragments_and_respects_filter() {
    // Human verification status: unverified
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
        // Filter out anything shorter than 5 (none here) and longer than 5 (none here)
        |f: &FragmentWithIndelCounts| f.len_indel_adjusted() <= 5,
        true,
    )
    .with_local_counters();

    // Act
    let frags: Vec<_> = iter.by_ref().map(|f| f.unwrap()).collect();

    // Assert
    assert_eq!(frags.len(), 2);
    let lengths: Vec<u32> = frags.iter().map(|f| f.len_indel_adjusted()).collect();
    assert_eq!(lengths, vec![5, 5]);
    let snap = iter.counters_snapshot();
    assert_eq!(snap.incoming_reads, 2);
    assert_eq!(snap.yielded_fragments, 2);
}

#[test]
fn pairs_reads_and_yields_single_fragment() {
    // Human verification status: unverified
    // Arrange
    let forward = Ok(make_record(b"r1", 10, 5, false)); // end 15
    let reverse = Ok(make_record(b"r1", 20, 5, true)); // end 25
    let records = vec![forward, reverse];
    let include_read = |_r: &bam::Record| true;
    let mut iter = fragments_with_indel_counts_from_bam(
        records.into_iter(),
        include_read,
        cfdnalab::shared::indel_mode::IndelMode::Ignore,
        |_f: &FragmentWithIndelCounts| true,
        false,
    )
    .with_local_counters();

    // Act
    let frags: Vec<_> = iter.by_ref().map(|f| f.unwrap()).collect();

    // Assert
    assert_eq!(frags.len(), 1);
    assert_eq!(frags[0].len_indel_adjusted(), 15); // end(reverse) - start(forward)
    let snap = iter.counters_snapshot();
    assert_eq!(snap.incoming_reads, 2);
    assert_eq!(snap.produced_fragments, 1);
    assert_eq!(snap.yielded_fragments, 1);
}
