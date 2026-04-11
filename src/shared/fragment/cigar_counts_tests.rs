use super::{CigarIndelInfo, InsertionAnchor, inspect_cigar_indels};
use crate::shared::interval::Interval;
use rust_htslib::bam::{
    Record,
    record::{Cigar, CigarString},
};

fn make_record(pos: i64, cigar_ops: Vec<Cigar>) -> Record {
    let mut record = Record::new();
    record.set_tid(0);
    record.set_pos(pos);
    record.set_flags(0);
    record.set_mapq(60);
    let seq_len: usize = cigar_ops
        .iter()
        .map(|op| match *op {
            Cigar::Match(bp)
            | Cigar::Equal(bp)
            | Cigar::Diff(bp)
            | Cigar::Ins(bp)
            | Cigar::SoftClip(bp) => bp as usize,
            Cigar::Del(_) | Cigar::RefSkip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => 0,
        })
        .sum();
    record.set(
        b"indels",
        Some(&CigarString(cigar_ops)),
        &vec![b'A'; seq_len],
        &vec![30_u8; seq_len],
    );
    record
}

#[test]
fn inspect_cigar_indels_merges_touching_deletion_like_ops_and_keeps_insert_anchor() {
    // Human verification status: unverified
    // Start at 100 with cigar 5M 2D 3N 4M 3I 2M.
    //
    // Reference walk:
    // - 5M  => ref 105
    // - 2D  => deletion [105,107), ref 107
    // - 3N  => deletion [107,110), ref 110
    // - 4M  => ref 114
    // - 3I  => insertion anchored at 114
    //
    // The touching deletion-like intervals [105,107) and [107,110) must merge to [105,110).
    let record = make_record(
        100,
        vec![
            Cigar::Match(5),
            Cigar::Del(2),
            Cigar::RefSkip(3),
            Cigar::Match(4),
            Cigar::Ins(3),
            Cigar::Match(2),
        ],
    );

    let indels = inspect_cigar_indels(&record);

    assert_eq!(
        indels,
        CigarIndelInfo {
            deletions: vec![Interval::new(105, 110).expect("expected merged deletion interval")],
            insertions: vec![InsertionAnchor {
                reference_position: 114,
                inserted_length: 3,
            }],
        }
    );
}

#[test]
fn inspect_cigar_indels_ignores_terminal_clips_and_pad_for_reference_walk() {
    // Human verification status: unverified
    // Start at 50 with cigar 2S 4M 2I 3M 1P 2H.
    //
    // Reference walk:
    // - 2S  => ignored, ref stays 50
    // - 4M  => ref 54
    // - 2I  => insertion anchored at 54
    // - 3M  => ref 57
    // - 1P  => ignored, ref stays 57
    // - 2H  => ignored, ref stays 57
    let record = make_record(
        50,
        vec![
            Cigar::SoftClip(2),
            Cigar::Match(4),
            Cigar::Ins(2),
            Cigar::Match(3),
            Cigar::Pad(1),
            Cigar::HardClip(2),
        ],
    );

    let indels = inspect_cigar_indels(&record);

    assert_eq!(
        indels,
        CigarIndelInfo {
            deletions: Vec::new(),
            insertions: vec![InsertionAnchor {
                reference_position: 54,
                inserted_length: 2,
            }],
        }
    );
}
