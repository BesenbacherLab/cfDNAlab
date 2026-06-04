#![allow(dead_code)]

// KEEP-IN-TESTS: bespoke fixture builders for integration tests that exercise public
// command/API behavior. Generic BAM, two-bit, BED, scaling, GC-package, and output-reader
// helpers live in `cfdnalab::testing`.

use anyhow::Result;
use cfdnalab::run_like_cli::common::{BaseSelectionArgs, FragmentPositionSelectionArgs};
use cfdnalab::run_like_cli::common::{BasesFrom, MismatchBasesFrom, ReferenceFrame};
use cfdnalab::testing::bam::{
    SAM_FLAG_FIRST_MATE, SAM_FLAG_MATE_REVERSE, SAM_FLAG_PROPER_PAIR, SAM_FLAG_SECOND_MATE,
};
use cfdnalab::testing::{
    Cigar, FragmentSpec, ReadSpec, TempBam, TempTwoBit, bam_from_fragments, twobit_from_sequences,
};

pub const LONG_FRAGMENT_LENGTH: i64 = 600;
pub const LONG_FRAGMENT_STARTS: [i64; 10] =
    [0, 400, 800, 1_200, 1_600, 2_000, 2_400, 2_800, 3_200, 3_600];

pub fn late_origin_gc_reference_sequence() -> String {
    let mut sequence = String::with_capacity(1_022);
    sequence.push_str(&"A".repeat(900));
    sequence.push_str(&"C".repeat(61));
    sequence.push_str(&"A".repeat(61));
    sequence
}

pub fn single_position_selection(
    frame: ReferenceFrame,
    positions: &str,
    step: usize,
) -> FragmentPositionSelectionArgs {
    FragmentPositionSelectionArgs {
        frame: vec![frame],
        positions: vec![positions.to_string()],
        step: vec![step],
    }
}

pub fn build_base_selection(
    bases_from: BasesFrom,
    mismatch_bases_from: MismatchBasesFrom,
) -> BaseSelectionArgs {
    BaseSelectionArgs {
        bases_from,
        mismatch_bases_from,
    }
}

fn repeat_pattern(pattern: &[u8], len: usize) -> String {
    let mut buffer = Vec::with_capacity(len);
    for index in 0..len {
        buffer.push(pattern[index % pattern.len()]);
    }
    String::from_utf8(buffer).expect("valid DNA pattern")
}

pub fn complex_reference_twobit() -> Result<TempTwoBit> {
    let chr1 = ("chr1".to_string(), repeat_pattern(b"ACGT", 500));
    let chr2 = ("chr2".to_string(), repeat_pattern(b"TGCA", 400));
    twobit_from_sequences("complex_reference", vec![chr1, chr2])
}

pub fn fragment_kmers_edge_reference() -> Result<TempTwoBit> {
    let chr1 = (
        "chr1".to_string(),
        "ACGTGACCTTAGGCTAACCGTACGTTAGCCGATTACAAGT".to_string(),
    );
    twobit_from_sequences("fragment_kmers_edge", vec![chr1])
}

pub fn fragment_kmers_edge_bam() -> Result<TempBam> {
    let chroms = vec![("chr1".to_string(), 40u32)];

    let fragments = vec![
        FragmentSpec {
            forward: ReadSpec {
                tid: 0,
                pos: 0,
                cigar: vec![Cigar::Match(10)],
                seq: seq(10, b'A'),
                base_quality: 40,
                is_reverse: false,
                mapq: 60,
                flags: SAM_FLAG_FIRST_MATE | SAM_FLAG_MATE_REVERSE | SAM_FLAG_PROPER_PAIR,
                mate_tid: Some(0),
                mate_pos: Some(14),
                insert_size: 24,
            },
            reverse: ReadSpec {
                tid: 0,
                pos: 14,
                cigar: vec![Cigar::Match(10)],
                seq: seq(10, b'T'),
                base_quality: 40,
                is_reverse: true,
                mapq: 60,
                flags: SAM_FLAG_SECOND_MATE | SAM_FLAG_PROPER_PAIR,
                mate_tid: Some(0),
                mate_pos: Some(0),
                insert_size: -24,
            },
        },
        FragmentSpec {
            forward: ReadSpec {
                tid: 0,
                pos: 5,
                cigar: vec![Cigar::Match(4), Cigar::Ins(1), Cigar::Match(4)],
                seq: seq(9, b'C'),
                base_quality: 35,
                is_reverse: false,
                mapq: 55,
                flags: SAM_FLAG_FIRST_MATE | SAM_FLAG_MATE_REVERSE | SAM_FLAG_PROPER_PAIR,
                mate_tid: Some(0),
                mate_pos: Some(13),
                insert_size: 16,
            },
            reverse: ReadSpec {
                tid: 0,
                pos: 13,
                cigar: vec![Cigar::Match(8)],
                seq: seq(8, b'G'),
                base_quality: 35,
                is_reverse: true,
                mapq: 55,
                flags: SAM_FLAG_SECOND_MATE | SAM_FLAG_PROPER_PAIR,
                mate_tid: Some(0),
                mate_pos: Some(5),
                insert_size: -16,
            },
        },
        FragmentSpec {
            forward: ReadSpec {
                tid: 0,
                pos: 16,
                cigar: vec![Cigar::Match(3), Cigar::Del(1), Cigar::Match(5)],
                seq: seq(8, b'A'),
                base_quality: 30,
                is_reverse: false,
                mapq: 50,
                flags: SAM_FLAG_FIRST_MATE | SAM_FLAG_MATE_REVERSE | SAM_FLAG_PROPER_PAIR,
                mate_tid: Some(0),
                mate_pos: Some(20),
                insert_size: 11,
            },
            reverse: ReadSpec {
                tid: 0,
                pos: 20,
                cigar: vec![Cigar::Match(7)],
                seq: seq(7, b'T'),
                base_quality: 30,
                is_reverse: true,
                mapq: 50,
                flags: SAM_FLAG_SECOND_MATE | SAM_FLAG_PROPER_PAIR,
                mate_tid: Some(0),
                mate_pos: Some(16),
                insert_size: -11,
            },
        },
    ];

    bam_from_fragments("fragment_kmers_edge", chroms, fragments, Vec::new())
}

fn seq(len: usize, base: u8) -> Vec<u8> {
    std::iter::repeat(base).take(len).collect()
}

pub fn complex_bam_fixture() -> Result<TempBam> {
    let chroms = vec![("chr1".to_string(), 500u32), ("chr2".to_string(), 400u32)];

    // Diverse fragments covering orientation, indels, skips, mismatched mates, etc.
    let fragments = vec![
        FragmentSpec {
            forward: ReadSpec {
                tid: 0,
                pos: 50,
                cigar: vec![Cigar::Match(40)],
                seq: seq(40, b'A'),
                base_quality: 30,
                is_reverse: false,
                mapq: 60,
                flags: SAM_FLAG_FIRST_MATE | SAM_FLAG_MATE_REVERSE | SAM_FLAG_PROPER_PAIR,
                mate_tid: Some(0),
                mate_pos: Some(120),
                insert_size: 120 - 50 + 40,
            },
            reverse: ReadSpec {
                tid: 0,
                pos: 120,
                cigar: vec![Cigar::Match(40)],
                seq: seq(40, b'T'),
                base_quality: 30,
                is_reverse: true,
                mapq: 60,
                flags: SAM_FLAG_SECOND_MATE | SAM_FLAG_PROPER_PAIR,
                mate_tid: Some(0),
                mate_pos: Some(50),
                insert_size: -(120 - 50 + 40) as i64,
            },
        },
        FragmentSpec {
            forward: ReadSpec {
                tid: 0,
                pos: 200,
                cigar: vec![
                    Cigar::Match(20),
                    Cigar::Ins(3),
                    Cigar::Match(10),
                    Cigar::Del(5),
                    Cigar::Match(12),
                ],
                seq: seq(45, b'C'),
                base_quality: 25,
                is_reverse: false,
                mapq: 50,
                flags: SAM_FLAG_FIRST_MATE | SAM_FLAG_MATE_REVERSE | SAM_FLAG_PROPER_PAIR,
                mate_tid: Some(0),
                mate_pos: Some(260),
                insert_size: 260 - 200 + 50,
            },
            reverse: ReadSpec {
                tid: 0,
                pos: 260,
                cigar: vec![
                    Cigar::SoftClip(5),
                    Cigar::Match(25),
                    Cigar::RefSkip(4),
                    Cigar::Match(16),
                ],
                seq: seq(46, b'G'),
                base_quality: 25,
                is_reverse: true,
                mapq: 40,
                flags: SAM_FLAG_SECOND_MATE | SAM_FLAG_PROPER_PAIR,
                mate_tid: Some(0),
                mate_pos: Some(200),
                insert_size: -(260 - 200 + 50) as i64,
            },
        },
        FragmentSpec {
            forward: ReadSpec {
                tid: 1,
                pos: 30,
                cigar: vec![Cigar::Match(25)],
                seq: seq(25, b'A'),
                base_quality: 30,
                is_reverse: false,
                mapq: 45,
                flags: SAM_FLAG_FIRST_MATE | SAM_FLAG_MATE_REVERSE | SAM_FLAG_PROPER_PAIR,
                mate_tid: Some(1),
                mate_pos: Some(80),
                insert_size: 80 - 30 + 25,
            },
            reverse: ReadSpec {
                tid: 1,
                pos: 80,
                cigar: vec![Cigar::Match(25)],
                seq: seq(25, b'T'),
                base_quality: 30,
                is_reverse: true,
                mapq: 45,
                flags: SAM_FLAG_SECOND_MATE | SAM_FLAG_PROPER_PAIR,
                mate_tid: Some(1),
                mate_pos: Some(30),
                insert_size: -(80 - 30 + 25) as i64,
            },
        },
        FragmentSpec {
            forward: ReadSpec {
                tid: 1,
                pos: 150,
                cigar: vec![Cigar::Match(20)],
                seq: seq(20, b'A'),
                base_quality: 20,
                is_reverse: false,
                mapq: 30,
                flags: SAM_FLAG_FIRST_MATE,
                mate_tid: Some(1),
                mate_pos: Some(180),
                insert_size: 180 - 150 + 20,
            },
            reverse: ReadSpec {
                tid: 1,
                pos: 180,
                cigar: vec![Cigar::Match(20)],
                seq: seq(20, b'C'),
                base_quality: 20,
                is_reverse: false,
                mapq: 30,
                flags: SAM_FLAG_SECOND_MATE,
                mate_tid: Some(1),
                mate_pos: Some(150),
                insert_size: -(180 - 150 + 20) as i64,
            },
        },
    ];

    let singles = vec![
        ReadSpec {
            tid: 0,
            pos: 320,
            cigar: vec![Cigar::Match(30)],
            seq: seq(30, b'A'),
            base_quality: 30,
            is_reverse: false,
            mapq: 10,
            flags: SAM_FLAG_FIRST_MATE | SAM_FLAG_MATE_REVERSE | SAM_FLAG_PROPER_PAIR,
            mate_tid: Some(1),
            mate_pos: Some(100),
            insert_size: 0,
        },
        ReadSpec {
            tid: 1,
            pos: 200,
            cigar: vec![
                Cigar::Match(10),
                Cigar::Pad(5),
                Cigar::Match(10),
                Cigar::HardClip(2),
            ],
            seq: seq(20, b'T'),
            base_quality: 30,
            is_reverse: true,
            mapq: 50,
            flags: SAM_FLAG_SECOND_MATE,
            mate_tid: Some(1),
            mate_pos: Some(210),
            insert_size: 0,
        },
    ];

    bam_from_fragments("complex", chroms, fragments, singles)
}
