#![allow(dead_code)]

use anyhow::{Context, Result};
use rust_htslib::bam::{self, header::HeaderRecord, record::Cigar, record::CigarString};
use std::{
    fs::{File, OpenOptions},
    io::{Read, Write},
    path::{Path, PathBuf},
};
use tempfile::TempDir;
use zstd::stream::read::Decoder as ZstdDecoder;

const FLAG_FIRST_MATE: u16 = 0x40;
const FLAG_SECOND_MATE: u16 = 0x80;
const FLAG_PROPER_PAIR: u16 = 0x2;
const FLAG_MATE_REVERSE: u16 = 0x20;

#[derive(Debug)]
pub struct BamFixture {
    _tempdir: TempDir,
    pub bam: PathBuf,
    pub bai: PathBuf,
}

impl BamFixture {
    fn new(tempdir: TempDir, bam: PathBuf, bai: PathBuf) -> Self {
        Self {
            _tempdir: tempdir,
            bam,
            bai,
        }
    }
}

#[derive(Clone)]
pub struct ReadSpec {
    pub tid: usize,
    pub pos: i64,
    pub cigar: Vec<(char, u32)>,
    pub seq: Vec<u8>,
    pub qual: u8,
    pub is_reverse: bool,
    pub mapq: u8,
    pub flags: u16,
    pub mate_tid: Option<usize>,
    pub mate_pos: Option<i64>,
    pub insert_size: i64,
}

impl ReadSpec {
    fn to_record(&self, qname: &[u8]) -> bam::Record {
        let mut rec = bam::Record::new();
        rec.set_tid(self.tid as i32);
        rec.set_pos(self.pos);
        rec.set_mapq(self.mapq);
        if let Some(mtid) = self.mate_tid {
            rec.set_mtid(mtid as i32);
        }
        if let Some(mpos) = self.mate_pos {
            rec.set_mpos(mpos);
        }
        rec.set_insert_size(self.insert_size);
        rec.set(
            qname,
            Some(&cigar(&self.cigar)),
            &self.seq,
            &vec![self.qual; self.seq.len()],
        );
        const FLAG_PAIRED: u16 = 0x1;
        const FLAG_REVERSE: u16 = 0x10;
        let mut flags = self.flags | FLAG_PAIRED;
        if self.is_reverse {
            flags |= FLAG_REVERSE;
        }
        rec.set_flags(flags);
        rec
    }
}

pub struct FragmentSpec {
    pub forward: ReadSpec,
    pub reverse: ReadSpec,
}

fn cigar(ops: &[(char, u32)]) -> CigarString {
    let mut v = Vec::with_capacity(ops.len());
    for (op, len) in ops {
        let c = match *op {
            'M' => Cigar::Match(*len),
            '=' => Cigar::Equal(*len),
            'X' => Cigar::Diff(*len),
            'I' => Cigar::Ins(*len),
            'D' => Cigar::Del(*len),
            'N' => Cigar::RefSkip(*len),
            'S' => Cigar::SoftClip(*len),
            'H' => Cigar::HardClip(*len),
            'P' => Cigar::Pad(*len),
            _ => panic!("Unsupported CIGAR op: {op}"),
        };
        v.push(c);
    }
    CigarString(v)
}

fn write_bam(
    chroms: &[(String, u32)],
    fragments: &[FragmentSpec],
    singles: &[ReadSpec],
    out_bam: &Path,
) -> Result<()> {
    let mut header = bam::Header::new();
    header.push_record(
        HeaderRecord::new(b"HD")
            .push_tag(b"VN", &"1.6")
            .push_tag(b"SO", &"coordinate"),
    );
    for (name, len) in chroms {
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", name)
                .push_tag(b"LN", len),
        );
    }

    let mut writer = bam::Writer::from_path(out_bam, &header, bam::Format::Bam)
        .with_context(|| format!("create bam at {}", out_bam.display()))?;

    let mut records: Vec<bam::Record> = Vec::new();

    for fragment in fragments {
        let qname = format!("frag{}_{}", fragment.forward.tid, fragment.forward.pos);
        records.push(fragment.forward.to_record(qname.as_bytes()));
        records.push(fragment.reverse.to_record(qname.as_bytes()));
    }

    for single in singles {
        let qname = format!("single{}_{}", single.tid, single.pos);
        records.push(single.to_record(qname.as_bytes()));
    }

    records.sort_by_key(|rec| (rec.tid(), rec.pos()));

    for rec in records {
        writer.write(&rec)?;
    }
    Ok(())
}

fn build_index(bam_path: &Path) -> Result<PathBuf> {
    let bai_path = bam_path.with_extension("bam.bai");
    bam::index::build(bam_path, None, bam::index::Type::Bai, 1)
        .with_context(|| format!("index bam {}", bam_path.display()))?;
    let target = bam_path.with_extension("bai");
    if bai_path.exists() {
        std::fs::rename(&bai_path, &target)?;
    }
    Ok(target)
}

fn seq(len: usize, base: u8) -> Vec<u8> {
    std::iter::repeat(base).take(len).collect()
}

pub fn bam_from_specs(
    chroms: Vec<(String, u32)>,
    fragments: Vec<FragmentSpec>,
    singles: Vec<ReadSpec>,
    name: &str,
) -> Result<BamFixture> {
    let tempdir = TempDir::new()?;
    let bam_path = tempdir.path().join(format!("{name}.bam"));

    write_bam(&chroms, &fragments, &singles, &bam_path)?;
    let bai = build_index(&bam_path)?;
    Ok(BamFixture::new(tempdir, bam_path, bai))
}

pub fn complex_bam_fixture() -> Result<BamFixture> {
    let chroms = vec![("chr1".to_string(), 500u32), ("chr2".to_string(), 400u32)];

    // Diverse fragments covering orientation, indels, skips, mismatched mates, etc.
    let fragments = vec![
        FragmentSpec {
            forward: ReadSpec {
                tid: 0,
                pos: 50,
                cigar: vec![('M', 40)],
                seq: seq(40, b'A'),
                qual: 30,
                is_reverse: false,
                mapq: 60,
                flags: FLAG_FIRST_MATE | FLAG_MATE_REVERSE | FLAG_PROPER_PAIR,
                mate_tid: Some(0),
                mate_pos: Some(120),
                insert_size: 120 - 50 + 40,
            },
            reverse: ReadSpec {
                tid: 0,
                pos: 120,
                cigar: vec![('M', 40)],
                seq: seq(40, b'T'),
                qual: 30,
                is_reverse: true,
                mapq: 60,
                flags: FLAG_SECOND_MATE | FLAG_PROPER_PAIR,
                mate_tid: Some(0),
                mate_pos: Some(50),
                insert_size: -(120 - 50 + 40) as i64,
            },
        },
        FragmentSpec {
            forward: ReadSpec {
                tid: 0,
                pos: 200,
                cigar: vec![('M', 20), ('I', 3), ('M', 10), ('D', 5), ('M', 12)],
                seq: seq(45, b'C'),
                qual: 25,
                is_reverse: false,
                mapq: 50,
                flags: FLAG_FIRST_MATE | FLAG_MATE_REVERSE | FLAG_PROPER_PAIR,
                mate_tid: Some(0),
                mate_pos: Some(260),
                insert_size: 260 - 200 + 50,
            },
            reverse: ReadSpec {
                tid: 0,
                pos: 260,
                cigar: vec![('S', 5), ('M', 25), ('N', 4), ('M', 16)],
                seq: seq(46, b'G'),
                qual: 25,
                is_reverse: true,
                mapq: 40,
                flags: FLAG_SECOND_MATE | FLAG_PROPER_PAIR,
                mate_tid: Some(0),
                mate_pos: Some(200),
                insert_size: -(260 - 200 + 50) as i64,
            },
        },
        FragmentSpec {
            forward: ReadSpec {
                tid: 1,
                pos: 30,
                cigar: vec![('M', 25)],
                seq: seq(25, b'A'),
                qual: 30,
                is_reverse: false,
                mapq: 45,
                flags: FLAG_FIRST_MATE | FLAG_MATE_REVERSE | FLAG_PROPER_PAIR,
                mate_tid: Some(1),
                mate_pos: Some(80),
                insert_size: 80 - 30 + 25,
            },
            reverse: ReadSpec {
                tid: 1,
                pos: 80,
                cigar: vec![('M', 25)],
                seq: seq(25, b'T'),
                qual: 30,
                is_reverse: true,
                mapq: 45,
                flags: FLAG_SECOND_MATE | FLAG_PROPER_PAIR,
                mate_tid: Some(1),
                mate_pos: Some(30),
                insert_size: -(80 - 30 + 25) as i64,
            },
        },
        FragmentSpec {
            forward: ReadSpec {
                tid: 1,
                pos: 150,
                cigar: vec![('M', 20)],
                seq: seq(20, b'A'),
                qual: 20,
                is_reverse: false,
                mapq: 30,
                flags: FLAG_FIRST_MATE,
                mate_tid: Some(1),
                mate_pos: Some(180),
                insert_size: 180 - 150 + 20,
            },
            reverse: ReadSpec {
                tid: 1,
                pos: 180,
                cigar: vec![('M', 20)],
                seq: seq(20, b'C'),
                qual: 20,
                is_reverse: false,
                mapq: 30,
                flags: FLAG_SECOND_MATE,
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
            cigar: vec![('M', 30)],
            seq: seq(30, b'A'),
            qual: 30,
            is_reverse: false,
            mapq: 10,
            flags: FLAG_FIRST_MATE | FLAG_MATE_REVERSE | FLAG_PROPER_PAIR,
            mate_tid: Some(1),
            mate_pos: Some(100),
            insert_size: 0,
        },
        ReadSpec {
            tid: 1,
            pos: 200,
            cigar: vec![('M', 10), ('P', 5), ('M', 10), ('H', 2)],
            seq: seq(20, b'T'),
            qual: 30,
            is_reverse: true,
            mapq: 50,
            flags: FLAG_SECOND_MATE,
            mate_tid: Some(1),
            mate_pos: Some(210),
            insert_size: 0,
        },
    ];

    bam_from_specs(chroms, fragments, singles, "complex")
}

pub fn simple_inward_bam() -> Result<BamFixture> {
    let chroms = vec![("chr1".to_string(), 200u32)];
    let fragments = vec![FragmentSpec {
        forward: ReadSpec {
            tid: 0,
            pos: 20,
            cigar: vec![('M', 20)],
            seq: seq(20, b'A'),
            qual: 35,
            is_reverse: false,
            mapq: 60,
            flags: FLAG_FIRST_MATE | FLAG_MATE_REVERSE | FLAG_PROPER_PAIR,
            mate_tid: Some(0),
            mate_pos: Some(60),
            insert_size: 60 - 20 + 20,
        },
        reverse: ReadSpec {
            tid: 0,
            pos: 60,
            cigar: vec![('M', 20)],
            seq: seq(20, b'T'),
            qual: 35,
            is_reverse: true,
            mapq: 60,
            flags: FLAG_SECOND_MATE | FLAG_PROPER_PAIR,
            mate_tid: Some(0),
            mate_pos: Some(20),
            insert_size: -(60 - 20 + 20) as i64,
        },
    }];
    bam_from_specs(chroms, fragments, Vec::new(), "simple_inward")
}

pub fn write_bed<P: AsRef<Path>>(path: P, rows: &[(&str, u64, u64, &str)]) -> Result<()> {
    let mut f = File::create(path)?;
    for (chr, start, end, name) in rows {
        writeln!(f, "{}\t{}\t{}\t{}", chr, start, end, name)?;
    }
    Ok(())
}

pub fn write_scaling_factors<P: AsRef<Path>>(
    path: P,
    rows: &[(&str, u64, u64, f32)],
) -> Result<()> {
    let mut f = File::create(path)?;
    for (chr, start, end, factor) in rows {
        writeln!(f, "{}\t{}\t{}\t{}", chr, start, end, factor)?;
    }
    Ok(())
}

pub fn read_zst_to_string(path: &Path) -> Result<String> {
    let reader = File::open(path)?;
    let mut decoder = ZstdDecoder::new(reader)?;
    let mut buf = String::new();
    decoder.read_to_string(&mut buf)?;
    Ok(buf)
}

pub fn read_binary_zst(path: &Path) -> Result<Vec<u8>> {
    let reader = File::open(path)?;
    let mut decoder = ZstdDecoder::new(reader)?;
    let mut buf = Vec::new();
    decoder.read_to_end(&mut buf)?;
    Ok(buf)
}

pub fn touch_file<P: AsRef<Path>>(path: P) -> Result<()> {
    OpenOptions::new().create(true).write(true).open(path)?;
    Ok(())
}
