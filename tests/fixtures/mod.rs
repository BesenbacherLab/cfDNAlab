#![allow(dead_code)]

use anyhow::{anyhow, ensure, Context, Result};
use cfdnalab::commands::cli_common::{BaseSelectionArgs, FragmentPositionSelectionArgs};
#[cfg(all(feature = "cmd_gc_bias", feature = "cmd_ref_gc_bias"))]
use cfdnalab::commands::cli_common::{
    ChromosomeArgs, GCWindowsArgs, IOCArgs, LoggingArgs, Ref2BitRequiredArgs,
};
#[cfg(all(feature = "cmd_gc_bias", feature = "cmd_ref_gc_bias"))]
use cfdnalab::commands::gc_bias::{config::GCConfig, gc_bias::run as run_gc_bias};
#[cfg(all(feature = "cmd_gc_bias", feature = "cmd_ref_gc_bias"))]
use cfdnalab::commands::ref_gc_bias::{
    config::RefGCBiasConfig, ref_gc_bias::run as run_ref_gc_bias,
};
use cfdnalab::shared::positioning::{BasesFrom, MismatchBasesFrom, ReferenceFrame};
#[cfg(all(feature = "cmd_gc_bias", feature = "cmd_ref_gc_bias"))]
use cfdnalab::shared::reference::twobit_contig_lengths;
use ndarray::{Array2, Array3};
use rust_htslib::bam::{self, header::HeaderRecord, record::Cigar, record::CigarString};
use std::{
    fs::{File, OpenOptions},
    io::{Read, Write},
    path::{Path, PathBuf},
    sync::Arc,
};
use tempfile::TempDir;
use twobit::convert::{fasta::FastaReader, to_2bit};
use zarrs::{array::Array, filesystem::FilesystemStore};
use zstd::stream::read::Decoder as ZstdDecoder;

const FLAG_FIRST_MATE: u16 = 0x40;
const FLAG_SECOND_MATE: u16 = 0x80;
const FLAG_PROPER_PAIR: u16 = 0x2;
const FLAG_MATE_REVERSE: u16 = 0x20;

pub const LONG_FRAGMENT_LENGTH: i64 = 600;
pub const LONG_FRAGMENT_READ_LENGTH: i64 = 100;
pub const LONG_FRAGMENT_STARTS: [i64; 10] =
    [0, 400, 800, 1_200, 1_600, 2_000, 2_400, 2_800, 3_200, 3_600];

/// Read the `counts` array from a midpoint Zarr output.
pub fn read_midpoint_zarr_counts<P: AsRef<Path>>(store_path: P) -> Result<Array3<f32>> {
    let array = open_zarr_array(store_path.as_ref(), "/counts")?;
    let shape = array.shape();
    ensure!(
        shape.len() == 3,
        "expected midpoint Zarr counts to be rank 3 but found rank {}",
        shape.len()
    );
    let values: Vec<f32> = array
        .retrieve_array_subset(&array.subset_all())
        .context("reading midpoint Zarr counts")?;
    let shape = (
        usize::try_from(shape[0]).context("group dimension exceeds usize")?,
        usize::try_from(shape[1]).context("length_bin dimension exceeds usize")?,
        usize::try_from(shape[2]).context("position dimension exceeds usize")?,
    );
    Array3::from_shape_vec(shape, values).context("building midpoint count array from Zarr values")
}

/// Read a one-dimensional signed-integer array from a midpoint Zarr output.
pub fn read_midpoint_zarr_i32_1d<P: AsRef<Path>>(
    store_path: P,
    array_path: &str,
) -> Result<Vec<i32>> {
    let array = open_zarr_array(store_path.as_ref(), array_path)?;
    ensure!(
        array.shape().len() == 1,
        "expected midpoint Zarr array {array_path} to be rank 1"
    );
    array
        .retrieve_array_subset(&array.subset_all())
        .with_context(|| format!("reading midpoint Zarr array {array_path}"))
}

/// Read a one-dimensional unsigned-integer array from a midpoint Zarr output.
pub fn read_midpoint_zarr_u32_1d<P: AsRef<Path>>(
    store_path: P,
    array_path: &str,
) -> Result<Vec<u32>> {
    let array = open_zarr_array(store_path.as_ref(), array_path)?;
    ensure!(
        array.shape().len() == 1,
        "expected midpoint Zarr array {array_path} to be rank 1"
    );
    array
        .retrieve_array_subset(&array.subset_all())
        .with_context(|| format!("reading midpoint Zarr array {array_path}"))
}

fn open_zarr_array(store_path: &Path, array_path: &str) -> Result<Array<FilesystemStore>> {
    let store = Arc::new(
        FilesystemStore::new(store_path)
            .with_context(|| format!("opening Zarr store {}", store_path.display()))?,
    );
    Array::open(store, array_path).with_context(|| {
        format!(
            "opening Zarr array {array_path} in {}",
            store_path.display()
        )
    })
}

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

/// Build a paired-end fragment from a start position plus fragment/read lengths.
pub fn paired_fragment(start: i64, fragment_len: i64, read_len: i64) -> FragmentSpec {
    let reverse_start = start + fragment_len - read_len;
    let insert_size = fragment_len;
    let forward = ReadSpec {
        tid: 0,
        pos: start,
        cigar: vec![('M', read_len as u32)],
        seq: vec![b'A'; read_len as usize],
        qual: 40,
        is_reverse: false,
        mapq: 60,
        flags: FLAG_FIRST_MATE | FLAG_MATE_REVERSE | FLAG_PROPER_PAIR,
        mate_tid: Some(0),
        mate_pos: Some(reverse_start),
        insert_size,
    };
    let reverse = ReadSpec {
        tid: 0,
        pos: reverse_start,
        cigar: vec![('M', read_len as u32)],
        seq: vec![b'T'; read_len as usize],
        qual: 40,
        is_reverse: true,
        mapq: 60,
        flags: FLAG_SECOND_MATE | FLAG_PROPER_PAIR,
        mate_tid: Some(0),
        mate_pos: Some(start),
        insert_size: -insert_size,
    };
    FragmentSpec { forward, reverse }
}

/// Create a BAM fixture from evenly described fragment starts.
pub fn bam_from_fragment_starts(
    name: &str,
    starts: &[i64],
    fragment_len: i64,
    read_len: i64,
) -> Result<BamFixture> {
    let max_end = starts
        .iter()
        .map(|start| start + fragment_len)
        .max()
        .unwrap_or(0);
    let chrom_len = (max_end + 500).max(1_000) as u32;
    let fragments: Vec<FragmentSpec> = starts
        .iter()
        .map(|&start| paired_fragment(start, fragment_len, read_len))
        .collect();
    bam_from_specs(
        vec![("chr1".to_string(), chrom_len)],
        fragments,
        Vec::new(),
        name,
    )
}

/// Convenience helper for the shared 10-fragment, 600bp scenario used by WPS tests.
pub fn long_fragment_bam(name: &str) -> Result<BamFixture> {
    bam_from_fragment_starts(
        name,
        &LONG_FRAGMENT_STARTS,
        LONG_FRAGMENT_LENGTH,
        LONG_FRAGMENT_READ_LENGTH,
    )
}

#[derive(Debug)]
pub struct TwoBitFixture {
    _tempdir: TempDir,
    pub path: PathBuf,
    sequences: Vec<(String, String)>,
}

impl TwoBitFixture {
    fn new(tempdir: TempDir, path: PathBuf, sequences: Vec<(String, String)>) -> Self {
        Self {
            _tempdir: tempdir,
            path,
            sequences,
        }
    }

    pub fn sequence(&self, chr: &str) -> Option<&str> {
        self.sequences
            .iter()
            .find(|(name, _)| name == chr)
            .map(|(_, seq)| seq.as_str())
    }

    pub fn sequences(&self) -> &[(String, String)] {
        &self.sequences
    }
}

pub fn twobit_from_sequences(
    name: &str,
    sequences: Vec<(String, String)>,
) -> Result<TwoBitFixture> {
    let normalized: Vec<(String, String)> = sequences
        .into_iter()
        .map(|(chr, seq)| (chr, seq.to_ascii_uppercase()))
        .collect();
    let tempdir = TempDir::new()?;
    let fasta_path = tempdir.path().join(format!("{name}.fasta"));
    write_fasta(&fasta_path, &normalized)?;
    let path = tempdir.path().join(format!("{name}.2bit"));
    {
        let reader = FastaReader::open(&fasta_path).map_err(|e| anyhow!(e.to_string()))?;
        let mut file = File::create(&path)?;
        to_2bit(&mut file, &reader).map_err(|e| anyhow!(e.to_string()))?;
    }
    Ok(TwoBitFixture::new(tempdir, path, normalized))
}

fn write_fasta<P: AsRef<Path>>(path: P, sequences: &[(String, String)]) -> Result<()> {
    let mut file = File::create(path)?;
    for (name, seq) in sequences {
        writeln!(file, ">{name}")?;
        for chunk in seq.as_bytes().chunks(60) {
            file.write_all(chunk)?;
            file.write_all(b"\n")?;
        }
    }
    Ok(())
}

pub fn simple_reference_twobit() -> Result<TwoBitFixture> {
    let chr1 = ("chr1".to_string(), "ACGT".repeat(64));
    twobit_from_sequences("simple_reference", vec![chr1])
}

#[cfg(all(feature = "cmd_gc_bias", feature = "cmd_ref_gc_bias"))]
fn base_chromosomes(chromosome_names: &[&str]) -> ChromosomeArgs {
    ChromosomeArgs {
        chromosomes: Some(
            chromosome_names
                .iter()
                .map(|name| name.to_string())
                .collect(),
        ),
        chromosomes_file: None,
    }
}

#[cfg(all(feature = "cmd_gc_bias", feature = "cmd_ref_gc_bias"))]
fn configure_gc_bias_common(gc_cfg: &mut GCConfig) {
    gc_cfg.set_min_mapq(0);
    gc_cfg.set_tile_size(1_000_000);
    gc_cfg.set_min_window_acgt_pct(0);
    gc_cfg.set_num_extreme_gc_bins(0);
    gc_cfg.set_num_short_length_bins(0);
    gc_cfg.outlier_method = cfdnalab::commands::gc_bias::config::OutlierMethodArg::None;
    gc_cfg.set_windows(GCWindowsArgs {
        by_size: None,
        by_bed: None,
        global: true,
    });
}

#[cfg(all(feature = "cmd_gc_bias", feature = "cmd_ref_gc_bias"))]
/// Build a "neutral" real GC-correction package from the full chromosome.
///
/// This helper runs the actual `ref-gc-bias` -> `gc-bias` producer chain, but keeps the reference
/// side intentionally broad and well-mixed:
/// - reference windows are global rather than BED-restricted
/// - `n_positions` is fixed to a moderate deterministic sample
/// - the resulting package is used in tests that want a valid end-to-end artifact without
///   depending on strong hand-derived GC skew
///
/// In other words, "neutral" here means "safe default real artifact for downstream consumers",
/// not "mathematically forced to exactly weight 1.0 everywhere".
pub fn build_real_neutral_gc_package_for_range(
    bam_path: &Path,
    reference_path: &Path,
    out_dir: &Path,
    min_fragment_length: u32,
    max_fragment_length: u32,
) -> Result<PathBuf> {
    let chromosomes = vec!["chr1".to_string()];
    let chrom_lengths = twobit_contig_lengths(reference_path, &chromosomes)?;
    let total_possible_starts: usize = chrom_lengths
        .values()
        .map(|&chrom_len| {
            chrom_len
                .checked_sub(max_fragment_length as usize)
                .map(|remaining| remaining + 1)
                .unwrap_or(0)
        })
        .sum();
    let n_positions = total_possible_starts.min(100);
    anyhow::ensure!(
        n_positions > 0,
        "neutral GC fixture has no valid reference start positions for fragment length range {}-{}",
        min_fragment_length,
        max_fragment_length
    );

    let ref_gc_dir = TempDir::new()?;
    let ref_cfg = RefGCBiasConfig {
        ref_genome: Ref2BitRequiredArgs {
            ref_2bit: reference_path.to_path_buf(),
        },
        output_dir: ref_gc_dir.path().to_path_buf(),
        output_prefix: String::new(),
        n_threads: 1,
        // Use a deterministic but valid number of sampled starts, capped at 100 while still
        // respecting the number of positions where the maximum fragment length fits in the
        // reference. Wider fixture length ranges can make the old fixed sample count invalid on
        // the tiny 256 bp test reference.
        n_positions,
        seed: Some(7),
        windows: Default::default(),
        chromosomes: base_chromosomes(&["chr1"]),
        blacklist: None,
        fragment_lengths: cfdnalab::commands::cli_common::FragmentLengthArgs {
            min_fragment_length,
            max_fragment_length,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 1_000_000,
        logging: LoggingArgs::default(),
    };
    run_ref_gc_bias(&ref_cfg)?;

    let gc_out_dir = out_dir.join(format!(
        "real_gc_bias_neutral_len_{}-{}",
        min_fragment_length, max_fragment_length
    ));
    std::fs::create_dir_all(&gc_out_dir)?;
    let mut gc_cfg = GCConfig::new(
        IOCArgs {
            bam: bam_path.to_path_buf(),
            output_dir: gc_out_dir.clone(),
            n_threads: 1,
        },
        reference_path.to_path_buf(),
        ref_gc_dir.path().join("ref_gc_package.npz"),
        base_chromosomes(&["chr1"]),
    );
    configure_gc_bias_common(&mut gc_cfg);
    // Keep sparse length ranges valid in tiny real-artifact fixtures. Tests that use this helper
    // often exercise only a few observed lengths but still need the package to honestly cover a
    // broader configured length range.
    gc_cfg.set_min_length_bin_mass(0.0);
    gc_cfg.set_min_length_bin_width(1);
    run_gc_bias(&gc_cfg)?;

    Ok(gc_out_dir.join("gc_bias_correction.npz"))
}

#[cfg(all(feature = "cmd_gc_bias", feature = "cmd_ref_gc_bias"))]
pub fn build_real_neutral_gc_package(
    bam_path: &Path,
    reference_path: &Path,
    out_dir: &Path,
    fragment_length: u32,
) -> Result<PathBuf> {
    build_real_neutral_gc_package_for_range(
        bam_path,
        reference_path,
        out_dir,
        fragment_length,
        fragment_length,
    )
}

#[cfg(feature = "cmd_gc_bias")]
pub fn write_constant_gc_package(path: &Path, fragment_length: u32, weight: f64) -> Result<()> {
    let package = cfdnalab::commands::gc_bias::package::GCCorrectionPackage {
        version: cfdnalab::shared::constants::GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![fragment_length, fragment_length + 1],
        gc_edges: vec![0, 101],
        length_bin_frequencies: ndarray::array![1.0_f64],
        reference_contig_footprint: Vec::new(),
        correction_matrix: ndarray::array![[weight]],
    };
    package.write_npz(path)?;
    Ok(())
}

#[cfg(feature = "cmd_gc_bias")]
pub fn write_two_bin_gc_package(
    path: &Path,
    fragment_length: u32,
    low_gc_weight: f64,
    high_gc_weight: f64,
    reference_contig_footprint: Vec<cfdnalab::shared::reference::ContigFootprintEntry>,
) -> Result<()> {
    let package = cfdnalab::commands::gc_bias::package::GCCorrectionPackage {
        version: cfdnalab::shared::constants::GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![fragment_length, fragment_length + 1],
        gc_edges: vec![0, 51, 101],
        length_bin_frequencies: ndarray::array![1.0_f64],
        reference_contig_footprint,
        correction_matrix: ndarray::array![[low_gc_weight, high_gc_weight]],
    };
    package.write_npz(path)?;
    Ok(())
}

pub fn late_origin_gc_reference_sequence() -> String {
    let mut sequence = String::with_capacity(1_022);
    sequence.push_str(&"A".repeat(900));
    sequence.push_str(&"C".repeat(61));
    sequence.push_str(&"A".repeat(61));
    sequence
}

#[cfg(all(feature = "cmd_gc_bias", feature = "cmd_ref_gc_bias"))]
/// Build a deliberately non-neutral real GC-correction package from caller-supplied BED windows.
///
/// Unlike `build_real_neutral_gc_package`, this helper is for tests that need the producer side to
/// create a specific, auditable GC bias profile:
/// - the reference support is restricted to the provided BED intervals
/// - callers choose `n_positions` explicitly so the sampled-start arithmetic stays derivable at the
///   test site
/// - `gc-bias` is configured with permissive bin-mass thresholds so sparse but intentionally skewed
///   fixture counts survive into the saved package
///
/// Use this when the test needs a real package whose non-unit weights can be reasoned about from
/// first principles, rather than just "some valid package".
pub fn build_real_non_neutral_gc_package(
    bam_path: &Path,
    reference_path: &Path,
    out_dir: &Path,
    fragment_length: u32,
    reference_windows_bed: &str,
    n_positions: usize,
) -> Result<PathBuf> {
    let ref_gc_dir = TempDir::new()?;
    let bed_path = ref_gc_dir.path().join("pure_windows.bed");
    std::fs::write(&bed_path, reference_windows_bed)?;
    let ref_cfg = RefGCBiasConfig {
        ref_genome: Ref2BitRequiredArgs {
            ref_2bit: reference_path.to_path_buf(),
        },
        output_dir: ref_gc_dir.path().to_path_buf(),
        output_prefix: String::new(),
        n_threads: 1,
        // Callers pass the exact number of sampled starts because the non-neutral tests rely on
        // fully hand-derived reference-side counts. `ref-gc-bias` only counts sampled starts that
        // both lie inside the BED interval and leave enough room for the full fragment before that
        // interval's right edge, so keeping the arithmetic at the call site makes the expected
        // producer/consumer weights easy to audit next to each test.
        n_positions,
        seed: Some(23),
        windows: cfdnalab::commands::ref_gc_bias::config::RefGCWindowsArgs {
            by_bed: Some(bed_path),
        },
        chromosomes: base_chromosomes(&["chr1"]),
        blacklist: None,
        fragment_lengths: cfdnalab::commands::cli_common::FragmentLengthArgs {
            min_fragment_length: fragment_length,
            max_fragment_length: fragment_length,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 1_000_000,
        logging: LoggingArgs::default(),
    };
    run_ref_gc_bias(&ref_cfg)?;

    let gc_out_dir = out_dir.join(format!("real_gc_bias_non_neutral_len_{fragment_length}"));
    std::fs::create_dir_all(&gc_out_dir)?;
    let mut gc_cfg = GCConfig::new(
        IOCArgs {
            bam: bam_path.to_path_buf(),
            output_dir: gc_out_dir.clone(),
            n_threads: 1,
        },
        reference_path.to_path_buf(),
        ref_gc_dir.path().join("ref_gc_package.npz"),
        base_chromosomes(&["chr1"]),
    );
    configure_gc_bias_common(&mut gc_cfg);
    gc_cfg.set_min_gc_bin_mass(1.0);
    gc_cfg.set_min_length_bin_mass(0.0);
    gc_cfg.set_min_length_bin_width(1);
    run_gc_bias(&gc_cfg)?;

    Ok(gc_out_dir.join("gc_bias_correction.npz"))
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
    let mut buf = Vec::with_capacity(len);
    for i in 0..len {
        buf.push(pattern[i % pattern.len()]);
    }
    String::from_utf8(buf).expect("valid DNA pattern")
}

pub fn complex_reference_twobit() -> Result<TwoBitFixture> {
    let chr1 = ("chr1".to_string(), repeat_pattern(b"ACGT", 500));
    let chr2 = ("chr2".to_string(), repeat_pattern(b"TGCA", 400));
    twobit_from_sequences("complex_reference", vec![chr1, chr2])
}

pub fn fragment_kmers_edge_reference() -> Result<TwoBitFixture> {
    let chr1 = (
        "chr1".to_string(),
        "ACGTGACCTTAGGCTAACCGTACGTTAGCCGATTACAAGT".to_string(),
    );
    twobit_from_sequences("fragment_kmers_edge", vec![chr1])
}

pub fn fragment_kmers_edge_bam() -> Result<BamFixture> {
    let chroms = vec![("chr1".to_string(), 40u32)];

    let fragments = vec![
        FragmentSpec {
            forward: ReadSpec {
                tid: 0,
                pos: 0,
                cigar: vec![('M', 10)],
                seq: seq(10, b'A'),
                qual: 40,
                is_reverse: false,
                mapq: 60,
                flags: FLAG_FIRST_MATE | FLAG_MATE_REVERSE | FLAG_PROPER_PAIR,
                mate_tid: Some(0),
                mate_pos: Some(14),
                insert_size: 24,
            },
            reverse: ReadSpec {
                tid: 0,
                pos: 14,
                cigar: vec![('M', 10)],
                seq: seq(10, b'T'),
                qual: 40,
                is_reverse: true,
                mapq: 60,
                flags: FLAG_SECOND_MATE | FLAG_PROPER_PAIR,
                mate_tid: Some(0),
                mate_pos: Some(0),
                insert_size: -24,
            },
        },
        FragmentSpec {
            forward: ReadSpec {
                tid: 0,
                pos: 5,
                cigar: vec![('M', 4), ('I', 1), ('M', 4)],
                seq: seq(9, b'C'),
                qual: 35,
                is_reverse: false,
                mapq: 55,
                flags: FLAG_FIRST_MATE | FLAG_MATE_REVERSE | FLAG_PROPER_PAIR,
                mate_tid: Some(0),
                mate_pos: Some(13),
                insert_size: 16,
            },
            reverse: ReadSpec {
                tid: 0,
                pos: 13,
                cigar: vec![('M', 8)],
                seq: seq(8, b'G'),
                qual: 35,
                is_reverse: true,
                mapq: 55,
                flags: FLAG_SECOND_MATE | FLAG_PROPER_PAIR,
                mate_tid: Some(0),
                mate_pos: Some(5),
                insert_size: -16,
            },
        },
        FragmentSpec {
            forward: ReadSpec {
                tid: 0,
                pos: 16,
                cigar: vec![('M', 3), ('D', 1), ('M', 5)],
                seq: seq(8, b'A'),
                qual: 30,
                is_reverse: false,
                mapq: 50,
                flags: FLAG_FIRST_MATE | FLAG_MATE_REVERSE | FLAG_PROPER_PAIR,
                mate_tid: Some(0),
                mate_pos: Some(20),
                insert_size: 11,
            },
            reverse: ReadSpec {
                tid: 0,
                pos: 20,
                cigar: vec![('M', 7)],
                seq: seq(7, b'T'),
                qual: 30,
                is_reverse: true,
                mapq: 50,
                flags: FLAG_SECOND_MATE | FLAG_PROPER_PAIR,
                mate_tid: Some(0),
                mate_pos: Some(16),
                insert_size: -11,
            },
        },
    ];

    bam_from_specs(chroms, fragments, Vec::new(), "fragment_kmers_edge")
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
        let qname = default_fragment_qname(fragment);
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

fn default_fragment_qname(fragment: &FragmentSpec) -> String {
    format!("frag{}_{}", fragment.forward.tid, fragment.forward.pos)
}

fn ensure_default_fragment_qnames_are_unique(fragments: &[FragmentSpec]) -> Result<()> {
    let mut seen = fxhash::FxHashSet::default();
    seen.reserve(fragments.len());
    for fragment in fragments {
        let qname = default_fragment_qname(fragment);
        ensure!(
            seen.insert(qname.clone()),
            "bam_from_specs would assign duplicate paired-read qname '{qname}'. \
             Duplicate qnames can make synthetic paired fragments overwrite each other in the \
             pairing stash. Use bam_from_specs_strict_identity when stacked fragments at the \
             same start are intentional."
        );
    }
    Ok(())
}

fn write_bam_with_strict_identity(
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

    // Give each synthetic fragment a unique qname so stacked fragments at the same genomic start
    // still represent distinct molecules to paired-end consumers.
    for (fragment_idx, fragment) in fragments.iter().enumerate() {
        let qname = format!(
            "frag{}_tid{}_pos{}",
            fragment_idx, fragment.forward.tid, fragment.forward.pos
        );
        records.push(fragment.forward.to_record(qname.as_bytes()));
        records.push(fragment.reverse.to_record(qname.as_bytes()));
    }

    // Singles need the same treatment for consistency if a test ever stacks reads at one start.
    for (single_idx, single) in singles.iter().enumerate() {
        let qname = format!("single{}_tid{}_pos{}", single_idx, single.tid, single.pos);
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
    ensure_default_fragment_qnames_are_unique(&fragments)?;

    let tempdir = TempDir::new()?;
    let bam_path = tempdir.path().join(format!("{name}.bam"));

    write_bam(&chroms, &fragments, &singles, &bam_path)?;
    let bai = build_index(&bam_path)?;
    Ok(BamFixture::new(tempdir, bam_path, bai))
}

pub fn single_read_bam_with_qualities(
    name: &str,
    pos: i64,
    cigar_ops: Vec<(char, u32)>,
    seq: &[u8],
    qualities: &[u8],
) -> Result<BamFixture> {
    if seq.len() != qualities.len() {
        return Err(anyhow!(
            "seq length ({}) must match qualities length ({})",
            seq.len(),
            qualities.len()
        ));
    }

    let tempdir = TempDir::new()?;
    let bam_path = tempdir.path().join(format!("{name}.bam"));
    let chrom_len = (pos.max(0) as u32)
        .saturating_add(seq.len() as u32)
        .saturating_add(100)
        .max(256);

    let mut header = bam::Header::new();
    header.push_record(
        HeaderRecord::new(b"HD")
            .push_tag(b"VN", &"1.6")
            .push_tag(b"SO", &"coordinate"),
    );
    header.push_record(
        HeaderRecord::new(b"SQ")
            .push_tag(b"SN", &"chr1")
            .push_tag(b"LN", chrom_len),
    );

    let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam)
        .with_context(|| format!("create bam at {}", bam_path.display()))?;

    let mut record = bam::Record::new();
    record.set_tid(0);
    record.set_pos(pos);
    record.set_mapq(60);
    record.set(
        b"single_custom_qualities",
        Some(&cigar(&cigar_ops)),
        seq,
        qualities,
    );
    record.set_flags(0);
    writer.write(&record)?;

    drop(writer);
    let bai = build_index(&bam_path)?;
    Ok(BamFixture::new(tempdir, bam_path, bai))
}

pub fn bam_from_specs_strict_identity(
    chroms: Vec<(String, u32)>,
    fragments: Vec<FragmentSpec>,
    singles: Vec<ReadSpec>,
    name: &str,
) -> Result<BamFixture> {
    let tempdir = TempDir::new()?;
    let bam_path = tempdir.path().join(format!("{name}.bam"));

    write_bam_with_strict_identity(&chroms, &fragments, &singles, &bam_path)?;
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
    writeln!(f, "chromosome\tstart\tend\tscaling_factor")?;
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

/// Read a zstd-compressed length-count TSV as text.
///
/// `cfdna lengths` writes zstd-compressed TSV output. Tests use this helper
/// when they need to inspect metadata columns directly.
pub fn read_length_counts_text<P: AsRef<Path>>(path: P) -> Result<String> {
    read_zst_to_string(path.as_ref())
}

/// Read a length-count TSV and return only the numeric `count_*` columns.
///
/// This keeps integration tests that previously read the old dense NPY matrix
/// focused on the same numeric contract while allowing metadata columns to
/// vary by output mode.
pub fn read_length_counts_tsv<P: AsRef<Path>>(path: P) -> Result<Array2<f64>> {
    let text = read_length_counts_text(path)?;
    let mut lines = text.lines();
    let header = lines
        .next()
        .context("length counts TSV must have a header")?;
    let headers: Vec<&str> = header.split('\t').collect();
    let first_count_column = headers
        .iter()
        .position(|column| column.starts_with("count_"))
        .context("length counts TSV must contain count columns")?;
    let count_column_count = headers.len() - first_count_column;
    ensure!(
        count_column_count > 0,
        "length counts TSV must contain at least one count column"
    );

    let mut values = Vec::new();
    let mut row_count = 0;
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        ensure!(
            fields.len() == headers.len(),
            "length counts row must match the header column count"
        );
        for value in &fields[first_count_column..] {
            values.push(value.parse::<f64>()?);
        }
        row_count += 1;
    }

    Ok(Array2::from_shape_vec(
        (row_count, count_column_count),
        values,
    )?)
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
