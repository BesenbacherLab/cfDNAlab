//! Builders for small synthetic BAM inputs.
//!
//! Use this module when a test needs an indexed BAM with known fragment spans.
//! The builders write coordinate-sorted BAM files and BAI indices in
//! a temporary directory. A returned `TempBam` owns that directory, so its paths
//! remain valid while the value is alive.
//!
//! Coordinates follow BAM conventions. `ReadSpec::pos` is zero-based. A paired
//! fragment span is directional and runs from the forward read `pos` to the
//! reverse read `reference_end`. This matches the fragment vocabulary used by
//! cfDNAlab commands and keeps expected fragment lengths derivable from the
//! fixture definition.
//!
//! The public surface favors explicit builders over vague canned fixtures.
//! Use `TempBamBuilder` for most cases. The named helpers are narrow scenarios
//! with documented coordinates and fragment spans.

use anyhow::{Context, Result, anyhow, ensure};
pub use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{self, header::HeaderRecord, record::CigarString};
use std::path::{Path, PathBuf};
use tempfile::TempDir;

/// SAM flag for the first read in a template.
pub const SAM_FLAG_FIRST_MATE: u16 = 0x40;
/// SAM flag for the second read in a template.
pub const SAM_FLAG_SECOND_MATE: u16 = 0x80;
/// SAM flag for a proper pair.
pub const SAM_FLAG_PROPER_PAIR: u16 = 0x2;
/// SAM flag indicating the mate is aligned to the reverse strand.
pub const SAM_FLAG_MATE_REVERSE: u16 = 0x20;

const SAM_FLAG_PAIRED: u16 = 0x1;
const SAM_FLAG_REVERSE: u16 = 0x10;
const MIN_FRAGMENT_LENGTH: i64 = 10;

/// Indexed BAM stored in an owned temporary directory.
///
/// `TempBam` is the return type for BAM builders in this module. It owns the
/// temporary directory containing both the BAM file and its BAI index. Keep the
/// value alive for as long as code under test needs to open either path.
///
/// The files are removed when the value is dropped. Clone the paths if another
/// API needs owned `PathBuf`s, but do not drop `TempBam` before those paths are
/// used.
#[derive(Debug)]
pub struct TempBam {
    _tempdir: TempDir,
    /// Path to the generated coordinate-sorted BAM file.
    pub bam: PathBuf,
    /// Path to the generated BAI index.
    pub bai: PathBuf,
}

impl TempBam {
    fn new(tempdir: TempDir, bam: PathBuf, bai: PathBuf) -> Self {
        Self {
            _tempdir: tempdir,
            bam,
            bai,
        }
    }

    /// Return the generated coordinate-sorted BAM path.
    pub fn bam_path(&self) -> &Path {
        &self.bam
    }

    /// Return the generated BAI path.
    pub fn bai_path(&self) -> &Path {
        &self.bai
    }
}

/// Explicit read record used to build a synthetic BAM.
///
/// Use `ReadSpec` when a test needs exact control over CIGAR operations, mate
/// fields, flags, mapping quality, or strand. For ordinary inward-facing pairs,
/// prefer `PairedFragmentSpec`, which derives both mates from fragment span
/// fields.
///
/// `tid` is the zero-based index into the contigs added to `TempBamBuilder`.
/// `pos` is the zero-based leftmost alignment position. `seq` must have the
/// same length as the query-consuming CIGAR operations, and all bases must be
/// `A`, `C`, `G`, `T`, or `N`.
///
/// When the record is written, the builder always adds the SAM paired flag
/// `0x1`. It also adds the SAM reverse-strand flag `0x10` when `is_reverse` is
/// true. Other flags come from `flags` unchanged, so callers should put mate
/// role, proper-pair, mate-reverse, duplicate, secondary, supplementary, or
/// other test-specific flags there.
#[derive(Clone, Debug)]
pub struct ReadSpec {
    /// Zero-based target id in the BAM header.
    pub tid: usize,
    /// Zero-based leftmost alignment position.
    pub pos: i64,
    /// CIGAR operations for the alignment.
    pub cigar: Vec<Cigar>,
    /// Read sequence.
    pub seq: Vec<u8>,
    /// Per-base quality value written uniformly for the read.
    pub base_quality: u8,
    /// Whether this read is aligned to the reverse strand.
    pub is_reverse: bool,
    /// Mapping quality.
    pub mapq: u8,
    /// Additional SAM flags. The builder adds paired and reverse-strand flags.
    pub flags: u16,
    /// Mate target id.
    pub mate_tid: Option<usize>,
    /// Mate leftmost alignment position.
    pub mate_pos: Option<i64>,
    /// Signed template length.
    pub insert_size: i64,
}

impl ReadSpec {
    fn to_record(&self, qname: &[u8]) -> bam::Record {
        let mut record = bam::Record::new();
        record.set_tid(self.tid as i32);
        record.set_pos(self.pos);
        record.set_mapq(self.mapq);
        if let Some(mate_tid) = self.mate_tid {
            record.set_mtid(mate_tid as i32);
        }
        if let Some(mate_pos) = self.mate_pos {
            record.set_mpos(mate_pos);
        }
        record.set_insert_size(self.insert_size);
        record.set(
            qname,
            Some(&cigar(&self.cigar)),
            &self.seq,
            &vec![self.base_quality; self.seq.len()],
        );
        let mut flags = self.flags | SAM_FLAG_PAIRED;
        if self.is_reverse {
            flags |= SAM_FLAG_REVERSE;
        }
        record.set_flags(flags);
        record
    }

    fn reference_end(&self) -> Result<i64> {
        let reference_len: u32 = self.cigar.iter().map(cigar_reference_len).sum();
        Ok(self.pos + i64::from(reference_len))
    }

    fn query_len(&self) -> usize {
        self.cigar
            .iter()
            .map(|op| cigar_query_len(op) as usize)
            .sum()
    }
}

/// Paired-end fragment represented by forward and reverse read records.
///
/// This type is useful when a test needs to control each mate directly. The
/// ordinary paired-end helper uses a forward read and a reverse read, but this
/// struct does not enforce orientation by itself. Validation happens when the
/// BAM is built.
///
/// `TempBamBuilder` validates each read independently: target id exists,
/// position is non-negative, CIGAR is non-empty, query-consuming CIGAR length
/// matches sequence length, sequence bases are supported, the read reference
/// span fits within its contig, and mate target/position fields point to valid
/// coordinates when present. It does not check that the two records are a
/// biologically consistent pair. Use `PairedFragmentSpec` when the test wants
/// the builder to derive a conventional inward-facing pair.
#[derive(Clone, Debug)]
pub struct FragmentSpec {
    /// Forward read record.
    pub forward: ReadSpec,
    /// Reverse read record.
    pub reverse: ReadSpec,
}

/// Compact specification for a conventional inward-facing paired fragment.
///
/// This is the preferred way to create ordinary paired fragments. It derives
/// both reads from a fragment start, a fragment length, and a read length:
///
/// - qname is assigned later by the BAM builder
/// - both reads use `tid`
/// - forward read `pos` is `start`
/// - forward CIGAR is `<read_length>M`
/// - forward sequence is `read_length` copies of `forward_base`
/// - forward MAPQ is `mapq`
/// - forward base quality is `base_quality`
/// - forward flags include first mate, mate reverse, proper pair, and paired
/// - forward mate position is `start + fragment_length - read_length`
/// - forward insert size is `fragment_length`
/// - reverse read `pos` is `start + fragment_length - read_length`
/// - reverse CIGAR is `<read_length>M`
/// - reverse sequence is `read_length` copies of `reverse_base`
/// - reverse MAPQ is `mapq`
/// - reverse base quality is `base_quality`
/// - reverse flags include second mate, proper pair, reverse strand, and paired
/// - reverse mate position is `start`
/// - reverse insert size is `-fragment_length`
///
/// Fragment length must be at least 10 bp, matching the minimum supported
/// fragment length used by cfDNAlab fixtures. Read length must be positive and
/// cannot exceed fragment length.
///
/// The directional fragment span is `[start, start + fragment_length)`. This
/// span is also `forward.pos` to `reverse.reference_end`.
#[derive(Clone, Debug)]
pub struct PairedFragmentSpec {
    /// Zero-based target id in the BAM header.
    pub tid: usize,
    /// Fragment start, corresponding to the forward read `pos`.
    pub start: i64,
    /// Directional fragment length from forward `pos` to reverse
    /// `reference_end`.
    pub fragment_length: i64,
    /// Aligned read length for both mates.
    pub read_length: i64,
    /// Mapping quality for both reads.
    pub mapq: u8,
    /// Per-base quality value for both reads.
    pub base_quality: u8,
    /// Base used for the forward-read sequence.
    pub forward_base: u8,
    /// Base used for the reverse-read sequence.
    pub reverse_base: u8,
}

impl PairedFragmentSpec {
    /// Create a paired-fragment spec with MAPQ 60, base quality 40, forward
    /// `A` bases, and reverse `T` bases.
    ///
    /// Use the setter methods to override mapping quality, base quality, or
    /// repeated read bases before passing the spec to a builder.
    pub fn new(tid: usize, start: i64, fragment_length: i64, read_length: i64) -> Self {
        Self {
            tid,
            start,
            fragment_length,
            read_length,
            mapq: 60,
            base_quality: 40,
            forward_base: b'A',
            reverse_base: b'T',
        }
    }

    /// Set the mapping quality for both reads.
    pub fn mapq(mut self, mapq: u8) -> Self {
        self.mapq = mapq;
        self
    }

    /// Set the per-base quality value for both reads.
    pub fn base_quality(mut self, base_quality: u8) -> Self {
        self.base_quality = base_quality;
        self
    }

    /// Set the repeated sequence bases for the forward and reverse reads.
    pub fn bases(mut self, forward_base: u8, reverse_base: u8) -> Self {
        self.forward_base = forward_base;
        self.reverse_base = reverse_base;
        self
    }

    /// Build the explicit paired-read records.
    pub fn build(&self) -> Result<FragmentSpec> {
        paired_fragment(self)
    }
}

/// Read-name behavior for generated paired fragments.
///
/// By default, generated paired fragments use a name derived from target id and
/// start position. That catches accidental duplicate molecules at the same
/// position. Use `UseRecordIndex` only when the test intentionally stacks
/// multiple fragments at the same start.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReadNamePolicy {
    /// Reject fragments that would receive the same synthetic read name.
    RejectDuplicates,
    /// Include the record index in synthetic read names so stacked fragments at
    /// the same position remain distinct molecules.
    UseRecordIndex,
}

/// Builder for indexed temporary BAM files.
///
/// `TempBamBuilder` writes a real BAM file and BAI index. It validates contig
/// definitions, read coordinates, CIGAR/query length agreement, sequence bases,
/// mate coordinates, and duplicate generated read names.
///
/// The builder starts empty. Add at least one contig, then add explicit
/// fragments, compact paired-fragment specs, or single reads. The result is a
/// `TempBam` that owns the generated files.
///
/// Records are sorted by `(tid, pos)` before writing. When
/// `ReadNamePolicy::RejectDuplicates` is active, paired-fragment qnames are
/// `frag{tid}_{pos}` from the forward read and duplicate generated qnames are
/// rejected. When `ReadNamePolicy::UseRecordIndex` is active, paired-fragment
/// qnames are `frag{record_index}_tid{tid}_pos{pos}`, allowing stacked
/// fragments at the same start. Single-read qnames are generated similarly,
/// using either `single{tid}_{pos}` or `single{record_index}_tid{tid}_pos{pos}`.
#[derive(Clone, Debug)]
pub struct TempBamBuilder {
    file_stem: String,
    contigs: Vec<(String, u32)>,
    fragments: Vec<FragmentSpec>,
    paired_fragments: Vec<PairedFragmentSpec>,
    single_reads: Vec<ReadSpec>,
    read_name_policy: ReadNamePolicy,
}

impl Default for TempBamBuilder {
    fn default() -> Self {
        Self {
            file_stem: "temp_bam".to_string(),
            contigs: Vec::new(),
            fragments: Vec::new(),
            paired_fragments: Vec::new(),
            single_reads: Vec::new(),
            read_name_policy: ReadNamePolicy::RejectDuplicates,
        }
    }
}

impl TempBamBuilder {
    /// Create an empty BAM builder.
    ///
    /// Add contigs and reads before calling `build`. The default generated
    /// file stem is `temp_bam`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the generated BAM file stem.
    ///
    /// The stem is used for the BAM file name inside the temporary directory.
    /// Passing `"sample"` creates `sample.bam` and a matching BAI path.
    pub fn name(mut self, file_stem: impl Into<String>) -> Self {
        self.file_stem = file_stem.into();
        self
    }

    /// Add a BAM contig with the given header name and length.
    ///
    /// Read `tid` values refer to the order contigs are added here. Contig
    /// length must be positive, and all reads must fit within their target
    /// contig after CIGAR reference-consuming operations are applied.
    pub fn contig(mut self, name: impl Into<String>, length: u32) -> Self {
        self.contigs.push((name.into(), length));
        self
    }

    /// Add an explicit paired fragment.
    ///
    /// Use this when both mates need custom CIGARs, flags, mate positions, or
    /// sequence content.
    pub fn fragment(mut self, fragment: FragmentSpec) -> Self {
        self.fragments.push(fragment);
        self
    }

    /// Add a conventional inward-facing paired fragment.
    ///
    /// This is the shortest path for ordinary paired-end molecules. The full
    /// read records are derived when `build` is called.
    pub fn paired_fragment(mut self, fragment: PairedFragmentSpec) -> Self {
        self.paired_fragments.push(fragment);
        self
    }

    /// Add a single read record.
    ///
    /// The read is written with the paired SAM flag because this builder is
    /// mainly used for cfDNAlab paired-end fixture inputs. Use explicit flags if
    /// a test needs unusual records.
    pub fn single_read(mut self, read: ReadSpec) -> Self {
        self.single_reads.push(read);
        self
    }

    /// Set how synthetic read names are assigned.
    ///
    /// `RejectDuplicates` is the default and is better for most tests because
    /// it catches accidental duplicate molecule names.
    pub fn read_name_policy(mut self, policy: ReadNamePolicy) -> Self {
        self.read_name_policy = policy;
        self
    }

    /// Use unique record-indexed read names.
    ///
    /// Choose this when stacked fragments at the same genomic start are
    /// intentional and should represent distinct molecules to paired-end
    /// consumers.
    pub fn use_record_indexed_read_names(mut self) -> Self {
        self.read_name_policy = ReadNamePolicy::UseRecordIndex;
        self
    }

    /// Write the BAM and BAI files into a temporary directory.
    ///
    /// Returns a `TempBam` that owns the temporary directory. The BAM is
    /// coordinate sorted before writing and indexed with a BAI index.
    pub fn build(self) -> Result<TempBam> {
        ensure!(
            !self.file_stem.is_empty(),
            "temporary BAM name must not be empty"
        );
        ensure!(
            !self.contigs.is_empty(),
            "temporary BAM must contain at least one contig"
        );
        for (name, length) in &self.contigs {
            ensure!(!name.is_empty(), "BAM contig names must not be empty");
            ensure!(*length > 0, "BAM contig {name} must have positive length");
        }

        let mut fragments = self.fragments;
        for fragment in &self.paired_fragments {
            fragments.push(fragment.build()?);
        }
        validate_records(&self.contigs, &fragments, &self.single_reads)?;
        if self.read_name_policy == ReadNamePolicy::RejectDuplicates {
            ensure_default_fragment_qnames_are_unique(&fragments)?;
        }

        let tempdir = TempDir::new()?;
        let bam_path = tempdir.path().join(format!("{}.bam", self.file_stem));
        match self.read_name_policy {
            ReadNamePolicy::RejectDuplicates => {
                write_bam(&self.contigs, &fragments, &self.single_reads, &bam_path)?
            }
            ReadNamePolicy::UseRecordIndex => write_bam_with_record_indexed_names(
                &self.contigs,
                &fragments,
                &self.single_reads,
                &bam_path,
            )?,
        }
        let bai_path = build_index(&bam_path)?;
        Ok(TempBam::new(tempdir, bam_path, bai_path))
    }
}

/// Build explicit read records for a conventional inward-facing pair.
///
/// Use this when a caller wants to inspect or modify the generated reads before
/// passing them to `TempBamBuilder`. Most tests can pass `PairedFragmentSpec`
/// directly to `TempBamBuilder::paired_fragment`.
///
/// The generated records follow the `PairedFragmentSpec` contract exactly:
/// the forward read starts at `start`, the reverse read starts at
/// `start + fragment_length - read_length`, both reads use `<read_length>M`,
/// and the directional fragment span is `[start, start + fragment_length)`.
/// The returned reads are not written to disk until they are passed to a BAM
/// builder.
pub fn paired_fragment(spec: &PairedFragmentSpec) -> Result<FragmentSpec> {
    ensure!(
        spec.start >= 0,
        "paired fragment start must be non-negative, got {}",
        spec.start
    );
    ensure!(
        spec.fragment_length >= MIN_FRAGMENT_LENGTH,
        "fragment length must be at least {MIN_FRAGMENT_LENGTH} bp, got {}",
        spec.fragment_length
    );
    ensure!(
        spec.read_length > 0,
        "read length must be positive, got {}",
        spec.read_length
    );
    ensure!(
        spec.read_length <= spec.fragment_length,
        "read length ({}) must not exceed fragment length ({})",
        spec.read_length,
        spec.fragment_length
    );
    ensure_valid_base(spec.forward_base)?;
    ensure_valid_base(spec.reverse_base)?;

    let reverse_start = spec.start + spec.fragment_length - spec.read_length;
    let forward = ReadSpec {
        tid: spec.tid,
        pos: spec.start,
        cigar: vec![Cigar::Match(spec.read_length as u32)],
        seq: repeated_base(spec.read_length, spec.forward_base)?,
        base_quality: spec.base_quality,
        is_reverse: false,
        mapq: spec.mapq,
        flags: SAM_FLAG_FIRST_MATE | SAM_FLAG_MATE_REVERSE | SAM_FLAG_PROPER_PAIR,
        mate_tid: Some(spec.tid),
        mate_pos: Some(reverse_start),
        insert_size: spec.fragment_length,
    };
    let reverse = ReadSpec {
        tid: spec.tid,
        pos: reverse_start,
        cigar: vec![Cigar::Match(spec.read_length as u32)],
        seq: repeated_base(spec.read_length, spec.reverse_base)?,
        base_quality: spec.base_quality,
        is_reverse: true,
        mapq: spec.mapq,
        flags: SAM_FLAG_SECOND_MATE | SAM_FLAG_PROPER_PAIR,
        mate_tid: Some(spec.tid),
        mate_pos: Some(spec.start),
        insert_size: -spec.fragment_length,
    };
    Ok(FragmentSpec { forward, reverse })
}

/// Create a temporary BAM from explicit paired fragments and single reads.
///
/// This is a convenience wrapper around `TempBamBuilder` for callers that
/// already have all records in vectors. It uses `ReadNamePolicy::RejectDuplicates`.
///
/// Records are coordinate sorted before writing. Paired-fragment qnames are
/// derived from the forward read as `frag{tid}_{pos}`. If two fragments would
/// receive the same generated qname, this function returns an error instead of
/// silently producing a BAM where paired-end consumers may collapse records.
///
/// The returned `TempBam` owns a temporary directory containing `<name>.bam`
/// and a matching BAI index.
pub fn bam_from_fragments(
    name: &str,
    contigs: Vec<(String, u32)>,
    fragments: Vec<FragmentSpec>,
    single_reads: Vec<ReadSpec>,
) -> Result<TempBam> {
    let mut builder = TempBamBuilder::new()
        .name(name)
        .read_name_policy(ReadNamePolicy::RejectDuplicates);
    for (contig_name, contig_length) in contigs {
        builder = builder.contig(contig_name, contig_length);
    }
    for fragment in fragments {
        builder = builder.fragment(fragment);
    }
    for read in single_reads {
        builder = builder.single_read(read);
    }
    builder.build()
}

/// Create a temporary BAM from explicit records using record-indexed read names.
///
/// Use this when intentionally stacking multiple fragments at the same start.
/// Record-indexed names keep those fragments distinct for paired-end consumers.
///
/// Records are coordinate sorted before writing. Paired-fragment qnames are
/// `frag{record_index}_tid{tid}_pos{pos}`, using the forward read's target id
/// and position. Single-read qnames are
/// `single{record_index}_tid{tid}_pos{pos}`.
///
/// The returned `TempBam` owns a temporary directory containing `<name>.bam`
/// and a matching BAI index.
pub fn bam_from_fragments_with_record_indexed_names(
    name: &str,
    contigs: Vec<(String, u32)>,
    fragments: Vec<FragmentSpec>,
    single_reads: Vec<ReadSpec>,
) -> Result<TempBam> {
    let mut builder = TempBamBuilder::new()
        .name(name)
        .use_record_indexed_read_names();
    for (contig_name, contig_length) in contigs {
        builder = builder.contig(contig_name, contig_length);
    }
    for fragment in fragments {
        builder = builder.fragment(fragment);
    }
    for read in single_reads {
        builder = builder.single_read(read);
    }
    builder.build()
}

/// Create a single-contig temporary BAM from paired fragment starts.
///
/// The generated BAM has one contig named `chr1`. Its length is derived from
/// the maximum fragment end plus padding, with a minimum contig length of 1000
/// bp. Each start creates one inward-facing pair with the supplied fragment and
/// read lengths.
///
/// For each `start`, this helper creates the same records as
/// `PairedFragmentSpec::new(0, start, fragment_length, read_length)`:
///
/// - qname: `frag0_{start}`
/// - contig: `chr1`
/// - forward read `pos`: `start`
/// - reverse read `pos`: `start + fragment_length - read_length`
/// - both CIGARs: `<read_length>M`
/// - forward sequence: `read_length` `A` bases
/// - reverse sequence: `read_length` `T` bases
/// - MAPQ: 60
/// - base quality: 40
/// - forward insert size: `fragment_length`
/// - reverse insert size: `-fragment_length`
/// - directional span: `[start, start + fragment_length)`
///
/// The contig length is `max(max_start + fragment_length + 500, 1000)`. If
/// `starts` is empty, the BAM contains no records and `chr1` length is 1000.
/// Read names use `ReadNamePolicy::RejectDuplicates`, so duplicate starts are
/// rejected.
pub fn bam_from_fragment_starts(
    name: &str,
    starts: &[i64],
    fragment_length: i64,
    read_length: i64,
) -> Result<TempBam> {
    let max_end = starts
        .iter()
        .map(|start| start + fragment_length)
        .max()
        .unwrap_or(0);
    let contig_length =
        u32::try_from((max_end + 500).max(1_000)).context("derived contig length exceeds u32")?;
    let mut builder = TempBamBuilder::new()
        .name(name)
        .read_name_policy(ReadNamePolicy::RejectDuplicates)
        .contig("chr1", contig_length);
    for &start in starts {
        builder = builder.paired_fragment(PairedFragmentSpec::new(
            0,
            start,
            fragment_length,
            read_length,
        ));
    }
    builder.build()
}

/// Create a temporary BAM with one `chr1` inward-facing pair spanning `[20, 80)`.
///
/// Use this when a test needs one clean paired-end molecule with easy
/// hand-derived coverage. The file contains one contig and one paired record
/// set:
///
/// - Contig:
///   - name: `chr1`
///   - length: 200 bp
///
/// - Forward read:
///   - qname: `frag0_20`
///   - `tid`: 0
///   - `pos`: 20
///   - CIGAR: `20M`
///   - query length: 20 bp
///   - reference-consuming length: 20 bp
///   - `reference_end`: 40
///   - sequence: 20 `A` bases
///   - MAPQ: 60
///   - base quality: 35
///   - flags: first mate, mate reverse, proper pair, paired
///   - mate position: 60
///   - insert size: 60
///
/// - Reverse read:
///   - qname: `frag0_20`
///   - `tid`: 0
///   - `pos`: 60
///   - CIGAR: `20M`
///   - query length: 20 bp
///   - reference-consuming length: 20 bp
///   - `reference_end`: 80
///   - sequence: 20 `T` bases
///   - MAPQ: 60
///   - base quality: 35
///   - flags: second mate, proper pair, reverse strand, paired
///   - mate position: 20
///   - insert size: -60
///
/// The directional span from the forward read `pos` to the reverse read
/// `reference_end` is `[20, 80)`, or 60 bp. Read names use
/// `ReadNamePolicy::RejectDuplicates`.
pub fn single_contig_inward_pair_bam() -> Result<TempBam> {
    TempBamBuilder::new()
        .name("single_contig_inward_pair")
        .read_name_policy(ReadNamePolicy::RejectDuplicates)
        .contig("chr1", 200)
        .paired_fragment(PairedFragmentSpec::new(0, 20, 60, 20).base_quality(35))
        .build()
}

/// Create a temporary BAM with ten 600 bp inward-facing fragments on one contig.
///
/// Fragment starts are fixed at `0, 400, ..., 3600`, read length is 100 bp,
/// and the contig is named `chr1`. This is useful for tests that need multiple
/// long fragments with easy hand-derived spans.
///
/// The generated `chr1` contig length is 4700 bp. For each start, the forward
/// read is at `start`, the reverse read is at `start + 500`, and the
/// directional span is `[start, start + 600)`. Both reads use `100M`, MAPQ 60,
/// and base quality 40. Forward reads contain 100 `A` bases and reverse reads
/// contain 100 `T` bases.
///
/// Qnames are `frag0_{start}`. Forward insert size is 600 and reverse insert
/// size is -600. Read names use `ReadNamePolicy::RejectDuplicates`.
pub fn long_inward_fragment_series_bam(name: &str) -> Result<TempBam> {
    const STARTS: [i64; 10] = [0, 400, 800, 1_200, 1_600, 2_000, 2_400, 2_800, 3_200, 3_600];
    bam_from_fragment_starts(name, &STARTS, 600, 100)
}

/// Create a temporary BAM whose paired reads include non-trivial CIGAR operations.
///
/// Use this when a test needs a compact indexed BAM that exercises
/// query-consuming and reference-consuming CIGAR operations. The file contains
/// one contig and one paired record set:
///
/// - Contig:
///   - name: `chr1`
///   - length: 500 bp
///
/// - Forward read:
///   - `tid`: 0
///   - `pos`: 200
///   - CIGAR: `20M3I10M5D12M`
///   - query length: 45 bp
///   - reference-consuming length: 47 bp
///   - `reference_end`: 247
///   - sequence: 45 `C` bases
///   - MAPQ: 50
///   - base quality: 25
///   - flags: first mate, mate reverse, proper pair
///   - mate position: 260
///   - insert size: 105
///
/// - Reverse read:
///   - `tid`: 0
///   - `pos`: 260
///   - CIGAR: `5S25M4N16M`
///   - query length: 46 bp
///   - reference-consuming length: 45 bp
///   - `reference_end`: 305
///   - sequence: 46 `G` bases
///   - MAPQ: 40
///   - base quality: 25
///   - flags: second mate, proper pair, reverse strand
///   - mate position: 200
///   - insert size: -105
///
/// The directional span from the forward read `pos` to the reverse read
/// `reference_end` is `[200, 305)`, or 105 bp. The SAM template length fields
/// use the same length.
pub fn bam_with_indel_and_softclip_reads() -> Result<TempBam> {
    let fragment = FragmentSpec {
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
            seq: repeated_base(45, b'C')?,
            base_quality: 25,
            is_reverse: false,
            mapq: 50,
            flags: SAM_FLAG_FIRST_MATE | SAM_FLAG_MATE_REVERSE | SAM_FLAG_PROPER_PAIR,
            mate_tid: Some(0),
            mate_pos: Some(260),
            insert_size: 105,
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
            seq: repeated_base(46, b'G')?,
            base_quality: 25,
            is_reverse: true,
            mapq: 40,
            flags: SAM_FLAG_SECOND_MATE | SAM_FLAG_PROPER_PAIR,
            mate_tid: Some(0),
            mate_pos: Some(200),
            insert_size: -105,
        },
    };
    bam_from_fragments(
        "indel_and_softclip_reads",
        vec![("chr1".to_string(), 500)],
        vec![fragment],
        Vec::new(),
    )
}

/// Create a temporary one-read BAM with caller-supplied base qualities.
///
/// The BAM has one contig named `chr1` and one unpaired read:
///
/// - qname: `single_custom_qualities`
/// - `tid`: 0
/// - `pos`: `pos`
/// - CIGAR: `cigar_ops`
/// - sequence: `seq`
/// - qualities: `qualities`
/// - MAPQ: 60
/// - flags: 0
///
/// The sequence length must match both the quality length and the
/// query-consuming CIGAR length.
///
/// Pass `Some(chr1_length)` when the surrounding test depends on the exact
/// chromosome extent, for example when deriving fixed-size output bins. The
/// supplied length must be positive and must contain the read's full
/// reference-consuming span. Pass `None` to derive a permissive contig length
/// as `max(pos + reference_consuming_cigar_length + 100, 256)`. Deletions and
/// reference skips are included in the reference-consuming length.
pub fn single_read_bam_with_qualities(
    name: &str,
    chr1_length: Option<u32>,
    pos: i64,
    cigar_ops: Vec<Cigar>,
    seq: &[u8],
    qualities: &[u8],
) -> Result<TempBam> {
    ensure!(
        seq.len() == qualities.len(),
        "seq length ({}) must match qualities length ({})",
        seq.len(),
        qualities.len()
    );
    ensure_valid_sequence(seq)?;
    ensure!(pos >= 0, "read position must be non-negative, got {pos}");
    ensure!(!cigar_ops.is_empty(), "read CIGAR must not be empty");
    let query_len: usize = cigar_ops
        .iter()
        .map(|op| cigar_query_len(op) as usize)
        .sum();
    ensure!(
        query_len == seq.len(),
        "read sequence length ({}) must match query-consuming CIGAR length ({query_len})",
        seq.len()
    );
    let reference_len: u32 = cigar_ops.iter().map(cigar_reference_len).sum();
    let pos_u32 = u32::try_from(pos).context("read position exceeds u32 coordinate range")?;
    let minimum_chr1_length = pos_u32
        .checked_add(reference_len)
        .context("read reference end exceeds u32 coordinate range")?;
    let chrom_len = match chr1_length {
        Some(length) => {
            ensure!(length > 0, "chr1 length must be positive, got {length}");
            ensure!(
                length >= minimum_chr1_length,
                "chr1 length ({length}) must contain read reference span ending at {minimum_chr1_length}"
            );
            length
        }
        None => minimum_chr1_length.saturating_add(100).max(256),
    };
    let tempdir = TempDir::new()?;
    let bam_path = tempdir.path().join(format!("{name}.bam"));

    let mut header = bam::Header::new();
    push_header(&mut header, &[("chr1".to_string(), chrom_len)]);

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

    let bai_path = build_index(&bam_path)?;
    Ok(TempBam::new(tempdir, bam_path, bai_path))
}

fn push_header(header: &mut bam::Header, contigs: &[(String, u32)]) {
    header.push_record(
        HeaderRecord::new(b"HD")
            .push_tag(b"VN", "1.6")
            .push_tag(b"SO", "coordinate"),
    );
    for (name, length) in contigs {
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", name)
                .push_tag(b"LN", *length),
        );
    }
}

fn cigar(ops: &[Cigar]) -> CigarString {
    CigarString(ops.to_vec())
}

fn cigar_reference_len(op: &Cigar) -> u32 {
    match *op {
        Cigar::Match(length)
        | Cigar::Equal(length)
        | Cigar::Diff(length)
        | Cigar::Del(length)
        | Cigar::RefSkip(length) => length,
        Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => 0,
    }
}

fn cigar_query_len(op: &Cigar) -> u32 {
    match *op {
        Cigar::Match(length)
        | Cigar::Equal(length)
        | Cigar::Diff(length)
        | Cigar::Ins(length)
        | Cigar::SoftClip(length) => length,
        Cigar::Del(_) | Cigar::RefSkip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => 0,
    }
}

fn validate_records(
    contigs: &[(String, u32)],
    fragments: &[FragmentSpec],
    single_reads: &[ReadSpec],
) -> Result<()> {
    for fragment in fragments {
        validate_read(contigs, &fragment.forward)?;
        validate_read(contigs, &fragment.reverse)?;
    }
    for read in single_reads {
        validate_read(contigs, read)?;
    }
    Ok(())
}

fn validate_read(contigs: &[(String, u32)], read: &ReadSpec) -> Result<()> {
    ensure!(
        read.tid < contigs.len(),
        "read tid {} has no matching BAM contig",
        read.tid
    );
    ensure!(
        read.pos >= 0,
        "read position must be non-negative, got {}",
        read.pos
    );
    ensure!(!read.cigar.is_empty(), "read CIGAR must not be empty");
    ensure!(
        read.query_len() == read.seq.len(),
        "read sequence length ({}) must match query-consuming CIGAR length ({})",
        read.seq.len(),
        read.query_len()
    );
    ensure_valid_sequence(&read.seq)?;
    let reference_end = read.reference_end()?;
    ensure!(
        reference_end <= i64::from(contigs[read.tid].1),
        "read reference_end ({reference_end}) exceeds contig {} length ({})",
        contigs[read.tid].0,
        contigs[read.tid].1
    );
    if let Some(mate_tid) = read.mate_tid {
        ensure!(
            mate_tid < contigs.len(),
            "read mate_tid {mate_tid} has no matching BAM contig"
        );
    }
    if let Some(mate_pos) = read.mate_pos {
        ensure!(
            mate_pos >= 0,
            "read mate position must be non-negative, got {mate_pos}"
        );
    }
    Ok(())
}

fn write_bam(
    contigs: &[(String, u32)],
    fragments: &[FragmentSpec],
    single_reads: &[ReadSpec],
    out_bam: &Path,
) -> Result<()> {
    let mut header = bam::Header::new();
    push_header(&mut header, contigs);
    let mut writer = bam::Writer::from_path(out_bam, &header, bam::Format::Bam)
        .with_context(|| format!("create bam at {}", out_bam.display()))?;

    let mut records = Vec::new();
    for fragment in fragments {
        let qname = default_fragment_qname(fragment);
        records.push(fragment.forward.to_record(qname.as_bytes()));
        records.push(fragment.reverse.to_record(qname.as_bytes()));
    }
    for read in single_reads {
        let qname = format!("single{}_{}", read.tid, read.pos);
        records.push(read.to_record(qname.as_bytes()));
    }
    write_sorted_records(&mut writer, records)
}

fn write_bam_with_record_indexed_names(
    contigs: &[(String, u32)],
    fragments: &[FragmentSpec],
    single_reads: &[ReadSpec],
    out_bam: &Path,
) -> Result<()> {
    let mut header = bam::Header::new();
    push_header(&mut header, contigs);
    let mut writer = bam::Writer::from_path(out_bam, &header, bam::Format::Bam)
        .with_context(|| format!("create bam at {}", out_bam.display()))?;

    let mut records = Vec::new();
    for (fragment_idx, fragment) in fragments.iter().enumerate() {
        let qname = format!(
            "frag{}_tid{}_pos{}",
            fragment_idx, fragment.forward.tid, fragment.forward.pos
        );
        records.push(fragment.forward.to_record(qname.as_bytes()));
        records.push(fragment.reverse.to_record(qname.as_bytes()));
    }
    for (read_idx, read) in single_reads.iter().enumerate() {
        let qname = format!("single{}_tid{}_pos{}", read_idx, read.tid, read.pos);
        records.push(read.to_record(qname.as_bytes()));
    }
    write_sorted_records(&mut writer, records)
}

fn write_sorted_records(writer: &mut bam::Writer, mut records: Vec<bam::Record>) -> Result<()> {
    records.sort_by_key(|record| (record.tid(), record.pos()));
    for record in records {
        writer.write(&record)?;
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
            "temporary BAM would assign duplicate paired-read qname '{qname}'. \
             Use ReadNamePolicy::UseRecordIndex when stacked fragments at the same start are intentional."
        );
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

fn repeated_base(length: i64, base: u8) -> Result<Vec<u8>> {
    ensure!(length >= 0, "sequence length must be non-negative");
    ensure_valid_base(base)?;
    Ok(std::iter::repeat_n(base, length as usize).collect())
}

fn ensure_valid_sequence(seq: &[u8]) -> Result<()> {
    ensure!(!seq.is_empty(), "read sequence must not be empty");
    for &base in seq {
        ensure_valid_base(base)?;
    }
    Ok(())
}

fn ensure_valid_base(base: u8) -> Result<()> {
    match base.to_ascii_uppercase() {
        b'A' | b'C' | b'G' | b'T' | b'N' => Ok(()),
        _ => Err(anyhow!(
            "synthetic sequence base must be A, C, G, T, or N, got {:?}",
            char::from(base)
        )),
    }
}
