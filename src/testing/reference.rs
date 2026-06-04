//! Builders for small two-bit reference inputs.
//!
//! Use this module when a test needs a real two-bit reference file with
//! compact, explicit sequence content. The builders write a temporary FASTA,
//! convert it to two-bit, and return a `TempTwoBit` that owns the generated
//! files.
//!
//! Sequence names and bases are validated before writing. Input sequences are
//! uppercased, and only `A`, `C`, `G`, `T`, and `N` are accepted. The returned
//! value also stores the normalized sequences so tests can derive expected GC
//! content, k-mers, or contig lengths without reopening the file.

use anyhow::{Result, anyhow, ensure};
use std::{
    fs::File,
    io::Write,
    path::{Path, PathBuf},
};
use tempfile::TempDir;
use twobit::convert::{fasta::FastaReader, to_2bit};

/// Two-bit reference stored in an owned temporary directory.
///
/// `TempTwoBit` owns the directory that contains the generated FASTA and
/// two-bit files. Keep the value alive while code under test needs the path.
/// The directory is removed when the value is dropped.
///
/// The `path` field points to the generated `<name>.2bit` file. The original
/// normalized contig sequences are available through `sequence` and
/// `sequences`, so tests can derive expected GC fractions, contig lengths, or
/// k-mers without parsing the two-bit file again.
#[derive(Debug)]
pub struct TempTwoBit {
    _tempdir: TempDir,
    /// Path to the generated two-bit file.
    pub path: PathBuf,
    sequences: Vec<(String, String)>,
}

impl TempTwoBit {
    fn new(tempdir: TempDir, path: PathBuf, sequences: Vec<(String, String)>) -> Self {
        Self {
            _tempdir: tempdir,
            path,
            sequences,
        }
    }

    /// Return the generated two-bit path.
    ///
    /// The path remains valid while this `TempTwoBit` value is alive.
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Return a generated contig sequence by name.
    ///
    /// The returned sequence is the normalized uppercase sequence that was
    /// written to the temporary FASTA before conversion.
    pub fn sequence(&self, contig_name: &str) -> Option<&str> {
        self.sequences
            .iter()
            .find(|(name, _)| name == contig_name)
            .map(|(_, sequence)| sequence.as_str())
    }

    /// Return all generated contig sequences.
    ///
    /// Sequences are returned in the order supplied to the builder.
    pub fn sequences(&self) -> &[(String, String)] {
        &self.sequences
    }
}

/// Specification for one repeating contig.
///
/// This is a compact way to generate references such as `ACGTACGT...` while
/// keeping the final contig length explicit. The repeated pattern is
/// uppercased and validated before sequence generation.
///
/// For example, `RepeatingContigSpec::new("chr1", "ACG", 8)` generates
/// `ACGACGAC`. The final repeat is truncated to exactly `length` bases.
#[derive(Clone, Debug)]
pub struct RepeatingContigSpec {
    /// Contig name.
    pub name: String,
    /// Repeated DNA pattern.
    pub pattern: String,
    /// Final contig length in bases.
    pub length: usize,
}

impl RepeatingContigSpec {
    /// Create a repeating contig spec.
    ///
    /// `length` is the final contig length, not the number of pattern repeats.
    /// If the requested length is not a multiple of the pattern length, the
    /// final repeat is truncated.
    pub fn new(name: impl Into<String>, pattern: impl Into<String>, length: usize) -> Self {
        Self {
            name: name.into(),
            pattern: pattern.into(),
            length,
        }
    }

    fn sequence(&self) -> Result<(String, String)> {
        ensure!(
            !self.name.is_empty(),
            "two-bit contig names must not be empty"
        );
        ensure!(
            !self.pattern.is_empty(),
            "repeating contig pattern must not be empty"
        );
        ensure!(
            self.length > 0,
            "repeating contig length must be greater than 0"
        );
        let pattern = normalize_sequence(&self.pattern)?;
        let pattern_bytes = pattern.as_bytes();
        let mut sequence = String::with_capacity(self.length);
        for position in 0..self.length {
            sequence.push(char::from(pattern_bytes[position % pattern_bytes.len()]));
        }
        Ok((self.name.clone(), sequence))
    }
}

/// Create a temporary two-bit reference from explicit contig sequences.
///
/// Sequences are uppercased and must contain only `A`, `C`, `G`, `T`, or `N`.
/// The returned `TempTwoBit` owns the temporary directory and stores the
/// normalized sequences for assertions.
///
/// The helper writes a temporary FASTA named `<name>.fasta`, then converts it
/// to `<name>.2bit`. FASTA sequence lines are wrapped at 60 bases. Contig order
/// is preserved exactly as supplied in `sequences`.
///
/// This helper does not invent contig sizes or padding. Each contig length is
/// exactly the length of the supplied sequence after uppercasing.
pub fn twobit_from_sequences(name: &str, sequences: Vec<(String, String)>) -> Result<TempTwoBit> {
    ensure!(!name.is_empty(), "temporary two-bit name must not be empty");
    ensure!(
        !sequences.is_empty(),
        "temporary two-bit reference must contain at least one contig"
    );
    let mut normalized = Vec::with_capacity(sequences.len());
    for (contig_name, sequence) in sequences {
        ensure!(
            !contig_name.is_empty(),
            "two-bit contig names must not be empty"
        );
        normalized.push((contig_name, normalize_sequence(&sequence)?));
    }

    let tempdir = TempDir::new()?;
    let fasta_path = tempdir.path().join(format!("{name}.fasta"));
    write_fasta(&fasta_path, &normalized)?;
    let path = tempdir.path().join(format!("{name}.2bit"));
    {
        let reader = FastaReader::open(&fasta_path).map_err(|error| anyhow!(error.to_string()))?;
        let mut file = File::create(&path)?;
        to_2bit(&mut file, &reader).map_err(|error| anyhow!(error.to_string()))?;
    }
    Ok(TempTwoBit::new(tempdir, path, normalized))
}

/// Create a temporary two-bit reference from repeating contig specs.
///
/// This helper is useful when expected GC content or k-mer content should be
/// obvious from a short pattern. Each `RepeatingContigSpec` defines one contig
/// and its final length.
///
/// Contigs are generated in the order supplied by `contigs`. Each pattern is
/// repeated from its first base and truncated to the requested final length.
/// The generated sequences are then passed through `twobit_from_sequences`, so
/// the same uppercase normalization, base validation, and temporary file
/// ownership rules apply.
pub fn twobit_with_repeating_contigs(
    name: &str,
    contigs: &[RepeatingContigSpec],
) -> Result<TempTwoBit> {
    ensure!(
        !contigs.is_empty(),
        "temporary two-bit reference must contain at least one repeating contig"
    );
    let sequences = contigs
        .iter()
        .map(RepeatingContigSpec::sequence)
        .collect::<Result<Vec<_>>>()?;
    twobit_from_sequences(name, sequences)
}

/// Create a temporary two-bit reference with one repeating contig.
///
/// The generated reference has exactly one contig. Use this for command tests
/// that only need a single `chr1`-style sequence and where the sequence pattern
/// should be visible at the call site.
///
/// The generated sequence is `pattern` repeated from its first base and
/// truncated to exactly `length` bases. For example, `pattern = "ACGT"` and
/// `length = 10` creates `ACGTACGTAC`. The returned `TempTwoBit` stores that
/// normalized sequence under `contig_name`.
pub fn twobit_with_single_repeating_contig(
    name: &str,
    contig_name: &str,
    pattern: &str,
    length: usize,
) -> Result<TempTwoBit> {
    twobit_with_repeating_contigs(
        name,
        &[RepeatingContigSpec::new(contig_name, pattern, length)],
    )
}

fn write_fasta(path: &Path, sequences: &[(String, String)]) -> Result<()> {
    let mut file = File::create(path)?;
    for (name, sequence) in sequences {
        writeln!(file, ">{name}")?;
        for chunk in sequence.as_bytes().chunks(60) {
            file.write_all(chunk)?;
            file.write_all(b"\n")?;
        }
    }
    Ok(())
}

fn normalize_sequence(sequence: &str) -> Result<String> {
    ensure!(
        !sequence.is_empty(),
        "two-bit contig sequences must not be empty"
    );
    let normalized = sequence.to_ascii_uppercase();
    for base in normalized.bytes() {
        match base {
            b'A' | b'C' | b'G' | b'T' | b'N' => {}
            _ => {
                return Err(anyhow!(
                    "two-bit sequence base must be A, C, G, T, or N, got {:?}",
                    char::from(base)
                ));
            }
        }
    }
    Ok(normalized)
}
