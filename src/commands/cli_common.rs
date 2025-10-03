use crate::shared::bam::bam_header_contigs;
use crate::shared::bam::{Contigs, bam_contigs_info};
use crate::shared::blacklist::load_blacklists;
use crate::shared::scale_genome::load_scaling_factors_tsv;
use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;
use std::path::Path;
use std::{path::PathBuf, str::FromStr};

/// Args for in-/output and core (threads).
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone)]
pub struct IOCArgs {
    /// Indexed, coordinate-sorted BAM input file `[path]`
    #[cfg_attr(
        feature = "cli",
        clap(
            short = 'i',
            long,
            value_parser,
            required = true,
            help_heading = "Core"
        )
    )]
    pub bam: PathBuf,

    /// Output directory for results `[path]`
    #[cfg_attr(
        feature = "cli",
        clap(
            short = 'o',
            long,
            value_parser,
            required = true,
            help_heading = "Core"
        )
    )]
    pub output_dir: PathBuf,

    /// Number of threads to use (increases RAM usage) `[integer]`
    ///
    /// Defaults to the number of available CPU cores (-1).
    #[cfg_attr(
        feature = "cli",
        clap(short = 't', long, default_value_t = (num_cpus::get()-1).max(1), help_heading = "Core")
    )]
    pub n_threads: usize,
}

#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone)]
pub struct Ref2BitRequiredArgs {
    /// 2bit reference genome file [path]
    ///
    /// E.g., "hg38.2bit" from UCSC ( https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit ).
    #[cfg_attr(
        feature = "cli",
        clap(
            short = 'r',
            long,
            value_parser,
            required = true,
            help_heading = "Core"
        )
    )]
    pub ref_2bit: PathBuf,
}

#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone)]
pub struct Ref2BitOptionalArgs {
    /// 2bit reference genome file [path]
    ///
    /// E.g., "hg38.2bit" from UCSC ( https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit ).
    #[cfg_attr(
        feature = "cli",
        clap(
            short = 'r',
            long,
            value_parser,
            required = false,
            help_heading = "Core"
        )
    )]
    pub ref_2bit: PathBuf,
}

/* Min/Max fragment lengths */

/// Args for setting minimum and maximum fragment length.
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone)]
pub struct FragmentLengthArgs {
    /// Minimum fragment length to include `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "20", value_parser = clap::value_parser!(u32).range(1..), help_heading="Filtering"))]
    pub min_fragment_length: u32,

    /// Maximum fragment length to include `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "1000", value_parser = clap::value_parser!(u32).range(1..), help_heading="Filtering"))]
    pub max_fragment_length: u32,
}

impl FragmentLengthArgs {
    /// Check whether a fragment length is within the configured inclusive range.
    pub fn contains(&self, len: u32) -> bool {
        len >= self.min_fragment_length && len <= self.max_fragment_length
    }
}

/* Window selection */

/// The windowing options `[ENUM]`
///
/// Whether to perform a command globally (1 overall genomic window)
/// or in windows specified with a BED file or a fixed window size.
#[derive(Debug, Clone)]
pub enum WindowSpec {
    Global,
    Size(u64),
    Bed(PathBuf),
}

#[cfg_attr(feature = "cli", derive(clap::Args))]
#[cfg_attr(
    feature = "cli",
    clap(
        // At most one of the two flags; if none -> Global in `resolve()`
        group = clap::ArgGroup::new("windows")
            .args(&["by_size", "by_bed"])
            .multiple(false)
    )
)]
#[derive(Debug, Clone, Default)]
pub struct WindowsArgs {
    /// Window definition: a fixed window size `[integer]`
    ///
    /// Default is one global window.
    #[cfg_attr(
        feature = "cli",
        clap(
            long = "by-size",
            alias = "by",
            value_parser,
            group = "windows",
            help_heading = "Windows (select max. one arg.)"
        )
    )]
    pub by_size: Option<u64>,

    /// Window definition: a BED file of windows `[path]`
    #[cfg_attr(
        feature = "cli",
        clap(
            long = "by-bed",
            value_parser,
            group = "windows",
            help_heading = "Windows (select max. one arg.)"
        )
    )]
    pub by_bed: Option<PathBuf>,
}

impl WindowsArgs {
    /// If neither flag is set, default to `Global`.
    pub fn resolve_windows(&self) -> WindowSpec {
        if let Some(n) = self.by_size {
            WindowSpec::Size(n)
        } else if let Some(p) = self.by_bed.clone() {
            WindowSpec::Bed(p)
        } else {
            WindowSpec::Global
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
/// How to assign a fragment to windows.
pub enum WindowAssigner {
    /// Count up the fraction of overlapping fragment bases.
    #[default]
    CountOverlap,
    /// Assign to windows overlapping any of the fragment bases.
    Any,
    /// Assign to windows overlapping all of the fragment bases.
    All,
    /// Assign to windows overlapping the fragment midpoint.
    Midpoint,
    /// Assign to windows overlapping a given percentage of the fragment bases.
    Proportion(f64),
}

impl FromStr for WindowAssigner {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s == "count-overlap" {
            Ok(WindowAssigner::CountOverlap)
        } else if s == "all" {
            Ok(WindowAssigner::All)
        } else if s == "any" {
            Ok(WindowAssigner::Any)
        } else if s == "midpoint" {
            Ok(WindowAssigner::Midpoint)
        } else if let Some(v) = s.strip_prefix("proportion=") {
            let thr: f64 = v
                .parse()
                .map_err(|e: std::num::ParseFloatError| e.to_string())?;
            if !(0.0..=1.0).contains(&thr) {
                Err("Proportion must be between 0.0 and 1.0".into())
            } else {
                Ok(WindowAssigner::Proportion(thr))
            }
        } else {
            Err("Use 'count-overlap', 'any', 'all', 'midpoint', or 'proportion=<0.0–1.0>'".into())
        }
    }
}

#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone, Default)]
pub struct AssignToWindowArgs {
    /// The fragment positions that should overlap a window for it to be counted in that window,
    /// OR the option to count the fraction of overlapping bases `[string]`
    ///
    /// Possible values:
    ///     "count-overlap", "any", "all", "midpoint", or "proportion=<threshold>"
    ///
    /// `'count-overlap'`: Count up the fraction of overlapping fragment bases.
    ///
    /// Example of proportion: `--assign-by proportion=0.2` (no space around `=`)
    ///
    /// Midpoints for even-sized fragments are randomly selected as either the left or right base
    /// to avoid bias.
    ///
    /// **NOTE**: Ignored when no windows are specified.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "count-overlap",
            ignore_case = true,
            help = "What to assign fragments to windows by (or count fragments as).",
            help_heading = "Window Assignment"
        )
    )]
    pub assign_by: WindowAssigner,
}

/* Chromosome selection */

#[cfg_attr(feature = "cli", derive(clap::Args))]
#[cfg_attr(
    feature = "cli",
    clap(
        group = clap::ArgGroup::new("chrom_select")
            .args(&["chromosomes", "chromosomes_file"])
            .multiple(false)))]
#[derive(Debug, Clone, Default)]
pub struct ChromosomeArgs {
    /// Names of chromosomes to process (comma-separated or repeated). E.g. `'chr1,chr2,chr3'`.
    ///
    /// When no chromosomes are specified, it defaults to `chr1..chr22`.
    ///
    /// Specify `"all"` *as the only string* to use all present chromosomes.
    /// Only works for tools where a BAM path is passed.
    #[cfg_attr(
        feature = "cli", clap(
            long, num_args = 1..,
            value_parser,
            value_delimiter = ',',
            group = "chrom_select", 
            help_heading="Chromosome Selection (select max. one arg.)"))]
    pub chromosomes: Option<Vec<String>>,

    /// File with chromosome names to process (one per line).
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser,
            group = "chrom_select",
            help_heading = "Chromosome Selection (select max. one arg.)"
        )
    )]
    pub chromosomes_file: Option<PathBuf>,
}

impl ChromosomeArgs {
    /// Returns the final chromosome list, in priority order:
    /// 1) from `--chromosomes-file`
    /// 2) from `--chromosomes`
    /// 3) default `chr1`..`chr22`
    pub fn resolve_chromosomes(
        &self,
        bam_path: Option<&std::path::Path>,
    ) -> anyhow::Result<Vec<String>> {
        if let Some(file) = &self.chromosomes_file {
            let text: String = std::fs::read_to_string(file)
                .context(format!("reading chromosome file {:?}", file))?;
            let list: Vec<String> = text
                .lines()
                .map(str::trim)
                .filter(|l| !l.is_empty() && !l.starts_with('#'))
                .map(String::from)
                .collect();
            Ok(list)
        } else if let Some(chrs) = &self.chromosomes {
            if chrs.len() == 1 && chrs[0].eq_ignore_ascii_case("all") {
                let Some(bam) = bam_path else {
                    bail!(
                        "`--chromosomes all` requires `--bam <file>` to read contigs from the BAM header"
                    );
                };
                return bam_header_contigs(bam);
            }
            Ok(chrs.clone())
        } else {
            Ok((1..=22).map(|i| format!("chr{}", i)).collect())
        }
    }
}

/* Genomic scaling (applying normalize_genome) */

#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone, Default)]
pub struct ScaleGenomeArgs {
    /// Optional path to *non-negative* scaling factors for normalizing/smoothing the genome `[path]`
    ///
    /// `.tsv` file as produced by `cfdna normalize-genome` containing a scaling factor to *multipy* by per **scaling-bin**.
    ///
    /// The scaling-bin-overlapping parts of the fragments are counted as the scaling factor of the bin (`w=sf`).
    ///
    /// File Requirements
    ///
    /// -----------------
    ///
    /// The TSV file **must** have a header. Column names are matched **case-insensitively**.
    ///
    /// Required columns: `chromosome`, `start`, `end`, `scaling_factor`.
    ///
    /// Coordinates are 0-based, half-open `[start, end)`.
    ///
    /// `scaling_factor` must be finite and strictly >= 0.
    ///
    /// Bins are filtered to the provided `chromosomes`.
    ///
    /// For every chromosome in `chromosomes`, bins must:
    ///
    ///   - start at 0
    ///
    ///   - be perfectly contiguous (no gaps, no overlaps)
    ///
    ///   - end exactly at that chromosome’s length (from `contigs`)
    #[cfg_attr(
        feature = "cli",
        clap(long, value_parser, help_heading = "Normalization")
    )]
    pub scaling_factors: Option<PathBuf>,
}

// Common loaders

/// Resolve chromosomes and BAM contig metadata once for a command.
///
/// Implementation details:
/// - Delegates to `ChromosomeArgs::resolve_chromosomes`, passing the BAM path so
///   aliases such as `--chromosomes all` work uniformly.
/// - Queries the BAM header via `bam_contigs_info` to obtain target lengths.
///
/// Parameters:
/// - `chrom_args`: Command-line chromosome selection configuration.
/// - `ioc`: Shared IO arguments providing the BAM path.
///
/// Returns:
/// - A tuple with the resolved chromosome names and their contig metadata.
///
/// Errors:
/// - Propagates IO and parsing failures when the BAM cannot be opened or lacks
///   the requested contigs.
pub fn resolve_chromosomes_and_contigs(
    chrom_args: &ChromosomeArgs,
    ioc: &IOCArgs,
) -> Result<(Vec<String>, Contigs)> {
    let chromosomes = chrom_args
        .resolve_chromosomes(Some(ioc.bam.as_path()))
        .context("resolve chromosomes")?;
    let contigs = bam_contigs_info(&ioc.bam, &chromosomes).context("fetch contig metadata")?;
    Ok((chromosomes, contigs))
}

/// Create the output directory if it does not exist.
///
/// Implementation details:
/// - Wraps `std::fs::create_dir_all` with an `anyhow` context to yield helpful
///   error messages tailored to the target path.
///
/// Parameters:
/// - `path`: Directory where the command should place its results.
///
/// Returns:
/// - `Ok(())` if the directory exists or was created successfully.
///
/// Errors:
/// - Returns an error when the directory cannot be created (for instance due to
///   missing permissions or an unwritable parent directory).
pub fn ensure_output_dir(path: &Path) -> Result<()> {
    std::fs::create_dir_all(path)
        .with_context(|| format!("cannot create output directory: {}", path.display()))
}

/// Load blacklist intervals when the user supplied one or more BED files.
///
/// Implementation details:
/// - Delegates to `load_blacklists`, which merges overlapping intervals per
///   chromosome and enforces `min_size` filtering.
/// - Returns an empty map when `beds` is `None` so callers can operate without
///   additional branching.
///
/// Parameters:
/// - `beds`: Optional list of BED paths.
/// - `min_size`: Minimum interval size (bp) to retain.
/// - `chromosomes`: Chromosomes the command intends to process.
///
/// Returns:
/// - A map keyed by chromosome name containing sorted blacklist intervals.
///
/// Errors:
/// - Propagates parsing errors if any BED file is malformed or unavailable.
pub fn load_blacklist_map(
    beds: Option<&Vec<std::path::PathBuf>>,
    min_size: u64,
    halo_bp: u64,
    chromosomes: &Vec<String>,
) -> Result<FxHashMap<String, Vec<(u64, u64)>>> {
    if let Some(paths) = beds {
        load_blacklists(paths, min_size, halo_bp, Some(chromosomes.as_slice()))
    } else {
        Ok(FxHashMap::default())
    }
}

/// Load per-chromosome scaling factors (if provided).
///
/// Implementation details:
/// - Uses `load_scaling_factors_tsv` to parse the command-line TSV into a
///   chromosome keyed map of `(start, end, factor)` tuples.
/// - Returns an empty map when no scaling factors were supplied, avoiding
///   unnecessary allocations inside the calling code.
///
/// Parameters:
/// - `scale_args`: Normalisation argument bundle.
/// - `chromosomes`: Chromosome ordering requested by the command.
/// - `contigs`: BAM target metadata, used to validate the TSV content.
///
/// Returns:
/// - A scaling factor map ready for lookups by chromosome.
///
/// Errors:
/// - Propagates IO or format errors when the TSV cannot be read or does not
///   match the BAM contigs.
pub fn load_scaling_map(
    scale_args: &ScaleGenomeArgs,
    chromosomes: &[String],
    contigs: &Contigs,
) -> Result<FxHashMap<String, Vec<(u64, u64, f32)>>> {
    if let Some(path) = &scale_args.scaling_factors {
        load_scaling_factors_tsv(path, chromosomes, contigs).context("load scaling factors")
    } else {
        Ok(FxHashMap::with_hasher(Default::default()))
    }
}
