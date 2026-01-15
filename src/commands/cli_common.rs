use crate::commands::fragment_kmers::positions::{BasesFrom, MismatchBasesFrom, ReferenceFrame};
use crate::shared::bam::bam_header_contigs;
use crate::shared::bam::{Contigs, bam_contigs_info};
use crate::shared::blacklist::load_blacklists;
use crate::shared::scale_genome::load_scaling_factors_tsv;
use anyhow::{Context, Result, bail, ensure};
use fxhash::FxHashMap;
use std::path::Path;
use std::{path::PathBuf, str::FromStr};

/// Minimum ACGT bases required when estimating GC fraction for sample reads.
pub const MIN_ACGT_BASES_FOR_GC_FRACTION: u32 = 10;

/// Args for in-/output and core (threads).
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone)]
pub struct IOCArgs {
    /// Indexed, coordinate-sorted BAM input file `[path]`
    ///
    /// Can be either **paired-end** or **single-end** (set `--single-end`).
    /// Single-end assumes the reads span their fragments exactly
    /// (so read size is fragment size).
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
pub struct SingleEndArgs {
    /// The input is single-end and the read spans exactly the full fragment `[flag]`
    ///
    /// Each aligned read is treated as a fragment spanning its aligned reference interval
    /// `[pos, reference_end)`. This uses the mapped span only
    /// (soft clips excluded).
    ///
    /// Cannot be combined with `--require-proper-pair` (when available).
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub single_end: bool,
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
pub struct Ref2BitOptionalForGCArgs {
    /// Optional 2bit reference genome file [path]
    ///
    /// NOTE: Required for GC correction, otherwise ignored.
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
    pub ref_2bit: Option<PathBuf>,
}

/* Min/Max fragment lengths */

/// Args for setting minimum and maximum fragment length.
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone)]
pub struct FragmentLengthArgs {
    /// Minimum fragment length to include `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "30", value_parser = clap::value_parser!(u32).range(MIN_ACGT_BASES_FOR_GC_FRACTION as i64..), help_heading="Filtering"))]
    pub min_fragment_length: u32,

    /// Maximum fragment length to include `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "1000", value_parser = clap::value_parser!(u32).range(MIN_ACGT_BASES_FOR_GC_FRACTION as i64..), help_heading="Filtering"))]
    pub max_fragment_length: u32,
}

impl FragmentLengthArgs {
    pub fn default() -> Self {
        Self {
            min_fragment_length: 30,
            max_fragment_length: 1000,
        }
    }
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
    ///     `"count-overlap"`, `"any"`, `"all"`, `"midpoint"`, or `"proportion=<threshold>"`
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

// TODO: Standardize whether lists should be comma-sep or space-sep

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
    /// For BAM-backed commands this uses the BAM header order.
    /// For commands that read chromosome order from their input,
    /// this may use the input order or some other order.
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
    /// `.tsv` file as produced by `cfdna coverage-weights` containing a scaling factor to *multipy* by per **scaling-bin**.
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
    ///   - end exactly at that chromosome’s length
    #[cfg_attr(
        feature = "cli",
        clap(long, value_parser, help_heading = "Normalization")
    )]
    pub scaling_factors: Option<PathBuf>,
}

#[cfg_attr(feature = "cli", derive(clap::Args))]
#[cfg_attr(
    feature = "cli",
    clap(
        group = clap::ArgGroup::new("gc_correction")
            .args(&["gc_file", "gc_tag"])
            .multiple(false)))]
#[derive(Debug, Clone, Default)]
pub struct ApplyGCArgs {
    /// Optional path to GC correction file *made from the same BAM file* with `gc-bias` `[path]`
    ///
    /// The file is usually called `gc_bias_correction.npz`.
    ///
    /// **NOTE**: Requires specifying the reference genome 2bit file as well.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser,
            group = "gc_correction",
            help_heading = "GC Correction (select max. one source)"
        )
    )]
    pub gc_file: Option<PathBuf>,

    /// Optional aux tag to get GC weight from when using external GC correction packages `[path]`
    ///
    /// Packages like `GCParagon` and `GCfix` allow saving GC weights directly to the reads
    /// in a BAM file. They often assign a "gc" aux tag.
    ///
    /// The average per-read weight is used to count the fragment. When any of the reads have a zero-weight,
    /// the fragment gets a zero-weight.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser,
            group = "gc_correction",
            help_heading = "GC Correction (select max. one source)"
        )
    )]
    pub gc_tag: Option<String>,

    /// Whether to drop fragments where the GC correction could no be calculated `[path]`
    ///
    /// If a GC correction weight could not be computed/retrieved for a fragment,
    /// the default is to weight it as `1.0` (no correction). If you prefer to
    /// exclude it instead, set this flag.
    #[cfg_attr(
        feature = "cli",
        clap(long, help_heading = "GC Correction (select max. one source)")
    )]
    pub drop_invalid_gc: bool,
}

#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone, Default)]
pub struct ApplyGCArgFileOnly {
    /// Optional path to GC correction file *made from the same BAM file* with `gc-bias` `[path]`
    ///
    /// The file is usually called `gc_bias_correction.npz`.
    ///
    /// **NOTE**: Requires specifying the reference genome 2bit file as well.
    #[cfg_attr(
        feature = "cli",
        clap(long, value_parser, help_heading = "GC Correction")
    )]
    pub gc_file: Option<PathBuf>,

    /// Whether to drop fragments where the GC correction could no be calculated `[path]`
    ///
    /// If a GC correction weight could not be computed for a fragment,
    /// the default is to weight it as `1.0` (no correction). If you prefer to
    /// exclude it instead, set this flag.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "GC Correction"))]
    pub drop_invalid_gc: bool,
}

// TODO: Is "nearest" clear enough in all usecases?
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone, Default)]
pub struct FragmentPositionSelectionArgs {
    /// Choose the reference frame that interprets every other region selection argument `[left|right|per-end|nearest|mid]`.
    ///
    /// When multiple frames are supplied (matching the same number of `positions` strings), the intersection of positions are used.
    /// The first frame+positions(+step) combination determines the output type.
    /// Using multiple specification allows selecting e.g., the `-N..N` bases around the fragment midpoint
    /// while limiting the distance to the fragment ends.
    ///
    /// Note: `--positions` describe positions to count at relative to the chosen frame.
    ///
    /// - **`left`** counts bases from the forward 5' end. Indices increase along the fragment and
    ///   k-mers are counted in the forward-orientation.
    ///
    /// - **`right`** counts bases from the reverse 5' end. Indices decrease along the fragment and
    ///   k-mers are counted in the reverse-orientation with **complemented** bases.
    ///
    /// - **`per-end`** counts both `left` and `right` simultaneously, producing two sets of k-mer counts.
    ///   The `step` start can differ per side.
    ///
    /// - **`nearest`** folds the fragment around the midpoint so distances grow away from the nearest end.
    ///   The positional keyword `half` represents the midpoint (and maximum position). For odd-sized fragments,
    ///   the single midpoint is not counted, as both sides count up-to it.
    ///   Bases contributed by the reverse 5' side are complemented.
    ///
    /// - **`mid`** centers the axis on the midpoint, allowing selections around zero with negative/positive offsets.
    ///   K-mers are counted in the forward-orientation.
    ///
    /// Pass multiple frames as e.g.: `--frame mid --frame left`.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            num_args = 1,
            action = clap::ArgAction::Append,
            default_values = ["left"],
            help_heading = "Region Selection"
        )
    )]
    pub frame: Vec<ReferenceFrame>,

    /// Describe which positions to count at relative to the selected frame `[string]`.
    ///
    /// Indices are **1-based inclusive**, why e.g. `1..10` would start at the first position and end at the tenth position (included).
    ///
    /// When multiple specifications are supplied (matching the same number of `frame`s), the intersection of positions are used.
    ///
    /// The allowed shapes depend on `--frame`:
    ///
    /// - **`left`**, **`right`**, **`per-end`**: use `..` for the full span or `A..B`, `A..`, `..B`, `A..-B`, `..half`, `A..half-K`.
    ///   For example, `1..10` keeps the first ten bases, `10..-10` trims both ends, and `..half-5`
    ///   includes bases from the start up to five before the fragment midpoint. Open intervals like `A..`
    ///   include every coordinate from `A` to the end of the frame.
    ///
    /// - **`nearest`** (folded 1..floor(length/2)): use `..` for every folded position or `A..B`, `A..`, `..B`, `..half`, `A..half-K`. Here, `half` expands to the
    ///   largest folded distance, ensuring the center base is maximally counted once. For odd-sized fragments, the central base remains uncounted.
    ///   Forms like `10..-10` are rejected for this frame.
    ///
    /// - **`mid`** (centered at 0): use `..` for the entire axis, `-M..N`, `-M..`, or `..N`. E.g. `-10..10` for the 20 bases around the midpoint.
    ///
    /// Pass multiple strings as e.g.: `--positions '-50..50' --positions '10..-10'`.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            help_heading = "Region Selection",
            num_args = 1,
            action = clap::ArgAction::Append,
            default_values = [".."],
            required = false,
            allow_hyphen_values = true
        )
    )]
    pub positions: Vec<String>,

    /// Downsample after selection by keeping every Nth index `[integer >= 1]`.
    ///
    /// When multiple `frame` and `positions` specifications are set, provide either a single step
    /// to use in all of them or a step per specification.
    ///
    /// Applied independently to each track in frame order (e.g., per-end left and right both stride through
    /// their own selections). Leave at 1 to keep every base.
    ///
    /// For the `mid` frame, zero is treated as the origin of the stride: when the chosen range includes the
    /// midpoint, it is always retained and every `step`th offset is kept symmetrically
    /// (`-2*step`, `-step`, `0`, `step`, `2*step`, ...). Ranges that exclude the origin fall back to the default stride.
    ///
    /// Pass multiple steps as e.g.: `--step 1 --step 2`.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_values_t = [1usize], num_args = 1, help_heading = "Region Selection")
    )]
    pub step: Vec<usize>,
}

pub struct UnparsedPositionalSelectionSpec {
    pub frame: ReferenceFrame,
    pub positions: String,
    pub step: usize,
}

impl std::fmt::Display for UnparsedPositionalSelectionSpec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "frame={}, positions=\"{}\", step={}",
            self.frame.as_str(),
            self.positions,
            self.step
        )
    }
}

impl FragmentPositionSelectionArgs {
    /// Check the selection args and convert to vec of specifications
    ///
    /// Each specification has one `frame`, `positions`, and `step`.
    pub fn into_positional_specs(self) -> Result<Vec<UnparsedPositionalSelectionSpec>> {
        // Destructure to get the fields as variables
        let FragmentPositionSelectionArgs {
            frame,
            positions,
            step,
        } = self;

        // Number of specifications
        let n = frame.len();

        ensure!(
            n == positions.len(),
            "--frame and --positions must have the same number of values (got {} vs {})",
            n,
            positions.len()
        );

        // Enforce each frame appears at most once without requiring Hash/Ord
        for i in 0..n {
            if frame[..i].contains(&frame[i]) {
                bail!("--frame contains a duplicate value: {:?}", frame[i]);
            }
        }

        // Resolve step: either one value reused or exactly n values
        let resolved_step = match step.len() {
            1 => vec![step[0]; n],
            len if len == n => step,
            other => bail!(
                "--step must be provided once or exactly {} times (got {})",
                n,
                other
            ),
        };

        // Basic sanity: steps must be >= 1
        if let Some(&bad) = resolved_step.iter().find(|&&s| s == 0) {
            bail!("--step must be >= 1 (found {})", bad);
        }

        Ok(frame
            .into_iter()
            .zip(positions.into_iter())
            .zip(resolved_step.into_iter())
            .map(
                |((frame, positions), step)| UnparsedPositionalSelectionSpec {
                    frame,
                    positions,
                    step,
                },
            )
            .collect())
    }
}

#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone, Default)]
pub struct BaseSelectionArgs {
    /// Choose which coordinate source defines the counted positions `[reference|prefer-reads|reads|nearest-read]`
    ///
    /// - `reference`: Always use the reference span, even when reads do not cover those bases.
    ///
    /// - `prefer-reads`: Use read-space coordinates whenever an observed base covers the requested position
    ///   and fall back to the reference span when reads don't cover the positions.
    ///
    /// - `reads`: Only count positions the reads cover.
    ///
    /// - `nearest-read`: Clamp the selection to the read that corresponds to the frame origin (e.g., the
    ///   left/forward read for the `left` frame).
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "reference",
            help_heading = "Region Selection"
        )
    )]
    pub bases_from: BasesFrom,

    /// Resolve overlapping read **mismatches** when preferring read bases `[nearest-read|base-quality|reference]`
    ///
    /// - `nearest-read`: Take the base from whichever read is closest to the frame origin. **NOTE**: Incompatible with `--frame mid`.
    ///
    /// - `base-quality`: Take the base with the highest quality score.
    ///
    /// - `reference`: Ignore the reads and fall back to the reference base for that coordinate.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "nearest-read",
            help_heading = "Region Selection"
        )
    )]
    pub mismatch_bases_from: MismatchBasesFrom,
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
    bam_path: &Path,
) -> Result<(Vec<String>, Contigs)> {
    let chromosomes = chrom_args
        .resolve_chromosomes(Some(bam_path))
        .context("resolve chromosomes")?;
    let contigs = bam_contigs_info(&bam_path, &chromosomes).context("fetch contig metadata")?;
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
