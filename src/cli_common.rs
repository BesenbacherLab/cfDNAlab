use crate::utils::bam::bam_header_contigs;
use anyhow::Context;
use anyhow::bail;
use std::{path::PathBuf, str::FromStr};

/// Args for in-/output and core (threads).
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug)]
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
    /// Defaults to the minimum of 22 (one thread per chromosome) and
    /// the number of available CPU cores (-1).
    #[cfg_attr(
        feature = "cli",
        clap(short = 't', long, default_value_t = (num_cpus::get()-1).max(1).min(22), help_heading = "Core")
    )]
    pub n_threads: usize,
}

/* Min/Max fragment lengths */

/// Args for setting minimum and maximum fragment length.
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug)]
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

/* Window selection */

// Windows option ENUM
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

// TODO: Consider allowing counting up the proportion?
#[derive(Debug, Clone, Copy, PartialEq, Default)]
/// How to assign a fragment to windows.
pub enum WindowAssigner {
    /// Assign to windows overlapping any of the fragment bases.
    #[default]
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
        if s == "all" {
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
            Err("Use 'any', 'all', 'midpoint', or 'proportion=<0.0–1.0>'".into())
        }
    }
}

#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone, Default)]
pub struct AssignToWindowArgs {
    /// The fragment positions that should overlap a window for it to be counted in that window `[string]`
    ///
    /// Possible values:
    ///     "any", "all", "midpoint", or "proportion=<threshold>"
    ///
    /// Example of proportion: `--assign-by proportion=0.2` (no space around `=`)
    ///
    /// **NOTE**: Ignored when no windows are specified.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "any",
            ignore_case = true,
            help = "What to assign fragments to windows by.",
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
    /// Optional scaling factors for normalizing/smoothing the genome `[path]`
    ///
    /// Coverage values are divided by the scaling factors of their bin.
    /// Scaling factors must thus be strictly positive.
    ///
    /// `.tsv` file with a header and the columns (names must match):
    ///     `chromosome, start, end, scaling_factor`.
    /// As produced by `cfdna normalize-genome`.
    ///
    /// Intervals are half-open `[start, end)` and must be non-overlapping per chromosome.
    /// We scale the coverage per position with `cov[i] / scaling_factor[bins[i]]` inside each bin.
    #[cfg_attr(
        feature = "cli",
        clap(
            short = 'i',
            long,
            value_parser,
            required = true,
            help_heading = "Normalization"
        )
    )]
    pub scaling_factors: Option<PathBuf>,
}
