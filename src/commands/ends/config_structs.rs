use std::str::FromStr;

// TODO these structs should use the format used by other cfDNAlab commands instead

#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
#[derive(Copy, Clone, Debug)]
pub enum KmerSource {
    /// Extract k-mer from the sequenced read
    Read,
    /// Extract k-mer from the reference genome
    Reference,
}

#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Copy, Clone)]
pub struct ClippingArgs {
    /// Soft-clip handling at fragment ends
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "raw",
            requires_if("FillWithRef", "ref_2bit")
        )
    )]
    pub clip_strategy: ClipStrategy,

    /// Skip a fragment if single-end (S+H) clipping exceeds this many bases (default= --kmer-size - 1)
    ///
    /// 0 = no limit, use `--kmer_source drop-clipped` to discard all clipped reads
    ///
    /// When total clipping reduces the remaining fragment's size to `< --min-length` we skip it as well
    #[clap(long)]
    pub max_end_clips: Option<usize>,

    /// (set internally after parse; not a CLI arg)
    #[clap(skip)]
    pub max_end_clips_resolved: usize,
}

#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
#[derive(Copy, Clone, Debug)]
pub enum ClipStrategy {
    /// Take the raw read bases, clipped or not
    Raw,

    /// Drop the fragment if that end is soft-clipped
    DropClipped,

    /// Slide into the aligned region: skip any soft-clips, take the first k aligned bases
    AlignStart,

    /// For reads only: when an end is soft-clipped, pull the missing bases from the reference
    FillWithRef,
}

#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone, Default)]
pub struct AssignMotifToWindowArgs {
    /// The **fragment positions** that should overlap a window for it to be counted in that window,
    /// OR the option to count the fraction of overlapping bases `[string]`
    ///
    /// Possible values:
    ///     `"endpoint"`, `"count-overlap"`, `"any"`, `"all"`, `"midpoint"`, or `"proportion=<threshold>"`
    ///
    /// `"endpoint"`: By default, the motif is counting in the windows overlapping the fragment end position.
    ///
    /// `"count-overlap"`: Count up the fraction of overlapping fragment bases.
    ///
    /// Example of proportion: `--assign-by proportion=0.2` (no space around `=`)
    ///
    /// Midpoints for even-sized fragments are randomly selected as either the left or right base
    /// to avoid bias.
    ///
    /// **NOTE**: In the rare case where windows are smaller than fragments, it's still
    /// the proportion of the fragment positions that overlap that is considered. If the window
    /// size is 30% of the fragment size, that fragment cannot overlap more than 30%.
    ///
    /// **NOTE**: Ignored when no windows are specified.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "endpoint",
            ignore_case = true,
            help = "What to assign fragments to windows by (or count fragments as).",
            help_heading = "Window Assignment"
        )
    )]
    pub assign_by: WindowMotifAssigner,
}

// TODO: In the future we might want to add window-based overlap variants (WindowProportion etc.). Not relevant yet.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
/// How to assign a fragment to windows.
///
/// NOTE: This only considers the proportion of **fragment positions**
/// overlapping the window. For window sizes smaller than fragments
/// this means a fragment could overlap a window fully but
/// have < 100% of fragment positions inside the window.
pub enum WindowMotifAssigner {
    /// Assign to windows overlapping the fragment midpoint.
    #[default]
    Endpoint,
    /// Count up the fraction of overlapping fragment bases.
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

impl FromStr for WindowMotifAssigner {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s == "count-overlap" {
            Ok(WindowMotifAssigner::CountOverlap)
        } else if s == "all" {
            Ok(WindowMotifAssigner::All)
        } else if s == "any" {
            Ok(WindowMotifAssigner::Any)
        } else if s == "endpoint" {
            Ok(WindowMotifAssigner::Endpoint)
        } else if s == "midpoint" {
            Ok(WindowMotifAssigner::Midpoint)
        } else if let Some(v) = s.strip_prefix("proportion=") {
            let thr: f64 = v
                .parse()
                .map_err(|e: std::num::ParseFloatError| e.to_string())?;
            if !(0.0..=1.0).contains(&thr) {
                Err("Proportion must be between 0.0 and 1.0".into())
            } else {
                Ok(WindowMotifAssigner::Proportion(thr))
            }
        } else {
            Err("Use 'endpoint', 'count-overlap', 'any', 'all', 'midpoint', or 'proportion=<0.0–1.0>'".into())
        }
    }
}
