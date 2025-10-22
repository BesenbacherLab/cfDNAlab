use crate::shared::blacklist::BlacklistStrategy;
#[cfg(feature = "cli")]
use clap::{ArgGroup, Parser, ValueEnum};
use std::path::PathBuf;

/// Clean and standardise genomic windows so downstream cfDNA tools receive a tidy BED.
///
/// `prep-windows` reads a delimited table with at least `chrom,start,end`, validates every row,
/// and writes a canonical BED-like file that downstream tools can reuse. The command keeps
/// your metadata columns during processing but emits only well-behaved coordinates plus an optional
/// `group` label for downstream tools.
///
/// A *group* is simply a tag that tells later analyses how to partition the windows (for example,
/// promoter vs enhancer sets, distance quartiles). You can supply it from existing
/// columns or instruct the command to derive it while reshaping the windows.
///
/// The command parses the TSV/CSV input and:
///
/// - Filters windows using score thresholds, blacklist overlap, deduplication, and distance to nearest same-group window.
///
/// - Adjust coordinates by resizing to a specific size or adding flanks to the current sizes (trimmed to chromosome limits).
///
/// - Build or refine groups by combining input columns or subdividing windows by their distance
///   to elements in the `near`-file. Windows can be merged when close to other windows in/across groups.
///
/// The output is minimal, headerless, sorted by `(chrom, start, end, group)`, and ready for commands
/// such as `profile-groups`.
///
/// ## Practical notes
///
/// - All coordinates are 0-based half-open `[start, end)`.
///
/// - Column indices are 0-based when you refer to them explicitly.
///
/// - Blacklist checks run on the final window span using the halo you configure.
///
/// - "Nearest distance" refers to the closest edge of the comparison interval. NOTE:
///   For point features (e.g., TSS), provide 1-bp intervals at the strand-specific coordinate.
///
/// - Output is sorted by `(chrom, start, end, group)`.
#[cfg_attr(feature = "cli", derive(Parser, Clone))]
#[cfg_attr(
    feature = "cli",
    command(
        group(
            ArgGroup::new("resize_group")
                .args(&["resize", "flank"])
                .multiple(false)
        )
    )
)]
pub struct PrepareConfig {
    // ─────────────────────────────────────────────────────────────────────────────
    // Core (I/O and schema)
    // ─────────────────────────────────────────────────────────────────────────────
    /// Input BED-like file `[path]`
    ///
    /// Compression inferred from file extension (.gz or .zst). Use '-' for stdin.
    #[cfg_attr(
        feature = "cli",
        clap(long, value_parser, required = true, help_heading = "Core")
    )]
    pub input: PathBuf,

    /// Output BED-like file `[path]`
    ///
    /// Compression inferred from file extension (.gz or .zst). Use '-' for stdout.
    #[cfg_attr(
        feature = "cli",
        clap(long, value_parser, default_value = "-", help_heading = "Core")
    )]
    pub output: PathBuf,

    /// Header presence in input `[string]`
    ///
    /// - "auto": Infer from first line.
    ///
    /// - "present": First line is a header line with column names.
    ///
    /// - "absent": No header; only indices allowed in `--cols` and related.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "auto",
            ignore_case = true,
            help_heading = "Core"
        )
    )]
    pub header: HeaderMode,

    /// Field separator for input and output `[char]`
    ///
    /// Common separators are `\t` (accepts "tab") for `.tsv` files and `,` or `;` for `.csv` files.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser = parse_sep,
            required = true,
            help_heading = "Core"
        )
    )]
    pub sep: char,

    /// Column mapping for the *required* first three columns `[string]`
    ///
    /// Format: `chrom=<idx>,start=<idx>,end=<idx>` (0-based indices only).
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser,
            default_value = "chrom=0,start=1,end=2",
            help_heading = "Core"
        )
    )]
    pub cols: String,

    // TODO: Change the other cols to accept u32
    /// Optional group columns `[strings]`
    ///
    /// Provide one or more column indices. When multiple are given, they will be
    /// concatenated using `__` into the single `group` output column.
    ///
    /// If omitted, the `group` will be derived from later subdivision steps (e.g. `--distance-bins`).
    ///
    /// If no subdivision occurs, no group column is written to the output.
    ///
    /// Example: `--group-cols 3`
    #[cfg_attr(feature = "cli", clap(long, value_parser, help_heading = "Core"))]
    pub group_cols: Vec<String>,

    /// Optional score column `[string]`
    ///
    /// Column index for a numeric score used by `--score-filter`.
    #[cfg_attr(
        feature = "cli",
        clap(long, value_parser, requires = "score-filter", help_heading = "Core")
    )]
    pub score_col: Option<String>,

    // ─────────────────────────────────────────────────────────────────────────────
    // Score filtering
    // ─────────────────────────────────────────────────────────────────────────────
    /// Score filter expression `[string]`
    ///
    /// Applies only if `--score-col` is set. Supported operators: `>=, >, <=, <, ==, !=`.
    /// Examples: `--score-filter ">=10"`, `--score-filter "<0.05"`.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser,
            requires = "score-col",
            help_heading = "Score filtering"
        )
    )]
    pub score_filter: Option<String>,

    /// Behavior for missing scores `[string]`
    ///
    /// - "keep": Keep records with missing/invalid scores.
    ///
    /// - "drop": Drop records with missing/invalid scores.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "keep",
            ignore_case = true,
            help_heading = "Score filtering"
        )
    )]
    pub score_missing: MissingScore,

    // ─────────────────────────────────────────────────────────────────────────────
    // Blacklist
    // ─────────────────────────────────────────────────────────────────────────────
    /// Optional BED file(s) with blacklisted regions `[path]`
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            short = 'b',
            value_parser,
            num_args = 1..,
            action = clap::ArgAction::Append,
            help_heading = "Blacklist"
        )
    )]
    pub blacklist: Option<Vec<PathBuf>>,

    /// Halo (bp) to expand blacklist intervals on both sides `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "0",
            value_parser = clap::value_parser!(u32),
            help_heading = "Blacklist"
        )
    )]
    pub blacklist_halo: u32,

    /// Strategy for determining when a window is blacklisted `[string]`
    ///
    /// Possible values:
    ///
    /// - `"any"`: Any overlap > 0 with a blacklist interval (after halo).
    ///
    /// - `"all"`: Window is fully contained within a blacklist interval.
    ///
    /// - `"midpoint"`: Window midpoint lies inside a blacklist interval.
    ///
    /// - `"proportion=<threshold>"`: Overlap proportion with respect to the window is ≥ threshold.
    ///
    /// Example: `--blacklist-strategy proportion=0.2`
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            alias = "bl-strategy",
            default_value = "any",
            ignore_case = true,
            help_heading = "Blacklist"
        )
    )]
    pub blacklist_strategy: BlacklistStrategy,

    // ─────────────────────────────────────────────────────────────────────────────
    // Distance to ‘near' intervals
    // ─────────────────────────────────────────────────────────────────────────────
    /// BED-like file with target intervals to compute nearest distance to `[path]`
    ///
    /// Expected columns:
    ///   `chromosome`, `start`, `end`, optional `strand`, optional `group`.
    ///
    /// - `strand` is one of `+`, `-`, or `.` (unknown). If absent, `+` is assumed.
    ///
    /// - Intervals must be half-open, non-overlapping, and have unique edges per chromosome.
    ///
    /// Distance:
    ///
    ///   - Overlap: Distance is `0`.
    ///
    ///   - Otherwise: Distance is the minimum from the window’s edges to the selected
    ///     target edge(s) (see `--near-edge`).
    ///
    /// Strand semantics:
    ///
    ///   - "Upstream/Downstream" are defined **relative to the near interval’s annotated strand**.
    ///
    ///   - If `strand` is unknown (`.`), upstream/downstream edge selection falls back to genomic-nearest.
    ///
    /// Header handling follows `--near-header`.
    ///
    /// If you require e.g., TSS distances, provide 1-bp intervals at the strand-specific TSS.
    ///
    /// When `--distance-bins` is used, the output group combines the original group (from
    /// `--group-cols`, if any), the nearest record’s group (if present in the near file),
    /// and the bin label:
    ///
    /// - With both: `{input_group}.{near_group}.{bin_label}`
    ///
    /// - Only input group: `{input_group}.{bin_label}`
    ///
    /// - Only near group: `{near_group}.{bin_label}`
    ///
    /// - Neither: `{bin_label}`
    #[cfg_attr(
        feature = "cli",
        clap(long, value_parser, help_heading = "Distance to near intervals")
    )]
    pub near: Option<PathBuf>,

    /// Header presence in the `--near` file `[string]`
    ///
    /// - "auto": Infer from first line.
    ///
    /// - "present": First line is a header line with column names.
    ///
    /// - "absent": No header; only indices allowed when referencing columns from the near file.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "auto",
            ignore_case = true,
            help_heading = "Distance to near intervals"
        )
    )]
    pub near_header: HeaderMode,

    /// Edge of near-intervals to consider in distance calculation `[string]`
    ///
    /// - "left": Use left genomic edge only.
    ///
    /// - "right": Use right genomic edge only.
    ///
    /// - "nearest": Use whichever genomic edge is closer (default).
    ///
    /// - "upstream": Use the edge that is upstream of each near interval given its strand (`+` uses left edge, `-` uses right edge).
    ///
    /// - "downstream": Use the edge that is downstream of each near interval given its strand (`+` uses right edge, `-` uses left edge).
    ///
    /// If a near interval’s strand is unknown (`.`), "upstream"/"downstream" behave like "nearest".
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "nearest",
            ignore_case = true,
            help_heading = "Distance to near intervals"
        )
    )]
    pub near_edge: NearEdge,

    /// Directionality of distance classification `[string]`
    ///
    /// - "upstream": Consider only near intervals that lie upstream (or overlap) relative to each near interval’s strand.
    ///
    /// - "downstream": Consider only near intervals that lie downstream (or overlap) relative to each near interval’s strand.
    ///
    /// - "both": Consider upstream and downstream (default).
    ///
    /// Overlaps are always allowed, returning zero distance.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "both",
            ignore_case = true,
            help_heading = "Distance to near intervals"
        )
    )]
    pub near_direction: NearDirection,

    /// How to respond when multiple near intervals tie for the minimum distance `[string]`
    ///
    /// - "annotate": keep the window and include both sides in the near label (e.g. `-A/+B`).
    ///
    /// - "drop": discard the window when a tie occurs.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "annotate",
            ignore_case = true,
            help_heading = "Distance to near intervals"
        )
    )]
    pub near_ties: NearTiePolicy,

    /// How to treat the computed distances when binning `[string]`
    ///
    /// - "absolute": Use `abs(distance)` for comparisons and thresholds.
    ///
    /// - "signed": Use signed distances.
    ///
    /// **Distance sign (when `--distance-sign signed`):**
    ///
    /// - Upstream of the near interval -> **negative** distance.
    ///
    /// - Downstream of the near interval -> **positive** distance.
    ///
    /// - Overlap/touch -> `0` distance.
    ///
    /// **Upstream/Downstream definition (strand-aware):**
    ///
    /// - For a `+` near interval: upstream is to the left (smaller genomic coordinates); downstream is to the right.
    ///
    /// - For a `-` near interval: upstream is to the right (larger genomic coordinates); downstream is to the left.
    ///
    /// - For an unknown strand (`.`): upstream/downstream are derived from genomic placement to the chosen target edge(s)
    ///   (falls back to genomic-nearest semantics).
    ///
    /// **Group label prefix (always emitted):**
    ///
    /// `-` = upstream, `+` = downstream, `=` = overlap.
    ///   
    /// Prefixes are included even when using `--distance-sign absolute`, so the side remains visible.
    ///
    /// **Examples:**
    ///
    /// Legend: '=' near interval, '#' window, '-' empty span. Signs are relative to the near interval’s strand.
    ///
    /// Case A: near is `+` strand
    ///
    ///     > coordinates:  100   120  140         200   220   240
    ///     >               |#####|----|===========|-----|#####|
    ///     >   upstream (-) ^^^^^         near           ^^^^^ downstream (+)
    ///     
    /// Case B: near is `-` strand
    ///
    ///     > coordinates:  100   120  140         200   220   240
    ///     >               |#####|----|===========|-----|#####|
    ///     > downstream (+) ^^^^^         near           ^^^^^ upstream (-)
    ///     
    /// Case C: overlap
    ///     
    ///     > coordinates:  100   120             200
    ///     >               |--###|===========###====|
    ///     >        touch (=) ^^^    near    ^^^ overlap (=)
    ///
    /// **Ties and overlaps:**
    ///
    /// - Overlap yields distance 0 and label prefix `=`.
    ///
    /// - If upstream/downstream tie and `--near-ties annotate`, both sides are reported,
    ///   e.g. `-GeneA/+GeneB`.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "absolute",
            ignore_case = true,
            help_heading = "Distance to near intervals"
        )
    )]
    pub distance_sign: DistSign,

    /// Distance bin specifications `[quoted strings]`
    ///
    /// Provide one or more `'<label>:<expr>'` rules. The first matching rule wins.
    /// Expression forms: `<N`, `<=N`, `A-B`, `>=N`, `>N` (N in bp).
    ///
    /// Examples:
    ///
    /// - `--distance-bins 'prox:<500' 'mid:500-2000' 'dist:>2000'`
    ///
    /// - `--distance-bins 'upstream:<0' 'at:0-0' 'downstream:>0'` (when using `--distance-sign signed`)
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser,
            num_args = 1..,
            help_heading = "Distance to near intervals"
        )
    )]
    pub distance_bins: Option<Vec<String>>,

    /// Maximum absolute distance to keep `[integer]`
    ///
    /// Windows with `|distance| > distance-max` are dropped prior to binning.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser = clap::value_parser!(u32),
            help_heading = "Distance to near intervals"
        )
    )]
    pub distance_max: Option<u32>,

    // ─────────────────────────────────────────────────────────────────────────────
    // Resizing / flanking (mutually exclusive)
    // ─────────────────────────────────────────────────────────────────────────────
    /// Resize window to a fixed size centered on midpoint (bp) `[integer]`
    ///
    /// For odd sizes, the midpoint base is centered; for even sizes, ties are resolved
    /// by randomly assigning either the left or right base (to avoid rounding bias).
    ///
    /// Only one of `--resize` and `--flank` can be specified at a time.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser = clap::value_parser!(u32).range(1..),
            requires = "chrom-sizes",
            help_heading = "Resizing / flanking"
        )
    )]
    pub resize: Option<u32>,

    /// Flank original window by the given left and right sizes (bp) `[integers]`
    ///
    /// Example: `--flank 0 100` to extend only to the right.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser = clap::value_parser!(i32),
            num_args = 2,
            requires = "chrom-sizes",
            help_heading = "Resizing / flanking"
        )
    )]
    pub flank: Option<Vec<i32>>,

    /// Chromosome sizes file (FAI or two-column sizes) `[path]`
    ///
    /// Required when either `--resize` or `--flank` are specified.
    #[cfg_attr(
        feature = "cli",
        clap(long, value_parser, help_heading = "Resizing / flanking")
    )]
    pub chrom_sizes: Option<PathBuf>,

    /// Policy for windows going out of bounds after transform `[string]`
    ///
    /// - "drop": Drop out-of-bounds windows (default).
    ///
    /// - "trim": Trim to chromosome bounds.
    ///
    /// - "allow": Allow out-of-bounds (unsafe).
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "drop",
            ignore_case = true,
            help_heading = "Resizing / flanking"
        )
    )]
    pub oob: OobPolicy,

    // ─────────────────────────────────────────────────────────────────────────────
    // Group-wise filters
    // ─────────────────────────────────────────────────────────────────────────────
    /// Minimum number of windows required per group `[integer]`
    ///
    /// Groups with fewer than this number of windows *after all other filtering and merging steps* are dropped.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser = clap::value_parser!(u32),
            help_heading = "Group-wise filters"
        )
    )]
    pub min_per_group: Option<u32>,

    /// Minimum spacing between windows within the same group (bp) `[integer]`
    ///
    /// Selection rule: Sort windows by `(chrom, start, end)`. Keep a window if its start
    /// is at least `min-distance-within-group` bp after the end of the last kept window on
    /// the same chromosome within the same group; otherwise, drop it. Ties are resolved by
    /// `--distance-ties`.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser = clap::value_parser!(u32),
            help_heading = "Group-wise filters"
        )
    )]
    pub min_distance_within_group: Option<u32>,

    /// How to resolve distance ties when enforcing `--min-distance-within-group` `[string]`
    ///
    /// - "keep-first": Keep the first; skip subsequent windows within distance.
    ///
    /// - "keep-highest-score": Prefer higher score (requires `--score-col`).
    ///
    /// - "keep-lowest-score": Prefer lower score (requires `--score-col`).
    ///
    /// - "keep-longest": Prefer longer windows.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "keep-first",
            ignore_case = true,
            help_heading = "Group-wise filters"
        )
    )]
    pub distance_ties: DistanceTiesPolicy,

    // TODO: Is this not itself a duplicate of min_distance_within_group=1 ?
    /// Deduplication policy for identical intervals within a group `[string]`
    ///
    /// Dedup acts only on **identical** `(chrom,start,end,group)` windows. It is different from
    /// `--min-distance-within-group` which considers physical spacing. Use dedup to collapse
    /// duplicated records; use min-distance to enforce spacing.
    ///
    /// - "none": No deduplication.
    ///
    /// - "keep-first": Keep the first occurrence.
    ///
    /// - "keep-highest-score": Prefer the window with the highest score (requires `--score-col`).
    ///
    /// - "keep-lowest-score": Prefer the window with the lowest score (requires `--score-col`).
    ///
    /// - "keep-longest": Keep the longest.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "none",
            ignore_case = true,
            help_heading = "Group-wise filters"
        )
    )]
    pub deduplicate: DedupKeep,

    /// Merging scope for nearby windows `[string]`
    ///
    /// - "none": Do not merge windows.
    ///
    /// - "within": Merge only windows belonging to the same group.
    ///
    /// - "across": Merge regardless of group (labels resolved by `--merge-label`).
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "none",
            ignore_case = true,
            help_heading = "Group-wise filters"
        )
    )]
    pub merge_scope: MergeScope,

    /// Maximum gap (bp) between windows to be merged `[integer]`
    ///
    /// Required when `--merge-scope` is "within" or "across".
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser = clap::value_parser!(u32).range(0..),
            requires = "merge-scope",
            help_heading = "Group-wise filters"
        )
    )]
    pub merge_gap: Option<u32>,

    /// Label policy when merging `[string]`
    ///
    /// - "join": Join labels with `__` (default).
    ///
    /// - "first": Keep the first label encountered.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "join",
            ignore_case = true,
            help_heading = "Group-wise filters"
        )
    )]
    pub merge_label: MergeLabel,

    // ─────────────────────────────────────────────────────────────────────────────
    // Reproducibility
    // ─────────────────────────────────────────────────────────────────────────────
    /// Seed for any randomized operations `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(long, value_parser = clap::value_parser!(u64), help_heading = "Reproducibility")
    )]
    pub seed: Option<u64>,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "cli", derive(ValueEnum))]
pub enum HeaderMode {
    Auto,
    Present,
    Absent,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "cli", derive(ValueEnum))]
pub enum MissingScore {
    Keep,
    Drop,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "cli", derive(ValueEnum))]
pub enum NearEdge {
    /// Use left genomic edge only.
    Left,
    /// Use right genomic edge only.
    Right,
    /// Use whichever genomic edge is closer.
    Nearest,
    /// Use the edge that is upstream of the near interval given its annotated strand orientation.
    Upstream,
    /// Use the edge that is downstream of the near interval given its annotated strand orientation.
    Downstream,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "cli", derive(ValueEnum))]
pub enum NearDirection {
    Upstream,
    Downstream,
    Both,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "cli", derive(ValueEnum))]
pub enum NearTiePolicy {
    #[cfg_attr(feature = "cli", value(name = "annotate"))]
    Annotate,
    #[cfg_attr(feature = "cli", value(name = "drop"))]
    Drop,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "cli", derive(ValueEnum))]
pub enum OobPolicy {
    Drop,
    Trim,
    Allow,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "cli", derive(ValueEnum))]
pub enum DistSign {
    Absolute,
    Signed,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "cli", derive(ValueEnum))]
pub enum DistanceTiesPolicy {
    #[cfg_attr(feature = "cli", value(name = "keep-first"))]
    KeepFirst,
    #[cfg_attr(feature = "cli", value(name = "keep-highest-score"))]
    KeepHighestScore,
    #[cfg_attr(feature = "cli", value(name = "keep-lowest-score"))]
    KeepLowestScore,
    #[cfg_attr(feature = "cli", value(name = "keep-longest"))]
    KeepLongest,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "cli", derive(ValueEnum))]
pub enum DedupKeep {
    None,
    #[cfg_attr(feature = "cli", value(name = "keep-first"))]
    KeepFirst,
    #[cfg_attr(feature = "cli", value(name = "keep-highest-score"))]
    KeepHighestScore,
    #[cfg_attr(feature = "cli", value(name = "keep-lowest-score"))]
    KeepLowestScore,
    #[cfg_attr(feature = "cli", value(name = "keep-longest"))]
    KeepLongest,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "cli", derive(ValueEnum))]
pub enum MergeScope {
    None,
    Within,
    Across,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "cli", derive(ValueEnum))]
pub enum MergeLabel {
    Join,
    First,
}

impl Default for PrepareConfig {
    fn default() -> Self {
        Self {
            input: "-".into(),
            output: "-".into(),
            header: HeaderMode::Auto,
            cols: "chrom=0,start=1,end=2".to_string(),
            group_cols: Vec::new(),
            score_col: None,
            sep: '\t',
            score_filter: None,
            score_missing: MissingScore::Keep,
            blacklist: None,
            blacklist_halo: 0,
            blacklist_strategy: BlacklistStrategy::Any,
            near: None,
            near_header: HeaderMode::Auto,
            near_edge: NearEdge::Nearest,
            near_direction: NearDirection::Both,
            near_ties: NearTiePolicy::Annotate,
            distance_sign: DistSign::Absolute,
            distance_bins: None,
            distance_max: None,
            resize: None,
            flank: None,
            chrom_sizes: None,
            oob: OobPolicy::Drop,
            min_per_group: None,
            min_distance_within_group: None,
            distance_ties: DistanceTiesPolicy::KeepFirst,
            deduplicate: DedupKeep::None,
            merge_scope: MergeScope::None,
            merge_gap: None,
            merge_label: MergeLabel::Join,
            seed: None,
        }
    }
}

fn parse_sep(input: &str) -> Result<char, String> {
    match input {
        r"\t" | "tab" => Ok('\t'),
        r"\n" | "nl" => Ok('\n'),
        r"\0" | "nul" => Ok('\0'),
        r"\r" => Ok('\r'),
        r"\x20" => Ok(' '),
        _ => {
            let mut it = input.chars();
            let ch = it.next().ok_or("separator cannot be empty")?;
            if it.next().is_some() {
                return Err("expected a single character or an escape like \\t".into());
            }
            Ok(ch)
        }
    }
}
