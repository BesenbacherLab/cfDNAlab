use crate::{
    commands::cli_common::ChromosomeArgs,
    commands::prepare_windows::{labels::validate_label_token, near_file::NearDuplicatesPolicy},
    shared::{blacklist::BlacklistStrategy, thread_pool::default_thread_count},
};
#[cfg(feature = "cli")]
use clap::{ArgGroup, Parser, ValueEnum};
use std::path::PathBuf;

/// Clean and standardise genomic windows for downstream cfDNAlab tools.
///
/// `prep-windows` reads a delimited table with `chrom,start,end` and optional
/// group and score columns, filters and validates each row, and writes a BED-like file for
/// downstream tools. The command writes validated coordinates plus specifiable
/// label columns.
///
/// A label is a tag that tells downstream analyses how to partition the windows. Labels can be
/// based on input columns, distance to a `--near` set of genomic intervals, overlap-based
/// "cluster" inclusion, or compositions of these parts that you define.
///
/// The command parses the BED/TSV/CSV input and:
///
/// - Filters windows using score thresholds, blacklist overlap, and label-based rules.
///
/// - Adjusts coordinates by resizing to a specific size or adding flanks to the
///   current sizes (trimmed to chromosome limits).
///
/// - Merges windows, either within groups or across groups.
///
/// - Builds labels by combining input columns and/or binning distances to elements
///   in the `near` file (e.g., proximal/distal TSS sites).
///
/// - Tags dense overlaps as "clusters".
///
/// The output is minimal, headerless, and sorted by `(chrom, start, end, labels)`,
/// where the label columns are specified via `--out-labels` and the chromosome
/// order is controlled by `--chromosomes`.
/// When `--chromosomes all` is specified, the output order follows the chromosome
/// order in the input.
///
/// ## Practical notes
///
/// - A temporary directory is created during processing and deleted afterwards.
///   When using `stdin`/`stdout`, ensure the working directory has read+write permissions.
///
/// - All coordinates are 0-based half-open `[start, end)`.
///
/// - Column indices are 0-based when you refer to them explicitly.
///
/// - Blacklist checks run on the resized window span using the halo (padding) you configure.
///   When no resize or flank is configured, "resized coordinates" match the original coordinates.
///
/// - "Nearest distance" refers to the closest considered edge of the comparison interval. NOTE:
///   For point features (e.g., TSS), provide 1-bp intervals at the strand-specific coordinate.
///
/// ## Examples
///
/// ### Case: Preparing binding sites for `cfdna midpoints`
///
/// Input: Transcription factor binding sites
///
/// Filters:
///
///   - Remove windows overlapping halo-expanded blacklisted regions ("halo": we flank blacklisted
///     regions with half the maximum fragment size (from downstream analysis) on both sides,
///     so midpoint profiles are not affected by problematic neighbouring regions.
///
///   - Deduplicate binding sites with same (chrom,start,end,transcription factor).
///
///   - Merge overlapping binding sites belonging to the same transcription factor
///     (using the original coordinates). These could be slightly shifted duplicates
///     from multiple experiments.
///
///   - Proximity to TSS sites (`--near` and `--distance-bins`). Keep only distal binding sites
///     within 100kb from a TSS site (`--distance-max`).
///     
///   - Min-per: At least 1000 windows per transcription factor (i.e., the `input` **group**).
///
///   - Detect "clusters" of binding sites when >15 binding sites overlap (after merging per
///     transcription factor). These are marked in the output group label, so we can count
///     midpoints separately for clustered and non-clustered binding sites.
///
/// We resize windows to 2001bp around the midpoint of the intervals.
///
/// We create an output label column that concatenates the input transcription factor
/// name and the clustering-status, so we can use that as group column in `cfdna midpoints`.
///
/// ```bash
///
/// cfdna prep-windows -i tf_coords.bed.gz -o new_tf_coords.tsv.zst \
///
///     --sep tab --group-cols 3  \
///
///     --compose out=input,cluster --out-labels out  \
///
///     --exclude-labels bin=prox --min-per input=1000  \
///
///     --blacklist encode_blacklist.bed --blacklist-halo 500  \
///
///     --near tss_coords.bed --near-strand-col 4 --near-duplicates keep-first  \
///
///     --distance-bins 'prox:<=2000' 'dist:>2000'  \
///
///     --distance-max 100000 --distance-from original  \
///
///     --resize 2001 --chrom-sizes hg38.chrom.sizes  \
///
///     --deduplicate keep-first \
///
///     --cluster-min-overlaps 15  \
///
///     --merge-scope within --merge-gap 0
///
/// ```
///
/// Alternatively, define proximal as e.g. 2.5kb upstream and 0.5kb downstream:
///
/// ```bash
///
///     --distance-bins 'prox:-2500-500' 'upDist:<-2500' 'downDist:>500'  \
///
///     --distance-sign signed
///
/// ```
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
    /// Input BED-like file `[path]`
    ///
    /// Compression inferred from file extension (`.gz` or `.zst`). Use `'-'` for stdin.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            short = 'i',
            value_parser,
            required = true,
            help_heading = "Core"
        )
    )]
    pub input: PathBuf,

    /// Output BED-like file `[path]`
    ///
    /// Compression inferred from file extension (`.gz` or `.zst`). Use `'-'` for stdout.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            short = 'o',
            value_parser,
            default_value = "-",
            help_heading = "Core"
        )
    )]
    pub output: PathBuf,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

    /// Header presence in input `[string]`
    ///
    /// - `"auto"`: Infer from first line.
    ///
    /// - `"present"`: First line is a header line with column names.
    ///
    /// - `"absent"`: No header. Only indices allowed in `--cols` and related.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "auto",
            ignore_case = true,
            hide_possible_values = true,
            help_heading = "Core"
        )
    )]
    pub header: HeaderMode,

    /// Field separator for input and output `[char]`
    ///
    /// Common separators are `\t` (accepts `tab`) for `.tsv` files and `,` or `;` for `.csv` files.
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

    /// Optional group columns `[strings]`
    ///
    /// Provide one or more column indices. When multiple are given, they will be
    /// concatenated using `__` into the `input` label.
    /// Values must be ASCII alphanumerics and cannot be `none`.
    /// Empty values are replaced with `[NA]` so the number of group segments stays fixed.
    ///
    /// If omitted, the `input` label is empty and can be composed with later labels.
    ///
    /// If no subdivision occurs and `--out-labels` is not set, the output still includes `input`.
    ///
    /// Example: `--group-cols 3`
    #[cfg_attr(feature = "cli", clap(long, value_parser, help_heading = "Core"))]
    pub group_cols: Vec<String>,

    /// Optional score column `[string]`
    ///
    /// Column index for a numeric score used by `--score-filter`.
    ///
    /// Example: `--score-col 4`
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser,
            requires = "score-filter",
            help_heading = "Score filtering"
        )
    )]
    pub score_col: Option<String>,

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
    /// - `"keep"`: Keep records with missing/invalid scores.
    ///
    /// - `"drop"`: Drop records with missing/invalid scores.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "keep",
            ignore_case = true,
            hide_possible_values = true,
            help_heading = "Score filtering"
        )
    )]
    pub score_missing: MissingScore,

    /// Define a named label composition `[string]`
    ///
    /// Use this to name a label that joins parts with dots in the order given.
    /// Parts can be atomic parts or earlier compositions.
    /// Names must be ASCII alphanumerics and cannot be `none`.
    ///
    /// **Format**:
    ///
    /// - `NAME=PART1,PART2,...`
    ///
    /// Example: `--compose core=input,bin`
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser,
            num_args = 1..,
            action = clap::ArgAction::Append,
            help_heading = "Labels and filtering"
        )
    )]
    pub compose: Vec<ComposeSpec>,

    /// Label columns to write after coordinates `[strings]`
    ///
    /// Use this to pick which labels are written, including atomic parts and named compositions.
    /// If omitted, the output includes `input` only.
    ///
    /// Rows are ordered by `chrom`, `start`, `end`, and then these label columns.
    ///
    /// Example: `--out-labels input bin`
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser,
            num_args = 1..,
            default_value = "input",
            help_heading = "Labels and filtering"
        )
    )]
    pub out_labels: Vec<String>,

    /// Minimum number of windows required per label key `[strings]`
    ///
    /// Use this to enforce minimum counts for atomic parts or named compositions.
    ///
    /// **Format**:
    ///
    /// - `KEY=COUNT`
    ///
    /// Example: `--min-per input=1000 core=250`
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser,
            num_args = 1..,
            action = clap::ArgAction::Append,
            help_heading = "Labels and filtering"
        )
    )]
    pub min_per: Vec<String>,

    /// Drop windows whose labels include these terms `[strings]`
    ///
    /// Use this to exclude windows based on atomic parts or compositions before
    /// any `--min-per` filtering.
    ///
    /// **Format**:
    ///
    /// - `KEY=TERM`
    ///
    /// **Examples**:
    ///
    /// - `--exclude-labels bin=prox cluster=cluster`
    ///
    /// - `--exclude-labels cluster=none` to keep only in-cluster windows
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser,
            num_args = 1..,
            action = clap::ArgAction::Append,
            help_heading = "Labels and filtering"
        )
    )]
    pub exclude_labels: Vec<String>,

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
    ///
    /// E.g., the maximum fragment length to avoid interval counts
    /// being affected by neighbouring blacklisted regions.
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
    /// - `"proportion=<threshold>"`: Overlap proportion with respect to the window is >= threshold.
    ///
    /// Example: `--blacklist-strategy proportion=0.2`
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            alias = "bl-strategy",
            default_value = "any",
            ignore_case = true,
            hide_possible_values = true,
            help_heading = "Blacklist"
        )
    )]
    pub blacklist_strategy: BlacklistStrategy,

    /// BED-like file with target intervals to compute nearest distance to `[path]`
    ///
    /// Each window can be annotated and filtered based on its distance to the nearest
    /// "near-interval", using specifiable distance bins.
    /// For signed distances, the distance is negative when the window lies upstream
    /// of the near-interval and positive when downstream.
    ///
    /// A common usecase is to annotate/filter windows based on the proximity
    /// to TSS sites (1-bp intervals with defined strands).
    ///
    /// Expected columns:
    ///   `chromosome`, `start`, `end`, optional `strand`, optional `group`.
    ///
    /// Use `--near-strand-col` to point to a strand column.
    /// When omitted, all intervals default to the `+` strand.
    ///
    /// Use `--near-group-cols` to point to group columns.
    ///
    /// - `strand` is one of `+`, `-`, or `.` (unknown). If absent, `+` is assumed.
    ///
    /// - Group names must be ASCII alphanumerics and cannot be the string `none`.
    ///
    /// - Intervals must be half-open (start=inclusive, end=exclusive).
    ///   Duplicate edges are resolved by `--near-duplicates`.
    ///
    /// **Distance**:
    ///
    ///   - Overlap: Distance is `0`.
    ///
    ///   - Otherwise: Distance is the minimum from the window's considered edges
    ///     to the selected target edge(s) of the near-intervals (see `--near-edge`).
    ///
    /// **Labels**:
    ///
    /// Supplying near-intervals makes the labels `win-direction` and `near-name`
    /// (when `--near-group-cols` is specified) available
    /// for use in `--compose`, `--out-labels`, and `--exclude-labels`.
    ///
    /// Win-direction label prefix (strand-aware):
    ///
    ///   `-` = window lies upstream of the near-interval
    ///
    ///   `+` = window lies downstream of the near-interval
    ///
    ///   `=` = window overlaps the near-interval
    ///
    /// Signed distances use the same orientation: the window lying upstream of
    /// near-interval has a negative distance, downstream is positive.
    ///
    /// The prefix is still included when `--distance-sign absolute` is selected.
    ///
    /// If the left and right intervals tie and `--near-ties annotate`, both sides are reported on the
    /// same window by expanding its labels. Output columns show paired comma-separated values
    /// (e.g., `win-direction: -,+` and `near-name: GeneA,GeneB`). When composing these labels,
    /// the side+name pairs are joined first, e.g.:
    /// `--compose near=win-direction,near-name` -> `-.GeneA,+.GeneB`. Note that the two
    /// intervals can have different strands, so `win-direction` can also be any combination of
    /// `+` and `-`, e.g. `+,+`.
    ///
    /// Output labels **compact identical values** after merges and tie
    /// expansion: if every surviving tuple has the same `win-direction`, it is written once
    /// (e.g., `+,+` -> `+`). If any tuple differs (e.g., a merged `-`), the list is kept. Use a
    /// composition such as `--compose near=win-direction,near-name --out-labels input near` when
    /// you need the side–name pairings preserved regardless of this compaction.
    ///
    /// **No near-interval on chromosome or in targeted direction**:
    ///
    /// When no near-intervals exist for a chromosome:
    ///
    ///   - If `--distance-max` is set, windows on that chromosome are dropped. Otherwise:
    ///
    ///   - If `--distance-bins` is set, `bin` is set to `[NO-NEAR]`.
    ///
    ///   - `win-direction` (and `near-name` when configured) are set to `[NONE]`.
    ///
    /// When a chromosome has near-intervals but none match the configured `--near-direction`
    /// or `--near-edge` for a given window:
    ///
    /// - If `--distance-max` is set, the window is dropped with a warning.
    ///
    /// - Otherwise, the window is kept (near labels remain empty).
    ///
    /// **Upstream/Downstream definition (strand-aware)**:
    ///
    /// `Upstream/Downstream` are defined **relative to the near-interval's annotated strand**:
    ///
    /// - For a `+`-stranded near-interval: upstream means the window is to the left
    ///   (smaller genomic coordinates) of the near-interval; downstream is to the right.
    ///
    /// - For a `-`-stranded near-interval: upstream is to the right (larger genomic coordinates); downstream is to the left.
    ///
    /// - For an unknown strand (`.`): upstream/downstream are derived from genomic placement
    ///   (falls back to genomic-nearest semantics).
    ///
    /// **Examples**
    ///
    /// Legend: '===' near-interval, '###' window, '---' empty span. Signs are relative to the near-interval's strand.
    ///
    /// Case A: near is on `+`-strand
    ///
    /// ```text
    ///
    /// | coordinates:  100   120  140         200   220   240
    ///
    /// |               |#####|----|===========|-----|#####|
    ///
    /// |   upstream (-) ^^^^^         near           ^^^^^ downstream (+)
    ///     
    /// ```
    ///
    /// Case B: near is on `-`-strand
    ///
    /// ```text
    ///
    /// | coordinates:  100   120  140         200   220   240
    ///
    /// |               |#####|----|===========|-----|#####|
    ///
    /// | downstream (+) ^^^^^         near           ^^^^^ upstream (-)
    ///     
    /// ```
    ///
    /// Case C: overlap
    ///     
    /// ```text
    ///
    /// | coordinates:  100   120             200
    ///
    /// |               |--###|===========###====|
    ///
    /// |        touch (=) ^^^    near    ^^^ overlap (=)
    ///
    /// ```
    #[cfg_attr(
        feature = "cli",
        clap(long, value_parser, help_heading = "Distance to near-intervals")
    )]
    pub near: Option<PathBuf>,

    /// Header presence in the `--near` file `[string]`
    ///
    /// - `"auto"`: Infer from first line.
    ///
    /// - `"present"`: First line is a header line with column names.
    ///
    /// - `"absent"`: No header. Only indices allowed when referencing columns from the near file.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "auto",
            ignore_case = true,
            hide_possible_values = true,
            help_heading = "Distance to near-intervals"
        )
    )]
    pub near_header: HeaderMode,

    /// Strand column index in the `--near` file `[string]`
    ///
    /// Use this when the near file includes a strand column.
    /// Index is 0-based.
    /// When omitted, all intervals default to the `+` strand.
    ///
    /// Example: `--near-strand-col 3`
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser,
            requires = "near",
            help_heading = "Distance to near-intervals"
        )
    )]
    pub near_strand_col: Option<String>,

    /// Group column indices in the `--near` file `[strings]`
    ///
    /// Use this when the near file includes group name columns.
    /// Indices are 0-based.
    /// When omitted, `near-name` is empty and not available for filtering.
    /// Empty values are replaced with `[NA]` so the number of group segments stays fixed.
    ///
    /// Example: `--near-group-cols 4 5`
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser,
            num_args = 1..,
            requires = "near",
            help_heading = "Distance to near-intervals"
        )
    )]
    pub near_group_cols: Vec<String>,

    /// Edge of near-intervals to consider in distance calculation `[string]`
    ///
    /// - `"left"`: Use left genomic edge only.
    ///
    /// - `"right"`: Use right genomic edge only.
    ///
    /// - `"nearest"`: Use whichever genomic edge is closer (default).
    ///
    /// - `"upstream"`: Use the edge towards the upstream side of the near-interval,
    ///   given its strand (`+`-strand uses left edge, `-` uses right edge).
    ///
    /// - `"downstream"`: Use the edge towards the downstream side of the near-interval,
    ///   given its strand (`+`-strand uses right edge, `-` uses left edge).
    ///
    /// If a near-interval's strand is unknown (`.`), "upstream"/"downstream" fall back to `"nearest"`.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "nearest",
            ignore_case = true,
            hide_possible_values = true,
            help_heading = "Distance to near-intervals"
        )
    )]
    pub near_edge: NearEdge,

    /// Which near-intervals are considered (upstream/downstream/both) `[string]`
    ///
    /// - `"upstream"`: Keep hits where the window is upstream of the near-interval (or overlaps it).
    ///
    /// - `"downstream"`: Keep hits where the window is downstream of the near-interval (or overlaps it).
    ///
    /// - `"both"`: Keep upstream, downstream, and overlaps (default).
    ///
    /// See `--near` for definitions of up-/downstream.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "both",
            ignore_case = true,
            hide_possible_values = true,
            help_heading = "Distance to near-intervals"
        )
    )]
    pub near_direction: NearDirection,

    /// How to respond when multiple near-intervals tie for the minimum distance `[string]`
    ///
    /// - `"annotate"`: keep the window and include both sides as (left interval, right interval).
    ///   `win-direction` becomes one of {`-,+`, `+,-`, `-` short for `-,-`, `+` short for `+,+`}
    ///   and `near-name` becomes `A,B`. When combined in a composition, they are paired
    ///   as e.g. '+A,+B'. **NOTE**: The `+,+ => +` compaction only happens in the final output annotation
    ///   when all values of `win-direction` are the same.
    ///
    /// - `"drop"`: discard the window when a tie occurs.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "annotate",
            ignore_case = true,
            hide_possible_values = true,
            help_heading = "Distance to near-intervals"
        )
    )]
    pub near_ties: NearTiePolicy,

    /// Policy for identical near-interval edges `[string]`
    ///
    /// Identical edges are records on the same chromosome with the same `(start, end)` and,
    /// when a strand column is provided, the same `strand`.
    ///
    /// Multiple groups at the exact same site create an ambiguous "nearest" unless resolved.
    ///
    ///  - `"error"`: Fail on identical edges with a descriptive message.
    ///
    ///  - `"keep-first"`: Keep the first record in each run of duplicates. Drop the rest.
    ///
    ///  - `"drop-all"`: Drop the entire set of duplicates.
    ///
    ///  - `"merge"`: Merge groups across identical edges (and strands when present) into one record.
    ///    Group names are joined with "`__`", with duplicates removed. Missing groups are ignored.
    ///
    /// Key used to detect “identical edges”:
    ///
    /// - When `--near-strand-col` is set, duplicates are keyed by `(start, end, strand)`.
    ///
    /// - When no strand column is present, duplicates are keyed by `(start, end)` while `strand` is ignored.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "error",
            ignore_case = true,
            hide_possible_values = true,
            help_heading = "Distance to near-intervals"
        )
    )]
    pub near_duplicates: NearDuplicatesPolicy,

    // TODO: Add "rest" definition
    // TODO: Ensure 'prox:-2500-500' works (when signed)
    /// Distance bin specifications `[quoted strings]`
    ///
    /// Provide one or more `'<label>:<expr>'` rules. The first matching rule wins.
    /// Expression forms: `<N`, `<=N`, `A-B`, `>=N`, `>N` (N in bp).
    /// Labels must be ASCII alphanumerics and cannot be the string `none`.
    ///
    /// When specified, the label `bin` becomes available.
    /// Chromosomes without near-intervals use the special bin `[NO-NEAR]`.
    ///
    /// **Examples**:
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
            requires = "near",
            help_heading = "Distance to near-intervals"
        )
    )]
    pub distance_bins: Option<Vec<String>>,

    /// How to treat the computed distances when binning `[string]`
    ///
    /// - `"absolute"`: Use `abs(distance)` for comparisons and thresholds.
    ///
    /// - `"signed"`: Use signed distances.
    ///
    /// **Distance sign (when `--distance-sign signed`):**
    ///
    /// - Window is **upstream** of the near-interval -> **negative** distance.
    ///
    /// - Window is **downstream** of the near-interval -> **positive** distance.
    ///
    /// - Overlap/touch -> `0` distance.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "absolute",
            ignore_case = true,
            requires = "near",
            help_heading = "Distance to near-intervals"
        )
    )]
    pub distance_sign: DistSign,

    /// Maximum absolute distance to keep `[integer]`
    ///
    /// Windows with `|distance| > --distance-max` are dropped prior to binning.
    /// When specified, windows without any near-intervals are dropped.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser = clap::value_parser!(u32),
            requires = "near",
            help_heading = "Distance to near-intervals"
        )
    )]
    pub distance_max: Option<u32>,

    /// Coordinates used for distance binning `[string]`
    ///
    /// Use this to choose which coordinates determine the near distance and bin label.
    /// When no resize or flank is configured, resized coordinates match the originals.
    ///
    /// One of:
    ///
    /// - `"resized"`: Use resized coordinates.
    ///
    /// - `"original"`: Use original coordinates.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "resized",
            ignore_case = true,
            hide_possible_values = true,
            help_heading = "Distance to near-intervals"
        )
    )]
    pub distance_from: CoordinateSet,

    /// Resize window to a fixed size (bp) centered on the midpoint `[integer]`
    ///
    /// Resizing centers the new window on the midpoint of the original interval.
    /// When the input length and target size have different parity (odd/even),
    /// there are two equally centered placements. In that case, the code chooses
    /// left or right with a deterministic hash based on the midpoint, input length,
    /// target size, and optional seed.
    ///
    /// The interval and resize combinations look like the following (2+3 leads to
    /// random selection of midpoint to reduce midpoint selection bias):
    ///
    /// ```text
    /// Interval size 6, resize 4: [011110] -> unique placement
    /// Interval size 6, resize 3: [001110] or [011100] -> left or right choice
    /// Interval size 5, resize 4: [11110] or [01111] -> left or right choice
    /// Interval size 5, resize 3: [01110] -> unique placement
    /// ```
    ///
    /// Only one of `--resize` and `--flank` can be specified at a time.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser = clap::value_parser!(u32).range(1..),
            requires = "chrom-sizes",
            help_heading = "Resizing / flanking (select max. one transformation)"
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
            help_heading = "Resizing / flanking (select max. one transformation)"
        )
    )]
    pub flank: Option<Vec<i32>>,

    /// Chromosome sizes file (FAI or two-column sizes) `[path]`
    ///
    /// Required when either `--resize` or `--flank` is specified. When no size transformation
    /// is specified, the windows are only checked for out-of-bounds when
    /// `--chrom-sizes` is specified (otherwise ignored).
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser,
            help_heading = "Resizing / flanking (select max. one transformation)"
        )
    )]
    pub chrom_sizes: Option<PathBuf>,

    /// Policy for windows going out of bounds `[string]`
    ///
    /// Negative coordinates are always caught.
    /// Coordinates beyond the chromosome length
    /// are checked when `--chrom-sizes` is specified.
    ///
    /// - `"drop"`: Drop out-of-bounds windows (default).
    ///
    /// - `"trim"`: Trim to chromosome bounds.
    ///
    /// - `"allow"`: Allow out-of-bounds (unsafe).
    ///   Windows with negative coordinates are dropped with a warning.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "drop",
            ignore_case = true,
            hide_possible_values = true,
            help_heading = "Resizing / flanking (select max. one transformation)"
        )
    )]
    pub oob: OobPolicy,

    /// Minimum distance between windows within the same group (bp) `[integer]`
    ///
    /// This drops windows that are too close to the previous kept window
    /// within the same **input group** on the same chromosome.
    ///
    /// It runs after within-group merging and before across-group merging.
    /// Use `--cluster-before-min-distance` to move clustering ahead of this step.
    /// It uses the same coordinate set as `--distance-from`.
    ///
    /// **Selection rule**:
    ///
    /// - Sort windows by `(input, chrom, start, end)` in the chosen coordinate set.
    ///
    /// - Walk windows in order per group and collect runs that violate the minimum distance.
    ///
    /// - Pick one window per run using `--distance-policy`, then continue from the chosen window's end.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser = clap::value_parser!(u32),
            help_heading = "Group-wise filters"
        )
    )]
    pub min_distance_within_group: Option<u32>,

    /// How to choose a window when a run of overlaps violates the minimum distance `[string]`
    ///
    /// - `"keep-first"`: Keep the first window and skip subsequent windows within distance.
    ///
    /// - `"keep-highest-score"`: Prefer higher score (requires `--score-col`).
    ///
    /// - `"keep-lowest-score"`: Prefer lower score (requires `--score-col`).
    ///
    /// - `"keep-longest"`: Prefer longer windows.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "keep-first",
            ignore_case = true,
            hide_possible_values = true,
            help_heading = "Group-wise filters"
        )
    )]
    pub distance_policy: DistancePolicy,

    /// Deduplication policy for identical intervals within a group `[string]`
    ///
    /// Windows are considered duplicates when they have identical `(chrom,start,end,input)`
    /// using the input group label.
    ///
    /// Deduplication runs before merging. It uses resized coordinates when resizing
    /// or flanking is enabled, otherwise original coordinates.
    ///
    /// - `"none"`: No deduplication.
    ///
    /// - `"keep-first"`: Keep the first occurrence.
    ///
    /// - `"keep-highest-score"`: Prefer the window with the highest score (requires `--score-col`).
    ///
    /// - `"keep-lowest-score"`: Prefer the window with the lowest score (requires `--score-col`).
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "none",
            ignore_case = true,
            hide_possible_values = true,
            help_heading = "Group-wise filters"
        )
    )]
    pub deduplicate: DedupKeep,

    /// Minimum overlapping windows to tag a window as a cluster `[integer]`
    ///
    /// A window is marked as a cluster when its average position-wise window overlap
    /// meets or exceeds this value after within-group merging, counting itself.
    /// The average is the total overlap depth across the window divided by its length.
    /// Overlap is evaluated across groups on the same chromosome.
    /// Use `--cluster-before-min-distance` to move clustering ahead of the
    /// minimum-distance filter. Across-group merging always happens after clustering.
    ///
    /// Use this to label dense regions so they can be filtered or stratified later.
    /// Non-cluster windows use the label value `none`.
    /// If omitted, cluster labels are not added.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser = clap::value_parser!(u32).range(1..),
            help_heading = "Labels and filtering"
        )
    )]
    pub cluster_min_overlaps: Option<u32>,

    /// Coordinates used for clustering overlap checks `[string]`
    ///
    /// One of:
    ///
    /// - `"original"`: Use original coordinates.
    ///
    /// - `"resized"`: Use resized coordinates. When no resize or flank is configured,
    ///   resized coordinates match the originals.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "original",
            ignore_case = true,
            hide_possible_values = true,
            help_heading = "Labels and filtering"
        )
    )]
    pub cluster_on: CoordinateSet,

    /// Compute cluster labels before `--min-distance-within-group`
    ///
    /// When unset, clustering runs after the minimum-distance filter.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            action = clap::ArgAction::SetTrue,
            help_heading = "Labels and filtering"
        )
    )]
    pub cluster_before_min_distance: bool,

    /// Merging scope for nearby windows `[string]`
    ///
    /// - `"none"`: Do not merge windows.
    ///
    /// - `"within"`: Merge only windows that share the same `--merge-key` value.
    ///
    /// - `"across"`: Merge regardless of labels (labels resolved by `--merge-label`).
    ///
    /// Within-group merges use `--merge-key` to decide which labels are grouped together.
    /// Across-group merges run after clustering and minimum-distance filtering.
    ///
    /// Merging requires specifying `--merge-gap`.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "none",
            ignore_case = true,
            hide_possible_values = true,
            help_heading = "Merging neighbours"
        )
    )]
    pub merge_scope: MergeScope,

    /// Label key used to define within-group merges `[string]`
    ///
    /// Use an atomic part or a named composition to decide which windows belong together.
    /// This applies only when `--merge-scope` is `"within"`.
    ///
    /// Example: `--merge-key input`
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser,
            default_value = "input",
            help_heading = "Merging neighbours"
        )
    )]
    pub merge_key: String,

    /// Coordinates used for merging `[string]`
    ///
    /// One of:
    ///
    /// - `"original"`: Merge using original coordinates. If resizing is configured,
    ///   the merged window is resized after merging.
    ///
    /// - `"resized"`: Merge using resized coordinates. **No resizing** is performed on the
    ///   merged windows. When no resize or flank is specified, the original coordinates
    ///   are used but no post-merge resizing is performed.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "original",
            ignore_case = true,
            hide_possible_values = true,
            help_heading = "Merging neighbours"
        )
    )]
    pub merge_on: CoordinateSet,

    /// Maximum gap (bp) between windows to be merged `[integer]`
    ///
    /// Use `--merge-gap 0` to merge overlapping/touching windows.
    ///
    /// Must be specified for merging occur. If omitted, merging is skipped even when `--merge-scope` is `"within"`
    /// or `"across"`.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser = clap::value_parser!(u32).range(0..),
            requires = "merge-scope",
            help_heading = "Merging neighbours"
        )
    )]
    pub merge_gap: Option<u32>,

    /// Label policy when merging `[string]`
    ///
    /// This controls how label tuples are combined when multiple windows merge.
    ///
    /// - `"join"`: Keep all label tuples. Output labels may become lists when tuples differ.
    ///
    /// - `"first"`: Keep only the first window's label tuple.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "join",
            ignore_case = true,
            hide_possible_values = true,
            help_heading = "Merging neighbours"
        )
    )]
    pub merge_label: MergeLabel,

    /// Seed for any randomized operations `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(long, value_parser = clap::value_parser!(u64), help_heading = "Reproducibility")
    )]
    pub seed: Option<u64>,

    /// Number of threads to use during sorting steps `[integer]`
    ///
    /// This command generally takes a long time and is hard to parallelize.
    ///
    /// Defaults to the number of available CPU cores (-1).
    #[cfg_attr(
        feature = "cli",
        clap(short = 't', long, default_value_t = default_thread_count(), help_heading = "Core")
    )]
    pub n_threads: usize,
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
    /// Use the edge that is upstream of the near-interval given its annotated strand orientation.
    Upstream,
    /// Use the edge that is downstream of the near-interval given its annotated strand orientation.
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

/// Coordinate set used for window operations.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "cli", derive(ValueEnum))]
pub enum CoordinateSet {
    #[cfg_attr(feature = "cli", value(name = "resized"))]
    Resized,
    #[cfg_attr(feature = "cli", value(name = "original"))]
    Original,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "cli", derive(ValueEnum))]
pub enum DistancePolicy {
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

/// Parsed `--compose` specification.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ComposeSpec {
    pub name: String,
    pub parts: Vec<String>,
}

impl std::str::FromStr for ComposeSpec {
    type Err = String;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        let (raw_name, raw_parts) = input
            .split_once('=')
            .ok_or_else(|| "compose spec must be NAME=PART1,PART2".to_string())?;

        let name = raw_name.trim();
        if name.is_empty() {
            return Err("compose name cannot be empty".to_string());
        }
        if let Err(message) = validate_label_token(name, "compose name") {
            return Err(message);
        }

        let parts: Vec<String> = raw_parts
            .split(',')
            .map(|part| part.trim())
            .filter(|part| !part.is_empty())
            .map(String::from)
            .collect();

        if parts.is_empty() {
            return Err("compose parts cannot be empty".to_string());
        }

        Ok(Self {
            name: name.to_string(),
            parts,
        })
    }
}

impl Default for PrepareConfig {
    fn default() -> Self {
        Self {
            input: "-".into(),
            output: "-".into(),
            chromosomes: ChromosomeArgs::default(),
            header: HeaderMode::Auto,
            cols: "chrom=0,start=1,end=2".to_string(),
            group_cols: Vec::new(),
            out_labels: vec!["input".to_string()],
            compose: Vec::new(),
            min_per: Vec::new(),
            exclude_labels: Vec::new(),
            score_col: None,
            sep: '\t',
            score_filter: None,
            score_missing: MissingScore::Keep,
            blacklist: None,
            blacklist_halo: 0,
            blacklist_strategy: BlacklistStrategy::Any,
            near: None,
            near_header: HeaderMode::Auto,
            near_strand_col: None,
            near_group_cols: Vec::new(),
            near_edge: NearEdge::Nearest,
            near_direction: NearDirection::Both,
            near_ties: NearTiePolicy::Annotate,
            near_duplicates: NearDuplicatesPolicy::Error,
            distance_sign: DistSign::Absolute,
            distance_bins: None,
            distance_max: None,
            distance_from: CoordinateSet::Resized,
            resize: None,
            flank: None,
            chrom_sizes: None,
            oob: OobPolicy::Drop,
            min_distance_within_group: None,
            distance_policy: DistancePolicy::KeepFirst,
            deduplicate: DedupKeep::None,
            cluster_min_overlaps: None,
            cluster_on: CoordinateSet::Original,
            cluster_before_min_distance: false,
            merge_scope: MergeScope::None,
            merge_key: "input".to_string(),
            merge_on: CoordinateSet::Original,
            merge_gap: None,
            merge_label: MergeLabel::Join,
            seed: None,
            n_threads: default_thread_count(),
        }
    }
}

#[cfg_attr(not(feature = "cli"), allow(dead_code))]
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
