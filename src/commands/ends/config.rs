use crate::{
    commands::{
        cli_common::{
            ApplyGCArgs, ChromosomeArgs, DistributionWindowsArgs, FragmentLengthArgs, IOCArgs,
            LoggingArgs, ScaleGenomeArgs, UnpairedArgs,
        },
        ends::config_structs::*,
    },
    shared::{
        blacklist::BlacklistStrategy, constants::DEFAULT_MAX_SOFT_CLIPS,
        indel_mode::IndelMotifFilterPolicy,
    },
};
use std::path::PathBuf;

const ENDS_ABOUT: &str = "Count fragment end- and breakpoint-motifs in a BAM-file.";

const ENDS_LONG_ABOUT: &str = concat!(
    "Count fragment end- and breakpoint-motifs in a BAM-file.\n\n",
    "For each fragment end, it extracts the `--k-outside` bases just outside the fragment and the ",
    "`--k-inside` bases just inside the fragment. For the right fragment end, these are ",
    "reverse-complemented together. Finally, they are combined to the reference 5'->3'-oriented ",
    "`\"<outside>_<inside>\"` motif.\n\n",
    "## Visualization of counting\n\n",
    "The following shows the counting for aligned fragment ends:\n\n",
    r#"For `--k-inside 2 --k-outside 2`:

```text
Reference 5' >>>>>>>>>>>>>>> 3'
             ATCGTTTTTTTCATC
Fragment     --|---------|--
Forward     5' |>>>>>>>| 3'
  Outside    AT
  Inside       CG
Reverse      3' |<<<<<<<<| 5'
  Inside                CA
  Outside                 TC
```

Reverse (`CATC`) is reverse complemented to `GATG`

Counts (`<outside>_<inside>`): `AT_CG: 1`, `GA_TG: 1`
"#,
    "\n",
    "## Output files\n\n",
    "Writes a self-contained `.end_motifs.zarr` store. The store contains either dense ",
    "`counts[row, motif]` when `--all-motifs` is enabled, or sparse COO arrays otherwise.\n\n",
    "Motif labels are saved as `<outside>_<inside>`.\n\n",
    "## GC correction\n\n",
    "Weight the contribution of each fragment based on their GC contents per fragment length.\n\n",
    "## Genomic smoothing (--scaling-factors)\n\n",
    "Weight how genomic regions contribute to the count distribution(s), e.g., to reduce the ",
    "influence of copy number alterations (if that is meaningful to your analysis). ",
    "This weights the contribution of each fragment by region-wise precomputed scaling factors.\n\n",
    "Can be precomputed with `cfdna fragment-count-weights` (recommended) or `cfdna coverage-weights`.\n\n",
    "## Window assignment\n\n",
    "By default, a motif is counted in the window the fragment end falls in with the weight 1.0 (before correction/scaling).\n\n",
    "With `--clip-strategy include-at-shifted-boundary`, that endpoint can move outside the aligned span by the soft-clipped length. ",
    "GC correction and scaling weights still use the aligned reference span.\n\n",
    "With `--clip-strategy include-at-aligned-boundary`, the inside motif includes soft-clipped read bases, but the endpoint assignment stays at the aligned boundary.\n\n",
    "Alternatively, we can weight the motif by how much the fragment overlaps the window or ",
    "we can count both end motifs of a fragment if the *fragment midpoint* or a given ",
    "*proportion* of positions overlaps the window.\n\n",
    "## Blacklisting\n\n",
    "1) Skips fragments that overlap blacklisted regions with a given proportion.\n\n",
    "2) Skips motifs overlapping blacklisted regions.\n\n",
    "Fragment-level blacklist filtering uses the same assignment coordinates as the selected clip strategy. ",
    "With `--clip-strategy include-at-shifted-boundary`, soft-clipped boundary shifts can therefore make a fragment overlap blacklisted regions outside its aligned span.\n\n",
    "With `--clip-strategy include-at-aligned-boundary`, motif-level blacklist validation only checks the part of the inside motif that still overlaps reference coordinates.\n\n",
    "## Always-on exclusion criteria\n\n",
    "The following criteria always exclude a read:\n\n",
    "The read is secondary, supplementary or duplicate. ",
    "The read failed quality check.\n\n",
    "**Paired-end input only**: ",
    "The read or mate read is unmapped. ",
    "The read is mapped to a different `tid` than the mate. ",
    "The paired reads are not inwardly directed (we require: `start(forward) <= start(reverse)`). ",
);

#[cfg_attr(feature = "cli", derive(clap::Args))]
#[cfg_attr(
    feature = "cli",
    clap(about = ENDS_ABOUT, long_about = ENDS_LONG_ABOUT)
)]
#[derive(Clone)]
pub struct EndsConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ioc: IOCArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub unpaired: UnpairedArgs,

    /// Optional prefix for output files (e.g., a sample name) `[string]`
    ///
    /// Leave empty to write filenames without a leading prefix.
    ///
    /// E.g., specify to enable writing to the same output directory from multiple calls to this software.
    ///
    /// Examples produce files like:
    ///   `<prefix>.end_motifs.zarr`
    ///   `<prefix>.end_settings.json`
    #[cfg_attr(
        feature = "cli",
        clap(long, short = 'x', default_value_t = String::new(), hide_default_value = true, value_parser = crate::commands::cli_common::parse_output_prefix, help_heading = "Core")
    )]
    pub output_prefix: String,

    /// 2bit reference genome file [path]
    ///
    /// NOTE: Required when using reference bases, blacklist filtering, or specifying `--gc-file`.
    ///
    /// E.g., "hg38.2bit" from UCSC (https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit).
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

    /// Number of bases to use from inside the fragment `[integer]`
    #[cfg_attr(feature = "cli", clap(long, required = true, help_heading = "Motifs"))]
    pub k_inside: usize,

    /// Number of bases to use from outside the fragment `[integer]`
    #[cfg_attr(feature = "cli", clap(long, required = true, help_heading = "Motifs"))]
    pub k_outside: usize,

    /// Whether to get the inside-fragment bases from the read or the reference `[string]`
    ///
    /// Possible values:
    ///
    /// - `"read"`:
    ///   Use the read sequence for bases inside the fragment.
    ///
    /// - `"reference"`:
    ///   Use the reference genome for bases inside the fragment.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            hide_possible_values = true,
            default_value = "read",
            requires_if("reference", "ref_2bit"),
            help_heading = "Motifs"
        )
    )]
    pub source_inside: KmerSource,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub clip: ClippingArgs,

    /// When to filter motifs due to indels.
    ///
    /// Deletions: Both 'D' and 'N' in the cigar string are considered deletions.
    ///
    /// Possible values:
    ///
    /// - `"auto"`:
    ///   Select the option based on the source.
    ///   
    ///   For **read**-sequence bases, allow indels in the alignment.
    ///
    ///   For **reference** bases, skip motifs with indels in the alignment.
    ///
    /// - `"skip-affected-end"`:
    /// Always skip motifs overlapping indels.
    ///
    /// - `"skip-affected-fragment"`:
    /// Skip **fragments** when either end overlap indels.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "auto",
            ignore_case = true,
            help_heading = "Filtering"
        )
    )]
    pub indel_filter: IndelMotifFilterPolicy,

    /// Size of tiles to parallelize over `[integer]`
    ///
    /// Chromosomes are processed in tiles of this size to reduce memory usage.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "20000000", value_parser = clap::value_parser!(u32).range(1000000..), help_heading="Core"))]
    pub tile_size: u32,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub windows: DistributionWindowsArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub window_assignment: AssignMotifToWindowArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub chromosomes: ChromosomeArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub scale_genome: ScaleGenomeArgs,

    /// Collapse each motif with its same-orientation complement [EXPERIMENTAL] [flag]
    ///
    /// This option is hidden by default and is only shown in CLI help when cfDNAlab is built
    /// with `--features ends_experimental`.
    ///
    /// **NOTE**: In many analyses, this may not be biologically meaningful.
    ///
    /// - (Always) Motifs are oriented so they run from the fragment end inward in 5'->3' direction.
    ///
    /// - Motifs are collapsed with their complement, using the lexicographically smaller motif representation.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            help_heading = "Motifs",
            hide = !cfg!(feature = "ends_experimental")
        )
    )]
    pub collapse_complement: bool,

    /// Include every possible motif in the output, even if its count is zero  [flag]
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Motifs"))]
    pub all_motifs: bool,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub fragment_lengths: FragmentLengthArgs,

    /// Minimum mapping quality to include `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(long, alias = "mq", default_value = "30", value_parser = clap::value_parser!(u8).range(0..), help_heading="Filtering"))]
    pub min_mapq: u8,

    /// Base-quality filter on the inside read bases `[string]`
    ///
    /// Filter either the whole fragment or individual ends based on the base qualities in the **inside** read bases of the motifs.
    ///
    /// Repeat `--bq-filter` to count only ends that pass all **end filters** and belong to fragments that pass all **fragment filters**.
    ///
    /// Examples:
    ///
    /// - `--bq-filter "min in end >= 30"` (for "all bases have decent quality")
    ///
    /// - `--bq-filter "mean in fragment >= 30"` (for "average bases have decent quality")
    ///
    /// - `--bq-filter "max in fragment < 20"` (for "no bases have decent quality")
    ///
    /// Each expression must use:
    ///
    /// - `<agg> in <scope> <op> <threshold>`
    ///
    /// With the following values:
    ///
    /// - with `<agg>` in `min`, `max`, or `mean`
    ///
    /// - with `<scope>` in `end` or `fragment`
    ///
    /// - with `<op>` in `>=`, `>`, `<=`, or `<`
    ///
    /// The keywords are parsed case-insensitively and ASCII whitespace is ignored.
    ///
    /// Scope semantics:
    ///
    /// - `end`: Score each fragment end independently and drop only the failing end.
    ///
    /// - `fragment`: Score the fragment from its two end scores and drop the full fragment when it fails.
    ///
    /// **NOTE**: `--bq-filter` requires `--k-inside > 0` and `--source-inside read`.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser,
            action = clap::ArgAction::Append,
            help_heading = "Filtering"
        )
    )]
    pub bq_filter: Vec<BaseQualityFilter>,

    /// Only count properly paired reads `[flag]`
    ///
    /// This is **NOT** recommended by default as it trims the tails of the length distribution.
    ///
    /// Note, that we only keep inward-directed fragments within the specified length range, so
    /// there's no real need for proper-pair filtering.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Filtering"))]
    pub require_proper_pair: bool,

    /// Optional BED file(s) with blacklisted regions `[path]`
    #[cfg_attr(
        feature = "cli", clap(short = 'b', long, value_parser, num_args = 1.., action = clap::ArgAction::Append, help_heading = "Filtering"))]
    pub blacklist: Option<Vec<PathBuf>>,

    /// Minimum size of blacklist intervals to load (bp) `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            alias = "bl-min-size",
            default_value = "1",
            help_heading = "Filtering"
        )
    )]
    pub blacklist_min_size: u64,

    /// The fragment positions that should overlap blacklisted regions for it to be excluded `[string]`
    ///
    /// **NOTE**: Motifs overlapping blacklisted regions are always skipped. This strategy is
    /// for further filtering of the full fragments. This is useful when you generally don't trust
    /// the reference sequences in blacklisted regions.
    ///
    /// Possible values:
    ///     `"any"`, `"all"`, `"midpoint"`, or `"proportion=<threshold>"`
    ///
    /// Example of proportion: `--blacklist-strategy proportion=0.2` (no space around `=`)
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            alias = "bl-strategy",
            default_value = "any",
            ignore_case = true,
            help_heading = "Filtering"
        )
    )]
    pub blacklist_strategy: BlacklistStrategy,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub gc: ApplyGCArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub logging: LoggingArgs,
}

impl EndsConfig {
    pub fn new(
        ioc: IOCArgs,
        chromosomes: ChromosomeArgs,
        k_inside: usize,
        k_outside: usize,
    ) -> Self {
        Self {
            ioc,
            output_prefix: String::new(),
            ref_2bit: None,
            k_inside,
            k_outside,
            source_inside: KmerSource::Read,
            clip: ClippingArgs {
                clip_strategy: ClipStrategy::Skip,
                max_soft_clips: DEFAULT_MAX_SOFT_CLIPS,
            },
            indel_filter: IndelMotifFilterPolicy::Auto,
            all_motifs: false,
            collapse_complement: false,
            windows: DistributionWindowsArgs::default(),
            window_assignment: AssignMotifToWindowArgs::default(),
            chromosomes,
            scale_genome: ScaleGenomeArgs::default(),
            fragment_lengths: FragmentLengthArgs::default(),
            unpaired: UnpairedArgs {
                reads_are_fragments: false,
            },
            tile_size: 20000000,
            min_mapq: 30,
            bq_filter: Vec::new(),
            require_proper_pair: false,
            blacklist: None,
            blacklist_min_size: 1,
            blacklist_strategy: BlacklistStrategy::default(),
            gc: ApplyGCArgs {
                gc_file: None,
                gc_tag: None,
                neutralize_invalid_gc: false,
            },
            logging: LoggingArgs::default(),
        }
    }

    pub fn set_indel_filter(&mut self, filter: IndelMotifFilterPolicy) {
        self.indel_filter = filter;
    }

    pub fn set_windows(&mut self, windows: DistributionWindowsArgs) {
        self.windows = windows;
    }

    pub fn set_window_assignment(&mut self, assign: AssignMotifToWindowArgs) {
        self.window_assignment = assign;
    }

    pub fn fragment_lengths_mut(&mut self) -> &mut FragmentLengthArgs {
        &mut self.fragment_lengths
    }

    pub fn set_unpaired(&mut self, unpaired: UnpairedArgs) {
        self.unpaired = unpaired;
    }

    pub fn set_scaling_factors(&mut self, scaling_factors: Option<PathBuf>) {
        self.scale_genome.scaling_factors = scaling_factors;
    }

    pub fn set_min_mapq(&mut self, min_mapq: u8) {
        self.min_mapq = min_mapq;
    }

    pub fn set_tile_size(&mut self, tile_size: u32) {
        self.tile_size = tile_size;
    }

    pub fn set_require_proper_pair(&mut self, require: bool) {
        self.require_proper_pair = require;
    }

    pub fn set_gc(&mut self, gc: ApplyGCArgs) {
        self.gc = gc;
    }

    pub fn set_ref_2bit(&mut self, ref_2bit: Option<PathBuf>) {
        self.ref_2bit = ref_2bit;
    }
}
