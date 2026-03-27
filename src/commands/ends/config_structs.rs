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

// TODO: for align-start, ensure 2xk_within is left after shifting - so (fragment length - total clipping) > 2xk_within
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Copy, Clone)]
pub struct ClippingArgs {
    /// How to extract a motif when its fragment end is clipped `[string]`
    ///
    /// Clipping means the read contains terminal bases that the aligner did not align normally.
    /// The choice here is thus what sequence object to count when that happens.
    ///
    /// We recommend either `"raw"` or `"drop"` unless you have a good reason for the other options.
    ///
    /// Possible values:
    ///
    /// - `"raw"`:
    ///   Use the raw read bases, clipped or not.
    ///
    ///   This treats the sequenced end as the object of interest.
    ///
    ///   Most faithful to the observed molecule end, but clipped technical sequence is kept.
    ///
    /// - `"drop"`:
    ///   Skip motifs where the end is clipped.
    ///
    ///   Use this when you only want ends the aligner fully trusted.
    ///
    ///   Conservative and simple, but discards potentially real biology at clipped ends.
    ///
    /// - `"align-start"`:
    ///   Skip clipped bases and take the first aligned bases instead.
    ///
    ///   This is useful when clipping mostly reflects adapters or technical junk,
    ///   but it changes the object from the terminal read bases to an interior aligned motif.
    ///
    ///   Can recover motifs when clipping is mostly technical, but the motif is shifted inward.
    ///
    /// - `"fill-with-ref"`:
    ///   Use read bases where aligned, and fill the clipped part from the reference.
    ///
    ///   This gives a hybrid observed/imputed motif. It assumes the clipped terminal bases
    ///   are better approximated by the reference than by the clipped read sequence.
    ///
    ///   Retains more motifs under clipping, but mixes observed read bases with imputed reference bases.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            default_value = "raw",
            requires_if("fill-with-ref", "ref_2bit"),
            help_heading = "Clipping"
        )
    )]
    pub clip_strategy: ClipStrategy,

    /// Skip motifs with a higher number of (S+H) clipped bases than this.
    ///
    /// Use `--clip-strategy drop` to discard all clipped motifs.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Clipping"))]
    pub max_clips: Option<usize>,
}

#[cfg_attr(feature = "cli", derive(clap::ValueEnum))]
#[derive(Copy, Clone, Debug)]
pub enum ClipStrategy {
    /// Use the raw read bases, clipped or not.
    Raw,

    /// Drop the motif if its end is clipped.
    Drop,

    /// Skip clipped bases and take the first aligned bases instead.
    AlignStart,

    /// Use read bases where aligned, and fill the clipped part from the reference.
    FillWithRef,
}

#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone, Default)]
pub struct AssignMotifToWindowArgs {
    /// When to assign motifs to windows `[string]`
    ///
    /// The default `"endpoint"` option assigns each motif by its own fragment-end position.
    ///
    /// The other modes ask which windows the **fragment** contributes to,
    /// and the fragment's motif(s) are then counted in those window(s).
    ///
    /// Possible values:
    ///     `"endpoint"`, `"count-overlap"`, `"any"`, `"all"`, `"midpoint"`, or `"proportion=<threshold>"`
    ///
    /// `"endpoint"`: Count each motif in the windows overlapping its fragment-end position.
    /// The two fragment ends may be counted in separate windows.
    ///
    /// `"count-overlap"`: Count up the fraction of fragment bases overlapping each window.
    ///
    /// `"any"`, `"all"`, or `"proportion=<threshold>"`:
    /// Assign motifs when a proportion of fragment bases overlap a window.
    ///
    /// Example of proportion: `--assign-by proportion=0.2` (no space around `=`)
    ///
    /// `"midpoint"`: Assign motifs when the fragment midpoint overlaps a window.
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
            help = "What to assign motifs to windows by (or count motifs as).",
            help_heading = "Window Assignment"
        )
    )]
    pub assign_by: WindowMotifAssigner,
}

// TODO: In the future we might want to add window-based overlap variants (WindowProportion etc.). Not relevant yet.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
/// How to assign a fragment to windows.
///
/// "Endpoint" is motif-end specific, the others are fragment-centric.
///
/// NOTE: This only considers the proportion of **fragment positions**
/// overlapping the window. For window sizes smaller than fragments
/// this means a fragment could overlap a window fully but
/// have < 100% of fragment positions inside the window.
pub enum WindowMotifAssigner {
    /// Assign the **motif** to windows overlapping its specific fragment end.
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
