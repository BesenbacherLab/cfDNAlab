use crate::shared::constants::{DEFAULT_MAX_SOFT_CLIPS, MAX_MAX_SOFT_CLIPS};
use std::str::FromStr;

#[cfg(feature = "cli")]
use clap::ValueEnum;

const BASE_QUALITY_FILTER_USAGE: &str = concat!(
    "Use '<agg> in <scope> <op> <threshold>' with ",
    "<agg> in {'min', 'mean', 'max'}, ",
    "<scope> in {'end', 'fragment'}, ",
    "and <op> in {'>=', '>', '<=', '<'}"
);

/// Aggregate base qualities across the inside bases of one counted unit.
///
/// These reductions each map to a clear thresholding question:
/// weakest point (`min`), overall level (`mean`), or best point (`max`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseQualityAggregation {
    Min,
    Mean,
    Max,
}

impl BaseQualityAggregation {
    #[inline]
    pub fn as_str(self) -> &'static str {
        match self {
            BaseQualityAggregation::Min => "min",
            BaseQualityAggregation::Mean => "mean",
            BaseQualityAggregation::Max => "max",
        }
    }
}

impl FromStr for BaseQualityAggregation {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "min" => Ok(BaseQualityAggregation::Min),
            "mean" => Ok(BaseQualityAggregation::Mean),
            "max" => Ok(BaseQualityAggregation::Max),
            _ => Err(format!(
                "Invalid base-quality aggregation '{s}'. {BASE_QUALITY_FILTER_USAGE}"
            )),
        }
    }
}

/// Decide whether a base-quality filter applies to one end or to the full fragment.
///
/// Fragment-level filters are intended for cases where the two fragment ends
/// should first be summarized into one score before thresholding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseQualityFilterScope {
    End,
    Fragment,
}

impl BaseQualityFilterScope {
    #[inline]
    pub fn as_str(self) -> &'static str {
        match self {
            BaseQualityFilterScope::End => "end",
            BaseQualityFilterScope::Fragment => "fragment",
        }
    }
}

impl FromStr for BaseQualityFilterScope {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "end" => Ok(BaseQualityFilterScope::End),
            "fragment" => Ok(BaseQualityFilterScope::Fragment),
            _ => Err(format!(
                "Invalid base-quality filter scope '{s}'. {BASE_QUALITY_FILTER_USAGE}"
            )),
        }
    }
}

/// Supported comparison operators for parsed base-quality filters.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseQualityComparisonOp {
    Gt,
    Ge,
    Lt,
    Le,
}

impl BaseQualityComparisonOp {
    #[inline]
    pub fn as_str(self) -> &'static str {
        match self {
            BaseQualityComparisonOp::Gt => ">",
            BaseQualityComparisonOp::Ge => ">=",
            BaseQualityComparisonOp::Lt => "<",
            BaseQualityComparisonOp::Le => "<=",
        }
    }
}

impl BaseQualityComparisonOp {
    #[inline]
    pub fn eval(self, value: f32, threshold: f32) -> bool {
        match self {
            BaseQualityComparisonOp::Gt => value > threshold,
            BaseQualityComparisonOp::Ge => value >= threshold,
            BaseQualityComparisonOp::Lt => value < threshold,
            BaseQualityComparisonOp::Le => value <= threshold,
        }
    }
}

/// One parsed `--bq-filter` expression.
///
/// The grammar is:
///
/// - `<agg> in <scope> <op> <threshold>`
///
/// Example expressions:
///
/// - `min in end >= 30`
///
/// - `mean in fragment < 25`
///
/// - `max in fragment < 20`
///
/// Repeating `--bq-filter` counts only ends that pass all end filters and belong to
/// fragments that pass all fragment filters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BaseQualityFilter {
    pub aggregation: BaseQualityAggregation,
    pub scope: BaseQualityFilterScope,
    pub op: BaseQualityComparisonOp,
    pub threshold: f32,
}

impl BaseQualityFilter {
    #[inline]
    pub fn as_cli_expr(self) -> String {
        format!(
            "{} in {} {} {}",
            self.aggregation.as_str(),
            self.scope.as_str(),
            self.op.as_str(),
            self.threshold
        )
    }

    #[inline]
    pub fn passes_value(self, value: f32) -> bool {
        self.op.eval(value, self.threshold)
    }
}

impl FromStr for BaseQualityFilter {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        let tokens: Vec<&str> = s.split_whitespace().collect();
        if tokens.len() < 4 {
            return Err(format!(
                "Invalid base-quality filter '{s}'. {BASE_QUALITY_FILTER_USAGE}"
            ));
        }
        if tokens.len() > 5 {
            return Err(format!(
                "Invalid base-quality filter '{s}'. {BASE_QUALITY_FILTER_USAGE}"
            ));
        }

        let aggregation = tokens[0].parse::<BaseQualityAggregation>()?;
        if !tokens[1].eq_ignore_ascii_case("in") {
            return Err(format!(
                "Invalid base-quality filter '{s}'. {BASE_QUALITY_FILTER_USAGE}"
            ));
        }
        let scope = tokens[2].parse::<BaseQualityFilterScope>()?;
        let comparison_expr = tokens[3..].join("");
        let (op, threshold_str) = parse_base_quality_comparison(&comparison_expr)
            .map_err(|msg| format!("Invalid base-quality filter '{s}'. {msg}"))?;
        let threshold = threshold_str
            .parse::<f32>()
            .map_err(|_| format!("Invalid base-quality threshold '{threshold_str}'"))?;
        if !threshold.is_finite() {
            return Err(format!(
                "Base-quality threshold must be finite, got '{threshold_str}'"
            ));
        }
        if threshold < 0.0 {
            return Err(format!(
                "Base-quality threshold must be >= 0, got '{threshold_str}'"
            ));
        }

        Ok(BaseQualityFilter {
            aggregation,
            scope,
            op,
            threshold,
        })
    }
}

fn parse_base_quality_comparison(
    s: &str,
) -> std::result::Result<(BaseQualityComparisonOp, &str), String> {
    if let Some(rest) = s.strip_prefix(">=") {
        Ok((BaseQualityComparisonOp::Ge, rest))
    } else if let Some(rest) = s.strip_prefix("<=") {
        Ok((BaseQualityComparisonOp::Le, rest))
    } else if let Some(rest) = s.strip_prefix('>') {
        Ok((BaseQualityComparisonOp::Gt, rest))
    } else if let Some(rest) = s.strip_prefix('<') {
        Ok((BaseQualityComparisonOp::Lt, rest))
    } else {
        Err(BASE_QUALITY_FILTER_USAGE.to_string())
    }
}

// TODO these structs should use the format used by other cfDNAlab commands instead

/// Select where the inside-fragment bases come from for end motifs.
///
/// The `ends` pipeline can either trust the read sequence itself or reconstruct
/// the inside-fragment half from the reference genome. The read-backed mode is
/// the default because it reflects what was actually sequenced, while the
/// reference-backed mode is useful when you want alignment-consistent sequence
/// context and are willing to skip indel-affected motifs.
#[cfg_attr(feature = "cli", derive(ValueEnum))]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum KmerSource {
    /// Extract the inside-fragment bases from the sequenced read.
    Read,
    /// Extract the inside-fragment bases from the reference genome.
    Reference,
}

/// Select how clipped fragment ends should be interpreted.
///
/// End motifs can either follow the aligned fragment span, include clipped read
/// bases at aligned genomic boundaries, include clipped read bases at shifted
/// genomic boundaries, or be skipped when clipping is present.
#[cfg_attr(feature = "cli", derive(ValueEnum))]
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub enum ClipStrategy {
    /// Use the aligned fragment ends.
    Aligned,

    /// Include clipped read bases, but keep the aligned genomic boundary.
    IncludeAtAlignedBoundary,

    /// Include clipped read bases and shift the genomic boundary outward.
    IncludeAtShiftedBoundary,

    /// Skip motifs whose end is soft-clipped.
    #[default]
    Skip,
}

impl ClipStrategy {
    #[inline]
    pub fn includes_clipped_inside_bases(self) -> bool {
        matches!(
            self,
            ClipStrategy::IncludeAtAlignedBoundary | ClipStrategy::IncludeAtShiftedBoundary
        )
    }

    #[inline]
    pub fn uses_shifted_boundary(self) -> bool {
        matches!(self, ClipStrategy::IncludeAtShiftedBoundary)
    }
}

#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct ClippingArgs {
    /// How to extract a motif when its fragment end is clipped `[string]`
    ///
    /// Clipping means the read contains terminal bases that the aligner did not align normally.
    /// The choice here is thus what positions to count when that happens.
    ///
    /// For extraction of **outside** bases, we suggest **skipping** fragments
    /// with soft clipping, as it is difficult to infer where on the
    /// reference genome the actual fragment end was. We do provide two
    /// include-at-boundary modes for this, but neither is perfect.
    ///
    /// **NOTE**: Fragments with **hard**-clipping are always discarded.
    ///
    /// Possible values:
    ///
    /// - `"skip"`:
    ///   Skip motifs when their fragment end is soft-clipped.
    ///
    /// - `"aligned"`:
    ///   Use the aligned start and end positions (the usual `cfDNAlab` fragment definition).
    ///   This ignores clipped bases in the read sequences.
    ///   
    ///   **NOTE**: If the aligner clipped the actual DNA molecule, these motifs may not reflect
    ///   the actual fragment ends.
    ///
    /// - `"include-at-aligned-boundary"`:
    ///   Include soft-clipped read bases, but keep the
    ///   **aligned** fragment-end genomic boundary for outside-base lookup,
    ///   window assignment, and motif-level blacklist validation.
    ///
    ///   This setting is only supported with `--source-inside read`.
    ///
    /// - `"include-at-shifted-boundary"`:
    ///   Include soft-clipped read bases, and **move** the
    ///   fragment-end boundary outside the aligned span by the clipped
    ///   length.
    ///
    ///   This shifted boundary is used for outside-base lookup,
    ///   window assignment, and blacklist filtering.
    ///
    ///   File-based GC correction and scaling-factor weighting still use the
    ///   aligned reference span. If the aligned length falls outside the GC
    ///   package range, the fragment is considered invalid and is included
    ///   in the GC correction failure statistics. When `--assign-by count-overlap`,
    ///   clipped-only window contributions use the nearest aligned reference base for scaling.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_enum,
            hide_possible_values = true,
            default_value = "skip",
            help_heading = "Clipping"
        )
    )]
    pub clip_strategy: ClipStrategy,

    /// Skip motifs whose relevant end has more soft-clipped bases than this `[integer]`
    ///
    /// This limit is applied independently to each fragment end.
    ///
    /// Fragment length filtering is applied after soft clip expansion.
    ///
    /// Use `--clip-strategy skip` to discard all soft-clipped motifs.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value_t = DEFAULT_MAX_SOFT_CLIPS,
            value_parser = clap::value_parser!(u16).range(0..=MAX_MAX_SOFT_CLIPS as i64),
            help_heading = "Clipping"
        )
    )]
    pub max_soft_clips: u16,
}

#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone, Default, PartialEq)]
pub struct AssignMotifToWindowArgs {
    /// When to assign motifs to windows `[string]`
    ///
    /// The default `"endpoint"` option assigns each motif separately by its own fragment-end position.
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
    /// Midpoints for even-sized fragments use a deterministic coordinate-derived random seed to
    /// select either the left or right base. Duplicate fragments with the same coordinates get the
    /// same choice. This avoids fixed rounding bias while keeping repeated runs reproducible.
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

#[cfg(test)]
mod tests {
    include!("config_structs_tests.rs");
}
