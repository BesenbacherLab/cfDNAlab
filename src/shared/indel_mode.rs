use std::str::FromStr;

/// Indel-handling strategy.
///
/// Possible values:
///     "ignore", "adjust", or "skip" [string]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum IndelMode {
    /// Ignore indels.
    #[default]
    Ignore,
    /// Adjust for observed insertions and deletions
    /// in the aligned bases. In read-overlap positions,
    /// the reads must agree about the presence of an indel.
    Adjust,
    /// Skip fragments with any insertions or deletions.
    Skip,
}

impl FromStr for IndelMode {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s == "ignore" {
            Ok(IndelMode::Ignore)
        } else if s == "adjust" {
            Ok(IndelMode::Adjust)
        } else if s == "skip" {
            Ok(IndelMode::Skip)
        } else {
            Err("Use 'ignore', 'adjust', or 'skip'".into())
        }
    }
}

/// Policy for when to filter motifs due to indels.
///
/// Possible values:
///     "auto", "skip-affected-end", or "skip-affected-fragment" [string]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum IndelMotifFilterPolicy {
    /// Select the option based on the source.
    ///
    /// - For read-sequence bases, allow indels in the alignment.
    ///
    /// - For reference bases, skip motifs with indels in the alignment.
    #[default]
    Auto,
    /// Always skip motifs overlapping indels.
    SkipAffectedEnd,
    /// Skip **fragments** when either end overlap indels.
    SkipAffectedFragment,
}

impl FromStr for IndelMotifFilterPolicy {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s == "auto" {
            Ok(IndelMotifFilterPolicy::Auto)
        } else if s == "skip-affected-end" {
            Ok(IndelMotifFilterPolicy::SkipAffectedEnd)
        } else if s == "skip-affected-fragment" {
            Ok(IndelMotifFilterPolicy::SkipAffectedFragment)
        } else {
            Err("Use 'auto', 'skip-affected-end', or 'skip-affected-fragment'".into())
        }
    }
}
