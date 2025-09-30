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
            Err("Use 'ignore', 'adjus', or 'skip'".into())
        }
    }
}
