use std::str::FromStr;

/// Length calculation mode.
///
/// Possible values:
///     "reference", "indel-adjusted", or "skip-indels" [string]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum LengthMode {
    /// Use reference length (`end(reverse) - start(forward)`)
    #[default]
    Reference,
    /// Adjust the reference length for observed insertions and deletions
    /// in the aligned bases. In a read-overlap, the reads must
    /// agree, with the minimum insertion selected.
    IndelAdjusted,
    /// Skip fragments with any insertions or deletions.
    SkipIndels,
}

impl FromStr for LengthMode {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s == "reference" {
            Ok(LengthMode::Reference)
        } else if s == "indel-adjusted" {
            Ok(LengthMode::IndelAdjusted)
        } else if s == "skip-indels" {
            Ok(LengthMode::SkipIndels)
        } else {
            Err("Use 'reference', 'indel-adjusted', or 'skip-indels'".into())
        }
    }
}
