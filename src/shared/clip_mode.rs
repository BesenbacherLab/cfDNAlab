use std::str::FromStr;

/// Soft clipping strategy.
///
/// Possible values:
///     "aligned", "adjust", or "skip" [string]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum ClipMode {
    /// Use the aligned fragment boundaries and ignore clipped bases in the length.
    #[default]
    Aligned,
    /// Adjust for observed soft clipping in the fragment ends.
    Adjust,
    /// Skip fragments with any clipping.
    Skip,
}

impl FromStr for ClipMode {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s == "aligned" {
            Ok(ClipMode::Aligned)
        } else if s == "adjust" {
            Ok(ClipMode::Adjust)
        } else if s == "skip" {
            Ok(ClipMode::Skip)
        } else {
            Err("Use 'aligned', 'adjust', or 'skip'".into())
        }
    }
}
