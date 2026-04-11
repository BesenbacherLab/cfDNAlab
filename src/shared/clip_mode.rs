use std::str::FromStr;

/// Soft clipping strategy.
///
/// Possible values:
///     "ignore", "adjust", or "skip" [string]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum ClipMode {
    /// Ignore soft clipped fragments.
    #[default]
    Ignore,
    /// Adjust for observed soft clipping in the fragment ends.
    Adjust,
    /// Skip fragments with any clipping.
    Skip,
}

impl FromStr for ClipMode {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s == "ignore" {
            Ok(ClipMode::Ignore)
        } else if s == "adjust" {
            Ok(ClipMode::Adjust)
        } else if s == "skip" {
            Ok(ClipMode::Skip)
        } else {
            Err("Use 'ignore', 'adjust', or 'skip'".into())
        }
    }
}
