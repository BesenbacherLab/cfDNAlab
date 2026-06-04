use std::str::FromStr;

/// What to do per window
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PeaksWindowAction {
    #[default]
    Stats,
    OnlyIncludeThesePositionsUnique,
    OnlyIncludeThesePositionsIndexed,
}

// For the CLI
impl FromStr for PeaksWindowAction {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s == "stats" {
            Ok(PeaksWindowAction::Stats)
        } else if s == "unique-positions" {
            Ok(PeaksWindowAction::OnlyIncludeThesePositionsUnique)
        } else if s == "indexed-positions" {
            Ok(PeaksWindowAction::OnlyIncludeThesePositionsIndexed)
        } else {
            Err("Use 'stats', 'indexed-positions', or 'unique-positions'".into())
        }
    }
}
