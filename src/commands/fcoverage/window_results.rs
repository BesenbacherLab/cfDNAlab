use std::str::FromStr;

/// What to do per window
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CoverageWindowAction {
    #[default]
    Average,
    Total,
    SummaryStats,
    AverageOnUniqueBases,
    TotalOnUniqueBases,
    SummaryStatsOnUniqueBases,
    OnlyIncludeThesePositionsUnique,
    OnlyIncludeThesePositionsIndexed,
}

impl CoverageWindowAction {
    pub fn is_positional(self) -> bool {
        matches!(
            self,
            Self::OnlyIncludeThesePositionsUnique | Self::OnlyIncludeThesePositionsIndexed
        )
    }

    pub fn is_summary_stats(self) -> bool {
        matches!(self, Self::SummaryStats | Self::SummaryStatsOnUniqueBases)
    }

    pub fn is_unique_base_grouped_action(self) -> bool {
        matches!(
            self,
            Self::AverageOnUniqueBases | Self::TotalOnUniqueBases | Self::SummaryStatsOnUniqueBases
        )
    }

    pub fn action_file_stem(self) -> &'static str {
        match self {
            Self::Average => "average",
            Self::Total => "total",
            Self::SummaryStats => "summary_stats",
            Self::AverageOnUniqueBases => "average_on_unique_bases",
            Self::TotalOnUniqueBases => "total_on_unique_bases",
            Self::SummaryStatsOnUniqueBases => "summary_stats_on_unique_bases",
            Self::OnlyIncludeThesePositionsUnique => "per_position",
            Self::OnlyIncludeThesePositionsIndexed => "per_position_per_window",
        }
    }
}

// For the CLI
impl FromStr for CoverageWindowAction {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s == "average" {
            Ok(CoverageWindowAction::Average)
        } else if s == "total" {
            Ok(CoverageWindowAction::Total)
        } else if s == "summary-stats" {
            Ok(CoverageWindowAction::SummaryStats)
        } else if s == "average-on-unique-bases" {
            Ok(CoverageWindowAction::AverageOnUniqueBases)
        } else if s == "total-on-unique-bases" {
            Ok(CoverageWindowAction::TotalOnUniqueBases)
        } else if s == "summary-stats-on-unique-bases" {
            Ok(CoverageWindowAction::SummaryStatsOnUniqueBases)
        } else if s == "unique-positions" {
            Ok(CoverageWindowAction::OnlyIncludeThesePositionsUnique)
        } else if s == "indexed-positions" {
            Ok(CoverageWindowAction::OnlyIncludeThesePositionsIndexed)
        } else {
            Err(
                "Use 'average', 'total', 'summary-stats', 'average-on-unique-bases', \
'total-on-unique-bases', 'summary-stats-on-unique-bases', 'indexed-positions', or \
'unique-positions'"
                    .into(),
            )
        }
    }
}
