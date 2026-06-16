use std::str::FromStr;

/// Blacklist strategy for fragment/read/interval filtering
///
/// Possible values:
///     `"any"`, `"all"`, `"midpoint"`, or `"proportion=<threshold>"` [string]
///
/// Example of proportion: `--blacklist-strategy proportion=0.2` (no space around `=`)
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum BlacklistStrategy {
    /// Any positions overlap with blacklisted regions.
    #[default]
    Any,
    /// All positions overlap with blacklisted regions.
    All,
    /// The central base overlaps with blacklisted regions.
    ///
    /// Even-length fragments have two central bases in discrete base coordinates.
    /// The fragment is blacklisted if either central base overlaps a blacklist interval.
    Midpoint,
    /// A given proportion of positions overlap with blacklisted regions (e.g. `proportion=0.2`).
    Proportion(f64),
}

impl FromStr for BlacklistStrategy {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s == "all" {
            Ok(BlacklistStrategy::All)
        } else if s == "any" {
            Ok(BlacklistStrategy::Any)
        } else if s == "midpoint" {
            Ok(BlacklistStrategy::Midpoint)
        } else if let Some(v) = s.strip_prefix("proportion=") {
            let thr: f64 = v
                .parse()
                .map_err(|e: std::num::ParseFloatError| e.to_string())?;
            if !(0.0..=1.0).contains(&thr) {
                Err("Proportion must be between 0.0 and 1.0".into())
            } else {
                Ok(BlacklistStrategy::Proportion(thr))
            }
        } else {
            Err("Use 'any', 'all', 'midpoint', or 'proportion=<0.0–1.0>'".into())
        }
    }
}
