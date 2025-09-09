use std::str::FromStr;

/// Blacklist strategy for fragment/read/interval filtering
#[derive(Debug, Clone)]
pub enum BlackStrategy {
    Full,
    Any,
    Midpoint,
    Proportion(f64),
}

impl FromStr for BlackStrategy {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s == "full" {
            Ok(BlackStrategy::Full)
        } else if s == "any" {
            Ok(BlackStrategy::Any)
        } else if s == "midpoint" {
            Ok(BlackStrategy::Midpoint)
        } else if let Some(v) = s.strip_prefix("proportion=") {
            let thr: f64 = v
                .parse()
                .map_err(|e: std::num::ParseFloatError| e.to_string())?;
            if !(0.0..=1.0).contains(&thr) {
                Err("Proportion must be between 0.0 and 1.0".into())
            } else {
                Ok(BlackStrategy::Proportion(thr))
            }
        } else {
            Err("Use 'full', 'midpoint', or 'proportion=<0.0–1.0>'".into())
        }
    }
}
