use std::str::FromStr;

/// How to write blacklisted windows / positions
/// (represented as NaN).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum NanPolicy {
    /// Skip rows where cov.is_nan()
    #[default]
    DropRow,
    /// Write the literal string "NaN"
    WriteLiteralNaN,
    /// Leave the field empty
    WriteEmptyCell,
}

// For the CLI
impl FromStr for NanPolicy {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s == "drop" {
            Ok(NanPolicy::DropRow)
        } else if s == "nan" {
            Ok(NanPolicy::WriteLiteralNaN)
        } else if s == "empty" {
            Ok(NanPolicy::WriteEmptyCell)
        } else {
            Err("Use 'drop', 'nan', or 'empty'".into())
        }
    }
}

impl NanPolicy {
    #[inline]
    pub fn drop_row(&self) -> bool {
        matches!(self, NanPolicy::DropRow)
    }

    /// Render the coverage cell for a blacklisted site
    /// - DropRow -> None (caller should skip the row)
    /// - WriteLiteralNaN -> Some("NaN")
    /// - WriteEmptyCell -> Some("")
    #[inline]
    pub fn render_masked_cell(&self) -> Option<&'static str> {
        match self {
            NanPolicy::DropRow => None,
            NanPolicy::WriteLiteralNaN => Some("NaN"),
            NanPolicy::WriteEmptyCell => Some(""),
        }
    }
}
