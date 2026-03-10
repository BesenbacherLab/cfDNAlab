use crate::commands::prepare_windows::labels::{MISSING_GROUP_LABEL, validate_label_token};
use anyhow::{Context, Result, bail};

/* Parse distance bins */

/// A single distance bin rule.
#[derive(Debug, Clone, Copy)]
pub enum DistanceExpr {
    Lt(i32),
    Le(i32),
    Gt(i32),
    Ge(i32),
    Range(i32, i32), // inclusive range [a,b]
}

#[derive(Debug, Clone)]
pub struct DistanceBin {
    label: String,
    expr: DistanceExpr,
}

/// Parsed distance-binning rules, checked in order; first match wins.
#[derive(Debug, Default)]
pub struct DistanceBins {
    rules: Vec<DistanceBin>,
}

impl DistanceBins {
    #[inline]
    pub fn match_label(&self, distance: i32) -> Option<&str> {
        for rule in &self.rules {
            let m = match rule.expr {
                DistanceExpr::Lt(x) => distance < x,
                DistanceExpr::Le(x) => distance <= x,
                DistanceExpr::Gt(x) => distance > x,
                DistanceExpr::Ge(x) => distance >= x,
                DistanceExpr::Range(a, b) => distance >= a && distance <= b,
            };
            if m {
                return Some(rule.label.as_str());
            }
        }
        None
    }
}

/// Parse distance bin specs like "prox:<500", "mid:500-2000", "dist:>2000".
///
/// Parameters
/// ----------
/// - specs:
///     List of bin specifications in the form "<label>:<expr>".
///
/// Returns
/// -------
/// - bins:
///     Parsed rules preserving input order.
pub fn parse_distance_bins(specs: &[String]) -> Result<DistanceBins> {
    let mut rules = Vec::with_capacity(specs.len());
    for spec in specs {
        let (label, expr) = spec
            .split_once(':')
            .with_context(|| format!("Invalid distance bin spec (missing ':'): '{}'", spec))?;

        let label = label.trim();
        if let Err(message) = validate_label_token(label, "distance bin label") {
            bail!(message);
        }

        let expr = expr.trim();
        let parsed = if let Some(num) = expr.strip_prefix("<=") {
            DistanceExpr::Le(num.parse::<i32>().context("Parsing <=N")?)
        } else if let Some(num) = expr.strip_prefix('<') {
            DistanceExpr::Lt(num.parse::<i32>().context("Parsing <N")?)
        } else if let Some(num) = expr.strip_prefix(">=") {
            DistanceExpr::Ge(num.parse::<i32>().context("Parsing >=N")?)
        } else if let Some(num) = expr.strip_prefix('>') {
            DistanceExpr::Gt(num.parse::<i32>().context("Parsing >N")?)
        } else if let Some((a, b)) = expr.split_once('-') {
            let a = a.trim().parse::<i32>().context("Parsing range A-B (A)")?;
            let b = b.trim().parse::<i32>().context("Parsing range A-B (B)")?;
            if b < a {
                bail!(
                    "Invalid distance range '{}': upper bound < lower bound",
                    expr
                );
            }
            DistanceExpr::Range(a, b)
        } else {
            bail!("Invalid distance expression '{}'", expr);
        };

        rules.push(DistanceBin {
            label: label.to_string(),
            expr: parsed,
        });
    }
    Ok(DistanceBins { rules })
}

/* Parse score filters */

/// Comparator used for score filter.
#[derive(Debug, Clone, Copy)]
pub enum CmpOp {
    Gt,
    Ge,
    Lt,
    Le,
    Eq,
    Ne,
}

#[derive(Debug, Clone)]
pub struct ScoreFilter {
    op: CmpOp,
    value: f32,
}

impl ScoreFilter {
    #[inline]
    pub fn eval(&self, score: f32) -> bool {
        match self.op {
            CmpOp::Gt => score > self.value,
            CmpOp::Ge => score >= self.value,
            CmpOp::Lt => score < self.value,
            CmpOp::Le => score <= self.value,
            CmpOp::Eq => score == self.value,
            CmpOp::Ne => score != self.value,
        }
    }
}

/// Parse a score filter string like ">=10" or "<0.05".
///
/// Parameters
/// ----------
/// - s:
///     The score filter expression.
///
/// Returns
/// -------
/// - filter:
///     Parsed filter ready to evaluate scores.
pub fn parse_score_filter(s: &str) -> Result<ScoreFilter> {
    let s = s.trim();
    let (op, rest) = if let Some(x) = s.strip_prefix(">=") {
        (CmpOp::Ge, x)
    } else if let Some(x) = s.strip_prefix("<=") {
        (CmpOp::Le, x)
    } else if let Some(x) = s.strip_prefix("==") {
        (CmpOp::Eq, x)
    } else if let Some(x) = s.strip_prefix("!=") {
        (CmpOp::Ne, x)
    } else if let Some(x) = s.strip_prefix('>') {
        (CmpOp::Gt, x)
    } else if let Some(x) = s.strip_prefix('<') {
        (CmpOp::Lt, x)
    } else {
        bail!("Invalid score filter '{}'", s);
    };
    let value = rest.trim().parse::<f32>().context("Parsing score value")?;
    Ok(ScoreFilter { op, value })
}

/* Parse bed */

/// Parse a single BED-like record from a line into (chrom, start, end, group?, score?).
///
/// This parser uses the config's `cols` mapping and `group_cols`.
///
/// Parameters
/// ----------
/// - line:
///     Input line (not including line terminator).
/// - separator:
///     Field separator character.
/// - cols_spec:
///     The `cols` mapping string (e.g., "chrom=0,start=1,end=2").
/// - group_cols:
///     Which columns to concatenate for group (may be empty).
/// - score_col:
///     Optional score column.
///
/// Returns
/// -------
/// - chrom:
///     Chromosome name.
/// - start:
///     Start coordinate.
/// - end:
///     End coordinate (exclusive).
/// - group:
///     Group string (empty if none).
/// - score:
///     Optional score.
pub struct ColumnIndices {
    chrom: usize,
    start: usize,
    end: usize,
    group: Vec<usize>,
    score: Option<usize>,
}

pub fn resolve_column_indices(
    cols_spec: &str,
    group_cols: &[String],
    score_col: Option<&str>,
) -> Result<ColumnIndices> {
    let (chrom, start, end) = parse_cols_indices(cols_spec)?;
    let mut group = Vec::with_capacity(group_cols.len());
    for spec in group_cols {
        group.push(parse_single_index(spec)?);
    }
    let score = if let Some(sc) = score_col {
        Some(parse_single_index(sc)?)
    } else {
        None
    };
    Ok(ColumnIndices {
        chrom,
        start,
        end,
        group,
        score,
    })
}

pub fn parse_record_line(
    line: &str,
    separator: char,
    cols: &ColumnIndices,
) -> Result<(String, u32, u32, String, Option<f32>)> {
    // Split all fields once
    let fields: Vec<&str> = line.split(separator).collect();

    let chrom = fields
        .get(cols.chrom)
        .context("Missing chrom field")?
        .trim()
        .to_string();

    let start: u32 = fields
        .get(cols.start)
        .context("Missing start field")?
        .trim()
        .parse()
        .context("Invalid start")?;

    let end: u32 = fields
        .get(cols.end)
        .context("Missing end field")?
        .trim()
        .parse()
        .context("Invalid end")?;

    if end <= start {
        bail!("End must be greater than start");
    }

    let group = if cols.group.is_empty() {
        String::new()
    } else {
        let mut parts: Vec<&str> = Vec::with_capacity(cols.group.len());
        for &idx in &cols.group {
            let val = fields.get(idx).unwrap_or(&"").trim();
            if val.is_empty() {
                parts.push(MISSING_GROUP_LABEL);
            } else {
                if let Err(message) = validate_label_token(val, "input group label") {
                    bail!(message);
                }
                parts.push(val);
            }
        }
        parts.join("__")
    };

    let score = if let Some(idx) = cols.score {
        let val = fields.get(idx).unwrap_or(&"").trim();
        if val.is_empty() {
            None
        } else {
            Some(val.parse::<f32>().context("Invalid score")?)
        }
    } else {
        None
    };

    Ok((chrom, start, end, group, score))
}

pub fn parse_cols_indices(cols_spec: &str) -> Result<(usize, usize, usize)> {
    let mut chrom_idx = None;
    let mut start_idx = None;
    let mut end_idx = None;
    for part in cols_spec.split(',') {
        let (key, val) = part
            .split_once('=')
            .with_context(|| format!("Invalid cols spec '{}'", part))?;
        let idx = parse_single_index(val.trim())?;
        match key.trim() {
            "chrom" => chrom_idx = Some(idx),
            "start" => start_idx = Some(idx),
            "end" => end_idx = Some(idx),
            other => bail!("Unknown cols key '{}'", other),
        }
    }
    Ok((
        chrom_idx.context("cols: missing chrom=")?,
        start_idx.context("cols: missing start=")?,
        end_idx.context("cols: missing end=")?,
    ))
}

#[inline]
pub(crate) fn parse_single_index(s: &str) -> Result<usize> {
    let idx = s
        .parse::<usize>()
        .with_context(|| format!("Expecting 0-based index, got '{}'", s))?;
    Ok(idx)
}
