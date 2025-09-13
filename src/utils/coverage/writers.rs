use crate::utils::coverage::window_results::{
    CoverageOutput, CoverageWindowAction, WindowResult, WindowValue,
};
use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;
use std::cmp::Ordering;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::path::Path;
use std::path::PathBuf;
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

/// Write per-window aggregates (Average or Total) to TSV
///
/// Parameters
/// ----------
/// - path
///     Output file path
/// - by_chr
///     Map from chromosome -> `CoverageOutput`
/// - action
///     Must be `CoverageWindowAction::Average` or `CoverageWindowAction::Total`
///
/// Output
/// ------
/// - Columns: `chromosome  start  end  <avg|total>  <blacklisted_positions|empty>`
pub fn write_aggregates_tsv<P: AsRef<std::path::Path>>(
    path: P,
    by_chr: &FxHashMap<String, CoverageOutput>,
    action: CoverageWindowAction,
) -> anyhow::Result<()> {
    use anyhow::{Context, bail};
    use std::cmp::Ordering;
    use std::fs::File;
    use std::io::{BufWriter, Write};

    if !matches!(
        action,
        CoverageWindowAction::Average | CoverageWindowAction::Total
    ) {
        bail!("write_aggregates_tsv requires action = Average or Total");
    }

    #[derive(Debug)]
    struct Row<'a> {
        chr: &'a str,
        start: u64,
        end: u64,
        value_f64: f64,
        num_blacklisted_pos: Option<u32>, // <- only written if present anywhere
    }
    let mut rows: Vec<Row> = Vec::new();

    for (chr, out) in by_chr {
        let results: &[WindowResult] = match out {
            CoverageOutput::PerWindow {
                action: act,
                results,
            } => {
                if *act != action {
                    bail!(
                        "chrom {} has PerWindow action {:?}, expected {:?}",
                        chr,
                        act,
                        action
                    );
                }
                results
            }
            CoverageOutput::WholePositional { .. } => {
                bail!(
                    "chrom {} has WholePositional output which is positional, not aggregate",
                    chr
                )
            }
        };

        for wr in results {
            match wr.value {
                WindowValue::Average(v) => rows.push(Row {
                    chr,
                    start: wr.start,
                    end: wr.end,
                    value_f64: v as f64,
                    num_blacklisted_pos: wr.num_blacklisted_pos, // use as-is; no fallback computation
                }),
                WindowValue::Total(v) => rows.push(Row {
                    chr,
                    start: wr.start,
                    end: wr.end,
                    value_f64: v,
                    num_blacklisted_pos: wr.num_blacklisted_pos, // use as-is; no fallback computation
                }),
                WindowValue::Positions(_) => bail!(
                    "chrom {} has positional WindowValue::Positions in an aggregate TSV",
                    chr
                ),
            }
        }
    }

    // Sort rows by chr, start, end
    rows.sort_by(|a, b| match a.chr.cmp(b.chr) {
        Ordering::Equal => match a.start.cmp(&b.start) {
            Ordering::Equal => a.end.cmp(&b.end),
            o => o,
        },
        o => o,
    });

    let any_bl = rows.iter().any(|r| r.num_blacklisted_pos.is_some());
    let mut w = BufWriter::new(
        File::create(&path).with_context(|| format!("Creating {:?}", path.as_ref()))?,
    );

    let value_col = match action {
        CoverageWindowAction::Average => "avg_coverage",
        CoverageWindowAction::Total => "total_coverage",
        _ => unreachable!(),
    };

    if any_bl {
        writeln!(
            w,
            "chromosome\tstart\tend\t{}\tblacklisted_positions",
            value_col
        )
        .context("Write header failed")?;
    } else {
        writeln!(w, "chromosome\tstart\tend\t{}", value_col).context("Write header failed")?;
    }

    for r in rows {
        if any_bl {
            if let Some(n) = r.num_blacklisted_pos {
                writeln!(
                    w,
                    "{}\t{}\t{}\t{}\t{}",
                    r.chr, r.start, r.end, r.value_f64, n
                )
                .context("Write data line failed")?;
            } else {
                // leave blank when not provided
                writeln!(w, "{}\t{}\t{}\t{}\t", r.chr, r.start, r.end, r.value_f64)
                    .context("Write data line failed")?;
            }
        } else {
            writeln!(w, "{}\t{}\t{}\t{}", r.chr, r.start, r.end, r.value_f64)
                .context("Write data line failed")?;
        }
    }
    Ok(())
}

/// Write positional TSV
///
/// Behavior depends on the payload:
/// - Whole-genome positional -> columns `chromosome  start  end  cov`
/// - Windowed positional     -> columns `chromosome  start  end  cov  orig_idx`
///
/// Parameters
/// ----------
/// - path:
///     Output file path
/// - by_chr:
///     Map from chromosome -> `CoverageOutput`
/// - policy:
///     What to write when coverage values are NaN.
///
/// Notes
/// -----
/// - All rows are collected and sorted by `chr, start, end, orig_idx`
/// - If you mix WholePositional and PerWindow::Positions across chromosomes, this will error
pub fn write_positions_tsv<P: AsRef<Path>>(
    path: P,
    by_chr: &FxHashMap<String, CoverageOutput>,
    nan_policy: NanPolicy,
) -> Result<()> {
    // Decide the mode across all chromosomes
    #[derive(Copy, Clone, Debug, Eq, PartialEq)]
    enum Mode {
        Whole,
        Windowed,
    }
    let mut mode: Option<Mode> = None;

    for out in by_chr.values() {
        let m = match out {
            CoverageOutput::WholePositional { .. } => Mode::Whole,
            CoverageOutput::PerWindow { results, .. } => {
                // Ensure all WindowValues are Positions
                for wr in results {
                    if !matches!(wr.value, WindowValue::Positions(_)) {
                        bail!(
                            "write_positions_tsv expects WindowValue::Positions, found {:?}",
                            wr.value
                        );
                    }
                }
                Mode::Windowed
            }
        };
        match mode {
            None => mode = Some(m),
            Some(prev) if prev != m => {
                bail!("Cannot mix WholePositional and windowed positional outputs in one file")
            }
            _ => {}
        }
    }
    let mode = mode.ok_or_else(|| anyhow::anyhow!("Empty output set"))?;

    // Collect rows
    #[derive(Debug)]
    struct Row<'a> {
        chr: &'a str,
        start: u64,
        end: u64,
        cov: f32,
        orig_idx: Option<u64>,
    }
    let mut rows: Vec<Row> = Vec::new();

    match mode {
        Mode::Whole => {
            for (chr, out) in by_chr {
                let (start, end, values) = match out {
                    CoverageOutput::WholePositional { start, end, values } => {
                        (*start, *end, values)
                    }
                    _ => unreachable!(),
                };
                // Emit one row per position [i, i+1)
                let mut pos = start;
                for &c in values {
                    rows.push(Row {
                        chr,
                        start: pos,
                        end: pos + 1,
                        cov: c,
                        orig_idx: None,
                    });
                    pos += 1;
                }
                // Sanity: pos should equal end
                debug_assert_eq!(pos, end);
            }
        }
        Mode::Windowed => {
            for (chr, out) in by_chr {
                let (action, results) = match out {
                    CoverageOutput::PerWindow { action, results } => (action, results),
                    _ => unreachable!(),
                };
                // Ensure action is the positional one
                if *action != CoverageWindowAction::OnlyIncludeThesePositions {
                    bail!(
                        "write_positions_tsv expects OnlyIncludeThesePositions, got {:?}",
                        action
                    );
                }
                for wr in results {
                    let WindowValue::Positions(vals) = &wr.value else {
                        unreachable!()
                    };
                    // Map each value to its genomic position
                    for (i, &c) in vals.iter().enumerate() {
                        let s = wr.start + i as u64;
                        rows.push(Row {
                            chr,
                            start: s,
                            end: s + 1,
                            cov: c,
                            orig_idx: Some(wr.original_idx),
                        });
                    }
                }
            }
        }
    }

    // Sort rows by chr, start, end, orig_idx
    rows.sort_by(|a, b| match a.chr.cmp(b.chr) {
        Ordering::Equal => match a.start.cmp(&b.start) {
            Ordering::Equal => match a.end.cmp(&b.end) {
                Ordering::Equal => a.orig_idx.cmp(&b.orig_idx),
                o => o,
            },
            o => o,
        },
        o => o,
    });

    // Write header
    let mut w = BufWriter::new(
        File::create(&path).with_context(|| format!("Creating {:?}", path.as_ref()))?,
    );

    match mode {
        Mode::Whole => {
            writeln!(w, "chromosome\tstart\tend\tcov").context("Write header failed")?;
        }
        Mode::Windowed => {
            writeln!(w, "chromosome\tstart\tend\tcov\torig_idx").context("Write header failed")?;
        }
    }

    // Emit rows
    for r in rows {
        // Skip this row if nan-policy is DropRow and value is NaN
        if matches!(nan_policy, NanPolicy::DropRow) && r.cov.is_nan() {
            continue;
        }
        match r.orig_idx {
            Some(idx) => {
                // Windowed positional
                match nan_policy {
                    NanPolicy::WriteLiteralNaN if r.cov.is_nan() => {
                        writeln!(w, "{}\t{}\t{}\t{}\t{}", r.chr, r.start, r.end, "NaN", idx)?
                    }
                    NanPolicy::WriteEmptyCell if r.cov.is_nan() => {
                        writeln!(w, "{}\t{}\t{}\t\t{}", r.chr, r.start, r.end, idx)?
                    }
                    _ => writeln!(w, "{}\t{}\t{}\t{}\t{}", r.chr, r.start, r.end, r.cov, idx)?,
                }
            }
            None => {
                // Whole positional
                match nan_policy {
                    NanPolicy::WriteLiteralNaN if r.cov.is_nan() => {
                        writeln!(w, "{}\t{}\t{}\t{}", r.chr, r.start, r.end, "NaN")?
                    }
                    NanPolicy::WriteEmptyCell if r.cov.is_nan() => {
                        writeln!(w, "{}\t{}\t{}\t{}", r.chr, r.start, r.end, "")?
                    }
                    _ => writeln!(w, "{}\t{}\t{}\t{}", r.chr, r.start, r.end, r.cov)?,
                }
            }
        }
    }

    Ok(())
}

/// Possible output types
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputKind {
    AggregateAverage,
    AggregateTotal,
    PosWhole,    // WholePositional
    PosWindowed, // PerWindow::OnlyIncludeThesePositions
}

#[inline]
fn output_suffix(kind: OutputKind) -> &'static str {
    match kind {
        OutputKind::AggregateAverage => ".avg.tsv",
        OutputKind::AggregateTotal => ".total.tsv",
        OutputKind::PosWhole => ".per_position.tsv", // Whole-genome positional
        OutputKind::PosWindowed => ".per_positions.tsv", // Windowed positional
    }
}

/// Detect the type of output based on outputs-per-chromosome hashmap
fn detect_output_kind(by_chr: &FxHashMap<String, CoverageOutput>) -> anyhow::Result<OutputKind> {
    use anyhow::{anyhow, bail};

    let mut kind: Option<OutputKind> = None;

    for (chr, out) in by_chr {
        let this = match out {
            CoverageOutput::WholePositional { .. } => OutputKind::PosWhole,
            CoverageOutput::PerWindow { action, results } => match action {
                CoverageWindowAction::Average => OutputKind::AggregateAverage,
                CoverageWindowAction::Total => OutputKind::AggregateTotal,
                CoverageWindowAction::OnlyIncludeThesePositions => {
                    if results
                        .iter()
                        .any(|wr| !matches!(wr.value, WindowValue::Positions(_)))
                    {
                        bail!(
                            "chrom {} has non-positional values under OnlyIncludeThesePositions",
                            chr
                        );
                    }
                    OutputKind::PosWindowed
                }
            },
        };

        match kind {
            None => kind = Some(this),
            Some(prev) if prev != this => {
                bail!(
                    "mixed output kinds across chromosomes: {:?} vs {:?}",
                    prev,
                    this
                )
            }
            _ => {}
        }
    }

    kind.ok_or_else(|| anyhow!("no outputs to write"))
}

/// Write a single TSV chosen by the detected output kind
///
/// - Aggregate (Average/Total) -> `write_aggregates_tsv`
/// - Positional (whole or windowed) -> `write_positions_tsv`
///
/// Parameters
/// ----------
/// - path: Where to write the TSV
/// - by_chr: Map `chr -> CoverageOutput`
/// - pos_nan_policy: How to handle masked positions in positional outputs
pub fn write_outputs_auto<P: AsRef<Path>>(
    path: P,
    by_chr: &FxHashMap<String, CoverageOutput>,
    pos_nan_policy: NanPolicy,
) -> anyhow::Result<()> {
    match detect_output_kind(by_chr)? {
        OutputKind::AggregateAverage => {
            write_aggregates_tsv(path, by_chr, CoverageWindowAction::Average)
        }
        OutputKind::AggregateTotal => {
            write_aggregates_tsv(path, by_chr, CoverageWindowAction::Total)
        }
        OutputKind::PosWhole | OutputKind::PosWindowed => {
            write_positions_tsv(path, by_chr, pos_nan_policy)
        }
    }
}

/// Write a single TSV chosen by the detected output kind
///
/// Accepts a path-prefix and appends the last part of the
/// path based on output type.
///
/// Uses these writers:
/// - Aggregate (Average/Total) -> `write_aggregates_tsv`
/// - Positional (whole or windowed) -> `write_positions_tsv`
///
/// Parameters
/// ----------
/// - prefix: Path prefix for where to write the TSV
/// - by_chr: Map `chr -> CoverageOutput`
/// - pos_nan_policy: How to handle masked positions in positional outputs
///
/// Returns
/// -------
/// final_path: The final path where the output was written
pub fn write_outputs_auto_with_prefix<P: AsRef<Path>>(
    prefix: P,
    by_chr: &FxHashMap<String, CoverageOutput>,
    pos_policy: NanPolicy,
) -> anyhow::Result<PathBuf> {
    let kind = detect_output_kind(by_chr)?;

    let final_path = {
        let ps = prefix.as_ref().to_string_lossy();
        PathBuf::from(format!("{}{}", ps, output_suffix(kind)))
    };

    match kind {
        OutputKind::AggregateAverage => {
            write_aggregates_tsv(&final_path, by_chr, CoverageWindowAction::Average)?
        }
        OutputKind::AggregateTotal => {
            write_aggregates_tsv(&final_path, by_chr, CoverageWindowAction::Total)?
        }
        OutputKind::PosWhole | OutputKind::PosWindowed => {
            write_positions_tsv(&final_path, by_chr, pos_policy)?
        }
    }

    Ok(final_path)
}
