use crate::commands::prepare_windows::labels::LabelTuple;
use crate::commands::prepare_windows::prepare_windows::FinalWindow;
use anyhow::{Context, Result, bail};
use std::io::Write;

const TUPLE_SEPARATOR: char = ';';
const TUPLE_SEPARATOR_STR: &str = ";";
const PART_SEPARATOR: char = '|';

/// Window representation stored in intermediate files for filtering passes.
#[derive(Clone, Debug)]
pub struct IntermediateWindow {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub label_tuples: Vec<LabelTuple>,
}

/// Write one intermediate window record.
pub fn write_intermediate_window<W: Write>(
    writer: &mut W,
    window: &IntermediateWindow,
    separator: char,
) -> Result<()> {
    let tuples = serialize_label_tuples(&window.label_tuples);
    writeln!(
        writer,
        "{}{sep}{}{sep}{}{sep}{}",
        window.chrom,
        window.start,
        window.end,
        tuples,
        sep = separator
    )?;
    Ok(())
}

/// Write intermediate windows as `chrom start end tuples`.
///
/// The tuple field stores `input|near-side|near-name|bin|cluster` values
/// joined by `;` for multiple tuples.
///
/// Parameters
/// ----------
/// - `writer`:
///     Output writer for the intermediate file.
/// - `windows`:
///     Window slice to serialize.
/// - `separator`:
///     Column separator character.
///
/// Returns
/// -------
/// `Ok(())` on success or an error if writing fails.
pub fn write_intermediate_windows<W: Write>(
    writer: &mut W,
    windows: &[FinalWindow],
    separator: char,
) -> Result<()> {
    for window in windows {
        let tuples = serialize_label_tuples(&window.label_tuples);
        writeln!(
            writer,
            "{}{sep}{}{sep}{}{sep}{}",
            window.chrom.as_ref(),
            window.resized_start,
            window.resized_end,
            tuples,
            sep = separator
        )?;
    }
    Ok(())
}

/// Parse a single intermediate line into a window.
pub fn parse_intermediate_line(line: &str, separator: char) -> Result<IntermediateWindow> {
    let mut fields = line.splitn(4, separator);
    let chrom = fields
        .next()
        .context("Missing chrom field")?
        .trim()
        .to_string();
    let start: u32 = fields
        .next()
        .context("Missing start field")?
        .trim()
        .parse()
        .context("Invalid start field")?;
    let end: u32 = fields
        .next()
        .context("Missing end field")?
        .trim()
        .parse()
        .context("Invalid end field")?;
    let tuples_raw = fields.next().context("Missing label tuples field")?.trim();
    let label_tuples = parse_label_tuples(tuples_raw)?;

    Ok(IntermediateWindow {
        chrom,
        start,
        end,
        label_tuples,
    })
}

fn serialize_label_tuples(tuples: &[LabelTuple]) -> String {
    let mut parts: Vec<String> = Vec::with_capacity(tuples.len());
    for tuple in tuples {
        let mut row = String::new();
        row.push_str(tuple.input.as_str());
        row.push(PART_SEPARATOR);
        row.push_str(tuple.near_side.as_deref().unwrap_or(""));
        row.push(PART_SEPARATOR);
        row.push_str(tuple.near_name.as_deref().unwrap_or(""));
        row.push(PART_SEPARATOR);
        row.push_str(tuple.bin.as_deref().unwrap_or(""));
        row.push(PART_SEPARATOR);
        row.push_str(tuple.cluster.as_deref().unwrap_or(""));
        parts.push(row);
    }
    parts.join(TUPLE_SEPARATOR_STR)
}

fn parse_label_tuples(raw: &str) -> Result<Vec<LabelTuple>> {
    if raw.is_empty() {
        return Ok(Vec::new());
    }

    let mut tuples: Vec<LabelTuple> = Vec::new();
    for tuple_raw in raw.split(TUPLE_SEPARATOR) {
        let fields: Vec<&str> = tuple_raw.split(PART_SEPARATOR).collect();
        if fields.len() != 5 {
            bail!("Invalid tuple field '{}'", tuple_raw);
        }
        let input = fields[0].to_string();
        let near_side = if fields[1].is_empty() {
            None
        } else {
            Some(fields[1].to_string())
        };
        let near_name = if fields[2].is_empty() {
            None
        } else {
            Some(fields[2].to_string())
        };
        let bin = if fields[3].is_empty() {
            None
        } else {
            Some(fields[3].to_string())
        };
        let cluster = if fields[4].is_empty() {
            None
        } else {
            Some(fields[4].to_string())
        };

        tuples.push(LabelTuple {
            input,
            near_side,
            near_name,
            bin,
            cluster,
        });
    }
    Ok(tuples)
}
