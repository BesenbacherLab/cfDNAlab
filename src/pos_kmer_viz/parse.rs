use std::error::Error;
use std::fmt;

use anyhow::{Context, Result, anyhow};

use super::model::{Anchor, LinearRange, MidRange, NearestRange, PositionsSpec};

const LINEAR_EXAMPLE: &str = "--positions 1..10";
const NEAREST_EXAMPLE: &str = "--positions ..half";
const MID_EXAMPLE: &str = "--positions -10..10";

/// Error type used when the range grammar does not match expectations.
#[derive(Debug)]
pub struct RangeParseError {
    message: String,
}

impl RangeParseError {
    pub fn new(message: impl Into<String>, example: &'static str) -> Self {
        let mut msg = message.into();
        if !msg.contains("Example:") {
            msg.push_str(" Example: ");
            msg.push_str(example);
        }
        Self { message: msg }
    }
}

impl fmt::Display for RangeParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl Error for RangeParseError {}

pub fn parse_positions(anchor: Anchor, input: &str) -> Result<PositionsSpec, RangeParseError> {
    match anchor {
        Anchor::Left | Anchor::Right | Anchor::PerEnd | Anchor::Span => {
            parse_linear_range(input).map(PositionsSpec::Linear)
        }
        Anchor::Nearest => parse_nearest_range(input).map(PositionsSpec::Nearest),
        Anchor::Mid => parse_mid_range(input).map(PositionsSpec::Mid),
    }
}

pub fn parse_lengths(list: Option<&str>, range: Option<&str>) -> Result<Vec<u32>> {
    if let Some(list) = list {
        let mut lengths = Vec::new();
        for item in list.split(',') {
            let trimmed = item.trim();
            if trimmed.is_empty() {
                continue;
            }
            let value = trimmed
                .parse::<u32>()
                .with_context(|| format!("\"{}\" is not a positive integer", trimmed))?;
            if value == 0 {
                return Err(anyhow!(
                    "fragment lengths must be positive (example: 120,150)"
                ));
            }
            lengths.push(value);
        }
        return Ok(lengths);
    }

    if let Some(range) = range {
        return parse_length_range(range);
    }

    Ok(Vec::new())
}

fn parse_length_range(spec: &str) -> Result<Vec<u32>> {
    let parts: Vec<&str> = spec.split(':').collect();
    if parts.len() < 2 || parts.len() > 3 {
        return Err(anyhow!(
            "length range must follow MIN:MAX[:STEP] (example: 100:200:25)"
        ));
    }
    let min = parts[0]
        .trim()
        .parse::<u32>()
        .with_context(|| format!("\"{}\" is not a positive integer", parts[0].trim()))?;
    let max = parts[1]
        .trim()
        .parse::<u32>()
        .with_context(|| format!("\"{}\" is not a positive integer", parts[1].trim()))?;
    if min == 0 || max == 0 {
        return Err(anyhow!(
            "length range values must be positive (example: 100:200:25)"
        ));
    }
    if min > max {
        return Err(anyhow!(
            "length range requires MIN <= MAX (example: 120:200:10)"
        ));
    }
    let step = if parts.len() == 3 {
        parts[2]
            .trim()
            .parse::<u32>()
            .with_context(|| format!("\"{}\" is not a positive integer", parts[2].trim()))?
    } else {
        10
    };
    if step == 0 {
        return Err(anyhow!(
            "length range step must be positive (example: 80:200:20)"
        ));
    }
    let mut result = Vec::new();
    let mut current = min;
    while current <= max {
        result.push(current);
        match current.checked_add(step) {
            Some(next) if next > current => current = next,
            _ => break,
        }
    }
    if *result.last().unwrap_or(&0) != max && (max - result.last().unwrap_or(&min)) % step != 0 {
        // Ensure the upper bound is included if it falls on the lattice.
        let mut last = *result.last().unwrap();
        while let Some(next) = last.checked_add(step) {
            if next > max {
                break;
            }
            result.push(next);
            last = next;
        }
        if *result.last().unwrap() != max {
            result.push(max);
        }
    }
    Ok(result)
}

fn parse_linear_range(input: &str) -> Result<LinearRange, RangeParseError> {
    if let Some((start_str, end_str)) = input.split_once("..") {
        if start_str.is_empty() && end_str.is_empty() {
            return Err(RangeParseError::new(
                "expected bounds around '..' (example: 1..10)",
                LINEAR_EXAMPLE,
            ));
        }
        if start_str.is_empty() {
            let end = parse_positive(end_str, "end", LINEAR_EXAMPLE)?;
            return Ok(LinearRange::To { end });
        }
        if end_str.is_empty() {
            let start = parse_positive(start_str, "start", LINEAR_EXAMPLE)?;
            return Ok(LinearRange::From { start });
        }
        if let Some(trim_str) = end_str.strip_prefix('-') {
            let start = parse_positive(start_str, "start", LINEAR_EXAMPLE)?;
            let trim = parse_positive(trim_str, "other-end trim", LINEAR_EXAMPLE)?;
            return Ok(LinearRange::TrimOtherEnd {
                start,
                other_end_trim: trim,
            });
        }
        let start = parse_positive(start_str, "start", LINEAR_EXAMPLE)?;
        let end = parse_positive(end_str, "end", LINEAR_EXAMPLE)?;
        if start > end {
            return Err(RangeParseError::new(
                "start must be <= end (example: 10..30)",
                LINEAR_EXAMPLE,
            ));
        }
        return Ok(LinearRange::Closed { start, end });
    }

    if input.contains(':') {
        return Err(RangeParseError::new(
            "colon ranges were replaced by '..' syntax (examples: 10.., ..25, 5..-5)",
            LINEAR_EXAMPLE,
        ));
    }

    Err(RangeParseError::new(
        "unsupported positions format for this anchor (examples: 1..10, 10.., ..25, 5..-5)",
        LINEAR_EXAMPLE,
    ))
}

fn parse_nearest_range(input: &str) -> Result<NearestRange, RangeParseError> {
    if let Some((start_str, tail)) = input.split_once("..half") {
        let minus = parse_optional_minus(tail, NEAREST_EXAMPLE)?;
        if start_str.is_empty() {
            return Ok(NearestRange::ToHalf { minus });
        }
        let start = parse_positive(start_str, "start", NEAREST_EXAMPLE)?;
        return Ok(NearestRange::FromToHalf { start, minus });
    }

    if input.contains("half") {
        return Err(RangeParseError::new(
            "half-relative ranges must follow A..half[-K] or ..half[-K]",
            NEAREST_EXAMPLE,
        ));
    }

    if let Some((start_str, end_str)) = input.split_once("..") {
        if start_str.is_empty() && end_str.is_empty() {
            return Err(RangeParseError::new(
                "expected bounds around '..' (example: 1..10)",
                NEAREST_EXAMPLE,
            ));
        }
        if start_str.is_empty() {
            if end_str.is_empty() {
                return Err(RangeParseError::new(
                    "expected digits after '..'",
                    NEAREST_EXAMPLE,
                ));
            }
            if let Some(tail) = end_str.strip_prefix("half") {
                let minus = parse_optional_minus(tail, NEAREST_EXAMPLE)?;
                return Ok(NearestRange::ToHalf { minus });
            }
            let end = parse_positive(end_str, "end", NEAREST_EXAMPLE)?;
            return Ok(NearestRange::Closed { start: 1, end });
        }
        if end_str.is_empty() {
            let start = parse_positive(start_str, "start", NEAREST_EXAMPLE)?;
            return Ok(NearestRange::From { start });
        }
        let start = parse_positive(start_str, "start", NEAREST_EXAMPLE)?;
        if let Some(tail) = end_str.strip_prefix("half") {
            let minus = parse_optional_minus(tail, NEAREST_EXAMPLE)?;
            return Ok(NearestRange::FromToHalf { start, minus });
        }
        let end = parse_positive(end_str, "end", NEAREST_EXAMPLE)?;
        if start > end {
            return Err(RangeParseError::new(
                "start must be <= end (example: 1..10)",
                NEAREST_EXAMPLE,
            ));
        }
        return Ok(NearestRange::Closed { start, end });
    }

    if input.contains(':') {
        return Err(RangeParseError::new(
            "colon ranges were replaced by '..' syntax (examples: 10.., ..10, ..half)",
            NEAREST_EXAMPLE,
        ));
    }

    Err(RangeParseError::new(
        "unsupported positions format for nearest anchor",
        NEAREST_EXAMPLE,
    ))
}

fn parse_mid_range(input: &str) -> Result<MidRange, RangeParseError> {
    if let Some(rest) = input.strip_prefix("..") {
        if rest.is_empty() {
            return Err(RangeParseError::new(
                "expected digits after '..'",
                MID_EXAMPLE,
            ));
        }
        let cleaned = rest.strip_prefix('+').unwrap_or(rest);
        if cleaned.is_empty() {
            return Err(RangeParseError::new(
                "expected digits after '..'",
                MID_EXAMPLE,
            ));
        }
        let pos = parse_positive(cleaned, "positive bound", MID_EXAMPLE)?;
        return Ok(MidRange::RightOpen { pos });
    }

    if let Some(rest) = input.strip_suffix("..") {
        if !rest.starts_with('-') || rest.len() <= 1 {
            return Err(RangeParseError::new("expected form '-M..'", MID_EXAMPLE));
        }
        let neg = parse_positive(&rest[1..], "negative bound", MID_EXAMPLE)?;
        return Ok(MidRange::LeftOpen { neg });
    }

    if let Some(idx) = input.find("..") {
        let left = &input[..idx];
        let right = &input[idx + 2..];
        if !left.starts_with('-') || left.len() <= 1 || right.is_empty() {
            return Err(RangeParseError::new("expected form '-M..N'", MID_EXAMPLE));
        }
        let neg = parse_positive(&left[1..], "negative bound", MID_EXAMPLE)?;
        let cleaned = right.strip_prefix('+').unwrap_or(right);
        if cleaned.is_empty() {
            return Err(RangeParseError::new(
                "expected digits after '..'",
                MID_EXAMPLE,
            ));
        }
        let pos = parse_positive(cleaned, "positive bound", MID_EXAMPLE)?;
        return Ok(MidRange::Closed { neg, pos });
    }

    Err(RangeParseError::new(
        "unsupported positions format for mid anchor",
        MID_EXAMPLE,
    ))
}

fn parse_positive(value: &str, field: &str, example: &'static str) -> Result<u32, RangeParseError> {
    value
        .trim()
        .parse::<u32>()
        .map_err(|_| RangeParseError::new(format!("{} must be a positive integer", field), example))
        .and_then(|v| {
            if v == 0 {
                Err(RangeParseError::new(
                    format!("{} must be positive", field),
                    example,
                ))
            } else {
                Ok(v)
            }
        })
}

fn parse_optional_minus(rest: &str, example: &'static str) -> Result<u32, RangeParseError> {
    if rest.is_empty() {
        return Ok(0);
    }
    if !rest.starts_with('-') {
        return Err(RangeParseError::new("expected '-K' after half", example));
    }
    parse_positive(&rest[1..], "offset", example)
}
