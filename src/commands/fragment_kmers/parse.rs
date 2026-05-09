use std::fmt;
use std::{error::Error, num::NonZeroUsize};

use anyhow::Result;

use crate::commands::{
    cli_common::UnparsedPositionalSelectionSpec,
    fragment_kmers::positions::{LinearRange, MidRange, NearestRange, PositionsSpec},
};
use crate::shared::positioning::ReferenceFrame;
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

#[derive(Debug, Clone)]
pub struct PositionalSelectionSpec {
    pub frame: ReferenceFrame,
    pub positions: PositionsSpec,
    pub positions_string: String,
    pub step: NonZeroUsize,
}

pub fn parse_positions(
    position_spec: &UnparsedPositionalSelectionSpec,
) -> Result<PositionalSelectionSpec, RangeParseError> {
    let step = NonZeroUsize::new(position_spec.step)
        .ok_or_else(|| RangeParseError::new("`--step` must be >= 1.", "--step 1"))?;

    let positions_spec = match position_spec.frame {
        ReferenceFrame::Left | ReferenceFrame::Right | ReferenceFrame::PerEnd => {
            parse_linear_range(&position_spec.positions).map(PositionsSpec::Linear)
        }
        ReferenceFrame::Nearest => {
            parse_nearest_range(&position_spec.positions).map(PositionsSpec::Nearest)
        }
        ReferenceFrame::Mid => parse_mid_range(&position_spec.positions).map(PositionsSpec::Mid),
    }?;

    Ok(PositionalSelectionSpec {
        frame: position_spec.frame,
        positions: positions_spec,
        positions_string: position_spec.positions.clone(),
        step,
    })
}

fn parse_linear_range(input: &str) -> Result<LinearRange, RangeParseError> {
    let input = input.trim();
    if input == ".." {
        return Ok(LinearRange::All);
    }
    if let Some((start_str, end_str)) = input.split_once("..") {
        if start_str.is_empty() && end_str.is_empty() {
            return Err(RangeParseError::new(
                "expected bounds around '..' (example: 1..10)",
                LINEAR_EXAMPLE,
            ));
        }
        if start_str.is_empty() {
            if let Some(tail) = end_str.strip_prefix("half") {
                let minus = parse_optional_minus(tail, LINEAR_EXAMPLE)?;
                return Ok(LinearRange::ToHalf { minus });
            }
            let end = parse_positive(end_str, "end", LINEAR_EXAMPLE)?;
            return Ok(LinearRange::To { end });
        }
        if end_str.is_empty() {
            let start = parse_positive(start_str, "start", LINEAR_EXAMPLE)?;
            return Ok(LinearRange::From { start });
        }
        if let Some(tail) = end_str.strip_prefix("half") {
            let start = parse_positive(start_str, "start", LINEAR_EXAMPLE)?;
            let minus = parse_optional_minus(tail, LINEAR_EXAMPLE)?;
            return Ok(LinearRange::FromToHalf { start, minus });
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
        "unsupported positions format for this frame (examples: 1..10, 10.., ..25, 5..-5, ..half)",
        LINEAR_EXAMPLE,
    ))
}

fn parse_nearest_range(input: &str) -> Result<NearestRange, RangeParseError> {
    let input = input.trim();
    if input == ".." {
        return Ok(NearestRange::All);
    }
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
        "unsupported positions format for nearest frame",
        NEAREST_EXAMPLE,
    ))
}

fn parse_mid_range(input: &str) -> Result<MidRange, RangeParseError> {
    let input = input.trim();
    if input == ".." {
        return Ok(MidRange::All);
    }
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
        "unsupported positions format for mid frame",
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
