use crate::Result;
use crate::shared::constants::MAX_SUPPORTED_FRAGMENT_LENGTH;
use crate::shared::fragment::segment_fragment::FragmentWithSegments;
use crate::shared::interval::{Interval, TouchingMergePolicy, merge_sorted_intervals};
use crate::shared::midpoint::midpoint_random_even_for_fragment;
use std::{fmt, str::FromStr};

const MIN_TRIM_TARGET_LENGTH: u32 = 1;

/// How `fcoverage` trims fragment coverage around the midpoint.
///
/// Both modes use an odd target length centered on the fragment midpoint. `AtMost` trims only
/// fragments longer than the target. `Exactly` additionally extends shorter fragments to the
/// target span.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FragmentSpanTrim {
    AtMost { target_length: u32 },
    Exactly { target_length: u32 },
}

impl FragmentSpanTrim {
    /// Return the requested odd outer span in base pairs.
    pub fn target_length(self) -> u32 {
        match self {
            Self::AtMost { target_length } | Self::Exactly { target_length } => target_length,
        }
    }

    /// Return whether shorter fragments should be extended to the target span.
    pub(crate) fn extends_shorter_fragments(self) -> bool {
        matches!(self, Self::Exactly { .. })
    }

    /// Validate a programmatically constructed trim rule.
    pub(crate) fn validate(self) -> std::result::Result<(), String> {
        let target_length = self.target_length();
        if target_length < MIN_TRIM_TARGET_LENGTH {
            return Err(format!(
                "trim target must be at least {} bp, got {target_length}",
                MIN_TRIM_TARGET_LENGTH
            ));
        }
        if target_length % 2 == 0 {
            return Err(format!(
                "trim target must be odd so it can be centered on one midpoint base, got {target_length}"
            ));
        }
        if target_length > MAX_SUPPORTED_FRAGMENT_LENGTH {
            return Err(format!(
                "trim target must be <= {MAX_SUPPORTED_FRAGMENT_LENGTH} bp, got {target_length}"
            ));
        }
        Ok(())
    }
}

impl fmt::Display for FragmentSpanTrim {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::AtMost { target_length } => write!(formatter, "at-most={target_length}"),
            Self::Exactly { target_length } => write!(formatter, "exactly={target_length}"),
        }
    }
}

impl FromStr for FragmentSpanTrim {
    type Err = String;

    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        let (mode, target_text) = value.trim().split_once('=').ok_or_else(|| {
            format!(
                "invalid trim specification '{value}'. Use 'at-most=<odd_bp>' or 'exactly=<odd_bp>'"
            )
        })?;
        let target_length = target_text.parse::<u32>().map_err(|_| {
            format!("invalid trim target '{target_text}'. Use an odd integer in base pairs")
        })?;
        let trim_to = if mode.eq_ignore_ascii_case("at-most") {
            Self::AtMost { target_length }
        } else if mode.eq_ignore_ascii_case("exactly") {
            Self::Exactly { target_length }
        } else {
            return Err(format!(
                "unsupported trim mode '{mode}'. Use 'at-most=<odd_bp>' or 'exactly=<odd_bp>'"
            ));
        };
        trim_to.validate()?;
        Ok(trim_to)
    }
}

/// Build the fragment segments used for coverage counting.
///
/// The fragment midpoint follows the shared deterministic midpoint rule. Because the requested
/// target is odd, the selected midpoint base has the same number of target bases on both sides.
/// Existing explicit segments remain authoritative. Trimming intersects them with the new span,
/// while extension adds only bases outside the original `pos` to `reference_end` span.
pub(crate) fn trimmed_counting_segments(
    fragment: &FragmentWithSegments,
    chromosome: &str,
    chromosome_length: u32,
    trim_to: FragmentSpanTrim,
) -> Result<Vec<Interval<u32>>> {
    let original_length = fragment.len();
    let target_length = trim_to.target_length();
    let fragment_needs_trimming = original_length > target_length;
    let fragment_needs_extension =
        original_length < target_length && trim_to.extends_shorter_fragments();

    if !fragment_needs_trimming && !fragment_needs_extension {
        return Ok(fragment.segments.as_ref().map_or_else(
            || vec![fragment.interval],
            |segments| segments.iter().copied().collect(),
        ));
    }

    let midpoint = midpoint_random_even_for_fragment(chromosome, fragment.start(), original_length);
    let half_target = target_length / 2;
    let trimmed_start = midpoint.saturating_sub(half_target);
    let trimmed_end =
        (midpoint as u64 + half_target as u64 + 1).min(chromosome_length as u64) as u32;
    let trimmed_interval = Interval::new(trimmed_start, trimmed_end)?;

    let Some(original_segments) = fragment.segments.as_ref() else {
        return Ok(vec![trimmed_interval]);
    };

    let mut countable_segments = Vec::with_capacity(original_segments.len() + 2);
    if fragment_needs_trimming {
        countable_segments.extend(
            original_segments
                .iter()
                .filter_map(|segment| segment.clip_to(trimmed_interval)),
        );
    } else {
        if trimmed_interval.start() < fragment.start() {
            countable_segments.push(Interval::new(trimmed_interval.start(), fragment.start())?);
        }
        countable_segments.extend(
            original_segments
                .iter()
                .filter_map(|segment| segment.clip_to(trimmed_interval)),
        );
        if trimmed_interval.end() > fragment.end() {
            countable_segments.push(Interval::new(fragment.end(), trimmed_interval.end())?);
        }
    }

    if countable_segments.is_empty() {
        return Ok(Vec::new());
    }

    countable_segments.sort_unstable_by_key(|segment| segment.start());
    let merged_segments =
        merge_sorted_intervals(countable_segments, TouchingMergePolicy::MergeTouching);
    Ok(merged_segments)
}

#[cfg(test)]
mod tests {
    include!("fragment_span_trim_tests.rs");
}
