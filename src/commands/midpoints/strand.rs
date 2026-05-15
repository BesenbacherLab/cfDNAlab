use crate::shared::bed::Strand;
use crate::shared::interval::Interval;
use anyhow::{Context, Result};

/// Map a genomic midpoint position to its profile-array index for one stranded window.
///
/// Forward and unstranded windows use the ordinary left-to-right genomic offset. Reverse windows
/// mirror the offset so the rightmost genomic base becomes profile position 0.
pub(crate) fn stranded_window_position(
    window: Interval<u64>,
    genomic_position: u64,
    strand: Strand,
) -> Result<usize> {
    let offset = match strand {
        Strand::Unstranded | Strand::Forward => genomic_position - window.start(),
        Strand::Reverse => (window.end() - 1) - genomic_position,
    };

    usize::try_from(offset).context("stranded window position does not fit in usize")
}

#[cfg(test)]
mod tests {
    include!("strand_tests.rs");
}
