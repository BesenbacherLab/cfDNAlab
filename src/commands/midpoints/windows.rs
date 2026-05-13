use crate::shared::{
    bam::Contigs,
    bed::GroupedWindows,
    interval::{IndexedInterval, Interval},
};
use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;

/// Get window length and ensure it's the same for ALL windows.
pub fn ensure_uniform_window_len(
    windows_by_chr: &FxHashMap<String, GroupedWindows>,
) -> Result<usize> {
    let mut reference_len: Option<usize> = None;

    for (chr, gw) in windows_by_chr {
        for window in &gw.windows {
            let start = window.start();
            let end = window.end();
            let len = end.checked_sub(start).with_context(|| {
                format!("Invalid window on {chr}: end ({end}) < start ({start})")
            })? as usize;

            match reference_len {
                None => reference_len = Some(len),
                Some(ref_len) if (len) != ref_len => {
                    anyhow::bail!(
                        "Non-uniform window length detected on {chr}: [{start},{end}) has len {}, expected {}",
                        len,
                        ref_len
                    );
                }
                _ => {}
            }
        }
    }

    reference_len.context("No windows found when checking uniform window length")
}

#[derive(Debug, Default, Clone, Copy, Eq, PartialEq)]
pub(crate) struct MidpointIntervalStats {
    pub loaded_after_chromosome_filtering: usize,
    pub dropped_by_blacklist_prefilter: usize,
    pub retained_for_counting: usize,
}

/// Prepare midpoint windows for counting without changing the public output span.
///
/// The input BED intervals define the output profile positions. Optional smoothing needs extra
/// context on both sides, so this helper expands retained windows by the derived smoothing flank
/// only for counting. Later post-processing trims the flank away again before writing output.
///
/// When blacklists are supplied, interval-level prefiltering removes windows whose output span plus
/// the safety margin touches a blacklist. Fragment-level blacklist filtering still runs later in
/// the normal counting loop.
///
/// This compacts each chromosome vector in place. That preserves the loader's group indices while
/// avoiding a second full allocation for large site sets.
///
/// The input window indices are transferred directly to the new intervals and remains intact.
pub(crate) fn prepare_count_windows(
    mut windows_by_chr: FxHashMap<String, Vec<IndexedInterval<u64>>>,
    contigs: &Contigs,
    blacklist_by_chr: &FxHashMap<String, Vec<Interval<u64>>>,
    smoothing_flank: u32,
    blacklist_margin: u64,
    use_blacklist_prefilter: bool,
) -> Result<(
    FxHashMap<String, Vec<IndexedInterval<u64>>>,
    MidpointIntervalStats,
)> {
    let mut stats = MidpointIntervalStats::default();
    let flank = u64::from(smoothing_flank);

    for (chromosome, windows) in windows_by_chr.iter_mut() {
        stats.loaded_after_chromosome_filtering += windows.len();
        let chrom_len = contigs
            .contigs
            .get(chromosome)
            .map(|(_, len)| u64::from(*len))
            .with_context(|| format!("missing contig length for {chromosome}"))?;
        let blacklist = blacklist_by_chr
            .get(chromosome)
            .map(Vec::as_slice)
            .unwrap_or(&[]);
        let mut bl_ptr = 0usize;
        let original_len = windows.len();
        // Compact the vector in place with separate read and write positions
        //
        // `window_idx` walks the original vector from left to right. `write_idx` marks the next
        // position in the retained prefix. Kept intervals are expanded and written to
        // `windows[write_idx]`, then `write_idx` advances
        //
        // Dropped intervals are not written. That is the actual filtering step: `write_idx` stays
        // put, so the next kept interval overwrites the dropped slot or another stale processed
        // slot. After the scan, every retained interval is in `windows[..write_idx]`
        let mut write_idx = 0usize;

        for window_idx in 0..original_len {
            let window = windows[window_idx];
            if window.end() > chrom_len {
                bail!(
                    "Invalid midpoint interval {chromosome}:{}-{} extends beyond chromosome length {}. Is it from the same assembly? Remove or correct this BED row.",
                    window.start(),
                    window.end(),
                    chrom_len
                );
            }

            if use_blacklist_prefilter
                && interval_with_margin_overlaps_blacklist(
                    window.interval,
                    blacklist_margin,
                    blacklist,
                    &mut bl_ptr,
                )
            {
                stats.dropped_by_blacklist_prefilter += 1;
                // Skip the write to remove this interval from the retained prefix
                continue;
            }

            let expanded_start = window.start().checked_sub(flank).ok_or_else(|| {
                anyhow::anyhow!(
                    "Cannot smooth interval {chromosome}:{}-{} because it is within {} bp of the chromosome start. Use a smaller --smooth window or remove this interval.",
                    window.start(),
                    window.end(),
                    flank
                )
            })?;
            let expanded_end = window.end().checked_add(flank).ok_or_else(|| {
                anyhow::anyhow!(
                    "expanded interval end overflow for {chromosome}:{}-{} with flank {}",
                    window.start(),
                    window.end(),
                    flank
                )
            })?;
            if expanded_end > chrom_len {
                bail!(
                    "Cannot smooth interval {chromosome}:{}-{} because it is within {} bp of the chromosome end. Use a smaller --smooth window or remove this interval.",
                    window.start(),
                    window.end(),
                    flank
                );
            }

            windows[write_idx] = IndexedInterval::new(expanded_start, expanded_end, window.idx())?;
            write_idx += 1;
        }

        // Drop the stale tail. It contains blacklisted intervals and old copies left behind after
        // retained intervals were moved forward
        windows.truncate(write_idx);
        stats.retained_for_counting += write_idx;
    }

    Ok((windows_by_chr, stats))
}

fn interval_with_margin_overlaps_blacklist(
    interval: Interval<u64>,
    margin: u64,
    blacklist: &[Interval<u64>],
    bl_ptr: &mut usize,
) -> bool {
    if blacklist.is_empty() {
        return false;
    }

    let start = interval.start().saturating_sub(margin);
    let end = interval.end().saturating_add(margin);

    // Windows arrive sorted by start and this command requires a uniform interval length
    // Because `margin` is constant for the run, the margin-expanded start coordinate is also
    // monotonic. Once a blacklist interval ends before that start, it cannot overlap this or any
    // later interval, so the linear-sweep pointer can advance without a lookback
    while *bl_ptr < blacklist.len() && blacklist[*bl_ptr].end() <= start {
        *bl_ptr += 1;
    }

    // `load_blacklist_map` returns sorted, merged blacklist intervals. If the first interval that
    // has not ended starts at or beyond the margin-expanded window end, later intervals start even
    // farther right and therefore cannot overlap either
    *bl_ptr < blacklist.len() && blacklist[*bl_ptr].start() < end
}

#[cfg(test)]
mod tests {
    include!("windows_tests.rs");
}
