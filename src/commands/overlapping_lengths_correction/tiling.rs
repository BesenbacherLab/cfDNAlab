//! Tile-local fragment scanning and positional sufficient-statistic accumulation.

use anyhow::{Context, Result, ensure};
use rust_htslib::bam::{Read, Record};

use crate::{
    commands::{
        counters::FCoverageCounters,
        gc_bias::{correct::GCCorrector, counting::build_gc_prefixes},
    },
    shared::{
        bam::create_chromosome_reader,
        fragment::segment_fragment::FragmentWithSegments,
        fragment_iterators::fragments_with_segments_from_bam,
        gc_tag::ClassifiedGCTagWeight,
        interval::Interval,
        read::{default_include_read_paired_end, default_include_read_unpaired},
        reference::ReferenceReader,
        scale_genome::ScalingBin,
        tiled_run::Tile,
    },
};

use super::{config::OverlappingLengthsCorrectionConfig, reducer::OverlappingLengthStatistics};

/// Additive output from a single tile core.
///
/// The dense difference arrays are consumed inside `process_tile`. Only these compact statistics
/// and command counters are returned to the parallel reducer.
#[derive(Debug)]
pub(crate) struct OverlappingLengthTileResult {
    /// Length-bin signal sums, base counts, and raw-depth frequencies for the tile core.
    pub(crate) statistics: OverlappingLengthStatistics,
    /// Read, fragment, and GC-filtering counters observed while processing the tile.
    pub(crate) counters: FCoverageCounters,
}

/// Collect sufficient statistics for one tile core.
///
/// Raw depth and fragment-length sums are discrete prefix sums. The optional signal prefix differs
/// only when GC correction is active. Genomic scaling is positional and is therefore applied while
/// scanning the finalized signal, never to the raw depth used by the mixture model.
///
/// The BAM fetch includes a maximum-fragment-length halo, while all returned statistics are
/// restricted to the non-overlapping tile core. This makes tile results additive and prevents
/// boundary bases from being counted twice.
///
/// Parameters
/// ----------
/// - `config`:
///   Fragment filters and optional observed-signal corrections.
/// - `tile`:
///   Core and halo coordinates for this worker.
/// - `length_bin_edges`:
///   Shared bin edges copied into the tile-local reducer.
/// - `blacklist`:
///   Sorted chromosome intervals excluded from the returned statistics.
/// - `scaling_bins`:
///   Sorted positional scaling intervals applied to observed coverage only.
/// - `gc_corrector`, `gc_tag`, `reference_reader`:
///   Mutually consistent GC correction inputs resolved by the command driver.
///
/// Returns
/// -------
/// - `OverlappingLengthTileResult`:
///   Compact sufficient statistics and counters for this tile core.
pub(crate) fn process_tile(
    config: &OverlappingLengthsCorrectionConfig,
    tile: &Tile,
    length_bin_edges: Vec<f64>,
    blacklist: &[Interval<u64>],
    scaling_bins: &[ScalingBin],
    gc_corrector: Option<GCCorrector>,
    gc_tag: Option<&str>,
    reference_reader: Option<&mut ReferenceReader>,
) -> Result<OverlappingLengthTileResult> {
    // Open a reader per worker because rust-htslib readers are stateful and not shared across tiles
    let (mut reader, bam_tid, _chromosome_length) =
        create_chromosome_reader(&config.ioc.bam, &tile.chr)?;
    tile.ensure_matches_bam_tid(bam_tid)?;
    let (fetch_start, fetch_end) = tile.fetch.as_tuple();
    reader
        .fetch((tile.tid, i64::from(fetch_start), i64::from(fetch_end)))
        .with_context(|| {
            format!(
                "fetch {} {}-{} for overlapping fragment length fitting",
                tile.chr, fetch_start, fetch_end
            )
        })?;

    // File-based GC correction needs reference prefixes over the complete fragment fetch halo
    let gc_prefixes = if gc_corrector.is_some() {
        let reference_reader = reference_reader
            .context("file-based GC correction requires a staged reference reader")?;
        let sequence = reference_reader
            .read_seq_in_range(&tile.chr, fetch_start as usize..fetch_end as usize)?;
        Some(build_gc_prefixes(&sequence))
    } else {
        None
    };

    // Integer arrays preserve exact raw depth and fragment length sums during model training
    let core_length = tile.core.len() as usize;
    let mut raw_depth_delta = vec![0_i64; core_length + 1];
    let mut fragment_length_sum_delta = vec![0_i64; core_length + 1];
    // Avoid a third dense array when the observed signal is identical to raw coverage
    let mut corrected_signal_delta =
        (gc_corrector.is_some() || gc_tag.is_some()).then(|| vec![0.0_f64; core_length + 1]);
    let lengths = config.fragment_lengths();
    let unpaired = config.unpaired.reads_are_fragments;
    // Apply the same read-level filtering vocabulary as fcoverage
    let include_read: Box<dyn Fn(&Record) -> bool + Send + Sync> = if unpaired {
        let minimum_mapq = config.min_mapq;
        Box::new(move |record| default_include_read_unpaired(record, minimum_mapq))
    } else {
        let minimum_mapq = config.min_mapq;
        let require_proper_pair = config.require_proper_pair;
        Box::new(move |record| {
            default_include_read_paired_end(record, require_proper_pair, minimum_mapq)
        })
    };
    let mut fragments = fragments_with_segments_from_bam(
        reader
            .records()
            .map(|result| result.map_err(anyhow::Error::from)),
        move |record| include_read(record),
        1,
        !config.ignore_gap,
        gc_tag.map(str::as_bytes),
        move |fragment: &FragmentWithSegments| lengths.contains(fragment.len()),
        unpaired,
    )
    .with_local_counters();
    let mut counters = FCoverageCounters::default();

    // Add complete directional fragments to core-clipped difference arrays
    for fragment_result in fragments.by_ref() {
        let fragment = fragment_result.context("reading fragment")?;
        if fragment.start() < fetch_start || fragment.end() > fetch_end {
            continue;
        }
        // GC affects the observed signal but never the raw overlap context or depth histogram
        let signal_weight = if let Some(corrector) = &gc_corrector {
            let relative_interval = fragment
                .interval
                .try_to_u64()?
                .shift_left(fetch_start as u64)?;
            match corrector.correct_fragment(
                relative_interval,
                gc_prefixes
                    .as_ref()
                    .context("GC prefixes are missing during correction")?,
            )? {
                Some(weight) => Some(weight),
                None => {
                    counters.gc_failed_fragments += 1;
                    if config.gc.neutralize_invalid_gc {
                        Some(1.0)
                    } else {
                        None
                    }
                }
            }
        } else if gc_tag.is_some() {
            match fragment.gc_tag.classify()? {
                ClassifiedGCTagWeight::Usable(weight) => Some(f64::from(weight)),
                ClassifiedGCTagWeight::Missing => {
                    counters.gc_failed_fragments += 1;
                    counters.gc_missing_tags += 1;
                    if config.gc.neutralize_invalid_gc {
                        Some(1.0)
                    } else {
                        None
                    }
                }
                ClassifiedGCTagWeight::Invalid { out_of_range } => {
                    counters.gc_failed_fragments += 1;
                    if out_of_range {
                        counters.gc_out_of_range_tags += 1;
                    }
                    if config.gc.neutralize_invalid_gc {
                        Some(1.0)
                    } else {
                        None
                    }
                }
            }
        } else {
            Some(1.0)
        };

        // Every counted segment receives the original full directional fragment length
        let fragment_length = i64::from(fragment.len());
        let mut counted = false;
        if let Some(segments) = &fragment.segments {
            for segment in segments {
                counted |= add_segment_to_core_deltas(
                    *segment,
                    tile.core,
                    fragment_length,
                    signal_weight,
                    &mut raw_depth_delta,
                    &mut fragment_length_sum_delta,
                    corrected_signal_delta.as_mut(),
                );
            }
        } else {
            counted = add_segment_to_core_deltas(
                fragment.interval,
                tile.core,
                fragment_length,
                signal_weight,
                &mut raw_depth_delta,
                &mut fragment_length_sum_delta,
                corrected_signal_delta.as_mut(),
            );
        }
        if counted {
            counters.base.counted_fragments += 1;
        }
    }
    counters.add_from_snapshot(fragments.counters_snapshot());

    // Resolve the three difference arrays in a single left-to-right scan of the tile core
    let mut statistics = OverlappingLengthStatistics::new(length_bin_edges)?;
    let mut raw_depth = 0_i64;
    let mut fragment_length_sum = 0_i64;
    let mut corrected_signal = 0.0_f64;
    let mut blacklist_index =
        blacklist.partition_point(|interval| interval.end() <= u64::from(tile.core_start()));
    let mut scaling_index =
        scaling_bins.partition_point(|bin| bin.interval.end() <= u64::from(tile.core_start()));
    for local_position in 0..core_length {
        raw_depth += raw_depth_delta[local_position];
        fragment_length_sum += fragment_length_sum_delta[local_position];
        if let Some(delta) = &corrected_signal_delta {
            corrected_signal += delta[local_position];
        }
        ensure!(
            raw_depth >= 0 && fragment_length_sum >= 0,
            "overlapping fragment length prefix sum became negative"
        );
        // Uncovered positions have no defined average overlapping fragment length
        if raw_depth == 0 {
            continue;
        }
        let genomic_position = u64::from(tile.core_start()) + local_position as u64;
        while blacklist_index < blacklist.len()
            && blacklist[blacklist_index].end() <= genomic_position
        {
            blacklist_index += 1;
        }
        // Mask before binning so blacklisted bases affect neither fit nor depth frequencies
        if blacklist_index < blacklist.len()
            && blacklist[blacklist_index].start() <= genomic_position
            && genomic_position < blacklist[blacklist_index].end()
        {
            continue;
        }
        // Reuse integer raw depth logically when no GC-weighted signal array was allocated
        let mut observed_signal = if corrected_signal_delta.is_some() {
            ensure!(
                corrected_signal >= -1.0e-9,
                "GC-corrected coverage prefix became materially negative: {}",
                corrected_signal
            );
            corrected_signal.max(0.0)
        } else {
            raw_depth as f64
        };
        while scaling_index < scaling_bins.len()
            && scaling_bins[scaling_index].interval.end() <= genomic_position
        {
            scaling_index += 1;
        }
        // Positional scaling changes only the signal whose length bias is being estimated
        if let Some(bin) = scaling_bins.get(scaling_index)
            && bin.interval.start() <= genomic_position
            && genomic_position < bin.interval.end()
        {
            observed_signal *= f64::from(bin.weight_per_base);
        }
        // The denominator remains raw integer fragment count under every correction mode
        let average_length = fragment_length_sum as f64 / raw_depth as f64;
        statistics.add_position(
            average_length,
            observed_signal,
            u32::try_from(raw_depth).context("raw positional fragment depth exceeds u32")?,
        )?;
    }

    Ok(OverlappingLengthTileResult {
        statistics,
        counters,
    })
}

#[allow(clippy::too_many_arguments)]
/// Add one covered segment to tile-core difference arrays.
///
/// The segment is clipped to the core before any write. Raw depth and fragment length sums are
/// always updated. The observed-signal array is updated only when it exists and the fragment has a
/// usable or neutralized GC weight.
///
/// Returning `true` means at least one base of this segment belongs to the tile core.
fn add_segment_to_core_deltas(
    segment: Interval<u32>,
    core: Interval<u32>,
    fragment_length: i64,
    signal_weight: Option<f64>,
    raw_depth_delta: &mut [i64],
    fragment_length_sum_delta: &mut [i64],
    corrected_signal_delta: Option<&mut Vec<f64>>,
) -> bool {
    let Some(clipped) = segment.clip_to(core) else {
        return false;
    };
    let start = (clipped.start() - core.start()) as usize;
    let end = (clipped.end() - core.start()) as usize;
    raw_depth_delta[start] += 1;
    raw_depth_delta[end] -= 1;
    fragment_length_sum_delta[start] += fragment_length;
    fragment_length_sum_delta[end] -= fragment_length;
    if let (Some(delta), Some(signal_weight)) = (corrected_signal_delta, signal_weight) {
        delta[start] += signal_weight;
        delta[end] -= signal_weight;
    }
    true
}

#[cfg(test)]
mod tests {
    include!("tiling_tests.rs");
}
