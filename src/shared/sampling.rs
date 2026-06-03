use fxhash::FxHashMap;
use rand::Rng;

/// Estimate the sampling density given a desired number of positions and the maximum window length.
///
/// The density is the ratio `n_samples / total_possible_starts` where `total_possible_starts`
/// counts all valid start indices across chromosomes (`len - max_window_len + 1` per chrom).
/// Returns `0.0` when no valid starts exist or when `n_samples` is zero.
pub(crate) fn sampling_density(
    chrom_sizes: &FxHashMap<String, usize>,
    max_window_len: u64,
    n_samples: usize,
) -> f64 {
    if n_samples == 0 {
        return 0.0;
    }

    let mut total_possible = 0u64;
    for &len in chrom_sizes.values() {
        if len as u64 >= max_window_len {
            total_possible += len as u64 - max_window_len + 1;
        }
    }

    if total_possible == 0 {
        0.0
    } else {
        n_samples as f64 / total_possible as f64
    }
}

/// Sample sorted start indices within a tile core using an approximate density.
///
/// The helper caps starts so the maximum fragment length fits inside the chromosome, then draws
/// `round(density * possible)` unique positions without replacement. When the estimate exceeds
/// the available positions, every valid start in the tile is returned.
pub(crate) fn sample_starts_in_core<R: Rng + ?Sized>(
    rng: &mut R,
    core_start: u64,
    core_end: u64,
    chrom_len: u64,
    max_fragment_length: u64,
    sampling_density: f64,
) -> Vec<usize> {
    if core_start >= core_end || chrom_len < max_fragment_length || max_fragment_length == 0 {
        return Vec::new();
    }

    // Restrict to starts that keep the maximum fragment length within the chromosome
    let max_start_exclusive = chrom_len - max_fragment_length + 1;
    let range_end = core_end.min(max_start_exclusive);
    if range_end <= core_start {
        return Vec::new();
    }

    let possible_in_tile = (range_end - core_start) as usize;
    let estimated = (sampling_density * possible_in_tile as f64).round() as usize;
    let n_to_sample = estimated.min(possible_in_tile);
    if n_to_sample == 0 {
        return Vec::new();
    }

    let idxs = rand::seq::index::sample(rng, possible_in_tile, n_to_sample);
    let mut starts: Vec<usize> = idxs
        .into_iter()
        .map(|offset| offset + core_start as usize)
        .collect();
    starts.sort_unstable();
    starts
}

#[cfg(test)]
mod tests {
    include!("sampling_tests.rs");
}
