use anyhow::Result;
use fxhash::FxHashMap;
use rand::Rng;
use std::cmp::Ordering;

/// Estimate the sampling density given a desired number of positions and the maximum window length.
///
/// The density is the ratio `n_samples / total_possible_starts` where `total_possible_starts`
/// counts all valid start indices across chromosomes (`len - max_window_len + 1` per chrom).
/// Returns `0.0` when no valid starts exist or when `n_samples` is zero.
pub fn sampling_density(
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
pub fn sample_starts_in_core<R: Rng + ?Sized>(
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

/// Sample unique start indices proportionally per chromosome, weighted by the number
/// of **valid** window starts (i.e., `chrom_len - max_window_len + 1`). Sampling is **without
/// replacement** within each chromosome and the resulting start indices are **sorted**.
///
/// Allocation of the per-chromosome quotas uses **Hamilton apportionment**:
/// floor each proportional share, then distribute the remaining samples to the chromosomes
/// with the largest fractional remainders. This ensures the final total equals `n_samples`.
///
/// Chromosomes with `chrom_len < max_window_len` receive **zero** samples and are **omitted** from
/// the output map.
///
/// Parameters
/// ----------
/// rng: Rng + ?Sized
///     Random number generator. Accepts concrete RNGs (e.g., `ThreadRng`, `StdRng`)
///     or trait objects via `&mut dyn RngCore`.
/// chrom_sizes: FxHashMap<String, usize>
///     Map from chromosome name -> chromosome length in bp.
/// n_samples: usize
///     Total number of start positions to sample across all chromosomes.
/// max_window_len: usize
///     Maximum window length you will use downstream. Valid start positions on a
///     chromosome of length `L` are in the inclusive range `0..=(L - max_window_len)`,
///     so the count is `L - max_window_len + 1` when `L >= max_window_len`, else `0`.
///
/// Returns
/// -------
/// starts_by_chrom: FxHashMap<String, Vec<usize>>
///     Per-chromosome sampled start positions. Each vector is **sorted ascending** and
///     contains **no duplicates**. Chromosomes with zero quota are **absent**.
pub fn sample_starts_per_chrom<R: rand::Rng + ?Sized>(
    rng: &mut R,
    chrom_sizes: &FxHashMap<String, usize>, // chrom -> length (bp)
    n_samples: usize,
    max_window_len: usize,
) -> Result<FxHashMap<String, Vec<usize>>> {
    // Compute valid start counts per chromosome and total
    let per_chrom_possible: Vec<(&str, usize)> = chrom_sizes
        .iter()
        .map(|(name, &len)| {
            // Number of valid starts for max window length: 0..=len-max_len  =>  (len - max_len + 1)
            let possible = if len >= max_window_len {
                len - max_window_len + 1
            } else {
                0
            };
            (name.as_str(), possible)
        })
        .filter(|&(_, possible)| possible > 0)
        .collect();

    let total_possible: usize = per_chrom_possible.iter().map(|&(_, p)| p).sum();
    if total_possible == 0 || n_samples == 0 {
        return Ok(FxHashMap::default());
    }

    // Hamilton apportionment: Floor + largest remainders to hit exactly n_samples
    let mut quotas: Vec<(&str, usize, f64)> = per_chrom_possible
        .iter()
        .map(|&(name, possible)| {
            let exact = (possible as f64) * (n_samples as f64) / (total_possible as f64);
            let base = exact.floor() as usize;
            (name, base, exact - base as f64)
        })
        .collect();

    let assigned_base: usize = quotas.iter().map(|&(_, b, _)| b).sum();
    let remaining = n_samples.saturating_sub(assigned_base);

    // Give remaining by largest fractional remainder
    quotas.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(Ordering::Equal));
    for q in quotas.iter_mut().take(remaining) {
        q.1 += 1;
    }

    // Build a quick lookup for possible starts per chrom
    let possible_map: FxHashMap<&str, usize> = per_chrom_possible.into_iter().collect();

    // Sample per chromosome (without replacement) and sort
    let mut out: FxHashMap<String, Vec<usize>> =
        FxHashMap::with_capacity_and_hasher(quotas.len(), Default::default());

    for (name, quota, _) in quotas {
        if quota == 0 {
            continue;
        }
        let possible = possible_map[&name];
        // Sample `quota` unique start positions from [0, possible)
        let idxs = rand::seq::index::sample(rng, possible, quota);
        let mut starts: Vec<usize> = idxs.into_iter().collect();
        starts.sort_unstable();
        out.insert(name.to_string(), starts);
    }

    Ok(out)
}
