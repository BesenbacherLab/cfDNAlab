use anyhow::Result;
use fxhash::FxHashMap;
use std::cmp::Ordering;

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
