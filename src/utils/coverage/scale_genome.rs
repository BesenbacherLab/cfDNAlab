use crate::utils::{bam::Contigs, overlaps::OverlappingWindows};
use anyhow::Context;
use fxhash::{FxHashMap, FxHashSet};
use std::io::{BufRead, BufReader};

/// Apply per-bin multiplicative scaling to the tile coverage in-place.
/// Assumes bins are sorted, non-overlapping, and fully cover the chromosome.
///
/// Parameters
/// ----------
/// - cov:
///     Mutable coverage slice for the tile core.
/// - core_start:
///     Absolute start of the tile core (0-based).
/// - bins:
///     Per-chromosome scaling bins `(start, end, sf)` with full coverage.
///
/// Returns
/// -------
/// - _:
///     Scales `cov[a..b]` by `1*sf` for each overlapping bin.
#[inline]
pub fn apply_scaling_to_coverage_in_place(
    cov: &mut [f32],
    core_start: u32,
    bins: &[(u64, u64, f32)],
) {
    if cov.is_empty() || bins.is_empty() {
        return;
    }
    let start_abs = core_start as u64;
    let end_abs = start_abs + cov.len() as u64;

    // Find first bin with end > start_abs (upper_bound on `end`)
    // First bin with end > start_abs (upper-bound on `end`)
    // Note: partition_point gets first element (from left-right) where the condition is false
    let mut i = bins.partition_point(|t| t.1 <= start_abs);

    // Linear sweep over bins until we pass end_abs
    while i < bins.len() {
        let (bs, be, sf) = bins[i];
        if bs >= end_abs {
            break;
        }
        let s = bs.max(start_abs); // Overlap start
        let e = be.min(end_abs); // Overlap end
        if e > s {
            let a = (s - start_abs) as usize; // Slice start in cov
            let b = (e - start_abs) as usize; // Slice end in cov
            for v in &mut cov[a..b] {
                *v *= sf; // Multiply by scaling factor
            }
        }
        if be >= end_abs {
            break; // Finished the tile
        }
        i += 1;
    }
}

/// Compute mass-conserving, scaling-aware weights for one fragment across overlapped count windows.
///
/// For each count window that overlaps the fragment, compute the fragment’s proportional
/// contribution and modulate it by the **average per-base scaling weight** over that *overlapped
/// span*. Here, the scaling weight stored in `scaling_chr[..].2` is already the value you want
/// to **multiply** with per base (i.e., if you conceptually “divide by coverage”, the stored
/// weight is the pre-inverted factor).
///
/// Details
/// -------
/// For fragment `[frag_start, frag_end)` and a count window `[win_start, win_end)`, the overlapped
/// span is `[overlap_start, overlap_end)`, where:
/// `overlap_start = max(frag_start, win_start)` and `overlap_end = min(frag_end, win_end)`.
///
/// Over this overlapped span, average the piecewise-constant scaling weights defined by
/// `scaling_chr` by summing `weight * overlapped_len` for each intersecting scaling bin and then
/// dividing by the total overlapped length. The final window weight is:
/// `weight = (overlap_len / fragment_len) * avg_scaling_weight_in_that_span`.
///
/// Parameters
/// ----------
/// - count_overlaps:
///     Overlaps of this fragment with the **count windows** from `find_overlapping_windows`.
///     Uses `interval_start` (i.e., fragment start), `interval_end` (i.e., fragment end),
///     and iterates `windows` with `win_start`, `win_end`, `idx`.
/// - scaling_bin_indices:
///     Sorted, non-empty indices into `scaling_chr` for scaling bins that overlap this fragment.
///     Passing only the bins that actually intersect the fragment reduces work per window.
/// - scaling_chr:
///     Slice of scaling bins as `(start, end, weight)`, half-open on end. `weight` is the
///     **already inverted per-base scaling factor to multiply** (e.g., `1.0 / coverage_norm`).
///
/// Returns
/// -------
/// - weights:
///     Vector of `(idx, weight)` per overlapped count window, where `idx` is the count-window
///     scan index and `weight` is the mass-conserving, scaling-aware contribution.
///
/// Errors
/// ------
/// - Returns an error if `count_overlaps.interval_start >= count_overlaps.interval_end`.
/// - Returns an error if `scaling_bin_indices` is empty.
#[inline]
pub fn compute_scaled_window_weights(
    count_overlaps: &OverlappingWindows,
    scaling_bin_indices: &[usize],
    scaling_chr: &[(u64, u64, f32)],
) -> anyhow::Result<Vec<(usize, f64)>> {
    let frag_start = count_overlaps.interval_start;
    let frag_end = count_overlaps.interval_end;

    if frag_end <= frag_start {
        anyhow::bail!(
            "count_overlaps.interval_start >= count_overlaps.interval_end. This should never happen, report please."
        );
    }
    if scaling_bin_indices.is_empty() {
        anyhow::bail!("scaling_bin_indices is empty but scaling was requested");
    }

    let fragment_len = (frag_end - frag_start) as f64;
    let mut weights = Vec::with_capacity(count_overlaps.windows.len());

    for w in &count_overlaps.windows {
        let win_start = w.win_start;
        let win_end = w.win_end;

        let overlap_start = frag_start.max(win_start);
        let overlap_end = frag_end.min(win_end);
        if overlap_end <= overlap_start {
            continue;
        }

        let overlap_len = (overlap_end - overlap_start) as f64;

        // Sum weight * overlapped_len across scaling bins intersecting [overlap_start, overlap_end)
        let mut weighted_bp_sum = 0.0f64;

        for &bin_idx in scaling_bin_indices {
            let (bin_start, bin_end, weight_per_base) = scaling_chr[bin_idx];

            if bin_end <= overlap_start {
                continue; // Bin entirely left of the overlap
            }
            if bin_start >= overlap_end {
                break; // Bins are sorted; we are past the overlap
            }

            let seg_start = overlap_start.max(bin_start);
            let seg_end = overlap_end.min(bin_end);
            if seg_end > seg_start {
                let seg_len = (seg_end - seg_start) as f64;
                weighted_bp_sum += seg_len * (weight_per_base as f64); // Multiply: weights are already inverted
            }
        }

        let avg_scaling_weight = weighted_bp_sum / overlap_len;
        let weight = (overlap_len / fragment_len) * avg_scaling_weight;

        if weight > 0.0 {
            weights.push((w.idx, weight));
        }
    }

    Ok(weights)
}

/// Load per-bin scaling factors from a TSV and validate that, for each requested
/// chromosome, the bins are **sorted**, **non-overlapping**, **contiguous**, and
/// provide **full coverage** from `0` up to the chromosome length taken from `contigs`.
///
/// Requirements
/// -----------
/// - The TSV **must** have a header. Column names are matched **case-insensitively**.
/// - Required columns: `chromosome`, `start`, `end`, `scaling_factor`.
/// - Coordinates are 0-based, half-open `[start, end)`.
/// - `scaling_factor` must be finite and strictly >= 0.
/// - Bins are filtered to the provided `chromosomes`.
/// - For every chromosome in `chromosomes`, bins must:
///   * start at 0,
///   * be perfectly contiguous (no gaps, no overlaps),
///   * end exactly at that chromosome’s length (from `contigs`).
///
/// Returns
/// -------
/// A map `chr -> Vec<(start, end, scaling_factor)>`, with entries **sorted by start**.
pub fn load_scaling_factors_tsv(
    path: &std::path::Path,
    chromosomes: &[String],
    contigs: &Contigs,
) -> anyhow::Result<FxHashMap<String, Vec<(u64, u64, f32)>>> {
    let f = std::fs::File::open(path)
        .with_context(|| format!("opening scaling TSV {}", path.display()))?;
    let mut r = BufReader::new(f);

    // Read and parse the header (required). We lower-case once to allow
    // case-insensitive matching, but keep the original names for error messages.
    let mut header = String::new();
    if r.read_line(&mut header)? == 0 {
        anyhow::bail!("{}: empty file; header required", path.display());
    }
    let cols: Vec<&str> = header.trim_end_matches('\n').split('\t').collect();
    let cols_lc: Vec<String> = cols.iter().map(|c| c.to_ascii_lowercase()).collect();
    let find = |name: &str| -> anyhow::Result<usize> {
        cols_lc.iter().position(|c| c == name).ok_or_else(|| {
            anyhow::anyhow!(
                "required column '{}' not found in header: {}",
                name,
                cols.join("\t")
            )
        })
    };
    // Column indices for required fields (case-insensitive lookup)
    let chr_i = find("chromosome")?;
    let s_i = find("start")?;
    let e_i = find("end")?;
    let sf_i = find("scaling_factor")?;

    // Find max column index
    let max_i = chr_i.max(s_i).max(e_i).max(sf_i);

    // Limit to the chromosomes the user requested (fast membership checks)
    let want: FxHashSet<&str> = chromosomes.iter().map(|s| s.as_str()).collect();

    // Accumulator: per-chromosome bins (unsorted initially; sorted below)
    let mut map: FxHashMap<String, Vec<(u64, u64, f32)>> =
        FxHashMap::with_hasher(Default::default());

    // Stream rows; reuse `line` to avoid allocations per row
    let mut line = String::new();
    let mut lineno = 1usize; // Header already read
    while {
        line.clear();
        r.read_line(&mut line)?
    } > 0
    {
        lineno += 1;
        let raw = line.trim_end_matches('\n');

        // Skip empty lines and comments.
        if raw.is_empty() || raw.starts_with('#') {
            continue;
        }

        // Split fields once; validate column count against the rightmost used index
        let fields: Vec<&str> = raw.split('\t').collect();
        if fields.len() <= max_i {
            anyhow::bail!(
                "{}:{}: not enough columns (have {}, need {})",
                path.display(),
                lineno,
                fields.len(),
                max_i + 1
            );
        }

        // Filter by requested chromosomes early to avoid storing irrelevant rows
        let chr = fields[chr_i];
        if !want.contains(chr) {
            continue;
        }

        // Parse coordinates and scaling factor with precise error context
        let s: u64 = fields[s_i].parse().with_context(|| {
            format!(
                "{}:{}: invalid start '{}'",
                path.display(),
                lineno,
                fields[s_i]
            )
        })?;
        let e: u64 = fields[e_i].parse().with_context(|| {
            format!(
                "{}:{}: invalid end '{}'",
                path.display(),
                lineno,
                fields[e_i]
            )
        })?;
        if s >= e {
            anyhow::bail!(
                "{}:{}: invalid interval [{}..{})",
                path.display(),
                lineno,
                s,
                e
            );
        }

        let sf: f32 = fields[sf_i].parse().with_context(|| {
            format!(
                "{}:{}: invalid scaling_factor '{}'",
                path.display(),
                lineno,
                fields[sf_i]
            )
        })?;
        if !sf.is_finite() || sf < 0.0 {
            anyhow::bail!(
                "{}:{}: scaling_factor must be finite and >= 0 (got {})",
                path.display(),
                lineno,
                sf
            );
        }

        // Stash; we’ll sort and validate contiguity/full coverage per chromosome below
        map.entry(chr.to_string()).or_default().push((s, e, sf));
    }

    // For each requested chromosome:
    //  - sort bins by start,
    //  - verify contiguity (no gaps or overlaps),
    //  - require full coverage [0..chrom_len).
    for chr in chromosomes {
        let v = map.get_mut(chr).ok_or_else(|| {
            anyhow::anyhow!("scaling TSV: no bins provided for chromosome '{}'", chr)
        })?;

        // Sort by start once; validation below assumes non-decreasing starts
        v.sort_unstable_by_key(|t| t.0);

        // Chromosome length from the BAM-derived contig map.
        let chrom_len = contigs
            .contigs
            .get(chr)
            .map(|&(_, len)| len as u64)
            .ok_or_else(|| anyhow::anyhow!("missing contig info for '{}'", chr))?;

        // Must start at 0; `unwrap_or(1)` ensures a clean error if somehow empty
        if v.first().map(|t| t.0).unwrap_or(1) != 0 {
            anyhow::bail!("scaling TSV: bins on '{}' must start at 0", chr);
        }

        // Sweep to ensure each bin begins exactly where the previous ended,
        // and that each bin has positive length
        let mut prev_end = 0u64;
        for &(s, e, _) in v.iter() {
            if s != prev_end {
                anyhow::bail!(
                    "scaling TSV: bins on '{}' are not contiguous at {}..{} (prev_end={})",
                    chr,
                    s,
                    e,
                    prev_end
                );
            }
            if s >= e {
                anyhow::bail!(
                    "scaling TSV: invalid empty/negative bin on '{}' at {}..{}",
                    chr,
                    s,
                    e
                );
            }
            prev_end = e;
        }

        // Require exact coverage up to the chromosome’s length
        if prev_end != chrom_len {
            anyhow::bail!(
                "scaling TSV: bins on '{}' must end at chrom_len={} (got end={})",
                chr,
                chrom_len,
                prev_end
            );
        }
    }

    Ok(map)
}
