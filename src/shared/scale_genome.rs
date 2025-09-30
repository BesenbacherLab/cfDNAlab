use crate::shared::{bam::Contigs, overlaps::OverlappingWindows};
use anyhow::{Context, Result, bail};
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

/// Compute per-window scaling factors averaged over the **overlapped span only**.
///
/// Returns
/// -------
/// For each count window that overlaps the fragment, returns:
///     `(window_idx, avg_scaling_over_overlap, overlap_fraction)`
///     where
///         `avg_scaling_over_overlap` is the average per-base
///         scaling evaluated over `overlap = window ∩ fragment`.
#[inline]
pub fn compute_window_scaling_over_overlap(
    count_overlaps: &OverlappingWindows,
    scaling_bin_indices: &[usize],
    scaling_chr: &[(u64, u64, f32)],
) -> Result<Vec<(usize, f64, f64)>> {
    let fragment_start_bp = count_overlaps.interval_start;
    let fragment_end_bp = count_overlaps.interval_end;

    if fragment_end_bp <= fragment_start_bp {
        bail!("count_overlaps.interval_start >= interval_end (empty fragment span)");
    }
    if scaling_bin_indices.is_empty() {
        bail!("scaling_bin_indices is empty but scaling was requested");
    }

    let mut per_window_scaling = Vec::with_capacity(count_overlaps.windows.len());

    for window in &count_overlaps.windows {
        let window_start_bp = window.win_start;
        let window_end_bp = window.win_end;

        let overlap_start_bp = fragment_start_bp.max(window_start_bp);
        let overlap_end_bp = fragment_end_bp.min(window_end_bp);
        if overlap_end_bp <= overlap_start_bp {
            continue; // No overlap with this window
        }

        let avg_scaling = avg_scaling_over_span(
            overlap_start_bp,
            overlap_end_bp,
            scaling_bin_indices,
            scaling_chr,
        )?;
        per_window_scaling.push((window.idx, avg_scaling, window.overlap_fraction as f64));
    }

    Ok(per_window_scaling)
}

/// Compute per-window scaling factors averaged over the **entire fragment** (treat fragment as fully included).
///
/// Returns
/// -------
/// For each count window that overlaps the fragment, returns:
/// `(window_idx, avg_scaling_over_fragment, full_overlap_fraction)` where `avg_scaling_over_fragment` is the average per-base
/// scaling evaluated over the **whole** fragment span `[fragment_start_bp, fragment_end_bp)`.
/// This value is identical for every returned window of the same fragment.
/// `full_overlap_fraction` is always `1.0`.
#[inline]
pub fn compute_window_scaling_over_fragment(
    count_overlaps: &OverlappingWindows,
    scaling_bin_indices: &[usize],
    scaling_chr: &[(u64, u64, f32)],
) -> Result<Vec<(usize, f64, f64)>> {
    let fragment_start_bp = count_overlaps.interval_start;
    let fragment_end_bp = count_overlaps.interval_end;

    if fragment_end_bp <= fragment_start_bp {
        bail!("count_overlaps.interval_start >= interval_end (empty fragment span)");
    }
    if scaling_bin_indices.is_empty() {
        bail!("scaling_bin_indices is empty but scaling was requested");
    }

    // Compute one average over the full fragment span.
    let avg_over_fragment = avg_scaling_over_span(
        fragment_start_bp,
        fragment_end_bp,
        scaling_bin_indices,
        scaling_chr,
    )?;

    // Emit the same value for every window that actually overlaps the fragment.
    let mut per_window_scaling = Vec::with_capacity(count_overlaps.windows.len());
    for window in &count_overlaps.windows {
        if window.win_end > fragment_start_bp && window.win_start < fragment_end_bp {
            per_window_scaling.push((window.idx, avg_over_fragment, 1.0));
        }
    }
    Ok(per_window_scaling)
}

/// Average per-base scaling over an arbitrary span `[span_start_bp, span_end_bp)`.
///
/// - Uses `scaling_bin_indices` (sorted, non-empty) to limit work to bins that touch the fragment.
/// - `scaling_chr` entries are `(bin_start_bp, bin_end_bp, weight_per_base)`.
/// - `weight_per_base` is already the factor you want to multiply (e.g., 1.0 / normalized_coverage).
#[inline]
fn avg_scaling_over_span(
    span_start_bp: u64,
    span_end_bp: u64,
    scaling_bin_indices: &[usize],
    scaling_chr: &[(u64, u64, f32)],
) -> Result<f64> {
    if span_end_bp <= span_start_bp {
        bail!("avg_scaling_over_span called with empty or inverted span");
    }
    if scaling_bin_indices.is_empty() {
        bail!("scaling_bin_indices is empty but scaling was requested");
    }

    let span_len_bp = (span_end_bp - span_start_bp) as f64;
    let mut weighted_sum_bp = 0.0f64;

    // Walk only the intersecting scaling bins and accumulate:
    // weighted_sum_bp += overlap_len_bp * weight_per_base
    for &bin_idx in scaling_bin_indices {
        let (bin_start_bp, bin_end_bp, weight_per_base) = scaling_chr[bin_idx];

        // Skip bins entirely left of the span
        if bin_end_bp <= span_start_bp {
            continue;
        }
        // Stop once bins start at or beyond the span end (bins must be sorted)
        if bin_start_bp >= span_end_bp {
            break;
        }

        let overlap_start_bp = span_start_bp.max(bin_start_bp);
        let overlap_end_bp = span_end_bp.min(bin_end_bp);
        if overlap_end_bp > overlap_start_bp {
            let overlap_len_bp = (overlap_end_bp - overlap_start_bp) as f64;
            weighted_sum_bp += overlap_len_bp * (weight_per_base as f64);
        }
    }

    Ok(weighted_sum_bp / span_len_bp)
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
