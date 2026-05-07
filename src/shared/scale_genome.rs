use crate::shared::{
    bam::Contigs,
    interval::Interval,
    overlaps::{OverlappingWindow, OverlappingWindows},
};
use anyhow::{Context, Result, bail, ensure};
use fxhash::{FxHashMap, FxHashSet};
use std::io::{BufRead, BufReader};
use std::path::Path;
use tracing::warn;

/// One bin in a per-chromosome scaling factor table.
///
/// Bins must be sorted, non-overlapping, and contiguous (validated on load).
/// `weight_per_base` is a multiplicative factor applied to each base of coverage
/// that falls within the bin (e.g., `1.0 / normalised_coverage`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScalingBin {
    pub interval: Interval<u64>,
    pub weight_per_base: f32,
}

impl ScalingBin {
    pub fn new(start: u64, end: u64, weight_per_base: f32) -> anyhow::Result<Self> {
        Ok(Self {
            interval: Interval::new(start, end)?,
            weight_per_base,
        })
    }

    /// Unpack the bin as `(start, end, weight_per_base)`.
    #[inline]
    pub fn as_tuple(&self) -> (u64, u64, f32) {
        let (start, end) = self.interval.as_tuple();
        (start, end, self.weight_per_base)
    }
}

/// Scaling weight assigned to one selected count window.
///
/// `window_interval` is the original count or assignment window. `scaling_interval`
/// is the aligned-reference span used to average `scaling_weight`; these differ
/// when clipped-only assignment windows are remapped to the nearest aligned base.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WindowScaling {
    pub window_idx: usize,
    pub window_interval: Interval<u64>,
    pub scaling_interval: Interval<u64>,
    pub scaling_weight: f64,
    pub overlap_fraction_to_count: f64,
}

/// GC correction mode for how scaling factors were built.
///
/// This is only used for compatibility checks. It must never affect how scaling
/// factors are applied numerically.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum ScalingGCMode {
    /// No metadata was present, so we trust the user and skip compatibility checks.
    #[default]
    Unknown,
    /// Scaling factors were built from uncorrected coverage.
    Uncorrected,
    /// Scaling factors were built with `--gc-file`.
    CorrectedFromFile,
    /// Scaling factors were built with `--gc-tag`.
    CorrectedFromTag,
}

impl ScalingGCMode {
    #[inline]
    pub fn as_metadata_value(self) -> &'static str {
        match self {
            ScalingGCMode::Unknown => "unknown",
            ScalingGCMode::Uncorrected => "uncorrected",
            ScalingGCMode::CorrectedFromFile => "corrected_file",
            ScalingGCMode::CorrectedFromTag => "corrected_tag",
        }
    }

    #[inline]
    pub fn describe(self) -> &'static str {
        match self {
            ScalingGCMode::Unknown => "unknown GC correction mode",
            ScalingGCMode::Uncorrected => "no GC correction",
            ScalingGCMode::CorrectedFromFile => "GC-corrected coverage via --gc-file",
            ScalingGCMode::CorrectedFromTag => "GC-corrected coverage via --gc-tag",
        }
    }
}

/// Metadata embedded in a scaling-factors TSV.
///
/// Files may include leading comment lines with `# key=value` metadata before
/// the header row. Unknown keys are ignored to keep the format extensible.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ScalingFactorsMetadata {
    /// Whether the source coverage used GC correction.
    pub gc_mode: ScalingGCMode,
    /// Whether coverage-based smoothing weights were built with `--ignore-gap`.
    ///
    /// Old files and count-based scaling files may omit this metadata.
    pub ignore_gap: Option<bool>,
}

/// Parsed scaling factors together with their file metadata.
#[derive(Debug, Clone)]
pub struct LoadedScalingFactors {
    pub bins_by_chromosome: FxHashMap<String, Vec<ScalingBin>>,
    pub metadata: ScalingFactorsMetadata,
}

#[inline]
pub fn scaling_gc_mode_for_run(gc_file_enabled: bool, gc_tag_enabled: bool) -> ScalingGCMode {
    if gc_file_enabled {
        ScalingGCMode::CorrectedFromFile
    } else if gc_tag_enabled {
        ScalingGCMode::CorrectedFromTag
    } else {
        ScalingGCMode::Uncorrected
    }
}

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
///     Per-chromosome scaling bins with full coverage.
///
/// Returns
/// -------
/// - _:
///     Scales `cov[a..b]` by `1*sf` for each overlapping bin.
#[inline]
pub fn apply_scaling_to_coverage_in_place(
    cov: &mut [f32],
    core_start: u32,
    bins: &[ScalingBin],
) {
    if cov.is_empty() || bins.is_empty() {
        return;
    }
    let start_abs = core_start as u64;
    let end_abs = start_abs + cov.len() as u64;

    // Find first bin whose end > start_abs (upper-bound on `end`).
    // partition_point returns the first index where the condition is false.
    let mut i = bins.partition_point(|b| b.interval.end() <= start_abs);

    // Linear sweep over bins until we pass end_abs
    while i < bins.len() {
        let bin_start = bins[i].interval.start();
        let bin_end = bins[i].interval.end();
        let scaling_factor = bins[i].weight_per_base;
        if bin_start >= end_abs {
            break;
        }
        let overlap_start = bin_start.max(start_abs);
        let overlap_end = bin_end.min(end_abs);
        if overlap_end > overlap_start {
            let slice_start = (overlap_start - start_abs) as usize;
            let slice_end = (overlap_end - start_abs) as usize;
            for v in &mut cov[slice_start..slice_end] {
                *v *= scaling_factor;
            }
        }
        if bin_end >= end_abs {
            break;
        }
        i += 1;
    }
}

/// Compute per-window scaling factors averaged over the **overlapped span only**.
///
/// By default, scaling is averaged over `count_overlaps`. Pass `scaling_overlaps`
/// only when the count-window coordinates and the scaling coordinates intentionally differ.
/// The two overlap collections must then have the same rows in the same order.
///
/// Returns
/// -------
/// For each count window that overlaps the fragment, returns:
///     `WindowScaling { window_idx, window_interval, scaling_interval, scaling_weight, overlap_fraction_to_count }`
///     where
///         `scaling_weight` is the average per-base
///         scaling evaluated over the overlap between the scaling interval and
///         the query interval used for scaling.
#[inline]
pub fn compute_per_window_scaling_over_overlap(
    count_overlaps: &OverlappingWindows,
    scaling_overlaps: Option<&OverlappingWindows>,
    scaling_bin_indices: &[usize],
    scaling_chr: &[ScalingBin],
) -> Result<Vec<WindowScaling>> {
    let scaling_overlaps = scaling_overlaps.unwrap_or(count_overlaps);
    ensure!(
        count_overlaps.windows.len() == scaling_overlaps.windows.len(),
        "count and scaling overlap row counts differ: {} vs {}",
        count_overlaps.windows.len(),
        scaling_overlaps.windows.len()
    );
    let scaling_query_start_bp = scaling_overlaps.query_start();
    let scaling_query_end_bp = scaling_overlaps.query_end();

    let mut per_window_scaling = Vec::with_capacity(count_overlaps.windows.len());

    for (count_window, scaling_window) in count_overlaps
        .windows
        .iter()
        .zip(scaling_overlaps.windows.iter())
    {
        ensure!(
            count_window.idx == scaling_window.idx,
            "count and scaling overlap window indices differ: {} vs {}",
            count_window.idx,
            scaling_window.idx
        );
        ensure!(
            count_window.overlap_fraction == scaling_window.overlap_fraction,
            "count and scaling overlap fractions differ for window {}: {} vs {}",
            count_window.idx,
            count_window.overlap_fraction,
            scaling_window.overlap_fraction
        );

        let scaling_window_start_bp = scaling_window.start();
        let scaling_window_end_bp = scaling_window.end();

        let overlap_start_bp = scaling_query_start_bp.max(scaling_window_start_bp);
        let overlap_end_bp = scaling_query_end_bp.min(scaling_window_end_bp);
        if overlap_end_bp <= overlap_start_bp {
            continue; // No overlap with this window
        }

        let avg_scaling = avg_scaling_over_span(
            overlap_start_bp,
            overlap_end_bp,
            scaling_bin_indices,
            scaling_chr,
        )?;
        let scaling_interval = Interval::new(overlap_start_bp, overlap_end_bp)?;
        per_window_scaling.push(WindowScaling {
            window_idx: count_window.idx,
            window_interval: count_window.interval,
            scaling_interval,
            scaling_weight: avg_scaling,
            overlap_fraction_to_count: count_window.overlap_fraction,
        });
    }

    Ok(per_window_scaling)
}

/// Compute per-window scaling factors averaged over the **entire fragment** (treat fragment as fully included).
///
/// Returns
/// -------
/// For each count window that overlaps the fragment, returns:
/// `WindowScaling { window_idx, window_interval, scaling_interval, scaling_weight, overlap_fraction_to_count }`
/// where `scaling_weight` is the average per-base
/// scaling evaluated over the **whole** fragment span `[fragment_start_bp, fragment_end_bp)`.
/// This value is identical for every returned window of the same fragment.
/// `overlap_fraction_to_count` is always `1.0`.
#[inline]
pub fn compute_per_window_scaling_over_fragment(
    fragment_interval: Interval<u64>,
    count_overlaps: &OverlappingWindows,
    scaling_bin_indices: &[usize],
    scaling_chr: &[ScalingBin],
) -> Result<Vec<WindowScaling>> {
    let fragment_start_bp = fragment_interval.start();
    let fragment_end_bp = fragment_interval.end();
    let avg_over_fragment = avg_scaling_over_span(
        fragment_start_bp,
        fragment_end_bp,
        scaling_bin_indices,
        scaling_chr,
    )?;

    // Emit the same value for every window that actually overlaps the fragment.
    let mut per_window_scaling = Vec::with_capacity(count_overlaps.windows.len());
    for window in &count_overlaps.windows {
        if window.end() > fragment_start_bp && window.start() < fragment_end_bp {
            per_window_scaling.push(WindowScaling {
                window_idx: window.idx,
                window_interval: window.interval,
                scaling_interval: fragment_interval,
                scaling_weight: avg_over_fragment,
                overlap_fraction_to_count: 1.0,
            });
        }
    }
    Ok(per_window_scaling)
}

/// Compute one fragment-level scaling factor and apply it to every selected count window.
///
/// Use this when a command has already selected candidate windows with assignment coordinates that
/// may differ from the aligned fragment span used for scaling. The scaling average is still
/// computed over `fragment_interval`, but every row in `count_overlaps.windows` is returned.
#[inline]
pub fn compute_per_window_scaling_over_fragment_for_selected_windows(
    fragment_interval: Interval<u64>,
    count_overlaps: &OverlappingWindows,
    scaling_bin_indices: &[usize],
    scaling_chr: &[ScalingBin],
) -> Result<Vec<WindowScaling>> {
    let avg_over_fragment = avg_scaling_over_span(
        fragment_interval.start(),
        fragment_interval.end(),
        scaling_bin_indices,
        scaling_chr,
    )?;

    let mut per_window_scaling = Vec::with_capacity(count_overlaps.windows.len());
    for window in &count_overlaps.windows {
        per_window_scaling.push(WindowScaling {
            window_idx: window.idx,
            window_interval: window.interval,
            scaling_interval: fragment_interval,
            scaling_weight: avg_over_fragment,
            overlap_fraction_to_count: 1.0,
        });
    }
    Ok(per_window_scaling)
}

/// Build a reference-based scaling view for assignment-selected overlap rows.
///
/// The returned rows keep each assignment window's `idx` and overlap fraction, but replace the
/// interval with the aligned-reference span used to average scaling. Windows with no aligned
/// overlap use the nearest aligned reference base.
pub fn build_reference_based_scaling_overlaps_for_assignment_overlaps(
    count_overlaps: &OverlappingWindows,
    aligned_fragment_interval: Interval<u64>,
) -> Result<OverlappingWindows> {
    let left_nearest_base_interval = Interval::new(
        aligned_fragment_interval.start(),
        aligned_fragment_interval.start() + 1,
    )?;
    let right_nearest_base_interval = Interval::new(
        aligned_fragment_interval.end() - 1,
        aligned_fragment_interval.end(),
    )?;

    let mut scaling_overlaps = OverlappingWindows::new(aligned_fragment_interval);
    for window in &count_overlaps.windows {
        let scaling_interval = if let Some(aligned_overlap_interval) =
            window.interval.clip_to(aligned_fragment_interval)
        {
            aligned_overlap_interval
        } else if window.end() <= aligned_fragment_interval.start() {
            left_nearest_base_interval
        } else if window.start() >= aligned_fragment_interval.end() {
            right_nearest_base_interval
        } else {
            bail!(
                "assignment overlap window {} has no aligned overlap but is not outside aligned interval",
                window.idx
            );
        };

        scaling_overlaps.windows.push(OverlappingWindow::new(
            window.idx,
            scaling_interval,
            window.overlap_fraction,
        )?);
    }

    Ok(scaling_overlaps)
}

/// Average per-base scaling over an arbitrary span `[span_start_bp, span_end_bp)`.
///
/// - Uses `scaling_bin_indices` (sorted, non-empty) to limit work to bins that touch the fragment.
/// - `weight_per_base` is already the factor you want to multiply (e.g., 1.0 / normalized_coverage).
#[inline]
fn avg_scaling_over_span(
    span_start_bp: u64,
    span_end_bp: u64,
    scaling_bin_indices: &[usize],
    scaling_chr: &[ScalingBin],
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
        let (bin_start_bp, bin_end_bp, weight_per_base) = scaling_chr[bin_idx].as_tuple();

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
/// - Leading comment lines starting with `#` are allowed before the header.
/// - Blank lines are not allowed before the header.
/// - Metadata comments can use `# key=value`, e.g. `# gc_mode=corrected_tag`.
/// - The TSV **must** have a header. Column names are matched **case-insensitively**.
/// - Required columns: `chromosome`, `start`, `end`, `scaling_factor`.
/// - Coordinates are 0-based, half-open `[start, end)`.
/// - `scaling_factor` must be finite and strictly >= 0.
/// - Bins are filtered to the provided `chromosomes`.
/// - For every chromosome in `chromosomes`, bins must:
///   * start at the 0-coordinate,
///   * be perfectly contiguous (no gaps, no overlaps),
///   * end exactly at that chromosome's length (from `contigs`).
///
/// Returns
/// -------
/// Parsed bins plus any metadata found in the leading comment block.
pub fn load_scaling_factors_tsv(
    path: &Path,
    chromosomes: &[String],
    contigs: &Contigs,
) -> anyhow::Result<LoadedScalingFactors> {
    let f = std::fs::File::open(path)
        .with_context(|| format!("opening scaling TSV {}", path.display()))?;
    let mut r = BufReader::new(f);

    // Read optional metadata comments, then parse the required header.
    // We lower-case the final header once to allow case-insensitive matching,
    // but keep the original names for error messages.
    let mut header = String::new();
    let mut metadata = ScalingFactorsMetadata::default();
    let mut lineno = 0usize;
    let mut saw_gc_mode_metadata = false;
    let mut saw_ignore_gap_metadata = false;
    loop {
        header.clear();
        if r.read_line(&mut header)? == 0 {
            anyhow::bail!("{}: empty file; header required", path.display());
        }
        lineno += 1;

        let trimmed = header.trim();
        if trimmed.is_empty() {
            anyhow::bail!(
                "{}:{}: blank lines are not allowed before the scaling TSV header",
                path.display(),
                lineno
            );
        }
        if let Some(comment_body) = trimmed.strip_prefix('#') {
            parse_scaling_metadata_comment(
                comment_body.trim(),
                path,
                lineno,
                &mut metadata,
                &mut saw_gc_mode_metadata,
                &mut saw_ignore_gap_metadata,
            )?;
            continue;
        }
        break;
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
    let mut map: FxHashMap<String, Vec<ScalingBin>> =
        FxHashMap::with_hasher(Default::default());

    // Stream rows; reuse `line` to avoid allocations per row
    let mut line = String::new();
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
        let start: u64 = fields[s_i].parse().with_context(|| {
            format!(
                "{}:{}: invalid start '{}'",
                path.display(),
                lineno,
                fields[s_i]
            )
        })?;
        let end: u64 = fields[e_i].parse().with_context(|| {
            format!(
                "{}:{}: invalid end '{}'",
                path.display(),
                lineno,
                fields[e_i]
            )
        })?;
        if start >= end {
            anyhow::bail!(
                "{}:{}: invalid interval [{}..{})",
                path.display(),
                lineno,
                start,
                end
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

        // Stash; we'll sort and validate contiguity/full coverage per chromosome below
        map.entry(chr.to_string())
            .or_default()
            .push(ScalingBin::new(start, end, sf)?);
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
        v.sort_unstable_by_key(|b| b.interval.start());

        // Chromosome length from the BAM-derived contig map.
        let chrom_len = contigs
            .contigs
            .get(chr)
            .map(|&(_, len)| len as u64)
            .ok_or_else(|| anyhow::anyhow!("missing contig info for '{}'", chr))?;

        // Must start at 0; `unwrap_or(1)` ensures a clean error if somehow empty
        if v.first().map(|b| b.interval.start()).unwrap_or(1) != 0 {
            anyhow::bail!("scaling TSV: bins on '{}' must start at 0", chr);
        }

        // Sweep to ensure each bin begins exactly where the previous ended.
        // (Positive length is already enforced by ScalingBin::new / Interval::new.)
        let mut prev_end = 0u64;
        for b in v.iter() {
            let (s, e) = b.interval.as_tuple();
            if s != prev_end {
                anyhow::bail!(
                    "scaling TSV: bins on '{}' are not contiguous at {}..{} (prev_end={})",
                    chr,
                    s,
                    e,
                    prev_end
                );
            }
            prev_end = e;
        }

        // Require exact coverage up to the chromosome's length
        if prev_end != chrom_len {
            anyhow::bail!(
                "scaling TSV: bins on '{}' must end at chrom_len={} (got end={})",
                chr,
                chrom_len,
                prev_end
            );
        }
    }

    Ok(LoadedScalingFactors {
        bins_by_chromosome: map,
        metadata,
    })
}

fn parse_scaling_metadata_comment(
    comment_body: &str,
    path: &Path,
    lineno: usize,
    metadata: &mut ScalingFactorsMetadata,
    saw_gc_mode_metadata: &mut bool,
    saw_ignore_gap_metadata: &mut bool,
) -> Result<()> {
    let Some((raw_key, raw_value)) = comment_body.split_once('=') else {
        return Ok(());
    };

    let key = raw_key.trim().to_ascii_lowercase();
    let value = raw_value.trim();
    match key.as_str() {
        "gc_mode" => {
            if *saw_gc_mode_metadata {
                bail!(
                    "{}:{}: duplicate scaling metadata key 'gc_mode'",
                    path.display(),
                    lineno
                );
            }
            metadata.gc_mode = parse_scaling_gc_mode(path, value)?;
            *saw_gc_mode_metadata = true;
        }
        "ignore_gap" => {
            if *saw_ignore_gap_metadata {
                bail!(
                    "{}:{}: duplicate scaling metadata key 'ignore_gap'",
                    path.display(),
                    lineno
                );
            }
            metadata.ignore_gap = Some(parse_scaling_bool(path, "ignore_gap", value)?);
            *saw_ignore_gap_metadata = true;
        }
        _ => {}
    }
    Ok(())
}

fn parse_scaling_gc_mode(path: &Path, value: &str) -> Result<ScalingGCMode> {
    match value.to_ascii_lowercase().as_str() {
        "unknown" => Ok(ScalingGCMode::Unknown),
        "uncorrected" => Ok(ScalingGCMode::Uncorrected),
        "corrected_file" => Ok(ScalingGCMode::CorrectedFromFile),
        "corrected_tag" => Ok(ScalingGCMode::CorrectedFromTag),
        _ => bail!(
            "{}: invalid value '{}' for scaling metadata key 'gc_mode'; expected one of 'unknown', 'uncorrected', 'corrected_file', or 'corrected_tag'",
            path.display(),
            value,
        ),
    }
}

fn parse_scaling_bool(path: &Path, key: &str, value: &str) -> Result<bool> {
    match value.to_ascii_lowercase().as_str() {
        "true" => Ok(true),
        "false" => Ok(false),
        _ => bail!(
            "{}: invalid value '{}' for scaling metadata key '{}'; expected 'true' or 'false'",
            path.display(),
            value,
            key,
        ),
    }
}

/// Error on known scaling GC-mode mismatches.
pub fn ensure_scaling_gc_compatibility(
    path: &Path,
    scaling_factor_metadata: ScalingFactorsMetadata,
    command_gc_mode: ScalingGCMode,
) -> Result<()> {
    match (scaling_factor_metadata.gc_mode, command_gc_mode) {
        // Unknown file metadata means we do not know how the scaling file was built.
        // Warn so users know the check was skipped, then trust them to decide.
        (ScalingGCMode::Unknown, _) => {
            warn!(
                target: "scale-genome",
                "warning: scaling factors file {} has no gc_mode metadata, so GC compatibility could not be checked",
                path.display()
            );
            Ok(())
        }

        // Known uncorrected/corrected mismatches should fail.
        (
            ScalingGCMode::Uncorrected,
            ScalingGCMode::CorrectedFromFile | ScalingGCMode::CorrectedFromTag,
        )
        | (
            ScalingGCMode::CorrectedFromFile | ScalingGCMode::CorrectedFromTag,
            ScalingGCMode::Uncorrected,
        ) => bail!(
            "Scaling factors file {} is incompatible with this run: file GC mode is '{}', but current run uses '{}'. Use scaling factors built with matching GC-correction status for this analysis.",
            path.display(),
            scaling_factor_metadata.gc_mode.describe(),
            command_gc_mode.describe(),
        ),

        // All other known combinations are acceptable:
        // - uncorrected with uncorrected
        // - file with file
        // - tag with tag
        // - file with tag
        // - tag with file
        _ => Ok(()),
    }
}

/// Warn when coverage-based smoothing factors were built with different gap handling.
///
/// Missing metadata is accepted silently. Older scaling TSVs and count-based scaling TSVs do not
/// necessarily know this setting.
pub fn warn_on_scaling_ignore_gap_mismatch(
    path: &Path,
    scaling_factor_metadata: ScalingFactorsMetadata,
    current_ignore_gap: Option<bool>,
) {
    let Some(file_ignore_gap) = scaling_factor_metadata.ignore_gap else {
        return;
    };
    let Some(current_ignore_gap) = current_ignore_gap else {
        return;
    };
    if file_ignore_gap != current_ignore_gap {
        warn!(
            target: "scale-genome",
            "warning: scaling factors file {} was built with ignore_gap={}, but current run uses ignore_gap={}. This can introduce a small gap-handling mismatch in coverage normalization.",
            path.display(),
            file_ignore_gap,
            current_ignore_gap
        );
    }
}
