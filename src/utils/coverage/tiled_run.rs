use crate::utils::coverage::coverage_prefix::CoveragePrefix;
use crate::utils::fragment::minimal_fragment::Fragment;
use crate::utils::fragment::segment_fragment::FragmentWithSegments;
use crate::utils::{bam::Contigs, coverage::window_results::CoverageWindowAction};
use anyhow::{Context, Result};
use rand::{Rng, distr::Alphanumeric};
use std::io::{BufWriter, Write};

/// A processing tile for one chromosome
#[derive(Debug, Clone)]
pub struct Tile {
    pub chr: String,
    pub tid: i32,
    pub index: u32,       // 0-based index within chromosome
    pub core_start: u32,  // inclusive
    pub core_end: u32,    // exclusive
    pub fetch_start: u32, // inclusive (core expanded by halo)
    pub fetch_end: u32,   // exclusive
}

/// Build non-overlapping core tiles with a fetch halo on each side.
/// Optionally align cores to multiples of `align_bp` (e.g. --by-size).
///
/// Alignment rule:
/// - If `align_bp.is_some()` and `tile_bp / align_bp >= 10`, we round the tile size
///   **down** to the nearest multiple of `align_bp`. This ensures tile core boundaries
///   land exactly on bin boundaries, eliminating cross-boundary bins for most of the genome.
/// - Otherwise, we keep `tile_bp` as-is.
///
/// Parameters
/// ----------
/// - tile_bp: target core size in bases (e.g. 20_000_000)
/// - halo_bp: fetch extension on both sides (e.g. max_fragment_length)
/// - align_bp: a fixed window size that tiles should align to
///
/// Returns
/// -------
///  - tiles:
///     Vec of `Tile`s
///
///  - guaranteed_aligned:
///     Whether the window edges are guaranteed to line up with the tile edges
pub fn build_tiles(
    chromosomes: &[String],
    contigs: &Contigs,
    tile_bp: u32,
    halo_bp: u32,
    align_bp: Option<u64>,
) -> anyhow::Result<(Vec<Tile>, bool)> {
    let mut tiles = Vec::new();

    // Decide the effective core size in bases (possibly aligned)
    let (effective_tile_bp, guaranteed_aligned) = match align_bp {
        Some(bin_size) if bin_size > 0 && bin_size <= tile_bp as u64 => {
            if (tile_bp as u64) % bin_size == 0 {
                (tile_bp, true)
            } else if (tile_bp as u64) / bin_size >= 10 {
                let k = (tile_bp as u64) / bin_size;
                ((k * bin_size) as u32, true) // Never drop below one bin
            } else {
                (tile_bp as u32, false)
            }
        }
        _ => (tile_bp as u32, false),
    };

    for chr in chromosomes {
        let &(tid, chrom_len_u32) = contigs
            .contigs
            .get(chr)
            .ok_or_else(|| anyhow::anyhow!("missing contig for '{}'", chr))?;
        let chrom_len = chrom_len_u32 as u32;

        let mut start = 0u32;
        let mut idx = 0u32;
        while start < chrom_len {
            let core_end = (start.saturating_add(effective_tile_bp)).min(chrom_len);

            // Halos do not need to be aligned; they are just fetch guards.
            let fetch_start = start.saturating_sub(halo_bp);
            let fetch_end = (core_end.saturating_add(halo_bp)).min(chrom_len);

            tiles.push(Tile {
                chr: chr.clone(),
                tid: tid as i32,
                index: idx,
                core_start: start,
                core_end,
                fetch_start,
                fetch_end,
            });

            idx += 1;
            start = core_end;
        }
    }

    // Just in case we decide to move start of cores in the future
    // We'll have an extensive debug test
    #[cfg(debug_assertions)]
    if let Some(bs) = align_bp {
        if guaranteed_aligned {
            for t in &tiles {
                // Starts/ends of *cores* (not final chromosome end) line up on the grid
                debug_assert_eq!((t.core_start as u64) % bs, 0);
            }
        }
    }

    Ok((tiles, guaranteed_aligned))
}

/// Add a possibly segmented fragment into a tile-local CoveragePrefix
///
/// CoveragePrefix must be initialized to length = core_end - core_start
/// This clips each segment (or full span) to the [core_start, core_end) interval
#[inline]
pub fn add_fragment_clipped_to_core(
    cp: &mut CoveragePrefix,
    fragment: &FragmentWithSegments,
    weight: f32,
    core_start: u32,
    core_end: u32,
) -> Result<()> {
    // Use explicit segments if present
    if let Some(segments) = &fragment.segments {
        for &(seg_start_abs, seg_end_abs) in segments {
            let s = seg_start_abs.max(core_start);
            let e = seg_end_abs.min(core_end);
            if s < e {
                // Skips fragments completely outside tile
                // Shift to tile-local coordinates
                let local = Fragment {
                    tid: fragment.tid,
                    start: s - core_start,
                    end: e - core_start,
                };
                cp.add_fragment_to_prefix_weighted(local, weight)?;
            }
        }
    } else {
        // No explicit segments -> treat as one span (this already encodes your include_inter_mate_gap policy)
        let s = fragment.start.max(core_start);
        let e = fragment.end.min(core_end);
        if s < e {
            // Skips fragments completely outside tile
            // Shift to tile-local coordinates
            let local = Fragment {
                tid: fragment.tid,
                start: s - core_start,
                end: e - core_start,
            };
            cp.add_fragment_to_prefix_weighted(local, weight)?;
        }
    }
    Ok(())
}

/// What the tile should write
pub enum TileMode<'w> {
    /// Whole positional coverage for the core,
    /// or windowed positional coverage without index (unique positions)
    Positional {
        windows: Option<&'w [(u64, u64, u64)]>, // Per-chr windows if provided
        out_path: std::path::PathBuf,           // Per-tile file path
        indexed: bool,                          // Whether to save index
    },
    AggregatesByBed {
        windows: &'w [(u64, u64, u64)],    // Per-chr windows
        masked: bool,                      // Use masked counts/sums
        partials_out: std::path::PathBuf, // Cross-boundary windows (idx, sum, allowed, blacklisted)
        cross_idx_out: std::path::PathBuf, // Sidecar listing crossers
    },
    AggregatesBySize {
        window_bp: u64,                    // Fixed window size in bases
        masked: bool,                      // Use masked counts/sums
        finals_out: std::path::PathBuf,    // Final windows that do not need reducing
        partials_out: std::path::PathBuf, // Cross-boundary windows (idx, sum, allowed, blacklisted)
        cross_idx_out: std::path::PathBuf, // Sidecar listing crossers
        guaranteed_aligned: bool, // Tiles and window_bp align, so write per-tile FINALS (no reducer needed)
    },
}

/// Restrict windows to those overlapping the tile core
#[inline]
pub fn windows_overlapping_core<'a>(
    windows_chr: &'a [(u64, u64, u64)],
    core_start: u32,
    core_end: u32,
) -> impl Iterator<Item = &'a (u64, u64, u64)> {
    let cs = core_start as u64;
    let ce = core_end as u64;
    windows_chr
        .iter()
        .filter(move |&&(ws, we, _idx)| we > cs && ws < ce)
}

/// Get the tile index from the filename by scanning segments from right to left
/// and picking the first purely-numeric token (before extensions).
/// Works for:
///   coverage.pos.chr1.000123.bedgraph.zst
///   coverage.pos.chr1.000123.pos.tsv.zst
///   coverage.part.chr1.000123.part.tsv.zst
pub fn parse_tile_index(file_name: &str) -> Option<u32> {
    for seg in file_name.rsplit('.') {
        if !seg.is_empty() && seg.chars().all(|c| c.is_ascii_digit()) {
            return seg.parse().ok();
        }
    }
    None
}

/// Merge compressed per-tile BedGraph chunks by *pure concatenation*
/// into a single `{final_name}` (also compressed).
///
/// Each tile file must be named `{per_tile_prefix}.{chr}.{index}.bedgraph.zst`.

/// Merge positional per-tile files (created in `temp_dir`) into one TSV.
/// Files must be named like: `{per_tile_prefix}.{chr}.{tile_index}.tsv`.
pub fn merge_positional_tiles(
    temp_dir: &std::path::Path,
    out_dir: &std::path::Path,
    chromosomes: &[String],
    per_tile_prefix: &str, // e.g. "coverage.pos" (whole-genome) or "coverage.pos.win" (windowed)
    final_name: &str,      // e.g. "coverage.per_position.tsv"
) -> Result<std::path::PathBuf> {
    let final_path = out_dir.join(final_name);
    let mut out = BufWriter::new(
        std::fs::File::create(&final_path)
            .with_context(|| format!("Creating merged output: {}", final_path.display()))?,
    );

    for chr in chromosomes {
        // Collect tile files for this chromosome from temp_dir
        let mut chr_files: Vec<(u32, std::path::PathBuf)> = Vec::new();
        for entry in std::fs::read_dir(temp_dir)
            .with_context(|| format!("Listing temp_dir: {}", temp_dir.display()))?
        {
            let path = entry?.path();
            if !path.is_file() {
                continue;
            }
            let fname = path.file_name().and_then(|s| s.to_str()).unwrap_or("");
            // Expect "{per_tile_prefix}.{chr}.{index}.tsv"
            if fname.starts_with(per_tile_prefix) && fname.contains(&format!(".{chr}.")) {
                if let Some(idx) = parse_tile_index(fname) {
                    chr_files.push((idx, path));
                }
            }
        }

        // Sort by tile index to preserve genomic order within chr
        chr_files.sort_by_key(|(i, _)| *i);

        // Stream copy each tile into the final file
        for (_idx, path) in chr_files {
            let mut f = std::fs::File::open(&path)
                .with_context(|| format!("Opening tile file: {}", path.display()))?;
            std::io::copy(&mut f, &mut out).with_context(|| {
                format!(
                    "Copying from {} into {}",
                    path.display(),
                    final_path.display()
                )
            })?;
        }
    }

    out.flush().context("Flushing merged output")?;
    Ok(final_path)
}

/// Concatenate compressed per-tile FINALS in chromosome/tile order.
/// A small compressed header frame is written first, then we stream-copy the tile frames.
/// This preserves full "concatenate zstd frames" behavior end-to-end.
pub fn concat_aligned_size_tile_finals(
    temp_dir: &std::path::Path,
    out_dir: &std::path::Path,
    chromosomes: &[String],
    per_tile_prefix: &str, // e.g., "<prefix>.fin"
    final_name: &str,      // e.g., "<prefix>.avg.tsv.zst"
    header_line: &str,     // single header line without trailing newline
) -> Result<std::path::PathBuf> {
    let final_path = out_dir.join(final_name);
    let mut out = BufWriter::new(
        std::fs::File::create(&final_path)
            .with_context(|| format!("Creating {}", final_path.display()))?,
    );

    // Write a compressed header frame first (so we never touch tile frames).
    let mut header_bytes = header_line.as_bytes().to_vec();
    header_bytes.push(b'\n');
    let header_frame =
        zstd::encode_all(&header_bytes[..], 3).context("Compressing header frame")?;
    out.write_all(&header_frame)?;

    // Then append each tile's compressed frame in genomic order
    for chr in chromosomes {
        // Collect tile files for this chromosome
        let mut chr_files: Vec<(u32, std::path::PathBuf)> = Vec::new();
        for entry in std::fs::read_dir(temp_dir)
            .with_context(|| format!("Listing {}", temp_dir.display()))?
        {
            let path = entry?.path();
            if !path.is_file() {
                continue;
            }
            let fname = path.file_name().and_then(|s| s.to_str()).unwrap_or("");
            if fname.starts_with(per_tile_prefix) && fname.contains(&format!(".{chr}.")) {
                if let Some(idx) = parse_tile_index(fname) {
                    chr_files.push((idx, path));
                }
            }
        }
        chr_files.sort_by_key(|(i, _)| *i);

        // Copy bytes verbatim (frame concatenation)
        for (_i, p) in chr_files {
            let mut f =
                std::fs::File::open(&p).with_context(|| format!("Opening {}", p.display()))?;
            std::io::copy(&mut f, &mut out).with_context(|| {
                format!("Copying {} into {}", p.display(), final_path.display())
            })?;
        }
    }

    out.flush()?;
    Ok(final_path)
}

pub fn adapt_fetch_to_extreme_windows(
    tile: &Tile,
    mode: &TileMode<'_>,
    chrom_len: u32,
) -> Option<(i64, i64)> {
    // Decide the fetch interval based on mode/windows.
    // For whole-genome positional: use the full tile fetch band.
    // For windowed runs: restrict to [min_overlapping_window, max_overlapping_window] ± halo,
    // intersected with the tile’s existing fetch band.
    let (fetch_from, fetch_to): (i64, i64) = match mode {
        // Whole positional coverage (no windows): keep the original tile fetch band
        TileMode::Positional { windows: None, .. } => {
            (tile.fetch_start as i64, tile.fetch_end as i64)
        }

        // Windowed positional coverage
        TileMode::Positional {
            windows: Some(wchr),
            ..
        } => {
            // Find the span of windows that overlap the tile core
            let mut found = false;
            let mut min_ws: u64 = u64::MAX;
            let mut max_we: u64 = 0;
            for &(ws, we, _) in windows_overlapping_core(wchr, tile.core_start, tile.core_end) {
                found = true;
                if ws < min_ws {
                    min_ws = ws;
                }
                if we > max_we {
                    max_we = we;
                }
            }
            // If nothing overlaps this core, skip this tile entirely
            if !found {
                return None;
            }

            // Use the tile's *actual* left/right halo (already edge-clamped)
            let left_halo = tile.core_start.saturating_sub(tile.fetch_start);
            let right_halo = tile.fetch_end.saturating_sub(tile.core_end);

            // Proposed narrower fetch band from window span ± halo
            let narrowed_start = (min_ws as u32).saturating_sub(left_halo);
            let narrowed_end = (max_we as u32).saturating_add(right_halo);

            // Intersect with the tile’s original fetch band, and clamp to chrom length
            let start_u32 = narrowed_start.max(tile.fetch_start);
            let end_u32 = narrowed_end.min(tile.fetch_end).min(chrom_len as u32);

            // It’s possible (though unlikely) numerical clamping collapses the band
            if start_u32 >= end_u32 {
                return None;
            }
            (start_u32 as i64, end_u32 as i64)
        }

        // Aggregates: same narrowing as windowed positional
        TileMode::AggregatesByBed { windows: wchr, .. } => {
            let mut found = false;
            let mut min_ws: u64 = u64::MAX;
            let mut max_we: u64 = 0;
            for &(ws, we, _) in windows_overlapping_core(wchr, tile.core_start, tile.core_end) {
                found = true;
                if ws < min_ws {
                    min_ws = ws;
                }
                if we > max_we {
                    max_we = we;
                }
            }
            if !found {
                return None;
            }

            let left_halo = tile.core_start.saturating_sub(tile.fetch_start);
            let right_halo = tile.fetch_end.saturating_sub(tile.core_end);

            let narrowed_start = (min_ws as u32).saturating_sub(left_halo);
            let narrowed_end = (max_we as u32).saturating_add(right_halo);

            let start_u32 = narrowed_start.max(tile.fetch_start);
            let end_u32 = narrowed_end.min(tile.fetch_end).min(chrom_len as u32);

            if start_u32 >= end_u32 {
                return None;
            }

            (start_u32 as i64, end_u32 as i64)
        }

        TileMode::AggregatesBySize { window_bp, .. } => {
            // Compute the extreme span of windows overlapping the tile core.
            let cs = tile.core_start as u64;
            let ce = tile.core_end as u64;
            let owned_window_bp = *window_bp;
            if cs >= chrom_len as u64 {
                return None;
            }

            // First window index that touches core_start, last that touches core_end-1.
            let k_lo = cs / owned_window_bp;
            let k_hi = (ce.saturating_sub(1)) / owned_window_bp;

            let min_ws = k_lo * owned_window_bp;
            let max_we = ((k_hi + 1) * owned_window_bp).min(chrom_len as u64);

            let left_halo = tile.core_start.saturating_sub(tile.fetch_start);
            let right_halo = tile.fetch_end.saturating_sub(tile.core_end);

            let narrowed_start = (min_ws as u32).saturating_sub(left_halo);
            let narrowed_end = (max_we as u32).saturating_add(right_halo);

            let start_u32 = narrowed_start.max(tile.fetch_start);
            let end_u32 = narrowed_end.min(tile.fetch_end).min(chrom_len as u32);

            if start_u32 >= end_u32 {
                return None;
            }

            (start_u32 as i64, end_u32 as i64)
        }
    };

    Some((fetch_from, fetch_to))
}

fn random_suffix(n: usize) -> String {
    rand::rng()
        .sample_iter(&Alphanumeric)
        .take(n)
        .map(char::from)
        .collect()
}

pub fn make_temp_dir(
    base_out: &std::path::Path,
    prefix: &str,
) -> anyhow::Result<std::path::PathBuf> {
    // Try a few times just in case
    for _ in 0..8 {
        let suffix = random_suffix(10);
        let p = base_out.join(format!("tmp.{prefix}.{suffix}"));
        if !p.exists() {
            std::fs::create_dir_all(&p)?;
            return Ok(p);
        }
    }
    // Fallback: timestamped
    let ts = chrono::Utc::now().timestamp_millis();
    let p = base_out.join(format!("tmp.{prefix}.{ts}"));
    std::fs::create_dir_all(&p)?;
    Ok(p)
}

// Round to number of decimals
pub fn round_to(x: f64, decimals: i32) -> f64 {
    if decimals <= 0 {
        return x.round();
    }
    let f = 10f64.powi(decimals);
    (x * f).round() / f
}

// Format as compactly as possible
pub fn format_number_simplify(v: f64, decimals: i32) -> String {
    // For integer formatting, don't trim trailing zeros from values like "30"
    if decimals <= 0 {
        let s = format!("{:.0}", v);
        return if s == "-0" { "0".to_string() } else { s };
    }
    // For decimal formatting, trim trailing zeros and a trailing dot
    let s = format!("{:.*}", decimals as usize, v);
    let s = s.trim_end_matches('0').trim_end_matches('.').to_string();
    if s == "-0" { "0".to_string() } else { s }
}

/// Emit BedGraph runs for cov[a..b), skipping masked bases.
/// - Standard 4-column BedGraph: chrom, start, end, value
/// - Masked bases (mask==1) are **not written** → gaps in the track.
/// - Values are rounded to `decimals` before comparing/printing to avoid run explosion.
pub fn emit_bedgraph_runs<W: Write>(
    chr: &str,
    cov: &[f32],
    mask: Option<&[u8]>,  // 1 = blacklisted(masked), 0 = allowed
    a: usize,             // Local start (inclusive)
    b: usize,             // Local end (exclusive)
    start_abs: u64,       // Absolute position of index 0 in `cov` (tile.core_start)
    decimals: i32,        // Decimals to round coverage
    keep_zero_runs: bool, // Whether to write zero-runs
    out: &mut W,
) -> Result<()> {
    if a >= b {
        return Ok(());
    }

    let m = mask.unwrap_or(&[]);
    let mut i = a;

    while i < b {
        // Skip masked stretch
        if !m.is_empty() && m[i] == 1 {
            i += 1;
            continue;
        }

        // Start unmasked run
        let run_start = i;
        let v0 = round_to(cov[i] as f64, decimals);

        let mut j = i + 1;
        while j < b {
            if !m.is_empty() && m[j] == 1 {
                break;
            }
            let vj = round_to(cov[j] as f64, decimals);
            if vj != v0 {
                break;
            }
            j += 1;
        }

        // Skip zero-runs unless specified otherwise
        if v0 == 0.0 && !keep_zero_runs {
            i = j;
            continue;
        }

        // Emit [run_start, j)
        let s_abs = start_abs + run_start as u64;
        let e_abs = start_abs + j as u64;
        writeln!(
            out,
            "{}\t{}\t{}\t{}",
            chr,
            s_abs,
            e_abs,
            format_number_simplify(v0, decimals)
        )?;

        i = j;
    }

    Ok(())
}

/// Emit either bedgraph or TSV runs for a single window cov[a..b):
/// chrom  start  end  value  (optional: orig_idx)
///
/// - Skips masked bases (creates gaps).
/// - Run-length encodes equal values inside the window to reduce size.
/// - Keeps `orig_idx` for downstream grouping.
/// - Use this for **windowed positional** outputs.
pub fn emit_windowed_runs<W: Write>(
    chr: &str,
    cov: &[f32],
    mask: Option<&[u8]>,   // 1 = blacklisted(masked), 0 = allowed
    a: usize,              // Local start (inclusive)
    b: usize,              // Local end (exclusive)
    start_abs: u64,        // Absolute position of index 0 in `cov` (tile.core_start)
    orig_idx: Option<u64>, // Window’s original index
    decimals: i32,         // Decimals to round coverage
    keep_zero_runs: bool,  // Whether to write zero-runs
    out: &mut W,
) -> Result<()> {
    if a >= b {
        return Ok(());
    }

    let m = mask.unwrap_or(&[]);
    let mut i = a;

    while i < b {
        // skip masked
        if !m.is_empty() && m[i] == 1 {
            i += 1;
            continue;
        }

        // start run
        let run_start = i;
        let v0 = round_to(cov[i] as f64, decimals);

        let mut j = i + 1;
        while j < b {
            if !m.is_empty() && m[j] == 1 {
                break;
            }
            let vj = round_to(cov[j] as f64, decimals);
            if vj != v0 {
                break;
            }
            j += 1;
        }

        // Skip zero-runs unless specified otherwise
        if v0 == 0.0 && !keep_zero_runs {
            i = j;
            continue;
        }

        // emit: chr  start  end  value  orig_idx
        let s_abs = start_abs + run_start as u64;
        let e_abs = start_abs + j as u64;

        if let Some(idx) = orig_idx {
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}",
                chr,
                s_abs,
                e_abs,
                format_number_simplify(v0, decimals),
                idx
            )?;
        } else {
            writeln!(
                out,
                "{}\t{}\t{}\t{}",
                chr,
                s_abs,
                e_abs,
                format_number_simplify(v0, decimals)
            )?;
        }

        i = j;
    }

    Ok(())
}

/// Sum coverage, respecting masked mode
#[inline]
pub fn tile_sum_and_counts(
    a: usize,
    b: usize,
    masked: bool,
    ps_all: &[f64],
    ps_allow: Option<&[f64]>,
    cnt_allow: Option<&[u32]>,
    mask: Option<&[u8]>,
) -> (f64, u64, u64) {
    let sum = if masked {
        if let Some(pa) = ps_allow {
            pa[b] - pa[a]
        } else {
            ps_all[b] - ps_all[a]
        }
    } else {
        ps_all[b] - ps_all[a]
    };

    let span = (b - a) as u64;

    let allowed = if masked {
        if let Some(cnt) = cnt_allow {
            (cnt[b] - cnt[a]) as u64
        } else if let Some(m) = mask {
            // Fall back to a tiny scan if count psum is absent
            let mut ok = 0u64;
            for i in a..b {
                if m[i] == 0 {
                    ok += 1;
                }
            }
            ok
        } else {
            span
        }
    } else {
        span
    };

    let blacklisted = span - allowed;
    (sum, allowed, blacklisted)
}

#[inline]
pub fn finalize_value(
    sum: f64,
    allowed_positions: u64,
    unmasked_span_bp: u64, // end-start when unmasked mode; ignored when masked mode
    masked: bool,
    mode: &CoverageWindowAction,
) -> f64 {
    match mode {
        CoverageWindowAction::Average => {
            if masked {
                if allowed_positions == 0 {
                    0.0
                } else {
                    sum / allowed_positions as f64
                }
            } else {
                if unmasked_span_bp == 0 {
                    0.0
                } else {
                    sum / unmasked_span_bp as f64
                }
            }
        }
        CoverageWindowAction::Total => sum,
        _ => unreachable!(),
    }
}
