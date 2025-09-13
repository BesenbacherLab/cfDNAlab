use crate::utils::bam::create_chromosome_reader;
use crate::utils::coverage::coverage_prefix::CoveragePrefix;
use crate::utils::fragment::minimal_fragment::Fragment;
use crate::utils::fragment::segment_fragment::FragmentWithSegments;
use anyhow::Result;
use std::io::Write;

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

/// Build non-overlapping core tiles with a fetch halo on each side
///
/// - tile_bp: target core size in bases (e.g. 20_000_000)
/// - halo_bp: fetch extension on both sides (e.g. max_fragment_length)
pub fn build_tiles(
    bam_path: &std::path::Path,
    chromosomes: &[String],
    tile_bp: u32,
    halo_bp: u32,
) -> Result<Vec<Tile>> {
    let mut tiles = Vec::new();

    for chr in chromosomes {
        // Reuse your helper to get tid and length once per chromosome
        let (_rdr, tid, chrom_len_u64) = create_chromosome_reader(bam_path, chr)?;
        let chrom_len = chrom_len_u64 as u32;

        let mut start = 0u32;
        let mut idx = 0u32;
        while start < chrom_len {
            let core_end = (start.saturating_add(tile_bp)).min(chrom_len);

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

    Ok(tiles)
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
    /// Whole positional coverage for the core, or windowed positional coverage
    Positional {
        windows: Option<&'w [(u64, u64, u64)]>, // per-chr windows if provided
        out_path: std::path::PathBuf,           // per-tile file path
    },
    /// Aggregate windows: write per-tile partials for later reduce
    Aggregates {
        windows: &'w [(u64, u64, u64)], // per-chr windows
        masked: bool,                   // use masked counts/sums
        out_path: std::path::PathBuf,   // per-tile partials
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

/// Merge positional per-tile files in (chr, index) order into one TSV
pub fn merge_positional_tiles(
    out_dir: &std::path::Path,
    chromosomes: &[String],
    per_tile_prefix: &str, // e.g. "coverage.pos"
    final_name: &str,      // e.g. "coverage.per_positions.tsv"
) -> anyhow::Result<std::path::PathBuf> {
    use std::io::{BufRead, BufReader, Write};

    let final_path = out_dir.join(final_name);
    let mut out = std::io::BufWriter::new(std::fs::File::create(&final_path)?);

    for chr in chromosomes {
        // List files for this chr and sort by index suffix
        let mut chr_files: Vec<(u32, std::path::PathBuf)> = Vec::new();
        for entry in std::fs::read_dir(out_dir)? {
            let p = entry?.path();
            if !p.is_file() {
                continue;
            }
            let fname = p.file_name().and_then(|s| s.to_str()).unwrap_or("");
            // Expect "{per_tile_prefix}.{chr}.{index}.tsv"
            if fname.starts_with(per_tile_prefix) && fname.contains(&format!(".{chr}.")) {
                if let Some(idx_str) = fname.rsplit('.').nth(2) {
                    if let Ok(idx) = idx_str.parse::<u32>() {
                        chr_files.push((idx, p));
                    }
                }
            }
        }
        chr_files.sort_by_key(|(i, _)| *i);

        // Stream copy in order
        for (_, path) in chr_files {
            let f = std::fs::File::open(&path)?;
            let mut r = BufReader::new(f);
            let mut line = String::new();
            while r.read_line(&mut line)? != 0 {
                out.write_all(line.as_bytes())?;
                line.clear();
            }
        }
    }

    out.flush()?;
    Ok(final_path)
}

/// Reduce aggregate partials for a chromosome to final windows in order
///
/// - windows_chr: same ordering used to assign original_idx
/// - masked: if true, output averages will use allowed_count; otherwise span length
pub fn reduce_aggregates_for_chr(
    chr: &str,
    out_dir: &std::path::Path,
    partial_prefix: &str, // e.g. "coverage.part"
    windows_chr: &[(u64, u64, u64)],
    masked: bool,
    final_writer: &mut std::io::BufWriter<std::fs::File>,
) -> anyhow::Result<()> {
    use std::io::{BufRead, BufReader};

    // Prepare accumulators sized by number of windows in this chr
    let n = windows_chr.len();
    let mut sum = vec![0.0_f64; n];
    let mut allowed = vec![0u64; n];
    let mut blacklisted = vec![0u64; n];

    // Collect & sort all partial files for this chr
    let mut chr_files: Vec<(u32, std::path::PathBuf)> = Vec::new();
    for entry in std::fs::read_dir(out_dir)? {
        let p = entry?.path();
        if !p.is_file() {
            continue;
        }
        let fname = p.file_name().and_then(|s| s.to_str()).unwrap_or("");
        // "{partial_prefix}.{chr}.{index}.tsv"
        if fname.starts_with(partial_prefix) && fname.contains(&format!(".{chr}.")) {
            if let Some(idx_str) = fname.rsplit('.').nth(2) {
                if let Ok(idx) = idx_str.parse::<u32>() {
                    chr_files.push((idx, p));
                }
            }
        }
    }
    chr_files.sort_by_key(|(i, _)| *i);

    // Accumulate
    for (_, path) in chr_files {
        let f = std::fs::File::open(&path)?;
        let r = BufReader::new(f);
        for line in r.lines() {
            let line = line?;
            // idx  sum  allowed  blacklisted
            let mut it = line.split('\t');
            let idx: usize = it.next().unwrap().parse()?;
            let s: f64 = it.next().unwrap().parse()?;
            let a: u64 = it.next().unwrap().parse()?;
            let b: u64 = it.next().unwrap().parse()?;
            sum[idx] += s;
            allowed[idx] += a;
            blacklisted[idx] += b;
        }
    }

    // Emit final rows in window order for this chromosome
    for (i, &(window_start, window_end, _orig_idx)) in windows_chr.iter().enumerate() {
        let span = (window_end - window_start) as f64;
        let a = allowed[i] as f64;
        let avg = if masked {
            if a == 0.0 { 0.0 } else { sum[i] / a }
        } else {
            if span == 0.0 { 0.0 } else { sum[i] / span }
        };
        // chr  start  end  avg  bl_pos
        writeln!(
            final_writer,
            "{}\t{}\t{}\t{}\t{}",
            chr, window_start, window_end, avg, blacklisted[i],
        )?;
    }

    Ok(())
}
