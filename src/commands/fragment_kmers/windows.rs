use crate::{
    commands::cli_common::WindowSpec,
    shared::{bam::Contigs, bed::Windows, blacklist::compute_blacklist_overlap},
};
use anyhow::{Context, Result};
use fxhash::FxHashMap;

/// Lightweight view into the window configuration for a given tile.
///
/// Stores the context needed to convert chromosome-local indices into global window ids for a tile.
pub struct WindowContext<'a> {
    pub spec: &'a WindowSpec,
    pub windows: Option<&'a [(u64, u64, u64)]>,
    pub chr_idx_offset: u64,
}

impl<'a> WindowContext<'a> {
    #[inline]
    /// Return the per-chromosome windows slice when operating in BED mode.
    pub fn windows_slice(&self) -> Option<&'a [(u64, u64, u64)]> {
        self.windows
    }

    #[inline]
    /// Map the provided chromosome-local window index to the global window identifier expected
    /// downstream.
    ///
    /// Parameters
    /// -----------
    /// `chrom_window_idx`: Index of the window relative to the *chromosome*, as supplied by
    /// [`find_overlapping_windows`] (it counts from the start of the chromosome, not the tile).
    ///
    pub fn original_idx(&self, chrom_window_idx: usize) -> u64 {
        match self.spec {
            WindowSpec::Global => 0,
            WindowSpec::Size(window_bp) => {
                debug_assert_ne!(*window_bp, 0);
                self.chr_idx_offset
                    .checked_add(chrom_window_idx as u64)
                    .expect("window index overflow for size-based windows")
            }
            WindowSpec::Bed(_) => {
                self.windows.expect("windows slice required for BED mode")[chrom_window_idx].2
            }
        }
    }
}

/// Compute the global window count together with per-chromosome index offsets.
///
/// Returns `(total_windows, chr_idx_offsets)` where `total_windows` is the expected length of the
/// window output vector and `chr_idx_offsets` maps each chromosome to the index of its first window
/// in that global vector.
///
/// Mode specifics:
/// - `Global`: there is exactly one window covering the whole genome; all chromosomes map to offset
///   `0`.
/// - `Size(window_bp)`: windows are created from the chromosome start in contiguous bins, so the
///   offset for chromosome *k* is the cumulative number of bins emitted for chromosomes `< k`.
/// - `Bed`: offsets remain `0` because each BED entry already carries its own globally unique
///   `original_idx`; consumers must use that `original_idx` when addressing global arrays.
pub fn compute_window_offsets(
    window_opt: &WindowSpec,
    chromosomes: &[String],
    contigs: &Contigs,
    windows_map: Option<&FxHashMap<String, Windows>>,
) -> Result<(u64, FxHashMap<String, u64>)> {
    let mut offsets: FxHashMap<String, u64> = FxHashMap::default();

    match window_opt {
        WindowSpec::Global => {
            for chr in chromosomes {
                offsets.insert(chr.clone(), 0);
            }
            Ok((1, offsets))
        }
        WindowSpec::Size(size) => {
            let mut running_window_idx_offset = 0u64; // number of windows seen so far
            for chr in chromosomes {
                offsets.insert(chr.clone(), running_window_idx_offset);
                let &(_, len_u32) = contigs
                    .contigs
                    .get(chr)
                    .with_context(|| format!("missing contig length for '{}'", chr))?;
                let len = len_u32 as u64;
                let bins = if len == 0 {
                    0
                } else {
                    (len + *size - 1) / *size
                };
                running_window_idx_offset = running_window_idx_offset.saturating_add(bins);
            }
            Ok((running_window_idx_offset, offsets))
        }
        WindowSpec::Bed(_) => {
            let win_map =
                windows_map.with_context(|| "window map required for --by-bed mode".to_string())?;
            let mut total = 0u64;
            for chr in chromosomes {
                // BED entries already encode their global `original_idx`, so the reducer should use
                // that value instead of a chromosome offset.
                offsets.insert(chr.clone(), 0);
                if let Some(windows) = win_map.get(chr) {
                    total = total.saturating_add(windows.len() as u64);
                }
            }
            Ok((total, offsets))
        }
    }
}

/// Build per-window metadata (coordinates, blacklist overlap, etc.) for downstream consumers.
///
/// When running in BED mode the `original_idx` embedded in the loaded windows is preserved so the
/// caller must continue using that identifier when addressing global vectors.
pub fn build_bin_info(
    window_opt: &WindowSpec,
    chromosomes: &[String],
    contigs: &Contigs,
    windows_map: Option<&FxHashMap<String, Windows>>,
    blacklist_map: &FxHashMap<String, Vec<(u64, u64)>>,
    chr_offsets: &FxHashMap<String, u64>,
) -> Result<Vec<(String, u64, u64, u64, f64)>> {
    let mut out = Vec::new();

    match window_opt {
        WindowSpec::Global => Ok(out),
        WindowSpec::Size(size) => {
            for chr in chromosomes {
                let &(_, len_u32) = contigs
                    .contigs
                    .get(chr)
                    .with_context(|| format!("missing contig length for '{}'", chr))?;
                let len = len_u32 as u64;
                let mut start = 0u64;
                let mut local_idx = 0u64;
                let mut bl_ptr = 0usize;
                let blacklist_intervals =
                    blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]);
                let chr_window_idx_offset = *chr_offsets.get(chr).unwrap_or(&0);

                while start < len {
                    let end = (start + *size).min(len);
                    let overlap =
                        compute_blacklist_overlap(blacklist_intervals, start, end, 0, &mut bl_ptr);
                    out.push((
                        chr.clone(),
                        start,
                        end,
                        chr_window_idx_offset + local_idx,
                        overlap,
                    ));
                    start += *size;
                    local_idx += 1;
                }
            }
            Ok(out)
        }
        WindowSpec::Bed(_) => {
            let win_map =
                windows_map.with_context(|| "window map required for --by-bed mode".to_string())?;
            for chr in chromosomes {
                let windows = win_map.get(chr).map(|w| w.as_slice()).unwrap_or(&[]);
                let mut bl_ptr = 0usize;
                let blacklist_intervals =
                    blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]);
                for &(start, end, original_idx) in windows {
                    let overlap =
                        compute_blacklist_overlap(blacklist_intervals, start, end, 0, &mut bl_ptr);
                    out.push((chr.clone(), start, end, original_idx, overlap));
                }
            }
            out.sort_unstable_by_key(|entry| entry.3);
            Ok(out)
        }
    }
}
