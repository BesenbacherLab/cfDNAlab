#[cfg(feature = "cmd_fcoverage")]
use crate::shared::bed::GroupedWindows;
use crate::shared::bed::Windows;
use anyhow::{Result, ensure};
use fxhash::FxHashMap;

#[cfg(any(feature = "cmd_ends", feature = "cmd_fragment_kmers"))]
pub(crate) use window_context::WindowContext;

#[cfg(any(
    feature = "cmd_ends",
    feature = "cmd_fragment_kmers",
    feature = "cmd_gc_bias",
    feature = "cmd_lengths",
    feature = "cmd_wps_peaks"
))]
pub(crate) use window_offsets::compute_window_offsets;

#[cfg(any(
    feature = "cmd_ends",
    feature = "cmd_fragment_kmers",
    feature = "cmd_lengths"
))]
pub(crate) use window_bin_info::{WindowBinInfo, build_bin_info};

#[cfg(any(feature = "cmd_ends", feature = "cmd_fragment_kmers"))]
pub(crate) use bin_info_tsv_writer::write_bin_info_tsv;

#[cfg(any(feature = "cmd_ends"))]
pub(crate) use group_index_writer::write_group_index_with_blacklist_tsv;

/// Validate that an ordinary BED map contains at least one selected window.
///
/// BED loaders keep empty chromosome entries when a chromosome whitelist is supplied. This helper
/// checks the actual surviving windows so commands can fail early when a BED file is empty after
/// chromosome filtering.
pub(crate) fn ensure_plain_bed_windows_not_empty(
    windows_map: &FxHashMap<String, Windows>,
) -> Result<()> {
    ensure!(
        windows_map.values().any(|windows| !windows.is_empty()),
        "BED file did not contain any valid windows on the selected chromosomes"
    );
    Ok(())
}

/// Validate that a grouped BED map contains at least one selected window.
///
/// Grouped BED loaders keep empty chromosome entries when a chromosome whitelist is supplied. This
/// checks the actual surviving grouped windows so commands fail early when every grouped BED row is
/// removed by chromosome filtering.
#[cfg(feature = "cmd_fcoverage")]
pub(crate) fn ensure_grouped_bed_windows_not_empty(
    windows_map: &FxHashMap<String, GroupedWindows>,
) -> Result<()> {
    ensure!(
        windows_map.values().any(|windows| !windows.is_empty()),
        "grouped BED file did not contain any valid windows on the selected chromosomes"
    );
    Ok(())
}

#[cfg(any(feature = "cmd_ends", feature = "cmd_fragment_kmers"))]
mod window_context {
    use crate::{commands::cli_common::WindowSpec, shared::interval::IndexedInterval};

    /// Lightweight view into the window configuration for a given tile.
    ///
    /// Stores the context needed to convert chromosome-local indices into global window ids for a tile.
    pub(crate) struct WindowContext<'a> {
        pub(crate) spec: &'a WindowSpec,
        pub(crate) windows: Option<&'a [IndexedInterval<u64>]>,
        pub(crate) chr_idx_offset: u64,
    }

    impl<'a> WindowContext<'a> {
        #[inline]
        /// Return the per-chromosome windows slice when operating in BED mode.
        #[cfg(feature = "cmd_fragment_kmers")]
        pub(crate) fn windows_slice(&self) -> Option<&'a [IndexedInterval<u64>]> {
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
        pub(crate) fn original_idx(&self, chrom_window_idx: usize) -> u64 {
            match self.spec {
                WindowSpec::Global => 0,
                WindowSpec::Size(window_bp) => {
                    debug_assert_ne!(*window_bp, 0);
                    self.chr_idx_offset
                        .checked_add(chrom_window_idx as u64)
                        .expect("window index overflow for size-based windows")
                }
                WindowSpec::Bed(_) => {
                    // In BED mode the stored idx is the original window index that
                    // downstream arrays and output rows are keyed by.
                    self.windows.expect("windows slice required for BED mode")[chrom_window_idx]
                        .idx()
                }
            }
        }
    }
}

#[cfg(any(
    feature = "cmd_ends",
    feature = "cmd_fragment_kmers",
    feature = "cmd_gc_bias",
    feature = "cmd_lengths",
    feature = "cmd_wps_peaks"
))]
mod window_offsets {
    use crate::{
        commands::cli_common::WindowSpec,
        shared::{bam::Contigs, bed::Windows},
    };
    use anyhow::{Context, Result};
    use fxhash::FxHashMap;

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
    pub(crate) fn compute_window_offsets(
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
                let mut running_window_idx_offset = 0u64;
                for chr in chromosomes {
                    offsets.insert(chr.clone(), running_window_idx_offset);
                    let &(_, len_u32) = contigs
                        .contigs
                        .get(chr)
                        .with_context(|| format!("missing contig length for '{}'", chr))?;
                    let len = len_u32 as u64;
                    let bins = if len == 0 { 0 } else { len.div_ceil(*size) };
                    running_window_idx_offset = running_window_idx_offset.saturating_add(bins);
                }
                Ok((running_window_idx_offset, offsets))
            }
            WindowSpec::Bed(_) => {
                let win_map = windows_map
                    .with_context(|| "window map required for --by-bed mode".to_string())?;
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
}

#[cfg(any(
    feature = "cmd_ends",
    feature = "cmd_fragment_kmers",
    feature = "cmd_lengths"
))]
mod window_bin_info {
    use crate::{
        commands::cli_common::WindowSpec,
        shared::{
            bam::Contigs, bed::Windows, blacklist::compute_blacklist_overlap, interval::Interval,
        },
    };
    use anyhow::{Context, Result};
    use fxhash::FxHashMap;

    /// Metadata for one ordinary output window.
    ///
    /// `output_index` is the global output row index for fixed-size windows. In BED mode, it is the
    /// original BED window index preserved by the loader.
    #[derive(Clone, Debug, PartialEq)]
    pub(crate) struct WindowBinInfo {
        pub(crate) chromosome: String,
        pub(crate) start: u64,
        pub(crate) end: u64,
        pub(crate) output_index: u64,
        pub(crate) blacklisted_fraction: f64,
    }

    /// Build per-window metadata (coordinates, blacklist overlap, etc.) for downstream consumers.
    ///
    /// When running in BED mode the `original_idx` embedded in the loaded windows is preserved so the
    /// caller must continue using that identifier when addressing global vectors.
    pub(crate) fn build_bin_info(
        window_opt: &WindowSpec,
        chromosomes: &[String],
        contigs: &Contigs,
        windows_map: Option<&FxHashMap<String, Windows>>,
        blacklist_map: &FxHashMap<String, Vec<Interval<u64>>>,
        chr_offsets: &FxHashMap<String, u64>,
    ) -> Result<Vec<WindowBinInfo>> {
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
                    let mut blacklist_ptr = 0usize;
                    let blacklist_intervals =
                        blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]);
                    let chr_window_idx_offset = *chr_offsets.get(chr).unwrap_or(&0);

                    while start < len {
                        let end = (start + *size).min(len);
                        let overlap = compute_blacklist_overlap(
                            blacklist_intervals,
                            Interval::new(start, end)?,
                            0,
                            &mut blacklist_ptr,
                        );
                        out.push(WindowBinInfo {
                            chromosome: chr.clone(),
                            start,
                            end,
                            output_index: chr_window_idx_offset + local_idx,
                            blacklisted_fraction: overlap,
                        });
                        start += *size;
                        local_idx += 1;
                    }
                }
                Ok(out)
            }
            WindowSpec::Bed(_) => {
                let win_map = windows_map
                    .with_context(|| "window map required for --by-bed mode".to_string())?;
                for chr in chromosomes {
                    let windows = win_map.get(chr).map(|w| w.as_slice()).unwrap_or(&[]);
                    let mut blacklist_ptr = 0usize;
                    let blacklist_intervals =
                        blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]);
                    for window in windows {
                        // Preserve the original BED window index rather than the chromosome-local loop
                        // position.
                        let (start, end, original_idx) = window.as_tuple();
                        let overlap = compute_blacklist_overlap(
                            blacklist_intervals,
                            Interval::new(start, end)?,
                            0,
                            &mut blacklist_ptr,
                        );
                        out.push(WindowBinInfo {
                            chromosome: chr.clone(),
                            start,
                            end,
                            output_index: original_idx,
                            blacklisted_fraction: overlap,
                        });
                    }
                }
                out.sort_unstable_by_key(|entry| entry.output_index);
                Ok(out)
            }
        }
    }
}

#[cfg(any(feature = "cmd_ends", feature = "cmd_fragment_kmers"))]
mod bin_info_tsv_writer {
    use super::WindowBinInfo;
    use crate::shared::io::create_text_writer;
    use anyhow::{Context, Result};
    use std::{io::Write, path::Path};

    /// Write ordinary window metadata to a TSV next to matrix outputs.
    ///
    /// Output has header `chrom\tstart\tend\tblacklisted_fraction`.
    pub(crate) fn write_bin_info_tsv(
        output_path: impl AsRef<Path>,
        bin_info: &[WindowBinInfo],
    ) -> Result<()> {
        let mut writer = create_text_writer(output_path.as_ref()).context("creating bins TSV")?;
        writeln!(writer, "chrom\tstart\tend\tblacklisted_fraction")
            .context("writing bins TSV header")?;
        for entry in bin_info {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}",
                entry.chromosome, entry.start, entry.end, entry.blacklisted_fraction
            )
            .context("writing bins TSV row")?;
        }
        writer.finish().context("finalizing bins.tsv writer")?;
        Ok(())
    }
}

#[cfg(any(feature = "cmd_ends"))]
mod group_index_writer {
    use crate::shared::{
        bed::GroupedWindows, blacklist::compute_blacklist_overlap, interval::Interval,
        io::create_text_writer,
    };
    use anyhow::{Context, Result};
    use fxhash::FxHashMap;
    use std::{io::Write, path::Path};

    /// Write grouped distribution row metadata to a TSV next to grouped outputs.
    ///
    /// Output always has `group_idx` and `group_name`.
    ///
    /// When `include_blacklisted_fraction` is true, the output also includes
    /// `blacklisted_fraction`, aggregated across all intervals assigned to the group and weighted by
    /// interval width:
    ///
    /// `sum(interval_blacklisted_bp) / sum(interval_bp)`
    ///
    /// Intervals are counted exactly as loaded, so overlapping intervals in the same group contribute
    /// separately to both the numerator and denominator.
    pub(crate) fn write_group_index_with_blacklist_tsv(
        output_path: impl AsRef<Path>,
        group_idx_to_name: &FxHashMap<u64, String>,
        chromosomes: &[String],
        grouped_windows_map: &FxHashMap<String, GroupedWindows>,
        blacklist_map: &FxHashMap<String, Vec<Interval<u64>>>,
        include_blacklisted_fraction: bool,
    ) -> Result<()> {
        let mut writer =
            create_text_writer(output_path.as_ref()).context("creating grouped group-index TSV")?;

        let mut entries: Vec<(u64, &str)> = group_idx_to_name
            .iter()
            .map(|(idx, name)| (*idx, name.as_str()))
            .collect();
        entries.sort_unstable_by_key(|(idx, _)| *idx);

        if !include_blacklisted_fraction {
            writeln!(writer, "group_idx\tgroup_name")
                .context("writing grouped group-index TSV header")?;
            for (group_idx, group_name) in entries {
                let group_name = group_name.replace('\t', "    ").replace('\n', " ");
                writeln!(writer, "{group_idx}\t{group_name}")
                    .context("writing grouped group-index TSV row")?;
            }
        } else {
            let mut total_group_bp: FxHashMap<u64, u64> = FxHashMap::default();
            let mut blacklisted_group_bp: FxHashMap<u64, f64> = FxHashMap::default();

            for chr in chromosomes {
                let windows = grouped_windows_map
                    .get(chr)
                    .map(|windows| windows.windows_as_slice())
                    .unwrap_or(&[]);
                let blacklist_intervals =
                    blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]);
                let mut blacklist_ptr = 0usize;
                for window in windows {
                    let (start, end, group_idx) = window.as_tuple();
                    let window_bp = end
                        .checked_sub(start)
                        .context("grouped window end must be >= start")?;
                    let blacklist_overlap_fraction = compute_blacklist_overlap(
                        blacklist_intervals,
                        Interval::new(start, end)?,
                        0,
                        &mut blacklist_ptr,
                    );
                    *total_group_bp.entry(group_idx).or_insert(0) += window_bp;
                    *blacklisted_group_bp.entry(group_idx).or_insert(0.0) +=
                        blacklist_overlap_fraction * window_bp as f64;
                }
            }

            writeln!(writer, "group_idx\tgroup_name\tblacklisted_fraction")
                .context("writing grouped group-index TSV header")?;
            for (group_idx, group_name) in entries {
                let total_bp = *total_group_bp.get(&group_idx).unwrap_or(&0);
                let blacklisted_fraction = if total_bp == 0 {
                    0.0
                } else {
                    blacklisted_group_bp.get(&group_idx).copied().unwrap_or(0.0) / total_bp as f64
                };
                let group_name = group_name.replace('\t', "    ").replace('\n', " ");
                writeln!(writer, "{group_idx}\t{group_name}\t{blacklisted_fraction}")
                    .context("writing grouped group-index TSV row")?;
            }
        }

        writer
            .finish()
            .context("finalizing grouped group-index TSV writer")?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    include!("windowing_tests.rs");
}
