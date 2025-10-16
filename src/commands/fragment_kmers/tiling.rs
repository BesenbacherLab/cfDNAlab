use crate::{
    commands::{
        cli_common::WindowSpec, counters::FragmentKmersCounters,
        fragment_kmers::windows::WindowContext,
    },
    shared::{
        kmers::{
            kmer_codec::{Kmer, KmerOrientation, KmerSpec},
            process_counts::{DecodedCounts, split_and_decode_counts},
        },
        tiled_run::{Tile, TileWindowSpan, clamp_fetch_to_window_span, tile_window_min_max},
    },
};
use anyhow::{Context, Result, bail};
use bincode::{
    config::standard,
    serde::{decode_from_std_read, encode_into_std_write},
};
use fxhash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::{
    fs::{self, File},
    io::{BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

#[cfg_attr(not(test), doc(hidden))]
#[derive(Debug, Serialize, Deserialize)]
pub struct TileKmerCountEntry {
    pub k: u8,
    pub code: u64,
    pub orientation: KmerOrientation,
    pub value: f64,
}

#[cfg_attr(not(test), doc(hidden))]
#[derive(Debug, Serialize, Deserialize)]
pub struct TileWindowCounts {
    pub original_idx: u64,
    pub entries: Vec<TileKmerCountEntry>,
}

/// Persist per-tile k-mer counts so they can be merged after parallel tile processing.
pub fn serialize_tile_counts(path: &Path, payload: &[TileWindowCounts]) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("creating tile counts file: {}", path.display()))?;
    let mut writer = BufWriter::with_capacity(512 * 1024, file);
    encode_into_std_write(payload, &mut writer, standard())
        .with_context(|| format!("serialising tile counts to {}", path.display()))?;
    writer.flush().with_context(|| {
        format!(
            "flushing tile counts file after serialisation: {}",
            path.display()
        )
    })
}

/// Load counts created by [`serialize_tile_counts`] during the reduction phase.
pub fn deserialize_tile_counts(path: &Path) -> Result<Vec<TileWindowCounts>> {
    let file = File::open(path)
        .with_context(|| format!("opening tile counts file: {}", path.display()))?;
    let mut reader = BufReader::with_capacity(512 * 1024, file);
    decode_from_std_read(&mut reader, standard())
        .with_context(|| format!("deserialising tile counts from {}", path.display()))
}

/// Per-tile bookkeeping for intermediate count files and fragment counters.
pub struct TileResult {
    pub chr: String,
    pub counts_path: Option<PathBuf>,
    pub counter: FragmentKmersCounters,
}

/// Determine the genomic span to request from the BAM reader for a tile.
pub fn determine_fetch_span(
    tile: &Tile,
    window_ctx: &WindowContext,
    tile_window_span: Option<&TileWindowSpan>,
    chrom_len: u64,
) -> Option<(i64, i64)> {
    let chrom_len_u32 = chrom_len.min(u32::MAX as u64) as u32;
    match window_ctx.spec {
        WindowSpec::Global | WindowSpec::Size(_) => {
            Some((tile.fetch_start as i64, tile.fetch_end as i64))
        }
        WindowSpec::Bed(_) => {
            let windows = window_ctx.windows_slice()?;
            let (min_ws, max_we) = tile_window_min_max(windows, tile, tile_window_span)?;
            clamp_fetch_to_window_span(tile, chrom_len.min(chrom_len_u32 as u64), min_ws, max_we)
        }
    }
}

/// Reduce per-tile count payloads into a dense vector aligned with the global window order.
#[cfg_attr(not(test), doc(hidden))]
pub fn merge_tile_counts<I>(
    payloads: I,
    total_windows: usize,
    kmer_specs: &FxHashMap<u8, KmerSpec>,
) -> Result<Vec<DecodedCounts>>
where
    I: IntoIterator<Item = Vec<TileWindowCounts>>,
{
    let mut aggregated_counts: FxHashMap<u64, FxHashMap<Kmer, f64>> = FxHashMap::default();

    for payload in payloads {
        for window_counts in payload {
            let entry = aggregated_counts
                .entry(window_counts.original_idx)
                .or_insert_with(FxHashMap::default);
            for count in window_counts.entries {
                let kmer = Kmer {
                    k: count.k,
                    code: count.code,
                    orientation: count.orientation,
                };
                *entry.entry(kmer).or_insert(0.0) += count.value;
            }
        }
    }

    let empty_counts: FxHashMap<Kmer, f64> = FxHashMap::default();
    let mut all_bins: Vec<DecodedCounts> = Vec::with_capacity(total_windows);
    for idx in 0..total_windows {
        if let Some(counts) = aggregated_counts.remove(&(idx as u64)) {
            all_bins.push(split_and_decode_counts(&counts, kmer_specs));
        } else {
            all_bins.push(split_and_decode_counts(&empty_counts, kmer_specs));
        }
    }

    if !aggregated_counts.is_empty() {
        bail!(
            "Received counts for unexpected window indices: {:?}",
            aggregated_counts.keys().collect::<Vec<&u64>>()
        );
    }

    Ok(all_bins)
}

pub fn reduce_chromosome_tile_results(
    tile_results: Vec<TileResult>,
) -> Result<Vec<TileWindowCounts>> {
    let mut aggregated: FxHashMap<u64, FxHashMap<Kmer, f64>> = FxHashMap::default();

    for tile_result in tile_results {
        let Some(path) = tile_result.counts_path else {
            continue;
        };
        let tile_payload = deserialize_tile_counts(&path)?;
        let _ = fs::remove_file(&path);

        for window_counts in tile_payload {
            let entry = aggregated
                .entry(window_counts.original_idx)
                .or_insert_with(FxHashMap::default);
            for count in window_counts.entries {
                let kmer = Kmer {
                    k: count.k,
                    code: count.code,
                    orientation: count.orientation,
                };
                *entry.entry(kmer).or_insert(0.0) += count.value;
            }
        }
    }

    let mut merged: Vec<TileWindowCounts> = aggregated
        .into_iter()
        .map(|(original_idx, counts_map)| {
            let mut entries: Vec<TileKmerCountEntry> = Vec::with_capacity(counts_map.len());
            for (kmer, value) in counts_map {
                entries.push(TileKmerCountEntry {
                    k: kmer.k,
                    code: kmer.code,
                    orientation: kmer.orientation,
                    value,
                });
            }
            TileWindowCounts {
                original_idx,
                entries,
            }
        })
        .collect();

    merged.sort_unstable_by_key(|window| window.original_idx);
    Ok(merged)
}
