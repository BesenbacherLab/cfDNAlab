use crate::commands::fragment_kmers::tiling::PositionDescriptor;
use crate::shared::io::dot_join;
use crate::shared::kmers::{kmer_codec::KmerSpec, process_counts::DecodedCounts};
use crate::shared::positioning::PositionGroup;
use anyhow::{Context, Result};
use fxhash::FxHashMap;
use ndarray::Array3;
use ndarray_npy::{WritableElement, WriteNpyExt, write_npy};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io::{Cursor, Write},
    path::{Path, PathBuf},
};
use zip::{ZipWriter, write::SimpleFileOptions};

/// Persist positional k-mer counts grouped by orientation/endpoint.
///
/// Dense matrices are written as `.npy` files with shape `(windows, positions, motifs)`.
/// Sparse outputs flatten the `(window, position)` axes into a single dimension; a companion
/// `*_grid.txt` file records the original `(windows, positions)` layout. Example for loading a
/// sparse matrix in Python:
///
/// ```python,ignore
/// import numpy as np
/// import scipy.sparse
///
/// coo = scipy.sparse.load_npz("PREFIX.k3_left_counts_sparse.npz")
/// grid = dict(line.strip().split() for line in open("PREFIX.k3_left_grid.txt"))
/// windows = int(grid["windows"])
/// positions = int(grid["positions"])
/// dense = coo.toarray().reshape(windows, positions, coo.shape[1])
/// motifs = [line.strip() for line in open("PREFIX.k3_left_motifs.txt")]
/// offsets = [int(line.strip()) for line in open("PREFIX.left_positions.txt")]
/// ```
pub(crate) fn write_positional_output(
    positional_counts: &[FxHashMap<PositionDescriptor, DecodedCounts>],
    motifs_by_k: &FxHashMap<u8, Vec<String>>,
    kmer_specs: &FxHashMap<u8, KmerSpec>,
    output_dir: &Path,
    prefix: &str,
    save_sparse: bool,
) -> Result<Vec<PathBuf>> {
    let mut written_paths = Vec::new();
    let n_windows = positional_counts.len();

    let positions_by_group = collect_positions_by_group(positional_counts);
    let position_paths = write_positions_metadata(prefix, output_dir, &positions_by_group)?;
    written_paths.extend(position_paths);

    // Iterate k's in a deterministic order
    let mut ks: Vec<u8> = kmer_specs.keys().copied().collect();
    ks.sort_unstable();

    for k in ks {
        let motifs = motifs_by_k
            .get(&k)
            .with_context(|| format!("missing motif list for k={k}"))?;
        let motif_index: FxHashMap<&str, usize> = motifs
            .iter()
            .enumerate()
            .map(|(idx, motif)| (motif.as_str(), idx))
            .collect();

        for (group, offsets) in &positions_by_group {
            let label = group_label(*group);
            if offsets.is_empty() {
                continue;
            }

            if save_sparse {
                let paths = write_positional_sparse_matrix(
                    k,
                    label,
                    offsets,
                    *group,
                    positional_counts,
                    &motif_index,
                    motifs,
                    output_dir,
                    prefix,
                )?;
                written_paths.extend(paths);
            } else {
                let paths = write_positional_dense_matrix(
                    k,
                    label,
                    offsets,
                    *group,
                    positional_counts,
                    &motif_index,
                    motifs,
                    output_dir,
                    prefix,
                    n_windows,
                )?;
                written_paths.extend(paths);
            }
        }
    }

    Ok(written_paths)
}

fn collect_positions_by_group(
    positional_counts: &[FxHashMap<PositionDescriptor, DecodedCounts>],
) -> BTreeMap<PositionGroup, Vec<i32>> {
    let mut positions: BTreeMap<PositionGroup, BTreeSet<i32>> = BTreeMap::new();

    for window in positional_counts {
        for descriptor in window.keys() {
            positions
                .entry(descriptor.group)
                .or_default()
                .insert(descriptor.offset);
        }
    }

    positions
        .into_iter()
        .map(|(group, set)| (group, set.into_iter().collect()))
        .collect()
}

fn write_positions_metadata(
    prefix: &str,
    output_dir: &Path,
    positions_by_group: &BTreeMap<PositionGroup, Vec<i32>>,
) -> Result<Vec<PathBuf>> {
    let mut written_paths = Vec::new();
    for (group, offsets) in positions_by_group {
        if offsets.is_empty() {
            continue;
        }
        let path = output_dir.join(format!("{}.{}_positions.txt", prefix, group_label(*group)));
        let mut file =
            File::create(&path).with_context(|| format!("creating {}", path.display()))?;
        for offset in offsets {
            writeln!(file, "{offset}")?;
        }
        written_paths.push(path);
    }
    Ok(written_paths)
}

fn write_positional_dense_matrix(
    k: u8,
    label: &str,
    offsets: &[i32],
    group: PositionGroup,
    positional_counts: &[FxHashMap<PositionDescriptor, DecodedCounts>],
    motif_index: &FxHashMap<&str, usize>,
    motifs: &[String],
    output_dir: &Path,
    prefix: &str,
    n_windows: usize,
) -> Result<Vec<PathBuf>> {
    let n_positions = offsets.len();
    let n_motifs = motifs.len();
    let mut array = Array3::<f64>::zeros((n_windows, n_positions, n_motifs));

    for (window_idx, window_counts) in positional_counts.iter().enumerate() {
        for (pos_idx, offset) in offsets.iter().enumerate() {
            let descriptor = PositionDescriptor {
                group,
                offset: *offset,
            };
            if let Some(decoded) = window_counts.get(&descriptor)
                && let Some(k_counts) = decoded.counts.get(&k)
            {
                for (motif, value) in k_counts {
                    if let Some(&col_idx) = motif_index.get(motif.as_str()) {
                        array[[window_idx, pos_idx, col_idx]] = *value;
                    }
                }
            }
        }
    }

    let counts_path = output_dir.join(dot_join(&[prefix, &format!("k{k}_{label}_counts.npy")]));
    write_npy(&counts_path, &array)?;

    let motifs_path = write_motif_list(prefix, k, label, motifs, output_dir)?;

    Ok(vec![counts_path, motifs_path])
}

fn write_positional_sparse_matrix(
    k: u8,
    label: &str,
    offsets: &[i32],
    group: PositionGroup,
    positional_counts: &[FxHashMap<PositionDescriptor, DecodedCounts>],
    motif_index: &FxHashMap<&str, usize>,
    motifs: &[String],
    output_dir: &Path,
    prefix: &str,
) -> Result<Vec<PathBuf>> {
    let n_windows = positional_counts.len();
    let n_positions = offsets.len();
    let n_motifs = motifs.len();

    let mut row = Vec::<u64>::new();
    let mut col = Vec::<u64>::new();
    let mut val = Vec::<f64>::new();

    for (w_idx, window_counts) in positional_counts.iter().enumerate() {
        for (p_idx, offset) in offsets.iter().enumerate() {
            let descriptor = PositionDescriptor {
                group,
                offset: *offset,
            };
            let Some(decoded) = window_counts.get(&descriptor) else {
                continue;
            };
            let Some(k_counts) = decoded.counts.get(&k) else {
                continue;
            };
            let row_index = (w_idx * n_positions + p_idx) as u64;
            for (motif, value) in k_counts {
                if let Some(&col_idx) = motif_index.get(motif.as_str()) {
                    row.push(row_index);
                    col.push(col_idx as u64);
                    val.push(*value);
                }
            }
        }
    }

    let row_npy = vec_to_npy(&row)?;
    let col_npy = vec_to_npy(&col)?;
    let val_npy = vec_to_npy(&val)?;
    let shape_arr = ndarray::arr1(&[(n_windows * n_positions) as i64, n_motifs as i64]);
    let mut shape_buf = Vec::<u8>::new();
    shape_arr.write_npy(Cursor::new(&mut shape_buf))?;
    let format_buf = numpy_string_scalar("coo")?;

    let npz_path = output_dir.join(dot_join(&[
        prefix,
        &format!("k{k}_{label}_counts_sparse.npz"),
    ]));
    let file =
        File::create(&npz_path).with_context(|| format!("creating {}", npz_path.display()))?;
    let mut npz = ZipWriter::new(file);
    let opts = SimpleFileOptions::default().compression_method(zip::CompressionMethod::Deflated);

    npz.start_file("row.npy", opts)?;
    npz.write_all(&row_npy)?;
    npz.start_file("col.npy", opts)?;
    npz.write_all(&col_npy)?;
    npz.start_file("data.npy", opts)?;
    npz.write_all(&val_npy)?;
    npz.start_file("shape.npy", opts)?;
    npz.write_all(&shape_buf)?;
    npz.start_file("format.npy", opts)?;
    npz.write_all(&format_buf)?;
    npz.finish()?;

    let motifs_path = write_motif_list(prefix, k, label, motifs, output_dir)?;
    let grid_path = write_grid_metadata(prefix, k, label, n_windows, n_positions, output_dir)?;

    Ok(vec![npz_path, motifs_path, grid_path])
}

fn write_motif_list(
    prefix: &str,
    k: u8,
    label: &str,
    motifs: &[String],
    output_dir: &Path,
) -> Result<PathBuf> {
    let path = output_dir.join(dot_join(&[prefix, &format!("k{k}_{label}_motifs.txt")]));
    let mut file = File::create(&path).with_context(|| format!("creating {}", path.display()))?;
    for motif in motifs {
        writeln!(file, "{motif}")?;
    }
    Ok(path)
}

fn write_grid_metadata(
    prefix: &str,
    k: u8,
    label: &str,
    windows: usize,
    positions: usize,
    output_dir: &Path,
) -> Result<PathBuf> {
    let path = output_dir.join(dot_join(&[prefix, &format!("k{k}_{label}_grid.txt")]));
    let mut file = File::create(&path).with_context(|| format!("creating {}", path.display()))?;
    writeln!(file, "windows\t{}", windows)?;
    writeln!(file, "positions\t{}", positions)?;
    Ok(path)
}

fn group_label(group: PositionGroup) -> &'static str {
    match group {
        PositionGroup::Left => "left",
        PositionGroup::Right => "right",
        PositionGroup::Mid => "mid",
    }
}

fn vec_to_npy<T: WritableElement>(v: &[T]) -> Result<Vec<u8>> {
    let view = ndarray::ArrayView1::from(v);
    let mut buf = Vec::<u8>::new();
    view.write_npy(Cursor::new(&mut buf))?;
    Ok(buf)
}

fn numpy_string_scalar(s: &str) -> Result<Vec<u8>> {
    let bytes = s.as_bytes();
    let len = bytes.len();
    let header_body = format!("{{'descr': '|S{len}', 'fortran_order': False, 'shape': (), }}");
    let mut header = header_body.into_bytes();
    header.push(b'\n');

    let mut header_len = header.len();
    let magic_len = 6 + 2 + 2;
    let pad = (16 - ((magic_len + header_len) % 16)) % 16;
    header.splice(header_len - 1..header_len - 1, vec![b' '; pad]);
    header_len += pad;

    let mut buf = Vec::<u8>::with_capacity(magic_len + header_len + len);
    buf.extend_from_slice(b"\x93NUMPY\x01\x00");
    buf.extend(&(header_len as u16).to_le_bytes());
    buf.extend_from_slice(&header);
    buf.extend_from_slice(bytes);
    Ok(buf)
}
