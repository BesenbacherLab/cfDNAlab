use crate::shared::kmers::kmer_codec::KmerSpec;
use crate::shared::kmers::process_counts::DecodedCounts;
use anyhow::{Context, Result};
use fxhash::FxHashMap;
use ndarray::{Array2, ArrayView1, arr1};
use ndarray_npy::WriteNpyExt; // trait brings .write_npy into scope
use ndarray_npy::{WritableElement, write_npy};
use num_traits::NumCast;
use std::fs::File;
use std::io::Cursor;
use std::io::Write;
use std::path::Path;
use zip::{ZipWriter, write::SimpleFileOptions};

/// Write one `.npy` matrix and a companion `*_motifs.txt` file for every
/// k present in `prepared_windows`.
///
/// * `prepared_windows` – windows of decoded counts.
/// * `kmer_specs`       – validated specs: the keys determine which k values
///                        will be written, and in which order.
/// * `output_dir`       – target directory.
///
/// * For reference windows the files are named  `k<k>_counts.npy`, e.g.
///   `k3_counts.npy`.  
///
/// The matrix dimensions are **windows × motifs** with the same column order
/// used across all windows of that k-mer size.
pub fn write_decoded_counts_matrix(
    prepared_windows: &[DecodedCounts],
    kmer_specs: &FxHashMap<u8, KmerSpec>,
    motifs_by_k: &FxHashMap<u8, Vec<String>>,
    output_dir: &Path,
    prefix: &str,
    save_sparse: bool,
) -> anyhow::Result<()> {
    let n_win = prepared_windows.len();

    for &k in kmer_specs.keys() {
        // Collect reference bins for this k
        let mut ref_bins: Vec<FxHashMap<String, f64>> = vec![FxHashMap::default(); n_win];
        for (idx, win) in prepared_windows.iter().enumerate() {
            if let Some(bin) = win.counts.get(&k) {
                ref_bins[idx] = bin.clone();
            }
        }
        let tag = format!("{}.k{}", prefix, k);
        if save_sparse {
            write_category_sparse(&ref_bins, &motifs_by_k[&k], &tag, output_dir)?;
        } else {
            write_category(&ref_bins, &motifs_by_k[&k], &tag, output_dir)?;
        }
    }

    Ok(())
}

/// Write <prefix>_counts.npy and <prefix>_motifs.txt
///
/// * `motifs`  - The motifs to include for all bins in the order you want it saved in.
fn write_category(
    bins: &[FxHashMap<String, f64>],
    motifs: &[String],
    prefix: &str,
    out_dir: &Path,
) -> anyhow::Result<()> {
    if bins.is_empty() {
        return Ok(()); // nothing to write
    }

    // Output matrix
    let n_rows = bins.len();
    let n_cols = motifs.len();
    let mut mat = Array2::<f64>::zeros((n_rows, n_cols));

    // Pre-compute motif -> column index once
    let col_of: FxHashMap<_, _> = motifs.iter().enumerate().map(|(c, m)| (m, c)).collect();

    for (row, hm) in bins.iter().enumerate() {
        for (motif, &cnt) in hm {
            if let Some(&col) = col_of.get(motif) {
                mat[(row, col)] = cnt; // Counts overwrite the zero
            }
        }
    }

    // Persist outputs
    write_npy(out_dir.join(format!("{prefix}_counts.npy")), &mat)?;

    let mut txt = File::create(out_dir.join(format!("{prefix}_motifs.txt")))?;
    for m in motifs {
        writeln!(txt, "{m}")?;
    }

    Ok(())
}

// Sparse version

type Idx = u64; // 64-bit row and column indices

/// Write SciPy-compatible COO matrix as <prefix>_counts_sparse.npz + <prefix>_motifs.txt
///
/// * `bins`   – Per-bin motif->count hash maps
/// * `motifs` – Full ordered motif list; defines column order
///
/// Examples
/// --------
///
/// Python: Loading with scipy.sparse and looping through each window:
///
/// ```python,ignore
/// coo = scipy.sparse.load_npz("k11_counts_sparse.npz")`
/// csr = coo.tocsr()
/// for i in range(csr.shape[0]):
///     sparse_row = csr.getrow(i)
///     dense_row = sparse_row.toarray().ravel()
///     ...
///
/// with open("my_prefix_motifs.txt") as f:
///     motifs = [line.strip() for line in f]
/// ```
pub fn write_category_sparse(
    bins: &[FxHashMap<String, f64>],
    motifs: &[String],
    prefix: &str,
    out_dir: &Path,
) -> Result<()> {
    if bins.is_empty() {
        return Ok(());
    }

    let n_rows = bins.len();
    let n_cols = motifs.len();

    // Motif --> column lookup
    let motif_index: FxHashMap<&str, Idx> = motifs
        .iter()
        .enumerate()
        .map(|(i, m)| (m.as_str(), i as Idx))
        .collect();

    // Collect triplets with one allocation
    let nnz: usize = bins.iter().map(|hm| hm.len()).sum();
    let mut row = Vec::<Idx>::with_capacity(nnz);
    let mut col = Vec::<Idx>::with_capacity(nnz);
    let mut val = Vec::<f64>::with_capacity(nnz);

    for (r, hm) in bins.iter().enumerate() {
        let ri: Idx = NumCast::from(r).context("row index overflow u64")?;
        for (motif, &count) in hm {
            if let Some(&ci) = motif_index.get(motif.as_str()) {
                row.push(ri);
                col.push(ci);
                val.push(count);
            }
        }
    }

    // Serialise numeric vectors
    let row_npy = vec_to_npy(&row)?;
    let col_npy = vec_to_npy(&col)?;
    let val_npy = vec_to_npy(&val)?;

    // shape = np.array([n_rows, n_cols], dtype=int64)
    let shape_arr = arr1(&[n_rows as i64, n_cols as i64]);
    let mut shape_buf = Vec::<u8>::new();
    shape_arr.write_npy(Cursor::new(&mut shape_buf))?;

    // format = np.array('coo', dtype='|S3')
    let format_buf = numpy_string_scalar("coo")?;

    // Pack everything into <prefix>_counts_sparse.npz
    let npz_path = out_dir.join(format!("{prefix}_counts_sparse.npz"));
    let file = File::create(&npz_path)?;
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

    // Plain-text motif list
    let mut txt = File::create(out_dir.join(format!("{prefix}_motifs.txt")))?;
    for m in motifs {
        writeln!(txt, "{m}")?;
    }

    Ok(())
}

// Vec --> .npy buffer helper
fn vec_to_npy<T: WritableElement>(v: &[T]) -> Result<Vec<u8>> {
    let view: ArrayView1<'_, T> = ArrayView1::from(v);
    let mut buf = Vec::<u8>::new();
    view.write_npy(Cursor::new(&mut buf))?;
    Ok(buf)
}

// Builds a scalar string .npy with dtype '|S{len}'
fn numpy_string_scalar(s: &str) -> Result<Vec<u8>> {
    let bytes = s.as_bytes();
    let len = bytes.len();
    let header_body = format!("{{'descr': '|S{len}', 'fortran_order': False, 'shape': (), }}",);
    let mut header = header_body.into_bytes();
    header.push(b'\n');

    // Pad header so that (10 + header_len) % 16 == 0
    let mut header_len = header.len();
    let magic_len = 6 + 2 + 2; // \x93NUMPY + ver + hdr_len field
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
