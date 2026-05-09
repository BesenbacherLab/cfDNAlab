use crate::shared::length_axis::LengthAxis;
use anyhow::{Context, Result, anyhow, bail, ensure};
use fxhash::FxHashMap;
use ndarray::{Array1, ArrayView1, ArrayView3};
use ndarray_npy::WriteNpyExt;
use ndarray_npy::{NpzReader, WritableElement};
use rayon::prelude::*;
use std::{
    fs::File,
    io::{Cursor, Write},
    ops::Range,
    path::Path,
    sync::{Arc, Mutex},
};
use zip::{CompressionMethod, ZipWriter, write::SimpleFileOptions};

/// Number of `f32` count entries protected by one lock during sparse partial file merging.
///
/// The final output is still one dense vector, but sparse tile partial files are merged in parallel.
/// Splitting the dense vector into lockable chunks lets different workers add into different
/// parts of the output at the same time. At four million entries, each lock protects about sixteen
/// MiB of dense output, which is coarse enough to keep lock overhead low while still
/// allowing useful parallelism when partial files touch different profile regions.
const SPARSE_MERGE_CHUNK_SIZE: usize = 4 * 1024 * 1024;

/// Convert validated profile coordinates to the single flat dense index.
///
/// This is the one canonical formula for midpoint profile storage. Both dense and sparse
/// accumulators call through this helper so a change to axis order cannot silently diverge between
/// in-memory counting and temporary file merging.
///
/// The function assumes that `position`, `group_idx`, and `length_bin_idx` have already been
/// checked against the profile shape. Keeping this helper unchecked avoids repeating bounds checks
/// in the hot path after `LengthAxis` has already mapped a fragment length to a valid bin.
#[inline]
fn profile_flat_index(
    position: usize,
    group_idx: usize,
    length_bin_idx: usize,
    window_size: usize,
    num_length_bins: usize,
) -> usize {
    group_idx * num_length_bins * window_size + length_bin_idx * window_size + position
}

/// Compute the total number of count entries needed for the full dense midpoint profile array.
///
/// This helper exists so every dense allocation uses the same checked
/// `group * length_bin * position` size calculation.
/// Axis order does not change the required length.
fn flattened_profile_size(
    window_size: usize,
    num_groups: usize,
    length_axis: &LengthAxis,
) -> usize {
    num_groups
        .checked_mul(window_size)
        .and_then(|size| size.checked_mul(length_axis.num_bins()))
        .expect("ProfileGroupsCounts shape overflow")
}

/// Validate raw profile coordinates and return the matching flat dense index.
///
/// `length` is a raw fragment length in base pairs, not a length-bin index. The shared
/// `LengthAxis` maps it to the configured half-open bin. The same helper is used by
/// `ProfileGroupsCounts` and `SparseProfileGroupsCounts`, which keeps indexing errors from
/// becoming format-dependent.
///
/// Errors are intentionally explicit because these paths usually indicate a bug in upstream
/// coordinate handling or a mismatch between fragment filtering and length-bin configuration.
fn checked_profile_flat_index(
    position: usize,
    group_idx: usize,
    length: usize,
    window_size: usize,
    num_groups: usize,
    length_axis: &LengthAxis,
) -> Result<usize> {
    if position >= window_size {
        bail!(
            "position {} out of bounds (0..{})",
            position,
            window_size.saturating_sub(1)
        );
    }
    if group_idx >= num_groups {
        bail!(
            "group_idx {} out of bounds (0..{})",
            group_idx,
            num_groups.saturating_sub(1)
        );
    }

    let Some(length_bin_idx) = length_axis.bin_index(length) else {
        bail!(
            "length {} out of allowed half-open range [{}..{})",
            length,
            length_axis.min_fragment_length(),
            length_axis.max_fragment_length() + 1
        );
    };

    Ok(profile_flat_index(
        position,
        group_idx,
        length_bin_idx,
        window_size,
        length_axis.num_bins(),
    ))
}

/// Sparse representation of a midpoint tile partial file after reading it from `.npz`.
///
/// `idx` stores flat dense indices in ascending order and `data` stores the corresponding `f32`
/// counts. The shape is validated while reading and is not stored here because merge code only
/// needs a checked sequence of destination positions and values.
#[derive(Debug)]
pub(crate) struct SparseProfilePartialFile {
    /// Flat indices into the final dense `ProfileGroupsCounts::counts` buffer.
    idx: Vec<u64>,
    /// Count values aligned one-to-one with `idx`.
    data: Vec<f32>,
}

/// Serialize a one-dimensional slice to an in-memory `.npy` payload.
///
/// NumPy `.npz` files are ZIP archives containing named `.npy` files. This helper creates each
/// named array payload before `write_sparse_profile_partial_file_npz` adds it to the ZIP container.
fn vec_to_npy<T: WritableElement>(values: &[T]) -> Result<Vec<u8>> {
    let view: ArrayView1<'_, T> = ArrayView1::from(values);
    let mut buffer = Vec::<u8>::new();
    view.write_npy(Cursor::new(&mut buffer))?;
    Ok(buffer)
}

/// Write a sparse midpoint tile partial file as a NumPy-compatible `.npz` archive.
///
/// The archive contains three arrays:
///
/// - `idx.npy`: `u64` flat dense indices in ascending order
/// - `data.npy`: `f32` counts aligned with `idx`
/// - `shape.npy`: `u64[3]` shape stored as `[group, length_bin, position]`
///
/// The temporary format is designed to reduce tile I/O and avoid materializing
/// dense per-tile arrays during merge, while still using ordinary NumPy containers
/// that can be inspected outside Rust if needed.
fn write_sparse_profile_partial_file_npz(
    path: &Path,
    idx: &[u64],
    data: &[f32],
    shape: &[u64],
) -> Result<()> {
    let idx_npy = vec_to_npy(idx)?;
    let data_npy = vec_to_npy(data)?;
    let shape_npy = vec_to_npy(shape)?;

    let file = File::create(path)
        .with_context(|| format!("creating sparse midpoint partial file {}", path.display()))?;
    let options = SimpleFileOptions::default().compression_method(CompressionMethod::Deflated);
    let mut npz = ZipWriter::new(file);

    npz.start_file("idx.npy", options)?;
    npz.write_all(&idx_npy)?;
    npz.start_file("data.npy", options)?;
    npz.write_all(&data_npy)?;
    npz.start_file("shape.npy", options)?;
    npz.write_all(&shape_npy)?;
    npz.finish()?;

    Ok(())
}

/// Dense midpoint profile accumulator used for final output.
///
/// This type owns the dense count vector that eventually becomes
/// `<prefix>.midpoint_profiles.npy`. During current midpoint runs, per-tile counting uses
/// `SparseProfileGroupsCounts` to reduce temporary file size, then sparse partial files are merged
/// into one `ProfileGroupsCounts` before writing the public dense output.
///
/// # Layout
/// The flattened buffer stores counts in group-major, then length-bin, then position order:
///
/// ```text
/// flat_idx = group_idx * (num_length_bins * window_size)
///     + length_bin_idx * window_size
///     + position
/// ```
///
/// where:
///
/// - `group_idx` is in `[0, num_groups)`
/// - `position` is in `[0, window_size)`
/// - `length_bin_idx` is in `[0, num_length_bins)`
///
/// # Length bins
///
/// Length bins are strictly increasing edges. For example, edges `[20, 50, 100]`
/// create bins `[20, 50)` and `[50, 100)`. A raw fragment length is accepted when it falls inside
/// the half-open range covered by the axis.
#[derive(Debug, Clone)]
pub struct ProfileGroupsCounts {
    /// Flattened dense counts in `(group, length_bin, position)` order.
    pub counts: Vec<f32>,
    /// Number of positions per group profile.
    pub window_size: usize,
    /// Number of grouped BED labels represented in the profile.
    pub num_groups: usize,
    /// Shared fragment length axis used to map raw lengths to output bins.
    length_axis: Arc<LengthAxis>,
}

impl ProfileGroupsCounts {
    /// Create a zero-initialized dense midpoint profile with the requested shape.
    ///
    /// The caller passes the already-resolved `LengthAxis` so counting, merging, and output shape
    /// all use the same bin lookup. The allocation size is
    /// `num_groups * window_size * length_axis.num_bins()`.
    ///
    /// Parameters
    /// ----------
    /// - `window_size`:
    ///     Number of positions per group profile.
    /// - `num_groups`:
    ///     Number of groups from the grouped BED input.
    /// - `length_axis`:
    ///     Shared length-bin lookup used to map raw fragment lengths to bins.
    pub fn new(window_size: usize, num_groups: usize, length_axis: Arc<LengthAxis>) -> Self {
        let flattened_size = flattened_profile_size(window_size, num_groups, &length_axis);

        let counts = vec![0f32; flattened_size];
        Self {
            counts,
            window_size,
            num_groups,
            length_axis,
        }
    }

    /// Return the minimum accepted fragment length.
    ///
    /// This is the first configured length-bin edge and is inclusive.
    #[inline]
    pub fn min_fragment_length(&self) -> u32 {
        self.length_axis.min_fragment_length()
    }

    /// Return the maximum accepted fragment length.
    ///
    /// Length bins are half-open internally, so this is one less than the final configured edge.
    #[inline]
    pub fn max_fragment_length(&self) -> u32 {
        self.length_axis.max_fragment_length()
    }

    /// Return the half-open length-bin edges used by this profile.
    #[inline]
    pub fn length_bins(&self) -> &[u32] {
        self.length_axis.edges()
    }

    /// Compute the flat dense index for a midpoint profile coordinate.
    ///
    /// `length` is an absolute fragment length in base pairs, not a length-bin index. The shared
    /// `LengthAxis` maps it to the configured bin before the flat index formula is applied.
    ///
    /// Returns an error if `position`, `group_idx`, or `length` is outside the profile shape.
    #[inline]
    pub fn index_of(&self, position: usize, group_idx: usize, length: usize) -> Result<usize> {
        checked_profile_flat_index(
            position,
            group_idx,
            length,
            self.window_size,
            self.num_groups,
            &self.length_axis,
        )
    }

    /// Return the current count for a midpoint profile coordinate.
    ///
    /// `length` is a raw fragment length in base pairs.
    #[inline]
    pub fn get(&self, position: usize, group_idx: usize, length: usize) -> Result<f32> {
        let flat_idx = self.index_of(position, group_idx, length)?;
        Ok(self.counts[flat_idx])
    }

    /// Return the number of configured fragment length bins.
    ///
    /// This is one less than the number of length-bin edges.
    #[inline]
    pub fn n_lengths(&self) -> usize {
        self.length_axis.num_bins()
    }

    /// Return the number of positions per group profile.
    #[inline]
    pub fn n_positions(&self) -> usize {
        self.window_size
    }

    /// Return the number of grouped BED labels represented.
    #[inline]
    pub fn n_groups(&self) -> usize {
        self.num_groups
    }

    /// View the dense count vector as a one-dimensional ndarray.
    ///
    /// This preserves the `(group, length_bin, position)` flattening. It is useful for low-level
    /// validation and inspection code that needs the raw vector shape.
    #[inline]
    pub fn as_ndarray1(&self) -> ArrayView1<'_, f32> {
        ArrayView1::from(&self.counts)
    }

    /// View the dense profile as `(group, length_bin, position)` without copying.
    ///
    /// The returned ndarray view is contiguous and has the same axis order as the final `.npy`
    /// output.
    #[inline]
    pub fn view_ndarray3_group_len_pos(&self) -> ArrayView3<'_, f32> {
        let num_groups = self.num_groups;
        let num_length_bins = self.n_lengths();
        let num_positions = self.window_size;
        ArrayView3::from_shape((num_groups, num_length_bins, num_positions), &self.counts)
            .expect("Shape mismatch for (group, length_bin, position)")
    }

    /// Merge sparse midpoint tile partial files into this dense final output buffer.
    ///
    /// Each worker reads one sparse `.npz` temp file, validates its shape, and adds only observed
    /// entries into locked chunks of the final dense vector. No dense per-tile array is allocated
    /// during this merge.
    ///
    /// The public output remains dense. Sparsity is only an internal optimization for temporary
    /// files and merge memory.
    pub fn add_from_sparse_npz_files_parallel<P>(&mut self, paths: Vec<P>) -> Result<()>
    where
        P: AsRef<Path> + Send + Sync,
    {
        self.add_from_sparse_npz_files_parallel_with_chunk_size(paths, SPARSE_MERGE_CHUNK_SIZE)
    }

    /// Merge sparse midpoint partial files with a caller-supplied dense chunk size.
    ///
    /// Production code uses `SPARSE_MERGE_CHUNK_SIZE`. Tests pass smaller chunk sizes to force
    /// multiple chunks in tiny fixtures. The merge owns exactly one final dense vector and borrows
    /// it as mutable chunks protected by separate mutexes.
    ///
    /// The flow is:
    ///
    /// 1. Validate the requested chunk size.
    /// 2. Move `self.counts` out so Rayon workers can borrow disjoint mutable chunks.
    /// 3. Read each sparse temp file in parallel.
    /// 4. Add sparse count entries into locked dense chunks.
    /// 5. Drop chunk locks and move the dense vector back into `self`.
    fn add_from_sparse_npz_files_parallel_with_chunk_size<P>(
        &mut self,
        paths: Vec<P>,
        chunk_size: usize,
    ) -> Result<()>
    where
        P: AsRef<Path> + Send + Sync,
    {
        ensure!(
            chunk_size > 0,
            "sparse midpoint merge chunk size must be greater than zero"
        );
        let expected_shape = self.sparse_shape()?;
        let destination_len = self.counts.len();

        // Borrow the final dense output as lockable chunks, without allocating another dense copy
        let mut dense_counts = std::mem::take(&mut self.counts);
        let chunks: Vec<Mutex<&mut [f32]>> = dense_counts
            .chunks_mut(chunk_size)
            .map(Mutex::new)
            .collect();

        // Spread workers across starting chunks so many broad partial files do not all contend for
        // chunk zero first
        let num_chunks = chunks.len();
        let merge_result =
            paths
                .par_iter()
                .enumerate()
                .try_for_each(|(path_idx, path)| -> Result<()> {
                    let partial_file = read_sparse_profile_partial_file(
                        path.as_ref(),
                        expected_shape,
                        destination_len,
                    )?;
                    let start_chunk = sparse_merge_start_chunk(path_idx, num_chunks);
                    merge_sparse_profile_partial_file(
                        &partial_file,
                        &chunks,
                        chunk_size,
                        start_chunk,
                    )
                    .with_context(|| {
                        format!("merging sparse midpoint partial file {:?}", path.as_ref())
                    })
                });

        // Drop chunk mutexes before moving the dense vector back into `self`
        drop(chunks);
        self.counts = dense_counts;
        merge_result
    }

    /// Return the sparse temp-file shape expected for this dense profile.
    ///
    /// Sparse midpoint partial files store shape as `[group, length_bin, position]`. This lets the
    /// reader reject partial files from a different run before any data are added.
    fn sparse_shape(&self) -> Result<[u64; 3]> {
        Ok([
            u64::try_from(self.num_groups).context("num_groups does not fit in u64")?,
            u64::try_from(self.n_lengths()).context("num length bins does not fit in u64")?,
            u64::try_from(self.window_size).context("window_size does not fit in u64")?,
        ])
    }
}

/// Sparse midpoint accumulator used while counting one tile.
///
/// The hashmap key is the same flat dense index used by `ProfileGroupsCounts`. That choice keeps
/// tile counting simple: the fragment loop can call `incr_weighted` with the same arguments as the
/// dense accumulator, and the merge step can add sparse values directly into the final dense vector
/// without converting coordinates back and forth.
///
/// The sparse accumulator is an internal representation only. It is written to a sparse `.npz`
/// temporary file, then merged into the final dense `.npy` output.
#[derive(Debug, Clone)]
pub(crate) struct SparseProfileGroupsCounts {
    /// Sparse counts keyed by final dense flat index.
    counts: FxHashMap<usize, f32>,
    /// Number of positions per group profile.
    window_size: usize,
    /// Number of grouped BED labels represented in the tile partial file.
    num_groups: usize,
    /// Shared fragment length axis used for raw length to bin lookup.
    length_axis: Arc<LengthAxis>,
}

impl SparseProfileGroupsCounts {
    /// Create an empty sparse accumulator for one tile.
    ///
    /// The accumulator starts with no allocated count entries beyond the hashmap itself. Entries are
    /// inserted only when a fragment midpoint contributes a nonzero weight to that profile
    /// coordinate.
    pub(crate) fn new(window_size: usize, num_groups: usize, length_axis: Arc<LengthAxis>) -> Self {
        Self {
            counts: FxHashMap::default(),
            window_size,
            num_groups,
            length_axis,
        }
    }

    /// Compute the final dense flat index for a sparse midpoint update.
    ///
    /// This uses the shared checked indexing helper so sparse and dense accumulators reject the
    /// same out-of-bounds coordinates.
    #[inline]
    pub(crate) fn index_of(
        &self,
        position: usize,
        group_idx: usize,
        length: usize,
    ) -> Result<usize> {
        checked_profile_flat_index(
            position,
            group_idx,
            length,
            self.window_size,
            self.num_groups,
            &self.length_axis,
        )
    }

    /// Add a weighted count to the sparse tile accumulator.
    ///
    /// Duplicate updates to the same `(group, length_bin, position)` coordinate are summed in
    /// place. The value remains sparse because only observed count entries are stored.
    #[inline]
    pub(crate) fn incr_weighted(
        &mut self,
        position: usize,
        group_idx: usize,
        length: usize,
        weight: f64,
    ) -> Result<()> {
        let flat_idx = self.index_of(position, group_idx, length)?;
        *self.counts.entry(flat_idx).or_insert(0.0) += weight as f32;
        Ok(())
    }

    /// Return whether this tile produced any midpoint counts.
    ///
    /// Empty sparse tiles do not write temp files. Skipping those files reduces merge work and
    /// avoids creating many tiny archives for tiles with no usable fragments or no overlapping
    /// windows.
    #[inline]
    pub(crate) fn is_empty(&self) -> bool {
        self.counts.is_empty()
    }

    /// Write the sparse tile accumulator to a sorted `.npz` partial file.
    ///
    /// Hashmap iteration order is deliberately not serialized. Entries are sorted by flat dense
    /// index before writing, which lets the merge code process each partial file in chunk runs.
    pub(crate) fn write_npz<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let mut entries: Vec<(usize, f32)> = self
            .counts
            .iter()
            .map(|(&flat_idx, &count)| (flat_idx, count))
            .collect();
        entries.sort_unstable_by_key(|(flat_idx, _)| *flat_idx);

        let mut idx = Vec::with_capacity(entries.len());
        let mut data = Vec::with_capacity(entries.len());
        // Store explicit portable indices, not hashmap internals
        for (flat_idx, count) in entries {
            idx.push(u64::try_from(flat_idx).context("sparse profile index does not fit in u64")?);
            data.push(count);
        }

        let shape = [
            u64::try_from(self.num_groups).context("num_groups does not fit in u64")?,
            u64::try_from(self.length_axis.num_bins())
                .context("num length bins does not fit in u64")?,
            u64::try_from(self.window_size).context("window_size does not fit in u64")?,
        ];

        write_sparse_profile_partial_file_npz(path.as_ref(), &idx, &data, &shape)?;

        Ok(())
    }
}

/// Read and validate one sparse midpoint partial file from a `.npz` temp file.
///
/// Validation happens before the partial file is returned so merge code can assume:
///
/// - `idx` and `data` have the same length
/// - `idx` is sorted ascending
/// - all indices fit the current platform `usize`
/// - all indices are inside the destination dense vector
/// - the stored shape matches the current run
///
/// Failing fast here gives a clear file-specific error for corrupt temp files or accidental shape
/// mismatches.
fn read_sparse_profile_partial_file(
    path: &Path,
    expected_shape: [u64; 3],
    destination_len: usize,
) -> Result<SparseProfilePartialFile> {
    let file = File::open(path)
        .with_context(|| format!("opening sparse midpoint partial file {}", path.display()))?;
    let mut npz = NpzReader::new(file)
        .with_context(|| format!("reading sparse midpoint partial file {}", path.display()))?;
    let idx_arr: Array1<u64> = npz
        .by_name("idx.npy")
        .with_context(|| format!("reading idx.npy from {}", path.display()))?;
    let data_arr: Array1<f32> = npz
        .by_name("data.npy")
        .with_context(|| format!("reading data.npy from {}", path.display()))?;
    let shape_arr: Array1<u64> = npz
        .by_name("shape.npy")
        .with_context(|| format!("reading shape.npy from {}", path.display()))?;

    let shape = shape_arr.to_vec();
    ensure!(
        shape.as_slice() == expected_shape,
        "Shape mismatch when merging {}: source shape {:?} vs destination shape {:?}",
        path.display(),
        shape,
        expected_shape
    );

    let idx = idx_arr.to_vec();
    let data = data_arr.to_vec();
    ensure!(
        idx.len() == data.len(),
        "Sparse midpoint partial file {} has {} indices but {} data values",
        path.display(),
        idx.len(),
        data.len()
    );

    // Validate sorted order and bounds once, before the parallel merge starts mutating output
    let mut previous_idx: Option<u64> = None;
    for &flat_idx_u64 in &idx {
        if let Some(previous) = previous_idx {
            ensure!(
                previous <= flat_idx_u64,
                "Sparse midpoint partial file {} indices must be sorted ascending",
                path.display()
            );
        }
        previous_idx = Some(flat_idx_u64);

        let flat_idx = usize::try_from(flat_idx_u64).with_context(|| {
            format!(
                "sparse midpoint partial file {} index {} does not fit in usize",
                path.display(),
                flat_idx_u64
            )
        })?;
        ensure!(
            flat_idx < destination_len,
            "Sparse midpoint partial file {} index {} out of bounds for dense length {}",
            path.display(),
            flat_idx,
            destination_len
        );
    }

    Ok(SparseProfilePartialFile { idx, data })
}

/// Merge one validated sparse partial file into the final dense chunks.
///
/// Sparse indices are sorted by flat dense index. To reduce lock convoys when many partial files
/// touch the same broad set of chunks, each partial file starts at a selected dense chunk and wraps
/// around to the beginning. The merge still visits every sparse entry exactly once.
///
/// `chunks` are mutable slices of the final dense output, each behind its own mutex. Only one chunk
/// lock is held at a time.
fn merge_sparse_profile_partial_file(
    partial_file: &SparseProfilePartialFile,
    chunks: &[Mutex<&mut [f32]>],
    chunk_size: usize,
    start_chunk: usize,
) -> Result<()> {
    ensure!(
        chunk_size > 0,
        "sparse midpoint merge chunk size must be greater than zero"
    );
    if partial_file.idx.is_empty() {
        return Ok(());
    }
    ensure!(
        !chunks.is_empty(),
        "sparse midpoint merge has entries but no dense chunks"
    );
    ensure!(
        start_chunk < chunks.len(),
        "sparse midpoint merge start chunk {} out of bounds for {} chunks",
        start_chunk,
        chunks.len()
    );

    if let Some(&last_flat_idx_u64) = partial_file.idx.last() {
        let last_flat_idx = usize::try_from(last_flat_idx_u64)
            .context("sparse midpoint index does not fit in usize")?;
        let required_chunks = last_flat_idx / chunk_size + 1;
        ensure!(
            required_chunks <= chunks.len(),
            "sparse midpoint index {} maps outside dense chunks",
            last_flat_idx
        );
    }

    // Jump to the first sparse entry at or after the chosen starting chunk
    let start_flat_idx = start_chunk
        .checked_mul(chunk_size)
        .context("sparse midpoint merge start chunk offset overflow")?;
    let start_flat_idx_u64 = u64::try_from(start_flat_idx)
        .context("sparse midpoint merge start chunk offset does not fit in u64")?;
    let split_idx = partial_file
        .idx
        .partition_point(|&flat_idx| flat_idx < start_flat_idx_u64);

    // Merge from the selected chunk to the end, then wrap to the early chunks
    merge_sparse_profile_entry_range(
        partial_file,
        chunks,
        chunk_size,
        split_idx..partial_file.idx.len(),
    )?;
    merge_sparse_profile_entry_range(partial_file, chunks, chunk_size, 0..split_idx)?;

    Ok(())
}

/// Pick the first dense chunk a sparse partial file should attempt to merge.
///
/// File order is used as a simple deterministic source of spread. Starts are restricted to a grid
/// of chunk positions based on the dense chunk count and Rayon thread count. This avoids opaque
/// randomization while still reducing the chance that every worker locks chunk zero first.
fn sparse_merge_start_chunk(path_idx: usize, num_chunks: usize) -> usize {
    if num_chunks <= 1 {
        return 0;
    }

    let num_threads = rayon::current_num_threads().max(1);
    let start_stride = (num_chunks / (num_threads * 2)).max(1);
    let num_candidate_starts = num_chunks.div_ceil(start_stride);
    let candidate_idx = path_idx % num_candidate_starts;

    (candidate_idx * start_stride).min(num_chunks - 1)
}

/// Merge a contiguous range of sparse entries into dense chunks.
///
/// `entry_range` is contiguous in sorted sparse-index order. The function locks the destination
/// chunk for the current entry, consumes all following entries that belong to the same dense chunk,
/// then releases the lock before moving to the next chunk.
fn merge_sparse_profile_entry_range(
    partial_file: &SparseProfilePartialFile,
    chunks: &[Mutex<&mut [f32]>],
    chunk_size: usize,
    entry_range: Range<usize>,
) -> Result<()> {
    let mut entry_idx = entry_range.start;
    while entry_idx < entry_range.end {
        // The first entry in the current run determines which dense chunk must be locked
        let flat_idx = usize::try_from(partial_file.idx[entry_idx])
            .context("sparse midpoint index does not fit in usize")?;
        let chunk_idx = flat_idx / chunk_size;

        let mut chunk = chunks[chunk_idx].lock().map_err(|poisoned| {
            anyhow!("sparse midpoint merge chunk mutex poisoned: {}", poisoned)
        })?;

        while entry_idx < entry_range.end {
            let run_flat_idx = usize::try_from(partial_file.idx[entry_idx])
                .context("sparse midpoint index does not fit in usize")?;
            // Stop before crossing into the next dense chunk, so one lock protects one run
            if run_flat_idx / chunk_size != chunk_idx {
                break;
            }
            let chunk_offset = run_flat_idx % chunk_size;
            chunk[chunk_offset] += partial_file.data[entry_idx];
            entry_idx += 1;
        }
        // `chunk` is dropped here, releasing the mutex before the next chunk run
    }

    Ok(())
}

impl std::fmt::Display for ProfileGroupsCounts {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let length_bins = self.length_bins();
        let len_str = if length_bins.len() > 2 {
            format!(
                "len:[{}..{}...={}]",
                self.min_fragment_length(),
                length_bins[1],
                self.max_fragment_length()
            )
        } else {
            format!(
                "len:[{}..={}]",
                self.min_fragment_length(),
                self.max_fragment_length()
            )
        };
        write!(
            f,
            "ProfileGroupsCounts(groups:[0..={}], {}, pos:[0..={}], size:({}) )",
            self.num_groups - 1,
            len_str,
            self.window_size - 1,
            self.counts.len()
        )
    }
}

#[cfg(test)]
mod tests {
    include!("counting_by_group_tests.rs");
}
