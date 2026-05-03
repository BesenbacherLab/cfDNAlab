use crate::shared::length_axis::LengthAxis;
use anyhow::{Context, Result, anyhow, bail, ensure};
use fxhash::FxHashMap;
use ndarray::{Array1, ArrayView1, ArrayView3, ShapeBuilder};
use ndarray_npy::WriteNpyExt;
use ndarray_npy::{NpzReader, WritableElement, read_npy};
use rayon::prelude::*;
use std::{
    fs::File,
    io::{Cursor, Write},
    ops::Range,
    path::Path,
    sync::{Arc, Mutex},
};
use zip::{CompressionMethod, ZipWriter, write::SimpleFileOptions};

const SPARSE_MERGE_CHUNK_SIZE: usize = 4 * 1024 * 1024;

/// Compute the flat index after shape construction and coordinate validation.
#[inline]
fn profile_flat_index(
    position: usize,
    group_idx: usize,
    length_bin_idx: usize,
    window_size: usize,
    num_length_bins: usize,
) -> usize {
    group_idx * window_size * num_length_bins + position * num_length_bins + length_bin_idx
}

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

#[derive(Debug)]
pub(crate) struct SparseProfilePartial {
    idx: Vec<u64>,
    data: Vec<f32>,
}

fn vec_to_npy<T: WritableElement>(values: &[T]) -> Result<Vec<u8>> {
    let view: ArrayView1<'_, T> = ArrayView1::from(values);
    let mut buffer = Vec::<u8>::new();
    view.write_npy(Cursor::new(&mut buffer))?;
    Ok(buffer)
}

fn write_sparse_profile_partial_npz(
    path: &Path,
    idx: &[u64],
    data: &[f32],
    shape: &[u64],
) -> Result<()> {
    let idx_npy = vec_to_npy(idx)?;
    let data_npy = vec_to_npy(data)?;
    let shape_npy = vec_to_npy(shape)?;

    let file = File::create(path)
        .with_context(|| format!("creating sparse midpoint partial {}", path.display()))?;
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

/// Count array for fragment coverage across fragment lengths.
///
/// # Layout
/// The flattened buffer stores counts in **group-major, then position, then length-bin** order:
///
/// ```text
/// flat_idx = group_idx * (window_size * num_length_bins)
///     + position  * (num_length_bins)
///     + length_bin_idx
/// ```
///
/// where:
/// - `group_idx ∈ [0, num_groups)`
/// - `position  ∈ [0, window_size)`
/// - `length_bin_idx ∈ [0, num_length_bins)`
///
/// # Length bins
/// `length_bins` are **edges**, sorted strictly increasing. For example,
/// `length_bins = [20, 50, 100]` creates two bins: `[20,50)` and `[50,100)`.  
/// Lengths are considered **in-bounds** iff `min_edge ≤ length < max_edge` (half-open at the top).
#[derive(Debug, Clone)]
pub struct ProfileGroupsCounts {
    /// 1D vector of flattened counts (see layout above).
    pub counts: Vec<f32>,
    pub window_size: usize,
    pub num_groups: usize,
    /// Shared fragment length axis used to map raw lengths to output bins.
    length_axis: Arc<LengthAxis>,
}

impl ProfileGroupsCounts {
    /// Create a new zero-initialized `ProfileGroupsCounts`.
    ///
    /// * `window_size`: number of positions per profile (e.g., 2001)
    /// * `num_groups`: number of groups (TFs)
    /// * `length_axis`: resolved fragment length axis shared with the counting path
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

    /// Get minimum allowed fragment length (first edge).
    #[inline]
    pub fn min_fragment_length(&self) -> u32 {
        self.length_axis.min_fragment_length()
    }

    /// Get maximum allowed fragment length (last edge `-1`, inclusive).
    #[inline]
    pub fn max_fragment_length(&self) -> u32 {
        self.length_axis.max_fragment_length()
    }

    /// Return the half-open length-bin edges used by this profile.
    #[inline]
    pub fn length_bins(&self) -> &[u32] {
        self.length_axis.edges()
    }

    /// Compute the flattened index for (`position`, `group_idx`, `length`) if all are in bounds.
    ///
    /// - `length` is an **absolute** fragment length (bp). It is binned into the unique
    ///   length-bin `i` such that `edges[i] ≤ length < edges[i+1]`.
    ///
    /// Returns an error if any argument is out of range.
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

    /// Compute the flattened index for (`position`, `group_idx`, `length_bin_idx`)
    /// WITHOUT bound checks.
    ///
    /// NOTE: Takes `length_bin_idx` NOT raw `length`.
    ///
    /// Assumes the arguments are valid.
    fn unchecked_index_of(
        &self,
        position: usize,
        group_idx: usize,
        length_bin_idx: usize,
    ) -> usize {
        let num_length_bins = self.n_lengths();

        profile_flat_index(
            position,
            group_idx,
            length_bin_idx,
            self.window_size,
            num_length_bins,
        )
    }

    /// Increment the counter by `1.0` at (`position`, `group_idx`, `length`), if in bounds.
    ///
    /// Errors if any argument is out of bounds.
    #[inline]
    pub fn incr(&mut self, position: usize, group_idx: usize, length: usize) -> Result<()> {
        let i = self.index_of(position, group_idx, length)?;
        self.counts[i] += 1.0;
        Ok(())
    }

    /// Increment by `weight` at (`position`, `group_idx`, `length`).
    ///
    /// Errors if any argument is out of bounds.
    #[inline]
    pub fn incr_weighted(
        &mut self,
        position: usize,
        group_idx: usize,
        length: usize,
        weight: f64,
    ) -> Result<()> {
        let i = self.index_of(position, group_idx, length)?;
        self.counts[i] += weight as f32;
        Ok(())
    }

    /// Get the count at (`position`, `group_idx`, `length`).
    ///
    /// Errors if any argument is out of bounds.
    #[inline]
    pub fn get(&self, position: usize, group_idx: usize, length: usize) -> Result<f32> {
        let i = self.index_of(position, group_idx, length)?;
        Ok(self.counts[i])
    }

    /// Number of **length bins** (rows in the length dimension).
    #[inline]
    pub fn n_lengths(&self) -> usize {
        self.length_axis.num_bins()
    }

    /// Number of **positions** per profile.
    #[inline]
    pub fn n_positions(&self) -> usize {
        self.window_size
    }

    /// Number of **groups**.
    #[inline]
    pub fn n_groups(&self) -> usize {
        self.num_groups
    }

    /// Create a zero-initialized copy with the same shape and binning.
    #[inline]
    pub fn zeroed_like(&self) -> Self {
        Self {
            counts: vec![0f32; self.counts.len()],
            window_size: self.window_size,
            num_groups: self.num_groups,
            length_axis: Arc::clone(&self.length_axis),
        }
    }

    /// Check if two `ProfileGroupsCounts` are compatible (same shape and bin edges).
    #[inline]
    pub fn is_compatible_with(&self, other: &Self) -> bool {
        self.window_size == other.window_size
            && self.num_groups == other.num_groups
            && self.length_axis.edges() == other.length_axis.edges()
            && self.counts.len() == other.counts.len()
    }

    /// Merge (sum) counts from `other` into `self`.
    ///
    /// Returns an error if shapes or binning are incompatible.
    pub fn merge_from(&mut self, other: &Self) -> anyhow::Result<()> {
        if !self.is_compatible_with(other) {
            anyhow::bail!(
                "incompatible ProfileGroupsCounts: self={} vs other={}",
                self,
                other
            );
        }
        for (dst, src) in self.counts.iter_mut().zip(other.counts.iter()) {
            *dst += *src;
        }
        Ok(())
    }

    /// Collapse (sum) an iterator of `ProfileGroupsCounts` into a single object.
    ///
    /// All inputs must be compatible (same shape and bin edges).
    pub fn collapse<'a, I>(iter: I) -> anyhow::Result<Self>
    where
        I: IntoIterator<Item = &'a ProfileGroupsCounts>,
    {
        let mut it = iter.into_iter();
        let first = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("collapse requires at least one ProfileGroupsCounts"))?;
        let mut acc = first.clone();
        for g in it {
            acc.merge_from(g)?;
        }
        Ok(acc)
    }

    /// Reshape into a 3D vector with axes `(group, length_bin, position)`.
    ///
    /// Note: This **reorders** data from the internal layout `(group, position, length_bin)`
    /// into `(group, length_bin, position)`, so it allocates and copies.
    pub fn to_3d_group_len_pos(&self) -> Vec<Vec<Vec<f32>>> {
        let g = self.num_groups;
        let p = self.window_size;
        let l = self.n_lengths();

        // Allocate [group][length_bin][position]
        let mut out = vec![vec![vec![0f32; p]; l]; g];

        // Read from (group, position, length_bin) and write to (group, length_bin, position)
        for group_idx in 0..g {
            for position in 0..p {
                for len_bin in 0..l {
                    out[group_idx][len_bin][position] =
                        self.counts[self.unchecked_index_of(position, group_idx, len_bin)];
                }
            }
        }
        out
    }

    /// Return a view over the flat counts as an `ndarray::ArrayView1<f32>`.
    ///
    /// Useful for saving a temporary flat vector with `ndarray_npy::write_npy`
    /// and reshaping on load.
    #[inline]
    pub fn as_ndarray1(&self) -> ArrayView1<'_, f32> {
        ArrayView1::from(&self.counts)
    }

    // TODO: Test!!

    /// Zero-copy `ndarray3` view with axes `(group, length_bin, position)`.
    ///
    /// Internally our flat layout is `(group, position, length_bin)` contiguous, so the strides for
    /// `(group, length_bin, position)` are:
    ///   - group stride      = window_size * n_lengths()
    ///   - length_bin stride = 1
    ///   - position stride   = n_lengths()
    #[inline]
    pub fn view_ndarray3_group_len_pos(&self) -> ArrayView3<'_, f32> {
        let g = self.num_groups;
        let l = self.n_lengths();
        let p = self.window_size;
        let group_stride = p
            .checked_mul(l)
            .expect("group stride overflow for (group, length_bin, position)");
        let strides = (group_stride, 1usize, l);
        ArrayView3::from_shape((g, l, p).strides(strides), &self.counts)
            .expect("Shape/stride mismatch for (group, length_bin, position)")
    }

    /// Add counts from a single `.npy` file (1D `f32`) into `self` in place.
    pub fn add_from_npy_1d_file<P: AsRef<Path>>(&mut self, npy_path: P) -> Result<()> {
        // Read as ndarray 1D array
        let src_arr: Array1<f32> = read_npy(&npy_path)
            .with_context(|| format!("Reading npy file {:?}", npy_path.as_ref()))?;
        // Ensure contiguous and get a slice view
        let src = src_arr
            .as_slice()
            .ok_or_else(|| anyhow::anyhow!("NPY array is not contiguous"))?;

        if src.len() != self.counts.len() {
            bail!(
                "Shape mismatch when merging {:?}: src len {} vs dest len {}",
                npy_path.as_ref(),
                src.len(),
                self.counts.len()
            );
        }
        for (dst, s) in self.counts.iter_mut().zip(src.iter()) {
            *dst += *s;
        }
        Ok(())
    }

    /// Sequentially merge a list of `.npy` files (1D `f32`) into `self`.
    pub fn add_from_npy_1d_files_sequential<I, P>(&mut self, paths: I) -> Result<()>
    where
        I: IntoIterator<Item = P>,
        P: AsRef<Path>,
    {
        for p in paths {
            self.add_from_npy_1d_file(p)?;
        }
        Ok(())
    }

    /// Parallel merge of `.npy` files (1D `f32`) into `self` using Rayon.
    ///
    /// Each worker loads one file into memory, then accumulates into a shared buffer under a mutex.
    pub fn add_from_npy_1d_files_parallel<P>(&mut self, paths: Vec<P>) -> Result<()>
    where
        P: AsRef<Path> + Send + Sync,
    {
        // Move destination counts behind a Mutex for shared accumulation.
        let acc = Arc::new(Mutex::new(std::mem::take(&mut self.counts)));
        let dest_len = acc
            .lock()
            .map_err(|poisoned| anyhow!("profile-group accumulator mutex poisoned: {}", poisoned))?
            .len();

        paths.par_iter().try_for_each(|p| -> Result<()> {
            // Read as 1D ndarray
            let src_arr: Array1<f32> =
                read_npy(p).with_context(|| format!("Reading npy file {:?}", p.as_ref()))?;
            let src = src_arr
                .as_slice()
                .ok_or_else(|| anyhow::anyhow!("NPY array is not contiguous"))?;
            if src.len() != dest_len {
                bail!(
                    "Shape mismatch when merging {:?}: src len {} vs dest len {}",
                    p.as_ref(),
                    src.len(),
                    dest_len
                );
            }
            // Accumulate while holding the lock
            let mut guard = acc.lock().map_err(|poisoned| {
                anyhow!(
                    "profile-group accumulator mutex poisoned while merging {:?}: {}",
                    p.as_ref(),
                    poisoned
                )
            })?;
            for (dst, s) in guard.iter_mut().zip(src.iter()) {
                *dst += *s;
            }
            Ok(())
        })?;

        // Restore the accumulated vector
        self.counts = Arc::try_unwrap(acc)
            .map_err(|_| anyhow!("profile-group accumulator still has multiple references"))?
            .into_inner()
            .map_err(|poisoned| {
                anyhow!("profile-group accumulator mutex poisoned: {}", poisoned)
            })?;

        Ok(())
    }

    /// Parallel merge of sparse midpoint tile partials into this dense final buffer.
    ///
    /// Each worker reads one sparse temp file and writes only observed cells into locked chunks of
    /// the dense output. This avoids materializing one dense array per temp file during merge.
    pub fn add_from_sparse_npz_files_parallel<P>(&mut self, paths: Vec<P>) -> Result<()>
    where
        P: AsRef<Path> + Send + Sync,
    {
        self.add_from_sparse_npz_files_parallel_with_chunk_size(paths, SPARSE_MERGE_CHUNK_SIZE)
    }

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
        let dest_len = self.counts.len();

        let mut dense_counts = std::mem::take(&mut self.counts);
        let chunks: Vec<Mutex<&mut [f32]>> = dense_counts
            .chunks_mut(chunk_size)
            .map(Mutex::new)
            .collect();

        let num_chunks = chunks.len();
        let merge_result =
            paths
                .par_iter()
                .enumerate()
                .try_for_each(|(path_idx, path)| -> Result<()> {
                    let partial =
                        read_sparse_profile_partial(path.as_ref(), expected_shape, dest_len)?;
                    let start_chunk = sparse_merge_start_chunk(path_idx, num_chunks);
                    merge_sparse_profile_partial(&partial, &chunks, chunk_size, start_chunk)
                        .with_context(|| {
                            format!("merging sparse midpoint partial {:?}", path.as_ref())
                        })
                });

        drop(chunks);
        self.counts = dense_counts;
        merge_result
    }

    fn sparse_shape(&self) -> Result<[u64; 3]> {
        Ok([
            u64::try_from(self.num_groups).context("num_groups does not fit in u64")?,
            u64::try_from(self.window_size).context("window_size does not fit in u64")?,
            u64::try_from(self.n_lengths()).context("num length bins does not fit in u64")?,
        ])
    }
}

#[derive(Debug, Clone)]
pub(crate) struct SparseProfileGroupsCounts {
    counts: FxHashMap<usize, f32>,
    window_size: usize,
    num_groups: usize,
    length_axis: Arc<LengthAxis>,
}

impl SparseProfileGroupsCounts {
    /// Create a sparse midpoint-profile accumulator for one tile.
    ///
    /// The key is the same flat index used by `ProfileGroupsCounts`, so sparse tile partials can
    /// be merged directly into the final dense buffer without coordinate conversion.
    pub(crate) fn new(window_size: usize, num_groups: usize, length_axis: Arc<LengthAxis>) -> Self {
        Self {
            counts: FxHashMap::default(),
            window_size,
            num_groups,
            length_axis,
        }
    }

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

    #[inline]
    pub(crate) fn is_empty(&self) -> bool {
        self.counts.is_empty()
    }

    pub(crate) fn write_npz<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let mut entries: Vec<(usize, f32)> = self
            .counts
            .iter()
            .map(|(&flat_idx, &count)| (flat_idx, count))
            .collect();
        entries.sort_unstable_by_key(|(flat_idx, _)| *flat_idx);

        let mut idx = Vec::with_capacity(entries.len());
        let mut data = Vec::with_capacity(entries.len());
        for (flat_idx, count) in entries {
            idx.push(u64::try_from(flat_idx).context("sparse profile index does not fit in u64")?);
            data.push(count);
        }

        let shape = [
            u64::try_from(self.num_groups).context("num_groups does not fit in u64")?,
            u64::try_from(self.window_size).context("window_size does not fit in u64")?,
            u64::try_from(self.length_axis.num_bins())
                .context("num length bins does not fit in u64")?,
        ];

        write_sparse_profile_partial_npz(path.as_ref(), &idx, &data, &shape)?;

        Ok(())
    }
}

fn read_sparse_profile_partial(
    path: &Path,
    expected_shape: [u64; 3],
    dest_len: usize,
) -> Result<SparseProfilePartial> {
    let file = File::open(path)
        .with_context(|| format!("opening sparse midpoint partial {}", path.display()))?;
    let mut npz = NpzReader::new(file)
        .with_context(|| format!("reading sparse midpoint partial {}", path.display()))?;
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
        "Shape mismatch when merging {}: src shape {:?} vs dest shape {:?}",
        path.display(),
        shape,
        expected_shape
    );

    let idx = idx_arr.to_vec();
    let data = data_arr.to_vec();
    ensure!(
        idx.len() == data.len(),
        "Sparse midpoint partial {} has {} indices but {} data values",
        path.display(),
        idx.len(),
        data.len()
    );

    let mut previous_idx: Option<u64> = None;
    for &flat_idx_u64 in &idx {
        if let Some(previous) = previous_idx {
            ensure!(
                previous <= flat_idx_u64,
                "Sparse midpoint partial {} indices must be sorted ascending",
                path.display()
            );
        }
        previous_idx = Some(flat_idx_u64);

        let flat_idx = usize::try_from(flat_idx_u64).with_context(|| {
            format!(
                "sparse midpoint partial {} index {} does not fit in usize",
                path.display(),
                flat_idx_u64
            )
        })?;
        ensure!(
            flat_idx < dest_len,
            "Sparse midpoint partial {} index {} out of bounds for dense length {}",
            path.display(),
            flat_idx,
            dest_len
        );
    }

    Ok(SparseProfilePartial { idx, data })
}

fn merge_sparse_profile_partial(
    partial: &SparseProfilePartial,
    chunks: &[Mutex<&mut [f32]>],
    chunk_size: usize,
    start_chunk: usize,
) -> Result<()> {
    ensure!(
        chunk_size > 0,
        "sparse midpoint merge chunk size must be greater than zero"
    );
    if partial.idx.is_empty() {
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

    if let Some(&last_flat_idx_u64) = partial.idx.last() {
        let last_flat_idx = usize::try_from(last_flat_idx_u64)
            .context("sparse midpoint index does not fit in usize")?;
        let required_chunks = last_flat_idx / chunk_size + 1;
        ensure!(
            required_chunks <= chunks.len(),
            "sparse midpoint index {} maps outside dense chunks",
            last_flat_idx
        );
    }

    let start_flat_idx = start_chunk
        .checked_mul(chunk_size)
        .context("sparse midpoint merge start chunk offset overflow")?;
    let start_flat_idx_u64 = u64::try_from(start_flat_idx)
        .context("sparse midpoint merge start chunk offset does not fit in u64")?;
    let split_idx = partial
        .idx
        .partition_point(|&flat_idx| flat_idx < start_flat_idx_u64);

    merge_sparse_profile_entry_range(partial, chunks, chunk_size, split_idx..partial.idx.len())?;
    merge_sparse_profile_entry_range(partial, chunks, chunk_size, 0..split_idx)?;

    Ok(())
}

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

fn merge_sparse_profile_entry_range(
    partial: &SparseProfilePartial,
    chunks: &[Mutex<&mut [f32]>],
    chunk_size: usize,
    entry_range: Range<usize>,
) -> Result<()> {
    let mut entry_idx = entry_range.start;
    while entry_idx < entry_range.end {
        let flat_idx = usize::try_from(partial.idx[entry_idx])
            .context("sparse midpoint index does not fit in usize")?;
        let chunk_idx = flat_idx / chunk_size;

        let mut chunk = chunks[chunk_idx].lock().map_err(|poisoned| {
            anyhow!("sparse midpoint merge chunk mutex poisoned: {}", poisoned)
        })?;

        while entry_idx < entry_range.end {
            let run_flat_idx = usize::try_from(partial.idx[entry_idx])
                .context("sparse midpoint index does not fit in usize")?;
            if run_flat_idx / chunk_size != chunk_idx {
                break;
            }
            let chunk_offset = run_flat_idx % chunk_size;
            chunk[chunk_offset] += partial.data[entry_idx];
            entry_idx += 1;
        }
        // `chunk` is dropped here, releasing the mutex before the next chunk run.
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
