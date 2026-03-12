use ndarray_npy::read_npy;
use rayon::prelude::*;
use std::{
    path::Path,
    sync::{Arc, Mutex},
};

use anyhow::{Context, Result, bail};
use ndarray::{Array1, ArrayView1, ArrayView3, ShapeBuilder};

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
    /// Sorted bin **edges** (len >= 2), strictly increasing.
    pub length_bins: Vec<u32>,
    /// Fast mapping from absolute length (bp) -> length_bin_idx,
    /// or usize::MAX if below lower bound.
    length_to_bin: Vec<usize>,
}

impl ProfileGroupsCounts {
    /// Create a new zero-initialized `ProfileGroupsCounts`.
    ///
    /// * `window_size`: number of positions per profile (e.g., 2001)
    /// * `num_groups`: number of groups (TFs)
    /// * `length_bins`: sorted length **edges** (len >= 2), strictly increasing
    pub fn new(window_size: usize, num_groups: usize, length_bins: Vec<u32>) -> Self {
        debug_assert!(
            length_bins.len() >= 2,
            "length_bins must have at least 2 edges"
        );
        debug_assert!(
            length_bins.windows(2).all(|w| w[0] < w[1]),
            "length_bins must be strictly increasing"
        );

        // Length -> Length bin lookup vector (no hashing)
        let length_to_bin = Self::precompute_length_bin_lookup(&length_bins);

        let num_length_bins = length_bins.len() - 1; // Edges -> bins
        let flattened_size = num_groups
            .checked_mul(window_size)
            .and_then(|x| x.checked_mul(num_length_bins))
            .expect("ProfileGroupsCounts shape overflow");

        let counts = vec![0f32; flattened_size];
        Self {
            counts,
            window_size,
            num_groups,
            length_bins,
            length_to_bin,
        }
    }

    #[inline]
    fn precompute_length_bin_lookup(length_bins: &[u32]) -> Vec<usize> {
        // Precompute length -> bin lookup for O(1) binning.
        // Fill for [min_edge, max_edge), others remain usize::MAX (invalid).
        let num_length_bins = length_bins.len() - 1; // Edges -> bins
        let max_edge = *length_bins.last().expect("length_bins non-empty") as usize; // Exclusive
        let mut length_to_bin = vec![usize::MAX; max_edge.max(1)]; // Avoid zero length vec

        for bin_idx in 0..num_length_bins {
            let lo = length_bins[bin_idx] as usize;
            let hi = length_bins[bin_idx + 1] as usize; // Exclusive
            for l in lo..hi {
                length_to_bin[l] = bin_idx;
            }
        }
        length_to_bin
    }

    /// Get minimum allowed fragment length (first edge).
    #[inline]
    pub fn min_fragment_length(&self) -> u32 {
        self.length_bins[0]
    }

    /// Get maximum allowed fragment length (last edge `-1`, inclusive).
    #[inline]
    pub fn max_fragment_length(&self) -> u32 {
        *self.length_bins.last().unwrap() - 1
    }

    /// Check bounds for `position`, `group_idx`, and `length`. Returns informative error.
    #[inline]
    fn check_bounds(&self, position: usize, group_idx: usize, length: usize) -> Result<()> {
        if position >= self.window_size {
            bail!(
                "position {} out of bounds (0..{})",
                position,
                self.window_size.saturating_sub(1)
            );
        }
        if group_idx >= self.num_groups {
            bail!(
                "group_idx {} out of bounds (0..{})",
                group_idx,
                self.num_groups.saturating_sub(1)
            );
        }
        let l = length as u32;
        if !(self.min_fragment_length() <= l && l <= self.max_fragment_length()) {
            bail!(
                "length {} out of allowed half-open range [{}..{})",
                length,
                self.min_fragment_length(),
                self.max_fragment_length() + 1
            );
        }
        Ok(())
    }

    /// Compute the flattened index for (`position`, `group_idx`, `length`) if all are in bounds.
    ///
    /// - `length` is an **absolute** fragment length (bp). It is binned into the unique
    ///   length-bin `i` such that `edges[i] ≤ length < edges[i+1]`.
    ///
    /// Returns `Some(idx)` on success, otherwise `None` if any argument is out of range.
    #[inline]
    pub fn index_of(&self, position: usize, group_idx: usize, length: usize) -> Result<usize> {
        self.check_bounds(position, group_idx, length)?;
        let len = length;

        // Map absolute length to bin in O(1) via LUT.
        if len >= self.length_to_bin.len() {
            bail!(
                "length {} out of LUT range [0..{}) (edges last={} exclusive)",
                length,
                self.length_to_bin.len(),
                self.max_fragment_length() + 1
            );
        }
        let length_bin_idx = self.length_to_bin[len];
        if length_bin_idx == usize::MAX {
            bail!(
                "length {} does not fall into any configured bin (edges are half-open; last is exclusive)",
                length
            );
        }

        let num_length_bins = self.n_lengths();
        let idx = group_idx * self.window_size * num_length_bins
            + position * num_length_bins
            + length_bin_idx;

        Ok(idx)
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

        group_idx * self.window_size * num_length_bins + position * num_length_bins + length_bin_idx
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

    //// Get the count at (`position`, `group_idx`, `length`).
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
        self.length_bins.len() - 1
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
            length_bins: self.length_bins.clone(),
            length_to_bin: self.length_to_bin.clone(),
        }
    }

    /// Check if two `ProfileGroupsCounts` are compatible (same shape and bin edges).
    #[inline]
    pub fn is_compatible_with(&self, other: &Self) -> bool {
        self.window_size == other.window_size
            && self.num_groups == other.num_groups
            && self.length_bins == other.length_bins
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
        let strides = ((p * l), 1usize, l);
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
        let dest_len = acc.lock().unwrap().len();

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
            let mut guard = acc.lock().unwrap();
            for (dst, s) in guard.iter_mut().zip(src.iter()) {
                *dst += *s;
            }
            Ok(())
        })?;

        // Restore the accumulated vector
        self.counts = Arc::try_unwrap(acc)
            .expect("Accumulator still has multiple references")
            .into_inner()
            .expect("Mutex poisoned");

        Ok(())
    }
}

impl std::fmt::Display for ProfileGroupsCounts {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let len_str = if self.length_bins.len() > 2 {
            format!(
                "len:[{}..{}...={}]",
                self.min_fragment_length(),
                self.length_bins[1],
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
