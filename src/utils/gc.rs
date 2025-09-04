use std::collections::HashMap;

use crate::cfdna_utils::fragment::{MinimalReadInfo, collect_fragment};
use anyhow::Context;
use rust_htslib::bam::Record;

/// Prefix sums (cumsum) to compute GC and AT fractions while excluding Ns.
/// `gc[i]`   = # of G/C in seq[0..i)
/// `acgt[i]`= # of A/T/G/C in seq[0..i)
pub struct GCPrefixes {
    pub gc: Vec<u32>,
    pub acgt: Vec<u32>,
}

/// Build prefix-sums (cumsum) for GC and ACGT (non-N) counts over a byte slice.
/// This lets you compute GC% on A/T/G/C only, so AT% = 1 - GC%.
///
/// Ignores Ns.
pub fn build_gc_prefixes(seq: &[u8]) -> GCPrefixes {
    let mut gc = Vec::with_capacity(seq.len() + 1);
    let mut acgt = Vec::with_capacity(seq.len() + 1);
    gc.push(0);
    acgt.push(0);

    for &b in seq {
        let is_gc = matches!(b, b'G' | b'g' | b'C' | b'c') as u32;
        let is_acgt = matches!(b, b'A' | b'a' | b'T' | b't' | b'G' | b'g' | b'C' | b'c') as u32;

        gc.push(gc.last().copied().unwrap() + is_gc);
        acgt.push(acgt.last().copied().unwrap() + is_acgt);
    }

    GCPrefixes { gc, acgt }
}

/// Compute the GC fraction for a window [start, end), excluding 'N's.
///
/// `min_acgt`: Minimum number of actual ACGT bases counted in the window.
///   E.g. if most of the window is blacklisted or Ns.
///
/// Returns None if the window has no A/T/G/C bases.
#[inline]
pub fn get_gc_fraction_in_window(
    prefixes: &GCPrefixes,
    start: usize,
    end: usize,
    min_acgt: u32,
) -> Option<f32> {
    debug_assert!(
        start < end && end <= prefixes.gc.len() - 1,
        "GC window [{}, {}) out of bounds (len={})",
        start,
        end,
        prefixes.gc.len() - 1
    );
    let gc = prefixes.gc[end] - prefixes.gc[start];
    let acgt = prefixes.acgt[end] - prefixes.acgt[start];
    if acgt == 0 || acgt < min_acgt as u32 {
        return None;
    }
    let gc_frac = gc as f32 / acgt as f32;
    Some(gc_frac)
}

// TODO: Generalize to any (Within-chromosome) region fetch? Requires knowing start and end? And user should submit subsetted ref_seq

/// reader.fetch(tid as u32, 0, ref_len)?;
/// let v = count_gc_in_chromosome(reader.records(), &ref_seq, by_len, |r| true)?;
///
// pub fn count_gc_in_chromosome<I, F>(
//     records: I,
//     ref_seq: &[u8],
//     length_range: Option<(u32, u32)>,
//     min_acgt: u32,
//     mut filter: F,
// ) -> anyhow::Result<Vec<Vec<u32>>>
// where
//     I: Iterator<Item = rust_htslib::errors::Result<Record>>,
//     F: FnMut(&Record) -> bool,
// {
//     let gc_prefix = build_gc_prefixes(ref_seq);

//     let (min_l, max_l, num_lengths) = if let Some((min_l, max_l)) = length_range {
//         (min_l, Some(max_l), (max_l - min_l + 1) as usize)
//     } else {
//         (0, None, 1usize)
//     };

//     let mut counts = vec![vec![0u32; 101]; num_lengths];

//     let mut stash: HashMap<Vec<u8>, MinimalReadInfo> = HashMap::new();

//     for rec_res in records {
//         let rec = rec_res.context("reading bam record")?;
//         if !filter(&rec) {
//             continue;
//         }

//         if let Some(mate) = stash.remove(rec.qname()) {
//             // Extract fragment
//             let fragment = if let Some(f) = collect_fragment(&MinimalReadInfo::from(&rec), &mate) {
//                 f
//             } else {
//                 continue;
//             };

//             debug_assert!(
//                 (fragment.end as usize) <= ref_seq.len(),
//                 "fragment end {} exceeds reference len {}",
//                 fragment.end,
//                 ref_seq.len()
//             );

//             // Check length is within allowed range
//             let fragment_length = fragment.len();
//             if fragment_length < min_l {
//                 continue;
//             }
//             if let Some(max_length) = max_l {
//                 if fragment_length > max_length {
//                     continue;
//                 }
//             }

//             let gc = get_gc_fraction_in_window(
//                 &gc_prefix,
//                 fragment.start as usize,
//                 fragment.end as usize,
//                 min_acgt,
//             );

//             let gc = match gc {
//                 Some(v) => v,
//                 None => continue,
//             };

//             assert!(gc.is_finite(), "GC not finite: {}", gc);
//             assert!(
//                 (0.0..=1.0).contains(&gc),
//                 "GC fraction out of [0,1]: {}",
//                 gc
//             );

//             // Choose length bin
//             let len_idx = if max_l.is_some() {
//                 (fragment_length - min_l) as usize
//             } else {
//                 0
//             };

//             // clamp GC bin to [0,100]
//             let gc_bin = (gc * 100.0).floor() as usize;

//             // Count!
//             counts[len_idx][gc_bin] += 1;
//         } else {
//             // Left-most read → stash
//             stash.insert(rec.qname().to_vec(), MinimalReadInfo::from(&rec));
//         }
//     }

//     Ok(counts)
// }

/// Count matrix for fragment coverage across GC fraction bins and fragment lengths.
///
/// The matrix is two-dimensional:
/// - Rows correspond to fragment lengths.
/// - Columns correspond to GC fraction bins (0–100).
pub struct GCCounts {
    counts: Vec<Vec<u32>>,
    gc_min: usize,
    gc_max: usize,
    length_min: usize,
    length_max: usize,
}

impl GCCounts {
    /// Create a new `GCCounts` with specified ranges and binning.
    ///
    /// Parameters
    /// ----------
    /// gc_min: usize
    ///     Minimum GC bin (inclusive).
    /// gc_max: usize
    ///     Maximum GC bin (inclusive).
    /// length_min: usize
    ///     Minimum fragment length (inclusive).
    /// length_max: usize
    ///     Maximum fragment length (inclusive).
    ///
    /// Returns
    /// -------
    /// counts: GCCounts
    ///     A `GCCounts` object with all counts initialized to zero.
    pub fn new(gc_min: usize, gc_max: usize, length_min: usize, length_max: usize) -> Self {
        let num_gc_bins = gc_max - gc_min + 1;
        let num_lengths = length_max - length_min + 1;
        let counts = vec![vec![0u32; num_gc_bins]; num_lengths];
        Self {
            counts,
            gc_min,
            gc_max,
            length_min,
            length_max,
        }
    }

    /// Check whether `(length, gc)` is within configured ranges.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length to test.
    /// gc: usize
    ///     GC bin to test.
    ///
    /// Returns
    /// -------
    /// ok: bool
    ///     True if both indices are in range.
    #[inline]
    fn in_bounds(&self, length: usize, gc: usize) -> bool {
        (self.length_min..=self.length_max).contains(&length)
            && (self.gc_min..=self.gc_max).contains(&gc)
    }

    /// Compute row/column indices for `(length, gc)` if in bounds.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length (absolute, not zero-based).
    /// gc: usize
    ///     GC bin (absolute, not zero-based).
    ///
    /// Returns
    /// -------
    /// idx: Option<(usize, usize)>
    ///     `(row, col)` zero-based indices if in range, otherwise `None`.
    #[inline]
    pub fn index_of(&self, length: usize, gc: usize) -> Option<(usize, usize)> {
        if self.in_bounds(length, gc) {
            Some((length - self.length_min, gc - self.gc_min))
        } else {
            None
        }
    }

    /// Increment the counter for a given fragment length and GC bin.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length (absolute).
    /// gc: usize
    ///     GC bin (absolute).
    pub fn incr(&mut self, length: usize, gc: usize) {
        if let Some((r, c)) = self.index_of(length, gc) {
            self.counts[r][c] = self.counts[r][c].saturating_add(1);
        }
    }

    // Get the count at a given fragment length and GC bin.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length (absolute).
    /// gc: usize
    ///     GC bin (absolute).
    ///
    /// Returns
    /// -------
    /// count: Option<u32>
    ///     The count if indices are in range, otherwise `None`.
    pub fn get(&self, length: usize, gc: usize) -> Option<u32> {
        self.index_of(length, gc).map(|(r, c)| self.counts[r][c])
    }

    /// Number of length rows.
    ///
    /// Returns
    /// -------
    /// n: usize
    ///     Count of length bins.
    #[inline]
    pub fn n_lengths(&self) -> usize {
        self.length_max - self.length_min + 1
    }

    /// Number of GC columns.
    ///
    /// Returns
    /// -------
    /// n: usize
    ///     Count of GC bins.
    #[inline]
    pub fn n_gc_bins(&self) -> usize {
        self.gc_max - self.gc_min + 1
    }
}

impl Default for GCCounts {
    /// Create an empty default `GCCounts` (0–100 GC, 20–600 length).
    fn default() -> Self {
        Self::new(0, 100, 20, 600)
    }
}
