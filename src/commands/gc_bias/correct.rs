use crate::commands::gc_bias::{
    binning::{BinnedAxis, bins_from_edges},
    counting::{GCPrefixes, get_gc_integer_percentage_for_window},
    package::GCCorrectionPackage,
};
use anyhow::{Context, Result, anyhow, ensure};
use ndarray::Array2;

#[derive(Debug, Clone)]
pub struct GCCorrector {
    correction_matrix: Array2<f64>,
    lengths_bins: BinnedAxis,
    gc_bins: BinnedAxis,
    length_min: usize,
    length_max: usize,
    gc_min: usize,
    gc_max: usize,
    end_offset: u64,
}

impl GCCorrector {
    pub fn from_package(package: &GCCorrectionPackage) -> Result<Self> {
        let length_bins = bins_from_edges(&package.length_edges)?;
        let gc_bins = bins_from_edges(&package.gc_edges)?;
        let length_min = *package
            .length_edges
            .first()
            .expect("GC correction package contained no length edges");
        let length_max = *package
            .length_edges
            .last()
            .expect("GC correction package contained no length edges");
        let gc_min = *package
            .gc_edges
            .first()
            .expect("GC correction package contained no GC edges");
        let gc_max = *package
            .gc_edges
            .last()
            .expect("GC correction package contained no GC edges");

        Ok(GCCorrector {
            correction_matrix: package.correction_matrix.clone(),
            lengths_bins: length_bins,
            gc_bins: gc_bins,
            length_min: length_min as usize,
            length_max: length_max as usize,
            gc_min: gc_min as usize,
            gc_max: gc_max as usize,
            end_offset: package.end_offset.clone(),
        })
    }
    /// Get GC correction weight from fragment coordinates and tile-/chromosome-wise prefix arrays
    ///
    /// **NOTE**: Coordinates must be relative to the prefix arrays.
    #[inline]
    pub fn correct_fragment(
        &self,
        start: u64,
        end: u64,
        gc_prefixes: &GCPrefixes,
    ) -> Result<Option<f64>> {
        let fragment_length = end.checked_sub(start).expect("fragment end precedes start") as usize;
        let offset_start = start.saturating_add(self.end_offset) as usize;
        let offset_end = end.saturating_sub(self.end_offset) as usize;
        ensure!(
            offset_end > offset_start,
            "GC correction: After applying end-offsets the fragment has no bases left to count GCs at.\
            Does the minimum fragment length match the one in the reference bias?"
        );
        let gc_bin = if let Some(gc_pct) =
            get_gc_integer_percentage_for_window(gc_prefixes, offset_start, offset_end, 0.0, 10)
        {
            gc_pct
        } else {
            return Ok(None);
        };
        Ok(Some(self.get_correction_weight(fragment_length, gc_bin)?))
    }

    /// Get the GC correction weight for a combination of fragment length and GC percentage.
    ///
    /// NOTE: The weight is **multiplicative**, so to correct a fragment's contribution,
    /// **multiply** its existing weight (e.g. `1.0`) with the correction weight.
    #[inline]
    pub fn get_correction_weight(&self, fragment_length: usize, gc_pct: usize) -> Result<f64> {
        // The length index and GC index has the minimum values assigned at 0
        // So we offset the values by their minimum to get the index, which
        // in turn will give us the bin (for which we have correction factors)
        let length_idx = fragment_length - self.length_min;
        let gc_idx = gc_pct - self.gc_min;
        let length_bin = self
            .lengths_bins
            .index_to_bin
            .get(&length_idx)
            .copied()
            .ok_or_else(|| {
                anyhow!(
                    "GC correction: unexpected fragment length {}",
                    fragment_length
                )
            })
            .with_context(|| format!("length range [{}-{}]", self.length_min, self.length_max))?;

        let gc_bin = self
            .gc_bins
            .index_to_bin
            .get(&gc_idx)
            .copied()
            .ok_or_else(|| anyhow!("GC correction: unexpected GC percentage {}", gc_pct))
            .with_context(|| format!("GC range [{}-{}]", self.gc_min, self.gc_max))?;

        Ok(self.correction_matrix[(length_bin, gc_bin)])
    }
}
