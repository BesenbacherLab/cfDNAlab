use crate::commands::gc_bias::{
    binning::{BinnedAxis, bins_from_edges},
    counting::{GCPrefixes, get_gc_integer_percentage_for_window},
    package::GCCorrectionPackage,
};
use anyhow::{Result, ensure};
use ndarray::Array2;

#[derive(Debug, Clone)]
pub struct GCCorrector {
    correction_matrix: Array2<f64>,
    lengths_bins: BinnedAxis,
    gc_bins: BinnedAxis,
    end_offset: u64,
}

impl GCCorrector {
    pub fn from_package(package: &GCCorrectionPackage) -> Result<Self> {
        let length_bins = bins_from_edges(&package.length_edges)?;
        let gc_bins = bins_from_edges(&package.gc_edges)?;

        Ok(GCCorrector {
            correction_matrix: package.correction_matrix.clone(),
            lengths_bins: length_bins,
            gc_bins: gc_bins,
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
        let length = offset_end - offset_start;
        Ok(Some(self.get_correction_weight(length, gc_bin)?))
    }

    #[inline]
    pub fn get_correction_weight(&self, fragment_length: usize, gc_pct: usize) -> Result<f64> {
        let length_bin = self
            .lengths_bins
            .index_to_bin
            .get(&fragment_length)
            .expect(&format!(
                "GC correction: Observed unexpected fragment length {}",
                fragment_length
            ));

        let gc_bin = self.gc_bins.index_to_bin.get(&gc_pct).expect(&format!(
            "GC correction: Observed unexpected GC percentage {}",
            gc_pct
        ));

        Ok(self.correction_matrix[(*length_bin, *gc_bin)])
    }
}
