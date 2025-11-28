use crate::commands::gc_bias::{
    GC_CORRECTION_SCHEMA_VERSION,
    binning::{BinnedAxis, compute_bin_edges},
    load_reference_bias::ReferenceGCMetadata,
};
use anyhow::{Context, Result, ensure};
use ndarray::{Array1, Array2};
use ndarray_npy::{NpzReader, NpzWriter};
use std::fs::File;

#[derive(Clone, Debug)]
pub struct GCCorrectionPackage {
    pub version: u32,
    pub end_offset: u64,
    pub length_edges: Vec<u32>,
    pub gc_edges: Vec<u32>,
    pub correction_matrix: Array2<f64>,
    pub length_bin_frequencies: Array1<f64>,
}

impl GCCorrectionPackage {
    pub fn from_components(
        version: u32,
        length_bins: &BinnedAxis,
        gc_bins: &BinnedAxis,
        correction_matrix: Array2<f64>,
        length_bin_frequencies: Array1<f64>,
        reference_metadata: &ReferenceGCMetadata,
    ) -> Result<Self> {
        let length_edges = compute_bin_edges(
            length_bins,
            reference_metadata.min_fragment_length as u32,
            reference_metadata.max_fragment_length as u32,
        )?;
        let gc_edges = compute_bin_edges(gc_bins, 0, 100)?;
        Ok(Self {
            version,
            end_offset: reference_metadata.end_offset as u64,
            length_edges,
            gc_edges,
            correction_matrix,
            length_bin_frequencies,
        })
    }

    pub fn write_npz<P: AsRef<std::path::Path>>(&self, path: P) -> Result<()> {
        let file = File::create(path)?;
        let mut npz = NpzWriter::new(file);
        npz.add_array("correction_matrix", &self.correction_matrix)?;
        npz.add_array("length_edges", &Array1::from(self.length_edges.clone()))?;
        npz.add_array("gc_edges", &Array1::from(self.gc_edges.clone()))?;
        npz.add_array("version", &Array1::from(vec![self.version]))?;
        npz.add_array("end_offset", &Array1::from(vec![self.end_offset]))?;
        npz.add_array(
            "length_bin_frequencies",
            &Array1::from(self.length_bin_frequencies.clone()),
        )?;
        npz.finish()?;
        Ok(())
    }

    pub fn from_file<P: AsRef<std::path::Path>>(path: P) -> Result<Self> {
        let file = File::open(path.as_ref())
            .with_context(|| format!("opening correction package {}", path.as_ref().display()))?;
        let mut reader = NpzReader::new(file)?;

        let correction_matrix: Array2<f64> = reader.by_name("correction_matrix")?;
        let length_edges_arr: Array1<u32> = reader.by_name("length_edges")?;
        let gc_edges_arr: Array1<u32> = reader.by_name("gc_edges")?;
        let version_arr: Array1<u32> = reader.by_name("version")?;
        let end_offset_arr: Array1<u64> = reader.by_name("end_offset")?;
        let length_bin_frequencies_arr: Array1<f64> = reader.by_name("length_bin_frequencies")?;

        let version = *version_arr
            .iter()
            .next()
            .context("version array in GC correction package is empty")?;
        ensure!(
            version == GC_CORRECTION_SCHEMA_VERSION,
            "GC correction package schema version mismatch: file={}, expected={}; \
            Incompatible with this version of cfDNAlab.",
            version,
            GC_CORRECTION_SCHEMA_VERSION
        );
        let end_offset = *end_offset_arr
            .iter()
            .next()
            .context("end_offset array in GC correction package is empty")?;

        let length_edges = length_edges_arr.to_vec();
        let gc_edges = gc_edges_arr.to_vec();

        ensure!(
            length_edges.len() == correction_matrix.dim().0 + 1,
            "Number of Length edges ({}) must match number of correction rows + 1 ({})",
            length_edges.len(),
            correction_matrix.dim().0 + 1
        );
        ensure!(
            gc_edges.len() == correction_matrix.dim().1 + 1,
            "Number of GC edges ({}) must match number of correction columns + 1 ({})",
            gc_edges.len(),
            correction_matrix.dim().1 + 1
        );
        ensure!(
            length_bin_frequencies_arr.len() == correction_matrix.dim().0,
            "Length frequency length ({}) must match number of correction rows ({})",
            length_bin_frequencies_arr.len(),
            correction_matrix.dim().0
        );

        Ok(Self {
            version,
            end_offset,
            length_edges,
            gc_edges,
            correction_matrix,
            length_bin_frequencies: length_bin_frequencies_arr,
        })
    }
}
