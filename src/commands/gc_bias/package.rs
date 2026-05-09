use crate::{
    commands::gc_bias::{
        binning::{BinnedAxis, compute_bin_edges},
        load_reference_bias::ReferenceGCMetadata,
    },
    shared::{constants::GC_CORRECTION_SCHEMA_VERSION, reference::ContigFootprintEntry},
};
use anyhow::{Context, Result, bail, ensure};
use ndarray::{Array1, Array2};
use ndarray_npy::{NpzReader, NpzWriter};
use std::{fs::File, path::Path};

#[derive(Clone, Debug)]
pub struct GCCorrectionPackage {
    pub version: u32,
    pub end_offset: u64,
    pub length_edges: Vec<u32>,
    pub gc_edges: Vec<u32>,
    pub correction_matrix: Array2<f64>,
    pub length_bin_frequencies: Array1<f64>,
    pub reference_contig_footprint: Vec<ContigFootprintEntry>,
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
        validate_correction_matrix_for_writing(&correction_matrix, &length_edges, &gc_edges)?;
        Ok(Self {
            version,
            end_offset: reference_metadata.end_offset as u64,
            length_edges,
            gc_edges,
            correction_matrix,
            length_bin_frequencies,
            reference_contig_footprint: reference_metadata.reference_contig_footprint.clone(),
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
        npz.add_array(
            "reference_contig_footprint_json",
            &Array1::from(serde_json::to_vec(&self.reference_contig_footprint)?),
        )?;
        npz.finish()?;
        Ok(())
    }

    pub fn from_file<P: AsRef<std::path::Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        validate_correction_package_path(path)?;

        let file = File::open(path)
            .with_context(|| format!("opening correction package {}", path.display()))?;
        let mut reader = NpzReader::new(file).with_context(|| {
            format!(
                "reading GC correction package {} as a .npz archive",
                path.display()
            )
        })?;
        let array_names = reader.names()?;

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
        ensure!(
            array_names
                .iter()
                .any(|name| name == "reference_contig_footprint_json"),
            "Missing reference_contig_footprint_json in GC correction package. Rebuild the package with the current schema."
        );
        let reference_contig_footprint_json: Array1<u8> =
            reader.by_name("reference_contig_footprint_json")?;
        let reference_contig_footprint: Vec<ContigFootprintEntry> = serde_json::from_slice(
            reference_contig_footprint_json
                .as_slice()
                .context("reference_contig_footprint_json should be contiguous")?,
        )
        .context("invalid reference_contig_footprint_json in GC correction package")?;

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
        validate_correction_matrix_values(&correction_matrix, &length_edges, &gc_edges)?;

        Ok(Self {
            version,
            end_offset,
            length_edges,
            gc_edges,
            correction_matrix,
            length_bin_frequencies: length_bin_frequencies_arr,
            reference_contig_footprint,
        })
    }
}

fn validate_correction_matrix_for_writing(
    correction_matrix: &Array2<f64>,
    length_edges: &[u32],
    gc_edges: &[u32],
) -> Result<()> {
    ensure!(
        length_edges.len() == correction_matrix.nrows() + 1,
        "Number of Length edges ({}) must match number of correction rows + 1 ({})",
        length_edges.len(),
        correction_matrix.nrows() + 1
    );
    ensure!(
        gc_edges.len() == correction_matrix.ncols() + 1,
        "Number of GC edges ({}) must match number of correction columns + 1 ({})",
        gc_edges.len(),
        correction_matrix.ncols() + 1
    );

    validate_correction_matrix_values(correction_matrix, length_edges, gc_edges)?;

    Ok(())
}

fn validate_correction_matrix_values(
    correction_matrix: &Array2<f64>,
    length_edges: &[u32],
    gc_edges: &[u32],
) -> Result<()> {
    for ((length_bin_idx, gc_bin_idx), &weight) in correction_matrix.indexed_iter() {
        ensure!(
            weight.is_finite() && weight >= 0.0,
            "GC correction matrix contains invalid weight {} at length bin {} [{}-{}], GC bin {} [{}-{}]. Correction weights must be finite and non-negative",
            weight,
            length_bin_idx,
            length_edges[length_bin_idx],
            length_edges[length_bin_idx + 1],
            gc_bin_idx,
            gc_edges[gc_bin_idx],
            gc_edges[gc_bin_idx + 1]
        );
    }

    Ok(())
}

fn validate_correction_package_path(path: &Path) -> Result<()> {
    if !path.is_file() {
        bail!(
            "GC correction package path must point to an existing .npz file: {}",
            path.display()
        );
    }

    let has_npz_extension = path
        .extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| extension.eq_ignore_ascii_case("npz"));

    ensure!(
        has_npz_extension,
        "GC correction package path must point to a .npz file: {}",
        path.display()
    );

    Ok(())
}
