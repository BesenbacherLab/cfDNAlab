use crate::{
    commands::gc_bias::{
        binning::{BinnedAxis, compute_bin_edges},
        load_reference_bias::ReferenceGCMetadata,
    },
    shared::{
        constants::{GC_CORRECTION_SCHEMA_VERSION, MIN_ACGT_BASES_FOR_GC_FRACTION},
        reference::ContigFootprintEntry,
        zarr::{
            create_zarr_store, ensure_zarr_schema, read_zarr_array1, read_zarr_array2,
            read_zarr_root_attributes, write_single_chunk_zarr_array, write_zarr_root_metadata,
        },
    },
};
use anyhow::{Context, Result, bail, ensure};
use ndarray::{Array1, Array2};
use serde_json::{Value, json};
use std::{path::Path, sync::Arc};
use zarrs::{array::data_type, filesystem::FilesystemStore};

const GC_CORRECTION_SCHEMA: &str = "gc_correction_package";

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
    pub(crate) fn from_components(
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

    pub fn write_zarr<P: AsRef<std::path::Path>>(&self, path: P) -> Result<()> {
        validate_correction_matrix_for_writing(
            &self.correction_matrix,
            &self.length_edges,
            &self.gc_edges,
        )?;
        ensure!(
            self.length_bin_frequencies.len() == self.correction_matrix.nrows(),
            "Length frequency length ({}) must match number of correction rows ({})",
            self.length_bin_frequencies.len(),
            self.correction_matrix.nrows()
        );

        let store = create_zarr_store(path.as_ref(), "GC correction package")?;
        write_zarr_root_metadata(
            store.clone(),
            "GC correction package",
            json!({
                "cfdnalab_schema": GC_CORRECTION_SCHEMA,
                "cfdnalab_schema_version": self.version,
                "package_role": "sample_gc_correction",
                "correction_units": "multiplicative_fragment_weight",
                "gc_percent_rounding": "integer_half_up",
                "minimum_acgt_bases_for_gc_fraction": MIN_ACGT_BASES_FOR_GC_FRACTION,
                "end_offset": self.end_offset,
            }),
        )?;

        let correction_matrix = self.correction_matrix.as_standard_layout();
        let correction_matrix_slice = correction_matrix
            .as_slice()
            .context("correction_matrix should be contiguous after layout conversion")?;
        write_single_chunk_zarr_array(
            store.clone(),
            "correction_matrix",
            &[
                self.correction_matrix.nrows(),
                self.correction_matrix.ncols(),
            ],
            &["length_bin", "gc_bin"],
            correction_matrix_slice,
            data_type::float64(),
            0.0,
            json!({"long_name": "multiplicative GC correction matrix"}),
        )?;
        write_single_chunk_zarr_array(
            store.clone(),
            "length_edges",
            &[self.length_edges.len()],
            &["length_edge"],
            &self.length_edges,
            data_type::uint32(),
            0u32,
            json!({"long_name": "half-open fragment length bin edges"}),
        )?;
        write_single_chunk_zarr_array(
            store.clone(),
            "gc_edges",
            &[self.gc_edges.len()],
            &["gc_edge"],
            &self.gc_edges,
            data_type::uint32(),
            0u32,
            json!({"long_name": "half-open GC percent bin edges"}),
        )?;
        let length_bin_frequencies = self.length_bin_frequencies.as_standard_layout();
        let length_bin_frequencies_slice = length_bin_frequencies
            .as_slice()
            .context("length_bin_frequencies should be contiguous after layout conversion")?;
        write_single_chunk_zarr_array(
            store.clone(),
            "length_bin_frequencies",
            &[self.length_bin_frequencies.len()],
            &["length_bin"],
            length_bin_frequencies_slice,
            data_type::float64(),
            0.0,
            json!({"long_name": "normalized length-bin frequency weights"}),
        )?;
        let reference_contig_footprint_json = serde_json::to_vec(&self.reference_contig_footprint)?;
        write_single_chunk_zarr_array(
            store,
            "reference_contig_footprint_json",
            &[reference_contig_footprint_json.len()],
            &["json_byte"],
            &reference_contig_footprint_json,
            data_type::uint8(),
            0u8,
            json!({"long_name": "JSON-encoded reference contig footprint"}),
        )?;
        Ok(())
    }

    pub fn from_file<P: AsRef<std::path::Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        validate_zarr_correction_package_path(path)?;
        let root = read_zarr_root_attributes(path)?;
        ensure_zarr_schema(
            &root,
            GC_CORRECTION_SCHEMA,
            GC_CORRECTION_SCHEMA_VERSION,
            "GC correction package",
        )?;

        let store = Arc::new(FilesystemStore::new(path)?);
        let correction_matrix = read_zarr_array2::<f64>(store.clone(), "/correction_matrix")?;
        let length_edges = read_zarr_array1::<u32>(store.clone(), "/length_edges")?;
        let gc_edges = read_zarr_array1::<u32>(store.clone(), "/gc_edges")?;
        let length_bin_frequencies_vec =
            read_zarr_array1::<f64>(store.clone(), "/length_bin_frequencies")?;
        let reference_contig_footprint_json =
            read_zarr_array1::<u8>(store, "/reference_contig_footprint_json")?;
        let reference_contig_footprint: Vec<ContigFootprintEntry> =
            serde_json::from_slice(&reference_contig_footprint_json)
                .context("invalid reference_contig_footprint_json in GC correction package")?;
        let version = root
            .get("cfdnalab_schema_version")
            .and_then(Value::as_u64)
            .context("GC correction package is missing cfdnalab_schema_version")?;
        let version = u32::try_from(version)
            .context("GC correction package schema version must fit in u32")?;
        let end_offset = root
            .get("end_offset")
            .and_then(Value::as_u64)
            .context("GC correction package is missing end_offset")?;
        let length_bin_frequencies = Array1::from(length_bin_frequencies_vec);

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
            length_bin_frequencies.len() == correction_matrix.dim().0,
            "Length frequency length ({}) must match number of correction rows ({})",
            length_bin_frequencies.len(),
            correction_matrix.dim().0
        );
        validate_correction_matrix_values(&correction_matrix, &length_edges, &gc_edges)?;

        Ok(Self {
            version,
            end_offset,
            length_edges,
            gc_edges,
            correction_matrix,
            length_bin_frequencies,
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

fn validate_zarr_correction_package_path(path: &Path) -> Result<()> {
    if !path.is_dir() {
        bail!(
            "GC correction package path must point to an existing .zarr directory: {}",
            path.display()
        );
    }

    let has_zarr_extension = path
        .extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| extension.eq_ignore_ascii_case("zarr"));

    ensure!(
        has_zarr_extension,
        "GC correction package path must point to a .zarr directory: {}",
        path.display()
    );

    Ok(())
}
