use super::*;
use crate::shared::{
    bed::GroupedWindows, kmers::kmer_codec::build_kmer_specs, reference::ContigFootprintEntry,
};
use serde_json::Value;
use std::sync::Arc;
use tempfile::TempDir;
use zarrs::{array::Array, filesystem::FilesystemStore};

fn kmer_spec(kmer_size: u8) -> KmerSpec {
    build_kmer_specs(&[kmer_size])
        .expect("k-mer spec should build")
        .remove(&kmer_size)
        .expect("requested k-mer spec should exist")
}

fn footprint() -> Vec<ContigFootprintEntry> {
    vec![ContigFootprintEntry {
        name: "chr1".to_string(),
        size: 100,
    }]
}

fn validate_sparse_global_package(
    motif_labels: &[String],
    column_kind: SelectedMotifColumnKind,
    kmer_size: u8,
    canonical: bool,
) -> Result<()> {
    let frequency_bins: Vec<RefKmerFrequencyBin> = vec![FxHashMap::default()];
    let scaling_factors = vec![0.0];
    let reference_contig_footprint = footprint();
    validate_ref_kmer_package(&RefKmerZarrPackage {
        frequency_bins: &frequency_bins,
        row_scaling_factors: &scaling_factors,
        motif_labels,
        column_kind,
        row_metadata: RefKmerRowMetadata::Global,
        write_dense_output: false,
        kmer_size,
        canonical,
        all_motifs: false,
        assign_by: WindowAssigner::Any,
        reference_contig_footprint: &reference_contig_footprint,
    })
}

fn assert_close(observed: f64, expected: f64) {
    assert!(
        (observed - expected).abs() < 1e-12,
        "observed {observed}, expected {expected}"
    );
}

#[test]
fn normalize_count_bins_writes_frequencies_and_scaling_factors() -> Result<()> {
    // Arrange: row 0 has total count 3, so its frequencies are 2/3 and 1/3. Row 1 has no
    // counts and therefore keeps scaling factor 0 with no stored frequencies.
    let count_bins = vec![
        FxHashMap::from_iter([(0, 2.0), (1, 1.0)]),
        FxHashMap::default(),
    ];

    // Act
    let normalized = normalize_count_bins_to_frequencies(count_bins)?;

    // Assert
    assert_eq!(normalized.row_scaling_factors, vec![3.0, 0.0]);
    assert_close(normalized.frequency_bins[0][&0], 2.0 / 3.0);
    assert_close(normalized.frequency_bins[0][&1], 1.0 / 3.0);
    assert!(normalized.frequency_bins[1].is_empty());
    assert_close(
        normalized.frequency_bins[0][&0] * normalized.row_scaling_factors[0],
        2.0,
    );
    assert_close(
        normalized.frequency_bins[0][&1] * normalized.row_scaling_factors[0],
        1.0,
    );
    Ok(())
}

#[test]
fn normalize_count_bins_rejects_negative_counts() {
    // Arrange: negative counts would make both frequencies and count reconstruction undefined for
    // this output contract.
    let count_bins = vec![FxHashMap::from_iter([(0, -1.0)])];

    // Act
    let error =
        normalize_count_bins_to_frequencies(count_bins).expect_err("negative count should fail");

    // Assert
    assert!(
        error.to_string().contains("negative"),
        "unexpected error: {error:#}"
    );
}

#[test]
fn normalize_count_bins_keeps_tiny_positive_reconstructable_frequencies() -> Result<()> {
    // Arrange: the rare motif has count 1 in a row with total count 1e18 + 1. Its frequency is far
    // below the sparse count tolerance, but dropping it would make count reconstruction incomplete.
    let count_bins = vec![FxHashMap::from_iter([(0, 1.0), (1, 1.0e18)])];

    // Act
    let normalized = normalize_count_bins_to_frequencies(count_bins)?;

    // Assert
    let scaling_factor = normalized.row_scaling_factors[0];
    let rare_frequency = normalized.frequency_bins[0][&0];
    assert!(rare_frequency > 0.0);
    assert_close(rare_frequency * scaling_factor, 1.0);
    Ok(())
}

#[test]
fn postprocess_ref_kmer_counts_decodes_and_canonicalizes_reverse_complements() -> Result<()> {
    // Arrange: for k = 2, GT reverse-complements to AC. Canonical output should collapse AC and GT
    // into one AC column with scaling factor 5 and frequency 1.
    let spec = kmer_spec(2);
    let code_ac = spec.encode_kmer_bytes(b"AC");
    let code_gt = spec.encode_kmer_bytes(b"GT");
    let mut counts_by_window = KmerCountsByWindow::default();
    counts_by_window.insert(
        0,
        KmerCounts {
            counts: FxHashMap::from_iter([
                (
                    Kmer {
                        k: 2,
                        code: code_ac,
                        orientation: KmerOrientation::Forward,
                    },
                    2.0,
                ),
                (
                    Kmer {
                        k: 2,
                        code: code_gt,
                        orientation: KmerOrientation::Forward,
                    },
                    3.0,
                ),
            ]),
        },
    );

    // Act
    let (frequencies, motif_order) =
        postprocess_ref_kmer_counts(counts_by_window, 1, &spec, true, false)?;

    // Assert
    assert_eq!(motif_order, vec!["AC"]);
    assert_eq!(frequencies.row_scaling_factors, vec![5.0]);
    assert_close(frequencies.frequency_bins[0][&0], 1.0);
    Ok(())
}

#[test]
fn postprocess_ref_kmer_counts_rejects_out_of_bounds_row() {
    // Arrange: total_windows = 1 allows only row 0.
    let spec = kmer_spec(2);
    let mut counts_by_window = KmerCountsByWindow::default();
    counts_by_window.insert(
        1,
        KmerCounts {
            counts: FxHashMap::from_iter([(
                Kmer {
                    k: 2,
                    code: spec.encode_kmer_bytes(b"AC"),
                    orientation: KmerOrientation::Forward,
                },
                1.0,
            )]),
        },
    );

    // Act
    let error = postprocess_ref_kmer_counts(counts_by_window, 1, &spec, false, false)
        .expect_err("out-of-bounds row should fail");

    // Assert
    assert!(
        error.to_string().contains("out of bounds"),
        "unexpected error: {error:#}"
    );
}

#[test]
fn build_all_ref_kmer_order_collapses_odd_canonical_universe() -> Result<()> {
    // Arrange: for k = 1, A and T collapse to A, while C and G collapse to C.
    let spec = kmer_spec(1);

    // Act
    let motif_order = build_all_ref_kmer_order(&spec, true)?;

    // Assert
    assert_eq!(motif_order, vec!["A", "C"]);
    Ok(())
}

#[test]
fn complete_ref_kmer_axis_len_counts_noncanonical_universe() {
    // Arrange: without canonicalization, every A/C/G/T string gets its own label.

    // Act
    let k1_axis_len = complete_ref_kmer_axis_len(1, false);
    let k2_axis_len = complete_ref_kmer_axis_len(2, false);
    let k3_axis_len = complete_ref_kmer_axis_len(3, false);

    // Assert
    assert_eq!(k1_axis_len, Some(4));
    assert_eq!(k2_axis_len, Some(16));
    assert_eq!(k3_axis_len, Some(64));
}

#[test]
fn complete_ref_kmer_axis_len_counts_odd_canonical_universe() {
    // Arrange: odd-length k-mers have no self reverse-complements, so canonicalization halves the
    // complete A/C/G/T universe.

    // Act
    let k1_axis_len = complete_ref_kmer_axis_len(1, true);
    let k3_axis_len = complete_ref_kmer_axis_len(3, true);

    // Assert
    assert_eq!(k1_axis_len, Some(2));
    assert_eq!(k3_axis_len, Some(32));
}

#[test]
fn complete_ref_kmer_axis_len_counts_even_canonical_universe() {
    // Arrange: even-length k-mers include self reverse-complements. For k = 2, the four fixed
    // points are AT, CG, GC, and TA. The remaining 12 k-mers form 6 reverse-complement pairs, so
    // k = 2 has 10 labels. For k = 4, there are 4^2 = 16 fixed points and 240 remaining k-mers,
    // giving 120 reverse-complement pairs plus 16 fixed points.

    // Act
    let k2_axis_len = complete_ref_kmer_axis_len(2, true);
    let k4_axis_len = complete_ref_kmer_axis_len(4, true);

    // Assert
    assert_eq!(k2_axis_len, Some(10));
    assert_eq!(k4_axis_len, Some(136));
}

#[test]
fn complete_ref_kmer_axis_len_returns_none_for_invalid_or_oversized_universes() {
    // Arrange: k = 0 is not a motif length, and 4^32 is one larger than u64::MAX.

    // Act
    let zero_axis_len = complete_ref_kmer_axis_len(0, false);
    let oversized_axis_len = complete_ref_kmer_axis_len(32, false);

    // Assert
    assert_eq!(zero_axis_len, None);
    assert_eq!(oversized_axis_len, None);
}

#[test]
fn ref_kmer_axis_is_complete_for_noncanonical_motif_axis() {
    // Arrange: for k = 1, the complete non-canonical A/C/G/T axis has exactly four labels.
    let motif_labels = ["A", "C", "G", "T"].map(str::to_string);

    // Act
    let is_complete = ref_kmer_axis_is_complete(1, false, SelectedMotifColumnKind::Motif, &motif_labels);

    // Assert
    assert!(is_complete);
}

#[test]
fn ref_kmer_axis_is_complete_for_canonical_motif_axis() {
    // Arrange: for k = 2, there are (4^2 + 4^1) / 2 = 10 reverse-complement classes.
    // The four self-reverse-complement motifs are AT, CG, GC, and TA.
    let motif_labels =
        ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "GA", "GC", "TA"].map(str::to_string);

    // Act
    let is_complete = ref_kmer_axis_is_complete(2, true, SelectedMotifColumnKind::Motif, &motif_labels);

    // Assert
    assert!(is_complete);
}

#[test]
fn ref_kmer_axis_is_not_complete_for_motif_group_axis() {
    // Arrange: group labels are not concrete k-mer labels, even if their count matches 4^k.
    let motif_labels = ["A", "C", "G", "T"].map(str::to_string);

    // Act
    let is_complete =
        ref_kmer_axis_is_complete(1, false, SelectedMotifColumnKind::MotifGroup, &motif_labels);

    // Assert
    assert!(!is_complete);
}

#[test]
fn ref_kmer_axis_is_not_complete_when_motif_count_is_short() {
    // Arrange: the complete non-canonical k = 1 axis has four A/C/G/T labels.
    let motif_labels = ["A", "C", "G"].map(str::to_string);

    // Act
    let is_complete = ref_kmer_axis_is_complete(1, false, SelectedMotifColumnKind::Motif, &motif_labels);

    // Assert
    assert!(!is_complete);
}

#[test]
fn validate_ref_kmer_package_rejects_invalid_concrete_motif_labels() {
    // Arrange: motif labels are concrete k-mers in public metadata, so an N base is invalid even
    // when sparse output would otherwise be able to write a string label.
    let motif_labels = vec!["AN".to_string()];

    // Act
    let error =
        validate_sparse_global_package(&motif_labels, SelectedMotifColumnKind::Motif, 2, false)
            .expect_err("invalid concrete motif label should fail");

    // Assert
    assert!(
        error.to_string().contains("contains invalid base `N`"),
        "unexpected error: {error:#}"
    );
}

#[test]
fn validate_ref_kmer_package_rejects_wrong_length_concrete_motif_labels() {
    // Arrange: public concrete motif labels must match the configured k-mer size exactly.
    let motif_labels = vec!["ACG".to_string()];

    // Act
    let error =
        validate_sparse_global_package(&motif_labels, SelectedMotifColumnKind::Motif, 2, false)
            .expect_err("wrong-length concrete motif label should fail");

    // Assert
    assert!(
        error.to_string().contains("has length 3, expected 2"),
        "unexpected error: {error:#}"
    );
}

#[test]
fn validate_ref_kmer_package_rejects_duplicate_concrete_motif_labels() {
    // Arrange: repeated concrete labels would make one motif column ambiguous.
    let motif_labels = vec!["A".to_string(), "A".to_string()];

    // Act
    let error =
        validate_sparse_global_package(&motif_labels, SelectedMotifColumnKind::Motif, 1, false)
            .expect_err("duplicate concrete motif label should fail");

    // Assert
    assert!(
        error
            .to_string()
            .contains("duplicate reference k-mer motif label `A`"),
        "unexpected error: {error:#}"
    );
}

#[test]
fn validate_ref_kmer_package_rejects_noncanonical_concrete_motif_labels() {
    // Arrange: canonical k = 2 output should write AA for the AA/TT reverse-complement class.
    let motif_labels = vec!["TT".to_string()];

    // Act
    let error =
        validate_sparse_global_package(&motif_labels, SelectedMotifColumnKind::Motif, 2, true)
            .expect_err("non-canonical motif label should fail");

    // Assert
    assert!(
        error.to_string().contains("should be represented as `AA`"),
        "unexpected error: {error:#}"
    );
}

#[test]
fn validate_ref_kmer_package_rejects_duplicate_motif_group_labels() {
    // Arrange: output group labels must still name distinct motif columns.
    let motif_labels = vec!["group".to_string(), "group".to_string()];

    // Act
    let error =
        validate_sparse_global_package(&motif_labels, SelectedMotifColumnKind::MotifGroup, 2, false)
            .expect_err("duplicate motif-group label should fail");

    // Assert
    assert!(
        error
            .to_string()
            .contains("duplicate reference k-mer motif-group label `group`"),
        "unexpected error: {error:#}"
    );
}

#[test]
fn dense_ref_kmer_zarr_writes_frequencies_scaling_and_window_metadata() -> Result<()> {
    // Arrange: row 0 frequencies reconstruct counts 2 and 1 with scaling factor 3. Row 1 has one
    // motif with scaling factor 4. BED row metadata is already in output-row order.
    let temp = TempDir::new().expect("temp dir should be created");
    let store_path = temp.path().join("sample.ref_kmer_counts.zarr");
    let motif_labels = vec!["AC".to_string(), "GT".to_string()];
    let frequency_bins = vec![
        FxHashMap::from_iter([(0, 2.0 / 3.0), (1, 1.0 / 3.0)]),
        FxHashMap::from_iter([(1, 1.0)]),
    ];
    let scaling_factors = vec![3.0, 4.0];
    let bin_info = vec![
        WindowBinInfo {
            chromosome: "chr2".to_string(),
            start: 10,
            end: 20,
            output_index: 0,
            blacklisted_fraction: 0.25,
        },
        WindowBinInfo {
            chromosome: "chr10".to_string(),
            start: 40,
            end: 60,
            output_index: 1,
            blacklisted_fraction: 0.0,
        },
    ];
    let reference_contig_footprint = footprint();

    // Act
    write_ref_kmer_zarr(
        &store_path,
        RefKmerZarrPackage {
            frequency_bins: &frequency_bins,
            row_scaling_factors: &scaling_factors,
            motif_labels: &motif_labels,
            column_kind: SelectedMotifColumnKind::Motif,
            row_metadata: RefKmerRowMetadata::Windows {
                bin_info: &bin_info,
                row_mode: RefKmerWindowRowMode::Bed,
            },
            write_dense_output: true,
            kmer_size: 2,
            canonical: false,
            all_motifs: true,
            assign_by: WindowAssigner::CountOverlap,
            reference_contig_footprint: &reference_contig_footprint,
        },
    )?;

    // Assert
    let root_metadata = read_json(&store_path.join("zarr.json"));
    assert_eq!(
        root_metadata["attributes"]["cfdnalab_schema"],
        "ref_kmer_frequencies"
    );
    assert_eq!(
        root_metadata["attributes"]["cfdnalab_schema_version"],
        serde_json::json!(1)
    );
    assert_eq!(root_metadata["attributes"]["storage_mode"], "dense");
    assert_eq!(root_metadata["attributes"]["row_mode"], "bed");
    assert_eq!(root_metadata["attributes"]["primary_array"], "frequencies");
    assert_eq!(
        root_metadata["attributes"]["count_reconstruction"],
        "reference_kmer_count = frequency * row_scaling_factor[row]"
    );
    assert_eq!(
        read_f64_array(&store_path, "/frequencies"),
        vec![2.0 / 3.0, 1.0 / 3.0, 0.0, 1.0]
    );
    assert_eq!(
        read_f64_array(&store_path, "/row_scaling_factor"),
        vec![3.0, 4.0]
    );
    assert_eq!(read_i32_array(&store_path, "/motif_index"), vec![0, 1]);
    assert_eq!(read_i32_array(&store_path, "/motif_byte"), vec![0, 1]);
    assert_eq!(read_u8_array(&store_path, "/motif_ascii"), b"ACGT".to_vec());
    assert_eq!(read_i32_array(&store_path, "/row_chromosome"), vec![0, 1]);
    assert_eq!(read_i64_array(&store_path, "/row_start_bp"), vec![10, 40]);
    assert_eq!(read_i64_array(&store_path, "/row_end_bp"), vec![20, 60]);
    assert_eq!(
        read_f64_array(&store_path, "/blacklisted_fraction"),
        vec![0.25, 0.0]
    );
    let chromosome_metadata = read_json(&store_path.join("chromosome/zarr.json"));
    assert_eq!(
        chromosome_metadata["attributes"]["labels"],
        serde_json::json!(["chr2", "chr10"])
    );
    let footprint_json = read_u8_array(&store_path, "/reference_contig_footprint_json");
    assert_eq!(
        serde_json::from_slice::<Vec<ContigFootprintEntry>>(&footprint_json)
            .expect("footprint JSON should parse"),
        reference_contig_footprint
    );
    Ok(())
}

#[test]
fn sparse_ref_kmer_zarr_writes_sorted_frequency_coo_arrays() -> Result<()> {
    // Arrange: sparse entries are inserted out of column order, but COO arrays should be sorted by
    // row and motif. The stored values are frequencies; row_scaling_factor reconstructs counts.
    let temp = TempDir::new().expect("temp dir should be created");
    let store_path = temp.path().join("sample.ref_kmer_counts.zarr");
    let motif_labels = vec!["AA".to_string(), "CC".to_string()];
    let frequency_bins = vec![FxHashMap::from_iter([(1, 0.25), (0, 0.75)])];
    let scaling_factors = vec![4.0];
    let reference_contig_footprint = footprint();

    // Act
    write_ref_kmer_zarr(
        &store_path,
        RefKmerZarrPackage {
            frequency_bins: &frequency_bins,
            row_scaling_factors: &scaling_factors,
            motif_labels: &motif_labels,
            column_kind: SelectedMotifColumnKind::Motif,
            row_metadata: RefKmerRowMetadata::Global,
            write_dense_output: false,
            kmer_size: 2,
            canonical: false,
            all_motifs: false,
            assign_by: WindowAssigner::Any,
            reference_contig_footprint: &reference_contig_footprint,
        },
    )?;

    // Assert
    let root_metadata = read_json(&store_path.join("zarr.json"));
    assert_eq!(root_metadata["attributes"]["storage_mode"], "sparse_coo");
    assert!(root_metadata["attributes"]["primary_array"].is_null());
    assert_eq!(root_metadata["attributes"]["primary_group"], "sparse");
    assert_eq!(read_i32_array(&store_path, "/sparse/row"), vec![0, 0]);
    assert_eq!(read_i32_array(&store_path, "/sparse/motif"), vec![0, 1]);
    assert_eq!(
        read_f64_array(&store_path, "/sparse/frequency"),
        vec![0.75, 0.25]
    );
    assert_eq!(read_i32_array(&store_path, "/sparse/shape"), vec![1, 2]);
    assert_eq!(
        read_f64_array(&store_path, "/row_scaling_factor"),
        vec![4.0]
    );
    assert!(!store_path.join("frequencies").exists());
    Ok(())
}

#[test]
fn sparse_ref_kmer_zarr_writes_row_major_coo_and_omits_zero_frequencies() -> Result<()> {
    // Arrange: sparse entries should be row-major by row and motif, regardless of hash-map
    // insertion order. Exact zero frequencies are omitted, but empty rows still contribute to
    // sparse/shape and row_scaling_factor.
    let temp = TempDir::new().expect("temp dir should be created");
    let store_path = temp.path().join("sample.ref_kmer_counts.zarr");
    let motif_labels = vec!["A".to_string(), "C".to_string(), "G".to_string()];
    let frequency_bins: Vec<RefKmerFrequencyBin> = vec![
        FxHashMap::from_iter([(2, 0.5), (1, 0.0), (0, 0.5)]),
        FxHashMap::default(),
        FxHashMap::from_iter([(2, 0.75), (1, 0.25)]),
    ];
    let scaling_factors = vec![2.0, 0.0, 4.0];
    let groups = vec![
        RefKmerGroupSummary {
            group_idx: 0,
            group_name: "row-zero",
            eligible_windows: 1,
            blacklisted_fraction: 0.0,
        },
        RefKmerGroupSummary {
            group_idx: 1,
            group_name: "row-one",
            eligible_windows: 0,
            blacklisted_fraction: 0.0,
        },
        RefKmerGroupSummary {
            group_idx: 2,
            group_name: "row-two",
            eligible_windows: 1,
            blacklisted_fraction: 0.0,
        },
    ];
    let reference_contig_footprint = footprint();

    // Act
    write_ref_kmer_zarr(
        &store_path,
        RefKmerZarrPackage {
            frequency_bins: &frequency_bins,
            row_scaling_factors: &scaling_factors,
            motif_labels: &motif_labels,
            column_kind: SelectedMotifColumnKind::Motif,
            row_metadata: RefKmerRowMetadata::Groups(groups),
            write_dense_output: false,
            kmer_size: 1,
            canonical: false,
            all_motifs: false,
            assign_by: WindowAssigner::Any,
            reference_contig_footprint: &reference_contig_footprint,
        },
    )?;

    // Assert
    assert_eq!(read_i32_array(&store_path, "/sparse/row"), vec![0, 0, 2, 2]);
    assert_eq!(read_i32_array(&store_path, "/sparse/motif"), vec![0, 2, 1, 2]);
    assert_eq!(
        read_f64_array(&store_path, "/sparse/frequency"),
        vec![0.5, 0.5, 0.25, 0.75]
    );
    assert_eq!(read_i32_array(&store_path, "/sparse/shape"), vec![3, 3]);
    assert_eq!(
        read_f64_array(&store_path, "/row_scaling_factor"),
        scaling_factors
    );
    Ok(())
}

#[test]
fn sparse_ref_kmer_zarr_writes_empty_coo_arrays_when_all_rows_are_empty() -> Result<()> {
    // Arrange: a sparse package may have no non-zero frequency entries. It should still write the
    // empty COO arrays plus shape metadata so loaders can distinguish an empty matrix from a
    // malformed package.
    let temp = TempDir::new().expect("temp dir should be created");
    let store_path = temp.path().join("sample.ref_kmer_counts.zarr");
    let motif_labels = vec!["A".to_string(), "C".to_string()];
    let frequency_bins: Vec<RefKmerFrequencyBin> = vec![FxHashMap::default()];
    let scaling_factors = vec![0.0];
    let reference_contig_footprint = footprint();

    // Act
    write_ref_kmer_zarr(
        &store_path,
        RefKmerZarrPackage {
            frequency_bins: &frequency_bins,
            row_scaling_factors: &scaling_factors,
            motif_labels: &motif_labels,
            column_kind: SelectedMotifColumnKind::Motif,
            row_metadata: RefKmerRowMetadata::Global,
            write_dense_output: false,
            kmer_size: 1,
            canonical: false,
            all_motifs: false,
            assign_by: WindowAssigner::Any,
            reference_contig_footprint: &reference_contig_footprint,
        },
    )?;

    // Assert
    assert_eq!(read_i32_array(&store_path, "/sparse/row"), Vec::<i32>::new());
    assert_eq!(
        read_i32_array(&store_path, "/sparse/motif"),
        Vec::<i32>::new()
    );
    assert_eq!(
        read_f64_array(&store_path, "/sparse/frequency"),
        Vec::<f64>::new()
    );
    assert_eq!(read_i32_array(&store_path, "/sparse/shape"), vec![1, 2]);
    assert_eq!(
        read_i32_array(&store_path, "/sparse/sparse_dimension"),
        vec![0, 1]
    );
    assert_eq!(
        read_f64_array(&store_path, "/row_scaling_factor"),
        vec![0.0]
    );
    Ok(())
}

#[test]
fn window_ref_kmer_zarr_rejects_out_of_order_output_indices() {
    // Arrange: window metadata rows are expected to already be sorted in output-row order. The writer
    // checks monotonic output_index values so swapped rows fail before a package can look valid.
    let temp = TempDir::new().expect("temp dir should be created");
    let store_path = temp.path().join("sample.ref_kmer_counts.zarr");
    let motif_labels = vec!["A".to_string()];
    let frequency_bins = vec![
        FxHashMap::from_iter([(0, 1.0)]),
        FxHashMap::from_iter([(0, 1.0)]),
    ];
    let scaling_factors = vec![1.0, 1.0];
    let bin_info = vec![
        WindowBinInfo {
            chromosome: "chr1".to_string(),
            start: 10,
            end: 20,
            output_index: 1,
            blacklisted_fraction: 0.0,
        },
        WindowBinInfo {
            chromosome: "chr1".to_string(),
            start: 0,
            end: 10,
            output_index: 0,
            blacklisted_fraction: 0.0,
        },
    ];
    let reference_contig_footprint = footprint();

    // Act
    let error = write_ref_kmer_zarr(
        &store_path,
        RefKmerZarrPackage {
            frequency_bins: &frequency_bins,
            row_scaling_factors: &scaling_factors,
            motif_labels: &motif_labels,
            column_kind: SelectedMotifColumnKind::Motif,
            row_metadata: RefKmerRowMetadata::Windows {
                bin_info: &bin_info,
                row_mode: RefKmerWindowRowMode::Bed,
            },
            write_dense_output: true,
            kmer_size: 1,
            canonical: false,
            all_motifs: false,
            assign_by: WindowAssigner::Any,
            reference_contig_footprint: &reference_contig_footprint,
        },
    )
    .expect_err("out-of-order output_index values should fail");

    // Assert
    assert!(
        error
            .to_string()
            .contains("sorted by increasing output_index"),
        "unexpected error: {error:#}"
    );
}

#[test]
fn motif_group_ref_kmer_zarr_writes_json_labels_without_motif_ascii() -> Result<()> {
    // Arrange: motifs-file group labels can have variable width, so they live as JSON labels on the
    // motif_index axis rather than in motif_ascii.
    let temp = TempDir::new().expect("temp dir should be created");
    let store_path = temp.path().join("sample.ref_kmer_counts.zarr");
    let motif_labels = vec!["short".to_string(), "group-two".to_string()];
    let frequency_bins = vec![FxHashMap::from_iter([(0, 0.25), (1, 0.75)])];
    let scaling_factors = vec![8.0];
    let reference_contig_footprint = footprint();

    // Act
    write_ref_kmer_zarr(
        &store_path,
        RefKmerZarrPackage {
            frequency_bins: &frequency_bins,
            row_scaling_factors: &scaling_factors,
            motif_labels: &motif_labels,
            column_kind: SelectedMotifColumnKind::MotifGroup,
            row_metadata: RefKmerRowMetadata::Global,
            write_dense_output: false,
            kmer_size: 4,
            canonical: false,
            all_motifs: false,
            assign_by: WindowAssigner::All,
            reference_contig_footprint: &reference_contig_footprint,
        },
    )?;

    // Assert
    let root_metadata = read_json(&store_path.join("zarr.json"));
    assert_eq!(
        root_metadata["attributes"]["motif_axis_kind"],
        serde_json::json!("motif_group")
    );
    let motif_metadata = read_json(&store_path.join("motif_index/zarr.json"));
    assert_eq!(motif_metadata["attributes"]["label_field"], "motif_group");
    assert_eq!(
        motif_metadata["attributes"]["labels"],
        serde_json::json!(["short", "group-two"])
    );
    assert!(!store_path.join("motif_byte").exists());
    assert!(!store_path.join("motif_ascii").exists());
    Ok(())
}

#[test]
fn grouped_ref_kmer_zarr_writes_group_metadata_and_dense_frequencies() -> Result<()> {
    // Arrange: group rows are already count-row indices. Group 1 has no nonzero frequencies but still
    // keeps its metadata and scaling factor 0.
    let temp = TempDir::new().expect("temp dir should be created");
    let store_path = temp.path().join("sample.ref_kmer_counts.zarr");
    let motif_labels = vec!["A".to_string(), "C".to_string()];
    let frequency_bins = vec![FxHashMap::from_iter([(0, 1.0)]), FxHashMap::default()];
    let scaling_factors = vec![5.0, 0.0];
    let groups = vec![
        RefKmerGroupSummary {
            group_idx: 0,
            group_name: "promoter",
            eligible_windows: 2,
            blacklisted_fraction: 0.125,
        },
        RefKmerGroupSummary {
            group_idx: 1,
            group_name: "enhancer",
            eligible_windows: 1,
            blacklisted_fraction: 0.0,
        },
    ];
    let reference_contig_footprint = footprint();

    // Act
    write_ref_kmer_zarr(
        &store_path,
        RefKmerZarrPackage {
            frequency_bins: &frequency_bins,
            row_scaling_factors: &scaling_factors,
            motif_labels: &motif_labels,
            column_kind: SelectedMotifColumnKind::Motif,
            row_metadata: RefKmerRowMetadata::Groups(groups),
            write_dense_output: true,
            kmer_size: 1,
            canonical: false,
            all_motifs: true,
            assign_by: WindowAssigner::Proportion(0.5),
            reference_contig_footprint: &reference_contig_footprint,
        },
    )?;

    // Assert
    let root_metadata = read_json(&store_path.join("zarr.json"));
    assert_eq!(root_metadata["attributes"]["row_mode"], "grouped_bed");
    assert_eq!(root_metadata["attributes"]["assign_by"], "proportion=0.5");
    assert_eq!(read_i32_array(&store_path, "/group"), vec![0, 1]);
    assert_eq!(read_i32_array(&store_path, "/eligible_windows"), vec![2, 1]);
    assert_eq!(
        read_f64_array(&store_path, "/blacklisted_fraction"),
        vec![0.125, 0.0]
    );
    assert_eq!(
        read_f64_array(&store_path, "/frequencies"),
        vec![1.0, 0.0, 0.0, 0.0]
    );
    assert_eq!(
        read_f64_array(&store_path, "/row_scaling_factor"),
        vec![5.0, 0.0]
    );
    let group_metadata = read_json(&store_path.join("group/zarr.json"));
    assert_eq!(
        group_metadata["attributes"]["labels"],
        serde_json::json!(["promoter", "enhancer"])
    );
    Ok(())
}

#[test]
fn grouped_ref_kmer_row_metadata_keeps_count_row_order_and_blacklist_fractions() -> Result<()> {
    // Arrange: group 1 has two windows totaling 30 bp. A 5 bp blacklist overlap gives 5/30.
    let mut group_idx_to_name = FxHashMap::default();
    group_idx_to_name.insert(0, "zero".to_string());
    group_idx_to_name.insert(1, "one".to_string());
    let mut grouped_windows_map = FxHashMap::default();
    grouped_windows_map.insert(
        "chr1".to_string(),
        GroupedWindows::from_tuples(&[(0, 10, 0), (10, 30, 1), (40, 50, 1)], None)?,
    );
    let mut blacklist_map = FxHashMap::default();
    blacklist_map.insert("chr1".to_string(), vec![Interval::new(15, 20)?]);

    // Act
    let summaries = grouped_ref_kmer_row_metadata(
        &group_idx_to_name,
        &["chr1".to_string()],
        &grouped_windows_map,
        &blacklist_map,
    )?;

    // Assert
    assert_eq!(summaries.len(), 2);
    assert_eq!(summaries[0].group_idx, 0);
    assert_eq!(summaries[0].group_name, "zero");
    assert_eq!(summaries[0].eligible_windows, 1);
    assert_close(summaries[0].blacklisted_fraction, 0.0);
    assert_eq!(summaries[1].group_idx, 1);
    assert_eq!(summaries[1].group_name, "one");
    assert_eq!(summaries[1].eligible_windows, 2);
    assert_close(summaries[1].blacklisted_fraction, 5.0 / 30.0);
    Ok(())
}

fn read_json(path: &std::path::Path) -> Value {
    serde_json::from_str(&std::fs::read_to_string(path).expect("JSON file should read"))
        .expect("JSON should parse")
}

fn read_f64_array(store_path: &std::path::Path, array_path: &str) -> Vec<f64> {
    read_array(store_path, array_path)
}

fn read_i32_array(store_path: &std::path::Path, array_path: &str) -> Vec<i32> {
    read_array(store_path, array_path)
}

fn read_i64_array(store_path: &std::path::Path, array_path: &str) -> Vec<i64> {
    read_array(store_path, array_path)
}

fn read_u8_array(store_path: &std::path::Path, array_path: &str) -> Vec<u8> {
    read_array(store_path, array_path)
}

fn read_array<T>(store_path: &std::path::Path, array_path: &str) -> Vec<T>
where
    T: zarrs::array::ElementOwned,
{
    let store = Arc::new(FilesystemStore::new(store_path).expect("Zarr store should open"));
    let array = Array::open(store, array_path).expect("Zarr array should open");
    array
        .retrieve_array_subset(&array.subset_all())
        .expect("Zarr array should read")
}
