#![cfg(all(feature = "cmd_ref_kmers", feature = "testing"))]

use anyhow::Result;
use cfdnalab::{
    RunOptions,
    output_loaders::{
        RefKmerFrequencyData, RefKmerMotifAxisKind, RefKmerRowMetadata, RefKmerRowMode,
        RefKmerStorageMode, RefKmerWindowMode, load_ref_kmers_output,
    },
    reference::twobit_contig_footprint,
    run_like_cli::{
        common::{ChromosomeArgs, DistributionWindowsArgs, WindowAssigner},
        ref_kmers::{RefKmersConfig, run_ref_kmers},
    },
    testing::{Bed4Row, twobit_from_sequences, write_bed4},
};
use serde_json::Value;
use std::{path::Path, sync::Arc};
use tempfile::{NamedTempFile, TempDir};
use zarrs::{array::Array, filesystem::FilesystemStore};

fn ref_kmers_config(reference_path: &Path, output_dir: &Path, kmer_size: u8) -> RefKmersConfig {
    let mut config = RefKmersConfig::new(
        reference_path.to_path_buf(),
        output_dir.to_path_buf(),
        kmer_size,
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
    );
    config.set_n_threads(1);
    config
}

fn run(config: &RefKmersConfig) -> Result<()> {
    run_ref_kmers(config, RunOptions::new_quiet()).map(|_| ())
}

fn write_motifs_file(contents: &str) -> Result<NamedTempFile> {
    let file = NamedTempFile::new()?;
    std::fs::write(file.path(), contents)?;
    Ok(file)
}

fn assert_close(observed: f64, expected: f64) {
    assert!(
        (observed - expected).abs() < 1e-12,
        "observed {observed}, expected {expected}"
    );
}

fn assert_slice_close(observed: &[f64], expected: &[f64]) {
    assert_eq!(observed.len(), expected.len());
    for (observed_value, expected_value) in observed.iter().zip(expected.iter()) {
        assert_close(*observed_value, *expected_value);
    }
}

#[test]
fn ref_kmers_by_size_writes_fractional_frequencies_and_scaling_factors() -> Result<()> {
    // Arrange:
    // Reference chr1 is ACGTACGT and k = 4. Fixed windows are [0,4) and [4,8).
    //
    // The five 4-mers are:
    //   ACGT [0,4), CGTA [1,5), GTAC [2,6), TACG [3,7), ACGT [4,8).
    //
    // Under count-overlap:
    //   row 0 counts are ACGT=1, CGTA=3/4, GTAC=1/2, TACG=1/4, total=2.5.
    //   row 1 counts are ACGT=1, CGTA=1/4, GTAC=1/2, TACG=3/4, total=2.5.
    let reference = twobit_from_sequences(
        "ref_kmers_by_size",
        vec![("chr1".to_string(), "ACGTACGT".to_string())],
    )?;
    let output_dir = TempDir::new()?;
    let mut config = ref_kmers_config(&reference.path, output_dir.path(), 4);
    config.set_output_prefix("unit_ref_kmers");
    config.set_windows(DistributionWindowsArgs {
        by_size: Some(4),
        by_bed: None,
        by_grouped_bed: None,
    });
    config.set_assign_by(WindowAssigner::CountOverlap);

    // Act
    run(&config)?;

    // Assert
    let package_path = output_dir
        .path()
        .join("unit_ref_kmers.ref_kmer_counts.zarr");
    let metadata = read_json(&package_path.join("zarr.json"));
    assert_eq!(
        metadata["attributes"]["cfdnalab_schema"],
        "ref_kmer_frequencies"
    );
    assert_eq!(metadata["attributes"]["storage_mode"], "sparse_coo");
    assert_eq!(metadata["attributes"]["row_mode"], "size");
    assert_eq!(metadata["attributes"]["assign_by"], "count-overlap");

    assert_eq!(
        read_u8_array(&package_path, "/motif_ascii"),
        b"ACGTCGTAGTACTACG".to_vec()
    );
    assert_eq!(
        read_f64_array(&package_path, "/row_scaling_factor"),
        vec![2.5, 2.5]
    );
    assert_eq!(read_i64_array(&package_path, "/row_start_bp"), vec![0, 4]);
    assert_eq!(read_i64_array(&package_path, "/row_end_bp"), vec![4, 8]);

    let sparse_rows = read_i32_array(&package_path, "/sparse/row");
    let sparse_motifs = read_i32_array(&package_path, "/sparse/motif");
    let sparse_frequencies = read_f64_array(&package_path, "/sparse/frequency");
    assert_eq!(sparse_rows, vec![0, 0, 0, 0, 1, 1, 1, 1]);
    assert_eq!(sparse_motifs, vec![0, 1, 2, 3, 0, 1, 2, 3]);
    assert_close(sparse_frequencies[0], 1.0 / 2.5);
    assert_close(sparse_frequencies[1], 0.75 / 2.5);
    assert_close(sparse_frequencies[2], 0.50 / 2.5);
    assert_close(sparse_frequencies[3], 0.25 / 2.5);
    assert_close(sparse_frequencies[4], 1.0 / 2.5);
    assert_close(sparse_frequencies[5], 0.25 / 2.5);
    assert_close(sparse_frequencies[6], 0.50 / 2.5);
    assert_close(sparse_frequencies[7], 0.75 / 2.5);

    // Count reconstruction uses frequency * row_scaling_factor[row].
    let scaling = read_f64_array(&package_path, "/row_scaling_factor");
    assert_close(
        sparse_frequencies[1] * scaling[sparse_rows[1] as usize],
        0.75,
    );
    assert_close(
        sparse_frequencies[7] * scaling[sparse_rows[7] as usize],
        0.75,
    );

    let loaded = load_ref_kmers_output(&package_path)?;
    assert_eq!(loaded.storage_mode(), RefKmerStorageMode::SparseCoo);
    assert_eq!(loaded.row_mode(), RefKmerRowMode::SizeWindows);
    assert_eq!(loaded.motif_axis_kind(), RefKmerMotifAxisKind::Motif);
    assert_eq!(loaded.kmer_size(), 4);
    assert_eq!(loaded.assign_by(), "count-overlap");
    assert_eq!(loaded.row_scaling_factors(), &[2.5, 2.5]);
    let expected_footprint = twobit_contig_footprint(&reference.path)?;
    assert_eq!(
        loaded.reference_contig_footprint(),
        expected_footprint.as_slice()
    );
    assert_eq!(
        loaded.output_metadata().reference_contig_footprint,
        expected_footprint
    );
    loaded.ensure_reference_2bit_matches(&reference.path)?;
    let mismatched_reference = twobit_from_sequences(
        "ref_kmers_mismatched_reference",
        vec![("chr1".to_string(), "ACGTACGTA".to_string())],
    )?;
    let mismatch_error = loaded
        .ensure_reference_2bit_matches(&mismatched_reference.path)
        .expect_err("different reference footprint should fail");
    assert!(
        mismatch_error
            .to_string()
            .contains("different reference contig footprint"),
        "unexpected error: {mismatch_error:#}"
    );
    assert_eq!(
        loaded.motif_labels(),
        &[
            "ACGT".to_string(),
            "CGTA".to_string(),
            "GTAC".to_string(),
            "TACG".to_string()
        ]
    );
    let cgta_index = loaded.motif_index("CGTA")?;
    assert_close(loaded.frequency(0, cgta_index).unwrap(), 0.75 / 2.5);
    assert_close(loaded.count(0, cgta_index).unwrap(), 0.75);
    let windows = loaded.window_metadata()?;
    assert_eq!(windows[0].chrom, "chr1");
    assert_eq!(windows[0].interval.as_tuple(), (0, 4));
    assert_eq!(windows[1].interval.as_tuple(), (4, 8));

    Ok(())
}

#[test]
fn ref_kmers_motifs_file_groups_selected_targets_end_to_end() -> Result<()> {
    // Arrange:
    // For ACGTACGT and k = 2, the valid starts are:
    //   AC, CG, GT, TA, AC, CG, GT.
    //
    // The motifs file selects AC and GT into group `edge`, and CG into group `middle`.
    // TA is unselected. The selected counts are edge=4, middle=2, total=6.
    let reference = twobit_from_sequences(
        "ref_kmers_grouped_motifs",
        vec![("chr1".to_string(), "ACGTACGT".to_string())],
    )?;
    let output_dir = TempDir::new()?;
    let motifs_file = write_motifs_file("AC\tedge\nCG\tmiddle\nGT\tedge\n")?;
    let mut config = ref_kmers_config(&reference.path, output_dir.path(), 2);
    config.set_output_prefix("unit_grouped_ref_kmers");
    config.set_motifs_file(Some(motifs_file.path().to_path_buf()));
    config.set_assign_by(WindowAssigner::Any);

    // Act
    run(&config)?;

    // Assert
    let package_path = output_dir
        .path()
        .join("unit_grouped_ref_kmers.ref_kmer_counts.zarr");
    let metadata = read_json(&package_path.join("zarr.json"));
    assert_eq!(metadata["attributes"]["storage_mode"], "sparse_coo");
    assert_eq!(metadata["attributes"]["row_mode"], "global");
    assert_eq!(metadata["attributes"]["motif_axis_kind"], "motif_group");

    let motif_metadata = read_json(&package_path.join("motif_index/zarr.json"));
    assert_eq!(
        motif_metadata["attributes"]["labels"],
        serde_json::json!(["edge", "middle"])
    );
    assert!(!package_path.join("motif_ascii").exists());

    assert_eq!(
        read_f64_array(&package_path, "/row_scaling_factor"),
        vec![6.0]
    );
    assert_eq!(read_i32_array(&package_path, "/sparse/row"), vec![0, 0]);
    assert_eq!(read_i32_array(&package_path, "/sparse/motif"), vec![0, 1]);
    let sparse_frequencies = read_f64_array(&package_path, "/sparse/frequency");
    assert_close(sparse_frequencies[0], 4.0 / 6.0);
    assert_close(sparse_frequencies[1], 2.0 / 6.0);

    let loaded = load_ref_kmers_output(&package_path)?;
    assert_eq!(loaded.storage_mode(), RefKmerStorageMode::SparseCoo);
    assert_eq!(loaded.row_mode(), RefKmerRowMode::Global);
    assert_eq!(loaded.motif_axis_kind(), RefKmerMotifAxisKind::MotifGroup);
    assert_eq!(
        loaded.motif_labels(),
        &["edge".to_string(), "middle".to_string()]
    );
    assert_eq!(loaded.row_scaling_factors(), &[6.0]);
    assert_close(loaded.frequency_for_motif(0, "edge")?.unwrap(), 4.0 / 6.0);
    assert_close(loaded.count_for_motif(0, "edge")?.unwrap(), 4.0);
    assert_close(
        loaded
            .count_for_motif(0, "middle")?
            .expect("middle should be in bounds"),
        2.0,
    );

    Ok(())
}

#[test]
fn ref_kmers_loader_exposes_sparse_windows_and_implicit_zero_cells() -> Result<()> {
    // Arrange:
    // Reference chr1 is AACC and k = 1. Fixed windows are [0,2) and [2,4).
    //
    // The observed k-mers are:
    //   row 0: A=2, C=0, total=2
    //   row 1: A=0, C=2, total=2
    //
    // The motif axis contains A and C because both are observed somewhere. Sparse storage only
    // writes the two nonzero cells, so the loader must expose missing in-bounds cells as zero.
    let reference = twobit_from_sequences(
        "ref_kmers_sparse_loader_api",
        vec![("chr1".to_string(), "AACC".to_string())],
    )?;
    let output_dir = TempDir::new()?;
    let mut config = ref_kmers_config(&reference.path, output_dir.path(), 1);
    config.set_output_prefix("unit_sparse_loader_api");
    config.set_windows(DistributionWindowsArgs {
        by_size: Some(2),
        by_bed: None,
        by_grouped_bed: None,
    });

    // Act
    run(&config)?;

    // Assert
    let package_path = output_dir
        .path()
        .join("unit_sparse_loader_api.ref_kmer_counts.zarr");
    let loaded = load_ref_kmers_output(&package_path)?;
    let metadata = loaded.output_metadata();
    assert_eq!(metadata.storage_mode, RefKmerStorageMode::SparseCoo);
    assert_eq!(metadata.row_mode, RefKmerRowMode::SizeWindows);
    assert_eq!(metadata.motif_axis_kind, RefKmerMotifAxisKind::Motif);
    assert_eq!(metadata.row_count, 2);
    assert_eq!(metadata.motif_count, 2);
    assert_eq!(metadata.kmer_size, 1);
    assert!(!metadata.canonical);
    assert!(!metadata.all_motifs);
    assert_eq!(metadata.assign_by, "count-overlap");

    assert_eq!(loaded.motif_labels(), &["A".to_string(), "C".to_string()]);
    assert!(loaded.has_motif("A"));
    assert!(!loaded.has_motif("G"));
    assert_eq!(loaded.row_scaling_factors(), &[2.0, 2.0]);
    assert_eq!(loaded.row_scaling_factor(2), None);

    match loaded.row_metadata() {
        RefKmerRowMetadata::Windows {
            window_mode,
            windows,
        } => {
            assert_eq!(*window_mode, RefKmerWindowMode::Size);
            assert_eq!(windows.len(), 2);
            assert_eq!(windows[0].interval.as_tuple(), (0, 2));
            assert_eq!(windows[1].interval.as_tuple(), (2, 4));
        }
        other => panic!("expected size-window metadata, got {other:?}"),
    }

    let sparse = loaded.sparse_frequencies()?;
    assert_eq!(sparse.shape(), (2, 2));
    assert_eq!(sparse.nnz(), 2);
    assert_eq!(sparse.row_indices(), &[0, 1]);
    assert_eq!(sparse.motif_indices(), &[0, 1]);
    assert_eq!(sparse.frequencies(), &[1.0, 1.0]);
    let sparse_entries: Vec<_> = sparse
        .entries()
        .map(|entry| (entry.row_index, entry.motif_index, entry.frequency))
        .collect();
    assert_eq!(sparse_entries, vec![(0, 0, 1.0), (1, 1, 1.0)]);
    let sparse_count_entries: Vec<_> = loaded
        .sparse_count_entries()?
        .into_iter()
        .map(|entry| (entry.row_index, entry.motif_index, entry.count))
        .collect();
    assert_eq!(sparse_count_entries, vec![(0, 0, 2.0), (1, 1, 2.0)]);

    assert!(loaded.dense_frequencies().is_err());
    match loaded.data() {
        RefKmerFrequencyData::Sparse(data) => assert_eq!(data.nnz(), 2),
        other => panic!("expected sparse frequency data, got {other:?}"),
    }

    let a_index = loaded.motif_index("A")?;
    let c_index = loaded.motif_index("C")?;
    assert_close(loaded.frequency(0, a_index).unwrap(), 1.0);
    assert_close(loaded.frequency(0, c_index).unwrap(), 0.0);
    assert_close(loaded.frequency(1, a_index).unwrap(), 0.0);
    assert_close(loaded.frequency(1, c_index).unwrap(), 1.0);
    assert_close(loaded.count(0, a_index).unwrap(), 2.0);
    assert_close(loaded.count(0, c_index).unwrap(), 0.0);
    assert_eq!(loaded.frequency(2, a_index), None);
    assert_eq!(loaded.count(0, 2), None);
    assert!(loaded.frequency_for_motif(0, "G").is_err());

    assert_eq!(
        loaded.to_dense_frequency_matrix()?.values_row_major(),
        &[1.0, 0.0, 0.0, 1.0]
    );
    assert_eq!(
        loaded.to_dense_count_matrix()?.values_row_major(),
        &[2.0, 0.0, 0.0, 2.0]
    );

    let selected = loaded
        .select()
        .windows(&[1, 0])
        .motifs_by_label(&["C", "A"])
        .read()?;
    assert_eq!(selected.storage_mode(), RefKmerStorageMode::SparseCoo);
    assert_eq!(selected.row_mode(), RefKmerRowMode::SizeWindows);
    assert_eq!(selected.motif_axis_kind(), RefKmerMotifAxisKind::Motif);
    assert_eq!(selected.kmer_size(), 1);
    assert!(!selected.canonical());
    assert!(!selected.source_all_motifs());
    assert_eq!(selected.assign_by(), "count-overlap");
    assert_eq!(selected.row_indices(), &[1, 0]);
    assert_eq!(selected.motif_indices(), &[1, 0]);
    assert_eq!(selected.motif_labels(), &["C".to_string(), "A".to_string()]);
    assert_eq!(selected.row_scaling_factors(), &[2.0, 2.0]);
    assert_eq!(selected.shape(), (2, 2));
    assert_eq!(selected.row_count(), 2);
    assert_eq!(selected.motif_count(), 2);
    assert_close(selected.frequency(0, 0).unwrap(), 1.0);
    assert_close(selected.frequency(0, 1).unwrap(), 0.0);
    assert_close(selected.frequency(1, 0).unwrap(), 0.0);
    assert_close(selected.frequency(1, 1).unwrap(), 1.0);
    assert_close(selected.count(0, 0).unwrap(), 2.0);
    assert_close(selected.count(1, 1).unwrap(), 2.0);
    assert_eq!(
        selected
            .window_metadata()?
            .iter()
            .map(|window| window.interval.as_tuple())
            .collect::<Vec<_>>(),
        vec![(2, 4), (0, 2)]
    );
    let selected_sparse_entries: Vec<_> = selected
        .sparse_frequencies()?
        .entries()
        .map(|entry| (entry.row_index, entry.motif_index, entry.frequency))
        .collect();
    assert_eq!(selected_sparse_entries, vec![(0, 0, 1.0), (1, 1, 1.0)]);
    let selected_count_entries: Vec<_> = selected
        .sparse_count_entries()?
        .into_iter()
        .map(|entry| (entry.row_index, entry.motif_index, entry.count))
        .collect();
    assert_eq!(selected_count_entries, vec![(0, 0, 2.0), (1, 1, 2.0)]);
    assert_eq!(
        selected.to_dense_count_matrix()?.values_row_major(),
        &[2.0, 0.0, 0.0, 2.0]
    );

    let duplicate_window_error = loaded
        .select()
        .windows(&[0, 0])
        .read()
        .expect_err("duplicate window selector should fail");
    assert!(
        duplicate_window_error
            .to_string()
            .contains("duplicate value 0")
    );
    let conflicting_motif_error = loaded
        .select()
        .motifs(&[0])
        .motifs_by_label(&["A"])
        .read()
        .expect_err("conflicting motif selectors should fail");
    assert!(
        conflicting_motif_error
            .to_string()
            .contains("cannot combine motifs() and motifs_by_label()")
    );
    Ok(())
}

#[test]
fn ref_kmers_loader_exposes_grouped_bed_rows_by_name() -> Result<()> {
    // Arrange:
    // Reference chr1 is ACGTACGT and k = 2. Grouped BED rows are [0,4) alpha and [4,8) beta.
    //
    // Under `all`, each group contains AC, CG, and GT exactly once. TA spans the group boundary
    // and is not fully contained by either group.
    let reference = twobit_from_sequences(
        "ref_kmers_grouped_bed_loader_api",
        vec![("chr1".to_string(), "ACGTACGT".to_string())],
    )?;
    let output_dir = TempDir::new()?;
    let grouped_bed = output_dir.path().join("groups.bed");
    write_bed4(
        &grouped_bed,
        &[
            Bed4Row::new("chr1", 0, 4, "alpha"),
            Bed4Row::new("chr1", 4, 8, "beta"),
        ],
    )?;
    let mut config = ref_kmers_config(&reference.path, output_dir.path(), 2);
    config.set_output_prefix("unit_grouped_bed_loader_api");
    config.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: None,
        by_grouped_bed: Some(grouped_bed),
    });
    config.set_assign_by(WindowAssigner::All);

    // Act
    run(&config)?;

    // Assert
    let package_path = output_dir
        .path()
        .join("unit_grouped_bed_loader_api.ref_kmer_counts.zarr");
    let loaded = load_ref_kmers_output(&package_path)?;
    assert_eq!(loaded.row_mode(), RefKmerRowMode::Groups);
    assert_eq!(loaded.motif_axis_kind(), RefKmerMotifAxisKind::Motif);
    assert!(loaded.window_metadata().is_err());

    let groups = loaded.group_metadata()?;
    assert_eq!(groups.len(), 2);
    assert_eq!(groups[0].index, 0);
    assert_eq!(groups[0].name, "alpha");
    assert_eq!(groups[0].eligible_windows, 1);
    assert_close(groups[0].blacklisted_fraction, 0.0);
    assert_eq!(groups[1].index, 1);
    assert_eq!(groups[1].name, "beta");
    assert_eq!(groups[1].eligible_windows, 1);
    assert_close(groups[1].blacklisted_fraction, 0.0);
    assert_eq!(loaded.group_index("alpha")?, 0);
    assert_eq!(loaded.group_index("beta")?, 1);
    assert_eq!(loaded.group(1)?.expect("beta group exists").name, "beta");
    assert!(loaded.group_index("missing").is_err());
    assert!(loaded.has_group("alpha"));
    assert!(!loaded.has_group("missing"));

    assert_eq!(
        loaded.motif_labels(),
        &["AC".to_string(), "CG".to_string(), "GT".to_string()]
    );
    for group_index in [loaded.group_index("alpha")?, loaded.group_index("beta")?] {
        assert_eq!(loaded.row_scaling_factor(group_index), Some(3.0));
        for motif_label in ["AC", "CG", "GT"] {
            assert_close(
                loaded.count_for_motif(group_index, motif_label)?.unwrap(),
                1.0,
            );
            assert_close(
                loaded
                    .frequency_for_motif(group_index, motif_label)?
                    .unwrap(),
                1.0 / 3.0,
            );
        }
    }

    let selected = loaded
        .select()
        .groups_by_name(&["beta", "alpha"])
        .motifs_by_label(&["GT", "AC"])
        .read()?;
    assert_eq!(selected.storage_mode(), RefKmerStorageMode::SparseCoo);
    assert_eq!(selected.row_indices(), &[1, 0]);
    assert_eq!(selected.motif_indices(), &[2, 0]);
    assert_eq!(
        selected.motif_labels(),
        &["GT".to_string(), "AC".to_string()]
    );
    assert_eq!(selected.row_scaling_factors(), &[3.0, 3.0]);
    assert_eq!(
        selected
            .group_metadata()?
            .iter()
            .map(|group| group.name.as_str())
            .collect::<Vec<_>>(),
        vec!["beta", "alpha"]
    );
    let selected_count_entries: Vec<_> = selected
        .sparse_count_entries()?
        .into_iter()
        .map(|entry| (entry.row_index, entry.motif_index, entry.count))
        .collect();
    assert_eq!(
        selected_count_entries,
        vec![(0, 0, 1.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 1.0)]
    );
    assert_eq!(
        selected.to_dense_count_matrix()?.values_row_major(),
        &[1.0, 1.0, 1.0, 1.0]
    );
    let selected_by_group_index = loaded.select().groups(&[0]).motifs(&[0]).read()?;
    assert_eq!(selected_by_group_index.row_indices(), &[0]);
    assert_eq!(selected_by_group_index.group_metadata()?[0].name, "alpha");
    assert_close(selected_by_group_index.count(0, 0).unwrap(), 1.0);
    let duplicate_group_name_error = loaded
        .select()
        .groups_by_name(&["alpha", "alpha"])
        .read()
        .expect_err("duplicate group names should fail");
    assert!(
        duplicate_group_name_error
            .to_string()
            .contains("duplicate value 'alpha'")
    );
    Ok(())
}

#[test]
fn ref_kmers_grouped_bed_count_overlap_uses_manual_overlap_mass() -> Result<()> {
    // Arrange:
    // Reference chr1 is fourteen A bases and k = 4, so all eleven starts are AAAA.
    //
    // Group alpha has [2,6) and [8,12):
    //   3.75 + 3.75 = 7.50 overlap-weighted AAAA counts.
    // Group beta has [6,9):
    //   3.00 overlap-weighted AAAA counts.
    let reference = twobit_from_sequences(
        "ref_kmers_grouped_bed_count_overlap",
        vec![("chr1".to_string(), "A".repeat(14))],
    )?;
    let output_dir = TempDir::new()?;
    let grouped_bed = output_dir.path().join("count_overlap_groups.bed");
    write_bed4(
        &grouped_bed,
        &[
            Bed4Row::new("chr1", 2, 6, "alpha"),
            Bed4Row::new("chr1", 6, 9, "beta"),
            Bed4Row::new("chr1", 8, 12, "alpha"),
        ],
    )?;
    let mut config = ref_kmers_config(&reference.path, output_dir.path(), 4);
    config.set_output_prefix("unit_grouped_bed_count_overlap_ref_kmers");
    config.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: None,
        by_grouped_bed: Some(grouped_bed),
    });
    config.set_assign_by(WindowAssigner::CountOverlap);

    // Act
    run(&config)?;

    // Assert
    let loaded = load_ref_kmers_output(
        output_dir
            .path()
            .join("unit_grouped_bed_count_overlap_ref_kmers.ref_kmer_counts.zarr"),
    )?;
    let alpha_idx = loaded.group_index("alpha")?;
    let beta_idx = loaded.group_index("beta")?;
    assert_close(loaded.row_scaling_factor(alpha_idx).unwrap(), 7.50);
    assert_close(loaded.row_scaling_factor(beta_idx).unwrap(), 3.00);
    assert_close(loaded.count_for_motif(alpha_idx, "AAAA")?.unwrap(), 7.50);
    assert_close(loaded.count_for_motif(beta_idx, "AAAA")?.unwrap(), 3.00);
    assert_close(loaded.frequency_for_motif(alpha_idx, "AAAA")?.unwrap(), 1.0);
    assert_close(loaded.frequency_for_motif(beta_idx, "AAAA")?.unwrap(), 1.0);

    let groups = loaded.group_metadata()?;
    let alpha = groups
        .iter()
        .find(|group| group.name == "alpha")
        .expect("alpha group should be present");
    let beta = groups
        .iter()
        .find(|group| group.name == "beta")
        .expect("beta group should be present");
    assert_eq!(alpha.eligible_windows, 2);
    assert_eq!(beta.eligible_windows, 1);

    Ok(())
}

#[test]
fn ref_kmers_loader_reconstructs_dense_all_motifs_counts() -> Result<()> {
    // Arrange:
    // Reference chr1 is ACGT and k = 1. The complete motif axis is A, C, G, T. Each motif occurs
    // once, so the global row has scaling factor 4 and frequency 1/4 for each motif.
    let reference = twobit_from_sequences(
        "ref_kmers_dense_loader",
        vec![("chr1".to_string(), "ACGT".to_string())],
    )?;
    let output_dir = TempDir::new()?;
    let mut config = ref_kmers_config(&reference.path, output_dir.path(), 1);
    config.set_output_prefix("unit_dense_ref_kmers");
    config.set_all_motifs(true);

    // Act
    run(&config)?;

    // Assert
    let package_path = output_dir
        .path()
        .join("unit_dense_ref_kmers.ref_kmer_counts.zarr");
    let loaded = load_ref_kmers_output(&package_path)?;
    assert_eq!(loaded.storage_mode(), RefKmerStorageMode::Dense);
    assert_eq!(loaded.row_mode(), RefKmerRowMode::Global);
    assert_eq!(
        loaded.output_metadata().storage_mode,
        RefKmerStorageMode::Dense
    );
    assert!(loaded.sparse_frequencies().is_err());
    assert!(loaded.sparse_count_entries().is_err());
    match loaded.data() {
        RefKmerFrequencyData::Dense(frequencies) => assert_eq!(frequencies.shape(), (1, 4)),
        other => panic!("expected dense frequency data, got {other:?}"),
    }
    assert_eq!(
        loaded.motif_labels(),
        &[
            "A".to_string(),
            "C".to_string(),
            "G".to_string(),
            "T".to_string()
        ]
    );
    assert_eq!(loaded.row_scaling_factors(), &[4.0]);

    let dense_frequencies = loaded.dense_frequencies()?;
    assert_eq!(dense_frequencies.shape(), (1, 4));
    for motif_label in ["A", "C", "G", "T"] {
        assert_close(
            loaded.frequency_for_motif(0, motif_label)?.unwrap(),
            1.0 / 4.0,
        );
        assert_close(loaded.count_for_motif(0, motif_label)?.unwrap(), 1.0);
    }
    assert_eq!(
        loaded.to_dense_count_matrix()?.values_row_major(),
        &[1.0, 1.0, 1.0, 1.0]
    );

    let selected = loaded.select().motifs_by_label(&["T", "A"]).read()?;
    assert_eq!(selected.storage_mode(), RefKmerStorageMode::Dense);
    assert_eq!(selected.row_mode(), RefKmerRowMode::Global);
    assert_eq!(selected.motif_axis_kind(), RefKmerMotifAxisKind::Motif);
    assert_eq!(selected.kmer_size(), 1);
    assert!(!selected.canonical());
    assert!(selected.source_all_motifs());
    assert_eq!(selected.assign_by(), "count-overlap");
    assert_eq!(selected.row_indices(), &[0]);
    assert_eq!(selected.motif_indices(), &[3, 0]);
    assert_eq!(selected.motif_labels(), &["T".to_string(), "A".to_string()]);
    assert_eq!(selected.row_scaling_factors(), &[4.0]);
    assert_eq!(
        selected.dense_frequencies()?.values_row_major(),
        &[1.0 / 4.0, 1.0 / 4.0]
    );
    assert_eq!(
        selected.to_dense_count_matrix()?.values_row_major(),
        &[1.0, 1.0]
    );
    let global_row_error = loaded
        .select()
        .rows(&[0])
        .read()
        .expect_err("global row selector should fail");
    assert!(global_row_error.to_string().contains("global"));

    Ok(())
}

#[test]
fn ref_kmers_small_tiles_match_single_tile_counts() -> Result<()> {
    // Arrange:
    // Reference chr1 is ACGTACGTAC and k = 4. Fixed windows are [0,5) and [5,10).
    // The seven 4-mers are ACGT, CGTA, GTAC, TACG, ACGT, CGTA, GTAC.
    //
    // Under count-overlap:
    //   row 0: ACGT=1+1/4, CGTA=1, GTAC=3/4, TACG=1/2, total=3.5.
    //   row 1: ACGT=3/4, CGTA=1, GTAC=1/4+1, TACG=1/2, total=3.5.
    //
    // The small-tile run cuts the reference at coordinate 5, so starts 4 and 5 exercise
    // k-mer spans that cross the processing boundary. The public counts should be independent
    // of that internal tiling.
    let reference = twobit_from_sequences(
        "ref_kmers_tiled_equivalence",
        vec![("chr1".to_string(), "ACGTACGTAC".to_string())],
    )?;
    let output_dir = TempDir::new()?;

    let run_by_size = |output_prefix: &str, tile_size: Option<u32>| -> Result<_> {
        let mut config = ref_kmers_config(&reference.path, output_dir.path(), 4);
        config.set_output_prefix(output_prefix);
        config.set_windows(DistributionWindowsArgs {
            by_size: Some(5),
            by_bed: None,
            by_grouped_bed: None,
        });
        if let Some(tile_size) = tile_size {
            config.set_tile_size(tile_size);
        }
        run(&config)?;
        load_ref_kmers_output(
            output_dir
                .path()
                .join(format!("{output_prefix}.ref_kmer_counts.zarr")),
        )
        .map_err(anyhow::Error::from)
    };

    // Act
    let single_tile = run_by_size("single_tile_ref_kmers", None)?;
    let small_tiles = run_by_size("small_tile_ref_kmers", Some(5))?;

    // Assert
    for output in [&single_tile, &small_tiles] {
        assert_eq!(output.row_scaling_factors(), &[3.5, 3.5]);
        assert_close(output.count_for_motif(0, "ACGT")?.unwrap(), 1.25);
        assert_close(output.count_for_motif(0, "CGTA")?.unwrap(), 1.0);
        assert_close(output.count_for_motif(0, "GTAC")?.unwrap(), 0.75);
        assert_close(output.count_for_motif(0, "TACG")?.unwrap(), 0.5);
        assert_close(output.count_for_motif(1, "ACGT")?.unwrap(), 0.75);
        assert_close(output.count_for_motif(1, "CGTA")?.unwrap(), 1.0);
        assert_close(output.count_for_motif(1, "GTAC")?.unwrap(), 1.25);
        assert_close(output.count_for_motif(1, "TACG")?.unwrap(), 0.5);
    }
    assert_eq!(
        single_tile.to_dense_count_matrix()?.values_row_major(),
        small_tiles.to_dense_count_matrix()?.values_row_major()
    );

    Ok(())
}

#[test]
fn ref_kmers_loader_count_conversions_leave_frequency_values_unchanged() -> Result<()> {
    // Arrange:
    // Sparse case: AACC with k = 1 and windows [0,2), [2,4) stores only A in row 0 and C in row 1.
    // Each row has count 2, so dense frequencies are [1, 0, 0, 1] and dense counts are [2, 0, 0, 2].
    let sparse_reference = twobit_from_sequences(
        "ref_kmers_sparse_conversion",
        vec![("chr1".to_string(), "AACC".to_string())],
    )?;
    let sparse_output_dir = TempDir::new()?;
    let mut sparse_config = ref_kmers_config(&sparse_reference.path, sparse_output_dir.path(), 1);
    sparse_config.set_output_prefix("unit_sparse_conversion");
    sparse_config.set_windows(DistributionWindowsArgs {
        by_size: Some(2),
        by_bed: None,
        by_grouped_bed: None,
    });
    run(&sparse_config)?;
    let sparse_loaded = load_ref_kmers_output(
        sparse_output_dir
            .path()
            .join("unit_sparse_conversion.ref_kmer_counts.zarr"),
    )?;
    let sparse_stored_frequencies = sparse_loaded.sparse_frequencies()?.frequencies().to_vec();

    // Act
    let sparse_counts = sparse_loaded.to_dense_count_matrix()?;

    // Assert
    assert_eq!(sparse_stored_frequencies, vec![1.0, 1.0]);
    assert_eq!(sparse_counts.values_row_major(), &[2.0, 0.0, 0.0, 2.0]);
    assert_eq!(
        sparse_loaded.sparse_frequencies()?.frequencies(),
        sparse_stored_frequencies
    );
    assert_eq!(
        sparse_loaded
            .to_dense_frequency_matrix()?
            .values_row_major(),
        &[1.0, 0.0, 0.0, 1.0]
    );

    // Arrange:
    // Dense case: ACGT with k = 1 and all motifs stores one global row. Each motif count is 1
    // out of a row total of 4, so dense frequencies are 1/4 and dense counts are 1.
    let dense_reference = twobit_from_sequences(
        "ref_kmers_dense_conversion",
        vec![("chr1".to_string(), "ACGT".to_string())],
    )?;
    let dense_output_dir = TempDir::new()?;
    let mut dense_config = ref_kmers_config(&dense_reference.path, dense_output_dir.path(), 1);
    dense_config.set_output_prefix("unit_dense_conversion");
    dense_config.set_all_motifs(true);
    run(&dense_config)?;
    let dense_loaded = load_ref_kmers_output(
        dense_output_dir
            .path()
            .join("unit_dense_conversion.ref_kmer_counts.zarr"),
    )?;
    let dense_stored_frequencies = dense_loaded
        .dense_frequencies()?
        .values_row_major()
        .to_vec();

    // Act
    let dense_counts = dense_loaded.to_dense_count_matrix()?;
    let dense_selected = dense_loaded.select().motifs_by_label(&["T", "A"]).read()?;
    let selected_stored_frequencies = dense_selected
        .dense_frequencies()?
        .values_row_major()
        .to_vec();
    let selected_counts = dense_selected.to_dense_count_matrix()?;

    // Assert
    assert_slice_close(
        &dense_stored_frequencies,
        &[1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0],
    );
    assert_eq!(dense_counts.values_row_major(), &[1.0, 1.0, 1.0, 1.0]);
    assert_slice_close(
        dense_loaded.dense_frequencies()?.values_row_major(),
        &dense_stored_frequencies,
    );
    assert_slice_close(
        dense_loaded.to_dense_frequency_matrix()?.values_row_major(),
        &dense_stored_frequencies,
    );
    assert_slice_close(&selected_stored_frequencies, &[1.0 / 4.0, 1.0 / 4.0]);
    assert_eq!(selected_counts.values_row_major(), &[1.0, 1.0]);
    assert_slice_close(
        dense_selected.dense_frequencies()?.values_row_major(),
        &selected_stored_frequencies,
    );
    assert_slice_close(
        dense_selected
            .to_dense_frequency_matrix()?
            .values_row_major(),
        &selected_stored_frequencies,
    );

    Ok(())
}

#[test]
fn ref_kmers_bed_count_overlap_matches_manual_counts_across_tiles() -> Result<()> {
    // Arrange:
    // Reference chr1 is fourteen A bases and k = 4, so all eleven starts are AAAA.
    //
    // Manual count-overlap row totals:
    //   [2,6):  2/4 + 3/4 + 4/4 + 3/4 + 2/4 + 1/4 = 3.75.
    //   [6,9):  1/4 + 2/4 + 3/4 + 3/4 + 2/4 + 1/4 = 3.00.
    //   [8,12): 1/4 + 2/4 + 3/4 + 4/4 + 3/4 + 2/4 = 3.75.
    let reference = twobit_from_sequences(
        "ref_kmers_bed_tiled_count_overlap",
        vec![("chr1".to_string(), "A".repeat(14))],
    )?;
    let output_dir = TempDir::new()?;
    let windows_bed = output_dir.path().join("count_overlap_windows.bed");
    write_bed4(
        &windows_bed,
        &[
            Bed4Row::new("chr1", 2, 6, "left"),
            Bed4Row::new("chr1", 6, 9, "middle"),
            Bed4Row::new("chr1", 8, 12, "right"),
        ],
    )?;
    let run_by_bed = |output_prefix: &str, tile_size: Option<u32>| -> Result<_> {
        let mut config = ref_kmers_config(&reference.path, output_dir.path(), 4);
        config.set_output_prefix(output_prefix);
        config.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(windows_bed.clone()),
            by_grouped_bed: None,
        });
        config.set_assign_by(WindowAssigner::CountOverlap);
        if let Some(tile_size) = tile_size {
            config.set_tile_size(tile_size);
        }
        run(&config)?;
        load_ref_kmers_output(
            output_dir
                .path()
                .join(format!("{output_prefix}.ref_kmer_counts.zarr")),
        )
        .map_err(anyhow::Error::from)
    };

    // Act
    let single_tile = run_by_bed("single_tile_bed_ref_kmers", None)?;
    let small_tiles = run_by_bed("small_tile_bed_ref_kmers", Some(5))?;

    // Assert
    for output in [&single_tile, &small_tiles] {
        assert_slice_close(output.row_scaling_factors(), &[3.75, 3.00, 3.75]);
        assert_close(output.count_for_motif(0, "AAAA")?.unwrap(), 3.75);
        assert_close(output.count_for_motif(1, "AAAA")?.unwrap(), 3.00);
        assert_close(output.count_for_motif(2, "AAAA")?.unwrap(), 3.75);
        for row_idx in 0..3 {
            assert_close(output.frequency_for_motif(row_idx, "AAAA")?.unwrap(), 1.0);
        }
        assert_eq!(
            output
                .window_metadata()?
                .iter()
                .map(|window| (window.chrom.as_str(), window.interval.as_tuple()))
                .collect::<Vec<_>>(),
            vec![("chr1", (2, 6)), ("chr1", (6, 9)), ("chr1", (8, 12))]
        );
    }
    assert_eq!(
        single_tile.to_dense_count_matrix()?.values_row_major(),
        small_tiles.to_dense_count_matrix()?.values_row_major()
    );

    Ok(())
}

#[test]
fn ref_kmers_proportion_assignment_counts_each_passing_kmer_once() -> Result<()> {
    // Arrange:
    // Reference chr1 is ACGTA and k = 3. The BED window is [1,3).
    //   ACG [0,3) overlaps 2/3 and passes proportion=2/3.
    //   CGT [1,4) overlaps 2/3 and passes proportion=2/3.
    //   GTA [2,5) overlaps 1/3 and fails.
    // Passing k-mers contribute 1.0 each, not their overlap fraction.
    let reference = twobit_from_sequences(
        "ref_kmers_proportion_assignment",
        vec![("chr1".to_string(), "ACGTA".to_string())],
    )?;
    let output_dir = TempDir::new()?;
    let windows_bed = output_dir.path().join("proportion_windows.bed");
    write_bed4(&windows_bed, &[Bed4Row::new("chr1", 1, 3, "target")])?;
    let mut config = ref_kmers_config(&reference.path, output_dir.path(), 3);
    config.set_output_prefix("unit_proportion_ref_kmers");
    config.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
        by_grouped_bed: None,
    });
    config.set_assign_by(WindowAssigner::Proportion(2.0 / 3.0));
    config.set_all_motifs(true);

    // Act
    run(&config)?;

    // Assert
    let loaded = load_ref_kmers_output(
        output_dir
            .path()
            .join("unit_proportion_ref_kmers.ref_kmer_counts.zarr"),
    )?;
    assert_eq!(loaded.row_scaling_factors(), &[2.0]);
    assert_close(loaded.count_for_motif(0, "ACG")?.unwrap(), 1.0);
    assert_close(loaded.count_for_motif(0, "CGT")?.unwrap(), 1.0);
    assert_close(loaded.count_for_motif(0, "GTA")?.unwrap(), 0.0);
    assert_close(loaded.frequency_for_motif(0, "ACG")?.unwrap(), 0.5);
    assert_close(loaded.frequency_for_motif(0, "CGT")?.unwrap(), 0.5);

    Ok(())
}

#[test]
fn ref_kmers_midpoint_assignment_uses_the_center_base() -> Result<()> {
    // Arrange:
    // Reference chr1 is ACGTA and k = 3. BED windows are [0,2) and [2,5).
    //   ACG [0,3) has center base at coordinate 1, so it goes to row 0.
    //   CGT [1,4) has center base at coordinate 2, so it goes to row 1.
    //   GTA [2,5) has center base at coordinate 3, so it goes to row 1.
    let reference = twobit_from_sequences(
        "ref_kmers_midpoint_assignment",
        vec![("chr1".to_string(), "ACGTA".to_string())],
    )?;
    let output_dir = TempDir::new()?;
    let windows_bed = output_dir.path().join("midpoint_windows.bed");
    write_bed4(
        &windows_bed,
        &[
            Bed4Row::new("chr1", 0, 2, "left"),
            Bed4Row::new("chr1", 2, 5, "right"),
        ],
    )?;
    let mut config = ref_kmers_config(&reference.path, output_dir.path(), 3);
    config.set_output_prefix("unit_midpoint_ref_kmers");
    config.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
        by_grouped_bed: None,
    });
    config.set_assign_by(WindowAssigner::Midpoint);
    config.set_all_motifs(true);

    // Act
    run(&config)?;

    // Assert
    let loaded = load_ref_kmers_output(
        output_dir
            .path()
            .join("unit_midpoint_ref_kmers.ref_kmer_counts.zarr"),
    )?;
    assert_eq!(loaded.row_scaling_factors(), &[1.0, 2.0]);
    assert_close(loaded.count_for_motif(0, "ACG")?.unwrap(), 1.0);
    assert_close(loaded.count_for_motif(0, "CGT")?.unwrap(), 0.0);
    assert_close(loaded.count_for_motif(0, "GTA")?.unwrap(), 0.0);
    assert_close(loaded.count_for_motif(1, "ACG")?.unwrap(), 0.0);
    assert_close(loaded.count_for_motif(1, "CGT")?.unwrap(), 1.0);
    assert_close(loaded.count_for_motif(1, "GTA")?.unwrap(), 1.0);

    Ok(())
}

#[test]
fn ref_kmers_large_k_motifs_file_counts_selected_subspace() -> Result<()> {
    // Arrange:
    // k = 30 is outside the full reference k-mer universe used without a motifs file. The motifs
    // file selects exactly two possible targets. The reference contains one A^30 k-mer and no C^30
    // k-mers, so all-motifs selected output should keep both columns with counts 1 and 0.
    let present_motif = "A".repeat(30);
    let absent_motif = "C".repeat(30);
    let reference = twobit_from_sequences(
        "ref_kmers_large_selected_subspace",
        vec![("chr1".to_string(), present_motif.clone())],
    )?;
    let output_dir = TempDir::new()?;
    let motifs_file = write_motifs_file(&format!("{present_motif}\n{absent_motif}\n"))?;
    let mut config = ref_kmers_config(&reference.path, output_dir.path(), 30);
    config.set_output_prefix("unit_large_k_ref_kmers");
    config.set_motifs_file(Some(motifs_file.path().to_path_buf()));
    config.set_all_motifs(true);

    // Act
    run(&config)?;

    // Assert
    let loaded = load_ref_kmers_output(
        output_dir
            .path()
            .join("unit_large_k_ref_kmers.ref_kmer_counts.zarr"),
    )?;
    assert_eq!(
        loaded.motif_labels(),
        &[present_motif.clone(), absent_motif.clone()]
    );
    assert_eq!(loaded.row_scaling_factors(), &[1.0]);
    assert_close(
        loaded.count_for_motif(0, present_motif.as_str())?.unwrap(),
        1.0,
    );
    assert_close(
        loaded.count_for_motif(0, absent_motif.as_str())?.unwrap(),
        0.0,
    );

    Ok(())
}

#[test]
fn ref_kmers_blacklist_excludes_kmers_touching_masked_bases() -> Result<()> {
    // Arrange:
    // Reference chr1 is ACGTAC and k = 2. The blacklist covers base [2,3).
    //   AC [0,2) is outside the blacklist and counts.
    //   CG [1,3) and GT [2,4) touch the blacklisted base and are excluded.
    //   TA [3,5) and AC [4,6) count.
    let reference = twobit_from_sequences(
        "ref_kmers_blacklist_exclusion",
        vec![("chr1".to_string(), "ACGTAC".to_string())],
    )?;
    let output_dir = TempDir::new()?;
    let windows_bed = output_dir.path().join("blacklist_window.bed");
    write_bed4(&windows_bed, &[Bed4Row::new("chr1", 0, 6, "window")])?;
    let blacklist_bed = output_dir.path().join("blacklist.bed");
    write_bed4(&blacklist_bed, &[Bed4Row::new("chr1", 2, 3, "masked_base")])?;
    let mut config = ref_kmers_config(&reference.path, output_dir.path(), 2);
    config.set_output_prefix("unit_blacklisted_ref_kmers");
    config.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
        by_grouped_bed: None,
    });
    config.set_blacklist(Some(vec![blacklist_bed]));
    config.set_all_motifs(true);

    // Act
    run(&config)?;

    // Assert
    let loaded = load_ref_kmers_output(
        output_dir
            .path()
            .join("unit_blacklisted_ref_kmers.ref_kmer_counts.zarr"),
    )?;
    assert_eq!(loaded.row_scaling_factors(), &[3.0]);
    assert_close(loaded.count_for_motif(0, "AC")?.unwrap(), 2.0);
    assert_close(loaded.count_for_motif(0, "CG")?.unwrap(), 0.0);
    assert_close(loaded.count_for_motif(0, "GT")?.unwrap(), 0.0);
    assert_close(loaded.count_for_motif(0, "TA")?.unwrap(), 1.0);
    let windows = loaded.window_metadata()?;
    assert_close(windows[0].blacklisted_fraction.unwrap(), 1.0 / 6.0);

    Ok(())
}

#[test]
fn ref_kmers_selected_motifs_keep_empty_rows_without_unselected_denominator() -> Result<()> {
    // Arrange:
    // Reference chr1 is ACGT and k = 2. BED windows are [0,2) and [2,4).
    // The motifs file selects only AC, and all assignment requires the complete k-mer span inside
    // a window. Row 0 contains AC once. Row 1 contains GT, but GT is unselected and must not create
    // a denominator for selected-motif frequencies.
    let reference = twobit_from_sequences(
        "ref_kmers_empty_selected_row",
        vec![("chr1".to_string(), "ACGT".to_string())],
    )?;
    let output_dir = TempDir::new()?;
    let windows_bed = output_dir.path().join("selected_empty_windows.bed");
    write_bed4(
        &windows_bed,
        &[
            Bed4Row::new("chr1", 0, 2, "has_ac"),
            Bed4Row::new("chr1", 2, 4, "no_selected_motif"),
        ],
    )?;
    let motifs_file = write_motifs_file("AC\n")?;
    let mut config = ref_kmers_config(&reference.path, output_dir.path(), 2);
    config.set_output_prefix("unit_empty_selected_row_ref_kmers");
    config.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
        by_grouped_bed: None,
    });
    config.set_assign_by(WindowAssigner::All);
    config.set_motifs_file(Some(motifs_file.path().to_path_buf()));
    config.set_all_motifs(true);

    // Act
    run(&config)?;

    // Assert
    let loaded = load_ref_kmers_output(
        output_dir
            .path()
            .join("unit_empty_selected_row_ref_kmers.ref_kmer_counts.zarr"),
    )?;
    assert_eq!(loaded.motif_labels(), &["AC".to_string()]);
    assert_eq!(loaded.row_scaling_factors(), &[1.0, 0.0]);
    assert_close(loaded.count_for_motif(0, "AC")?.unwrap(), 1.0);
    assert_close(loaded.frequency_for_motif(0, "AC")?.unwrap(), 1.0);
    assert_close(loaded.count_for_motif(1, "AC")?.unwrap(), 0.0);
    assert_close(loaded.frequency_for_motif(1, "AC")?.unwrap(), 0.0);
    assert_eq!(
        loaded
            .window_metadata()?
            .iter()
            .map(|window| (window.chrom.as_str(), window.interval.as_tuple()))
            .collect::<Vec<_>>(),
        vec![("chr1", (0, 2)), ("chr1", (2, 4))]
    );

    Ok(())
}

#[test]
fn ref_kmers_fixed_size_rows_are_offset_across_chromosomes() -> Result<()> {
    // Arrange:
    // k = 1 and fixed windows of width 2.
    //   chr1=AAAA gives row 0 A=2 and row 1 A=2.
    //   chr2=CCCC gives row 2 C=2 and row 3 C=2.
    // Rows should follow the selected chromosome order without reusing row indices per chromosome.
    let reference = twobit_from_sequences(
        "ref_kmers_multi_contig_size_rows",
        vec![
            ("chr1".to_string(), "AAAA".to_string()),
            ("chr2".to_string(), "CCCC".to_string()),
        ],
    )?;
    let output_dir = TempDir::new()?;
    let chromosome_args = ChromosomeArgs {
        chromosomes: Some(vec!["chr1".to_string(), "chr2".to_string()]),
        chromosomes_file: None,
    };
    let mut config = RefKmersConfig::new(
        reference.path.clone(),
        output_dir.path().to_path_buf(),
        1,
        chromosome_args,
    );
    config.set_output_prefix("unit_multi_contig_ref_kmers");
    config.set_n_threads(1);
    config.set_windows(DistributionWindowsArgs {
        by_size: Some(2),
        by_bed: None,
        by_grouped_bed: None,
    });
    config.set_all_motifs(true);

    // Act
    run(&config)?;

    // Assert
    let loaded = load_ref_kmers_output(
        output_dir
            .path()
            .join("unit_multi_contig_ref_kmers.ref_kmer_counts.zarr"),
    )?;
    assert_eq!(loaded.row_scaling_factors(), &[2.0, 2.0, 2.0, 2.0]);
    assert_close(loaded.count_for_motif(0, "A")?.unwrap(), 2.0);
    assert_close(loaded.count_for_motif(1, "A")?.unwrap(), 2.0);
    assert_close(loaded.count_for_motif(2, "C")?.unwrap(), 2.0);
    assert_close(loaded.count_for_motif(3, "C")?.unwrap(), 2.0);
    assert_close(loaded.count_for_motif(0, "C")?.unwrap(), 0.0);
    assert_close(loaded.count_for_motif(2, "A")?.unwrap(), 0.0);
    assert_eq!(
        loaded
            .window_metadata()?
            .iter()
            .map(|window| (window.chrom.as_str(), window.interval.as_tuple()))
            .collect::<Vec<_>>(),
        vec![
            ("chr1", (0, 2)),
            ("chr1", (2, 4)),
            ("chr2", (0, 2)),
            ("chr2", (2, 4)),
        ]
    );

    Ok(())
}

fn read_json(path: &Path) -> Value {
    serde_json::from_str(&std::fs::read_to_string(path).expect("JSON file should read"))
        .expect("JSON should parse")
}

fn read_f64_array(store_path: &Path, array_path: &str) -> Vec<f64> {
    read_array(store_path, array_path)
}

fn read_i32_array(store_path: &Path, array_path: &str) -> Vec<i32> {
    read_array(store_path, array_path)
}

fn read_i64_array(store_path: &Path, array_path: &str) -> Vec<i64> {
    read_array(store_path, array_path)
}

fn read_u8_array(store_path: &Path, array_path: &str) -> Vec<u8> {
    read_array(store_path, array_path)
}

fn read_array<T>(store_path: &Path, array_path: &str) -> Vec<T>
where
    T: zarrs::array::ElementOwned,
{
    let store = Arc::new(FilesystemStore::new(store_path).expect("Zarr store should open"));
    let array = Array::open(store, array_path).expect("Zarr array should open");
    array
        .retrieve_array_subset(&array.subset_all())
        .expect("Zarr array should read")
}
