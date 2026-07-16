use super::*;
use crate::{
    command_run::RunOptions,
    commands::cli_common::{ChromosomeArgs, DistributionWindowsArgs, WindowAssigner},
    output_loaders::{
        load_ref_kmers_output, RefKmerMotifAxisKind, RefKmerRowMode, RefKmerStorageMode,
    },
    shared::kmers::kmer_codec::MAX_RADIX5_KMER_SIZE,
    testing::{twobit_from_sequences, write_bed4, Bed4Row},
};
use std::path::PathBuf;

fn base_config(kmer_size: u8, output_dir: PathBuf) -> RefKmersConfig {
    RefKmersConfig::new(
        PathBuf::from("missing-reference.2bit"),
        output_dir,
        kmer_size,
        ChromosomeArgs::default(),
    )
}

#[test]
fn rejects_zero_kmer_size_before_reading_reference() {
    // Arrange: the missing 2bit path would fail later, so seeing the k-mer size error proves
    // validation happened first.
    let output_dir = tempfile::tempdir().expect("temp output directory should be created");
    let config = base_config(0, output_dir.path().to_path_buf());

    // Act
    let error =
        run_ref_kmers(&config, RunOptions::new_quiet()).expect_err("zero k-mer size should fail");

    // Assert
    assert!(
        error
            .to_string()
            .contains("`--kmer-size` must be greater than 0"),
        "unexpected error: {error:#}"
    );
}

#[test]
fn rejects_zero_by_size_before_reading_reference() {
    // Arrange: programmatic configs bypass clap's positive range parser, so command validation must
    // still reject fixed windows that would otherwise divide by zero.
    let output_dir = tempfile::tempdir().expect("temp output directory should be created");
    let mut config = base_config(3, output_dir.path().to_path_buf());
    config.set_windows(DistributionWindowsArgs {
        by_size: Some(0),
        by_bed: None,
        by_grouped_bed: None,
    });

    // Act
    let error = run_ref_kmers(&config, RunOptions::new_quiet())
        .expect_err("zero fixed window size should fail");

    // Assert
    assert!(
        error
            .to_string()
            .contains("`--by-size` must be greater than 0"),
        "unexpected error: {error:#}"
    );
}

#[test]
fn rejects_large_full_space_kmer_size_before_reading_reference() {
    // Arrange: full-space radix-5 encoding is only supported through k = 27. The missing 2bit path
    // would fail later, so seeing the k-mer error proves validation happened first.
    let output_dir = tempfile::tempdir().expect("temp output directory should be created");
    let config = base_config(
        (MAX_RADIX5_KMER_SIZE + 1) as u8,
        output_dir.path().to_path_buf(),
    );

    // Act
    let error = run_ref_kmers(&config, RunOptions::new_quiet())
        .expect_err("large full-space k-mer size should fail");

    // Assert
    assert!(
        error.to_string().contains("requires `--motifs-file`"),
        "unexpected error: {error:#}"
    );
}

#[test]
fn allows_large_kmer_size_to_reach_motifs_file_parser() {
    // Arrange: a motifs file switches large k-mer counting to the selected-subspace path. The file
    // is intentionally missing so the expected failure point is unambiguous and before reference IO.
    let output_dir = tempfile::tempdir().expect("temp output directory should be created");
    let mut config = base_config(
        (MAX_RADIX5_KMER_SIZE + 1) as u8,
        output_dir.path().to_path_buf(),
    );
    config.set_motifs_file(Some(output_dir.path().join("missing-ref-kmers-motifs.tsv")));

    // Act
    let error = run_ref_kmers(&config, RunOptions::new_quiet())
        .expect_err("missing motifs file should fail after k-mer validation passes");

    // Assert
    let message = error.to_string();
    assert!(
        message.contains("reading ref-kmers motifs file"),
        "unexpected error: {error:#}"
    );
    assert!(
        !message.contains("requires `--motifs-file`"),
        "large k-mer validation should not reject runs with --motifs-file: {error:#}"
    );
}

#[test]
fn manualish_example() -> anyhow::Result<()> {
    let ref_seq_chrom1 = "ACGTGCAACCGGTTGGCCAGAGATATATCGCTCGTAACCAGGGTTTAAACCCAAAACCCCTTTTGGGG"; // 68 long
    let ref_seq_chrom2 = "GTGCAACCGGTTGGCCAGAGATATATCGCTCGTAACCAGGGTNTAAACCCAAAACCNCTTTTGNGGCA"; // 68 long, shifted 2 left with Ns after the first two windows
    let bed_starts = vec![2, 15, 40, 50];
    let bed_ends = vec![13, 32, 50, 68]; // exclusive
    let bed_groups = vec!["A", "A", "B", "B"];

    // This fixture checks BED, grouped BED, fixed-size, global, canonical, motifs-file, and
    // blacklist output from the same small reference. The expected counts below are hand-derived in
    // motif-axis order.
    let assert_close = |observed: f64, expected: f64| {
        assert!(
            (observed - expected).abs() < 1e-12,
            "observed {observed}, expected {expected}"
        );
    };
    let assert_counts_close = |observed: &[f64], expected: &[f64]| {
        assert_eq!(observed.len(), expected.len());
        for (motif_index, (observed_count, expected_count)) in
            observed.iter().zip(expected).enumerate()
        {
            assert!(
                (*observed_count - *expected_count).abs() < 1e-12,
                "motif index {motif_index}: observed {observed_count}, expected {expected_count}"
            );
        }
    };

    let reference = twobit_from_sequences(
        "ref_kmers_manual_example",
        vec![
            ("chr1".to_string(), ref_seq_chrom1.to_string()),
            ("chr2".to_string(), ref_seq_chrom2.to_string()),
        ],
    )?;
    let output_dir = tempfile::tempdir()?;
    let windows_bed = output_dir.path().join("manual_windows.bed");
    let mut bed_rows = Vec::new();
    for chromosome in ["chr1", "chr2"] {
        for window_index in 0..bed_starts.len() {
            bed_rows.push(Bed4Row::new(
                chromosome,
                bed_starts[window_index],
                bed_ends[window_index],
                bed_groups[window_index],
            ));
        }
    }
    write_bed4(&windows_bed, &bed_rows)?;

    let chromosome_args = ChromosomeArgs {
        chromosomes: Some(vec!["chr1".to_string(), "chr2".to_string()]),
        chromosomes_file: None,
    };
    let run_output = |output_prefix: &str, windows: DistributionWindowsArgs| {
        let mut config = RefKmersConfig::new(
            reference.path.clone(),
            output_dir.path().to_path_buf(),
            2,
            chromosome_args.clone(),
        );
        config.set_output_prefix(output_prefix);
        config.set_n_threads(1);
        config.set_windows(windows);
        config.set_assign_by(WindowAssigner::CountOverlap);
        config.set_all_motifs(true);
        let result = run_ref_kmers(&config, RunOptions::new_quiet())?;
        load_ref_kmers_output(result.ref_kmer_counts_path).map_err(anyhow::Error::from)
    };

    let expected_motif_labels = vec![
        "AA".to_string(),
        "AC".to_string(),
        "AG".to_string(),
        "AT".to_string(),
        "CA".to_string(),
        "CC".to_string(),
        "CG".to_string(),
        "CT".to_string(),
        "GA".to_string(),
        "GC".to_string(),
        "GG".to_string(),
        "GT".to_string(),
        "TA".to_string(),
        "TC".to_string(),
        "TG".to_string(),
        "TT".to_string(),
    ];
    // Count-overlap uses half-open k-mer intervals. In the final chr1 row [50,68), CC gets 0.5
    // from start 49 and 1.0 from starts 50, 56, 57, and 58. No start 67 is counted because
    // [67,69) crosses the chromosome end.
    let expected_bed_counts: [[f64; 16]; 8] = [
        [
            1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.5, 0.0, 0.0, 1.0, 1.0, 2.0, 0.0, 0.0, 1.0, 0.5,
        ],
        [
            0.0, 0.0, 2.0, 3.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 0.5, 0.0, 2.0, 1.5, 0.0, 0.0,
        ],
        [
            2.0, 1.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0, 1.0, 0.0, 0.0, 2.0,
        ],
        [
            3.0, 1.0, 0.0, 0.0, 1.0, 4.5, 0.0, 1.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 1.0, 3.0,
        ],
        [
            1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.5, 1.0, 0.0, 0.0, 1.5, 1.0,
        ],
        [
            0.0, 0.0, 2.0, 3.0, 1.0, 0.5, 2.0, 1.0, 2.0, 1.0, 0.0, 0.5, 2.0, 2.0, 0.0, 0.0,
        ],
        [
            2.0, 1.0, 0.0, 0.0, 0.5, 2.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 0.0, 0.0, 0.0,
        ],
        [
            3.0, 1.0, 0.0, 0.0, 1.5, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 3.0,
        ],
    ];
    let expected_bed_scaling = [11.0, 17.0, 10.0, 17.5, 11.0, 17.0, 8.0, 13.5];

    /* --by-bed */

    let bed_output = run_output(
        "manual_bed",
        DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(windows_bed.clone()),
            by_grouped_bed: None,
        },
    )?;
    assert_eq!(bed_output.storage_mode(), RefKmerStorageMode::Dense);
    assert_eq!(bed_output.row_mode(), RefKmerRowMode::BedWindows);
    assert_eq!(bed_output.motif_labels(), expected_motif_labels.as_slice());
    assert_eq!(
        bed_output.row_scaling_factors().len(),
        expected_bed_scaling.len()
    );
    for (observed_scaling, expected_scaling) in bed_output
        .row_scaling_factors()
        .iter()
        .zip(expected_bed_scaling)
    {
        assert_close(*observed_scaling, expected_scaling);
    }
    let expected_windows = [
        ("chr1", 2, 13),
        ("chr1", 15, 32),
        ("chr1", 40, 50),
        ("chr1", 50, 68),
        ("chr2", 2, 13),
        ("chr2", 15, 32),
        ("chr2", 40, 50),
        ("chr2", 50, 68),
    ];
    let bed_windows = bed_output.window_metadata()?;
    assert_eq!(bed_windows.len(), expected_windows.len());
    for (window, (chromosome, start, end)) in bed_windows.iter().zip(expected_windows) {
        assert_eq!(window.chrom.as_str(), chromosome);
        assert_eq!(window.interval.as_tuple(), (start, end));
        assert_close(window.blacklisted_fraction.unwrap(), 0.0);
    }
    let bed_counts = bed_output.to_dense_count_matrix()?;
    assert_eq!(bed_counts.shape(), (8, 16));
    for (row_index, expected_counts) in expected_bed_counts.iter().enumerate() {
        let observed_counts = bed_counts.row(row_index).expect("BED count row exists");
        assert_counts_close(observed_counts, expected_counts);
    }

    /* --by-grouped-bed */

    let sum_expected_rows = |row_indices: &[usize]| -> [f64; 16] {
        let mut totals = [0.0; 16];
        for &row_index in row_indices {
            for (motif_index, count) in expected_bed_counts[row_index].iter().enumerate() {
                totals[motif_index] += count;
            }
        }
        totals
    };
    let expected_group_a_counts = sum_expected_rows(&[0, 1, 4, 5]);
    let expected_group_b_counts = sum_expected_rows(&[2, 3, 6, 7]);
    let grouped_output = run_output(
        "manual_grouped_bed",
        DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(windows_bed.clone()),
        },
    )?;
    assert_eq!(grouped_output.storage_mode(), RefKmerStorageMode::Dense);
    assert_eq!(grouped_output.row_mode(), RefKmerRowMode::Groups);
    assert_eq!(grouped_output.row_scaling_factors(), &[56.0, 49.0]);
    let groups = grouped_output.group_metadata()?;
    assert_eq!(groups.len(), 2);
    assert_eq!(groups[0].name, "A");
    assert_eq!(groups[0].eligible_windows, 4);
    assert_close(groups[0].blacklisted_fraction, 0.0);
    assert_eq!(groups[1].name, "B");
    assert_eq!(groups[1].eligible_windows, 4);
    assert_close(groups[1].blacklisted_fraction, 0.0);
    let grouped_counts = grouped_output.to_dense_count_matrix()?;
    assert_eq!(grouped_counts.shape(), (2, 16));
    assert_counts_close(
        grouped_counts.row(0).expect("group A count row exists"),
        &expected_group_a_counts,
    );
    assert_counts_close(
        grouped_counts.row(1).expect("group B count row exists"),
        &expected_group_b_counts,
    );

    /* --by-size vs global */

    // Global counts skip only 2-mers that include an N. For chr2 this removes starts 41, 42, 55,
    // 56, 62, and 63. Starts 54=CC, 61=TG, and 66=CA remain valid.
    let expected_global_counts = [
        14.0, 9.0, 6.0, 6.0, 9.0, 14.0, 7.0, 4.0, 4.0, 7.0, 12.0, 8.0, 8.0, 4.0, 6.0, 10.0,
    ];
    let sum_count_rows = |counts: &crate::output_loaders::DenseMatrix<f64>| -> Vec<f64> {
        let mut totals = vec![0.0; counts.column_count()];
        for row in counts.rows() {
            for (motif_index, count) in row.iter().enumerate() {
                totals[motif_index] += count;
            }
        }
        totals
    };
    let by_size_output = run_output(
        "manual_by_size",
        DistributionWindowsArgs {
            by_size: Some(10),
            by_bed: None,
            by_grouped_bed: None,
        },
    )?;
    assert_eq!(by_size_output.row_mode(), RefKmerRowMode::SizeWindows);
    assert_eq!(by_size_output.row_scaling_factors().len(), 14);
    assert_close(
        by_size_output.row_scaling_factors().iter().sum::<f64>(),
        128.0,
    );
    let by_size_counts = by_size_output.to_dense_count_matrix()?;
    let by_size_total_counts = sum_count_rows(&by_size_counts);
    assert_counts_close(&by_size_total_counts, &expected_global_counts);

    let global_output = run_output(
        "manual_global",
        DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: None,
        },
    )?;
    assert_eq!(global_output.row_mode(), RefKmerRowMode::Global);
    assert_eq!(global_output.row_scaling_factors(), &[128.0]);
    let global_counts = global_output.to_dense_count_matrix()?;
    assert_eq!(global_counts.shape(), (1, 16));
    let observed_global_counts = global_counts.row(0).expect("global count row exists");
    assert_counts_close(observed_global_counts, &expected_global_counts);
    assert_counts_close(&by_size_total_counts, observed_global_counts);

    /* Canonical */

    let mut canonical_config = RefKmersConfig::new(
        reference.path.clone(),
        output_dir.path().to_path_buf(),
        2,
        chromosome_args.clone(),
    );
    canonical_config.set_output_prefix("manual_canonical");
    canonical_config.set_n_threads(1);
    canonical_config.set_all_motifs(true);
    canonical_config.set_canonical(true);
    let canonical_result = run_ref_kmers(&canonical_config, RunOptions::new_quiet())?;
    let canonical_output = load_ref_kmers_output(canonical_result.ref_kmer_counts_path)?;
    assert_eq!(canonical_output.row_mode(), RefKmerRowMode::Global);
    assert_eq!(canonical_output.storage_mode(), RefKmerStorageMode::Dense);
    assert!(canonical_output.canonical());
    // Canonical 2-mers are the lexicographically smaller member of each reverse-complement pair.
    // The expected counts come from `expected_global_counts` in `expected_motif_labels` order:
    // AA = AA + TT = 14 + 10 = 24
    // AC = AC + GT = 9 + 8 = 17
    // AG = AG + CT = 6 + 4 = 10
    // AT = AT = 6
    // CA = CA + TG = 9 + 6 = 15
    // CC = CC + GG = 14 + 12 = 26
    // CG = CG = 7
    // GA = GA + TC = 4 + 4 = 8
    // GC = GC = 7
    // TA = TA = 8
    assert_eq!(
        canonical_output.motif_labels(),
        &[
            "AA".to_string(),
            "AC".to_string(),
            "AG".to_string(),
            "AT".to_string(),
            "CA".to_string(),
            "CC".to_string(),
            "CG".to_string(),
            "GA".to_string(),
            "GC".to_string(),
            "TA".to_string()
        ]
    );
    assert_eq!(canonical_output.row_scaling_factors(), &[128.0]);
    let canonical_counts = canonical_output.to_dense_count_matrix()?;
    let expected_canonical_counts = [24.0, 17.0, 10.0, 6.0, 15.0, 26.0, 7.0, 8.0, 7.0, 8.0];
    assert_counts_close(
        canonical_counts
            .row(0)
            .expect("canonical global count row exists"),
        &expected_canonical_counts,
    );

    /* --motifs-file */

    let motifs_file = output_dir.path().join("manual_motifs.tsv");
    std::fs::write(
        &motifs_file,
        "AC\tedge\nGT\tedge\nCC\tgc_rich\nGG\tgc_rich\nAA\thomopolymer\nTT\thomopolymer\n",
    )?;
    let mut motifs_config = RefKmersConfig::new(
        reference.path.clone(),
        output_dir.path().to_path_buf(),
        2,
        chromosome_args.clone(),
    );
    motifs_config.set_output_prefix("manual_motifs_file");
    motifs_config.set_n_threads(1);
    motifs_config.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed.clone()),
        by_grouped_bed: None,
    });
    motifs_config.set_assign_by(WindowAssigner::CountOverlap);
    motifs_config.set_all_motifs(true);
    motifs_config.set_motifs_file(Some(motifs_file));
    let motifs_result = run_ref_kmers(&motifs_config, RunOptions::new_quiet())?;
    let motifs_output = load_ref_kmers_output(motifs_result.ref_kmer_counts_path)?;
    assert_eq!(motifs_output.row_mode(), RefKmerRowMode::BedWindows);
    assert_eq!(motifs_output.storage_mode(), RefKmerStorageMode::Dense);
    assert_eq!(
        motifs_output.motif_axis_kind(),
        RefKmerMotifAxisKind::MotifGroup
    );
    assert_eq!(
        motifs_output.motif_labels(),
        &[
            "edge".to_string(),
            "gc_rich".to_string(),
            "homopolymer".to_string()
        ]
    );
    let expected_motifs_file_scaling = [6.5, 1.5, 8.5, 14.5, 6.5, 1.0, 6.5, 9.0];
    assert_eq!(
        motifs_output.row_scaling_factors().len(),
        expected_motifs_file_scaling.len()
    );
    for (observed_scaling, expected_scaling) in motifs_output
        .row_scaling_factors()
        .iter()
        .zip(expected_motifs_file_scaling)
    {
        assert_close(*observed_scaling, expected_scaling);
    }
    let motifs_counts = motifs_output.to_dense_count_matrix()?;
    assert_eq!(motifs_counts.shape(), (8, 3));
    let expected_motifs_file_counts: [[f64; 3]; 8] = [
        [3.0, 2.0, 1.5],
        [0.0, 1.5, 0.0],
        [2.0, 2.5, 4.0],
        [1.0, 7.5, 6.0],
        [2.0, 2.5, 2.0],
        [0.5, 0.5, 0.0],
        [2.0, 2.5, 2.0],
        [1.0, 2.0, 6.0],
    ];
    for (row_index, expected_counts) in expected_motifs_file_counts.iter().enumerate() {
        let observed_counts = motifs_counts
            .row(row_index)
            .expect("motifs-file BED count row exists");
        assert_counts_close(observed_counts, expected_counts);
    }

    let blacklist_bed = output_dir.path().join("manual_blacklist.bed");
    write_bed4(
        &blacklist_bed,
        &[
            Bed4Row::new("chr1", 6, 7, "mask_chr1_first_a"),
            Bed4Row::new("chr1", 40, 41, "mask_chr1_g"),
            Bed4Row::new("chr2", 66, 67, "mask_chr2_c"),
        ],
    )?;
    let mut blacklist_config = RefKmersConfig::new(
        reference.path.clone(),
        output_dir.path().to_path_buf(),
        2,
        chromosome_args,
    );
    blacklist_config.set_output_prefix("manual_blacklist");
    blacklist_config.set_n_threads(1);
    blacklist_config.set_all_motifs(true);
    blacklist_config.set_blacklist(Some(vec![blacklist_bed]));
    let blacklist_result = run_ref_kmers(&blacklist_config, RunOptions::new_quiet())?;
    let blacklist_output = load_ref_kmers_output(blacklist_result.ref_kmer_counts_path)?;
    assert_eq!(blacklist_output.row_mode(), RefKmerRowMode::Global);
    assert_eq!(blacklist_output.row_scaling_factors(), &[122.0]);
    let expected_blacklisted_global_counts = [
        13.0, 9.0, 5.0, 6.0, 7.0, 14.0, 7.0, 4.0, 4.0, 6.0, 11.0, 8.0, 8.0, 4.0, 6.0, 10.0,
    ];
    let blacklist_counts = blacklist_output.to_dense_count_matrix()?;
    assert_counts_close(
        blacklist_counts
            .row(0)
            .expect("blacklisted global count row exists"),
        &expected_blacklisted_global_counts,
    );

    Ok(())
}
