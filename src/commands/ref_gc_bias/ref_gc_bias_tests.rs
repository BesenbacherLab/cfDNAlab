use super::*;

use crate::commands::cli_common::{FragmentLengthArgs, LoggingArgs, Ref2BitRequiredArgs};
use anyhow::anyhow;
use std::{
    fs::File,
    io::Write,
    path::{Path, PathBuf},
};
use tempfile::TempDir;
use twobit::convert::{fasta::FastaReader, to_2bit};

struct TwoBitFixture {
    _tempdir: TempDir,
    path: PathBuf,
}

fn write_fasta<P: AsRef<Path>>(path: P, sequences: &[(String, String)]) -> Result<()> {
    let mut file = File::create(path)?;
    for (name, seq) in sequences {
        writeln!(file, ">{name}")?;
        for chunk in seq.as_bytes().chunks(60) {
            file.write_all(chunk)?;
            file.write_all(b"\n")?;
        }
    }
    Ok(())
}

fn twobit_from_sequences(name: &str, sequences: Vec<(String, String)>) -> Result<TwoBitFixture> {
    let normalized: Vec<(String, String)> = sequences
        .into_iter()
        .map(|(chr, seq)| (chr, seq.to_ascii_uppercase()))
        .collect();
    let tempdir = TempDir::new()?;
    let fasta_path = tempdir.path().join(format!("{name}.fasta"));
    write_fasta(&fasta_path, &normalized)?;
    let path = tempdir.path().join(format!("{name}.2bit"));
    {
        let reader = FastaReader::open(&fasta_path).map_err(|e| anyhow!(e.to_string()))?;
        let mut file = File::create(&path)?;
        to_2bit(&mut file, &reader).map_err(|e| anyhow!(e.to_string()))?;
    }
    Ok(TwoBitFixture {
        _tempdir: tempdir,
        path,
    })
}

fn base_test_config(ref_2bit: PathBuf) -> RefGCBiasConfig {
    RefGCBiasConfig {
        ref_genome: Ref2BitRequiredArgs { ref_2bit },
        output_dir: PathBuf::from("out"),
        output_prefix: String::new(),
        n_threads: 1,
        n_positions: 100,
        seed: Some(1),
        windows: crate::commands::ref_gc_bias::config::RefGCWindowsArgs {
            by_bed: Some(PathBuf::from("windows.bed")),
        },
        chromosomes: ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        blacklist: None,
        fragment_lengths: FragmentLengthArgs {
            min_fragment_length: 10,
            max_fragment_length: 10,
        },
        end_offset: 0,
        skip_interpolation: true,
        smoothing_sigma: 0.55,
        smoothing_radius: 2,
        skip_smoothing: true,
        tile_size: 10,
        logging: LoggingArgs::default(),
    }
}

#[test]
fn process_tile_skips_a_halo_only_bed_window_in_core_overlap_mode() -> Result<()> {
    let reference = twobit_from_sequences(
        "ref_gc_bias_halo_only_tile_local",
        vec![("chr1".into(), "A".repeat(40))],
    )?;
    let cfg = base_test_config(reference.path.clone());
    let tile = Tile::from_coords("chr1".to_string(), 0, 1, 10, 20, 10, 20)?;
    let windows = vec![IndexedInterval::new(22, 32, 0)?];
    let span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 1,
    };
    let start_positions: Vec<usize> = (10..20).collect();

    let (counts, total_acgt_in_core) = process_tile(
        &tile,
        Some(&span),
        40,
        Some(&windows),
        &start_positions,
        &[],
        &cfg,
    )?;

    assert_eq!(counts.sum(), 0.0);
    assert_eq!(total_acgt_in_core, 0);
    Ok(())
}

#[test]
fn process_tile_counts_boundary_crossing_bed_window_from_core_start_through_fetch_context() -> Result<()> {
    // Manual derivation:
    // - tile core is [10,20) and loaded sequence is [10,30)
    // - BED window [18,30) overlaps the core and remains countable for this tile
    // - `ref_gc_bias` owns sampled starts by tile core, and those starts may still count through
    //   the full overlapping window instead of stopping at `core_end`
    // - starts 18 and 19 are therefore owned by this tile and each fit one length-10 fragment
    //   fully inside the counting window [18,30)
    // - counts.sum() is thus 2.0 for the single configured fragment length
    // - support bookkeeping still clips to the tile-owned core overlap [18,20), which has 2 ACGT
    //   bases on an all-A reference
    let reference =
        twobit_from_sequences("ref_gc_bias_boundary_clip", vec![("chr1".into(), "A".repeat(40))])?;
    let cfg = base_test_config(reference.path.clone());
    let tile = Tile::from_coords("chr1".to_string(), 0, 1, 10, 20, 10, 30)?;
    let windows = vec![IndexedInterval::new(18, 30, 0)?];
    let span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 1,
    };
    let start_positions: Vec<usize> = (10..20).collect();

    let (counts, total_acgt_in_core) = process_tile(
        &tile,
        Some(&span),
        40,
        Some(&windows),
        &start_positions,
        &[],
        &cfg,
    )?;

    assert_eq!(counts.sum(), 2.0);
    assert_eq!(total_acgt_in_core, 2);
    Ok(())
}

#[test]
fn process_tile_blacklist_uses_reference_coordinates_with_nonzero_sequence_origin() -> Result<()> {
    // Manual derivation:
    // - The tile loads reference slice [900,920), so prefix-local coordinate 0 maps to reference
    //   coordinate 900.
    // - The blacklist interval [905,910) should mask local [5,10), leaving only five ACGT bases
    //   in the tile-owned core [900,910).
    // - With fragment length 10 and min ACGT fraction 1.0, every possible core-owned start
    //   900..909 overlaps the masked [905,910) segment. No fragment has full ACGT support, so
    //   no reference GC count survives.
    // - If masking incorrectly used prefix-local origin 0, [905,910) would miss this 20 bp slice,
    //   `total_acgt_in_core` would be 10, and starts would be counted.
    let mut sequence = "A".repeat(900);
    sequence.push_str(&"C".repeat(20));
    let reference = twobit_from_sequences(
        "ref_gc_bias_blacklist_nonzero_origin",
        vec![("chr1".into(), sequence)],
    )?;
    let cfg = base_test_config(reference.path.clone());
    let tile = Tile::from_coords("chr1".to_string(), 0, 90, 900, 910, 900, 920)?;
    let windows = vec![IndexedInterval::new(900, 920, 0)?];
    let span = TileWindowSpan {
        first_idx: 0,
        last_idx_exclusive: 1,
    };
    let start_positions: Vec<usize> = (900..910).collect();
    let blacklist_intervals = vec![Interval::new(905_u64, 910_u64)?];

    let (counts, total_acgt_in_core) = process_tile(
        &tile,
        Some(&span),
        920,
        Some(&windows),
        &start_positions,
        &blacklist_intervals,
        &cfg,
    )?;

    assert_eq!(counts.sum(), 0.0);
    assert_eq!(total_acgt_in_core, 5);
    Ok(())
}
