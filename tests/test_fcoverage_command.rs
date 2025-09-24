mod fixtures;

use anyhow::Result;
use cfdnalab::cli_common::{ChromosomeArgs, IOCArgs, WindowsArgs};
use cfdnalab::fcoverage::{FCoverageConfig, run};
use cfdnalab::utils::coverage::window_results::CoverageWindowAction;
use cfdnalab::utils::fragment::minimal_fragment::collect_fragment_from_records;
use cfdnalab::utils::read::default_include_read;
use fixtures::{read_zst_to_string, simple_inward_bam};
use rust_htslib::bam::{Read, Reader};
use tempfile::TempDir;

fn base_chromosomes(chrs: &[&str]) -> ChromosomeArgs {
    ChromosomeArgs {
        chromosomes: Some(chrs.iter().map(|c| c.to_string()).collect()),
        chromosomes_file: None,
    }
}

#[test]
fn per_position_outputs_basic_fragment() -> Result<()> {
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;

    let mut cfg = FCoverageConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.path().to_path_buf(),
            n_threads: 2,
        },
        base_chromosomes(&["chr1"]),
    );
    cfg.set_output_prefix("testcov");
    cfg.set_decimals(2);
    cfg.set_tile_size(1_000);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_ignore_gap(false);
    cfg.set_keep_zero_runs(false);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 200;
    }

    let mut reader = Reader::from_path(&bam.bam)?;
    let mut accepted = (0, 0);
    let mut pair_store = Vec::new();
    for (idx, result) in reader.records().enumerate() {
        let rec = result?;
        if default_include_read(&rec, cfg.require_proper_pair, cfg.min_mapq) {
            if rec.is_reverse() {
                accepted.1 += 1;
            } else {
                accepted.0 += 1;
            }
            pair_store.push(rec);
        } else if idx == 0 {
            eprintln!(
                "forward read filtered: flags={:#x}, is_reverse={}, mate_reverse={}, mapq={}",
                rec.flags(),
                rec.is_reverse(),
                rec.is_mate_reverse(),
                rec.mapq()
            );
        }
    }
    assert_eq!(
        accepted,
        (1, 1),
        "expected both mates accepted, got {accepted:?}"
    );
    assert_eq!(pair_store.len(), 2);
    let frag = collect_fragment_from_records(&pair_store[0], &pair_store[1]);
    assert!(frag.is_some(), "expected fragment collection to succeed");

    run(cfg)?;

    let bedgraph = out_dir.path().join("testcov.per_position.bedgraph.zst");
    assert!(bedgraph.exists(), "expected positional bedgraph output");
    let text = read_zst_to_string(&bedgraph)?;
    assert!(
        text.contains("chr1\t20\t80\t1"),
        "expected contiguous coverage run, got: {text}"
    );

    Ok(())
}

#[test]
fn by_size_total_and_average_outputs() -> Result<()> {
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;

    let mut windows = WindowsArgs::default();
    windows.by_size = Some(40);

    let mut cfg = FCoverageConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.path().to_path_buf(),
            n_threads: 2,
        },
        base_chromosomes(&["chr1"]),
    );
    cfg.set_output_prefix("testcov");
    cfg.set_decimals(2);
    cfg.set_tile_size(1_000);
    cfg.set_per_window(CoverageWindowAction::Total);
    cfg.set_keep_zero_runs(true);
    cfg.set_ignore_gap(false);
    cfg.set_windows(windows);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 200;
    }

    run(cfg)?;

    let totals = out_dir.path().join("testcov.total.tsv.zst");
    assert!(totals.exists(), "expected per-window totals output");
    let text = read_zst_to_string(&totals)?;
    let mut lines = text.lines();
    let _header = lines.next().unwrap_or("");
    let first = lines.next().unwrap_or("");
    assert!(first.starts_with("chr1\t0\t40"));
    assert!(
        first.ends_with("20\t0"),
        "expected total coverage 20, got: {first}"
    );
    let second = lines.next().unwrap_or("");
    assert!(second.starts_with("chr1\t40\t80"));
    assert!(
        second.ends_with("40\t0"),
        "expected total coverage 40, got: {second}"
    );

    Ok(())
}
