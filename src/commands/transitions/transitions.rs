use crate::{
    commands::{
        cli_common::ensure_output_dir,
        fragment_kmers::{config::*, fragment_kmers},
        transitions::config::TransitionsConfig,
    },
    shared::tiled_run::make_temp_dir,
};
use anyhow::{Context, Result};
use std::time::Instant;

/// Execute the base transition probability counting pipeline end-to-end.
///
/// Wraps fragment-kmers and calculates frequencies.
///
/// Parameters:
/// - `opt`: Fully resolved configuration for the `transitions` command.
///
/// Returns:
/// - `Ok(())` when the counts and accompanying metadata files are written successfully.
///
/// Errors:
/// - Propagates IO and parsing errors when reading inputs or writing results, aborting the run on
///   the first failure.
pub fn run(opt: &TransitionsConfig) -> Result<()> {
    let start_time = Instant::now();

    let prefix = opt.shared_args.output_prefix.trim();

    // Create output directory
    ensure_output_dir(&opt.shared_args.ioc.output_dir)?;

    // Build temporary directory
    let temp_dir = make_temp_dir(&opt.shared_args.ioc.output_dir, prefix)
        .context("create per-run temp dir")?;

    let mut tmp_ioc = opt.shared_args.ioc.clone();
    tmp_ioc.output_dir = temp_dir.join("/kmer_counts");

    let kmer_sizes = opt.orders.iter().map(|o| o + 1).collect();

    let fk_cfg = FragmentKmersConfig {
        shared_args: FragmentKmersSharedArgs::new(
            tmp_ioc,
            opt.shared_args.ref_genome.clone(),
            opt.shared_args.chromosomes.clone(),
            "fragment_kmers".to_string(),
        ),
        kmer_sizes: kmer_sizes,
        canonical: false,
        positional_counts: true,
        save_sparse: false, // TODO: Might be necessary later?
    };

    let global_counter = fragment_kmers::run_inner(&fk_cfg)?;

    // TODO: Load k-mer positional counts and calculate transition probabilities (frequencies)

    println!("");
    println!("Statistics");
    println!("----------");

    // Print summary statistics and execution time
    let elapsed = start_time.elapsed();
    println!("  Total reads: {}", global_counter.base.total_reads);
    println!(
        "  Initially accepted reads: {} ({:.2}%, forward: {}, reverse: {})",
        global_counter.base.accepted_forward + global_counter.base.accepted_reverse,
        (global_counter.base.accepted_forward + global_counter.base.accepted_reverse) as f64
            / global_counter.base.total_reads as f64
            * 100.0,
        global_counter.base.accepted_forward,
        global_counter.base.accepted_reverse
    );
    println!(
        "  Blacklist-excluded fragments: {}",
        global_counter.blacklisted_fragments
    );
    // if opt.gc.bin_by_gc {
    //     println!("GC-excluded reads: {}", global_counter.base.gc_excl);
    // }
    println!(
        "  Fragments counted one or more times: {}",
        global_counter.base.counted_fragments
    );
    println!("----------");
    println!("Elapsed time: {:.2?}", elapsed);
    Ok(())
}
