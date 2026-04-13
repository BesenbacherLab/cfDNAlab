use crate::{
    commands::{
        cli_common::{ensure_output_dir, resolve_chromosomes_and_contigs},
        coverage_weights::scaling_weights_config::ScalingWeightsArgs,
        coverage_weights::striding::{
            StrideBin, fill_triangular_overlap, normalize_avg_overlap_by_global_mean,
        },
        fcoverage::{
            config::FCoverageConfig,
            fcoverage::{FCoverageRunResult, run_inner as fcoverage_run_inner},
            window_results::CoverageWindowAction,
        },
    },
    shared::{
        interval::Interval,
        io::{dot_join, open_text_reader},
        tiled_run::make_temp_dir,
    },
};
use anyhow::{Context, Result, bail, ensure};
use fxhash::FxHashMap;
use std::{
    fs::File,
    io::{BufRead, BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};

const FCOVERAGE_INTERMEDIATE_DECIMALS: u8 = 12;

/// Calculates weights for genomic smoothing using large bins and a stride.
///
/// Technical details:
/// - Reuses `fcoverage --by-size <stride> --per-window average` as the raw counting step so
///   fragment handling, GC correction, blacklisting, and tiling stay consistent with `fcoverage`.
/// - Reads the resulting stride-bin averages back from disk, smooths them with a triangular
///   kernel, and writes the final scaling factors as TSV.
/// - Tracks the internal `fcoverage` counters so the printed summary reflects the fragments
///   that contributed to the scaling factors.
///
/// Parameters
/// ----------
/// - `opt`:
///     Fully resolved configuration for the `coverage-weights` command.
///
/// Returns
/// -------
/// - `Ok(())`:
///     Scaling factors were written successfully.
///
/// Errors
/// ------
/// - Returns an error if internal `fcoverage` fails, the intermediate TSV is malformed, or the
///   final scaling output cannot be written.
pub fn run(opt: &crate::commands::coverage_weights::config::CoverageWeightsConfig) -> Result<()> {
    run_with_fcoverage(&opt.shared, false, "coverage.scaling_factors.tsv")
}

/// Shared implementation for coverage-based and count-normalized scaling weights.
///
/// The future count-based command can call this with `normalize_by_length = true`
/// to reuse the exact same `fcoverage`-based counting path.
pub(crate) fn run_with_fcoverage(
    opt: &ScalingWeightsArgs,
    normalize_by_length: bool,
    output_file_name: &str,
) -> Result<()> {
    let start_time = Instant::now();
    let (chromosomes, _contigs) =
        resolve_chromosomes_and_contigs(&opt.chromosomes, opt.ioc.bam.as_path())?;
    opt.check_bin_sizes()?;

    // Keep all intermediate files under the user-chosen output directory so disk usage stays
    // within the filesystem location the user already selected for results.
    ensure_output_dir(&opt.ioc.output_dir)?;
    let fcoverage_output_dir = make_temp_dir(
        &opt.ioc.output_dir,
        &dot_join(&[opt.output_prefix.as_str(), "coverage_weights_source"]),
    )
    .context("creating internal fcoverage output directory")?;
    let _fcoverage_output_cleanup = RemoveDirOnDrop::new(fcoverage_output_dir.clone());

    let fcoverage_cfg =
        build_fcoverage_average_config(opt, &fcoverage_output_dir, normalize_by_length);

    println!("Calling fcoverage");
    println!("-----------------");
    let fcoverage_result =
        fcoverage_run_inner(&fcoverage_cfg).context("running internal fcoverage")?;
    println!("-----------------");
    println!("Reading fcoverage output");

    let mut bins_by_chr =
        load_stride_bins_from_fcoverage_average_tsv(&fcoverage_result, chromosomes.as_slice())?;

    for chromosome in &chromosomes {
        let bins = bins_by_chr
            .get_mut(chromosome)
            .with_context(|| format!("missing stride bins for chromosome '{}'", chromosome))?;
        fill_triangular_overlap(bins, opt.bin_size, opt.stride);
    }

    let global_avg_overlap_coverage =
        normalize_avg_overlap_by_global_mean(&mut bins_by_chr, true, true)?;

    println!(
        "Calculated the global average overlapping position-coverage: {}",
        global_avg_overlap_coverage
    );

    println!("Start: Writing stride-bin coordinates and scaling factors to disk");
    let file_name = dot_join(&[opt.output_prefix.as_str(), output_file_name]);
    let mut tsv_writer = BufWriter::new(
        File::create(opt.ioc.output_dir.join(file_name)).context("creating scaling-factors TSV")?,
    );
    writeln!(
        tsv_writer,
        "# gc_mode={}",
        crate::shared::scale_genome::scaling_gc_mode_for_run(
            opt.gc.gc_file.is_some(),
            opt.gc.gc_tag.is_some(),
        )
        .as_metadata_value()
    )
    .context("writing TSV metadata")?;
    writeln!(
        tsv_writer,
        "chromosome\tstart\tend\tavg_pos_cov\tavg_overlapping_pos_cov\tscaling_factor"
    )
    .context("writing TSV header")?;

    for chromosome in chromosomes {
        let bins = bins_by_chr
            .get(&chromosome)
            .with_context(|| format!("missing bins for chromosome: {}", chromosome))?;

        for bin in bins {
            writeln!(
                tsv_writer,
                "{}\t{}\t{}\t{}\t{}\t{}",
                chromosome,
                bin.start(),
                bin.end(),
                bin.avg_coverage,
                bin.avg_overlap_coverage,
                bin.scaling_factor
            )
            .context("writing TSV row")?;
        }
    }

    let global_counter = fcoverage_result.counters;
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
    if opt.gc.gc_file.is_some() || opt.gc.gc_tag.is_some() {
        let gc_fail_action = if opt.gc.skip_invalid_gc {
            "fragment skipped"
        } else {
            "fragment counted with weight 1.0"
        };
        println!(
            "  GC correction failures ({}): {}",
            gc_fail_action, global_counter.gc_failed_fragments
        );
    }
    if opt.gc.gc_tag.is_some() && global_counter.gc_out_of_range_tags > 0 {
        println!(
            "  GC tag values outside [0, {:.0}] treated as invalid: {}",
            crate::shared::gc_tag::MAX_REASONABLE_GC_WEIGHT,
            global_counter.gc_out_of_range_tags
        );
    }
    println!(
        "  Fragments counted one or more times: {}",
        global_counter.base.counted_fragments
    );
    println!("Elapsed time: {:.2?}", elapsed);
    Ok(())
}

fn build_fcoverage_average_config(
    opt: &ScalingWeightsArgs,
    output_dir: &Path,
    normalize_by_length: bool,
) -> FCoverageConfig {
    let mut cfg = FCoverageConfig::new(
        crate::commands::cli_common::IOCArgs {
            bam: opt.ioc.bam.clone(),
            output_dir: output_dir.to_path_buf(),
            n_threads: opt.ioc.n_threads,
        },
        opt.chromosomes.clone(),
    );

    cfg.set_unpaired(opt.unpaired.clone());
    cfg.set_normalize_by_length(normalize_by_length);
    cfg.set_output_prefix(opt.output_prefix.clone());
    cfg.set_decimals(FCOVERAGE_INTERMEDIATE_DECIMALS);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_windows(crate::commands::cli_common::WindowsArgs {
        by_size: Some(opt.stride as u64),
        by_bed: None,
    });
    cfg.set_fragment_lengths(opt.fragment_lengths.clone());
    cfg.set_min_mapq(opt.min_mapq);
    cfg.set_require_proper_pair(opt.require_proper_pair);
    cfg.set_blacklist(opt.blacklist.clone());
    cfg.set_gc(opt.gc.clone());
    cfg.set_ref_2bit(opt.ref_2bit.clone());
    cfg
}

fn load_stride_bins_from_fcoverage_average_tsv(
    fcoverage_result: &FCoverageRunResult,
    chromosomes: &[String],
) -> Result<FxHashMap<String, Vec<StrideBin>>> {
    let path = &fcoverage_result.final_out_path;
    let mut reader = open_text_reader(path)?;
    let mut line = String::new();

    line.clear();
    if reader.read_line(&mut line)? == 0 {
        bail!("{}: empty file; header required", path.display());
    }
    let header = line.trim_end();
    ensure!(
        header == "chromosome\tstart\tend\tavg_coverage\tblacklisted_positions",
        "{}: unexpected fcoverage header: '{}'",
        path.display(),
        header
    );

    let mut bins_by_chr =
        FxHashMap::with_capacity_and_hasher(chromosomes.len(), Default::default());

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            continue;
        }

        let cols: Vec<&str> = trimmed.split('\t').collect();
        ensure!(
            cols.len() == 5,
            "{}: expected 5 tab-separated columns, got {} in line '{}'",
            path.display(),
            cols.len(),
            trimmed
        );

        let chromosome = cols[0].to_string();
        let start: u32 = cols[1]
            .parse()
            .with_context(|| format!("{}: invalid start '{}'", path.display(), cols[1]))?;
        let end: u32 = cols[2]
            .parse()
            .with_context(|| format!("{}: invalid end '{}'", path.display(), cols[2]))?;
        let avg_coverage: f32 = cols[3]
            .parse()
            .with_context(|| format!("{}: invalid avg_coverage '{}'", path.display(), cols[3]))?;
        let _: u64 = cols[4].parse().with_context(|| {
            format!(
                "{}: invalid blacklisted_positions '{}'",
                path.display(),
                cols[4]
            )
        })?;

        bins_by_chr
            .entry(chromosome)
            .or_insert_with(Vec::new)
            .push(StrideBin {
                interval: Interval::new(start, end)?,
                avg_coverage,
                avg_overlap_coverage: 0.0,
                scaling_factor: 0.0,
            });
    }

    for chromosome in chromosomes {
        let bins = bins_by_chr
            .get(chromosome)
            .with_context(|| format!("{}: missing chromosome '{}'", path.display(), chromosome))?;
        ensure!(
            !bins.is_empty(),
            "{}: chromosome '{}' had no stride bins",
            path.display(),
            chromosome
        );
        ensure!(
            bins[0].start() == 0,
            "{}: chromosome '{}' did not start at 0",
            path.display(),
            chromosome
        );
        for pair in bins.windows(2) {
            ensure!(
                pair[0].end() == pair[1].start(),
                "{}: chromosome '{}' had non-contiguous stride bins at {}..{} and {}..{}",
                path.display(),
                chromosome,
                pair[0].start(),
                pair[0].end(),
                pair[1].start(),
                pair[1].end()
            );
        }
    }

    Ok(bins_by_chr)
}

struct RemoveDirOnDrop {
    path: PathBuf,
}

impl RemoveDirOnDrop {
    fn new(path: PathBuf) -> Self {
        Self { path }
    }
}

impl Drop for RemoveDirOnDrop {
    fn drop(&mut self) {
        if let Err(err) = std::fs::remove_dir_all(&self.path) {
            eprintln!(
                "warning: failed to remove temp dir {}: {}",
                self.path.display(),
                err
            );
        }
    }
}

#[cfg(test)]
mod tests {
    include!("coverage_weights_tests.rs");
}
