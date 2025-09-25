use std::path::Path;

use anyhow::{Context, Result};
use fxhash::FxHashMap;

use crate::cli_common::{ChromosomeArgs, IOCArgs, ScaleGenomeArgs};
use crate::utils::bam::{Contigs, bam_contigs_info};
use crate::utils::blacklist::load_blacklists;
use crate::utils::coverage::scale_genome::load_scaling_factors_tsv;

/// Resolve chromosomes and BAM contig metadata once for a command.
///
/// Implementation details:
/// - Delegates to `ChromosomeArgs::resolve_chromosomes`, passing the BAM path so
///   aliases such as `--chromosomes all` work uniformly.
/// - Queries the BAM header via `bam_contigs_info` to obtain target lengths.
///
/// Parameters:
/// - `chrom_args`: Command-line chromosome selection configuration.
/// - `ioc`: Shared IO arguments providing the BAM path.
///
/// Returns:
/// - A tuple with the resolved chromosome names and their contig metadata.
///
/// Errors:
/// - Propagates IO and parsing failures when the BAM cannot be opened or lacks
///   the requested contigs.
pub fn resolve_chromosomes_and_contigs(
    chrom_args: &ChromosomeArgs,
    ioc: &IOCArgs,
) -> Result<(Vec<String>, Contigs)> {
    let chromosomes = chrom_args
        .resolve_chromosomes(Some(ioc.bam.as_path()))
        .context("resolve chromosomes")?;
    let contigs = bam_contigs_info(&ioc.bam, &chromosomes).context("fetch contig metadata")?;
    Ok((chromosomes, contigs))
}

/// Create the output directory if it does not exist.
///
/// Implementation details:
/// - Wraps `std::fs::create_dir_all` with an `anyhow` context to yield helpful
///   error messages tailored to the target path.
///
/// Parameters:
/// - `path`: Directory where the command should place its results.
///
/// Returns:
/// - `Ok(())` if the directory exists or was created successfully.
///
/// Errors:
/// - Returns an error when the directory cannot be created (for instance due to
///   missing permissions or an unwritable parent directory).
pub fn ensure_output_dir(path: &Path) -> Result<()> {
    std::fs::create_dir_all(path)
        .with_context(|| format!("cannot create output directory: {}", path.display()))
}

/// Load blacklist intervals when the user supplied one or more BED files.
///
/// Implementation details:
/// - Delegates to `load_blacklists`, which merges overlapping intervals per
///   chromosome and enforces `min_size` filtering.
/// - Returns an empty map when `beds` is `None` so callers can operate without
///   additional branching.
///
/// Parameters:
/// - `beds`: Optional list of BED paths.
/// - `min_size`: Minimum interval size (bp) to retain.
/// - `chromosomes`: Chromosomes the command intends to process.
///
/// Returns:
/// - A map keyed by chromosome name containing sorted blacklist intervals.
///
/// Errors:
/// - Propagates parsing errors if any BED file is malformed or unavailable.
pub fn load_blacklist_map(
    beds: Option<&Vec<std::path::PathBuf>>,
    min_size: u64,
    chromosomes: &Vec<String>,
) -> Result<FxHashMap<String, Vec<(u64, u64)>>> {
    if let Some(paths) = beds {
        load_blacklists(paths, min_size, chromosomes)
    } else {
        Ok(FxHashMap::default())
    }
}

/// Load per-chromosome scaling factors (if provided).
///
/// Implementation details:
/// - Uses `load_scaling_factors_tsv` to parse the command-line TSV into a
///   chromosome keyed map of `(start, end, factor)` tuples.
/// - Returns an empty map when no scaling factors were supplied, avoiding
///   unnecessary allocations inside the calling code.
///
/// Parameters:
/// - `scale_args`: Normalisation argument bundle.
/// - `chromosomes`: Chromosome ordering requested by the command.
/// - `contigs`: BAM target metadata, used to validate the TSV content.
///
/// Returns:
/// - A scaling factor map ready for lookups by chromosome.
///
/// Errors:
/// - Propagates IO or format errors when the TSV cannot be read or does not
///   match the BAM contigs.
pub fn load_scaling_map(
    scale_args: &ScaleGenomeArgs,
    chromosomes: &[String],
    contigs: &Contigs,
) -> Result<FxHashMap<String, Vec<(u64, u64, f32)>>> {
    if let Some(path) = &scale_args.scaling_factors {
        load_scaling_factors_tsv(path, chromosomes, contigs).context("load scaling factors")
    } else {
        Ok(FxHashMap::with_hasher(Default::default()))
    }
}
