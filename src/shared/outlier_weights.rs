use crate::shared::bam::Contigs;
use crate::shared::interval::IndexedInterval;
use crate::shared::overlaps::OverlappingWindows;
use anyhow::{Context, Result, bail, ensure};
use fxhash::{FxHashMap, FxHashSet};
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Sparse outlier intervals and their fragment keep weights for a chromosome.
///
/// `intervals` and `keep_weights` have matching indices. The interval vector deliberately uses
/// the same representation consumed by the project's existing linear overlap sweep.
#[derive(Debug, Clone, Default)]
pub(crate) struct ChromosomeOutlierWeights {
    pub(crate) intervals: Vec<IndexedInterval<u64>>,
    keep_weights: Vec<f64>,
}

/// Validated sparse outlier weights for the selected chromosomes.
#[derive(Debug, Clone)]
pub(crate) struct LoadedOutlierWeights {
    pub(crate) by_chromosome: FxHashMap<String, ChromosomeOutlierWeights>,
    /// Smallest strictly positive listed weight, including the implicit outside weight of `1.0`.
    pub(crate) minimum_positive_keep_weight: f64,
}

impl LoadedOutlierWeights {
    /// Construct the identity case used when no outlier-weight file was supplied.
    pub(crate) fn identity() -> Self {
        Self {
            by_chromosome: FxHashMap::with_hasher(Default::default()),
            minimum_positive_keep_weight: 1.0,
        }
    }
}

/// Load the sparse fragment keep weights written by `cfdna outliers`.
///
/// Listed intervals are retained in file order because downstream commands pass them directly to
/// the existing linear overlap sweep. Missing intervals and missing selected chromosomes are
/// valid and mean keep weight `1.0`.
pub(crate) fn load_outlier_weights_tsv(
    path: &Path,
    chromosomes: &[String],
    contigs: &Contigs,
) -> Result<LoadedOutlierWeights> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("opening outlier keep-weight TSV {}", path.display()))?;
    let mut reader = BufReader::new(file);
    let mut header = String::new();
    let mut line_number = 0usize;

    loop {
        header.clear();
        if reader.read_line(&mut header)? == 0 {
            bail!("{}: empty file, header required", path.display());
        }
        line_number += 1;

        let trimmed = header.trim();
        if trimmed.is_empty() {
            bail!(
                "{}:{}: blank lines are not allowed before the outlier keep-weight header",
                path.display(),
                line_number
            );
        }
        if trimmed.starts_with('#') {
            continue;
        }
        break;
    }

    let header_fields: Vec<&str> = header
        .trim_end_matches(&['\r', '\n'][..])
        .split('\t')
        .collect();
    let lowercase_header: Vec<String> = header_fields
        .iter()
        .map(|field| field.to_ascii_lowercase())
        .collect();
    let required_column = |name: &str| -> Result<usize> {
        lowercase_header
            .iter()
            .position(|field| field == name)
            .with_context(|| {
                format!(
                    "required column '{}' not found in outlier keep-weight header: {}",
                    name,
                    header_fields.join("\t")
                )
            })
    };
    let chromosome_column = required_column("chromosome")?;
    let start_column = required_column("start")?;
    let end_column = required_column("end")?;
    let keep_weight_column = required_column("keep_weight")?;
    let maximum_required_column = chromosome_column
        .max(start_column)
        .max(end_column)
        .max(keep_weight_column);

    let selected_chromosomes: FxHashSet<&str> = chromosomes.iter().map(String::as_str).collect();
    let mut by_chromosome: FxHashMap<String, ChromosomeOutlierWeights> =
        FxHashMap::with_hasher(Default::default());
    let mut previous_end_by_chromosome: FxHashMap<String, u64> =
        FxHashMap::with_hasher(Default::default());
    let mut minimum_positive_keep_weight = 1.0_f64;
    let mut line = String::new();

    while {
        line.clear();
        reader.read_line(&mut line)?
    } > 0
    {
        line_number += 1;
        let raw_line = line.trim_end_matches(&['\r', '\n'][..]);
        if raw_line.is_empty() || raw_line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = raw_line.split('\t').collect();
        if fields.len() <= maximum_required_column {
            bail!(
                "{}:{}: not enough columns in outlier keep-weight row, found {} and need at least {}",
                path.display(),
                line_number,
                fields.len(),
                maximum_required_column + 1
            );
        }

        let chromosome = fields[chromosome_column];
        if !selected_chromosomes.contains(chromosome) {
            continue;
        }

        let start: u64 = fields[start_column].parse().with_context(|| {
            format!(
                "{}:{}: invalid start '{}'",
                path.display(),
                line_number,
                fields[start_column]
            )
        })?;
        let end: u64 = fields[end_column].parse().with_context(|| {
            format!(
                "{}:{}: invalid end '{}'",
                path.display(),
                line_number,
                fields[end_column]
            )
        })?;
        let interval = IndexedInterval::new(start, end, line_number as u64).with_context(|| {
            format!(
                "{}:{}: invalid outlier interval [{}..{})",
                path.display(),
                line_number,
                start,
                end
            )
        })?;

        let keep_weight: f64 = fields[keep_weight_column].parse().with_context(|| {
            format!(
                "{}:{}: invalid keep_weight '{}'",
                path.display(),
                line_number,
                fields[keep_weight_column]
            )
        })?;
        ensure!(
            keep_weight.is_finite() && (0.0..=1.0).contains(&keep_weight),
            "{}:{}: keep_weight must be finite and between 0.0 and 1.0 inclusive, got {}",
            path.display(),
            line_number,
            keep_weight
        );

        let chromosome_length = contigs
            .contigs
            .get(chromosome)
            .map(|&(_, length)| length as u64)
            .with_context(|| format!("missing BAM contig information for '{chromosome}'"))?;
        ensure!(
            end <= chromosome_length,
            "{}:{}: outlier interval on '{}' ends at {}, beyond chromosome length {}",
            path.display(),
            line_number,
            chromosome,
            end,
            chromosome_length
        );

        if let Some(previous_end) = previous_end_by_chromosome.get(chromosome) {
            ensure!(
                start >= *previous_end,
                "{}:{}: outlier intervals on '{}' must be sorted and non-overlapping, previous end is {} but next start is {}",
                path.display(),
                line_number,
                chromosome,
                previous_end,
                start
            );
        }
        previous_end_by_chromosome.insert(chromosome.to_string(), end);

        let chromosome_weights = by_chromosome.entry(chromosome.to_string()).or_default();
        chromosome_weights.intervals.push(interval);
        chromosome_weights.keep_weights.push(keep_weight);
        if keep_weight > 0.0 {
            minimum_positive_keep_weight = minimum_positive_keep_weight.min(keep_weight);
        }
    }

    Ok(LoadedOutlierWeights {
        by_chromosome,
        minimum_positive_keep_weight,
    })
}

/// Resolve the minimum keep weight from overlap rows produced by the common overlap sweep.
///
/// A missing overlap has implicit keep weight `1.0`. Positive overlap fractions are not used to
/// dilute a regional weight because the correction applies to the complete fragment.
pub(crate) fn minimum_overlapping_keep_weight(
    overlaps: Option<&OverlappingWindows>,
    chromosome_weights: &ChromosomeOutlierWeights,
) -> Result<f64> {
    let Some(overlaps) = overlaps else {
        return Ok(1.0);
    };

    let mut minimum_keep_weight = 1.0_f64;
    for overlap in &overlaps.windows {
        let keep_weight = chromosome_weights
            .keep_weights
            .get(overlap.idx)
            .with_context(|| {
                format!(
                    "outlier overlap index {} has no matching keep weight among {} intervals",
                    overlap.idx,
                    chromosome_weights.keep_weights.len()
                )
            })?;
        minimum_keep_weight = minimum_keep_weight.min(*keep_weight);
    }
    Ok(minimum_keep_weight)
}

#[cfg(test)]
mod tests {
    include!("outlier_weights_tests.rs");
}
