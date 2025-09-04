use anyhow::{Context, Result};
use std::fs::File;
use std::{
    collections::HashMap,
    io::{BufRead, BufReader},
    path::Path,
};

/// Load windows from a BED file into a per-chromosome map.
/// 
/// Parameters
/// ----------
///  - bed: Path to BED file.
///  - chromosomes: Names of chromosomes to include in output, 
///   even when not present in the BED file.
/// 
/// Returns
/// -------
///  - Mapping of 'chromosome -> sorted window coordinates (start, end, original window index)'.
pub fn load_windows_from_bed(
    bed: &Path,
    chromosomes: &Vec<String>,
) -> Result<HashMap<String, Vec<(u64, u64, u64)>>> {
    let f = File::open(bed).context("Opening BED file with windows/intervals")?;
    let reader = BufReader::new(f);
    let mut mapping: HashMap<String, Vec<(u64, u64, u64)>> = HashMap::new();
    // Ensure all chromosomes are added
    chromosomes.iter().for_each(|chr| {
        mapping.entry(chr.to_string()).or_default();
    });
    // Original interval index for reconstructing order
    let mut win_idx = 0u64;
    for line in reader.lines() {
        let l = line?;
        if l.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = l.split_whitespace().collect();
        let chr = cols[0];
        if !chromosomes.contains(&chr.to_owned()) {
            continue;
        }
        let start: u64 = cols[1].parse().context("Parsing window start")?;
        let end: u64 = cols[2].parse().context("Parsing window end")?;
        mapping
            .entry(chr.to_string())
            .or_default()
            .push((start, end, win_idx));
        win_idx += 1;
    }
    for v in mapping.values_mut() {
        // Ensure sorted windows
        v.sort_unstable_by_key(|&(s, e, _)| (s, e));
    }
    Ok(mapping)
}
