use crate::{
    commands::{
        cli_common::{ensure_output_dir, load_blacklist_map},
        frag_to_bam::config::FragToBamConfig,
    },
    shared::{
        blacklist::is_blacklisted,
        cli_output,
        interval::Interval,
        io::{dot_join, open_text_reader},
        reference::load_chrom_sizes_with_order,
        tiled_run::make_temp_dir,
    },
};
use anyhow::{Context, Result, anyhow, bail};
use fxhash::{FxHashMap, FxHashSet};
use rust_htslib::bam::{
    self, Format, Header,
    header::HeaderRecord,
    record::{Aux, Cigar, CigarString, Record},
};
use std::collections::hash_map::Entry;
use std::{
    fs,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};
use tracing::warn;

const COMMAND_TARGET: &str = "frag-to-bam";

#[derive(Debug, Default)]
struct FragToBamCounters {
    lines: u64,
    parsed_fragments: u64,
    rejected_chromosome: u64,
    rejected_length: u64,
    rejected_mapq: u64,
    rejected_blacklist: u64,
    written: u64,
}

#[derive(Debug)]
struct ParsedFragment {
    chrom: String,
    start: u64,
    end: u64,
    mapq: u8,
    strand: char,
    gc_weight: Option<f32>,
    coverage_scaling_weight: Option<f32>,
    count_scaling_weight: Option<f32>,
    flen: Option<u32>,
}

#[derive(Debug, Clone)]
struct FragColumnIndices {
    chromosome: usize,
    start: usize,
    end: usize,
    mapq: usize,
    strand: usize,
    gc_weight: Option<usize>,
    coverage_scaling_weight: Option<usize>,
    count_scaling_weight: Option<usize>,
    flen: Option<usize>,
}

#[derive(Debug, Clone)]
struct FragColumnLayout {
    indices: FragColumnIndices,
    skip_first_non_empty_line: bool,
}

/// Execute the frag-to-bam conversion.
///
/// Parameters:
/// - `opt`: Fully resolved configuration for the `frag-to-bam` command.
///
/// Returns:
/// - `Ok(())` when the BAM is written successfully.
///
/// Errors:
/// - Propagates IO and parsing errors when reading inputs or writing results, aborting the run on
///   the first failure.
pub fn run(opt: &FragToBamConfig) -> Result<()> {
    let start_time = Instant::now();
    let (counters, output_path) = run_inner(opt)?;

    let elapsed = start_time.elapsed();
    cli_output::write_primary_line("");
    cli_output::write_primary_line("Statistics");
    cli_output::write_primary_line("----------");
    cli_output::write_primary_line(&format!("  Input lines: {}", counters.lines));
    cli_output::write_primary_line(&format!(
        "  Parsed fragments: {}",
        counters.parsed_fragments
    ));
    cli_output::write_primary_line(&format!(
        "  Rejected (chromosome filter): {}",
        counters.rejected_chromosome
    ));
    cli_output::write_primary_line(&format!(
        "  Rejected (length): {}",
        counters.rejected_length
    ));
    cli_output::write_primary_line(&format!("  Rejected (mapq): {}", counters.rejected_mapq));
    cli_output::write_primary_line(&format!(
        "  Rejected (blacklist): {}",
        counters.rejected_blacklist
    ));
    cli_output::write_primary_line(&format!("  Written to BAM: {}", counters.written));
    cli_output::write_primary_line("----------");
    cli_output::write_primary_line(&format!("Output BAM: {}", output_path.display()));
    cli_output::write_primary_line(&format!("Elapsed time: {:.2?}", elapsed));

    Ok(())
}

fn run_inner(opt: &FragToBamConfig) -> Result<(FragToBamCounters, PathBuf)> {
    ensure_output_dir(&opt.output_dir)?;
    let column_layout = resolve_frag_column_layout(opt)?;

    let (chrom_sizes_order, chrom_sizes) = load_chrom_sizes_with_order(&opt.chrom_sizes)
        .context("Loading chromosome sizes for BAM header")?;

    let chromosomes = {
        let want_all = opt
            .chromosomes
            .chromosomes
            .as_ref()
            .map(|chrs| chrs.len() == 1 && chrs[0].eq_ignore_ascii_case("all"))
            .unwrap_or(false);

        if want_all {
            chrom_sizes_order.clone()
        } else {
            opt.chromosomes.resolve_chromosomes(None)?
        }
    };

    if chromosomes.is_empty() {
        bail!("No chromosomes configured to read");
    }

    for chr in &chromosomes {
        if !chrom_sizes.contains_key(chr) {
            bail!("Chromosome '{}' missing from chrom sizes file", chr);
        }
    }

    // Chromosome membership, used to ensure inputs only contain expected chromosomes
    let allowed_chromosomes: FxHashSet<String> = chromosomes.iter().cloned().collect();

    let blacklist_map = load_blacklist_map(
        opt.blacklist.as_ref(),
        opt.blacklist_min_size,
        0,
        &chromosomes,
    )
    .context("Loading blacklist intervals")?;

    let temp_dir = make_temp_dir(&opt.output_dir, opt.output_prefix.trim())
        .context("Creating temp directory for frag-to-bam")?;

    let reader = open_text_reader(&opt.frag)
        .with_context(|| format!("Opening fragment file {}", opt.frag.display()))?;

    let mut counters = FragToBamCounters::default();
    let mut current_chr: Option<String> = None;
    let mut finished_chromosomes: FxHashSet<String> = FxHashSet::default();
    let mut last_start: Option<u64> = None;
    let mut current_chrom_len: u64 = 0;

    // Pointer into the current chromosome's blacklist intervals for streaming overlap checks
    let mut bl_ptr: usize = 0;
    let mut chroms_observed: Vec<String> = Vec::new();
    let mut temp_paths: FxHashMap<String, PathBuf> = FxHashMap::default();
    let mut temp_writers: FxHashMap<String, BufWriter<File>> = FxHashMap::default();

    /* First pass - validate and filter fragments */

    let mut non_empty_lines_seen = 0_u64;
    for (line_idx, line_res) in reader.lines().enumerate() {
        let line_number = line_idx as u64 + 1;
        counters.lines += 1;
        let line = line_res.with_context(|| format!("Reading line {}", line_number))?;
        if line.trim().is_empty() {
            continue;
        }
        non_empty_lines_seen += 1;
        if column_layout.skip_first_non_empty_line && non_empty_lines_seen == 1 {
            continue;
        }
        let frag = parse_frag_line(&line, line_number, &column_layout.indices)?;

        if !allowed_chromosomes.contains(&frag.chrom) {
            counters.rejected_chromosome += 1;
            continue;
        }
        match current_chr.as_deref() {
            None => {
                // First chromosome encountered
                current_chr = Some(frag.chrom.clone());
                current_chrom_len = *chrom_sizes.get(&frag.chrom).ok_or_else(|| {
                    anyhow!(
                        "Chromosome '{}' from the fragment file was not found in --chrom-sizes",
                        frag.chrom
                    )
                })? as u64;
                last_start = Some(frag.start);
                bl_ptr = 0;
                chroms_observed.push(frag.chrom.clone());
            }
            Some(chr_name) if chr_name == frag.chrom.as_str() => {
                // Enforce coordinate monotonicity within the chromosome
                if let Some(prev_start) = last_start
                    && frag.start < prev_start
                {
                    bail!(
                        "Order error: Fragment out of order on {} at line {} (saw {}-{}, previous start {})",
                        frag.chrom,
                        line_number,
                        frag.start,
                        frag.end,
                        prev_start
                    );
                }
                last_start = Some(frag.start);
            }
            Some(chr_name) => {
                // New chromosome encountered. Previous chromosome is finished
                finished_chromosomes.insert(chr_name.to_string());
                if finished_chromosomes.contains(&frag.chrom) {
                    bail!(
                        "Order error: Chromosome '{}' appears after moving past it (line {})",
                        frag.chrom,
                        line_number
                    );
                }
                current_chr = Some(frag.chrom.clone());
                current_chrom_len = *chrom_sizes.get(&frag.chrom).ok_or_else(|| {
                    anyhow!(
                        "Chromosome '{}' from the fragment file was not found in --chrom-sizes",
                        frag.chrom
                    )
                })? as u64;
                last_start = Some(frag.start);
                bl_ptr = 0;
                chroms_observed.push(frag.chrom.clone());
            }
        }

        debug_assert!(current_chrom_len > 0, "chromosome length not set");
        let chrom_len = current_chrom_len;
        if frag.end > chrom_len {
            bail!(
                "Fragment exceeds chromosome bounds on {} ({}-{}, chrom len {}) at line {}",
                frag.chrom,
                frag.start,
                frag.end,
                chrom_len,
                line_number
            );
        }
        let frag_len = frag.end - frag.start;
        counters.parsed_fragments += 1;
        if frag_len < opt.fragment_lengths.min_fragment_length as u64
            || frag_len > opt.fragment_lengths.max_fragment_length as u64
        {
            counters.rejected_length += 1;
            continue;
        }
        if frag.mapq < opt.min_mapq {
            counters.rejected_mapq += 1;
            continue;
        }

        // Check overlap with blacklisted regions
        let chrom_blacklist = blacklist_map
            .get(&frag.chrom)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);
        let is_in_blacklist = !chrom_blacklist.is_empty()
            && is_blacklisted(
                chrom_blacklist,
                opt.blacklist_strategy,
                Interval::new(frag.start, frag.end)?,
                opt.fragment_lengths.max_fragment_length as u64,
                &mut bl_ptr,
            );
        if is_in_blacklist {
            counters.rejected_blacklist += 1;
            continue;
        }

        let writer = match temp_writers.entry(frag.chrom.clone()) {
            Entry::Occupied(o) => o.into_mut(),
            Entry::Vacant(v) => {
                let path = temp_dir.join(format!("{}.frag.tmp", frag.chrom));
                temp_paths.insert(frag.chrom.clone(), path.clone());
                let file = File::create(&path)
                    .with_context(|| format!("Creating temp frag file {}", path.display()))?;
                v.insert(BufWriter::with_capacity(1 << 20, file))
            }
        };
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            frag.chrom,
            frag.start,
            frag.end,
            frag.mapq,
            frag.strand,
            format_optional_f32(frag.gc_weight),
            format_optional_f32(frag.coverage_scaling_weight),
            format_optional_f32(frag.count_scaling_weight),
            format_optional_u32(frag.flen),
        )
        .with_context(|| format!("Writing temp fragment for line {}", line_number))?;
        counters.written += 1;
    }

    // Flush temp writers to ensure contents are readable in second pass
    for writer in temp_writers.values_mut() {
        writer.flush().context("Flushing temp fragment files")?;
    }
    drop(temp_writers);

    /* Second pass (from temps) - Write to BAM */

    if chroms_observed.is_empty() {
        fs::remove_dir_all(&temp_dir).context("Cleaning up temp directory")?;
        bail!("No fragments passed filters; no BAM to write");
    }

    // Build header using the chrom_sizes file order so tids match the reference order
    let header_chroms: Vec<String> = chrom_sizes_order.clone();
    let (header, tid_lookup) = build_header(&header_chroms, &chrom_sizes)?;
    let output_path = opt
        .output_dir
        .join(dot_join(&[opt.output_prefix.trim(), "fragments.bam"]));
    let mut writer = bam::Writer::from_path(&output_path, &header, Format::Bam)
        .context("Creating BAM writer")?;

    // Second pass: write BAM in chrom_sizes order for predictable tid ordering
    let mut write_idx: u64 = 0;
    let qname_prefix = "fragment";
    for chr in &chrom_sizes_order {
        let path = match temp_paths.get(chr) {
            Some(p) => p,
            None => continue,
        };
        let file =
            File::open(path).with_context(|| format!("Opening temp fragment file for {}", chr))?;
        let reader = BufReader::with_capacity(1 << 20, file);
        for (line_idx, line_res) in reader.lines().enumerate() {
            let line_number = line_idx as u64 + 1;
            let line = line_res.with_context(|| {
                format!("Reading temp fragment line {} for {}", line_number, chr)
            })?;
            let frag = parse_temp_fragment_line(&line, line_number)?;
            let tid = *tid_lookup
                .get(&frag.chrom)
                .expect("tid lookup constructed for all chromosomes");
            write_idx = write_idx.saturating_add(1);
            let record = make_record(&frag, tid, qname_prefix, write_idx)?;
            writer
                .write(&record)
                .with_context(|| format!("Writing BAM record for {}", chr))?;
        }
    }

    fs::remove_dir_all(&temp_dir).context("Cleaning up temp directory")?;

    Ok((counters, output_path))
}

fn build_header(
    chromosomes: &[String],
    chrom_sizes: &FxHashMap<String, u32>,
) -> Result<(Header, FxHashMap<String, i32>)> {
    let mut header = Header::new();
    header.push_record(
        HeaderRecord::new(b"HD")
            // SAM format version per hts-specs (currently 1.6)
            .push_tag(b"VN", "1.6")
            // Chromosome order in the input frag file defines the reference order we emit
            .push_tag(b"SO", "coordinate"),
    );

    let mut tid_lookup: FxHashMap<String, i32> =
        FxHashMap::with_capacity_and_hasher(chromosomes.len(), Default::default());
    for (idx, chr) in chromosomes.iter().enumerate() {
        let len = *chrom_sizes
            .get(chr)
            .ok_or_else(|| anyhow::anyhow!("Chromosome '{}' missing length", chr))?;
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", chr.as_str())
                .push_tag(b"LN", len),
        );
        tid_lookup.insert(chr.clone(), idx as i32);
    }

    Ok((header, tid_lookup))
}

fn parse_frag_line(
    line: &str,
    line_number: u64,
    indices: &FragColumnIndices,
) -> Result<ParsedFragment> {
    let columns: Vec<&str> = line.split('\t').collect();

    let chrom =
        get_required_column(&columns, indices.chromosome, "chromosome", line_number)?.to_string();
    let start: u64 = get_required_column(&columns, indices.start, "start", line_number)?
        .parse()
        .with_context(|| format!("Invalid start coordinate on line {}", line_number))?;
    let end: u64 = get_required_column(&columns, indices.end, "end", line_number)?
        .parse()
        .with_context(|| format!("Invalid end coordinate on line {}", line_number))?;
    let mapq: u8 = get_required_column(&columns, indices.mapq, "mapq", line_number)?
        .parse()
        .with_context(|| format!("Invalid mapq on line {}", line_number))?;
    let strand = get_required_column(&columns, indices.strand, "strand", line_number)?
        .as_bytes()
        .first()
        .copied()
        .map(|base| base as char)
        .ok_or_else(|| anyhow::anyhow!("Missing strand on line {}", line_number))?;
    if strand != '+' && strand != '-' {
        bail!("Strand must be '+' or '-' on line {}", line_number);
    }
    if end <= start {
        bail!(
            "Fragment end must be greater than its start on line {} ({}-{})",
            line_number,
            start,
            end
        );
    }
    Ok(ParsedFragment {
        chrom: chrom.to_string(),
        start,
        end,
        mapq,
        strand,
        gc_weight: parse_optional_f32_column(
            &columns,
            indices.gc_weight,
            "gc_weight",
            line_number,
        )?,
        coverage_scaling_weight: parse_optional_f32_column(
            &columns,
            indices.coverage_scaling_weight,
            "coverage_scaling_weight",
            line_number,
        )?,
        count_scaling_weight: parse_optional_f32_column(
            &columns,
            indices.count_scaling_weight,
            "count_scaling_weight",
            line_number,
        )?,
        flen: parse_optional_u32_column(&columns, indices.flen, "flen", line_number)?,
    })
}

fn make_record(frag: &ParsedFragment, tid: i32, prefix: &str, idx: u64) -> Result<Record> {
    let frag_len: u32 = (frag
        .end
        .checked_sub(frag.start)
        .ok_or_else(|| anyhow::anyhow!("Negative fragment length"))?)
        as u32;

    let qname = format!("{}_{}", prefix, idx);
    let seq = vec![b'N'; frag_len as usize];
    let qual = vec![40u8; frag_len as usize];
    let cigar = CigarString(vec![Cigar::Match(frag_len)]);

    let mut record = Record::new();
    record.set_tid(tid);
    record.set_pos(frag.start as i64);
    record.set_insert_size(0);
    record.set_mapq(frag.mapq);
    // Flag only the reverse strand
    // These are unpaired records with no mate information
    let flags = if frag.strand == '-' { 0x10 } else { 0 };
    record.set_flags(flags);
    record.set(qname.as_bytes(), Some(&cigar), &seq, &qual);

    if let Some(gc_weight) = frag.gc_weight {
        record
            .push_aux(b"GC", Aux::Float(gc_weight))
            .with_context(|| {
                format!(
                    "Failed writing GC aux tag for fragment {}:{}-{}",
                    frag.chrom, frag.start, frag.end
                )
            })?;
    }
    if let Some(coverage_scaling_weight) = frag.coverage_scaling_weight {
        record
            .push_aux(b"COV", Aux::Float(coverage_scaling_weight))
            .with_context(|| {
                format!(
                    "Failed writing COV aux tag for fragment {}:{}-{}",
                    frag.chrom, frag.start, frag.end
                )
            })?;
    }
    if let Some(count_scaling_weight) = frag.count_scaling_weight {
        record
            .push_aux(b"CNT", Aux::Float(count_scaling_weight))
            .with_context(|| {
                format!(
                    "Failed writing CNT aux tag for fragment {}:{}-{}",
                    frag.chrom, frag.start, frag.end
                )
            })?;
    }
    if let Some(fragment_length_tag) = frag.flen {
        record
            .push_aux(b"FLEN", Aux::U32(fragment_length_tag))
            .with_context(|| {
                format!(
                    "Failed writing FLEN aux tag for fragment {}:{}-{}",
                    frag.chrom, frag.start, frag.end
                )
            })?;
    }

    Ok(record)
}

fn parse_temp_fragment_line(line: &str, line_number: u64) -> Result<ParsedFragment> {
    let columns: Vec<&str> = line.split('\t').collect();
    if columns.len() != 9 {
        bail!(
            "Invalid temporary fragment row at line {}. Expected 9 columns, got {}",
            line_number,
            columns.len()
        );
    }
    let indices = FragColumnIndices {
        chromosome: 0,
        start: 1,
        end: 2,
        mapq: 3,
        strand: 4,
        gc_weight: Some(5),
        coverage_scaling_weight: Some(6),
        count_scaling_weight: Some(7),
        flen: Some(8),
    };
    parse_frag_line(line, line_number, &indices)
}

fn get_required_column<'a>(
    columns: &'a [&str],
    index: usize,
    column_name: &str,
    line_number: u64,
) -> Result<&'a str> {
    columns
        .get(index)
        .map(|value| value.trim())
        .filter(|value| !value.is_empty())
        .ok_or_else(|| {
            anyhow::anyhow!(
                "Missing {} value on line {} (expected column index {})",
                column_name,
                line_number,
                index
            )
        })
}

fn parse_optional_f32_column(
    columns: &[&str],
    index: Option<usize>,
    column_name: &str,
    line_number: u64,
) -> Result<Option<f32>> {
    let Some(column_index) = index else {
        return Ok(None);
    };
    let value = columns
        .get(column_index)
        .map(|value| value.trim())
        .unwrap_or("");
    if value.is_empty() || value == "." || value.eq_ignore_ascii_case("na") {
        return Ok(None);
    }
    let parsed: f32 = value
        .parse()
        .with_context(|| format!("Invalid {} value on line {}", column_name, line_number))?;
    if !parsed.is_finite() {
        bail!("{} must be finite on line {}", column_name, line_number);
    }
    Ok(Some(parsed))
}

fn parse_optional_u32_column(
    columns: &[&str],
    index: Option<usize>,
    column_name: &str,
    line_number: u64,
) -> Result<Option<u32>> {
    let Some(column_index) = index else {
        return Ok(None);
    };
    let value = columns
        .get(column_index)
        .map(|value| value.trim())
        .unwrap_or("");
    if value.is_empty() || value == "." || value.eq_ignore_ascii_case("na") {
        return Ok(None);
    }
    let parsed: u32 = value
        .parse()
        .with_context(|| format!("Invalid {} value on line {}", column_name, line_number))?;
    Ok(Some(parsed))
}

fn format_optional_f32(value: Option<f32>) -> String {
    value
        .map(|v| v.to_string())
        .unwrap_or_else(|| ".".to_string())
}

fn format_optional_u32(value: Option<u32>) -> String {
    value
        .map(|v| v.to_string())
        .unwrap_or_else(|| ".".to_string())
}

fn resolve_frag_column_layout(opt: &FragToBamConfig) -> Result<FragColumnLayout> {
    let first_non_empty_line = read_first_non_empty_line(&opt.frag)?;
    let inline_header_columns = first_non_empty_line
        .as_deref()
        .and_then(detect_inline_header_columns);

    let explicit_header = if let Some(path) = &opt.frag_header {
        Some((path.clone(), read_header_columns(path)?))
    } else {
        None
    };

    let companion_header = if explicit_header.is_none() {
        if let Some(path) = infer_companion_header_path(&opt.frag) {
            if path.exists() {
                Some((path.clone(), read_header_columns(&path)?))
            } else {
                None
            }
        } else {
            None
        }
    } else {
        None
    };

    if inline_header_columns.is_some() {
        if let Some((explicit_path, _)) = &explicit_header {
            bail!(
                "Conflicting headers detected: both --frag-header ({}) and an inline header row in {}. Use only one header source",
                explicit_path.display(),
                opt.frag.display()
            );
        }
        if let Some((companion_path, _)) = &companion_header {
            bail!(
                "Conflicting headers detected: both companion header file ({}) and an inline header row in {}. Use only one header source",
                companion_path.display(),
                opt.frag.display()
            );
        }
    }

    let use_inline_header =
        explicit_header.is_none() && companion_header.is_none() && inline_header_columns.is_some();
    let header_columns = explicit_header
        .map(|(_, columns)| columns)
        .or_else(|| companion_header.map(|(_, columns)| columns))
        .or(inline_header_columns);

    let indices = if let Some(columns) = header_columns {
        resolve_indices_from_header(&columns, opt.ignore_extras, opt.allow_unknown_extras)?
    } else {
        resolve_default_indices(opt.ignore_extras)
    };

    Ok(FragColumnLayout {
        indices,
        skip_first_non_empty_line: use_inline_header,
    })
}

fn read_first_non_empty_line(path: &Path) -> Result<Option<String>> {
    let reader = open_text_reader(path)
        .with_context(|| format!("Opening fragment file {}", path.display()))?;
    for line_result in reader.lines() {
        let line = line_result?;
        if !line.trim().is_empty() {
            return Ok(Some(line));
        }
    }
    Ok(None)
}

fn detect_inline_header_columns(line: &str) -> Option<Vec<String>> {
    let columns: Vec<String> = line.split('\t').map(normalize_column_name).collect();
    let has_required_names = find_column_index(&columns, &["chromosome", "chrom"]).is_some()
        && find_column_index(&columns, &["start"]).is_some()
        && find_column_index(&columns, &["end"]).is_some()
        && find_column_index(&columns, &["mapq", "min_mapq"]).is_some()
        && find_column_index(&columns, &["strand", "read1_strand"]).is_some();
    if has_required_names {
        Some(columns)
    } else {
        None
    }
}

fn read_header_columns(path: &Path) -> Result<Vec<String>> {
    let reader = open_text_reader(path)
        .with_context(|| format!("Opening header file {}", path.display()))?;
    for line_result in reader.lines() {
        let line = line_result?;
        if line.trim().is_empty() {
            continue;
        }
        let columns: Vec<String> = line.split('\t').map(normalize_column_name).collect();
        if columns.is_empty() {
            continue;
        }
        return Ok(columns);
    }
    bail!("Header file {} was empty", path.display());
}

fn normalize_column_name(name: &str) -> String {
    name.trim().to_ascii_lowercase()
}

fn resolve_indices_from_header(
    columns: &[String],
    ignore_extras: bool,
    allow_unknown_extras: bool,
) -> Result<FragColumnIndices> {
    let chromosome_index = find_column_index(columns, &["chromosome", "chrom"])
        .ok_or_else(|| anyhow::anyhow!("Could not find chromosome column in header"))?;
    let start_index = find_column_index(columns, &["start"])
        .ok_or_else(|| anyhow::anyhow!("Could not find start column in header"))?;
    let end_index = find_column_index(columns, &["end"])
        .ok_or_else(|| anyhow::anyhow!("Could not find end column in header"))?;
    let mapq_index = find_column_index(columns, &["mapq", "min_mapq"])
        .ok_or_else(|| anyhow::anyhow!("Could not find mapq column in header"))?;
    let strand_index = find_column_index(columns, &["strand", "read1_strand"])
        .ok_or_else(|| anyhow::anyhow!("Could not find strand column in header"))?;

    if !ignore_extras {
        validate_extra_column_names(columns, allow_unknown_extras)?;
    }

    let gc_weight_index = if ignore_extras {
        None
    } else {
        find_column_index(columns, &["gc_weight"])
    };
    let coverage_scaling_weight_index = if ignore_extras {
        None
    } else {
        find_column_index(columns, &["coverage_scaling_weight"])
    };
    let count_scaling_weight_index = if ignore_extras {
        None
    } else {
        find_column_index(columns, &["count_scaling_weight"])
    };
    let flen_index = if ignore_extras {
        None
    } else {
        find_column_index(columns, &["flen"])
    };

    Ok(FragColumnIndices {
        chromosome: chromosome_index,
        start: start_index,
        end: end_index,
        mapq: mapq_index,
        strand: strand_index,
        gc_weight: gc_weight_index,
        coverage_scaling_weight: coverage_scaling_weight_index,
        count_scaling_weight: count_scaling_weight_index,
        flen: flen_index,
    })
}

fn resolve_default_indices(_ignore_extras: bool) -> FragColumnIndices {
    FragColumnIndices {
        chromosome: 0,
        start: 1,
        end: 2,
        mapq: 3,
        strand: 4,
        gc_weight: None,
        coverage_scaling_weight: None,
        count_scaling_weight: None,
        flen: None,
    }
}

fn find_column_index(columns: &[String], aliases: &[&str]) -> Option<usize> {
    columns.iter().position(|column| {
        aliases
            .iter()
            .any(|alias| column == &alias.trim().to_ascii_lowercase())
    })
}

fn validate_extra_column_names(columns: &[String], allow_unknown_extras: bool) -> Result<()> {
    let unsupported_columns = collect_unsupported_extra_columns(columns);
    if unsupported_columns.is_empty() {
        return Ok(());
    }

    if allow_unknown_extras {
        warn!(
            target: COMMAND_TARGET,
            "Warning: Ignoring unsupported frag header column name(s): {}. Recognized extra columns are gc_weight, coverage_scaling_weight, count_scaling_weight, and flen",
            unsupported_columns.join(", ")
        );
        Ok(())
    } else {
        bail!(
            "Unsupported frag header column name(s): {}. Extra columns must be named exactly gc_weight, coverage_scaling_weight, count_scaling_weight, or flen. Use --ignore-extras to ignore all extra columns or --allow-unknown-extras to ignore only unknown names",
            unsupported_columns.join(", ")
        );
    }
}

fn collect_unsupported_extra_columns(columns: &[String]) -> Vec<String> {
    let mut unsupported_columns = Vec::new();
    for column_name in columns {
        let is_core_column = matches!(
            column_name.as_str(),
            "chromosome"
                | "chrom"
                | "start"
                | "end"
                | "mapq"
                | "min_mapq"
                | "strand"
                | "read1_strand"
        );
        let is_supported_extra = matches!(
            column_name.as_str(),
            "gc_weight" | "coverage_scaling_weight" | "count_scaling_weight" | "flen"
        );
        if !is_core_column && !is_supported_extra {
            unsupported_columns.push(column_name.clone());
        }
    }
    unsupported_columns.sort();
    unsupported_columns.dedup();
    unsupported_columns
}

fn infer_companion_header_path(frag_path: &Path) -> Option<PathBuf> {
    let file_name = frag_path.file_name()?.to_str()?;
    const KNOWN_SUFFIXES: [&str; 4] = [
        ".frag.tsv.gz",
        ".frag.tsv.zst",
        ".frag.tsv.bgz",
        ".frag.tsv",
    ];
    for suffix in KNOWN_SUFFIXES {
        if let Some(prefix) = file_name.strip_suffix(suffix) {
            let header_name = dot_join(&[prefix, "frag.header.tsv"]);
            return Some(frag_path.with_file_name(header_name));
        }
    }
    None
}
