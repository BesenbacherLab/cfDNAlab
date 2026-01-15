use crate::{
    commands::{
        cli_common::{ensure_output_dir, load_blacklist_map},
        frag_to_bam::config::FragToBamConfig,
    },
    shared::{
        blacklist::is_blacklisted, io::open_text_reader, reference::load_chrom_sizes_with_order,
        tiled_run::make_temp_dir,
    },
};
use anyhow::{Context, Result, bail};
use fxhash::{FxHashMap, FxHashSet};
use rust_htslib::bam::{
    self, Format, Header,
    header::HeaderRecord,
    record::{Cigar, CigarString, Record},
};
use std::collections::hash_map::Entry;
use std::{
    fs,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::PathBuf,
    time::Instant,
};

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

    println!();
    println!("Statistics");
    println!("----------");
    let elapsed = start_time.elapsed();
    println!("  Input lines: {}", counters.lines);
    println!("  Parsed fragments: {}", counters.parsed_fragments);
    println!(
        "  Rejected (chromosome filter): {}",
        counters.rejected_chromosome
    );
    println!("  Rejected (length): {}", counters.rejected_length);
    println!("  Rejected (mapq): {}", counters.rejected_mapq);
    println!("  Rejected (blacklist): {}", counters.rejected_blacklist);
    println!("  Written to BAM: {}", counters.written);
    println!("----------");
    println!("Output BAM: {}", output_path.display());
    println!("Elapsed time: {:.2?}", elapsed);

    Ok(())
}

fn run_inner(opt: &FragToBamConfig) -> Result<(FragToBamCounters, PathBuf)> {
    ensure_output_dir(&opt.output_dir)?;

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

    for (line_idx, line_res) in reader.lines().enumerate() {
        let line_number = line_idx as u64 + 1;
        counters.lines += 1;
        let line = line_res.with_context(|| format!("Reading line {}", line_number))?;
        if line.trim().is_empty() {
            continue;
        }
        let frag = parse_frag_line(&line, line_number)?;

        if !allowed_chromosomes.contains(&frag.chrom) {
            counters.rejected_chromosome += 1;
            continue;
        }
        match current_chr.as_ref().map(|s| s.as_str()) {
            None => {
                // First chromosome encountered
                current_chr = Some(frag.chrom.clone());
                current_chrom_len = *chrom_sizes
                    .get(&frag.chrom)
                    .expect("chromosome length available for first chromosome")
                    as u64;
                last_start = Some(frag.start);
                bl_ptr = 0;
                chroms_observed.push(frag.chrom.clone());
            }
            Some(chr_name) if chr_name == frag.chrom.as_str() => {
                // Enforce coordinate monotonicity within the chromosome
                if let Some(prev_start) = last_start {
                    if frag.start < prev_start {
                        bail!(
                            "Order error: Fragment out of order on {} at line {} (saw {}-{}, previous start {})",
                            frag.chrom,
                            line_number,
                            frag.start,
                            frag.end,
                            prev_start
                        );
                    }
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
                current_chrom_len = *chrom_sizes
                    .get(&frag.chrom)
                    .expect("chromosome length available for next chromosome")
                    as u64;
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
                opt.blacklist_strategy.clone(),
                frag.start,
                frag.end,
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
            "{}\t{}\t{}\t{}\t{}",
            frag.chrom, frag.start, frag.end, frag.mapq, frag.strand
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
        .join(format!("{}.bam", opt.output_prefix.trim()));
    let mut writer = bam::Writer::from_path(&output_path, &header, Format::Bam)
        .context("Creating BAM writer")?;

    // Second pass: write BAM in chrom_sizes order for predictable tid ordering
    let mut write_idx: u64 = 0;
    let qname_prefix = "fragment_";
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
            let frag = parse_frag_line(&line, line_number)?;
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
            .push_tag(b"VN", &"1.6")
            // Chromosome order in the input frag file defines the reference order we emit
            .push_tag(b"SO", &"coordinate"),
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
                .push_tag(b"LN", &len),
        );
        tid_lookup.insert(chr.clone(), idx as i32);
    }

    Ok((header, tid_lookup))
}

fn parse_frag_line(line: &str, line_number: u64) -> Result<ParsedFragment> {
    let mut parts = line.split('\t');
    let chrom = parts
        .next()
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| anyhow::anyhow!("Missing chromosome on line {}", line_number))?;
    let start: u64 = parts
        .next()
        .ok_or_else(|| anyhow::anyhow!("Missing start on line {}", line_number))?
        .trim()
        .parse()
        .with_context(|| format!("Invalid start coordinate on line {}", line_number))?;
    let end: u64 = parts
        .next()
        .ok_or_else(|| anyhow::anyhow!("Missing end on line {}", line_number))?
        .trim()
        .parse()
        .with_context(|| format!("Invalid end coordinate on line {}", line_number))?;
    let mapq: u8 = parts
        .next()
        .ok_or_else(|| anyhow::anyhow!("Missing mapq on line {}", line_number))?
        .trim()
        .parse()
        .with_context(|| format!("Invalid mapq on line {}", line_number))?;
    let strand = parts
        .next()
        .ok_or_else(|| anyhow::anyhow!("Missing strand on line {}", line_number))?
        .trim()
        .as_bytes()
        .get(0)
        .copied()
        .map(|b| b as char)
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
    })
}

fn make_record(frag: &ParsedFragment, tid: i32, prefix: &str, idx: u64) -> Result<Record> {
    let frag_len: u32 = (frag
        .end
        .checked_sub(frag.start)
        .ok_or_else(|| anyhow::anyhow!("Negative fragment length"))?)
        as u32;

    let qname = format!("{}:{}", prefix, idx);
    let seq = vec![b'N'; frag_len as usize];
    let qual = vec![40u8; frag_len as usize];
    let cigar = CigarString(vec![Cigar::Match(frag_len)]);

    let mut record = Record::new();
    record.set_tid(tid);
    record.set_pos(frag.start as i64);
    record.set_insert_size(0);
    record.set_mapq(frag.mapq);
    // Flag only the reverse strand; these are single-end records with no mate information
    let flags = if frag.strand == '-' { 0x10 } else { 0 };
    record.set_flags(flags);
    record.set(qname.as_bytes(), Some(&cigar), &seq, &qual);
    Ok(record)
}
