use crate::command_run::{CommandRunResult, RunOptions};
use crate::commands::cli_common::{
    ChromosomeArgs, FragmentPositionSelectionArgs, IOCArgs, Ref2BitRequiredArgs, WindowsArgs,
};
use crate::commands::fragment_kmers::config::FragmentKmersConfig;
use crate::commands::visualize_positions::config::VisualizePositionsConfig;
use crate::commands::visualize_positions::model::{LengthVisualization, Style, VizConfig};
use crate::commands::visualize_positions::select::ReadClamp;
use crate::commands::visualize_positions::{render_ascii, render_svg};
use crate::shared::interval::Interval;
use crate::shared::positioning::{BasesFrom, PositionGroup, ReferenceFrame};
use crate::shared::visualization::{AxisBounds, Track};
use anyhow::{Context, Result, anyhow, bail, ensure};
use ndarray::Array3;
use ndarray_npy::read_npy;
use rust_htslib::bam::{self, Writer, header::HeaderRecord, record::Cigar, record::CigarString};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::{self, File};
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use tempfile::Builder as TempDirBuilder;
use twobit::convert::{fasta::FastaReader, to_2bit};

const CHROM_NAME: &str = "chr1";
const FRAGMENT_GAP: u32 = 20;

/// Result from `visualize-positions`.
///
/// The command writes or prints a visualization for selected fragment positions. The result records
/// the output path when one exists.
#[derive(Debug)]
pub struct VisualizePositionsRunResult {
    /// Empty counter placeholder for the shared command result interface.
    pub counters: (),
    /// Output path when the visualization is written to a file.
    pub output_path: Option<PathBuf>,
    /// Final output files produced by the command.
    pub output_files: Vec<PathBuf>,
}

impl CommandRunResult for VisualizePositionsRunResult {
    type Counters = ();

    fn counters(&self) -> &Self::Counters {
        &self.counters
    }

    fn output_files(&self) -> &[PathBuf] {
        &self.output_files
    }

    fn primary_output(&self) -> Option<&Path> {
        self.output_path.as_deref()
    }
}

/// Run the `visualize-positions` command.
///
/// This command builds a small reference and fragment fixture, computes selected positional k-mer
/// tracks, and renders them as ASCII or SVG. It is intended for visual inspection of positional
/// selection behavior rather than high-throughput analysis.
///
/// The current implementation does not use reporting options.
///
/// Parameters
/// ----------
/// - `cfg`:
///     Fully resolved configuration for the `visualize-positions` command.
/// - `_options`:
///     Reserved reporting controls for consistency with other command runners.
///
/// Returns
/// -------
/// - `Ok(VisualizePositionsRunResult)`:
///     Output path information for the completed run.
///
/// Errors
/// ------
/// Returns an error when the fixture cannot be created, k-mer counts cannot be computed, or the
/// visualization cannot be written.
pub fn run_visualize_positions(
    cfg: &VisualizePositionsConfig,
    options: RunOptions,
) -> Result<VisualizePositionsRunResult> {
    let viz_cfg = cfg.build()?;

    if options.log_equivalent_cli {
        let command = crate::ToCliCommand::to_cli_string(cfg)?;
        let message = crate::command_run::equivalent_cli_log_message(&command);
        tracing::info!(target: "visualize-positions", "{message}");
    }

    fs::create_dir_all(&cfg.work_dir)
        .with_context(|| format!("creating work directory {}", cfg.work_dir.display()))?;

    let results = compute_visualizations(
        &viz_cfg,
        &cfg.position_selection,
        &cfg.work_dir,
        options.log_equivalent_cli,
    )?;

    let rendered = match viz_cfg.style {
        Style::Ascii => render_ascii(&results, &viz_cfg),
        Style::Svg => render_svg(&results, &viz_cfg),
    };

    if let Some(path) = &viz_cfg.output {
        fs::write(path, rendered)
            .with_context(|| format!("writing visualization to {}", path.display()))?;
    } else {
        let mut stdout = io::stdout().lock();
        stdout
            .write_all(rendered.as_bytes())
            .context("writing visualization to stdout")?;
    }

    let output_path = viz_cfg.output.clone();
    let output_files = output_path.iter().cloned().collect();
    Ok(VisualizePositionsRunResult {
        counters: (),
        output_path,
        output_files,
    })
}

fn compute_visualizations(
    viz_cfg: &VizConfig,
    position_args: &FragmentPositionSelectionArgs,
    work_dir: &Path,
    log_equivalent_cli: bool,
) -> Result<Vec<LengthVisualization>> {
    let temp_dir = TempDirBuilder::new()
        .prefix("cfdna_viz")
        .tempdir_in(work_dir)
        .context("creating visualize-positions temp dir")?;

    let synthetic = synthesize_inputs(viz_cfg, temp_dir.path())?;

    let run_k_sizes: Vec<u8> = match viz_cfg.kmer_sizes.as_ref() {
        Some(sizes) if !sizes.is_empty() => sizes.clone(),
        _ => vec![1],
    };

    let prefix = "viz_counts";
    run_fragment_kmers(
        viz_cfg,
        position_args,
        &synthetic,
        temp_dir.path(),
        &run_k_sizes,
        prefix,
        log_equivalent_cli,
    )?;

    let main_positional_spec = viz_cfg.position_specs[0].clone();

    let counts = collect_counts(temp_dir.path(), prefix, &run_k_sizes, &synthetic.windows)?;

    let clamp_mode = match viz_cfg.bases {
        BasesFrom::NearestRead => ReadClamp::Nearest,
        BasesFrom::Reads => ReadClamp::Both,
        _ => ReadClamp::None,
    };

    let mut results = Vec::with_capacity(synthetic.windows.len());
    for (window, window_counts) in synthetic.windows.iter().zip(counts.into_iter()) {
        let mut viz = build_tracks_from_counts(
            main_positional_spec.frame,
            window.len(),
            clamp_mode,
            &window_counts.offsets,
            &window_counts.coverage,
        );

        if let Some(kmer_sizes) = viz_cfg.kmer_sizes.as_ref()
            && !kmer_sizes.is_empty()
        {
            let overlays = build_overlays_from_counts(
                main_positional_spec.frame,
                window.len(),
                &viz.tracks,
                kmer_sizes,
                &window_counts.offsets_by_k,
            );
            viz.tracks.extend(overlays);
        }

        results.push(viz);
    }

    Ok(results)
}

struct FragmentWindow {
    interval: Interval<u32>,
}

impl FragmentWindow {
    #[inline]
    fn len(&self) -> u32 {
        self.interval.len()
    }

    #[inline]
    fn start(&self) -> u32 {
        self.interval.start()
    }

    #[inline]
    fn end(&self) -> u32 {
        self.interval.end()
    }
}

struct SyntheticInputs {
    bam: PathBuf,
    reference: PathBuf,
    bed: PathBuf,
    windows: Vec<FragmentWindow>,
}

fn synthesize_inputs(viz_cfg: &VizConfig, temp_dir: &Path) -> Result<SyntheticInputs> {
    if viz_cfg.fragment_lengths.is_empty() {
        bail!("no fragment lengths provided");
    }

    let mut fragments: Vec<FragmentSpec> = Vec::with_capacity(viz_cfg.fragment_lengths.len());
    let mut windows: Vec<FragmentWindow> = Vec::with_capacity(viz_cfg.fragment_lengths.len());

    let mut cursor: u32 = FRAGMENT_GAP;
    for &length in &viz_cfg.fragment_lengths {
        let start = cursor;
        let end = start
            .checked_add(length)
            .ok_or_else(|| anyhow!("fragment length overflow"))?;
        let interval = Interval::new(start, end)?;
        fragments.push(FragmentSpec { interval });
        windows.push(FragmentWindow { interval });
        cursor = end
            .checked_add(FRAGMENT_GAP)
            .ok_or_else(|| anyhow!("fragment length overflow"))?;
    }
    let reference_length = cursor as usize;

    let reference_sequence = build_reference_sequence(reference_length);
    let reference = write_reference(temp_dir, &reference_sequence)?;
    let bam = write_bam(temp_dir, &reference_sequence, &fragments)?;
    let bed = write_windows_bed(temp_dir, &windows)?;

    Ok(SyntheticInputs {
        bam,
        reference,
        bed,
        windows,
    })
}

struct FragmentSpec {
    interval: Interval<u32>,
}

impl FragmentSpec {
    #[inline]
    fn start(&self) -> u32 {
        self.interval.start()
    }

    #[inline]
    fn len(&self) -> u32 {
        self.interval.len()
    }
}

fn write_bam(temp_dir: &Path, reference: &[u8], fragments: &[FragmentSpec]) -> Result<PathBuf> {
    let bam_path = temp_dir.join("synthetic.bam");

    let mut header = bam::Header::new();
    header.push_record(
        HeaderRecord::new(b"HD")
            .push_tag(b"VN", "1.6")
            .push_tag(b"SO", "coordinate"),
    );
    header.push_record(
        HeaderRecord::new(b"SQ")
            .push_tag(b"SN", CHROM_NAME)
            .push_tag(b"LN", reference.len() as u32),
    );

    let mut writer = Writer::from_path(&bam_path, &header, bam::Format::Bam)
        .with_context(|| format!("creating {}", bam_path.display()))?;

    let mut records = Vec::with_capacity(fragments.len() * 2);
    for (idx, fragment) in fragments.iter().enumerate() {
        let fragment_length = fragment.len();
        let forward_len = (fragment_length / 2).max(1);
        let mut reverse_len = fragment_length.saturating_sub(forward_len);
        if reverse_len == 0 {
            reverse_len = 1;
        }

        let forward_start = fragment.start();
        let reverse_start = fragment
            .start()
            .checked_add(fragment_length)
            .and_then(|end| end.checked_sub(reverse_len))
            .ok_or_else(|| anyhow!("invalid fragment coordinates"))?;

        let fragment_size = (fragment_length as i64).max(1);
        let qname = format!("frag{}_{}", idx, fragment.start());

        let forward = build_read(
            &qname,
            forward_start,
            forward_len,
            reverse_start,
            fragment_size,
            false,
            true,
            reference,
        )?;
        let reverse = build_read(
            &qname,
            reverse_start,
            reverse_len,
            forward_start,
            -fragment_size,
            true,
            false,
            reference,
        )?;

        records.push(forward);
        records.push(reverse);
    }

    records.sort_by_key(|rec| (rec.tid(), rec.pos()));
    for record in &records {
        writer
            .write(record)
            .context("writing synthetic fragment record")?;
    }

    drop(writer);

    bam::index::build(&bam_path, None, bam::index::Type::Bai, 1)
        .with_context(|| format!("indexing {}", bam_path.display()))?;

    let bai_candidate = bam_path.with_extension("bam.bai");
    let bai_path = bam_path.with_extension("bai");
    let final_bai = if bai_candidate.exists() {
        fs::rename(&bai_candidate, &bai_path).with_context(|| {
            format!(
                "renaming {} -> {}",
                bai_candidate.display(),
                bai_path.display()
            )
        })?;
        bai_path
    } else {
        bai_path
    };

    if !final_bai.exists() {
        bail!("failed to locate BAM index for {}", bam_path.display());
    }

    Ok(bam_path)
}

fn build_read(
    qname: &str,
    start: u32,
    read_len: u32,
    mate_start: u32,
    insert_size: i64,
    is_reverse: bool,
    mate_is_reverse: bool,
    reference: &[u8],
) -> Result<bam::Record> {
    use rust_htslib::bam::record::Record;

    const FLAG_PAIRED: u16 = 0x1;
    const FLAG_PROPER_PAIR: u16 = 0x2;
    const FLAG_MATE_REVERSE: u16 = 0x20;
    const FLAG_REVERSE: u16 = 0x10;
    const FLAG_FIRST_MATE: u16 = 0x40;
    const FLAG_SECOND_MATE: u16 = 0x80;

    let mut record = Record::new();
    record.set_tid(0);
    record.set_pos(start as i64);
    record.set_mtid(0);
    record.set_mpos(mate_start as i64);
    record.set_insert_size(insert_size);
    record.set_mapq(60);

    let mut flags = FLAG_PAIRED | FLAG_PROPER_PAIR;
    if mate_is_reverse {
        flags |= FLAG_MATE_REVERSE;
    }
    if is_reverse {
        flags |= FLAG_SECOND_MATE | FLAG_REVERSE;
    } else {
        flags |= FLAG_FIRST_MATE;
    }
    record.set_flags(flags);

    let end = start
        .checked_add(read_len)
        .ok_or_else(|| anyhow!("read length overflow"))?;
    let end = end.min(reference.len() as u32);
    let seq_slice = &reference[start as usize..end as usize];
    let sequence: Vec<u8> = if is_reverse {
        seq_slice.iter().rev().map(complement_base).collect()
    } else {
        seq_slice.to_vec()
    };
    let qualities = vec![40u8; sequence.len()];
    let cigar = CigarString(vec![Cigar::Match(sequence.len() as u32)]);

    record.set(qname.as_bytes(), Some(&cigar), &sequence, &qualities);
    Ok(record)
}

fn complement_base(base: &u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => b'N',
    }
}

fn write_reference(temp_dir: &Path, sequence: &[u8]) -> Result<PathBuf> {
    let fasta_path = temp_dir.join("reference.fa");
    let mut fasta =
        File::create(&fasta_path).with_context(|| format!("creating {}", fasta_path.display()))?;
    writeln!(fasta, ">{}", CHROM_NAME)?;
    for chunk in sequence.chunks(60) {
        fasta.write_all(chunk)?;
        fasta.write_all(b"\n")?;
    }
    drop(fasta);

    let reference_path = temp_dir.join("reference.2bit");
    let reader = FastaReader::open(&fasta_path)
        .with_context(|| format!("opening {}", fasta_path.display()))?;
    let mut file = File::create(&reference_path)
        .with_context(|| format!("creating {}", reference_path.display()))?;
    to_2bit(&mut file, &reader)
        .map_err(|err| anyhow::anyhow!("converting FASTA to 2bit: {err}"))?;
    Ok(reference_path)
}

fn write_windows_bed(temp_dir: &Path, windows: &[FragmentWindow]) -> Result<PathBuf> {
    let bed_path = temp_dir.join("windows.bed");
    let mut file =
        File::create(&bed_path).with_context(|| format!("creating {}", bed_path.display()))?;
    for (index, window) in windows.iter().enumerate() {
        writeln!(
            file,
            "{}\t{}\t{}\twin{}",
            CHROM_NAME,
            window.start(),
            window.end(),
            index
        )?;
    }
    Ok(bed_path)
}

fn build_reference_sequence(length: usize) -> Vec<u8> {
    const PATTERN: &[u8] = b"ACGT";
    (0..length)
        .map(|idx| PATTERN[idx % PATTERN.len()])
        .collect()
}

fn run_fragment_kmers(
    viz_cfg: &VizConfig,
    position_args: &FragmentPositionSelectionArgs,
    inputs: &SyntheticInputs,
    temp_dir: &Path,
    kmer_sizes: &[u8],
    prefix: &str,
    log_equivalent_cli: bool,
) -> Result<()> {
    let mut cfg = FragmentKmersConfig::new(
        IOCArgs {
            bam: inputs.bam.clone(),
            output_dir: temp_dir.to_path_buf(),
            n_threads: 1,
        },
        Ref2BitRequiredArgs {
            ref_2bit: inputs.reference.clone(),
        },
        ChromosomeArgs {
            chromosomes: Some(vec![CHROM_NAME.to_string()]),
            chromosomes_file: None,
        },
    );

    cfg.set_output_prefix(prefix.to_string());
    cfg.set_kmer_sizes(kmer_sizes.to_vec());
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_ignore_gap(false);
    cfg.set_canonical(false);
    cfg.set_positional_counts(true);
    cfg.set_save_sparse(false);
    cfg.set_position_selection(position_args.clone());

    if let Some(&min_len) = viz_cfg.fragment_lengths.iter().min() {
        cfg.fragment_lengths_mut().min_fragment_length = min_len;
    }
    if let Some(&max_len) = viz_cfg.fragment_lengths.iter().max() {
        cfg.fragment_lengths_mut().max_fragment_length = max_len;
    }

    let mut windows_args = WindowsArgs::default();
    windows_args.by_bed = Some(inputs.bed.clone());
    cfg.set_windows(windows_args);
    crate::commands::fragment_kmers::fragment_kmers::run_fragment_kmers(
        &cfg,
        RunOptions {
            log_equivalent_cli,
            ..RunOptions::new_quiet()
        },
    )
    .context("running fragment-kmers for visualize-positions")?;
    Ok(())
}

struct WindowCounts {
    offsets: HashMap<PositionGroup, BTreeSet<i32>>,
    offsets_by_k: BTreeMap<u8, HashMap<PositionGroup, BTreeSet<i32>>>,
    coverage: HashMap<PositionGroup, BTreeSet<i32>>,
}

impl WindowCounts {
    fn new() -> Self {
        Self {
            offsets: HashMap::new(),
            offsets_by_k: BTreeMap::new(),
            coverage: HashMap::new(),
        }
    }
}

fn collect_counts(
    temp_dir: &Path,
    prefix: &str,
    kmer_sizes: &[u8],
    windows: &[FragmentWindow],
) -> Result<Vec<WindowCounts>> {
    ensure!(!kmer_sizes.is_empty(), "k-mer size list must be non-empty");
    let base_k = kmer_sizes[0];
    let mut results: Vec<WindowCounts> = windows.iter().map(|_| WindowCounts::new()).collect();

    let groups = [
        (PositionGroup::Left, "left"),
        (PositionGroup::Right, "right"),
        (PositionGroup::Mid, "mid"),
    ];

    for (group, label) in groups {
        let positions_path = temp_dir.join(crate::shared::io::dot_join(&[
            prefix,
            &format!("{label}_positions.txt"),
        ]));
        if !positions_path.exists() {
            continue;
        }

        let raw_positions = fs::read_to_string(&positions_path)
            .with_context(|| format!("reading {}", positions_path.display()))?;
        let offsets: Vec<i32> = raw_positions
            .lines()
            .filter_map(|line| line.trim().parse::<i32>().ok())
            .collect();

        if offsets.is_empty() {
            continue;
        }

        for window_counts in &mut results {
            window_counts.offsets.entry(group).or_default();
            for &k in kmer_sizes {
                window_counts
                    .offsets_by_k
                    .entry(k)
                    .or_insert_with(HashMap::new)
                    .entry(group)
                    .or_default();
            }
        }

        for &k in kmer_sizes {
            let counts_path = temp_dir.join(crate::shared::io::dot_join(&[
                prefix,
                &format!("k{k}_{label}_counts.npy"),
            ]));
            if !counts_path.exists() {
                continue;
            }

            let counts: Array3<f64> = read_npy(&counts_path)
                .with_context(|| format!("reading {}", counts_path.display()))?;
            let shape = counts.shape();
            if shape.len() != 3 || shape[0] != windows.len() || shape[1] != offsets.len() {
                bail!(
                    "unexpected shape {:?} in {} (expected {} windows × {} positions × motifs)",
                    shape,
                    counts_path.display(),
                    windows.len(),
                    offsets.len()
                );
            }

            for window_idx in 0..windows.len() {
                for (pos_idx, offset) in offsets.iter().enumerate() {
                    let mut has_signal = false;
                    for motif_idx in 0..shape[2] {
                        if counts[[window_idx, pos_idx, motif_idx]] > 0.0 {
                            has_signal = true;
                            if k == base_k {
                                let cov_set =
                                    results[window_idx].coverage.entry(group).or_default();
                                let span = coverage_span(
                                    windows[window_idx].len(),
                                    *offset,
                                    k as i32,
                                    group,
                                );
                                for idx in span {
                                    cov_set.insert(idx);
                                }
                            }
                            break;
                        }
                    }
                    if has_signal {
                        results[window_idx]
                            .offsets
                            .entry(group)
                            .or_default()
                            .insert(*offset);
                        results[window_idx]
                            .offsets_by_k
                            .entry(k)
                            .or_default()
                            .entry(group)
                            .or_default()
                            .insert(*offset);
                    }
                }
            }
        }
    }

    Ok(results)
}

fn build_tracks_from_counts(
    frame: ReferenceFrame,
    length: u32,
    clamp: ReadClamp,
    offsets: &HashMap<PositionGroup, BTreeSet<i32>>,
    coverage: &HashMap<PositionGroup, BTreeSet<i32>>,
) -> LengthVisualization {
    let mut tracks = match frame {
        ReferenceFrame::Left => {
            let indices = coverage
                .get(&PositionGroup::Left)
                .map(|set| set.iter().copied().collect())
                .unwrap_or_else(|| {
                    map_linear_positions(
                        length,
                        offsets.get(&PositionGroup::Left),
                        PositionGroup::Left,
                    )
                });
            vec![Track {
                name: "left".to_string(),
                axis: AxisBounds::new(1, length as i32),
                selected_indices: indices,
            }]
        }
        ReferenceFrame::Right => {
            let indices = coverage
                .get(&PositionGroup::Right)
                .map(|set| set.iter().copied().collect())
                .unwrap_or_else(|| {
                    map_linear_positions(
                        length,
                        offsets.get(&PositionGroup::Right),
                        PositionGroup::Right,
                    )
                });
            vec![Track {
                name: "right".to_string(),
                axis: AxisBounds::new(1, length as i32),
                selected_indices: indices,
            }]
        }
        ReferenceFrame::PerEnd => {
            let left_indices = coverage
                .get(&PositionGroup::Left)
                .map(|set| set.iter().copied().collect())
                .unwrap_or_else(|| {
                    map_linear_positions(
                        length,
                        offsets.get(&PositionGroup::Left),
                        PositionGroup::Left,
                    )
                });
            let right_indices = coverage
                .get(&PositionGroup::Right)
                .map(|set| set.iter().copied().collect())
                .unwrap_or_else(|| {
                    map_linear_positions(
                        length,
                        offsets.get(&PositionGroup::Right),
                        PositionGroup::Right,
                    )
                });
            vec![
                Track {
                    name: "left".to_string(),
                    axis: AxisBounds::new(1, length as i32),
                    selected_indices: left_indices,
                },
                Track {
                    name: "right".to_string(),
                    axis: AxisBounds::new(1, length as i32),
                    selected_indices: right_indices,
                },
            ]
        }
        ReferenceFrame::Nearest => {
            let left_positions = coverage
                .get(&PositionGroup::Left)
                .map(|set| set.iter().copied().collect::<Vec<_>>())
                .unwrap_or_else(|| {
                    map_linear_positions(
                        length,
                        offsets.get(&PositionGroup::Left),
                        PositionGroup::Left,
                    )
                });
            let right_positions = coverage
                .get(&PositionGroup::Right)
                .map(|set| set.iter().copied().collect::<Vec<_>>())
                .unwrap_or_else(|| {
                    map_linear_positions(
                        length,
                        offsets.get(&PositionGroup::Right),
                        PositionGroup::Right,
                    )
                });
            let mut fragment_positions = left_positions.clone();
            fragment_positions.extend(right_positions.iter().copied());
            fragment_positions.sort_unstable();
            fragment_positions.dedup();
            let distances = fold_fragment_positions(length, &fragment_positions);
            let mut tracks = vec![
                Track {
                    name: "fragment".to_string(),
                    axis: AxisBounds::new(1, length as i32),
                    selected_indices: fragment_positions,
                },
                Track {
                    name: "nearest".to_string(),
                    axis: AxisBounds::new(1, (length / 2).max(1) as i32),
                    selected_indices: distances,
                },
            ];
            if !left_positions.is_empty() {
                let half = length.div_ceil(2) as i32;
                if let Some(&idx) = left_positions.iter().find(|&&idx| idx > half) {
                    panic!(
                        "Nearest left track received index {} beyond half {}. \
This indicates fragment-kmers emitted starts past the nearest-read boundary.",
                        idx, half
                    );
                }
                tracks.push(Track {
                    name: "left".to_string(),
                    axis: AxisBounds::new(1, length as i32),
                    selected_indices: left_positions,
                });
            }
            if !right_positions.is_empty() {
                let half = length.div_ceil(2) as i32;
                let right_start = (length as i32 + 1) - half;
                if let Some(&idx) = right_positions.iter().find(|&&idx| idx < right_start) {
                    panic!(
                        "Nearest right track received index {} below start {}. \
This indicates fragment-kmers emitted starts past the nearest-read boundary.",
                        idx, right_start
                    );
                }
                tracks.push(Track {
                    name: "right".to_string(),
                    axis: AxisBounds::new(1, length as i32),
                    selected_indices: right_positions,
                });
            }
            tracks
        }
        ReferenceFrame::Mid => {
            let center = (length as i64) / 2;
            let mut indices: Vec<i32> = offsets
                .get(&PositionGroup::Mid)
                .map(|set| {
                    set.iter()
                        .filter_map(|offset| {
                            let relative = (*offset as i64) - center;
                            (relative >= i64::from(i32::MIN) && relative <= i64::from(i32::MAX))
                                .then_some(relative as i32)
                        })
                        .collect()
                })
                .unwrap_or_default();
            indices.sort_unstable();
            vec![Track {
                name: "mid".to_string(),
                axis: mid_axis_bounds(length),
                selected_indices: indices,
            }]
        }
    };

    apply_read_clamp_local(&mut tracks, frame, length, clamp);

    LengthVisualization {
        fragment_length: length,
        tracks,
    }
}

fn build_overlays_from_counts(
    frame: ReferenceFrame,
    length: u32,
    base_tracks: &[Track],
    overlay_k_sizes: &[u8],
    offsets_by_k: &BTreeMap<u8, HashMap<PositionGroup, BTreeSet<i32>>>,
) -> Vec<Track> {
    if overlay_k_sizes.is_empty() {
        return Vec::new();
    }

    match frame {
        ReferenceFrame::Left => {
            let base = match base_tracks.iter().find(|track| track.name == "left") {
                Some(track) => track,
                None => return Vec::new(),
            };
            let mut overlays = Vec::new();
            for &k in overlay_k_sizes {
                let Some(group_map) = offsets_by_k.get(&k) else {
                    continue;
                };
                let Some(offsets) = group_map.get(&PositionGroup::Left) else {
                    continue;
                };
                if offsets.is_empty() {
                    continue;
                }
                let indices: Vec<i32> = offsets.iter().map(|offset| offset + 1).collect();
                let mut overlay = base.clone();
                overlay.name = format!("{} k-mer starts (k={})", base.name, k);
                overlay.selected_indices = indices;
                clamp_overlay_axis(&mut overlay, length, k);
                overlays.push(overlay);
            }
            overlays
        }
        ReferenceFrame::Right => {
            let base = match base_tracks.iter().find(|track| track.name == "right") {
                Some(track) => track,
                None => return Vec::new(),
            };
            let mut overlays = Vec::new();
            for &k in overlay_k_sizes {
                let Some(group_map) = offsets_by_k.get(&k) else {
                    continue;
                };
                let Some(offsets) = group_map.get(&PositionGroup::Right) else {
                    continue;
                };
                if offsets.is_empty() {
                    continue;
                }
                let mut indices: Vec<i32> = offsets
                    .iter()
                    .map(|offset| (length as i32) - *offset)
                    .collect();
                indices.sort_unstable();
                let mut overlay = base.clone();
                overlay.name = format!("{} k-mer starts (k={})", base.name, k);
                overlay.selected_indices = indices;
                clamp_overlay_axis(&mut overlay, length, k);
                overlays.push(overlay);
            }
            overlays
        }
        ReferenceFrame::PerEnd => {
            let left_base = match base_tracks.iter().find(|track| track.name == "left") {
                Some(track) => track,
                None => return Vec::new(),
            };
            let right_base = match base_tracks.iter().find(|track| track.name == "right") {
                Some(track) => track,
                None => return Vec::new(),
            };
            let mut overlays = Vec::new();
            for &k in overlay_k_sizes {
                let Some(group_map) = offsets_by_k.get(&k) else {
                    continue;
                };

                if let Some(offsets) = group_map.get(&PositionGroup::Left)
                    && !offsets.is_empty()
                {
                    let mut indices: Vec<i32> = offsets.iter().map(|offset| offset + 1).collect();
                    indices.sort_unstable();
                    let mut overlay = left_base.clone();
                    overlay.name = format!("{} k-mer starts (k={})", left_base.name, k);
                    overlay.selected_indices = indices;
                    clamp_overlay_axis(&mut overlay, length, k);
                    overlays.push(overlay);
                }

                if let Some(offsets) = group_map.get(&PositionGroup::Right)
                    && !offsets.is_empty()
                {
                    let mut indices: Vec<i32> = offsets
                        .iter()
                        .map(|offset| (length as i32) - *offset)
                        .collect();
                    indices.sort_unstable();
                    let mut overlay = right_base.clone();
                    overlay.name = format!("{} k-mer starts (k={})", right_base.name, k);
                    overlay.selected_indices = indices;
                    clamp_overlay_axis(&mut overlay, length, k);
                    overlays.push(overlay);
                }
            }
            overlays
        }
        ReferenceFrame::Nearest => {
            let fragment_base = match base_tracks.iter().find(|track| track.name == "fragment") {
                Some(track) => track,
                None => return Vec::new(),
            };
            let nearest_base = match base_tracks.iter().find(|track| track.name == "nearest") {
                Some(track) => track,
                None => return Vec::new(),
            };
            let left_base = base_tracks.iter().find(|track| track.name == "left");
            let right_base = base_tracks.iter().find(|track| track.name == "right");
            let mut overlays = Vec::new();
            for &k in overlay_k_sizes {
                let Some(group_map) = offsets_by_k.get(&k) else {
                    continue;
                };
                let left_positions = map_linear_positions(
                    length,
                    group_map.get(&PositionGroup::Left),
                    PositionGroup::Left,
                );
                let right_positions = map_linear_positions(
                    length,
                    group_map.get(&PositionGroup::Right),
                    PositionGroup::Right,
                );
                if left_positions.is_empty() && right_positions.is_empty() {
                    continue;
                }
                let mut fragment_positions = left_positions.clone();
                fragment_positions.extend(right_positions.iter().copied());
                fragment_positions.sort_unstable();
                fragment_positions.dedup();

                let mut fragment_overlay = fragment_base.clone();
                fragment_overlay.name = format!("{} k-mer starts (k={})", fragment_base.name, k);
                fragment_overlay.selected_indices = fragment_positions.clone();
                overlays.push(fragment_overlay);

                let mut nearest_overlay = nearest_base.clone();
                nearest_overlay.name = format!("{} k-mer starts (k={})", nearest_base.name, k);
                nearest_overlay.selected_indices =
                    fold_fragment_positions(length, &fragment_positions);
                overlays.push(nearest_overlay);

                if let Some(base) = left_base
                    && !left_positions.is_empty()
                {
                    let mut overlay = base.clone();
                    overlay.name = format!("{} k-mer starts (k={})", base.name, k);
                    overlay.selected_indices = left_positions.clone();
                    overlays.push(overlay);
                }
                if let Some(base) = right_base
                    && !right_positions.is_empty()
                {
                    let mut overlay = base.clone();
                    overlay.name = format!("{} k-mer starts (k={})", base.name, k);
                    overlay.selected_indices = right_positions.clone();
                    overlays.push(overlay);
                }
            }
            overlays
        }
        ReferenceFrame::Mid => {
            let base = match base_tracks.iter().find(|track| track.name == "mid") {
                Some(track) => track,
                None => return Vec::new(),
            };
            let mut overlays = Vec::new();
            for &k in overlay_k_sizes {
                let Some(group_map) = offsets_by_k.get(&k) else {
                    continue;
                };
                let Some(offsets) = group_map.get(&PositionGroup::Mid) else {
                    continue;
                };
                if offsets.is_empty() {
                    continue;
                }
                let center = (length as i64) / 2;
                let mut indices: Vec<i32> = offsets
                    .iter()
                    .filter_map(|offset| {
                        let relative = (*offset as i64) - center;
                        (relative >= i64::from(i32::MIN) && relative <= i64::from(i32::MAX))
                            .then_some(relative as i32)
                    })
                    .collect();
                indices.sort_unstable();
                let mut overlay = base.clone();
                overlay.name = format!("{} k-mer starts (k={})", base.name, k);
                overlay.selected_indices = indices;
                overlays.push(overlay);
            }
            overlays
        }
    }
}

fn fold_fragment_positions(length: u32, starts: &[i32]) -> Vec<i32> {
    if length == 0 {
        return Vec::new();
    }
    let half = length / 2;
    let mut distances = Vec::with_capacity(starts.len());
    for &start in starts {
        if start <= 0 {
            continue;
        }
        let start_u32 = start as u32;
        let distance = if start_u32 <= half {
            start_u32
        } else {
            length - start_u32 + 1
        };
        if distance > 0 {
            distances.push(distance as i32);
        }
    }
    distances.sort_unstable();
    distances.dedup();
    distances
}

fn map_linear_positions(
    length: u32,
    offsets: Option<&BTreeSet<i32>>,
    group: PositionGroup,
) -> Vec<i32> {
    let Some(offsets) = offsets else {
        return Vec::new();
    };
    let mut values: Vec<i32> = offsets
        .iter()
        .filter_map(|offset| match group {
            PositionGroup::Left => {
                let value = offset + 1;
                (value > 0 && value <= length as i32).then_some(value)
            }
            PositionGroup::Right => {
                let value = length as i32 - offset;
                (value > 0 && value <= length as i32).then_some(value)
            }
            PositionGroup::Mid => None,
        })
        .collect();
    values.sort_unstable();
    values.dedup();
    values
}

fn coverage_span(length: u32, offset: i32, k_len: i32, group: PositionGroup) -> Vec<i32> {
    assert!(k_len > 0, "k-mer length must be positive");
    let fragment_len = length as i32;
    match group {
        PositionGroup::Left | PositionGroup::Mid => {
            let start = offset + 1;
            if start < 1 {
                panic!(
                    "{} coverage start {} fell below 1 for fragment length {}",
                    match group {
                        PositionGroup::Left => "Left",
                        PositionGroup::Mid => "Mid",
                        _ => unreachable!(),
                    },
                    start,
                    length
                );
            }
            let end = start + k_len - 1;
            if end > fragment_len {
                panic!(
                    "{} coverage index {} exceeds fragment length {}",
                    match group {
                        PositionGroup::Left => "Left",
                        PositionGroup::Mid => "Mid",
                        _ => unreachable!(),
                    },
                    end,
                    length
                );
            }
            (0..k_len).map(|delta| start + delta).collect()
        }
        PositionGroup::Right => {
            let start = fragment_len - offset;
            if start < 1 {
                panic!(
                    "Right coverage start {} fell below 1 for fragment length {}",
                    start, length
                );
            }
            let end = start - (k_len - 1);
            if end < 1 {
                panic!(
                    "Right coverage index {} fell below 1 for fragment length {}",
                    end, length
                );
            }
            (0..k_len).map(|delta| start - delta).collect()
        }
    }
}

fn apply_read_clamp_local(
    tracks: &mut [Track],
    frame: ReferenceFrame,
    length: u32,
    clamp: ReadClamp,
) {
    if matches!(clamp, ReadClamp::None) || length == 0 {
        return;
    }

    let half = length.div_ceil(2) as i32;
    let right_start = (length as i32 + 1) - half;

    for track in tracks {
        match clamp {
            ReadClamp::None => {}
            ReadClamp::Nearest => clamp_track_nearest(track, frame, half, right_start),
            ReadClamp::Both => clamp_track_both_reads(track, frame, half, right_start),
        }
    }
}

fn clamp_track_nearest(track: &mut Track, frame: ReferenceFrame, half: i32, right_start: i32) {
    match frame {
        ReferenceFrame::Nearest => {
            if track.name == "fragment" {
                if let Some(&idx) = track
                    .selected_indices
                    .iter()
                    .find(|&&idx| !(idx <= half || idx >= right_start))
                {
                    panic!(
                        "Nearest-read clamp detected nearest fragment track index {} outside <= {} or >= {}.",
                        idx, half, right_start
                    );
                }
                track
                    .selected_indices
                    .retain(|&idx| idx <= half || idx >= right_start);
            } else if track.name == "left" {
                if let Some(&idx) = track.selected_indices.iter().find(|&&idx| idx > half) {
                    panic!(
                        "Nearest-read clamp detected nearest left track index {} outside <= {}.",
                        idx, half
                    );
                }
                track.selected_indices.retain(|&idx| idx <= half);
            } else if track.name == "right" {
                if let Some(&idx) = track
                    .selected_indices
                    .iter()
                    .find(|&&idx| idx < right_start)
                {
                    panic!(
                        "Nearest-read clamp detected nearest right track index {} outside >= {}.",
                        idx, right_start
                    );
                }
                track.selected_indices.retain(|&idx| idx >= right_start);
            }
        }
        ReferenceFrame::Mid => {
            track.selected_indices.retain(|&idx| idx.abs() <= half);
        }
        ReferenceFrame::Left => {
            track.selected_indices.retain(|&idx| idx <= half);
        }
        ReferenceFrame::Right => {
            track.selected_indices.retain(|&idx| idx >= right_start);
        }
        ReferenceFrame::PerEnd => {
            if track.name == "left" {
                track.selected_indices.retain(|&idx| idx <= half);
            } else if track.name == "right" {
                track.selected_indices.retain(|&idx| idx >= right_start);
            }
        }
    }
}

fn clamp_track_both_reads(track: &mut Track, frame: ReferenceFrame, half: i32, right_start: i32) {
    match frame {
        ReferenceFrame::Nearest => {
            if track.name == "fragment" {
                if let Some(&idx) = track
                    .selected_indices
                    .iter()
                    .find(|&&idx| !(idx <= half || idx >= right_start))
                {
                    panic!(
                        "Both-read clamp detected nearest fragment track index {} outside <= {} or >= {}.",
                        idx, half, right_start
                    );
                }
                track
                    .selected_indices
                    .retain(|&idx| idx <= half || idx >= right_start);
            } else if track.name == "left" {
                if let Some(&idx) = track.selected_indices.iter().find(|&&idx| idx > half) {
                    panic!(
                        "Both-read clamp detected nearest left track index {} outside <= {}.",
                        idx, half
                    );
                }
                track.selected_indices.retain(|&idx| idx <= half);
            } else if track.name == "right" {
                if let Some(&idx) = track
                    .selected_indices
                    .iter()
                    .find(|&&idx| idx < right_start)
                {
                    panic!(
                        "Both-read clamp detected nearest right track index {} outside >= {}.",
                        idx, right_start
                    );
                }
                track.selected_indices.retain(|&idx| idx >= right_start);
            }
        }
        ReferenceFrame::Mid => {
            track.selected_indices.retain(|&idx| idx.abs() <= half);
        }
        ReferenceFrame::Left => {
            track.selected_indices.retain(|&idx| idx <= half);
        }
        ReferenceFrame::Right => {
            track.selected_indices.retain(|&idx| idx >= right_start);
        }
        ReferenceFrame::PerEnd => {
            if track.name == "left" {
                track.selected_indices.retain(|&idx| idx <= half);
            } else if track.name == "right" {
                track.selected_indices.retain(|&idx| idx >= right_start);
            }
        }
    }
}

fn mid_axis_bounds(length: u32) -> AxisBounds {
    let half = (length / 2) as i32;
    if length.is_multiple_of(2) {
        AxisBounds::new(-half, half - 1)
    } else {
        AxisBounds::new(-half, half)
    }
}

fn clamp_overlay_axis(overlay: &mut Track, length: u32, k: u8) {
    if let Some(max_start) = length
        .checked_sub(u32::from(k))
        .and_then(|value| value.checked_add(1))
    {
        let max_start_i32 = if max_start > i32::MAX as u32 {
            i32::MAX
        } else {
            max_start as i32
        };
        overlay.axis.end = overlay.axis.end.min(max_start_i32);
    }
}

#[cfg(test)]
mod coverage_tests {
    use super::*;

    fn sorted(span: Vec<i32>) -> Vec<i32> {
        let mut s = span;
        s.sort_unstable();
        s
    }

    #[test]
    fn expands_left_span_from_offset_zero() {
        assert_eq!(
            sorted(coverage_span(10, 0, 3, PositionGroup::Left)),
            vec![1, 2, 3]
        );
    }

    #[test]
    fn expands_left_span_at_tail() {
        assert_eq!(
            sorted(coverage_span(10, 7, 3, PositionGroup::Left)),
            vec![8, 9, 10]
        );
    }

    #[test]
    fn expands_right_span_near_end() {
        assert_eq!(
            sorted(coverage_span(10, 0, 3, PositionGroup::Right)),
            vec![8, 9, 10]
        );
    }

    #[test]
    fn expands_right_span_away_from_end() {
        assert_eq!(
            sorted(coverage_span(10, 5, 2, PositionGroup::Right)),
            vec![4, 5]
        );
    }

    #[test]
    #[should_panic]
    fn panics_when_right_span_underflows() {
        let _ = coverage_span(5, 4, 3, PositionGroup::Right);
    }

    #[test]
    #[should_panic]
    fn panics_when_left_span_overflows() {
        let _ = coverage_span(5, 4, 3, PositionGroup::Left);
    }
}
