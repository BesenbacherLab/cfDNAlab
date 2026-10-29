# End-motif concordance and speed

Compare runnable tools on the same prepared, coordinate-sorted hg38 BAM, retaining non-autosomes and all alignment flags. Each tool applies its own filters. Use the exact alignment reference as a `.2bit` file with matching contig names and lengths. Disable GC correction, genomic scaling, blacklists, base-quality filtering and complement collapsing. No sample analyses or timings have been run while preparing these scripts.

The main comparison uses all 256 inward-oriented 5' 4-mers, MAPQ >=30 and fragment lengths 30-1000 bp inclusive. FrEIA needs a separate comparison of all 64 3-mers. Normalize each profile by its total counted ends, then compare motifs by their labels.

## Tools and matching settings

| CLI tool | Native command and settings | Closest cfDNAlab settings | Expected differences |
| --- | --- | --- | --- |
| FinaleToolkit 1.1.0 | `end-motifs BAM REF -k 4 --strand both --min-length 30 --max-length 1000 --min-mapq 30 --threads C` | `ends --k-inside 4 --k-outside 0 --source-inside reference --clip-strategy aligned --all-motifs --chromosomes all`, matching requested MAPQ/length bounds | Reference bases at both ends, proper pairs and read1/TLEN-based selection. **This release does not apply the requested end-motif length bounds.** Its window queries can also count boundary-crossing read1 alignments twice. Keep these limitations visible. |
| FrEIA, commit `2c32fb9` | Native BAM filter, `1_extract_fragment_ends.py`, then `2_data_transformation.py --both -fsmin 30 -fsmax 1001 -subs 1 --bootstrap_sample 0 --subchroff` | `ends --k-inside 3 --k-outside 0 --source-inside read --clip-strategy aligned --all-motifs --chromosomes all`, matching MAPQ/length bounds | Native aggregation produces 3-mers, adds coordinate deduplication and drops pairs with N in either extracted 10-base sequence. Right-end labels need complementing before pooling the two ends. These differences remain in concordance results. |
| cfDNAFE, commit `90c6979` | `cfDNAFE bam` | No equivalent configuration for the inspected counting code | Apparent end-motif counting error, detailed below. Excluded from equivalent-output timing pending a correction. |

The original Lo-group FRAGMA paper's [Methods and availability section](https://pmc.ncbi.nlm.nih.gov/articles/PMC9636966/) did not identify a runnable original software release. FRAGMA is therefore excluded. The [epifluidlab repository](https://github.com/epifluidlab/FRAGMA) calls itself a reimplementation for FinaleMe review and is not substituted for original software.

### FinaleToolkit differences

cfDNAlab defines paired fragment spans as `forward.pos` to `reverse.reference_end`, requires both mates to pass filtering, and keeps its normal pairing policy unless `--require-proper-pair` is requested. FinaleToolkit selects read1 and reconstructs spans using TLEN. Its reader does not check the mate's `MQ` tag. A separately labeled cfDNAlab `--require-proper-pair` run isolates the proper-pair selection difference.

cfDNAlab excludes fragments with hard clipping on either mate. Reference-derived motifs affected by terminal indels are excluded at the affected end, and outer soft clipping is limited to 256 bp per end. FinaleToolkit has no corresponding CIGAR filters or unique fragment ownership check across its 1 Mb queries. These differences prevent a general promise of identical counts even when proper-pair selection is matched.

The ignored length bounds are documented from source inspection. Benchmark only the current release, without patches or an older-version comparison.

Sources: [cfDNAlab pairing and end handling](../../src/shared/fragment/ends_fragment.rs), [motif encoding](../../src/commands/ends/motifs.rs), [FinaleToolkit end motifs](https://github.com/epifluidlab/FinaleToolkit/blob/v1.1.0/src/finaletoolkit/frag/_end_motifs.py), [alignment reader](https://github.com/epifluidlab/FinaleToolkit/blob/v1.1.0/src/finaletoolkit/io/alignment.py), [window construction](https://github.com/epifluidlab/FinaleToolkit/blob/v1.1.0/src/finaletoolkit/frag/_motif_common.py).

### FrEIA differences and workload

The runner follows the authors' [BAM filtering and extraction rule](https://github.com/mouliere-lab/FrEIA/blob/2c32fb9d9fab774f4ae57931df7fd7d75654d343/workflow/rules/2_FrEIA/1_preprocessing.smk), skipping remapping. Its [extraction script](https://github.com/mouliere-lab/FrEIA/blob/2c32fb9d9fab774f4ae57931df7fd7d75654d343/workflow/scripts/1_FrEIA/1_extract_fragment_ends.py) uses aligned read sequences and TLEN. It caches records by QNAME without cfDNAlab's orientation checks, and clears the cache when it exceeds 10,000 entries. Its preprocessing excludes unmapped, secondary, supplementary and duplicate records but has no explicit QC-fail or proper-pair requirement.

The [aggregation script](https://github.com/mouliere-lab/FrEIA/blob/2c32fb9d9fab774f4ae57931df7fd7d75654d343/workflow/scripts/1_FrEIA/2_data_transformation.py) uses an exclusive maximum length, hence `1001` for inclusive 1000 bp. `--both` is necessary because it otherwise reports only left ends. No bootstrap sampling or intentional downsampling is requested. The comparator sums native `CountLen` values, complements the already reversed right-end labels, pools both ends and supplies zeros for absent ACGT motifs. It does not remove FrEIA's deduplication or N filtering.

FrEIA forces extraction to **one worker when its filtered BAM exceeds 8000 MiB**, regardless of `-t`. The timing TSV records that actual extraction-worker count. Filtering and indexing still use the allocated core budget, and aggregation is serial. Its timed workload also writes an intermediate BAM/index, a per-fragment Parquet table and mononucleotide summaries. Label this as the complete native BAM-to-profile workflow. The per-fragment table is held in memory, so pilot it before scheduling the full repetitions.

### cfDNAFE source audit

At commit `90c6979a98ef9424a7659ac88522b2a8f6e0578f`, [the counting code](https://github.com/Cuiwanxin1998/cfDNAFE/blob/90c6979a98ef9424a7659ac88522b2a8f6e0578f/src/cfDNAFE/utils.py#L195-L237) calls `get_End_motif(End_motif, forward_end3, forward_end3)`. It supplies the right-end sequence twice, then complements one copy without reversing it. This is not the conventional left plus reverse-complemented right profile. No upstream source is patched for the benchmark.

Other differences include proper-pair selection, rejecting both mates when either contains S/H/P CIGAR operations, strict length inequalities and a default chromosome list omitting chr17. Parallelism distributes BAM files rather than splitting a single BAM. Its `bam` command also writes fragment BED, GC values, breakpoint motifs and MDS. Motif counts are accumulated before its BED blacklist subtraction.

## Install and run

[../install_tools.sh](../install_tools.sh) installs all coverage and end-motif benchmark tools, with separate environments for each tool and for reading outputs. It builds the current cfDNAlab checkout. FrEIA uses the authors' dependency pins. Only the newest tool releases are compared, pinned for reproducibility after verification on 2026-10-28. FrEIA has no published release, so its latest source commit is pinned. Installation is separate from timing. Both benchmarks use the same FinaleToolkit archive.

```bash
bash benchmarks/install_tools.sh /path/to/benchmark-tools

export BENCH_BAM=/path/to/BH01.hg38.paired.15x.bam
export BENCH_REF_2BIT=/path/to/the-alignment-reference.2bit
export BENCH_CFDNA=/path/to/benchmark-tools/envs/cfdnalab/bin/cfdna
export BENCH_FINALE_ARCHIVE=/path/to/benchmark-tools/finaletoolkit.tar.gz
export BENCH_FREIA_ARCHIVE=/path/to/benchmark-tools/freia.tar.gz
export BENCH_LOCAL_ROOT=/path/to/node-local-scratch

# 4-mer pilot: native tool outputs and cfDNAlab configurations for concordance
export BENCH_RESULTS=/path/to/end-motif-pilot
export BENCH_K=4
export BENCH_CORES="1"
export BENCH_WARMUPS=0
export BENCH_REPETITIONS=1
export BENCH_METHODS="cfdnalab-read cfdnalab-read-aligned cfdnalab-read-include cfdnalab-reference cfdnalab-reference-proper cfdnalab-reference-autosomes finaletoolkit"
sbatch --constraint=gen6 --time=12:00:00 benchmarks/ends/benchmark_end_motifs.sh

# After the pilot finishes, compare its saved native motif profiles
/path/to/benchmark-tools/envs/comparison/bin/python \
    benchmarks/ends/compare_end_motifs.py "$BENCH_RESULTS/profiles" \
    > "$BENCH_RESULTS/concordance.tsv"
```

Run a separate FrEIA pilot with `BENCH_K=3`, a new `BENCH_RESULTS` directory and `BENCH_METHODS="cfdnalab-read-aligned cfdnalab-reference finaletoolkit freia"`. Pass `--k 3` to the comparison script for those outputs. Other timing settings work the same way.

The cfDNAlab run labels select these configurations. They are not additional competing tools:

| Run label | Configuration |
| --- | --- |
| `cfdnalab-reference` | Reference bases, aligned boundaries, normal pairing, all BAM contigs. |
| `cfdnalab-reference-proper` | Same, adding `--require-proper-pair` as a diagnostic comparison. |
| `cfdnalab-read` | Read bases, `--clip-strategy skip`. Omit a soft-clipped end while retaining the other usable end. |
| `cfdnalab-read-aligned` | Read bases, `--clip-strategy aligned`. Compare against reference/aligned to quantify read/reference differences. |
| `cfdnalab-read-include` | Read bases, `--clip-strategy include-at-aligned-boundary`. Include terminal soft-clipped bases, keeping aligned spans for length filtering. |
| `cfdnalab-reference-autosomes` | Reference/aligned, selecting the BAM's actual autosome names. Its output and runtime have a smaller scope than whole-reference calls. |

For timing, choose the tool calls and cfDNAlab configuration supported by the concordance pilot. For example:

```bash
export BENCH_RESULTS=/path/to/end-motif-timing
export BENCH_K=4
export BENCH_CORES="1 4 8 12 32"
export BENCH_WARMUPS=3
export BENCH_REPETITIONS=10
export BENCH_METHODS="cfdnalab-reference cfdnalab-reference-proper finaletoolkit"
sbatch --constraint=gen6 --time=24:00:00 benchmarks/ends/benchmark_end_motifs.sh
```

Replace the generation selector and wall time with the cluster's settings. The runner stages input copies and environments locally, reserves an exclusive node and all its memory, and uses sequential `srun --jobid` steps with the requested binding. It records the actual hostname and logical/physical CPU IDs and rejects unexpected placement. Hard CPU confinement also depends on the cluster's Slurm cgroup configuration. Source files are not modified and no sample is downloaded.

Each repetition block shuffles tool/core combinations and records their actual order. All runs, including warmups, remain in `timings.tsv`. Native command startup, internal reference staging, counting and native output writing are timed. Scheduler startup, input staging, profile preservation and cleanup are outside the timer. No cache eviction is performed. One native profile per run label is retained from repetition 1 at the first requested core budget. Other temporary outputs and logs are removed. Failed tool messages go to the Slurm log, and a failed timing row is retained before stopping.

The comparator reports Pearson r, frequency MAE/RMSE, maximum absolute difference, total variation distance and normalized Shannon MDS differences. Native counted-end totals are included when available. FinaleToolkit's global files only contain frequencies, so their raw totals are unavailable. Metrics use all `4^k` motifs, including zeros. All pairwise rows are written. Interpret deliberately different scopes and configurations using the tables above. Frequencies use fractions, not percentages.

## Supplementary result placeholders

| Tool comparison | k | Pearson r | MAE / RMSE | Max. diff. / TV | Expected differences |
| --- | --- | --- | --- | --- | --- |
| cfDNAlab reference/proper vs FinaleToolkit 1.1.0 | 4 | [ ] | [ ] / [ ] | [ ] / [ ] | Ignored length bounds, both-mate filtering, actual span vs TLEN, CIGAR and query boundaries. |
| cfDNAlab reference/native vs FinaleToolkit 1.1.0 | 4 | [ ] | [ ] / [ ] | [ ] / [ ] | Above differences plus proper-pair selection. |
| cfDNAlab read/aligned vs FrEIA | 3 | [ ] | [ ] / [ ] | [ ] / [ ] | Native deduplication, N filtering, pair cache and TLEN. |

All rows request all contigs, MAPQ >=30 and inclusive lengths 30-1000 bp, with the exceptions noted. Report clipping/read-source diagnostic results separately if useful for explaining differences.

| Tool / configuration | Cores | Warmups / measured | Mean seconds | SD seconds | Max. RSS, KiB | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| cfDNAlab reference/native | 1 / 4 / 8 / 12 / 32 | 3 / 10 each | [ ] | [ ] | [ ] | 4-mers, normal pairing. |
| cfDNAlab reference/proper | 1 / 4 / 8 / 12 / 32 | 3 / 10 each | [ ] | [ ] | [ ] | 4-mers, flag 0x2 required. |
| FinaleToolkit 1.1.0 | 1 / 4 / 8 / 12 / 32 | 3 / 10 each | [ ] | [ ] | [ ] | 4-mers, length bounds ignored. |
| cfDNAlab read/aligned | 1 / 4 / 8 / 12 / 32 | 3 / 10 each | [ ] | [ ] | [ ] | Separate 3-mer series. |
| FrEIA | 1 / 4 / 8 / 12 / 32 | 3 / 10 each | [ ] | [ ] | [ ] | 3-mers, full BAM-to-profile workflow. Extraction forced serial above 8000 MiB. |

Expand timing rows to a core budget per row for the final A4 supplementary table. Keep long notes below the table if needed. GNU time's maximum RSS is not the simultaneous total memory of a multiprocessing tree. Retain the native CLI versions and dependency environment archives with the benchmark record.

Timing rows also include sampled peak and time-weighted mean process-tree RSS, with sample counts and interval. See the [shared timing TSV documentation](../README.md) for units and measurement limits. Submit from the repository root so helper scripts can be staged before warmups.
