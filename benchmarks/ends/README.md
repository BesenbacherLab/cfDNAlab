# End-motif concordance and speed

Compare runnable tools on the same prepared, coordinate-sorted hg38 BAM, retaining non-autosomes and all alignment flags. Each tool applies its own filters. Use the exact alignment reference as a `.2bit` file with matching contig names and lengths. Disable GC correction, genomic scaling, blacklists, base-quality filtering and complement collapsing.

The main comparison uses all 256 inward-oriented 5' 4-mers, MAPQ >=30 and fragment lengths 30-1000 bp inclusive. cfDNAPro supplies native reference-derived left ends. Normalize each profile by its native counted-end total and compare motifs by label. Left-end versus both-end frequency comparisons are explicitly approximate under the side-mirroring assumption. Native left counts are retained without doubling.

## Tools and matching settings

| CLI tool | Native command and settings | Closest cfDNAlab settings | Expected differences |
| --- | --- | --- | --- |
| FinaleToolkit `d16632ec` | `end-motifs BAM REF -k 4 --strand both --min-length 30 --max-length 1000 --min-mapq 30 --threads C` | Reference/aligned cfDNAlab with matching bounds | Applies inclusive length bounds and midpoint window ownership. Proper read1/TLEN selection remains native. |
| cfDNAPro 1.7.4 | Curated `readBam()` with explicit MAPQ/flag selection, user-side inclusive R width subsetting, then `callMotif(motif_type="s", motif_length=4)` | Reference/aligned cfDNAlab, normal pairs/all contigs, matching MAPQ/length bounds except cfDNAPro reference compatibility | Leftmost reference bases only. Uses the selected bundled BSgenome reference and native pair/orientation/reference-bound checks. The reference manifest is retained. Comparison against both ends assumes mirrored side frequencies. |

### FinaleToolkit differences

cfDNAlab defines paired fragment spans as `forward.pos` to `reverse.reference_end`, requires both mates to pass filtering, and keeps its normal pairing policy unless `--require-proper-pair` is requested. FinaleToolkit selects read1 and reconstructs spans using TLEN. Its reader does not check the mate's `MQ` tag. A separately labeled cfDNAlab `--require-proper-pair` run isolates the proper-pair selection difference.

cfDNAlab excludes fragments with hard clipping on either mate. Reference-derived motifs affected by terminal indels are excluded at the affected end, and outer soft clipping is limited to 256 bp per end. FinaleToolkit has no corresponding CIGAR filters. The pinned fix assigns motifs by midpoint to avoid counting a fragment in two adjacent windows. Counts can therefore differ even when proper-pair selection is matched.

The pinned FinaleToolkit commit applies both length bounds and counts each fragment in a single window. A fragment longer than twice its read length can still be lost when its midpoint lies in a different window from its read1 alignment. The job refuses contigs with lengths exactly divisible by 1 Mb, because the native window builder omits their last full window.

Sources: [cfDNAlab pairing and end handling](../../src/shared/fragment/ends_fragment.rs), [motif encoding](../../src/commands/ends/motifs.rs), [FinaleToolkit end motifs](https://github.com/epifluidlab/FinaleToolkit/blob/d16632ec5ae8999e774234457f6b0f9405c8f590/src/finaletoolkit/frag/_end_motifs.py), [alignment reader](https://github.com/epifluidlab/FinaleToolkit/blob/d16632ec5ae8999e774234457f6b0f9405c8f590/src/finaletoolkit/io/alignment.py), [window construction](https://github.com/epifluidlab/FinaleToolkit/blob/d16632ec5ae8999e774234457f6b0f9405c8f590/src/finaletoolkit/frag/_motif_common.py).

### cfDNAPro motifs

The [R entry point](cfdnapro_end_motifs.R) uses the BAM contigs present in its selected BSgenome reference, computed from both contig lists and passed through the documented `chromosome_to_keep` argument. Dropped contigs are logged and retained in `contig_scope.tsv`. For this hg38 no-alt BAM, the missing contig is chrEBV (13 MAPQ >=30 alignments reported for the input). Other tools retain whole-BAM scope. Inclusive fragment widths are filtered by user-side R subsetting inside timing because cfDNAPro has no motif length-filter option. Native type `s` extracts the leftmost K forward-reference bases. BAM reading, curation, reference loading, motif extraction and table writing remain timed. There is no native thread option.

`BENCH_GENOME=hg38` selects the image's UCSC hg38 BSgenome. The result retains the package name/version and sequence metadata in `reference.tsv`. Verify its sequence identity with the alignment reference before treating reference-source differences as resolved. cfDNAPro does not use the supplied `.2bit` directly. Its native pairing/orientation/reference-bound checks and cfDNAlab's CIGAR/end handling can still differ.

Sources: pinned [callMotif](https://github.com/hw538/cfDNAPro/blob/ab9eeb166b8614b9ec2ce4ceda190dc295c65a7c/R/callMotif.R), [sequence extraction](https://github.com/hw538/cfDNAPro/blob/ab9eeb166b8614b9ec2ce4ceda190dc295c65a7c/R/internalFunctions.R#L676-L734), [reference provenance](https://github.com/hw538/cfDNAPro/blob/ab9eeb166b8614b9ec2ce4ceda190dc295c65a7c/R/genomeReference.R).

Each successful cfDNAPro invocation saves `runs/<run_order>/cfdnapro-stage-timings.tsv`. It reports elapsed seconds and the fraction of R elapsed time for package/reference setup, native `readBam()` and curation, the benchmark's inclusive length filter, native `callMotif()`, and output writing. The clocks surround the existing operations without repeating them or forcing garbage collection. Native `readBam()` materializes the paired alignments and curated fragment object in memory, as in the authors' [quick start](https://github.com/hw538/cfDNAPro#quick-start-1). The benchmark filter then subsets that object. The complete container/R process remains the end-to-end timing result, including startup and the small stage-timing table. Stage fractions exclude R/container startup and timing-table writing, so they describe the cost within R rather than a second timing benchmark.

## Install and run

[../install_tools.sh](../install_tools.sh) installs all coverage and end-motif benchmark tools, with separate environments for each tool and for reading outputs. It builds the current cfDNAlab checkout. Tool versions and source revisions are pinned for reproducibility, as listed in the [installation reference](../README.md#installation). Installation is separate from timing. Both benchmarks use the same FinaleToolkit archive.

Use the commented [production run procedure](../README.md#runs-and-budgets) for paths, submissions at 1/4/12/32 cores, fresh result directories and accounting. It times `cfdnalab-reference`, `finaletoolkit` and `cfdnapro` with `BENCH_K=4`. A separate `cfdnalab-reference-proper` run at 32 cores checks concordance with FinaleToolkit once and is excluded from timing summaries.

Full runs automatically write `concordance.tsv` after all tool attempts finish, using the job's `BENCH_K`. This comparison is outside every tool's timing and step memory measurements. Single-method jobs retain their profile and explain that there is no pair to compare. A comparison failure fails the job, with details in `concordance-error.log`, while preserving completed timings and profiles.

The cfDNAlab run labels select the timing configuration and its proper-pair concordance check:

| Run label | Configuration |
| --- | --- |
| `cfdnalab-reference` | Reference bases, aligned boundaries, normal pairing, all BAM contigs. |
| `cfdnalab-reference-proper` | Same, adding `--require-proper-pair` as a diagnostic comparison. |

For full runs, the submission launcher reserves an exclusive node and all its memory. The runner stages input copies and environments locally and uses sequential `srun --jobid` steps with the requested binding. It records the actual hostname and logical/physical CPU IDs and rejects unexpected placement. Hard CPU confinement also depends on the cluster's Slurm cgroup configuration. Source files are not modified and no sample is downloaded.

Submit separate jobs for the other core settings, each with a new results directory. Methods run in contiguous blocks with warmups before measurements at the job's single core count. All attempted runs, including warmups, remain in `timings.tsv`. Native command startup, internal reference staging, counting and native output writing are timed. Scheduler startup, input staging, profile preservation and cleanup are outside the timer. No cache eviction is performed. One native profile per run label is retained from the first measured repetition. Other temporary outputs and logs are removed. Failed diagnostics are retained, the remaining repetitions of that method are skipped and other selected methods continue. The job exits nonzero after those methods finish. See the [shared failure and timeout behavior](../README.md#timing-and-slurm-memory).

The comparator reports Pearson/Spearman r, maximum absolute log2 ratio for motifs reaching frequency >=1e-4 in either profile (infinite when the other is zero), frequency MAE/RMSE, maximum absolute difference, total variation distance and normalized Shannon MDS differences. Native counted-end totals are included when available, with `left_end_selection` and `right_end_selection` identifying `left` or `both`. `side_comparison=approximate-side-mirroring` marks different end selections. This assumes sufficiently similar inward frequencies on the two sides. Counts are never doubled. FinaleToolkit's global files only contain frequencies, so their raw totals are unavailable. Metrics use all `4^k` motifs, including zeros. All pairwise rows are written. Interpret deliberately different scopes and configurations using the tables above. Frequencies use fractions, not percentages.

## Supplementary result templates

The tables below contain no measured results. Fill them from the saved concordance metrics and completed timing/accounting records.

| Tool comparison | k | Pearson r | MAE / RMSE | Max. diff. / TV | Expected differences |
| --- | --- | --- | --- | --- | --- |
| cfDNAlab reference/proper vs FinaleToolkit `d16632ec` | 4 | [ ] | [ ] / [ ] | [ ] / [ ] | Both-mate filtering, actual span vs TLEN, CIGAR and query boundaries. |
| cfDNAlab reference/native vs FinaleToolkit `d16632ec` | 4 | [ ] | [ ] / [ ] | [ ] / [ ] | Above differences plus proper-pair selection. |
| cfDNAlab reference/native vs cfDNAPro left | 4 | [ ] | [ ] / [ ] | [ ] / [ ] | Approximate side mirroring, reference identity and native selection. |

All rows request all contigs, MAPQ >=30 and inclusive lengths 30-1000 bp, with the exceptions noted.

| Tool / configuration | Cores | Warmups / measured | Mean process+flush, s | SD seconds | Slurm peak, KiB | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| cfDNAlab reference/native | 1 / 4 / 12 / 32, separate jobs | 1 / 3 each | [ ] | [ ] | [ ] | 4-mers, normal pairing. |
| cfDNAlab reference/proper | 32 | 0 / 1 | Not a timing result | Not a timing result | [ ] | Concordance diagnostic only, flag 0x2 required. |
| FinaleToolkit `d16632ec` | 1 / 4 / 12 / 32, separate jobs | 1 / 3 each | [ ] | [ ] | [ ] | 4-mers, native inclusive length bounds. |
| cfDNAPro 1.7.4 | 1 / 4 / 12 / 32, separate jobs | 1 / 3 each | [ ] | [ ] | [ ] | Reference-derived left motifs, full BAM-to-profile process. No native thread argument. |

Expand timing rows to a core budget per row for the final A4 supplementary table. Keep long notes below the table if needed. GNU time's maximum RSS is not the simultaneous total memory of a multiprocessing tree. Retain the native CLI versions and dependency environment archives with the benchmark record.

Join completed Slurm step memory after the job. Timing rows also record filesystem flush time, process-plus-flush time and output bytes. See the [shared timing TSV documentation](../README.md) for units and measurement limits. Submit from the repository root so helper scripts can be staged before warmups.

FinaleToolkit is installed from upstream merge commit `d16632ec5ae8999e774234457f6b0f9405c8f590` ([PR #194](https://github.com/epifluidlab/FinaleToolkit/pull/194)). Timing rows include that commit. Reinstall only FinaleToolkit with `bash benchmarks/install_tools.sh ROOT --only finaletoolkit` before these runs.
