# Remaining benchmark work

Implemented behavior is specified in [benchmarks_spec.md](../../../../benchmarks/benchmarks_spec.md) and [benchmarks_concordance_spec.md](../../../../benchmarks/benchmarks_concordance_spec.md). This plan retains only unfinished implementation, verification and experimental work.

Preserve the supplied production sites and the existing remapped BAM. Do not restore automatic production-site selection or regenerate that BAM to obtain extra preparation statistics.

## Remaining experimental choices

Unresolved parts of D2/D4/D6 and D7-D10 remain in [new_review_plan.md](new_review_plan.md#decisions-and-current-user-direction). Do not treat those choices as settled by the implementation.

## Planned regression fixtures still to prepare

The existing [fixture file](../../../../benchmarks/tests/test_benchmarks.py) covers native coverage formats, base weighting, residual distributions including disk batches, base-resolution Pearson, contig order, concatenated zstd frames, incomplete/invalid tracks, regional coordinate identity, scope intersection, timing summaries, Slurm joins, single-core-setting jobs without automatic tool exclusions, Picard version status, failed installer retries, reconstructed TLEN/flags and native left-motif comparison. [Review regressions](../../../../benchmarks/tests/test_review_regressions.py) prepare explicit-method/duplicate rejection, early scope errors, native signal/timeout statuses, failure diagnostics and continuation, partial/all-failed budgets, Griffin demo isolation from measured runs, samtools tail bins, missing Picard FR and midpoint normalization/empty-profile checks. Neither file has been executed during implementation.

These cases from the original plan are not fully represented in those fixtures:

| Area | Remaining fixture and manually derived expectation |
| --- | --- |
| Selected-tool setup | A cfDNAlab-only selection must not probe/install Picard, cfDNAPro, FinaleToolkit or Griffin. Shared samtools/setup-Python dependencies remain explicit. |
| Installer controls | Exercise the actual `--mark-outdated` option, a missing packed archive and successful resumption after interruption, preserving completed unrelated tools. |
| Native gap/D/N coverage | Forward `pos=100`, `30M`, reverse `pos=170`, `30M` spans 100 bp with a 40 bp inter-mate gap. Read depth covers 60 bp. Replacing forward CIGAR with `15M2D15M` leaves `[115,117)` unsupported in cfDNAlab fragment coverage, while mosdepth fragment mode fills the 100 bp span. |
| Length/orientation | Check ordinary inward bwa TLEN against directional spans, reverse-start-first overlap selection differences and integer truncation when samtools halves an odd passing-record count. |
| Motif range/window | A valid 20 bp pair is excluded by the normal 30 bp minimum. A 1 Mb contig triggers the selected FinaleToolkit window guard. Test rank/log-ratio metrics, including qualifying zero frequencies. |
| Midpoint input | Small synthetic parser fixtures preserve `+`, `-`, `.`, repeats and overlaps, and reject invalid strands/numeric-only Griffin groups. No biological site selection or preparation helper is needed. |
| Failure rows/measurement | Native failure and continuation are prepared through the lengths runner. Add coverage/ends/midpoints integration cases and explicit placement, missing-output and failed-flush cases, checking TSV field counts, statuses, missing measurements and available diagnostics. |
| Preparation edges | Unmapped merged records, unchanged existing pairs, deterministic sampling and the new unpaired-name guard before sampling. |

Fragment fixtures must be at least 10 bp. Fixtures using normal filters must respect their default 30 bp minimum. Derive expectations before execution. The repository owner runs tests. Do not run them automatically.

## Linux and cluster verification not yet performed

- [ ] Verify complete revised Linux Conda solves, including FinaleToolkit pins and Griffin Python 3.7.12/PyYAML. Save resolved environments rather than treating package availability as a successful solve.
- [ ] Verify interrupted installation/resumption, force/outdated reinstall and conda-pack 0.9.2 relocation including its sibling progress module. Exercise the actual Picard version probe.
- [ ] Use a small input on the compute node to verify native CLI output, linkage, Python/container isolation, Griffin's demo warmup and default cfDNAlab QC plots/fonts. An `ldd` pass alone does not establish Fontconfig/font availability.
- [ ] Run the prepared fixtures in the appropriate environment and inspect numerical results. Static Python/Bash and Cargo checks passed during implementation, but do not establish runtime correctness.
- [ ] Establish Slurm plugin/frequency, descendant coverage and file-cache treatment once on the cluster. Verify actual `sacct` job/step joins, missing/delayed records, failures and the primary memory label.

## Experimental records and measurements still required

- [ ] Record the existing final BAM/BAI and exact FASTA/FAI/2bit identities/checksums, header, flagstat, idxstats, MarkDuplicates metrics and actual remapping versions/commands. Retain original preparation logs if available. Do not rebuild the BAM solely to fill missing counters.
- [ ] Record the supplied midpoint BED identity/checksum and its derivation from the pinned public Griffin set. Keep its agreed coordinates, groups, strands, site count and selection unchanged during benchmarking.
- [ ] Run the proper-pair concordance diagnostics once at 32 cores on the actual remap, outside the normal timing matrix. Measure filtered count effects from outputs rather than reusing the original BAM's sampled non-proper percentage. Keep biological differences/scopes explicit.
- [ ] Record runtime, memory and scratch use from the selected production jobs. Preserve resource failures as outcomes. The [production procedure](../../../../benchmarks/README.md#runs-and-budgets) uses no pilot jobs or illustrative scaling estimates.
- [ ] Use separate timing jobs for 1, 4, 12 and 32 cores, with 1+3 normally and three measurements plus the short demo for Griffin. Budget setup, preservation, concordance and cleanup within each allocation. Use matched hardware and retain node identity with each result.
- [ ] Collect production timings, completed accounting, retained profiles and measured points. Fill supplementary tables with sample sizes, mean/SD, median/quartiles, memory policy and workload notes. No final timing, capacity or concordance conclusion has been measured here.

## Optional later scope, still deferred

- deepTools fragment coverage, only with a scientific comparison justifying its selection/CIGAR differences.
- A 1x/5x/15x depth series after the primary fixed-input experiment succeeds.
- Bootstrap intervals, extra plots or drift-monitoring anchors if observed results justify them.
