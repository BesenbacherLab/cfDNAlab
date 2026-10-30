# Remaining benchmark work

Implemented behavior is specified in [benchmarks_spec.md](benchmarks_spec.md) and [benchmarks_concordance_spec.md](benchmarks_concordance_spec.md). This plan retains only unfinished implementation, verification and experimental work.

The implemented production defaults are still **two warmups and five measured repetitions**. The user has now agreed **one warmup plus three measurements**, with no full Griffin warmup. FrEIA's existing broken both-end workflow remains excluded from concordance until the accepted left-only replacement is implemented. Midpoint sites will come from Griffin's public paper set through the supplied BED. Do not restore automatic site selection or regenerate the existing remapped BAM to obtain extra preparation statistics.

Coverage comparisons use the already implemented base-resolution MAE as the primary similarity measure. Signed mean, RMSE, identical fraction, base-resolution Pearson and the residual distribution are implemented secondary diagnostics. Per-base comparisons do not use genomic binning.

## Follow-up direction, implementation pending

Keep [D1-D10](new_review_plan.md#decisions-and-current-user-direction) for Monday's discussion. D1-D6 now record the user's current direction rather than the review's earlier objections.

- [ ] D1: Update normal defaults, estimates, examples, budgets and relevant fixtures to 1+3. Treat Griffin's full warmups separately.
- [ ] D2: Put Griffin in its own jobs, permit runs longer than 24 hours and omit a full Griffin warmup. Decide repetition splitting and whether cache preparation is needed. Running Griffin on its 1000-interval demo is a potential short warmup, as suggested by the user. It reads data around those intervals, so its benefit for caching the full production input is uncertain. No decision to use it has been made. Other packing and serial-tool settings remain open.
- [ ] D3: Apply the settled simple resource policy: allocate N cores and configure native thread/core/worker options for that budget, respecting additional-thread semantics. Keep no-thread-option tools in their selected core comparisons. Remove inconsistent feature-specific library overrides without adding controls for hypothetical nesting.
- [ ] D4: Retain the pinned public Griffin site-set provenance and verify its actual size. Discuss window width, normalization and any subsets before changing them.
- [ ] D5: Implement FrEIA's intact left-only workflow and comparison. Use native frequencies and/or explicitly labeled doubled-count estimates under the documented mirrored-side assumption. Retain native left counts.
- [ ] D6: Add cfDNAPro end motifs. Its inclusion and the use of clearly disclosed side mirroring for approximate similarity are accepted. Matching settings, reference identity and reporting details remain open.

## Planned regression fixtures still to prepare

The existing [fixture file](tests/test_benchmarks.py) covers native coverage formats, base weighting, residual distributions including disk batches, base-resolution Pearson, contig order, concatenated zstd frames, incomplete/invalid tracks, regional coordinate identity, scope intersection, timing summaries, Slurm joins, single-core-setting jobs without automatic tool exclusions, Picard version status, failed installer retries, reconstructed TLEN/flags and FrEIA exclusion. [Review regressions](tests/test_review_regressions.py) now prepare explicit-method/duplicate rejection, early scope errors, native signal/timeout statuses, failure diagnostics and continuation, partial/all-failed budgets, samtools tail bins, missing Picard FR and midpoint normalization/empty-profile checks. Neither file has been executed during implementation.

These cases from the original plan are not fully represented in those fixtures:

| Area | Remaining fixture and hand-derived expectation |
| --- | --- |
| Selected-tool setup | A cfDNAlab-only selection must not probe/install Picard, cfDNAPro, FinaleToolkit or Griffin. Shared samtools/setup-Python dependencies remain explicit. |
| Installer controls | Exercise the actual `--mark-outdated` option, a missing packed archive and successful resumption after interruption, preserving completed unrelated tools. |
| Native gap/D/N coverage | Forward `pos=100`, `30M`, reverse `pos=170`, `30M` spans 100 bp with a 40 bp inter-mate gap. Read depth covers 60 bp. Replacing forward CIGAR with `15M2D15M` leaves `[115,117)` unsupported in cfDNAlab fragment coverage, while mosdepth fragment mode fills the 100 bp span. |
| Length/orientation | Check ordinary inward bwa TLEN against directional spans, reverse-start-first overlap selection differences and integer truncation when samtools halves an odd passing-record count. |
| Motif range/window | A valid 20 bp pair is excluded by the normal 30 bp minimum and included by the wide diagnostic. A 1 Mb contig triggers the selected FinaleToolkit window guard. Test rank/log-ratio metrics, including qualifying zero frequencies. |
| Midpoint input | Small synthetic parser fixtures preserve `+`, `-`, `.`, repeats and overlaps, and reject invalid strands/numeric-only Griffin groups. No biological site selection or preparation helper is needed. |
| Failure rows/measurement | Native failure and continuation are prepared through the lengths runner. Add coverage/ends/midpoints integration cases and explicit placement, missing-output and failed-flush cases, checking TSV field counts, statuses, missing measurements and available diagnostics. |
| Preparation edges | Unmapped merged records, unchanged existing pairs, deterministic sampling and the new unpaired-name guard before sampling. |

Fragment fixtures must be at least 10 bp. Fixtures using normal filters must respect their default 30 bp minimum. Derive expectations before execution. The repository owner runs tests. Do not run them automatically.

## Linux and cluster verification not yet performed

- [ ] Verify complete revised Linux Conda solves, including FinaleToolkit pins and Griffin Python 3.7.12/PyYAML. Save resolved environments rather than treating package availability as a successful solve.
- [ ] Verify interrupted installation/resumption, force/outdated reinstall and conda-pack 0.9.2 relocation including its sibling progress module. Exercise the actual Picard version probe.
- [ ] Use a tiny wiring input on the compute node to verify native CLI output, linkage, Python/container isolation and default cfDNAlab QC plots/fonts. An `ldd` pass alone does not establish Fontconfig/font availability.
- [ ] Run the prepared fixtures in the appropriate environment and inspect numerical results. Static Python/Bash and Cargo checks passed during implementation, but do not establish runtime correctness.
- [ ] Establish Slurm plugin/frequency, descendant coverage and file-cache treatment once on the cluster. Verify actual `sacct` job/step joins, missing/delayed records, failures and the primary memory label.

## Experimental records and measurements still required

- [ ] Record the existing final BAM/BAI and exact FASTA/FAI/2bit identities/checksums, header, flagstat, idxstats, MarkDuplicates metrics and actual remapping versions/commands. Retain original preparation logs if available. Do not rebuild the BAM solely to fill missing counters.
- [ ] Record the supplied midpoint BED identity/checksum and its derivation from the pinned public Griffin set. Keep its agreed coordinates, groups, strands, site count and selection unchanged during benchmarking.
- [ ] Run native/proper-pair and wide-bound concordance on the actual remap. Measure filtered count effects from outputs rather than reusing the original BAM's sampled non-proper percentage. Keep biological differences/scopes explicit.
- [ ] Pilot cfDNAPro progressively from about 0.1x to choose a physical node with demonstrated headroom for all methods. Run full-input method/core pilots, splitting when necessary, and establish scratch capacity/filesystem behavior.
- [ ] Replace the illustrative [scaling scenarios](README.md#runs-and-budgets) with full-input pilot budgets using the agreed 1+3 and no full Griffin warmup, plus margin and staging/preservation/cleanup. Use separate jobs for 1, 4, 12 and 32 cores. Griffin gets its own jobs and is not subject to a 24-hour cap. Other packing and serial-tool settings remain open in D2. Use matched hardware, preferring the same physical node when available, and retain resource failures as outcomes.
- [ ] Collect production timings, completed accounting, retained profiles and measured points. Fill supplementary tables with sample sizes, mean/SD, median/quartiles, memory policy and workload notes. No final timing, capacity or concordance conclusion has been measured here.

## Optional later scope, still deferred

- deepTools fragment coverage, only with a scientific comparison justifying its selection/CIGAR differences.
- A 1x/5x/15x depth series after the primary fixed-input experiment succeeds.
- Bootstrap intervals, extra plots or drift-monitoring anchors if observed results justify them.

FrEIA's broken both-end aggregate remains excluded. The accepted replacement is the intact left-only workflow in D5, with native frequencies and/or clearly labeled doubled-count estimates. Recovery from the damaged aggregate remains excluded.
