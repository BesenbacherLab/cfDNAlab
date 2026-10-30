# Remaining benchmark work

Implemented behavior is specified in [benchmarks_spec.md](benchmarks_spec.md) and [benchmarks_concordance_spec.md](benchmarks_concordance_spec.md). This plan retains only unfinished implementation, verification and experimental work. Review history remains in [claude_review.md](claude_review.md).

Production defaults are **two warmups and five measured repetitions**. FrEIA is excluded from concordance. Midpoint sites are supplied by the user. Do not restore site generation or regenerate the existing remapped BAM to obtain extra preparation statistics.

## Unimplemented comparison output

- [ ] Add the full base-weighted signed coverage-difference distribution requested by the original plan. The current comparator gives the signed mean and below/equal/above fractions, which do not describe the full distribution of difference magnitudes. Keep memory bounded and label the evaluated scope.

## Planned regression fixtures still to prepare

The existing [fixture file](tests/test_benchmarks.py) covers native coverage formats, base weighting, contig order, concatenated zstd frames, incomplete/invalid tracks, regional coordinate identity, scope intersection, timing summaries, Slurm joins, serial endpoints, Picard version status, failed installer retries, reconstructed TLEN/flags and FrEIA exclusion. It has not been executed during implementation.

These cases from the original plan are not fully represented in those fixtures:

| Area | Remaining fixture and hand-derived expectation |
| --- | --- |
| Selected-tool setup | A cfDNAlab-only selection must not probe/install Picard, cfDNAPro, FinaleToolkit or Griffin. Shared samtools/setup-Python dependencies remain explicit. |
| Installer controls | Exercise the actual `--mark-outdated` option, a missing packed archive and successful resumption after interruption, preserving completed unrelated tools. |
| Native gap/D/N coverage | Forward `pos=100`, `30M`, reverse `pos=170`, `30M` spans 100 bp with a 40 bp inter-mate gap. Read depth covers 60 bp. Replacing forward CIGAR with `15M2D15M` leaves `[115,117)` unsupported in cfDNAlab fragment coverage, while mosdepth fragment mode fills the 100 bp span. |
| Length/orientation | Check ordinary inward bwa TLEN against directional spans, reverse-start-first overlap selection differences and integer truncation when samtools halves an odd passing-record count. |
| Motif range/window | A valid 20 bp pair is excluded by the normal 30 bp minimum and included by the wide diagnostic. A 1 Mb contig triggers the selected FinaleToolkit window guard. Test rank/log-ratio metrics, including qualifying zero frequencies. |
| Midpoint input | Small synthetic parser fixtures preserve `+`, `-`, `.`, repeats and overlaps, and reject invalid strands/numeric-only Griffin groups. No biological site selection or preparation helper is needed. |
| Failure rows/measurement | Placement failure, native failure, missing output and failed filesystem flush preserve the correct TSV field counts, statuses, missing measurements and available failed diagnostics. |
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
- [ ] Record the user's midpoint BED identity/checksum and provenance without changing its coordinates, groups, strands, site count or selection.
- [ ] Run native/proper-pair and wide-bound concordance on the actual remap. Measure filtered count effects from outputs rather than reusing the original BAM's sampled non-proper percentage. Keep biological differences/scopes explicit.
- [ ] Pilot cfDNAPro progressively from about 0.1x to choose a physical node with demonstrated headroom for all methods. Run full-input method/core pilots, splitting when necessary, and establish scratch capacity/filesystem behavior.
- [ ] Budget the full matrix with two warmups plus five measurements, margin and staging/cleanup. Pin production blocks to the same node and retain resource failures as outcomes.
- [ ] Collect production timings, completed accounting, retained profiles and measured points. Fill supplementary tables with sample sizes, mean/SD, median/quartiles, memory policy and workload notes. No final timing, capacity or concordance conclusion has been measured here.

## Optional later scope, still deferred

- deepTools fragment coverage, only with a scientific comparison justifying its selection/CIGAR differences.
- A 1x/5x/15x depth series after the primary fixed-input experiment succeeds.
- Bootstrap intervals, extra plots or drift-monitoring anchors if observed results justify them.

FrEIA concordance can be reconsidered only after native aggregation is fixed and checked. Recovery from its damaged aggregate remains excluded.
