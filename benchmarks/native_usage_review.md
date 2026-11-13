# Native usage review of the benchmark competitors

This review checks every competitor command in [experiment_commands.md](experiment_commands.md) against that tool's own documentation and source at the pinned version. The question is whether each command runs the tool the way its authors present it, or whether the benchmark changed it in ways the documentation does not propose. For midpoints, it also checks whether cfDNAlab is configured to reproduce Griffin's profiles. Reviewed at commit `a84c9a69` on 2026-11-12. No code was changed.

## Rules applied

- Shared biological settings may differ from a tool's default: MAPQ ≥30, fragment length 30-1000 bp (100-200 bp for midpoints), k, flag exclusions and contig scope. What matters is how a setting reaches the tool. A documented option of that subcommand is fine. A workaround the tool does not present, such as a prefilter pipeline, a generic htslib option on a subcommand that does not list it, or JVM flags, is flagged together with the documented alternative.
- Every native thread or worker option receives the allocated core count N. Parallelism the tool does not offer is not added.
- A documented optimization may be used. A mode the documentation recommends for normal use, but which the benchmark skips, is flagged.
- Added work inside a competitor's timer and work removed from it are both flagged, as are method labels that misdescribe the timed work.

Statements marked *inference* were derived from reading code or docs and were not confirmed by running anything.

## Summary

**Threads.** No N-1 remains. Every native thread or worker option receives N. Two settings add parallelism the tools do not offer: the `nthreads=` option injected into `samtools bedcov` and the JVM flag `-XX:ActiveProcessorCount` for Picard.

**Griffin.** The four timed Griffin scripts reproduce the official Snakemake rules argument by argument, apart from MAPQ, which is now decided as Griffin's 20 for all midpoint jobs. The problems are around those calls: a timed workload Griffin does not offer, the window used for concordance, and the undocumented preparation of the site BED.

**Status.** All 19 findings are decided as of 2026-11-13, including both High severity findings: `bedcov` runs without threads (finding 2), and Picard is timed without the samtools prefilter (finding 3). The agreed changes are implemented in the runners, comparators and documentation. Syntax and Cargo checks passed. Tests, installations and native cluster runs have not been executed during this implementation.

| # | Workload | Finding | Favours | Severity | Action |
| --- | --- | --- | --- | --- | --- |
| 1 | `mosdepth`, `mosdepth-regions` | `--fast-mode`, recommended by mosdepth for most use-cases, is never used. Decided: exact mode computes the matched quantity | cfDNAlab | Low | State in the methods why exact mode was used |
| 2 | `samtools-regions` | `nthreads=` is injected into `bedcov`, which disables thread options. samtools(1) says `nthreads=` is encoding-only for BAM. Decided: no threads | samtools | High | Run without threads and record `tool_threads=NA` |
| 3 | `picard` | The timed work is `samtools view \| picard`. Picard has no MAPQ option and its docs do not present a prefilter. Decided: time plain Picard | Mixed | High | Time the documented Picard call. Take the concordance histogram from the pipeline, untimed |
| 4 | `griffin` | The unit-GC-table run is not a Griffin workflow. Its uncorrected profile is identical to the one `griffin-gc` already writes. Decided: drop it | Griffin, in timing | Medium | Time `griffin-gc` against `cfdnalab-gc`. Take both Griffin tables from `griffin-gc` |
| 5 | Midpoint concordance | Comparing the full ±4995 bp window includes Griffin's edge bins, which its official output crops. Decided: official saved window | Against cfDNAlab | Medium | Compare on Griffin's official saved window after full-window normalization |
| 6 | Midpoint sites | The `prep-windows` step was undocumented. Decided: its 1 bp tie-break and the shared pre-filter are acceptable | Neither | Low | Done: documented in [experiment_commands.md](experiment_commands.md#site-preparation), with the sources of the two blacklists and `hg38.chrom.sizes` |
| 7 | `samtools-depth` | `-aa` writes zero-depth positions, while cfDNAlab's per-position track omits zero runs by default (`--keep-zero-runs` is not passed). Decided: drop `-aa` | cfDNAlab | Medium | Drop `-aa`, let the comparator accept sparse samtools output, and add the per-chromosome Jaccard truncation diagnostic |
| 8 | `pandepth-regions` | A converted BED with `-b` replaces PanDepth's documented fixed-window mode `-w`. Both modes are documented, and cfDNAlab and mosdepth use their own window options. Decided: `-w` | Unclear | Low | Use `-w 500000`. Adapt the comparator to its output file, and let it tolerate a missing 1 bp final window |
| 9 | `samtools` lengths | `samtools stats` times a full QC report, not a length histogram. Decided: label it | cfDNAlab | Medium | Label it "samtools stats (full QC report)" |
| 10 | `cfdnapro` motifs | `callMotif()` checks fragment contigs against the bundled UCSC hg38 BSgenome, which lacks the no-alt analysis set's `chrEBV`. The BAM has 13 MAPQ ≥30 reads there. `readBam()` does not check, so lengths are unaffected | Neither (failure) | Medium | Decided: pass every BAM contig present in cfDNAPro's reference to `chromosome_to_keep`, i.e. all but `chrEBV`. Other tools keep whole-BAM scope |
| 11 | samtools `depth`, `bedcov`, `stats` | Length and MAPQ filters use the generic htslib `filter=` option, while PanDepth gets no length filter. Decided: keep the samtools filters | samtools, mildly | Low-Medium | Describe them as samtools' global input option. PanDepth exposes no such option, so it keeps `length_filter=unsupported` |
| 12 | `picard` JVM | `-Xms512m` and `-XX:ActiveProcessorCount` are not presented. Picard's docs recommend `-Xmx2g`. Decided with finding 3 | Unclear | Low-Medium | Use `-Xmx2g` only, as in the finding 3 command |
| 13 | `freia` | `--subchroff` is not in the workflow rule and does nothing at this commit. Snakemake is bypassed and the run stops before group comparison and scoring. Decided: remove the flag | FrEIA, slightly | Low-Medium | Remove the flag and relabel the workload |
| 14 | `cfdnapro` lengths | `read_bam_insert_metrics(fragment_obj = )` instead of the `callLength()` route used in the paper and tutorial. Both run the same counting code. Decided: switch | Neither | Low | Switch to `callLength()` so the call matches the tutorial |
| 15 | `cfdnapro` motifs | The length filter is R subsetting, because cfDNAPro has no length option for motifs. Decided: keep it | Neither | Low | Label it as user-side subsetting |
| 16 | Griffin, minor | MAPQ 30 instead of Griffin's 20, `-4995 4995` instead of `-5000 5000`, untimed plotting rule, "eligible" Griffin site count. Decided: MAPQ 20 for all midpoint jobs, disclose the untimed plotting, relabel the site count | Mixed | Low | See the Griffin section |
| 17 | `mosdepth-regions` | Regional means keep the default two decimals. Decided: 5 decimals | Neither | Low | Set the documented `MOSDEPTH_PRECISION=5` |
| 18 | FrEIA `--both` | The right-end sequence is reversed but not complemented. The runner already uses left ends only | Neither | Low | None. The upstream FrEIA project has not responded to the first report |
| 19 | `picard`, minor | `MINIMUM_PCT=0` is unnecessary. The documented example renders a histogram PDF, which this call skips. Decided with finding 3 | Picard, slightly | Low | The finding 3 command omits `MINIMUM_PCT`. Disclose the missing PDF |

## Midpoints: Griffin v0.2.0

Sources: the [profiling config](https://github.com/adoebley/Griffin/blob/v0.2.0/snakemakes/griffin_nucleosome_profiling/config/config.yaml), [profiling Snakefile](https://github.com/adoebley/Griffin/blob/v0.2.0/snakemakes/griffin_nucleosome_profiling/griffin_nucleosome_profiling.snakefile), [cluster config](https://github.com/adoebley/Griffin/blob/v0.2.0/snakemakes/griffin_nucleosome_profiling/config/cluster_slurm.yaml), [GC config](https://github.com/adoebley/Griffin/blob/v0.2.0/snakemakes/griffin_GC_and_mappability_correction/config/config.yaml), [GC Snakefile](https://github.com/adoebley/Griffin/blob/v0.2.0/snakemakes/griffin_GC_and_mappability_correction/griffin_GC_and_mappability_correction.snakefile), [griffin_coverage.py](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_coverage.py), [griffin_merge_sites.py](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_merge_sites.py) and the [pinned TFBS site lists](https://github.com/adoebley/Griffin_analyses/tree/75467e438a88029206d8d778ddce05eb0ab6efa8/sites/TFBS/10000_unfiltered_sites_CIS_BP_v2).

### Timed calls against the official rules

| Official rule argument | Official value | Benchmark | Verdict |
| --- | --- | --- | --- |
| `calc_cov --GC_bias` | Fitted table from the GC workflow | Unit table for `griffin`, fitted table for `griffin-gc` | Unit table is not an upstream configuration. See finding 4 |
| `--mappability_correction`, `--mappability_bias` | `False`, `none` | Same | As presented |
| `--reference_genome`, `--mappability_bw` | `hg38.fa`, `k100.Umap.MultiTrackMappability.bw` | Staged FASTA, same bigWig | As presented |
| `--chrom_sizes_path` | `hg38.standard.chrom.sizes` | BAM sizes for the BED's chromosomes, lexicographic order | No effect on output |
| `--chroms` | `chr1` to `chr22` | Chromosomes found in the BED | Same if the BED is autosomal. See finding 6 |
| `--norm_window` | `-5000 5000` | `-4995 4995` | Same output. See the minor items |
| `--size_range` | `100 200` | `100 200` | As presented |
| `--map_quality` | `20` | `30` | Decided: 20 for all midpoint jobs. See the minor items |
| `--number_of_sites`, `--sort_by`, `--ascending` | `none` | Script defaults | Equivalent: all sites |
| `--CPU` (both scripts) | 8 in the cluster config, "for parallelization" | N | As presented |
| `merge_sites --save_window` | `-1000 1000` | `-4995 4995` | Changes the compared bins and the `mean_coverage` feature. See finding 5 |
| `--center_window`, `--fft_window`, `--fft_index` | `-30 30`, `-960 960`, `10` | Same | As presented |
| `--smoothing`, `--smoothing_length`, `--step` | `True`, `165`, `15` | Same | As presented |
| `--exclude_paths` | Five Griffin `Ref/` BEDs | Same five | As presented |
| `--exclude_zero_mappability`, `--exclude_outliers` | `True`, `True` | Same | As presented |
| `--CNA_normalization`, `--individual` | `False`, `False` | Same | As presented |
| `generate_plots` rule | Part of `rule all` | Not run | See the minor items |
| `GC_counts` | Mappable-regions BED, reference, MAPQ 20, `15 500`, CPU 8 | Same BED and range, MAPQ 30, CPU N | As presented apart from MAPQ, decided as 20 |
| `GC_bias` | `--mappable_name`, `--genome_GC_frequency`, `15 500` | Same | As presented |
| `mappability_bias` rule | Only runs with mappability correction, which the GC config does not recommend | Not run | As presented |

Calling the scripts directly instead of through Snakemake is reasonable. Snakemake only schedules these commands.

### Finding 4: the uncorrected `griffin` workload is not a Griffin workflow

The profiling Snakefile loads `config/samples.GC.yaml` with the comment "this is generated by griffin_GC_correction", and its `rule all` requires `GC_corrected.coverage.tsv`. Griffin's presented pipeline is therefore GC estimation followed by profiling, which writes both an uncorrected and a GC-corrected table. The unit GC table is a construct of the runner's setup.

The uncorrected table does not depend on the GC table. `griffin_coverage.py` counts every passing midpoint with `cov_dict[midpoint]+=1` ([L395](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_coverage.py#L395)). The bias lookup only feeds the separate `GC_cov_dict`. The outlier mask in `griffin_merge_sites.py` is built from the uncorrected coverage and applied to both. A `griffin-gc` run therefore writes the same `sample.uncorrected.coverage.tsv` as the `griffin` run, for the same BAM, sites and settings.

Consequences:

- **Timing.** `griffin` times Griffin's profiling scripts with compulsory GC weighting but no GC estimation. That is neither the workflow Griffin presents nor work matched to uncorrected cfDNAlab.
- **Compute.** The README budgets roughly 30 hours per full Griffin run, for three measurements at each of four core settings, for both `griffin` and `griffin-gc`. Dropping the `griffin` jobs halves Griffin's compute without losing a profile.

**Decision (2026-11-13).** Drop the `griffin` jobs. Report Griffin's timing from `griffin-gc` against `cfdnalab-gc`. Preserve both Griffin tables from the first measured `griffin-gc` repetition, and compare its uncorrected table with `cfdnalab`. [compare_midpoints.py](midpoints/compare_midpoints.py) currently reads the uncorrected table from the `griffin` directory, so it must read it from `griffin-gc` instead.

### Finding 5: the compared window includes Griffin's edge bins

Griffin smooths the whole normalization window with `savgol_filter(results[norm_columns], savgol_window, 3)` ([L477](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_merge_sites.py#L477)). SciPy's default `mode='interp'`, also in the pinned SciPy 1.7.1, fills the outer 5 bins on each side from a polynomial fitted to the last 11 bins. Griffin then normalizes over the full window ([L481](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_merge_sites.py#L481)) and keeps only `save_columns` ([L485](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_merge_sites.py#L485)). With the official `-1000 1000`, rounded inward to `[-990, 990)` or 132 bins, Griffin's output never contains those edge bins. cfDNAlab counts flank positions and smooths with real data before trimming.

Setting `--save_window` equal to the normalization window puts 10 of the 666 compared bins, the outer 75 bp on each side, under different edge handling. *Inference:* those bins are likely to set `max_abs_mean_one` and to raise MAE and RMSE for reasons unrelated to midpoint counting.

The saved window never affects Griffin's normalization, because Griffin always normalizes over the normalization window before cropping. The README's reason for saving the full window, "so cropping does not change the normalization denominator", is therefore unnecessary. The saved window does change the `mean_coverage` feature ([L497](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_merge_sites.py#L497)), which becomes about 1 by construction. This matters only if Griffin's features are reported.

**Decision (2026-11-13).** Run Griffin with the official `--save_window -1000 1000`. In the comparator, normalize cfDNAlab over all 666 bins as now, then select the 132 bins `[-990, 990)`. Two comparator checks change:

- [compare_midpoints.py L85-86](midpoints/compare_midpoints.py#L85-L86) require Griffin's position labels to cover the whole cfDNAlab window. They should require exactly the saved bins instead.
- The mean-one check on Griffin's profile ([L114](midpoints/compare_midpoints.py#L114)) no longer holds, because Griffin normalizes before cropping. Drop it, or check the mean of the uncropped window if Griffin is also asked to save it.

### Finding 6: the site BED preparation is undocumented and shifts anchors

[experiment_commands.md](experiment_commands.md#midpoint-profiles) and [midpoints/README.md](midpoints/README.md) start from "supplied" 9,990 bp windows. Nothing in `benchmarks/` records how they were made. The preparation command open in the editor is:

```bash
cfdna prep-windows -i griffin_tfbs.cleaned.bed -o griffin_tfbs.prepared_9990bp.bed \
    --sep tab --group-cols 3 \
    --blacklist hg38-blacklist.v2.bed --blacklist k100.umap.exclusion_intervals.bed \
    --blacklist-halo 200 --resize 9990 --chrom-sizes hg38.chrom.sizes
```

**Anchor shifts.** Griffin's TFBS files have the columns `Chrom Start End uniprotId Gene exp.count peak.count position` and no strand column. Across all 10,000 rows of the CTCF list, `position = floor((Start + End) / 2)`. Griffin centers each site at `position`. `prep-windows --resize` ([resizers.rs L54-82](../src/commands/prepare_windows/resizers.rs#L54-L82)) centers on `start + floor(length / 2)`, which equals `position`. When the interval length is odd and the target length is even, it picks between two placements with a hash. The right-hand placement starts the window at `position - 4994`, so the runner's Griffin anchor becomes `position + 1`. The first two CTCF intervals are 37 and 39 bp long, so odd lengths are common. *Inference:* about half of the odd-length sites move by 1 bp.

Both tools receive the same shifted anchor, so tool-against-tool concordance is unaffected. The sites are no longer Griffin's sites at Griffin's positions, though. Writing windows directly as `[position - 4995, position + 4995)` reproduces Griffin's anchoring exactly. Passing 2 bp intervals `[position - 1, position + 1)` to `prep-windows` does the same, because even lengths have a unique placement. This assumes the cleaned BED keeps Griffin's Start and End. A cleaned BED of 1 bp intervals at `position` has the same problem.

**Pre-filtering.** The command removes windows within 200 bp of `hg38-blacklist.v2` and a k100 Umap exclusion BED. Griffin's own analyses use all 10,000 sites per TF and handle problem regions by bin masking, using the ENCODE unified exclusion list, centromeres, gaps, fix patches, alternative haplotypes, zero mappability and outliers. Neither pre-filter list is one Griffin uses. Griffin still runs natively on the sites it receives. However, the pre-filter removes most of the cases where the tools' native masking differs, because cfDNAlab drops whole sites while Griffin masks bins. High concordance would then partly reflect a cfDNAlab-style site selection.

**Chromosomes.** Griffin's config restricts sites to `chr1`-`chr22`. The runner passes the chromosomes found in the BED. The CTCF list contains only `chr1`-`chr22`, so this is identical if the cleaned BED is autosomal.

**Decision (2026-11-12).** The shifted anchors and the pre-filter stay, because both tools receive the same windows and the benchmark tests whether they produce the same profiles. For even-length intervals, `prep-windows` reproduces Griffin's placement exactly. For odd-length intervals, Griffin always places the motif center on the right of the window's two central bases, an offset of half a base. The random tie-break removes that systematic offset on average.

Done: the derivation is now in [experiment_commands.md](experiment_commands.md#site-preparation), including a reconstruction of the lost name-cleaning command and the sources of `hg38-blacklist.v2.bed`, `k100.umap.exclusion_intervals.bed` and `hg38.chrom.sizes`. Check `cut -f 1 | sort -u` on the final BED.

### cfDNAlab configuration against Griffin

The matched settings are the 100-200 bp inclusive length range (`--length-bins 100 201`), 15 bp bins aligned with Griffin's signed bins, the 165 bp order-3 Savitzky-Golay filter (cfDNAlab's default was taken from Griffin), the five exclusion BEDs, the chromosomes and MAPQ, and mean-one group normalization in the comparator.

These native differences cannot be removed through cfDNAlab options. All four are already disclosed in [midpoints/README.md](midpoints/README.md):

1. **Even-length midpoints.** Griffin uses `floor((start + start + TLEN) / 2)` ([L373](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_coverage.py#L373)), which is the right central base. cfDNAlab picks either central base reproducibly. *Inference:* the mean offset is about 0.25 bp, which is negligible at 15 bp resolution.
2. **Fragment selection.** Griffin counts forward records with positive TLEN, and checks MAPQ, duplicate, QC-fail and paired status on that record only ([L367-370](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_coverage.py#L367-L370)). cfDNAlab checks both mates.
3. **Masking.** Griffin masks a site's bins that overlap the exclusion BEDs, bins with any zero-mappability base, and bins with a z-score of at least 10 across the site-by-bin matrix ([L428](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_merge_sites.py#L428)). It then averages the remaining sites per bin. cfDNAlab removes intervals within ceil(200 / 2) + 82 = 182 bp of the exclusion BEDs and drops fragments that overlap them. It has no mappability or outlier mask.
4. **Smoothing.** Griffin smooths 15 bp bin sums with an 11-bin filter. cfDNAlab smooths base-resolution counts with a 165 bp filter, then averages bins.

To show that the remaining differences come from these rules rather than from counting, add an untimed diagnostic. Run cfDNAlab with `--smoothing none --bin-size 15`, then apply Griffin's own `savgol_filter(profile, 11, 3)` and mean-one normalization in the comparator. This isolates difference 4 and needs only cheap cfDNAlab reruns.

Reporting fix: [compare_midpoints.py](midpoints/compare_midpoints.py) and the [concordance spec](benchmarks_concordance_spec.md) describe `griffin_sites` as an eligible site count. Griffin writes `number_of_sites = len(current_sites)` ([L572](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_merge_sites.py#L572)), which counts every site on the selected chromosomes, including fully masked sites. The site-preparation pre-filter has probably removed most or all fully masked sites, so the numbers may barely differ. **Decision (2026-11-13):** relabel the column as the number of sites in Griffin's site list.

### GC-corrected Griffin

The four timed calls match the `GC_counts`, `GC_bias`, `calc_cov` and `merge_sites` rules. They use the mappable-regions BED, the reference GC-frequency directory, the 15-500 bp estimation range and CPU N. MAPQ 30 against the config's 20 is the only difference, and it is decided as 20. `--chrom_sizes` receives the BAM sizes instead of `hg38.standard.chrom.sizes`, and `griffin_GC_counts.py` reads that argument but never uses it.

The two GC models differ by design. Griffin models per-length GC counts over the whole fragment. cfDNAlab models GC percentages with the reference package's end offset over 30-1000 bp. Expect lower agreement for corrected profiles than for uncorrected ones. That is a fair native comparison and needs no change.

### Minor Griffin items

- **MAPQ.** Both Griffin workflows are configured with 20. **Decision (2026-11-13):** every midpoint job uses MAPQ 20, Griffin's own setting, for both tools and for the GC estimation steps. The rest of the benchmark keeps the shared 30.
- **Normalization window.** `--norm_window -4995 4995` replaces the official `-5000 5000`, to match the 15 bp bins of the 9,990 bp windows. Griffin rounds the official window inward to the same `[-4995, 4995)` itself ([L225](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_merge_sites.py#L225)). The 100 bp fetch padding covers the counted positions in both cases, so the profiles are identical. No change is needed. Passing the official values instead would let a reader match the command to Griffin's config without knowing about the rounding.
- **Plotting.** The official profiling `rule all` includes `generate_plots` (`griffin_plot.py`), which the benchmark neither runs nor times. cfDNAlab's default QC plot is timed. This favours Griffin. **Decision (2026-11-13):** leave both as they are and disclose it.
- **Parallelism.** `griffin_merge_sites.py` parallelizes only across site lists, so a group cannot use more than a worker. This is native and already noted.

## Coverage

### mosdepth 0.3.14 (finding 1)

The v0.3.14 help ([mosdepth.nim L905](https://github.com/brentp/mosdepth/blob/v0.3.14/mosdepth.nim#L905)) says: "-x --fast-mode dont look at internal cigar operations or correct mate overlaps (recommended for most use-cases)." The README's WGS window example is `mosdepth -n --fast-mode --by 500 sample.wgs $sample.wgs.cram`, explained as: "--fast-mode avoids the extra calculations of mate pair overlap and cigar operations, ... providing a substantial speed improvement."

Neither `mosdepth` nor `mosdepth-regions` uses it. Fast mode cannot be combined with `--fragment-mode` (mosdepth exits with an error), so `mosdepth-fragment` is correct as it is.

Fast mode changes the output. It counts overlapping mates twice and counts D and N positions as covered. The exact mode is therefore the right concordance partner for cfDNAlab's overlap-merged read depth. Timing only the exact mode, though, shows mosdepth in a mode its authors do not recommend for most use. A reviewer who knows mosdepth is likely to notice.

For regions, the regional peers `samtools bedcov` and PanDepth already count overlapping mates twice. `mosdepth -n --fast-mode --by 500000` is both mosdepth's own recipe and closer to those peers.

**Decision (2026-11-12).** No fast-mode rows. Exact mode computes the same quantity as cfDNAlab, with mate overlap counted once and deletions excluded, so it is the matched workload for both per-base and regional read depth. Add a sentence to the methods saying so, because mosdepth's own help recommends fast mode.

The other mosdepth options are documented and used as presented: `--threads N` ("number of BAM decompression threads". The README notes no benefit beyond about 4 threads, without advising fewer), `--include-flag 1`, `--flag 3852`, `--mapq`, `--min-frag-len`/`--max-frag-len` (inclusive on absolute TLEN), `--no-per-base` and `--by`.

Finding 17: the README documents `MOSDEPTH_PRECISION=5` for "5 decimals of precision rather than the default of 2". **Decision (2026-11-13):** set `MOSDEPTH_PRECISION=5`. The comparator multiplies each rounded mean by the window length, so with two decimals a 500 kb window's total can be off by up to 2,500, or about 0.02% of the window's total at 30x. Five decimals reduce that to 2.5. The default-precision notes in [experiment_commands.md](experiment_commands.md) and [coverage/README.md](coverage/README.md) change with it. The README also presents `--d4` output (separate `mosdepth_d4` binary) as a large per-base write speed-up. That is an optional extra row, since its output format would need a separate reader.

### samtools 1.25 `depth` (finding 7)

- `-@ N` appears in `depth`'s usage text ("Number of additional threads to use"), not in the man page. Fine.
- `--require-flags 1`, `-G 3852`, `-Q 30`, `-q 0` and `-s` are documented. No `-J` correctly excludes deletions.
- `-aa` is documented but not the default. It writes a line for every zero-depth base and every base of contigs without reads. cfDNAlab's per-position track omits zero-coverage runs by default ("By default, only covered positions are written to the output", [config.rs](../src/commands/fcoverage/config.rs)), and the benchmark does not pass `--keep-zero-runs`. Both default outputs contain the same information, because absent positions on a processed contig are zero. `-aa` therefore makes samtools write content that cfDNAlab does not write. *Inference:* the cost scales with the zero-depth fraction, which is small for deep WGS and large for low-pass data.
- **Decision (2026-11-12).** Drop `-aa`. samtools then writes every position with depth from passing reads, as cfDNAlab does. The comparator must accept that: [compare_coverage.py L86](coverage/compare_coverage.py#L86) currently treats only `cfdnalab*` tracks as sparse and rejects a samtools table that skips positions. The comparator already treats absent positions in sparse tracks as zero through `zero_filled_runs()`, so only that line needs to change.

Allowing sparse samtools output gives up the dense-format truncation check. Two checks already remain: a non-zero samtools exit fails the run, and the comparator rejects cut-off lines. A truncation that ends on a line boundary would only show up as missing coverage. Replacement, decided on 2026-11-13: count bases covered by both tools in the same pass that already counts `nonzero_union_bases`, and report the Jaccard index (both / either) genome-wide and per chromosome. samtools writes in BAM order, so truncation has a distinctive pattern: one chromosome with a partial Jaccard and every later chromosome at 0. Native differences instead lower all chromosomes similarly. Report these as diagnostics and flag that pattern, rather than failing on an arbitrary threshold, because tools legitimately differ in which bases they cover.
- The TLEN filter uses `--input-fmt-option filter=...`. `depth` has no fragment length option. The generic filter is documented in [samtools(1)](https://www.htslib.org/doc/samtools.html#GLOBAL_COMMAND_OPTIONS) and listed in `depth`'s usage text, but not in the `depth` man page. It saves an extra pass. See finding 11.

### samtools 1.25 `bedcov` (finding 2)

- `bedcov` disables the global thread option: `SAM_OPT_GLOBAL_OPTIONS('-', 0, '-', '-', 0, '-')` ([bedcov.c L143](https://github.com/samtools/samtools/blob/1.25/bedcov.c#L143)), where the last position is `--threads`. Its usage and man page list no thread option. [samtools(1)](https://www.htslib.org/doc/samtools.html#GLOBAL_COMMAND_OPTIONS) says of `nthreads=`: "For BAM this will be encoding only."
- *Inference from htslib source:* the injected `nthreads=` still starts BAM decoding threads, so it speeds `bedcov` up in a way samtools does not present.
- `timings.tsv` records `tool_threads=$cores` for this method, and the runner comment and experiment_commands.md call it a native thread option. That metadata is wrong.
- **Decision (2026-11-13).** Run `bedcov` without threads and record `tool_threads=NA`, as for cfDNAPro and `frag-length-bins`. Correct the runner comment and [experiment_commands.md](experiment_commands.md) to match.
- The paired and TLEN filter is also the generic htslib option. `bedcov` requires an index, so a `view | bedcov` pipeline is not possible. **Decision (2026-11-13):** keep it, described as samtools' global input option rather than a `bedcov` option. The same applies to `depth` and `stats`. PanDepth exposes no htslib input options, so it keeps `length_filter=unsupported`.
- `-j`, `-Q 30` and `-G 3852` are documented.

### PanDepth 2.26 (finding 8)

- `-t N` (default 3), `-q 30` and `-x 3852` are documented. The FAQ explicitly presents `-x` for custom flag filters. *Inference from source:* PanDepth caps its workers at the number of contigs longer than 10% of the longest contig, so 32 cores gain little over about 24.
- No length or paired filter exists, and the runner correctly records `length_filter=unsupported`.
- `-a` (per-site output) is in the help but not among the README examples. It is PanDepth's only per-base route, so keep it. Disclose that its per-base time includes a single-threaded gzip of every base, which happens after the workers finish (source).
- `-b` with a 1-based copy of the shared BED replaces the fixed-window mode the README presents in example 3.4: `pandepth -i test.bam -w 100 -o test4`. Both modes are documented, so `-b` is not a hack. The case for `-w 500000` is consistency: cfDNAlab and mosdepth make their 500 kb windows with their own options (`--by-size`, `--by`), while only FinaleToolkit and `samtools bedcov`, which have no window option, read the shared BED. `-w` also removes the coordinate conversion. **Decision (2026-11-13):** use `-w 500000`. *Inference from source:* a 1 bp tail may be dropped when the contig length is 1 more than a multiple of 500,000. Its effect on the values is negligible. The comparator, though, requires PanDepth's windows to match the shared BED exactly ([compare_coverage.py L361-364](coverage/compare_coverage.py#L361-L364)), so it would stop with "Regional row identity differs". When implementing, compare PanDepth on the windows it reports and report how many shared windows it lacks. The comparator also parses the `-b` output file (`pandepth.bed.stat.gz`), so it needs a reader for the `-w` output.

### FinaleToolkit 1.1.0 `coverage`

As presented. `--threads N` is the documented worker-process option. The default midpoint policy is what the [features page](https://epifluidlab.github.io/FinaleToolkit/documentation/user_guide/features.html) recommends for binned, genome-wide features. Normalization is off by default, with a scale factor of 1. The read1-overlap loss described in experiment_commands.md matches the source.

## Fragment lengths

### samtools 1.25 `stats` (finding 9)

- `-@ N`, `-f 1`, `-F 3852`, `-m 1` and `-i 1001` are documented. The overflow-bin behavior of `-i` is source-only.
- `stats` has no MAPQ option. Its `-q` is "The BWA trimming parameter". The MAPQ filter uses the generic htslib `filter=` option, which is kept (finding 11). The alternative would be `samtools view -u -q 30 -f 1 -F 3852 sample.bam | samtools stats -@ N -m 1 -i 1001 -`.
- `stats` always computes its full report: checksums, per-cycle qualities, GC, read lengths, indels, the coverage distribution and more. That disadvantages samtools against a length-only tool. It is the only native samtools route, so keep it, and label it "samtools stats (full QC report, insert-size section extracted)" in figures and text. Decided on 2026-11-12.

### Picard 3.5.0 `CollectInsertSizeMetrics` (findings 3, 12 and 19)

The documented example is `java -jar picard.jar CollectInsertSizeMetrics I=input.bam O=insert_size_metrics.txt H=insert_size_histogram.pdf M=0.5` ([source L77-83](https://github.com/broadinstitute/picard/blob/3.5.0/src/main/java/picard/analysis/CollectInsertSizeMetrics.java#L77-L83)).

**The prefilter pipeline.** Picard has no MAPQ argument, and its docs never suggest prefiltering. Picard already skips first-of-pair, unpaired, unmapped, mate-unmapped, secondary, supplementary and duplicate records and records with TLEN 0. The samtools prefilter therefore only adds MAPQ ≥30, effectively on the second-of-pair record Picard counts, and QC-fail exclusion. In return, the timed `picard` method becomes a two-process pipeline:

- At 4 to 32 cores, samtools decodes the BAM with N threads, so most of the "Picard" scaling curve is samtools. *Inference.*
- At 1 core, the extra decode and uncompressed re-encode pass slows Picard down.
- Natively, Picard stops reading at the unplaced, unmapped tail of the BAM, while `samtools view` reads the whole file.
- *Inference:* samtools runs N+1 threads and the JVM is sized for N, all on N cores.

Decision (2026-11-12): time the documented call, `picard -Xmx2g CollectInsertSizeMetrics I=sample.bam O=metrics.txt W=1000`, and disclose that it keeps MAPQ<30 and QC-fail pairs. Produce the concordance histogram from the samtools pipeline in an untimed step.

**JVM flags.** Picard's [command-line overview](https://broadinstitute.github.io/picard/command-line-overview.html) says: "We do not provide guidance for specifying JVM (Java Virtual Machine) arguments except for specifying memory: most of the commands are designed to run in 2GB of JVM, so we recommend using the JVM argument -Xmx2g."

- `-Xms512m` makes the Bioconda launcher drop its default `-Xms512m -Xmx1g`, so Java picks a quarter of available memory.
- `-XX:ActiveProcessorCount=N` is a JVM flag, not a Picard option, and should not be described alongside native thread options. *Inference:* since the runner already verifies an N-CPU set, it probably has no effect.
- `-Djava.io.tmpdir` corresponds to Picard's documented `TMP_DIR=`. *Inference:* this tool probably writes no temporary files.
- **Decision (2026-11-13), with finding 3:** use `-Xmx2g` and drop the other two flags. `CollectInsertSizeMetrics` streams the BAM and keeps only histograms in memory, so a 2 GB heap does not limit its speed.

**Minor settings.**

- `HISTOGRAM_WIDTH=1000` is the documented way to stop tail truncation and is needed.
- `INCLUDE_DUPLICATES=false` is the default.
- `INPUT=/dev/stdin` is documented in the [FAQ](https://broadinstitute.github.io/picard/faq.html).
- `MINIMUM_PCT=0` differs from both the default (0.05) and the example (0.5). It only discards orientation categories below the threshold, and FR dominates cfDNA, so omit it.
- `HISTOGRAM_FILE` became optional in 3.5.0. The documented example renders the PDF, which costs R startup and plotting. picard-slim has no R, so skipping it is necessary. Disclose that cfDNAlab's default plot is timed while Picard's is not.
- Not used: `VALIDATION_STRINGENCY=SILENT`, which the docs say "can improve performance". The default is defensible. Mention it only if asked whether Picard was tuned.

### FinaleToolkit 1.1.0 `frag-length-bins`

As presented. `--bin-size 1`, `--min-length`, `--max-length` and `--min-mapq 30` are documented. They match the authors' own Figure 2a command, which was run on a fragment file. The command has no worker option at 1.1.0, so `tool_threads=NA` is correct.

### cfDNAPro 1.7.4 lengths (finding 14)

- The container route with `--cleanenv` follows the README's recommended installation. `--containall` is not shown there, but no effect is expected.
- No cfDNAPro function used here has a thread or BiocParallel argument, so there is nothing to set to N.
- `readBam()` settings:
  - `genome_label`, `chromosome_to_keep = FALSE` (default `chr1`-`chr22`), `curate_start_and_end` and `galp_flag` are documented arguments used to harmonize.
  - `galp_mapqFilter = 30` matches the README, tutorial and paper ("mapping qualities less than 30 were discarded"), although the code default is 40.
  - As a result, `readBam(bamfile)` as written in the Quick Start runs at MAPQ 40, hg19 and `chr1`-`chr22`.
- `read_bam_insert_metrics(fragment_obj = ...)` uses a documented argument. The paper's Methods and readthedocs tutorial 1 present `readBam(...) |> callLength(...)` instead. With a fragment object, both functions run the same counting code line for line: `width()`, inclusive bounds and the same zero-filling of missing lengths. The only difference is that `callLength()` adds a `prop` column when it zero-fills. The runner comment explains why `fragment_obj` is passed (the `bamfile` shortcut would drop the MAPQ and flag settings) but not why this function was chosen. **Decision (2026-11-13):** switch to `callLength(fragments, isize_min = 30L, isize_max = 1000L)` so the call matches the published route. Timing and counts do not change. Check that the comparator tolerates the extra `prop` column.

## End motifs

### FinaleToolkit 1.1.0 `end-motifs`

The command is as presented: `-k 4` is the default, `--strand both` is the default, `--threads N` is documented, and `--min-mapq 30` replaces the CLI default of 20 to harmonize. The ignored length bounds are confirmed in source, and the maintainers confirm them in [PR #194](https://github.com/epifluidlab/FinaleToolkit/pull/194): "`region_end_motifs` took `fraction_low`/`fraction_high` but never used them ... v0.12.0 through v1.1.0 are affected". The same PR confirms double counting of fragments that straddle 1 Mb window boundaries.

**Decision (2026-11-13): use FinaleToolkit `main` at the PR #194 merge commit.** PR #194 ("fix(motifs): restore fragment-length filter and midpoint membership", fixing [issue #192](https://github.com/epifluidlab/FinaleToolkit/issues/192)) was merged into `main` on 2026-11-12 as `d16632ec`. No release contains it yet; the latest release is still v1.1.0. At that commit:

- `--min-length` and `--max-length` apply as inclusive bounds, so the end-motif length filter becomes native.
- Each fragment is counted only in the window that contains its midpoint, so the double counting is gone.
- A much smaller boundary loss remains. The window query fetches read1 alignments that overlap the window and then keeps fragments whose midpoint lies inside it. A fragment longer than twice the read length can have its midpoint in one window and its read1 entirely in the neighbouring window, so neither window counts it. FinaleToolkit `coverage` already has the same behaviour, which [experiment_commands.md](experiment_commands.md) discloses.
- The window builder (`_genome_window_args`) is unchanged, so contigs whose length is an exact multiple of 1 Mb still lose their last full window. Keep the runner's guard.
- Between v1.1.0 and `d16632ec`, only motif, DELFI and release files changed. `coverage` and `frag-length-bins` run the same code as v1.1.0, so a single installation from that commit serves all three FinaleToolkit commands.

FinaleToolkit's wiki presents `filter-file --min-length --max-length` as a prefilter. Its BAM length filter compares signed TLEN, though, so it removes negative-TLEN records, which include about half of the read1 records this command uses. It is not a usable workaround, and with the fix it is no longer needed.

### cfDNAPro 1.7.4 motifs (findings 10 and 15)

`callMotif(motif_type = "s", motif_length = k)` is the documented entry point, and the paper uses `s` motifs with `motif_length = 3`. cfDNAPro has no native both-end 5′ motif. Type `e` returns the rightmost forward-strand bases without reverse complementing. Left-end counting with the disclosed side-mirroring approximation is therefore justified.

**Length filter.** Neither `readBam()` nor `callMotif()` has a length argument. Only the length functions take `isize_min` and `isize_max`. The runner therefore subsets fragment widths in R. The paper says only 50-450 bp fragments were kept but does not show how. **Decision (2026-11-13):** keep the subsetting inside the timed script, and describe it in [ends/README.md](ends/README.md) as user-side R subsetting because cfDNAPro offers no length filter for motifs. It is a single vectorized comparison, so its cost is negligible next to `readBam()`.

**Contig risk.** [data/README.md](data/README.md) uses the GRCh38 no-alt analysis set with UCSC names, which contains `chrEBV`. `BSgenome.Hsapiens.UCSC.hg38` does not. `readBam()` only attaches the reference, so the lengths workload is unaffected. `callMotif()` calls `get_genome_reference()`, which stops with "The supplied genome is missing fragment chromosomes" (`genomeReference.R` L82-92). This is not cfDNAPro mishandling BAM files. It follows from the harmonized `chromosome_to_keep = FALSE`, because the default `chr1`-`chr22` never reaches `chrEBV`. Whether it fails depends on whether any MAPQ ≥30 pair maps there. `samtools view -c -q 30 sample.bam chrEBV` returns 13, so the whole-BAM motif call is expected to fail.

**Decision (2026-11-13).** An autosome-only cfDNAPro run would make the all-against-all motif table mix scopes, and restricting every motif tool to autosomes is not attractive. Instead, keep whole-BAM scope and drop only the contigs cfDNAPro's reference lacks. `chromosome_to_keep` is a documented argument that accepts any vector of contig names. `curate_bam_seqinfo()` prunes the other contigs with `keepSeqlevels(..., pruning.mode = "coarse")` before `callMotif()` validates the remaining contigs. For this BAM that removes only `chrEBV`, at most 13 MAPQ ≥30 reads, against whole-BAM counts in the other tools. No extra table is needed, and the other tools are unchanged. Implementation, when code changes are allowed:

- Pass the BAM contigs that exist in `BSgenome.Hsapiens.UCSC.hg38` as `chromosome_to_keep` in [cfdnapro_end_motifs.R](ends/cfdnapro_end_motifs.R), computed from the two contig lists rather than by hard-coding `chrEBV`.
- Record the scope as "BAM contigs present in cfDNAPro's reference" for `cfdnapro` and disclose the dropped contigs.

The lengths workload keeps `chromosome_to_keep = FALSE`, because `readBam()` does not validate contigs.

Two documented fixes exist:

- `genome_label = "user_define"` with a BSgenome built from the exact alignment FASTA. This also settles the sequence-identity question in D6.
- An explicit `chromosome_to_keep` limited to contigs present in both references, with the scope disclosed.

### FrEIA `2c32fb9` (findings 13 and 18)

Sources: [README](https://github.com/mouliere-lab/FrEIA/blob/2c32fb9d9fab774f4ae57931df7fd7d75654d343/README.md), [config](https://github.com/mouliere-lab/FrEIA/blob/2c32fb9d9fab774f4ae57931df7fd7d75654d343/config/config.yaml), [preprocessing rules](https://github.com/mouliere-lab/FrEIA/blob/2c32fb9d9fab774f4ae57931df7fd7d75654d343/workflow/rules/2_FrEIA/1_preprocessing.smk) and [FrEIA rules](https://github.com/mouliere-lab/FrEIA/blob/2c32fb9d9fab774f4ae57931df7fd7d75654d343/workflow/rules/2_FrEIA/2_FrEIA.smk).

**As presented:**

- The filter, `tee` and index pipeline matches rule `filter_and_index_FrEIA` exactly, including both `-@ N`.
- `-q 30` replaces the config's `FiltQual: 5` to harmonize.
- QC-fail records are not excluded, which is faithful to the rule. Leave it.
- `-c <all contigs>` replaces the config's `chr1`-`chr22` plus `chrM` to harmonize.
- `-n 3` is the config's `nBase`.
- `-t N` matches the rule's threads.
- `-fsmin 30 -fsmax 1001` use the config keys, and the +1 is correct because FrEIA uses `range()`.
- `-subs 1` and `--bootstrap_sample 0` are the config defaults that the rule passes.

**`--subchroff`.** Its help text says "Turn off subchromosomal analysis. Subchromosomal analysis can use a lot of memory". The `data_transformation` rule does not pass this flag. At this commit the subchromosomal code is commented out and both branches select `Genome_Lvl` (`2_data_transformation.py` L352-358 and L396-399), so it has no effect. [ends/README.md](ends/README.md) lists it as part of the native call. **Decision (2026-11-13):** remove it.

**Snakemake.** The README presents FrEIA as `snakemake --use-conda --cores N -s 1_preprocessing.smk`, then `-s 2_FrEIA.smk`. The benchmark runs the rules' shell commands directly and stops after `data_transformation`, before `compare_groups` and `FrEIA_score`. That is a reasonable way to isolate per-sample work, and `3_Abundances/Genome_Lvl` is an output the README lists. It leaves out Snakemake and conda startup, which slightly favours FrEIA. Relabel the workload as "FrEIA rules filter_and_index_FrEIA, FragmEndSeq and data_transformation, run without Snakemake" wherever it is called a complete native workflow.

**Additional `--both` defect (finding 18).** `getReadSeq` takes the right end as `query_alignment_sequence[-nrBase:][::-1]` ([1_extract_fragment_ends.py L44-50](https://github.com/mouliere-lab/FrEIA/blob/2c32fb9d9fab774f4ae57931df7fd7d75654d343/workflow/scripts/1_FrEIA/1_extract_fragment_ends.py#L44-L50)). That reverses the forward-strand bases without complementing them, so a true 5′ `CCA` at the right end is reported as `GGT`. Even with the index bug already reported upstream fixed, `--both` would not give 5′ motifs at both ends. No benchmark change is needed, because the runner uses left ends. **Decision (2026-11-13):** no upstream report, because the project has not responded to the first one.

**Expected concordance difference (native).** The extractor clears its whole mate cache when it exceeds 10,000 reads (L90-91), which drops pending pairs. *Inference:* MAPQ 30 filtering leaves more reads whose mates never arrive, so the cache fills and clears more often than at FrEIA's default MAPQ 5. List this among the expected differences.

## cfDNAlab commands

Every cfDNAlab call uses the subcommand its label names, and every flag exists with the meaning the runners describe:

- `fcoverage --ignore-gap`, `--per-window total`, `--normalize-by-length` and `--by-size`.
- `lengths --indel-mode ignore` and `--clip-mode aligned`, both defaults passed explicitly.
- `ends --source-inside`, `--clip-strategy`, `--k-outside 0` and `--all-motifs`.
- `midpoints --length-bins`, `--bin-size`, `--blacklist` and the default `savgol=165` smoothing.

## Verified as presented

These need no action beyond existing disclosures:

- **mosdepth.** `--fragment-mode` and its proper-pair requirement.
- **samtools.** All `-@` options except `bedcov`.
- **FinaleToolkit.** All three commands. The end-motif release bugs are handled by using `main` at `d16632ec`.
- **cfDNAPro.** `readBam()` harmonization and `callMotif` type `s`.
- **FrEIA.** Preprocessing, extraction and transformation arguments, except `--subchroff`.
- **Griffin.** Every profiling and GC-estimation argument, apart from the items above.

## Implementation checklist

The checkmarks below record implemented code/documentation and prepared test updates, not executed tests or cluster validation. Every item in the first five groups is decided. The reasoning is in the finding sections above. Line numbers are from the working tree on 2026-11-13 and will shift as items are done. Update tests that assert the old behaviour as part of each item.

### Coverage

- [x] **Finding 2.** Remove `--input-fmt-option "nthreads=$tool_threads"` from `samtools bedcov` ([benchmark_coverage.sh L360](coverage/benchmark_coverage.sh#L360)). Record `tool_threads=NA` for `samtools-regions` (L356), and correct the comment at L353-355, which presents the option as decoding threads.
- [x] **Finding 7.** Drop `-aa` from `samtools depth` ([benchmark_coverage.sh L350](coverage/benchmark_coverage.sh#L350)) and from the comment at L338.
- [x] **Finding 7.** In [compare_coverage.py L86](coverage/compare_coverage.py#L86), treat `samtools-depth` as sparse like `cfdnalab*`, so its absent positions count as zero through `zero_filled_runs()`. Update the docstring at L80-81.
- [x] **Finding 7.** In the pass that counts `nonzero_union_bases`, also count bases covered by both tracks. For every positional pair, report the Jaccard index (both / either) genome-wide and per chromosome. Flag the truncation pattern: one chromosome with a partial Jaccard and every later chromosome, in BAM order, at 0. Report it as a diagnostic, not a pass/fail threshold.
- [x] **Finding 8.** Replace PanDepth's `-b "$BENCH_WORK/pandepth-regions.bed"` with `-w 500000` ([benchmark_coverage.sh L378-381](coverage/benchmark_coverage.sh#L378-L381)), point `track` at the `-w` output file, and remove the 1-based BED copy (L189-194).
- [x] **Finding 8.** In [compare_coverage.py L335-341](coverage/compare_coverage.py#L335-L341), read the `-w` output. The v2.26 writer confirms 1-based inclusive starts/ends and TotalDepth in column 6. This was checked from pinned source, not from a cluster output file. For PanDepth only, relax the exact row-identity check (L361-364): compare the windows it reports and record how many shared windows it lacks, because a 1 bp final window may be dropped.
- [x] **Finding 17.** Run `mosdepth-regions` with `MOSDEPTH_PRECISION=5` in its environment ([benchmark_coverage.sh L330-332](coverage/benchmark_coverage.sh#L330-L332)).
- [x] Docs: [coverage/README.md](coverage/README.md) (L12 `-aa`; L15 and L58 precision; L58 PanDepth BED conversion), [experiment_commands.md](experiment_commands.md) (L86-90 `-aa`; L143 and L185 PanDepth BED; L174 `nthreads=`; L188 precision) and [benchmarks_concordance_spec.md](benchmarks_concordance_spec.md) (L31 `-aa`; L35 `nthreads=`).
- [x] Tests: [test_benchmarks.py L527](tests/test_benchmarks.py#L527) asserts `nthreads=` in the `bedcov` command.

### Fragment lengths

- [x] **Findings 3, 12 and 19.** Time the documented Picard call alone: `picard -Xmx2g CollectInsertSizeMetrics I=sample.bam O=metrics.txt HISTOGRAM_WIDTH=$BENCH_MAX_LENGTH` ([benchmark_lengths.sh L272-292](lengths/benchmark_lengths.sh#L272-L292)). Drop the `samtools view` prefilter, `-Xms512m`, `-XX:ActiveProcessorCount`, `MINIMUM_PCT=0` and `INPUT=/dev/stdin`. `INCLUDE_DUPLICATES=false` is the default and can go. The temporary directory can stay, preferably as Picard's documented `TMP_DIR=`. Rewrite the comment block, which explains the pipeline and the old heap settings. `tool_threads` stays `NA`.
- [x] **Finding 3.** Produce the concordance histogram in an untimed step from the current `samtools view -u -q MAPQ -f 1 -F 3852 | picard` pipeline. [compare_lengths.py L94-95](lengths/compare_lengths.py#L94-L95) must read that histogram, not the timed output.
- [x] **Finding 14.** Replace `read_bam_insert_metrics(fragment_obj = ...)` with `callLength(frag_obj = fragments, isize_min = ..., isize_max = ...)` ([cfdnapro_lengths.R L47-53](lengths/cfdnapro_lengths.R#L47-L53)). Make [compare_lengths.py](lengths/compare_lengths.py) accept the `prop` column that `callLength()` adds when it zero-fills missing lengths.
- [x] **Finding 9.** Label the samtools lengths method "samtools stats (full QC report, insert-size section extracted)" in tables, figures and [lengths/README.md](lengths/README.md).
- [x] Docs: [lengths/README.md](lengths/README.md) (L9, L20 and L47 cfDNAPro call; L12 and L75 Picard pipeline), [experiment_commands.md](experiment_commands.md) (L282-288 Picard; L436 cfDNAPro) and [benchmarks_concordance_spec.md](benchmarks_concordance_spec.md) (L113 cfDNAPro; L117 Picard).

### End motifs

- [x] **Finding 10.** In [cfdnapro_end_motifs.R L27](ends/cfdnapro_end_motifs.R#L27), replace `chromosome_to_keep = FALSE` with the BAM's contigs that exist in the BSgenome selected by `genome_label`, computed from the two contig lists rather than by naming `chrEBV`. Log the dropped contigs. Record the `cfdnapro` scope as "BAM contigs present in cfDNAPro's reference" ([benchmark_end_motifs.sh L326-337](ends/benchmark_end_motifs.sh#L326-L337)). The lengths workload keeps `chromosome_to_keep = FALSE`.
- [x] **Finding 13.** Remove `--subchroff` ([benchmark_end_motifs.sh L322](ends/benchmark_end_motifs.sh#L322), [ends/README.md L12](ends/README.md#L12), [experiment_commands.md L387](experiment_commands.md#L387)). Wherever FrEIA is presented as a complete native workflow, relabel it "FrEIA rules filter_and_index_FrEIA, FragmEndSeq and data_transformation, run without Snakemake".
- [x] **Finding 15.** In [ends/README.md](ends/README.md), describe cfDNAPro's length filter as user-side R subsetting, because cfDNAPro has no length option for motifs. It stays inside the timed script.
- [x] Docs: record the cfDNAPro motif scope and the dropped contigs in [ends/README.md](ends/README.md) and [experiment_commands.md](experiment_commands.md).
- [x] **FinaleToolkit `main`.** Install FinaleToolkit from GitHub at the PR #194 merge commit `d16632ec5ae8999e774234457f6b0f9405c8f590` instead of the bioconda 1.1.0 package ([install_tools.sh L66 and L241-263](install_tools.sh#L241-L263)). Record the commit with the version in timing rows ([benchmark_end_motifs.sh L180](ends/benchmark_end_motifs.sh#L180), [benchmark_coverage.sh L162](coverage/benchmark_coverage.sh#L162), [benchmark_lengths.sh L125](lengths/benchmark_lengths.sh#L125)), because the package may still report 1.1.0. If a release containing the fix appears before the runs, use that release instead.
- [x] **FinaleToolkit `main`.** In [benchmark_end_motifs.sh L338-345](ends/benchmark_end_motifs.sh#L338-L345), record the requested length bounds as native and effective, replacing `length_filter=ignored-by-release` and the `NA` effective bounds, and update the comments at L115 and L340. Keep the guard for contigs whose length is an exact multiple of 1 Mb (L174), because the window builder is unchanged. Reword its message so it no longer names 1.1.0 only.
- [x] Docs: update FinaleToolkit's version, length filter and double-counting text in [ends/README.md](ends/README.md) (L11, L22-24, L122-123, L133 and the source links pinned to v1.1.0), [experiment_commands.md](experiment_commands.md) (L48 and L343), [benchmarks_concordance_spec.md L152](benchmarks_concordance_spec.md#L152), and the FinaleToolkit version in [lengths/README.md](lengths/README.md) (L10, L136 and L154). Add the remaining boundary loss for fragments longer than twice the read length to the end-motif differences, as already disclosed for FinaleToolkit `coverage`.

### Midpoints

- [x] **Finding 16.** Default `BENCH_MIN_MAPQ` to 20 for midpoints ([benchmark_midpoints.sh L32](midpoints/benchmark_midpoints.sh#L32)). That covers cfDNAlab, Griffin profiling and both GC estimation paths ([gc_corrected_profile.sh](midpoints/gc_corrected_profile.sh) L22, L27 and L45). The other features keep 30.
- [x] **Finding 4.** Remove the `griffin` method ([benchmark_midpoints.sh L67 and L92](midpoints/benchmark_midpoints.sh#L67)) and its run accounting in [summarize_timings.py L202](summarize_timings.py#L202). Keep the unit GC table (L306-308) only if the untimed Griffin demo still needs it. Preserve both Griffin tables, `sample.uncorrected.coverage.tsv` and `sample.GC_corrected.coverage.tsv`, from the first measured `griffin-gc` repetition.
- [x] **Finding 4.** In [compare_midpoints.py L137-143](midpoints/compare_midpoints.py#L137-L143), read the uncorrected Griffin table from `griffin-gc` for the uncorrected comparison with `cfdnalab`.
- [x] **Finding 5.** Pass Griffin's official `--save_window -1000 1000` ([griffin_profile.sh L49](midpoints/griffin_profile.sh#L49)). Keep `--norm_window` at the full window (L39).
- [x] **Finding 5.** In [compare_midpoints.py](midpoints/compare_midpoints.py), normalize cfDNAlab over all its bins as now (L207), then compare only the 132 bins whose signed start lies in `[-990, 990)`. Change the position-label check (L79-86) to expect exactly those labels. Remove the mean-one check on Griffin's profile (L110-115), which no longer holds after cropping. Report the compared bins and their width in `positions` and `window_bp`.
- [x] **Finding 16.** Describe `griffin_sites` as the number of sites in Griffin's site list, not an eligible count ([compare_midpoints.py L176](midpoints/compare_midpoints.py#L176), [benchmarks_concordance_spec.md L194](benchmarks_concordance_spec.md#L194), [midpoints/README.md L194](midpoints/README.md#L194)). cfDNAlab's `eligible_intervals` is named correctly.
- [x] Docs: [midpoints/README.md](midpoints/README.md) (L11 MAPQ; L14 saved window and its normalization rationale; L143, L164 and L168 `BENCH_METHODS=griffin`), [experiment_commands.md](experiment_commands.md) (L563-576 Griffin calls, MAPQ and unit GC setup), [benchmarks_concordance_spec.md](benchmarks_concordance_spec.md) (L13, L174 and L182-194), [benchmarks_spec.md](benchmarks_spec.md) (L128 MAPQ default; L148 `BENCH_METHODS=griffin`) and the Griffin compute budget in [README.md L111](README.md#L111).
- [x] Tests: [test_benchmarks.py L654-662](tests/test_benchmarks.py#L654-L662), [test_review_regressions.py](tests/test_review_regressions.py) (L906-1028 and L1288-1435 use the `griffin` method and directory) and [test_midpoint_gc.py L198](tests/test_midpoint_gc.py#L198) (`--save_window`).

- [ ] Confirm the PanDepth fixed-window reader against the first real cluster output. No cluster output was available during implementation.

### Disclosures for the methods text

These disclosures are now included in the benchmark README files, specifications and command report. Transferring them to the manuscript remains the authors' task.

- [ ] **Finding 1.** mosdepth runs in exact mode, although its help recommends `--fast-mode`, because exact mode computes the same quantity as cfDNAlab: mate overlap counted once and deletions excluded.
- [ ] **Finding 3.** Timed Picard keeps MAPQ <30 and QC-fail pairs, because it has no MAPQ option. Its concordance histogram comes from an untimed, prefiltered run.
- [ ] **Finding 10.** cfDNAPro motifs exclude contigs missing from its reference, which is `chrEBV` for this BAM.
- [ ] **Finding 11.** samtools applies the length and MAPQ filters through its global htslib input option. PanDepth has no equivalent, so it runs without a length filter.
- [ ] **Findings 16 and 19.** Griffin's plotting step and Picard's PDF are not run, while cfDNAlab's default QC plot is timed.
- [ ] **PanDepth per-base.** Its time includes a single-threaded gzip of every base, which runs after its worker threads finish. PanDepth has no uncompressed output, so its per-base time scales poorly with cores by design. Also add this sentence to [coverage/README.md](coverage/README.md).
- [ ] **FinaleToolkit.** End motifs use unreleased `main` at `d16632ec`, which contains the maintainers' fixes for the ignored length bounds and the double counting ([PR #194](https://github.com/epifluidlab/FinaleToolkit/pull/194), [issue #192](https://github.com/epifluidlab/FinaleToolkit/issues/192)).

### Dropped or deferred

Do not implement these.

- **FrEIA mate-cache clearing.** Dropped. It is how FrEIA's code works, and each clear loses only the few pairs still waiting for their mate.
- **Finding 18 upstream report.** Dropped. The project has not responded to the first report.
- **Midpoint smoothing diagnostic.** Deferred. Use it only if cfDNAlab's and Griffin's profiles differ a lot. Run cfDNAlab untimed with `--smoothing none --bin-size 15`, then apply Griffin's `savgol_filter(profile, 11, 3)` and mean-one normalization in the comparator. If that matches Griffin much better than default cfDNAlab, the difference comes from the order of smoothing and binning, not from midpoint counting. Without smoothing, cfDNAlab's exclusion margin shrinks from 182 to 100 bp, so a few more sites may be kept. Compare site counts too.
- [ ] Finding 18: optionally add the right-end complement bug to the upstream FrEIA issue. This is for the authors, not Codex.
