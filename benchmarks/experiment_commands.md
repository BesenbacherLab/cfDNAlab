# Benchmark experiment commands

Each measured repetition runs the native commands below on the same staged BAM. Commands are shown once, with variables for the allocated cores, paths and experimental settings. Executable aliases resolve to the staged paths used by the runners. Shell argument arrays are written as ordinary command lines, and timed shell pipelines are expanded into their constituent commands. The executable definitions are in each feature's `commands.sh` ([coverage](coverage/commands.sh), [lengths](lengths/commands.sh), [end motifs](ends/commands.sh), [midpoints](midpoints/commands.sh)), which uses the same variable names. Every repetition runs its command through [measure_command.sh](measure_command.sh) in its own Slurm step, and `runs/<run_order>/command.sh` in the results records the exact command. For `cfdnalab-gc` and `griffin-gc`, which run a helper script, `runs/<run_order>/command-trace.log` also records each native call the helper made, with its variables filled in.

The configurations describe the current runners. For an existing result, `run.tsv` identifies the benchmark revision and `timings.tsv` records the tool version, filters and core allocation. Installation and cfDNAlab build records identify the installed tools independently of the benchmark revision.

## Shared settings and timing

| Setting | Coverage and lengths | End motifs | Midpoints |
| --- | --- | --- | --- |
| Allocated cores | 1, 4, 12 or 32, in separate jobs | Same | Same |
| Minimum MAPQ | 30, except native timed Picard | 30 | 20 for profiling and GC estimation |
| Included fragment lengths | 30-1000 bp inclusive | Requested 30-1000 bp inclusive | 100-200 bp inclusive for profiling |
| Chromosome selection | All BAM contigs, except named autosome configurations | All BAM contigs, except cfDNAPro reference compatibility and named autosome configuration | Chromosomes in the supplied site set |
| Production repetitions | 1 warmup, then 3 measurements | Same | Same for cfDNAlab. Griffin uses an untimed 1000-site demo, then 3 measurements |

The variables in the calls correspond to the job settings:

```bash
cores=$BENCH_CORES
min_mapq=$BENCH_MIN_MAPQ
min_length=$BENCH_MIN_LENGTH
max_length=$BENCH_MAX_LENGTH
length_end=$((max_length + 1))

bam="$BENCH_WORK/sample.bam"
output="$BENCH_WORK/output"
regions="$BENCH_WORK/regions.bed"
reference_2bit="$BENCH_WORK/reference.2bit"
genome=$BENCH_GENOME
autosomes=...  # Comma-separated BAM contigs named chr1..chr22 or 1..22

cfdna="$BENCH_WORK/bin/cfdna"
mosdepth="$BENCH_WORK/bin/mosdepth"
pandepth="$BENCH_WORK/bin/pandepth"
samtools="$BENCH_WORK/samtools-env/bin/samtools"
finale_python="$BENCH_WORK/finale-env/bin/python"
finale="$BENCH_WORK/finale-env/bin/finaletoolkit"
picard="$BENCH_WORK/picard-env/bin/picard"
cfdnapro_image="$BENCH_WORK/cfdnapro.sif"
apptainer=$(command -v apptainer)
```

`BENCH_WORK` is the job's temporary directory under the compute node's `$TMPDIR`. The runners copy the BAM, index, selected tools and required references there before timing. `output` is created empty for every repetition. Container paths under `/work` refer to the same directory through an Apptainer bind mount.

Every available native thread/core/worker option receives `cores` directly. Commands without a native thread option still run within the allocated CPU set. This also applies to samtools inside the Picard pipeline.

The native options have different meanings. For example, the [mosdepth CLI](https://github.com/brentp/mosdepth/blob/v0.3.14/README.md#usage) defines `--threads` as BAM decompression threads, while coverage calculation runs separately. A 4-core allocation uses `--threads 4`, and all of the tool's threads share those allocated CPUs. `timings.tsv` records the core allocation and the native thread/worker setting separately.

The installer requests mosdepth 0.3.14, samtools 1.25, PanDepth 2.26, FinaleToolkit commit `d16632ec5ae8999e774234457f6b0f9405c8f590`, cfDNAPro 1.7.4, Picard 3.5.0 and Griffin v0.2.0. cfDNAlab is built from the checkout recorded in `cfdnalab-build.tsv`.

The exclusion mask `3852`, or `0xF0C`, combines unmapped (`4`), mate unmapped (`8`), secondary (`256`), QC-fail (`512`), duplicate (`1024`) and supplementary (`2048`) flags. It does not require the proper-pair flag.

GNU time surrounds each complete native command or shell pipeline. `wall_seconds` includes startup, BAM reading, native calculations, output writing/compression and enabled default plots. The full R process and Apptainer startup are included for cfDNAPro. Pipeline shells use `set -euo pipefail`.

After native success, a separate timer measures `sync -f` on the output filesystem. `process_plus_flush_seconds` adds that duration to `wall_seconds`. Staging, Slurm startup, CPU-placement checks, saved-profile copying, cleanup and concordance are outside the native command timer. Per-command timeouts are disabled when `BENCH_RUN_TIMEOUT_SECONDS=0`, the default. Slurm's allocation limit applies to the whole job.

## Coverage

The [coverage runner](coverage/benchmark_coverage.sh) contains three comparisons: positional read depth, positional fragment coverage, and regional summaries. Regional read depth and regional fragment counts are compared separately.

### Read depth at individual bases

**cfDNAlab, method `cfdnalab`:**

```bash
"$cfdna" fcoverage \
    --bam "$bam" --output-dir "$output" \
    --min-mapq "$min_mapq" \
    --min-fragment-length "$min_length" --max-fragment-length "$max_length" \
    --n-threads "$cores" --chromosomes all --ignore-gap
```

**mosdepth, method `mosdepth`:**

```bash
"$mosdepth" --threads "$cores" \
    --include-flag 1 --flag 3852 --mapq "$min_mapq" \
    --min-frag-len "$min_length" --max-frag-len "$max_length" \
    "$output/coverage" "$bam"
```

**samtools, method `samtools-depth`:**

```bash
"$samtools" depth -@ "$cores" \
    --input-fmt-option "filter=(tlen >= $min_length && tlen <= $max_length) || (tlen <= -$min_length && tlen >= -$max_length)" \
    -s -q 0 -Q "$min_mapq" -G 3852 --require-flags 1 \
    -o "$output/samtools.depth.tsv" "$bam"
```

Here `-q 0` is the base-quality threshold and `-Q` is MAPQ. `-s` removes mate overlap. Native sparse output omits zero-depth positions, which the comparator treats as zero. samtools' global HTSlib input expression applies inclusive absolute TLEN bounds to both signs during BAM reading.

**PanDepth, method `pandepth`:**

```bash
"$pandepth" -i "$bam" -o "$output/pandepth" \
    -q "$min_mapq" -x 3852 -t "$cores" -a
```

PanDepth counts overlapping mates twice and has no fragment length filter. cfDNAlab merges mate overlap. mosdepth's native overlap correction depends on proper-pair status. The [concordance specification](benchmarks_concordance_spec.md#native-differences) describes the remaining native pairing and CIGAR differences.

### Fragment coverage at individual bases

**cfDNAlab, method `cfdnalab-fragment`:**

```bash
"$cfdna" fcoverage \
    --bam "$bam" --output-dir "$output" \
    --min-mapq "$min_mapq" \
    --min-fragment-length "$min_length" --max-fragment-length "$max_length" \
    --n-threads "$cores" --chromosomes all
```

**mosdepth, method `mosdepth-fragment`:**

```bash
"$mosdepth" --threads "$cores" \
    --include-flag 1 --flag 3852 --mapq "$min_mapq" \
    --min-frag-len "$min_length" --max-frag-len "$max_length" \
    --fragment-mode "$output/coverage" "$bam"
```

Both include the inter-mate gap. mosdepth fragment mode requires proper pairs and fills the complete TLEN span. cfDNAlab uses directional paired spans and excludes D/N positions unsupported by the other mate.

### Shared regional windows

Before timing, the runner creates disjoint 500 kb BED4 windows over every BAM contig, with shorter final windows:

```bash
"$samtools" idxstats "$bam" |
    awk 'BEGIN {OFS="\t"} $1 != "*" {
        for (window_start=0; window_start<$2; window_start+=500000) {
            window_end=window_start+500000
            if (window_end>$2) window_end=$2
            print $1, window_start, window_end, "."
        }
    }' > "$regions"
```

cfDNAlab, mosdepth and PanDepth generate 500 kb windows with their native fixed-window options. samtools and FinaleToolkit read the shared BED. PanDepth coordinates are converted only when its native output is read for concordance.

### Regional read depth

These four methods measure summed read depth, or covered-base totals, in each 500 kb window. No fragment length normalization is applied.

**cfDNAlab, method `cfdnalab-regions`:**

```bash
"$cfdna" fcoverage \
    --bam "$bam" --output-dir "$output" \
    --min-mapq "$min_mapq" \
    --min-fragment-length "$min_length" --max-fragment-length "$max_length" \
    --n-threads "$cores" --chromosomes all --ignore-gap \
    --by-size 500000 --per-window total
```

**mosdepth, method `mosdepth-regions`:**

```bash
env MOSDEPTH_PRECISION=5 "$mosdepth" --threads "$cores" \
    --include-flag 1 --flag 3852 --mapq "$min_mapq" \
    --min-frag-len "$min_length" --max-frag-len "$max_length" \
    --by 500000 --no-per-base "$output/coverage" "$bam"
```

**samtools, method `samtools-regions`:**

```bash
"$samtools" bedcov \
    --input-fmt-option "filter=flag.paired && ((tlen >= $min_length && tlen <= $max_length) || (tlen <= -$min_length && tlen >= -$max_length))" \
    -j -Q "$min_mapq" -G 3852 "$regions" "$bam" \
    > "$output/samtools.bedcov.tsv"
```

**PanDepth, method `pandepth-regions`:**

```bash
"$pandepth" -i "$bam" -o "$output/pandepth" \
    -q "$min_mapq" -x 3852 -t "$cores" \
    -w 500000
```

These commands write regional output directly. mosdepth uses its documented `MOSDEPTH_PRECISION=5`, and concordance multiplies each mean by its actual window width. PanDepth writes `pandepth.win.stat.gz` with 1-based inclusive coordinates and `TotalDepth`. Its missing final 1 bp windows are counted in the comparison metadata. samtools bedcov has no native thread option (`tool_threads=NA`). Its paired/TLEN filter is samtools' global input option. `-j` excludes D/N positions but does not remove mate overlap. PanDepth also counts overlapping mates twice and has no length filter.

### Regional fragment counts

**cfDNAlab, method `cfdnalab-region-counts`:**

```bash
"$cfdna" fcoverage \
    --bam "$bam" --output-dir "$output" \
    --chromosomes all --by-size 500000 --per-window total --normalize-by-length \
    --min-fragment-length "$min_length" --max-fragment-length "$max_length" \
    --min-mapq "$min_mapq" --n-threads "$cores"
```

**FinaleToolkit, method `finaletoolkit-coverage`:**

```bash
"$finale_python" "$finale" coverage "$bam" "$regions" \
    --output "$output/finaletoolkit.coverage.bed.gz" \
    --min-length "$min_length" --max-length "$max_length" \
    --min-mapq "$min_mapq" --threads "$cores"
```

cfDNAlab gives each fragment total weight 1 before output rounding, split across windows when it crosses a boundary. FinaleToolkit uses its default midpoint assignment and requires proper read1/TLEN selection. Its indexed read1 query can miss fragments whose midpoint is inside the window while read1 lies outside. These counts are compared separately from regional read-depth totals. The end-motif length-filter bug described below does not apply to this coverage call.

### Additional coverage configurations

| Method | Exact change to the corresponding call above |
| --- | --- |
| `cfdnalab-autosomes` | Read-depth call, replacing `--chromosomes all` with `--chromosomes "$autosomes"` |
| `cfdnalab-fragment-proper` | Fragment-coverage call plus `--require-proper-pair` |
| `cfdnalab-region-counts-proper` | Normalized regional-count call plus `--require-proper-pair` |

`autosomes` is the comma-separated set of BAM contigs named `chr1` through `chr22`, or `1` through `22`.

## Fragment lengths

The [lengths runner](lengths/benchmark_lengths.sh) measures global histograms with one bin for each integer fragment length from `min_length` through `max_length`, inclusive.

### cfDNAlab

Method `cfdnalab`:

```bash
"$cfdna" lengths --bam "$bam" --output-dir "$output" \
    --length-bins "$min_length:$length_end:1" \
    --indel-mode ignore --clip-mode aligned --min-mapq "$min_mapq" \
    --n-threads "$cores" --chromosomes all
```

### cfDNAPro

Method `cfdnapro`, with the final argument selecting native coordinate curation:

```bash
"$apptainer" exec --cleanenv --containall \
    --bind "$BENCH_WORK:/work" --pwd /work \
    "$cfdnapro_image" Rscript --vanilla /work/cfdnapro_lengths.R \
    /work/sample.bam /work/output/lengths.tsv \
    "$min_length" "$max_length" "$min_mapq" "$genome" TRUE
```

The native R calculation is shown under [cfDNAPro BAM reading and calculations](#cfdnapro-bam-reading-and-calculations). This function has no native thread option. The allocated CPU set still applies to the complete container/R process.

### FinaleToolkit

Method `finaletoolkit`:

```bash
"$finale_python" "$finale" frag-length-bins "$bam" \
    --bin-size 1 --min-length "$min_length" --max-length "$max_length" \
    --min-mapq "$min_mapq" --output "$output/lengths.tsv"
```

This command has no native worker-count option. Its length filtering is separate from FinaleToolkit's end-motif bug.

### samtools

Method `samtools`:

```bash
"$samtools" stats -@ "$cores" \
    --input-fmt-option "filter=mapq >= $min_mapq" \
    -f 1 -F 3852 -m 1 -i "$length_end" "$bam" > "$output/stats.txt"
```

The timer includes the complete native `stats` command. The comparator later selects its inward-pair `IS` column and the requested length range. `-m 1` prevents tail truncation, and the overflow bin is placed above the requested maximum.

### Picard

Method `picard` times the standalone native command:

```bash
"$picard" -Xmx2g CollectInsertSizeMetrics \
    "I=$bam" "O=$output/metrics.txt" "HISTOGRAM_WIDTH=$max_length" \
    "TMP_DIR=$BENCH_WORK/tmp"
```

Picard has no MAPQ option, so this timed call retains low-MAPQ and QC-fail pairs. Its timing rows record `min_mapq=NA`. The 2 GB heap is Picard's documented recommendation. No PDF is requested because picard-slim lacks R, while cfDNAlab's default plot stays timed.

Once after all timed measurements, the following untimed pipeline in [picard_filtered_histogram.sh](lengths/picard_filtered_histogram.sh) supplies the common-filter concordance histogram. It runs in its own Slurm step with the same CPU binding, adds job wall time and has no benchmark timing row:

```bash
"$samtools" view -@ "$cores" -u -q "$min_mapq" -f 1 -F 3852 "$bam" |
    "$picard" -Xmx2g CollectInsertSizeMetrics I=/dev/stdin \
        "O=$BENCH_RESULTS/profiles/picard/concordance.metrics.txt" \
        "HISTOGRAM_WIDTH=$max_length" "TMP_DIR=$BENCH_WORK/tmp"
```

The comparator reads only `concordance.metrics.txt`, selects the all-reads FR histogram, and limits it to the requested length range. The timed `metrics.txt` is retained separately.

### Additional length configurations

| Method | Exact change |
| --- | --- |
| `cfdnalab-autosomes` | cfDNAlab call with `--chromosomes "$autosomes"` instead of `all` |
| `cfdnalab-proper` | cfDNAlab all-contig call plus `--require-proper-pair` |
| `cfdnapro-outer-span` | cfDNAPro container call with its final argument changed from `TRUE` to `FALSE` |

cfDNAlab uses directional paired spans. Curated cfDNAPro also uses directional spans, while its uncurated diagnostic uses outer alignment spans. FinaleToolkit, samtools and Picard use native TLEN-based selection and can retain different fragment populations.

## End motifs

The [end-motif runner](ends/benchmark_end_motifs.sh) uses `k=$BENCH_K`. The four-mer comparison selects `cfdnalab-reference`, `cfdnalab-reference-proper`, `finaletoolkit` and `cfdnapro`. Every selected tool can be allocated each requested core setting, including tools without a native thread option.

### cfDNAlab reference motifs

Method `cfdnalab-reference`:

```bash
"$cfdna" ends --bam "$bam" --output-dir "$output" \
    --k-inside "$k" --k-outside 0 --all-motifs \
    --min-mapq "$min_mapq" \
    --min-fragment-length "$min_length" --max-fragment-length "$max_length" \
    --n-threads "$cores" --source-inside reference --clip-strategy aligned \
    --ref-2bit "$reference_2bit" --chromosomes all
```

Method `cfdnalab-reference-proper` runs this same call with `--require-proper-pair` appended.

### cfDNAlab aligned read motifs

Method `cfdnalab-read-aligned`:

```bash
"$cfdna" ends --bam "$bam" --output-dir "$output" \
    --k-inside "$k" --k-outside 0 --all-motifs \
    --min-mapq "$min_mapq" \
    --min-fragment-length "$min_length" --max-fragment-length "$max_length" \
    --n-threads "$cores" --source-inside read --clip-strategy aligned \
    --chromosomes all
```

### FinaleToolkit

Method `finaletoolkit`:

```bash
"$finale_python" "$finale" end-motifs "$bam" "$reference_2bit" -k "$k" \
    --min-length "$min_length" --max-length "$max_length" \
    --min-mapq "$min_mapq" --threads "$cores" --strand both \
    --output "$output/motifs.tsv"
```

FinaleToolkit uses the upstream merge commit `d16632ec5ae8999e774234457f6b0f9405c8f590`, which applies inclusive length bounds and unique midpoint window membership. Requested and effective bounds match (`length_filter=applied`). A fragment longer than twice its read length can still be lost when read1 does not overlap the window containing its midpoint. The unchanged 1 Mb window-builder omission is guarded before timing. Native proper read1/TLEN selection remains.

### cfDNAPro

Method `cfdnapro`:

```bash
"$apptainer" exec --cleanenv --containall \
    --bind "$BENCH_WORK:/work" --pwd /work \
    "$cfdnapro_image" Rscript --vanilla /work/cfdnapro_end_motifs.R \
    /work/sample.bam /work/output \
    "$min_length" "$max_length" "$min_mapq" "$genome" "$k"
```

The native R calculation below selects reference-derived left ends. Reference loading, BAM reading, curation, motif calling and output/provenance writing are timed. The package's bundled reference is selected by `BENCH_GENOME`, normally `hg38`. Sequence identity with the supplied alignment reference remains to be established before interpreting differences solely as tool differences.

### Additional motif configurations

| Method | Exact change |
| --- | --- |
| `cfdnalab-read` | Aligned-read call with `--clip-strategy skip` |
| `cfdnalab-read-include` | Aligned-read call with `--clip-strategy include-at-aligned-boundary` |
| `cfdnalab-reference-autosomes` | Reference call with `--chromosomes "$autosomes"` instead of `all` |

These configurations are only run when explicitly selected in `BENCH_METHODS`.

## cfDNAPro BAM reading and calculations

The complete R helpers are [cfdnapro_lengths.R](lengths/cfdnapro_lengths.R) and [cfdnapro_end_motifs.R](ends/cfdnapro_end_motifs.R). Their computational calls are expanded here so the container invocation does not hide the filtering or calculation. Argument parsing and error checks remain in the linked helpers.

The length helper reads the BAM with the following options. The motif helper uses the same call with `chromosome_to_keep = retained_contigs`, computed from the BAM and reference headers as shown below. `bam_path`, `genome_label`, `minimum_mapq`, `minimum_length` and `maximum_length` are their positional command-line arguments. `curate_coordinates` is `TRUE` for normal lengths and motifs, and `FALSE` only for `cfdnapro-outer-span`.

```r
suppressPackageStartupMessages(library(cfDNAPro))

fragments <- readBam(
    bamfile = bam_path,
    genome_label = genome_label,
    chromosome_to_keep = FALSE,
    curate_start_and_end = curate_coordinates,
    galp_mapqFilter = minimum_mapq,
    galp_flag = Rsamtools::scanBamFlag(
        isPaired = TRUE,
        isUnmappedQuery = FALSE,
        hasUnmappedMate = FALSE,
        isSecondaryAlignment = FALSE,
        isSupplementaryAlignment = FALSE,
        isDuplicate = FALSE,
        isNotPassingQualityControls = FALSE
    ),
    outdir = NA
)
```

The **length** calculation is:

```r
counts <- callLength(
    frag_obj = fragments,
    isize_min = minimum_length,
    isize_max = maximum_length
)
write.table(counts, output_path, sep = "\t", quote = FALSE, row.names = FALSE)
```

The **motif** calculation applies user-side R subsetting for inclusive length bounds to the curated fragments, then calls the native left-end function. `motif_length` is the supplied `k` argument:

```r
fragment_lengths <- BiocGenerics::width(fragments)
fragments <- fragments[fragment_lengths >= minimum_length & fragment_lengths <= maximum_length]
counts <- callMotif(
    frag_obj = fragments,
    motif_type = "s",
    motif_length = motif_length,
    integrate_mut = FALSE,
    downsample_ref = FALSE
)
write.table(counts, file.path(output_directory, "cfdnapro-motifs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
```

The BAM-reading call above uses `chromosome_to_keep = FALSE` for lengths. For motifs it instead uses `chromosome_to_keep = retained_contigs`, computed inside the timed R process:

```r
reference_package <- switch(
    genome_label,
    hg19 = "BSgenome.Hsapiens.UCSC.hg19",
    hg38 = "BSgenome.Hsapiens.UCSC.hg38",
    `hg38-NCBI` = "BSgenome.Hsapiens.NCBI.GRCh38"
)
reference_contigs <- GenomeInfoDb::seqlevels(BSgenome::getBSgenome(reference_package))
bam_contigs <- names(Rsamtools::scanBamHeader(bam_path)[[1L]]$targets)
retained_contigs <- intersect(bam_contigs, reference_contigs)
dropped_contigs <- setdiff(bam_contigs, retained_contigs)
```

The dropped contigs are logged and saved in `contig_scope.tsv`. The motif timing scope is `bam-contigs-in-cfdnapro-reference`, which excludes chrEBV for this BAM. Other tools keep whole-BAM scope. The length benchmark keeps every BAM contig.

The motif helper also writes `reference.tsv` from the native fragment object's reference metadata. Both reads and writes occur inside the R process timer.

## Midpoint profiles

The [midpoint runner](midpoints/benchmark_midpoints.sh) uses the supplied grouped 9,990 bp windows (`center - 4995` to `center + 4995`) and profiling length bounds 100-200 bp. `BENCH_CHROMOSOMES` is derived from the chromosomes in those windows. The runner preserves the sites, including repeated and overlapping rows, and creates the native Griffin anchor tables outside timing with [prepare_sites.py](midpoints/prepare_sites.py). Each workload runs in its own job. cfDNAlab retains 666 bins of 15 bp, with native 165 bp smoothing. Griffin normalizes over that full window and saves its official central window. Concordance normalizes cfDNAlab over all 666 bins before selecting the same 132 central bins. The following shared exclusion paths name the staged release BEDs:

```bash
exclusions=(
    "$BENCH_WORK/exclusions/encode_unified_GRCh38_exclusion_list.bed"
    "$BENCH_WORK/exclusions/hg38_centromeres.bed"
    "$BENCH_WORK/exclusions/hg38_gaps.bed"
    "$BENCH_WORK/exclusions/hg38_fix_patches.bed"
    "$BENCH_WORK/exclusions/hg38_alternative_haplotypes.bed"
)
```

### Site preparation

The supplied windows are prepared once, before any job, from Griffin's [paper TFBS set](https://github.com/adoebley/Griffin_analyses/tree/75467e438a88029206d8d778ddce05eb0ab6efa8/sites/TFBS/10000_unfiltered_sites_CIS_BP_v2) of 377 transcription factors with up to 10,000 sites each. The repository's `main` branch has pointed to that pinned revision since December 2022.

```bash
git clone --depth 1 --filter=blob:none --sparse --branch main \
    https://github.com/adoebley/Griffin_analyses.git griffin-sites
site_folder=sites/TFBS/10000_unfiltered_sites_CIS_BP_v2
git -C griffin-sites sparse-checkout set "$site_folder"

awk -F '\t' 'BEGIN {OFS="\t"} FNR>1 {print $1,$2,$3,$5}' \
    "griffin-sites/$site_folder/"*.txt |
    LC_ALL=C sort -t $'\t' -k1,1 -k2,2n -k3,3n \
    > griffin_tfbs.bed
```

This keeps each site's `Chrom`, `Start`, `End` and `Gene` columns. The site files have no strand column, so every window is unstranded. Group names are then cleaned by removing the `.hg38` suffix and all non-alphanumeric characters, so `NKX2-1.hg38` becomes `NKX21`. The original cleaning command was not kept. This reconstruction applies the same rule and gives no name collisions across the 377 site lists:

```bash
awk 'BEGIN {FS=OFS="\t"} {sub(/\.[Hh][Gg]38$/, "", $4); gsub(/[^A-Za-z0-9]/, "", $4); print}' \
    griffin_tfbs.bed > griffin_tfbs.cleaned.bed
```

The final windows are 9,990 bp long:

```bash
cfdna prep-windows -i griffin_tfbs.cleaned.bed -o griffin_tfbs.prepared_9990bp.bed \
    --sep tab --group-cols 3 \
    --blacklist hg38-blacklist.v2.bed --blacklist k100.umap.exclusion_intervals.bed \
    --blacklist-halo 200 --resize 9990 --chrom-sizes hg38.chrom.sizes
```

The inputs are:

- `hg38-blacklist.v2.bed` is the ENCODE Blacklist v2 file `lists/hg38-blacklist.v2.bed.gz` from [Boyle-Lab/Blacklist](https://github.com/Boyle-Lab/Blacklist), decompressed.
- `k100.umap.exclusion_intervals.bed` was made by the benchmark author from the [Umap](https://bismap.hoffmanlab.org/) k100 hg38 mappability intervals.
- `hg38.chrom.sizes` is a standard hg38 chromosome sizes file. The primary chromosomes have the same lengths in all GRCh38 builds, so the file's source does not change the prepared windows.

`--resize` centers each window on `start + floor(length / 2)`, which matches Griffin's `position` column. For odd-length site intervals, the two equally centered placements are chosen by a deterministic hash, while Griffin always uses the placement starting at `position - 4995`. Windows within 200 bp of the two blacklists are removed before benchmarking. These lists are separate from the shared exclusion BEDs above, which both tools apply natively. Both tools receive the identical prepared BED as `BENCH_INTERVALS`.

### cfDNAlab without GC correction

Method `cfdnalab` supplies cfDNAlab's uncorrected profiles, which are compared with the uncorrected table that `griffin-gc` also writes. With a GC file, cfDNAlab writes only corrected profiles. Its timing has no Griffin counterpart:

```bash
"$cfdna" midpoints --bam "$bam" --output-dir "$output" \
    --intervals "$BENCH_WORK/sites.bed" --chromosomes "$BENCH_CHROMOSOMES" \
    --length-bins "$min_length" "$length_end" \
    --min-mapq "$min_mapq" --bin-size 15 --n-threads "$cores" \
    --blacklist "${exclusions[@]}"
```

### Griffin profiling calls

Method `griffin-gc` times GC estimation followed by both native profiling scripts in [griffin_profile.sh](midpoints/griffin_profile.sh). The following variables express that helper's production arguments:

```bash
griffin_python="$BENCH_WORK/griffin-env/bin/python"
griffin_scripts="$BENCH_WORK/griffin-env/share/griffin"
sites_yaml="$BENCH_WORK/sites.yaml"
chromosome_sizes="$BENCH_WORK/chrom.sizes"
reference_fasta="$BENCH_WORK/reference.fa"
half_width=$BENCH_HALF_WIDTH
window_end=$((BENCH_WINDOW_WIDTH - half_width))
gc_bias="$output/gc-bias/GC_bias/sample.GC_bias.txt"
IFS=, read -r -a chromosomes <<< "$BENCH_CHROMOSOMES"
```

The 9,990 bp windows use normalization and saved bounds `-4995 4995`. Ends are exclusive. Both bounds are multiples of 15, so Griffin's native inward rounding leaves them unchanged and the two tools' bins align. The runner does not resize the supplied BED. Auxiliary feature windows follow the official configuration. The arguments shared by the two scripts are:

```bash
griffin_common=(
    --sample_name sample --tmp_dir "$output/tmp"
    --sites_yaml "$sites_yaml" --griffin_scripts_dir "$griffin_scripts"
    --chrom_sizes_path "$chromosome_sizes" --chroms "${chromosomes[@]}"
    --chrom_column Chrom --position_column position --strand_column Strand
    --norm_window "-$half_width" "$window_end"
    --mappability_correction False --mappability_bw "$BENCH_WORK/mappability.bw" --CPU "$cores"
)
```

The two timed native calls are:

```bash
"$griffin_python" "$griffin_scripts/griffin_coverage.py" "${griffin_common[@]}" \
    --bam "$bam" --GC_bias "$gc_bias" \
    --reference_genome "$reference_fasta" --mappability_bias none \
    --size_range "$min_length" "$max_length" --map_quality "$min_mapq"

"$griffin_python" "$griffin_scripts/griffin_merge_sites.py" "${griffin_common[@]}" \
    --uncorrected_bw_path "$output/tmp/sample/tmp_bigWig/sample.uncorrected.bw" \
    --GC_corrected_bw_path "$output/tmp/sample/tmp_bigWig/sample.GC_corrected.bw" \
    --GC_map_corrected_bw_path none --results_dir "$output/results" \
    --save_window -1000 1000 \
    --center_window -30 30 --fft_window -960 960 --fft_index 10 \
    --step 15 --exclude_paths "${exclusions[@]}" --exclude_zero_mappability True \
    --exclude_outliers True --CNA_normalization False --individual False \
    --smoothing True --smoothing_length 165
```

The `griffin-gc` workload retains both `sample.uncorrected.coverage.tsv` and `sample.GC_corrected.coverage.tsv`. There is no separate timed `griffin` method. Corrected timings compare `griffin-gc` with `cfdnalab-gc`, and uncorrected concordance compares `cfdnalab` with the raw table from `griffin-gc`. Griffin normalizes before saving `[-1000, 1000)`, rounded to `[-990, 990)`. The comparator normalizes cfDNAlab over its complete 9,990 bp output, then selects the same 132 bins. Griffin's cropped profile need not have mean one. No per-site CNA correction or fragment-count weights are added. `cfdnalab_sites` counts eligible intervals, while `griffin_sites` counts sites in its input site list, including masked sites.

### cfDNAlab with GC correction

Method `cfdnalab-gc` times the following two calls together:

```bash
bias_directory="$output/gc-bias"

"$cfdna" gc-bias --bam "$bam" --output-dir "$bias_directory" \
    --ref-2bit "$reference_2bit" \
    --ref-gc-file "$BENCH_WORK/reference-gc.zarr" \
    --min-mapq "$min_mapq" --n-threads "$cores"

"$cfdna" midpoints --bam "$bam" --output-dir "$output" \
    --intervals "$BENCH_WORK/sites.bed" --chromosomes "$BENCH_CHROMOSOMES" \
    --length-bins "$min_length" "$length_end" \
    --min-mapq "$min_mapq" --bin-size 15 --n-threads "$cores" \
    --blacklist "${exclusions[@]}" \
    --ref-2bit "$reference_2bit" \
    --gc-file "$bias_directory/gc_bias_correction.zarr"
```

The reference GC package is prepared once beforehand and supplied through `BENCH_CFDNA_REF_GC`. The runner stages it as `reference-gc.zarr`. It does not calculate reference bias during repetitions. Sample bias is recalculated for every repetition with native settings, including its native windows/chromosome selection and the reference package's default 30-1000 bp range. Profiling still uses 100-200 bp. Native default diagnostic plots remain timed.

### Griffin with GC correction

Method `griffin-gc` first runs these two native estimation calls, then the same two profiling calls shown above. All four native calls share one timer:

```bash
bias_directory="$output/gc-bias"

"$griffin_python" "$griffin_scripts/griffin_GC_counts.py" \
    --bam_file "$bam" --bam_file_name sample \
    --mappable_regions_path "$BENCH_WORK/gc-regions.bed" \
    --ref_seq "$reference_fasta" --chrom_sizes "$BENCH_WORK/bam.chrom.sizes" \
    --out_dir "$bias_directory" --map_q "$min_mapq" \
    --size_range 15 500 --CPU "$cores"

"$griffin_python" "$griffin_scripts/griffin_GC_bias.py" \
    --bam_file_name sample --mappable_name "$BENCH_GC_REFERENCE_NAME" \
    --genome_GC_frequency "$BENCH_WORK/reference-gc" \
    --out_dir "$bias_directory" --size_range 15 500

gc_bias="$bias_directory/GC_bias/sample.GC_bias.txt"
bash "$BENCH_WORK/griffin_profile.sh" "$cores" "$sites_yaml" \
    "$chromosome_sizes" "$BENCH_CHROMOSOMES" "$output" "$gc_bias"
```

The final helper invocation expands to the two Griffin profiling calls above, with the newly fitted `gc_bias`. Estimation uses Griffin's native 15-500 bp range. Profiling uses the supplied 100-200 bp range. MAPQ is 20 for both tools and both estimation paths. The estimator's region BED and reference GC frequencies come from Griffin's bundled inputs, independently of cfDNAlab's reference package. `BENCH_GC_REFERENCE_NAME` is the estimation BED's basename without its final extension. The retained profile is `sample.GC_corrected.coverage.tsv`.

### Griffin demo warmup outside timing

Before the first measured Griffin repetition, both corrected and uncorrected jobs run the release's 1000-site CTCF demo through the same profiling helper:

```bash
bash "$BENCH_WORK/griffin_profile.sh" "$cores" \
    "$BENCH_WORK/warmup-sites.yaml" "$BENCH_WORK/warmup.chrom.sizes" \
    "$warmup_chromosomes" "$BENCH_WORK/warmup-output"
```

`warmup_chromosomes` lists the demo's chromosomes. The helper uses the same BAM, reference, profile width, length bounds, MAPQ and core allocation as production, with unit GC bias and the demo's own sites. No sample GC estimate is fitted for this warmup. The demo has a separate untimed Slurm step, and its outputs are removed before production measurements. There is no full-input Griffin warmup.

mosdepth stays in exact mode to match cfDNAlab's mate-overlap and deletion handling, although its help recommends fast mode for many uses. PanDepth per-base timing includes its native single-threaded gzip after worker calculations. samtools lengths is labeled "samtools stats (full QC report, insert-size section extracted)". Griffin's final plotting rule is omitted, while its GC-bias plots and compulsory FFT features remain timed.
