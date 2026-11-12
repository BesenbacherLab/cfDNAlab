# Benchmark experiment commands

Each measured repetition runs the native commands below on the same staged BAM. Commands are shown once, with variables for the allocated cores, paths and experimental settings. Executable aliases resolve to the staged paths used by the runners. Shell argument arrays are written as ordinary command lines, and timed shell pipelines are expanded into their constituent commands.

The configurations describe the current runners. For an existing result, `run.tsv` identifies the benchmark revision and `timings.tsv` records the tool version, filters and core allocation. Installation and cfDNAlab build records identify the installed tools independently of the benchmark revision.

## Shared settings and timing

| Setting | Coverage and lengths | End motifs | Midpoints |
| --- | --- | --- | --- |
| Allocated cores | 1, 4, 12 or 32, in separate jobs | Same | Same |
| Minimum MAPQ | 30 | 30 | 30 |
| Included fragment lengths | 30-1000 bp inclusive | Requested 30-1000 bp inclusive | 100-200 bp inclusive for profiling |
| Chromosome selection | All BAM contigs, except named autosome configurations | All BAM contigs, except named autosome configuration | Chromosomes in the supplied site set |
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
coverage_output="$BENCH_WORK/coverage"
regions="$BENCH_WORK/regions.bed"
reference_2bit="$BENCH_WORK/reference.2bit"

cfdna="$BENCH_WORK/bin/cfdna"
mosdepth="$BENCH_WORK/bin/mosdepth"
pandepth="$BENCH_WORK/bin/pandepth"
samtools="$BENCH_WORK/samtools-env/bin/samtools"
finale_python="$BENCH_WORK/finale-env/bin/python"
finale="$BENCH_WORK/finale-env/bin/finaletoolkit"
picard="$BENCH_WORK/picard-env/bin/picard"
cfdnapro_image="$BENCH_WORK/cfdnapro.sif"
```

`BENCH_WORK` is the job's temporary directory under the compute node's `$TMPDIR`. The runners copy the BAM, index, selected tools and required references there before timing. `output` and `coverage_output` are recreated between repetitions. Container paths under `/work` refer to the same directory through an Apptainer bind mount.

Every available native thread/core/worker option receives `cores` directly. Commands without a native thread option still run within the allocated CPU set. This also applies to samtools inside the Picard and FrEIA pipelines.

The native options have different meanings. For example, the [mosdepth CLI](https://github.com/brentp/mosdepth/blob/v0.3.14/README.md#usage) defines `--threads` as BAM decompression threads, while coverage calculation runs separately. A 4-core allocation uses `--threads 4`, and all of the tool's threads share those allocated CPUs. `timings.tsv` records the core allocation and the native thread/worker setting separately.

The installer requests mosdepth 0.3.14, samtools 1.25, PanDepth 2.26, FinaleToolkit 1.1.0, cfDNAPro 1.7.4, Picard 3.5.0 and Griffin v0.2.0. FrEIA is pinned to `2c32fb9d9fab774f4ae57931df7fd7d75654d343`. cfDNAlab is built from the checkout recorded in `cfdnalab-build.tsv`.

The exclusion mask `3852`, or `0xF0C`, combines unmapped (`4`), mate unmapped (`8`), secondary (`256`), QC-fail (`512`), duplicate (`1024`) and supplementary (`2048`) flags. It does not require the proper-pair flag. FrEIA's native preprocessing uses its own exclusions, shown explicitly below.

GNU time surrounds each complete native command or shell pipeline. `wall_seconds` includes startup, BAM reading, native calculations, output writing/compression and enabled default plots. The full R process and Apptainer startup are included for cfDNAPro. Pipeline shells use `set -euo pipefail`.

After native success, a separate timer measures `sync -f` on the output filesystem. `process_plus_flush_seconds` adds that duration to `wall_seconds`. Staging, Slurm startup, CPU-placement checks, saved-profile copying, cleanup and concordance are outside the native command timer. Per-command timeouts are disabled when `BENCH_RUN_TIMEOUT_SECONDS=0`, the default. Slurm's allocation limit applies to the whole job.

## Coverage

The [coverage runner](coverage/benchmark_coverage.sh) contains three comparisons: positional read depth, positional fragment coverage, and regional summaries. Regional read depth and regional fragment counts are compared separately.

### Read depth at individual bases

**cfDNAlab, method `cfdnalab`:**

```bash
"$cfdna" fcoverage \
    --bam "$bam" --output-dir "$coverage_output" \
    --min-mapq "$min_mapq" \
    --min-fragment-length "$min_length" --max-fragment-length "$max_length" \
    --n-threads "$cores" --chromosomes all --ignore-gap
```

**mosdepth, method `mosdepth`:**

```bash
"$mosdepth" --threads "$cores" \
    --include-flag 1 --flag 3852 --mapq "$min_mapq" \
    --min-frag-len "$min_length" --max-frag-len "$max_length" \
    "$coverage_output/coverage" "$bam"
```

**samtools, method `samtools-depth`:**

```bash
"$samtools" depth -@ "$cores" \
    --input-fmt-option "filter=(tlen >= $min_length && tlen <= $max_length) || (tlen <= -$min_length && tlen >= -$max_length)" \
    -aa -s -q 0 -Q "$min_mapq" -G 3852 --require-flags 1 \
    -o "$coverage_output/samtools.depth.tsv" "$bam"
```

Here `-q 0` is the base-quality threshold and `-Q` is MAPQ. `-s` removes mate overlap, and `-aa` includes zero-depth positions. The input expression applies inclusive absolute TLEN bounds to both signs during BAM reading.

**PanDepth, method `pandepth`:**

```bash
"$pandepth" -i "$bam" -o "$coverage_output/pandepth" \
    -q "$min_mapq" -x 3852 -t "$cores" -a
```

PanDepth counts overlapping mates twice and has no fragment length filter. cfDNAlab merges mate overlap. mosdepth's native overlap correction depends on proper-pair status. The [concordance specification](benchmarks_concordance_spec.md#native-differences) describes the remaining native pairing and CIGAR differences.

### Fragment coverage at individual bases

**cfDNAlab, method `cfdnalab-fragment`:**

```bash
"$cfdna" fcoverage \
    --bam "$bam" --output-dir "$coverage_output" \
    --min-mapq "$min_mapq" \
    --min-fragment-length "$min_length" --max-fragment-length "$max_length" \
    --n-threads "$cores" --chromosomes all
```

**mosdepth, method `mosdepth-fragment`:**

```bash
"$mosdepth" --threads "$cores" \
    --include-flag 1 --flag 3852 --mapq "$min_mapq" \
    --min-frag-len "$min_length" --max-frag-len "$max_length" \
    --fragment-mode "$coverage_output/coverage" "$bam"
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

cfDNAlab and mosdepth generate the same window coordinates through their fixed-size options. samtools and FinaleToolkit read the shared BED. PanDepth 2.26 interprets its BED input as 1-based inclusive, so its scratch-only BED3 is prepared outside timing:

```bash
awk 'BEGIN {OFS="\t"} {print $1, $2+1, $3}' \
    "$regions" > "$BENCH_WORK/pandepth-regions.bed"
```

### Regional read depth

These four methods measure summed read depth, or covered-base totals, in each 500 kb window. No fragment length normalization is applied.

**cfDNAlab, method `cfdnalab-regions`:**

```bash
"$cfdna" fcoverage \
    --bam "$bam" --output-dir "$coverage_output" \
    --min-mapq "$min_mapq" \
    --min-fragment-length "$min_length" --max-fragment-length "$max_length" \
    --n-threads "$cores" --chromosomes all --ignore-gap \
    --by-size 500000 --per-window total
```

**mosdepth, method `mosdepth-regions`:**

```bash
"$mosdepth" --threads "$cores" \
    --include-flag 1 --flag 3852 --mapq "$min_mapq" \
    --min-frag-len "$min_length" --max-frag-len "$max_length" \
    --by 500000 --no-per-base "$coverage_output/coverage" "$bam"
```

**samtools, method `samtools-regions`:**

```bash
"$samtools" bedcov \
    --input-fmt-option "nthreads=$cores" \
    --input-fmt-option "filter=flag.paired && ((tlen >= $min_length && tlen <= $max_length) || (tlen <= -$min_length && tlen >= -$max_length))" \
    -j -Q "$min_mapq" -G 3852 "$regions" "$bam" \
    > "$coverage_output/samtools.bedcov.tsv"
```

**PanDepth, method `pandepth-regions`:**

```bash
"$pandepth" -i "$bam" -o "$coverage_output/pandepth" \
    -q "$min_mapq" -x 3852 -t "$cores" \
    -b "$BENCH_WORK/pandepth-regions.bed"
```

These commands write regional output directly. mosdepth's native regional means are multiplied by the actual window widths during concordance, outside timing. Their default two-decimal precision can introduce rounding differences in the resulting totals. PanDepth supplies `TotalDepth`, and the comparator converts its coordinates back to BED coordinates. samtools `-j` excludes D/N positions, but `bedcov` does not remove mate overlap. PanDepth also retains overlapping mates and has no length filter. Regional summaries do not replace the separate base-resolution concordance calculation.

### Regional fragment counts

**cfDNAlab, method `cfdnalab-region-counts`:**

```bash
"$cfdna" fcoverage \
    --bam "$bam" --output-dir "$coverage_output" \
    --chromosomes all --by-size 500000 --per-window total --normalize-by-length \
    --min-fragment-length "$min_length" --max-fragment-length "$max_length" \
    --min-mapq "$min_mapq" --n-threads "$cores"
```

**FinaleToolkit, method `finaletoolkit-coverage`:**

```bash
"$finale_python" "$finale" coverage "$bam" "$regions" \
    --output "$coverage_output/finaletoolkit.coverage.bed.gz" \
    --intersect-policy any --scale-factor 1 \
    --min-length "$min_length" --max-length "$max_length" \
    --min-mapq "$min_mapq" --threads "$cores"
```

cfDNAlab gives each fragment total weight 1 before output rounding, split across windows when it crosses a boundary. FinaleToolkit counts a fragment once in every overlapping window and requires proper read1/TLEN selection. Its indexed read1 query can miss fragments whose span alone overlaps a boundary. These counts are compared separately from regional read-depth totals. The end-motif length-filter bug described below does not apply to this coverage call.

### Additional coverage configurations

| Method | Exact change to the corresponding call above |
| --- | --- |
| `cfdnalab-autosomes` | Read-depth call, replacing `--chromosomes all` with `--chromosomes "$BENCH_AUTOSOMES"` |
| `cfdnalab-fragment-proper` | Fragment-coverage call plus `--require-proper-pair` |
| `cfdnalab-region-counts-proper` | Normalized regional-count call plus `--require-proper-pair` |

`BENCH_AUTOSOMES` is the comma-separated set of BAM contigs named `chr1` through `chr22`, or `1` through `22`.

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
"$BENCH_APPTAINER" exec --cleanenv --containall \
    --bind "$BENCH_WORK:/work" --pwd /work \
    "$cfdnapro_image" Rscript --vanilla /work/cfdnapro_lengths.R \
    /work/sample.bam /work/output/lengths.tsv \
    "$min_length" "$max_length" "$min_mapq" "$BENCH_GENOME" TRUE
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

Method `picard`. Both programs below share one timer and the same allocated CPU set:

```bash
"$samtools" view -@ "$cores" -u -q "$min_mapq" -f 1 -F 3852 "$bam" |
    "$picard" -Xms512m -XX:ActiveProcessorCount="$cores" \
        -Djava.io.tmpdir="$BENCH_WORK/tmp" CollectInsertSizeMetrics \
        INPUT=/dev/stdin OUTPUT="$output/metrics.txt" \
        HISTOGRAM_WIDTH="$max_length" MINIMUM_PCT=0 INCLUDE_DUPLICATES=false
```

The input filter supplies MAPQ selection absent from this Picard command. `-Xms512m` preserves the launcher's initial heap while preventing its fixed maximum-heap fallback. Java chooses the maximum heap. No PDF is requested. The comparator later selects the all-reads FR histogram and requested length range.

### Additional length configurations

| Method | Exact change |
| --- | --- |
| `cfdnalab-autosomes` | cfDNAlab call with `--chromosomes "$BENCH_AUTOSOMES"` instead of `all` |
| `cfdnalab-proper` | cfDNAlab all-contig call plus `--require-proper-pair` |
| `cfdnapro-outer-span` | cfDNAPro container call with its final argument changed from `TRUE` to `FALSE` |

cfDNAlab uses directional paired spans. Curated cfDNAPro also uses directional spans, while its uncurated diagnostic uses outer alignment spans. FinaleToolkit, samtools and Picard use native TLEN-based selection and can retain different fragment populations.

## End motifs

The [end-motif runner](ends/benchmark_end_motifs.sh) uses `k=$BENCH_K`. The four-mer comparison selects `cfdnalab-reference`, `cfdnalab-reference-proper`, `finaletoolkit` and `cfdnapro`. The three-mer comparison selects `cfdnalab-read-aligned` and `freia`. Every selected tool can be allocated each requested core setting, including tools without a native thread option.

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

FinaleToolkit 1.1.0 ignores the supplied motif length bounds. Timing rows record `length_filter=ignored-by-release` and unknown effective bounds. It also has a native last-window omission when a contig length is exactly divisible by 1,000,000, which the runner rejects before timing. This call uses reference-derived motifs and proper read1/TLEN selection.

### cfDNAPro

Method `cfdnapro`:

```bash
"$BENCH_APPTAINER" exec --cleanenv --containall \
    --bind "$BENCH_WORK:/work" --pwd /work \
    "$cfdnapro_image" Rscript --vanilla /work/cfdnapro_end_motifs.R \
    /work/sample.bam /work/output \
    "$min_length" "$max_length" "$min_mapq" "$BENCH_GENOME" "$k"
```

The native R calculation below selects reference-derived left ends. Reference loading, BAM reading, curation, motif calling and output/provenance writing are timed. The package's bundled reference is selected by `BENCH_GENOME`, normally `hg38`. Sequence identity with the supplied alignment reference remains to be established before interpreting differences solely as tool differences.

### FrEIA

Method `freia` requires `k=3`. Its extraction and aggregation use the pinned environment's own Python and samtools:

```bash
freia_environment="$BENCH_WORK/freia"
freia_scripts="$freia_environment/share/freia"
```

The samples table is prepared outside the timer:

```bash
printf 'sample_name group\nsample benchmark\n' > "$output/samples.tsv"
```

All four native calls below share one timer:

```bash
"$freia_environment/bin/samtools" view -@ "$cores" -q "$min_mapq" \
    -F 4 -F 256 -F 1024 -F 2048 -b \
    -o "$output/filtered.bam" "$bam"
"$freia_environment/bin/samtools" index -@ "$cores" "$output/filtered.bam"
"$freia_environment/bin/python" "$freia_scripts/1_extract_fragment_ends.py" \
    -b "$output/filtered.bam" -o "$output/sample.pq" \
    -c "$BENCH_CONTIGS" -t "$cores" -n 3
"$freia_environment/bin/python" "$freia_scripts/2_data_transformation.py" \
    -i "$output/sample.pq" -o "$output/" -st "$output/samples.tsv" \
    -fsmin "$min_length" -fsmax "$length_end" \
    -subs 1 --bootstrap_sample 0 --subchroff
```

`BENCH_CONTIGS` lists all BAM contigs. The maximum is supplied as `max_length + 1` because FrEIA's upper bound is exclusive. The native transformation retains coordinate deduplication, rejects fragments containing N in either extracted end, and calculates mononucleotide as well as trinucleotide abundances. There is no subsampling or bootstrap resampling in this call. The extractor itself forces one worker when its filtered BAM exceeds 8000 MiB, and the timing row records that worker count.

The broken native `--both` aggregate is not requested. The intact left-only table is retained as `left-motifs.pq`. FrEIA and cfDNAPro left-end frequencies are compared with cfDNAlab/FinaleToolkit both-end frequencies using the documented side-mirroring approximation. Counts are not doubled.

### Additional motif configurations

| Method | Exact change |
| --- | --- |
| `cfdnalab-read` | Aligned-read call with `--clip-strategy skip` |
| `cfdnalab-read-include` | Aligned-read call with `--clip-strategy include-at-aligned-boundary` |
| `cfdnalab-reference-autosomes` | Reference call with `--chromosomes "$BENCH_AUTOSOMES"` instead of `all` |
| `cfdnalab-reference-proper-wide` | Reference call with `--require-proper-pair`, minimum length 10 and maximum length 50000 |

These configurations are only run when explicitly selected in `BENCH_METHODS`.

## cfDNAPro BAM reading and calculations

The complete R helpers are [cfdnapro_lengths.R](lengths/cfdnapro_lengths.R) and [cfdnapro_end_motifs.R](ends/cfdnapro_end_motifs.R). Their computational calls are expanded here so the container invocation does not hide the filtering or calculation. Argument parsing and error checks remain in the linked helpers.

Both helpers read the BAM explicitly with the following options. `bam_path`, `genome_label`, `minimum_mapq`, `minimum_length` and `maximum_length` are their positional command-line arguments. `curate_coordinates` is `TRUE` for normal lengths and motifs, and `FALSE` only for `cfdnapro-outer-span`.

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
    galp_tag = character(),
    outdir = NA
)
```

The **length** calculation is:

```r
counts <- read_bam_insert_metrics(
    fragment_obj = fragments,
    isize_min = minimum_length,
    isize_max = maximum_length
)
write.table(counts, output_path, sep = "\t", quote = FALSE, row.names = FALSE)
```

The **motif** calculation applies inclusive length bounds to the curated fragments, then calls the native left-end function. `motif_length` is the supplied `k` argument:

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

The motif helper also writes `reference.tsv` from the native fragment object's reference metadata. Both reads and writes occur inside the R process timer.

## Midpoint profiles

The [midpoint runner](midpoints/benchmark_midpoints.sh) uses the supplied grouped windows, normally 2001 bp wide, and profiling length bounds 100-200 bp. `BENCH_CHROMOSOMES` is derived from the chromosomes in those windows. The runner preserves the sites, including repeated and overlapping rows, and creates the native Griffin anchor tables outside timing. Each workload runs in its own job.

### cfDNAlab without GC correction

Method `cfdnalab`:

```bash
"$cfdna" midpoints --bam "$bam" --output-dir "$output" \
    --intervals "$BENCH_WORK/sites.bed" --chromosomes "$BENCH_CHROMOSOMES" \
    --length-bins "$min_length" "$length_end" \
    --min-mapq "$min_mapq" --smoothing none --bin-size 1 --n-threads "$cores"
```

### Griffin profiling calls

Method `griffin` times both native profiling scripts in [griffin_profile.sh](midpoints/griffin_profile.sh). The following variables express that helper's production arguments:

```bash
griffin_python="$BENCH_WORK/griffin-env/bin/python"
griffin_scripts="$BENCH_WORK/griffin-env/share/griffin"
sites_yaml="$BENCH_WORK/sites.yaml"
chromosome_sizes="$BENCH_WORK/chrom.sizes"
reference_fasta="$BENCH_WORK/reference.fa"
half_width=$BENCH_HALF_WIDTH
window_end=$((BENCH_WINDOW_WIDTH - half_width))
gc_bias="$BENCH_WORK/unit_gc_bias.tsv"
IFS=, read -r -a chromosomes <<< "$BENCH_CHROMOSOMES"
```

For a 2001 bp window, every relative window argument below is `-1000 1001`, with an exclusive end. The arguments shared by the two scripts are:

```bash
griffin_common=(
    --sample_name sample --tmp_dir "$output/tmp"
    --sites_yaml "$sites_yaml" --griffin_scripts_dir "$griffin_scripts"
    --chrom_sizes_path "$chromosome_sizes" --chroms "${chromosomes[@]}"
    --chrom_column Chrom --position_column position --strand_column Strand
    --norm_window "-$half_width" "$window_end"
    --mappability_correction False --mappability_bw none --CPU "$cores"
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
    --save_window "-$half_width" "$window_end" \
    --center_window "-$half_width" "$window_end" \
    --fft_window "-$half_width" "$window_end" --fft_index 1 \
    --step 1 --exclude_paths none --exclude_zero_mappability False \
    --exclude_outliers False --CNA_normalization False --individual False --smoothing False
```

For the uncorrected experiment, setup supplies unit GC bias for all requested lengths and GC counts. Griffin still performs its compulsory reference/GC work and writes uncorrected and weighted bigWigs inside timing. The retained profile is `sample.uncorrected.coverage.tsv`. Its native group normalization gives mean-one profiles. Concordance normalizes cfDNAlab's group profile over the complete window outside timing, without shifting, smoothing or rebinning.

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
    --min-mapq "$min_mapq" --smoothing none --bin-size 1 --n-threads "$cores" \
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

The final helper invocation expands to the two Griffin profiling calls above, with the newly fitted `gc_bias`. Estimation uses Griffin's native 15-500 bp range. Profiling uses the supplied 100-200 bp range. MAPQ is 30 for both. The estimator's region BED and reference GC frequencies come from Griffin's bundled inputs, independently of cfDNAlab's reference package. `BENCH_GC_REFERENCE_NAME` is the estimation BED's basename without its final extension. The retained profile is `sample.GC_corrected.coverage.tsv`.

### Griffin demo warmup outside timing

Before the first measured Griffin repetition, both corrected and uncorrected jobs run the release's 1000-site CTCF demo through the same profiling helper:

```bash
bash "$BENCH_WORK/griffin_profile.sh" "$cores" \
    "$BENCH_WORK/warmup-sites.yaml" "$BENCH_WORK/warmup.chrom.sizes" \
    "$warmup_chromosomes" "$BENCH_WORK/warmup-output"
```

`warmup_chromosomes` lists the demo's chromosomes. The helper uses the same BAM, reference, profile width, length bounds, MAPQ and core allocation as production, with unit GC bias and the demo's own sites. No sample GC estimate is fitted for this warmup. The demo has a separate untimed Slurm step, and its outputs are removed before production measurements. There is no full-input Griffin warmup.
