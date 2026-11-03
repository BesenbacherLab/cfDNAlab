# Benchmark inputs

The BH01 preparation workflow starts with the published Snyder BAM, applies deterministic downsampling and merged-read reconstruction, then remaps to hg38 with the lionheart workflow and marks duplicates. Benchmark runners consume the final BAM and midpoint site BED as external inputs. Record file checksums, observed depth, reconstruction counts and QC results for each prepared input.

## BAM preparation

Download the original [BH01 BAM](https://kircherlab.bihealth.org/download/cfDNA/BH01.bam) and retain its checksum. The [Snyder processing description](https://shendurelab.github.io/cfDNA/#primary-data-processing-of-sequencing-data) uses GRCh37/hs37d5, rather than the hg38 reference used below.

```bash
"$BENCH_PREPARE_PYTHON" benchmarks/data/prepare_snyder.py \
    BH01.bam BH01.reconstructed.downsampled.bam \
    --keep-probability 0.15625 --seed 42 --threads 4 \
    2> preparation.log
```

The probability targets 15x from the published approximately 96x sequenced-base depth. Treat that as a target until measuring the remapped BAM. Sequenced-base depth and fragment coverage are different quantities. The script samples by read name, keeps existing paired records, and reconstructs only recognized `M_` merged reads. Unknown unpaired names fail. It preserves filtering flags and gives the forward reconstructed mate positive TLEN. These synthetic mates preserve the observed merged fragment boundaries, without recovering the original separate read sequences or qualities. The log records input, merged, retained, reconstructed and output counts.

Use the pinned [lionheart remap instructions](https://github.com/BesenbacherLab/lionheart/blob/085324a/remap/README.md) and [workflow](https://github.com/BesenbacherLab/lionheart/blob/085324a/remap/workflow.py): collate and convert to FASTQ, adapter/quality trimming, `bwa mem -p -Y -K 100000000`, coordinate sorting and Picard MarkDuplicates. Retain the actual invocation, workflow configuration, resolved program versions and duplicate metrics. Its installation example does not freeze every dependency. Index the final coordinate-sorted BAM with BAI.

Use the workflow's GRCh38 no-alt analysis reference with UCSC contig names. Record its exact download URL, release/basename and checksum, along with checksums of the derived uncompressed FASTA, `.fai` and `.2bit`. Derive the 2bit from the same FASTA. Matching contig names and lengths alone cannot verify matching reference bases. Remapping recomputes TLEN and proper-pair flags, so the final alignments, rather than the reconstruction convention, determine tool inclusion.

Save provenance once, outside timing:

```bash
# Set these paths to the actual final files
sha256sum "$BENCH_BAM" "$BENCH_INDEX" "$BENCH_REF_FASTA" \
    "$BENCH_FASTA_INDEX" "$BENCH_REF_2BIT" > inputs.sha256
md5sum "$BENCH_BAM" > final-bam.md5
samtools view -H "$BENCH_BAM" > final-bam.header.sam
samtools flagstat "$BENCH_BAM" > final-bam.flagstat.txt
samtools idxstats "$BENCH_BAM" > final-bam.idxstats.tsv
samtools --version > samtools-version.txt
```

Keep the MarkDuplicates metrics beside these records. Measure proper-pair effects from the benchmark's native/proper-pair output totals after the actual filters. `flagstat` alone does not give that conditional fragment fraction. Proper-pair rates from the original BH01 BAM do not describe the fragment population retained after remapping and benchmark filtering.

## Midpoint sites

Supply the midpoint site BED through `BENCH_INTERVALS`. The benchmark uses its coordinates, groups and strands, retaining repeated and overlapping rows. Window width comes from the supplied intervals, with no additional site selection. See the [midpoint input format](../midpoints/README.md#inputs-and-run) for the equal, even-width intervals required by the comparison. The runner saves the normalized rows and source checksum with the results.

Use Griffin's [public paper TFBS set](https://github.com/adoebley/Griffin_analyses/tree/75467e438a88029206d8d778ddce05eb0ab6efa8/sites/TFBS/10000_unfiltered_sites_CIS_BP_v2), pinned at `75467e438a88029206d8d778ddce05eb0ab6efa8`, as the production source. Record the supplied BED's derivation, window width and actual site/group counts. The runner does not choose or resize those production sites. Griffin's separate 1000-site release demo is used only for its short warmup.

## Check the setup with a small input

```bash
samtools view -b -s 42.01 -o check.bam "$BENCH_BAM"
samtools index check.bam
```

This keeps approximately 1% by read-name hash, so both mates share a decision. Reduce the fraction if needed. Use `BENCH_PILOT=1`, one core and a small site set to check relocation, native CLI outputs and cfDNAlab's default QC plotting/fonts. Then pilot cfDNAPro on approximately 0.1x and increase input size to establish memory headroom. Use the chosen physical node for every tool. A small input checks the setup but cannot predict full-input runtime.
