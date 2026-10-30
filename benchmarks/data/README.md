# Benchmark inputs

The selected route is the published Snyder BH01 BAM, deterministic downsampling and merged-read reconstruction, then the lionheart hg38 remap with duplicates marked. The final BAM and biological site lists are external inputs. No checksums, observed depth, reconstruction counts or QC results are claimed here.

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

Keep the MarkDuplicates metrics beside these records. Measure proper-pair effects from the benchmark's native/proper-pair output totals after the actual filters. `flagstat` alone does not give that conditional fragment fraction. Do not reuse the original BAM's sampled non-proper percentage for the remap.

## Midpoint sites

Supply your midpoint site BED through `BENCH_INTERVALS`. The benchmark uses the supplied coordinates, groups and strands, retaining repeated and overlapping rows. It does not generate sites, select a subset or impose a site count or 2 kb window size. See the [midpoint input format](../midpoints/README.md#inputs-and-run) for the equal, even-width intervals required by the comparison. The runner saves the normalized rows and source checksum with the results.

## Small wiring input

```bash
samtools view -b -s 42.01 -o wiring.bam "$BENCH_BAM"
samtools index wiring.bam
```

This keeps approximately 1% by read-name hash, so both mates share a decision. Reduce the fraction if needed. Use `BENCH_PILOT=1`, one core and a small site set to check relocation, native CLI outputs and cfDNAlab's default QC plotting/fonts. Then pilot cfDNAPro on approximately 0.1x and increase input size to establish memory headroom. Use the chosen physical node for every tool. A wiring input cannot predict full-input runtime.
