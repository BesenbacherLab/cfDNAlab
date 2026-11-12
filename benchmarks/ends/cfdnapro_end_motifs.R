#!/usr/bin/env Rscript

# Read the BAM and call cfDNAPro's native reference-derived left-end motifs
# The runner times the whole process, including loading R and the reference
# Usage: Rscript --vanilla cfdnapro_end_motifs.R BAM OUTPUT MIN MAX MAPQ GENOME K

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 7L) {
    stop("Expected BAM OUTPUT MIN MAX MAPQ GENOME K")
}
bam_path <- arguments[[1L]]
output_directory <- arguments[[2L]]
minimum_length <- as.integer(arguments[[3L]])
maximum_length <- as.integer(arguments[[4L]])
minimum_mapq <- as.integer(arguments[[5L]])
genome_label <- arguments[[6L]]
motif_length <- as.integer(arguments[[7L]])

suppressPackageStartupMessages(library(cfDNAPro))

# Match MAPQ and flag exclusions without requiring proper pairs or autosomes
# Native curation uses directional fragment spans. Both mates must survive
# Native alignment fields and NM/MD tag loading are left at their defaults
fragments <- readBam(
    bamfile = bam_path,
    genome_label = genome_label,
    chromosome_to_keep = FALSE,
    curate_start_and_end = TRUE,
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
fragment_lengths <- BiocGenerics::width(fragments)
fragments <- fragments[fragment_lengths >= minimum_length & fragment_lengths <= maximum_length]
if (length(fragments) == 0L) {
    stop("No fragments remain within the requested length bounds")
}

# Native type "s" extracts the leftmost k reference bases on the forward strand
# Do not double counts or substitute a both-end calculation
counts <- callMotif(
    frag_obj = fragments,
    motif_type = "s",
    motif_length = motif_length,
    integrate_mut = FALSE,
    downsample_ref = FALSE
)
if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
}
write.table(counts, file.path(output_directory, "cfdnapro-motifs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Retain the native reference manifest so package identity and scope are visible
# Sequence identity with the alignment reference still needs separate verification
reference <- S4Vectors::metadata(fragments)$cfDNAPro_reference
if (is.null(reference) || is.null(reference$seqinfo)) {
    stop("cfDNAPro did not retain its reference provenance")
}
reference_table <- reference$seqinfo
reference_table$label <- reference$label
reference_table$package <- reference$package
reference_table$package_version <- reference$package_version
write.table(reference_table, file.path(output_directory, "reference.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
