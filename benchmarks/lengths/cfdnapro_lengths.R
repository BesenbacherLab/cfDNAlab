#!/usr/bin/env Rscript

# Read the BAM and write cfDNAPro's native fragment length histogram
# The runner times this complete process, including loading R packages
# Usage: Rscript --vanilla cfdnapro_lengths.R BAM OUTPUT MIN MAX MAPQ GENOME CURATE
# CURATE is TRUE for directional spans and FALSE for outer alignment spans

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 7L) {
    stop("Expected BAM OUTPUT MIN MAX MAPQ GENOME CURATE")
}
bam_path <- arguments[[1L]]
output_path <- arguments[[2L]]
minimum_length <- as.integer(arguments[[3L]])
maximum_length <- as.integer(arguments[[4L]])
minimum_mapq <- as.integer(arguments[[5L]])
genome_label <- arguments[[6L]]
if (!arguments[[7L]] %in% c("TRUE", "FALSE")) {
    stop("CURATE must be TRUE or FALSE")
}
curate_coordinates <- as.logical(arguments[[7L]])

suppressPackageStartupMessages(library(cfDNAPro))

# Match cfDNAlab's MAPQ and flag exclusions without requiring proper pairs
# Both mates must survive reading before cfDNAPro can form a fragment
# Keep read names, which its native outward-pair curation uses for ordering
# NM and MD tags are unnecessary for a length-only calculation
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

# Supplying fragment_obj is necessary because the BAM shortcut does not forward
# the MAPQ/flag options above. It would silently use readBam's defaults instead
counts <- read_bam_insert_metrics(
    fragment_obj = fragments,
    isize_min = minimum_length,
    isize_max = maximum_length
)
write.table(counts, output_path, sep = "\t", quote = FALSE, row.names = FALSE)
