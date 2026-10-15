#!/usr/bin/env Rscript
# Plot selected loci from matched raw and outlier-weighted fcoverage runs, using base R
# Usage: Rscript scripts/plot_outliers.R raw.bedgraph.zst weighted.bedgraph.zst models.tsv audit.tsv loci.bed output_prefix

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 6L) {
    stop("Usage: plot_outliers.R raw.bedgraph[.zst] weighted.bedgraph[.zst] models.tsv audit.tsv loci.bed output_prefix", call. = FALSE)
}

# Read selected-locus outputs. Compressed bedGraphs require the zstd command-line executable
read_lines <- function(path) {
    if (!file.exists(path)) stop("Input does not exist: ", path)
    if (grepl("\\.zst$", path)) {
        executable <- Sys.which("zstd")
        if (!nzchar(executable)) stop("Reading .zst files requires zstd, or supply decompressed files")
        connection <- pipe(paste(shQuote(executable), "-q -d -c --", shQuote(path)), "r")
    } else {
        connection <- file(path, "r")
    }
    contents <- readLines(connection, warn = FALSE)
    status <- close(connection)
    if (is.numeric(status) && length(status) && status != 0) stop("Could not decompress ", path)
    contents[nzchar(trimws(contents)) & !grepl("^#", contents)]
}

read_table <- function(path, header) {
    contents <- read_lines(path)
    if (!length(contents)) return(data.frame())
    read.table(text = contents, header = header, sep = "\t", quote = "", comment.char = "",
               stringsAsFactors = FALSE, check.names = FALSE)
}

read_coverage <- function(path) {
    rows <- read_table(path, FALSE)
    if (!nrow(rows)) return(data.frame(chromosome = character(), start = numeric(), end = numeric(), value = numeric()))
    if (ncol(rows) != 4L) stop("Expected a four-column positional bedGraph: ", path)
    names(rows) <- c("chromosome", "start", "end", "value")
    if (any(!is.finite(rows$value)) || any(rows$value < 0) || any(rows$end <= rows$start)) {
        stop("Coverage must be finite and nonnegative, with nonempty intervals: ", path)
    }
    rows
}

raw <- read_coverage(arguments[[1]])
weighted <- read_coverage(arguments[[2]])
models <- read_table(arguments[[3]], TRUE)
audit <- read_table(arguments[[4]], TRUE)
loci <- read_table(arguments[[5]], FALSE)
if (ncol(loci) < 3L || !nrow(loci)) stop("loci.bed must contain chromosome, start, and end")
names(loci)[1:3] <- c("chromosome", "start", "end")
if (any(!is.finite(loci$start)) || any(!is.finite(loci$end)) || any(loci$start < 0) || any(loci$end <= loci$start)) {
    stop("Locus coordinates must be nonnegative, 0-based, and half-open")
}
required_models <- c("chromosome", "start", "end", "underlying_component_mean", "underlying_mean", "final_threshold")
required_audit <- c("chromosome", "start", "end", "action", "window_score")
if (!all(required_models %in% names(models))) stop("Model table lacks required columns")
if (!all(required_audit %in% names(audit))) stop("Audit table lacks required columns")
dir.create(dirname(arguments[[6]]), recursive = TRUE, showWarnings = FALSE)

subset_runs <- function(rows, locus) {
    selected <- rows[rows$chromosome == locus$chromosome & rows$end > locus$start & rows$start < locus$end, , drop = FALSE]
    selected <- selected[order(selected$start), , drop = FALSE]
    if (nrow(selected) > 1L && any(selected$start[-1L] < head(selected$end, -1L))) {
        stop("Coverage intervals overlap. Use fcoverage --per-window unique-positions")
    }
    selected
}

# Missing bedGraph runs represent zero coverage. The combined knots preserve exact run boundaries
coverage_at <- function(rows, positions) {
    values <- numeric(length(positions))
    if (!nrow(rows)) return(values)
    indices <- findInterval(positions, rows$start)
    selected <- which(indices > 0L)
    selected <- selected[positions[selected] < rows$end[indices[selected]]]
    values[selected] <- rows$value[indices[selected]]
    values
}

draw_locus <- function(locus) {
    raw_rows <- subset_runs(raw, locus)
    weighted_rows <- subset_runs(weighted, locus)
    knots <- sort(unique(c(locus$start, locus$end,
                           pmax(locus$start, raw_rows$start), pmin(locus$end, raw_rows$end),
                           pmax(locus$start, weighted_rows$start), pmin(locus$end, weighted_rows$end))))
    raw_values <- coverage_at(raw_rows, knots)
    weighted_values <- coverage_at(weighted_rows, knots)
    removed <- raw_values - weighted_values
    if (any(removed < -1e-5 * pmax(1, raw_values))) {
        warning("Weighted coverage exceeds raw coverage at ", locus$chromosome, ":", locus$start,
                ". Check that filters and other corrections match")
    }
    widths <- diff(knots)
    raw_mass <- sum(head(raw_values, -1L) * widths)
    removed_mass <- sum(head(removed, -1L) * widths)
    removed_percentage <- if (raw_mass > 0) 100 * removed_mass / raw_mass else NA_real_
    local_models <- subset_runs(models, locus)
    local_calls <- audit[audit$chromosome == locus$chromosome & audit$end > locus$start & audit$start < locus$end, , drop = FALSE]
    horizontal <- (knots - locus$start) / 1000
    x_limits <- c(0, (locus$end - locus$start) / 1000)
    layout(matrix(1:4, ncol = 1), heights = c(1, 1, 0.8, 0.5))
    par(mar = c(3.4, 4.5, 2.4, 1))

    for (logarithmic in c(FALSE, TRUE)) {
        transform <- if (logarithmic) log1p else identity
        model_values <- c(local_models$final_threshold, local_models$underlying_component_mean, local_models$underlying_mean)
        limits <- range(transform(c(0, raw_values, weighted_values, model_values)), finite = TRUE)
        if (diff(limits) == 0) limits[[2]] <- limits[[1]] + 1
        plot(horizontal, transform(raw_values), type = "s", col = "#777777", xlim = x_limits,
             ylim = limits, xlab = "Position from locus start (kb)",
             ylab = if (logarithmic) "log(1 + coverage)" else "Fragment coverage",
             main = if (logarithmic) "Coverage and fitted baseline" else
                 sprintf("%s:%s-%s | %.3f%% of raw coverage mass removed", locus$chromosome, locus$start, locus$end, removed_percentage))
        lines(horizontal, transform(weighted_values), type = "s", col = "#2166ac", lwd = 1.2)
        if (nrow(local_models)) {
            starts <- (pmax(locus$start, local_models$start) - locus$start) / 1000
            ends <- (pmin(locus$end, local_models$end) - locus$start) / 1000
            segments(starts, transform(local_models$underlying_component_mean), ends,
                     transform(local_models$underlying_component_mean), col = "#1b9e77", lwd = 1.4)
            segments(starts, transform(local_models$underlying_mean), ends,
                     transform(local_models$underlying_mean), col = "#7570b3", lty = 3)
            segments(starts, transform(local_models$final_threshold), ends,
                     transform(local_models$final_threshold), col = "#d95f02", lty = 2)
        }
        if (logarithmic) legend("topright", c("Raw", "Weighted", "Active mean", "Unconditional mean", "T2"),
                                col = c("#777777", "#2166ac", "#1b9e77", "#7570b3", "#d95f02"),
                                lty = c(1, 1, 1, 3, 2), bty = "n", cex = 0.75, ncol = 3)
    }
    plot(horizontal, removed, type = "s", xlim = x_limits, col = "#b2182b",
         xlab = "Position from locus start (kb)", ylab = "Coverage removed", main = "Complete-fragment correction, including outside exact calls")
    plot(NA, xlim = x_limits, ylim = c(0.5, 2.5), yaxt = "n", xlab = "Position from locus start (kb)",
         ylab = "", main = "Exact candidates and final merged-region decisions")
    axis(2, at = c(1, 2), labels = c("Preserved", "Actionable"), las = 1, cex.axis = 0.8)
    if (nrow(local_calls)) {
        levels <- ifelse(local_calls$action == "preserve", 1, 2)
        segments((pmax(locus$start, local_calls$start) - locus$start) / 1000, levels,
                 (pmin(locus$end, local_calls$end) - locus$start) / 1000, levels,
                 col = ifelse(levels == 1, "#e6ab02", "#b2182b"), lwd = 5)
    }
}

for (locus_index in seq_len(nrow(loci))) {
    locus <- loci[locus_index, , drop = FALSE]
    suffix <- paste(gsub("[^A-Za-z0-9_.-]", "_", locus$chromosome), locus$start, locus$end, sep = "_")
    prefix <- paste(arguments[[6]], suffix, sep = ".")
    png(paste0(prefix, ".png"), width = 1800, height = 1600, res = 150)
    draw_locus(locus)
    dev.off()
    pdf(paste0(prefix, ".pdf"), width = 12, height = 10.7)
    draw_locus(locus)
    dev.off()
    message("Saved ", prefix, ".{png,pdf}")
}
