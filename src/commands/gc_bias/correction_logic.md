# GC Correction Logic

This document captures the working plan for our GC-bias correction pass. It mirrors the GCparagon strategy where it helps, and simplifies where our pipeline can process whole genomes quickly. Every fragment-length/GC bin ends up with a weight centered around 1.0 so downstream coverage tools can treat values above 1 as down-weighted and below 1 as up-weighted.

## 1. Build the interval grid
- Default to 1 Mb windows so examples stay manageable. Users with more compute can switch to genome-wide or custom window sizes.
- Always apply the blacklist up front. A window with all bases masked simply contributes zero fragments and drops out naturally.
- Precompute and store the reference GC matrix (`reference_gc`) for every allowed window so we never resimulate reference fragments at runtime.

## 2. Count observed fragments
- For each window, iterate all fragments once (`bam_to_frag` can do this in seconds) and fill a 2D histogram indexed by fragment length and absolute GC count.
- Keep per-window metadata: total fragments, mean cell value, any QC flags (e.g., high-N content, too few fragments).

## 3. Derive per-window correction matrices
- Fetch the matching reference matrix and compute `window_weight = observed / reference` element-wise.
- Use the window’s mean observed count as a weight when combining windows later; this lets high-coverage regions drive the final correction without ever resorting to unbounded ratios.
- If a window fails QC (missing reference bases, too many Ns, no fragments), drop it and optionally log it for diagnostics.
- Keep a per-window mask of fragment-length/GC bins with at least `min_frag_occurs` observations on both the observed and reference sides. Bins that miss this threshold stay at the neutral weight (set them back to 1.0 after computing the window-level mean-scaled correction, then apply the mean coverage weight) so sparse regions cannot explode.
- Discard windows that fall below a minimum usable-fragment threshold (e.g., <5k fragments after masking). This protects against windows where blacklisting removed nearly everything.
- Also drop windows whose mask covers more than a set percentage (for example 80%) of bins; if almost every bin reverted to 1.0 there is no meaningful signal to contribute.

## 4. Aggregate windows
- Multiply each per-window correction matrix by its mean observed count (or store that scalar as the weight). This ensures the final average reflects fragment availability—high-coverage windows steer the combined matrix, while sparse windows cannot inject outlier ratios.
- Sum the weighted matrices and divide by the sum of weights. With full-genome processing this is equivalent to a global computation, but the per-window path lets us drop or inspect windows independently.
- Divide the combined matrix by its own mean so weights are roughly 1-centered before post-processing.

## 5. Post-processing (borrowed from GCparagon defaults)
- **Outlier trimming:** run a configurable pass (e.g., IQR-based) that clips extreme cells to within `(10 - stringency)` standard deviations of the median. This prevents bins with very low coverage from exploding.
- **Smoothing:** apply a Gaussian blur across fragment length and GC axes. Intensity defaults to 5 (same as GCparagon’s preset) but remains user-configurable. Cells that were exactly 1.0 remain untouched so unknown bins stay neutral.
- After each transformation, re-center by dividing through the matrix mean. This keeps the global correction neutral even if smoothing or clipping shifts the average.

## 6. Output and tagging
- Persist both the raw averaged matrix and the final post-processed matrix; the raw version helps debugging.
- When tagging BAMs, look up the fragment-length/GC bin. Values outside the matrix bounds fall back to 1.0.
- Any optional genomic smoothing on downstream counts should reuse the same window grid so users can toggle it without recomputing GC weights.

This plan keeps the code straightforward (no interval selection heuristics) while preserving the proven bits from GCparagon: per-window normalization, weighted averaging, outlier limitation, smoothing, and consistent mean-centering.
