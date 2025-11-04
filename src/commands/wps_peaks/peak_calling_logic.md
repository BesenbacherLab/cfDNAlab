# Peak Calling Logic

This document explains how the `wps-peaks` command detects nucleosome peaks. It is written for future maintainers, reviewers, and anyone who needs to reason about the implementation without reverse engineering the original Snyder script. The goal is to keep the default behavior aligned with Snyder while making it straightforward to experiment with improved heuristics.

## Overview

We analyze windowed protection scores (WPS) computed at single base resolution. The pipeline converts raw per-base WPS into peak intervals by applying the following stages:

1. **Preprocessing**: optional smoothing, then subtracting a local (pre-smoothing) baseline.

2. **Positive run detection**: grouping consecutive positive residuals into candidate regions.

3. **Peak selection**: pruning candidate regions to a set of peak intervals with scores.

4. **Window aggregation**: collect peaks per processing tile and per requested genomic window.

All stages are designed to mirror Snyder et al. (2016) in default behavior while avoiding their hard-to-follow branching.

## Detailed Steps

### 1. Input Data and Definitions

We start with a dense in-memory array of raw WPS values for the tile segment under analysis. Each position is represented by its chromosome name plus a zero-based start offset and an exclusive end offset within that chromosome. The peak calling function also receives a mask that marks bases to ignore (blacklist or incomplete context) and optional genomic windows.

Terminology used throughout the code:

- `center_index`: zero-based index of the genomic position whose WPS we are processing.

- `raw_wps`: raw WPS value at `center_index`.

- `baseline_window`: the 1,000 bp neighborhood centered on `center_index`, covering 500 bp upstream and 500 bp downstream (matching Snyder's fixed window size).

- `residual_wps`: adjusted value after smoothing (if enabled) and baseline subtraction.

### 2. Optional Savitzky-Golay Smoothing

Snyder applies a Savitzky-Golay filter with a 21 bp window and a quadratic polynomial. They also mirror the signal at the edges so that every position can be smoothed without shrinking the output. We reproduce this behavior exactly:

1. Expand the signal by mirroring the necessary prefix and suffix (length equals half the SG window, meaning 10 bases on each side).

2. Convolve with the coefficient vector derived from the 21 bp, order 2 design.

3. Replace each center's numerator with the smoothed value while keeping the set of raw values used for the baseline window unchanged (the baseline still inspects the raw 1,000 bp neighborhood).

The smoothing flag is enabled by default and can be disabled (`--no-smoothing`) for experiments that require raw WPS.

### 3. Baseline Subtraction via Rolling Median (1,000 bp window)

To highlight local enrichments and suppress long-range trends, we subtract the median WPS from the surrounding 1,000 bp baseline window (500 bp on each side of the center). Snyder used the median instead of the mean because it resists outliers. We keep the same choice to stay compatible with Snyder's behavior. It works even when future versions introduce floating point corrections, such as GC scaling.

Implementation details:

1. Maintain a sliding multiset (e.g., two heaps or a balanced tree) over the 1,000 bp window.

2. For each `center_index`, compute the median of the raw WPS values in the window.

3. Store `residual_wps = (smoothed_or_raw_value - window_median)`.

The window slides one base at a time. When the mask marks any member of the window, we exclude that base from the multiset. If the window ends up with too few unmasked values (for example, fewer than 400 bases), we mark the residual as `NaN`. This approach avoids dilating blacklist regions beyond the mask we already computed during WPS collection, and it naturally handles mask entries that originate from tile edges rather than true blacklists.

### 4. Positive Run Detection

Snyder converts the residual stream into candidate regions by grouping consecutive positive values. They also allow short gaps (up to 5 bp) to bridge small dips inside a peak. We follow the same approach while keeping all work tile-local:

1. Iterate `center_index` in ascending order.

2. Whenever `residual_wps > 0`, append it to the current run.

3. If the residual is non-positive but the gap size is <= 5 bp, append zeros to preserve continuity.

4. If the gap exceeds the threshold or the masked region starts, close the run and evaluate it.

Each run tracks its genomic start and end, the sum of positive residuals, and the maximum residual encountered. Zero padding never increases the peak height but keeps the run connected.

### 5. Peak Selection Within Runs

Snyder enforces domain-specific constraints on run lengths (50-150 bp, or 50-150 bp segments carved out of longer runs) and keeps only the strongest sub-interval. We replicate this logic:

1. If the run length is < 50 bp or > 150 bp x 3, discard it.

2. If the run length is between 150 bp and 450 bp, split it into overlapping 150 bp windows sliding by one base, evaluate each, and keep sub-runs that satisfy the minimum length requirement. These sub-windows can overlap, but a later step picks at most one representative interval per parent run.

3. For each sub-run, compute its median residual. Only positions whose residual >= median are kept; the rest of the sub-run is ignored. If nothing survives the threshold, the sub-run is discarded and the parent run does not produce a peak at that offset.

4. Collapse adjacent retained positions into peak intervals.

5. For each interval compute:

   - `peak_height`: maximum residual inside the interval.

   - `peak_sum`: sum of residuals. This is the discrete area under the adjusted curve that Snyder maximizes, and it favors peaks that are both tall and wide without averaging away intensity.

   - Genomic `start` and `end` coordinates (inclusive start, exclusive end).

6. Keep the interval with the highest `peak_sum` as the representative peak for the parent run. This matches Snyder's final filter that emits the strongest sub-window and prevents duplicate peaks from a single run.

### 6. Outputs

For each chosen peak we output:

- Chromosome name.

- Start coordinate (0-based).

- End coordinate (exclusive).

- Peak height (float).

Additional statistics, such as peak width (computed as `end - start`) or peak area, are retained in memory so window-level summaries can be generated. Peaks outside the requested windows are ignored. Chromosomes without windows produce no output so windowed analyses are restricted to the provided intervals.

## Reproducibility and Flexibility

The implementation should default to Snyder-equivalent behavior, meaning the same inputs yield roughly the same peak intervals and heights (blacklisting is allowed to differ). At the same time, the code must expose extension points:

- Smoothing can be disabled or swapped for other filters.

- Baseline windows can be adjusted (e.g., switch to an 800 bp weighted rolling mean) or replaced with an alternative detrending function (e.g., subtraction of triangularly-weighted average).

- Gap thresholds, run length limits, and scoring functions can be adjusted.

To keep experimentation safe, we isolate each stage behind functions with clear arguments and ensure variable names describe domain concepts (`baseline_window`, `residual_wps`, `peak_runs`) instead of terse counters. These choices make the pipeline easy to read and modify without reintroducing Snyder's dense branching.

## Implementation Notes

- Operate on tile-sized buffers. Reuse the dilation already applied during WPS collection and enlarge it further (e.g., an extra 500 bp halo) whenever the baseline window would extend beyond the available WPS values, so both smoothing and baseline subtraction have the necessary neighborhood.

- Use descriptive structs and enums rather than tuples so intent is obvious at call sites.

- Keep helper functions small and pure where possible to simplify unit testing.

- Document any deviations from Snyder's defaults in-line and in this file.

- Guard experiments behind explicit configuration flags so reproducible defaults stay intact.

By following these guidelines we can stay close to the Snyder peak calls today, maintain them long term, and iterate on improvements without sacrificing clarity.
