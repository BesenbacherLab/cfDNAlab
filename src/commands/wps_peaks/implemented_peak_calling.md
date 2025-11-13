# Implemented WPS Peak Calling

This generated report captures the behavior that currently ships in `cfdna wps-peaks`. It is generated from code review rather than intent, so it reflects what the binary does today regardless of what `peak_calling_logic.md` specifies.

## Scope and entry points

- The command is compiled only when both `cmd_wps_peaks` and `cmd_wps` features are enabled. The CLI exposes it through `cfdna wps-peaks`.
- The configuration struct is `WPSPeaksConfig` (`src/commands/wps_peaks/config.rs`). It embeds `WPSSharedConfig`, so every filter that applies to `cfdna wps` (paired-end requirements, fragment length bounds, blacklist handling, scaling, tile sizing) is also in effect here.
- Peaks can be emitted for entire chromosomes (no windows) or confined to BED / fixed-size windows. Window modes drive the writer behavior (unique positions, indexed positions, stats).

## End-to-end data flow

1. **Resolve inputs**
   - Chromosomes are resolved once (`resolve_chromosomes_and_contigs`).
   - The output directory is created up front. Temporary tile files live under `tmp.<prefix>.*` and are deleted at the end of the run.
   - A blacklist map is mandatory whenever `--blacklist` is supplied. Each blacklist interval is dilated by `max_fragment_length + ceil(window_size / 2)` (code uses integer arithmetic). The halo prevents partially overlapping fragments from contaminating the baseline.
   - Optional scaling factors are loaded and later used to rescale per-position WPS after the core calculation (identical to `cfdna wps`).

2. **Tile construction**
  - `build_tiles` splits chromosomes into overlapping tiles. Each tile’s core is surrounded by a symmetric halo large enough to cover the WPS window. When `wps_for_tile` runs it extends both sides by `window_left + EXTRA_PEAK_HALO_BP` on the start and `window_right + EXTRA_PEAK_HALO_BP` on the end. The left extension also drives `initial_segment_marker`, ensuring the first peak in the tile is aware of the last masked region upstream.
   - For fixed-size windows, `build_tiles` also reports whether tile boundaries align with window boundaries. When aligned, downstream code can treat stats as pure bin aggregates without streaming per-peak data.

3. **Per-tile WPS computation**
  - `peaks_for_tile` calls `wps_for_tile` in positional mode. It ignores per-window aggregation at this stage and requests both the WPS vector and the associated mask.
  - `wps_for_tile` enforces fragment-level filters declared in `WPSSharedConfig`: fragment length bounds, MAPQ cutoff, inward orientation, required proper pair (depending on config), and blacklist checks.
  - The mask returned by `wps_for_tile` marks chromosome edges, dilated blacklist regions, and any bases whose baseline window would leave the fetch halo. `peaks_for_tile` pads/truncates the mask to match the WPS length before sending both downstream.

4. **Signal conditioning** (`PeakSignalProcessingOptions`)
   - **Smoothing:** Enabled unless `--no-smoothing` is passed. The Savitzky Golay filter uses a 21 bp, order-2 kernel (`SG_HALF_WINDOW = 10`). Masked positions split the signal into independent segments; each segment is mirrored on both sides before convolution.
   - **Normalization:** Always enabled in the CLI path. The smoothed numerator minus the sliding median of the *raw* WPS reference yields the residual. The median window size is `--normalize-bp` (default 1000) and the window slides by one base. If fewer than `--min-unmasked` bases remain inside the window (default 400), the residual is `NaN`.
   - **Initial segment markers:** Each tile inherits the coordinate of the last blacklist end before the dilated tile start (`last_mask_end_before`). The value is passed to peak calling as `initial_segment_marker` so that the first unmasked region in a tile starts a fresh segment ID even if the tile begins mid run.

5. **Peak detection** (`call_peaks.rs`)
   - Residuals (i.e. median-subtracted smoothed WPS) are scanned left to right. Masked bases and `NaN`s finalize the current run.
  - Short gaps (`MAX_GAP = 5 bp`) are bridged by inserting zeros, so minor dips do not break a run. Snyder’s script does the same by pushing `0` values while `ipos <= cend + 5` (see `snyder_code.md`).
   - Runs shorter than 50 bp or longer than 450 bp are discarded (`MIN_LENGTH = 50`, `MAX_LENGTH = 150`, `MAX_RUN_LENGTH = 3 * MAX_LENGTH`).
   - A run is filtered through its own median: only positions with residuals >= run median survive.
   - Surviving positions are collapsed into sub windows. For runs up to 150 bp the densest sub window (by summed residual) is emitted once. For longer runs (<= 450 bp) every contiguous sub window whose length sits within 50-150 bp and whose max residual exceeds `--min-peak-height` becomes a peak.
   - Each emitted `PeakCall` records: chromosome, inclusive start, exclusive end, `peak_position` (where the residual reaches its max), `height` (that max), and `segment_id`. The `segment_id` equals the barrier marker that was active when the run started and is incremented every time a mask or run reset occurs. Stats use it to avoid crossing masked gaps.
   - Peaks are clipped to the tile core before being persisted to ensure no peaks are included in more than one tile.

6. **Tile outputs and persistence**
  - Every tile writes its peaks to `tmp.<prefix>.tile.<chrom>.<idx>.peaks.tmp`. When stats are requested the tile also returns `WindowStatsContribution` structs in-memory (not on disk) so the writer can stitch aligned windows without re-reading peaks.
   - Non-windowed runs call `GlobalWriter` to concatenate every tile file into `<prefix>.wps.peaks.tsv.zst` with the header `chromosome start end peak_position height`.
   - When windows are requested (`--by-bed` or `--by-size`), `WindowOutputWriter` orchestrates the merge:
     - **Window sources**
       - `WindowSource::Bed`: the BED entries are read once. If `--per-window unique-positions` is selected, overlapping BED windows are merged up front, otherwise they are kept as-is.
       - `WindowSource::FixedSizeBuffered`: used when tile and window boundaries do not align. `FixedSizeWindows` streams windows per chromosome, maintains `next_start` and `next_idx`, and can spawn windows that begin before the current tile so long as they overlap it.
       - `WindowSource::FixedSizeAligned`: used when tiles align with bins. Because windows are deterministic, stats do not need to examine peaks directly; the writer only needs the per-window contributions that `peaks_for_tile` precomputed.
     - **Output modes**
       - `Unique`: per-window merge of peaks by genomic position. Peaks sharing the same base keep only the highest height. Outputs `chromosome start end height` (start=end-1).
       - `Indexed`: all peaks inside a window are emitted individually along with the zero-based window index. Peaks are not deduplicated.
       - `Stats`: per-window peak counts plus average and median inter-peak distances.
     - `WindowAccumulator` keeps state for all windows overlapping the running tile. It stores either a `BTreeMap` (unique) or a `Vec<PeakCall>` (indexed). For stats it stores counts, first/last peaks, first/last segments, and a histogram of observed distances.
     - Completed windows are flushed when their end coordinate drops behind the end of the current tile. A chromosome change flushes any remaining windows.

7. **Stats computation**
   - For each tile the code optionally pre-computes `WindowStatsContribution` objects. Each contribution covers one window and contains: peak count, first/last peak position, first/last segment ID, sum of distances, and a histogram keyed by distance (bp).
   - `compute_window_stats_contributions` enforces that distances are only measured within the same `segment_id`, so masked gaps or tile boundaries do not inflate inter-peak spacing.
   - When the writer processes tiles in stats mode: 
     - If contributions exist (window-aligned fixed-size runs), it merges them through `apply_stats_contribution`. This function also stitches the distance between the last peak of the previous tile and the first peak of the new tile **if** both belong to the same window **and** share the same `segment_id`.
     - If no contributions exist (BED windows or buffered fixed-size windows), peaks are streamed through `WindowAccumulator.push_peak` to build the histogram directly.
   - Average distances are `distance_sum / total_distances`. Medians are computed from the histogram via cumulative rank.
   - Windows with fewer than two peaks report `NaN` for both metrics.

8. **Temp files and cleanup**
   - After all tiles are processed, the temp directory is deleted. Failures only emit a warning; stale temp directories may persist if the process is interrupted.

## Configuration knobs and their real effects

| Option                        | Effect in code                                                                                                                                                                                             |
| ----------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--normalize-bp`              | Sliding-median window size inside `normalize_wps`. Values below 100 are rejected at CLI parsing.                                                                                                           |
| `--min-unmasked`              | Minimum usable bases for the median window. When unmet, the residual becomes `NaN` so peaks cannot start there. Lowering it lets sparse regions survive normalization but makes the median noisier.         |
| `--no-smoothing`              | Bypasses `smoothe_wps`. Residuals are then the raw WPS minus the median.                                                                                                                                   |
| `--min-peak-height`           | Final threshold applied to every peak. Snyder's `variCutoff` was implicitly 5; here it defaults to 3.0 (our placeholder) and is fully user-facing.                                                       |
| `--per-window`                | Required whenever `--by-bed` or `--by-size` is used. Drives the writer mode (unique, indexed, stats).                                                                                                      |
| `--by-size`                   | When tile and bin boundaries align, the code uses `FixedSizeAlignedWindows` and does not stream peaks to reconstruct stats. Misaligned runs fall back to buffered windows that accumulate peaks in-memory. |
| Blacklist options             | Dilated blacklist slices are passed both to `wps_for_tile` and to peak segmentation. They determine the `initial_segment_marker` and set mask bytes that terminate runs.                                   |
| Scaling (`--scaling-factors`) | Applied to the positional WPS before smoothing/normalization, same as in `cfdna wps`. Peak detection consequently operates on the scaled coverage.                                                         |

## Known limitations and TODO references

- `normalize_halo` adjustments were removed earlier; code comments still ask for better documentation (`wps_peaks.rs` around `peaks_for_tile`). The current implementation assumes that the halo provided by `wps_for_tile` plus `EXTRA_PEAK_HALO_BP` is sufficient. Edge masking near tile boundaries can still zero out potential peaks; this is noted in TODO comments.
- `WindowAccumulator::push_peak` relies on `scan_start` linear scans. With many overlapping windows this can become O(n^2). There is a TODO about making the function more readable and possibly more efficient.
- Stats currently ignore distances between peaks that fall in adjacent windows, even if the windows overlap. Each window’s histogram only sees peaks fully contained in that window, so the gap between the last peak of window A and the first peak of window B is never counted. `WindowAccumulator` has TODO comments about threading “first/last peaks” across windows to remedy this, but nothing is implemented yet.
- `FixedSizeWindows` uses in-memory `Vec`s per chromosome. If a chromosome has millions of small bins, the vector persists until the chromosome finishes. There is no streaming eviction beyond flushing once the chromosome changes.
- There is no guard against windows that are smaller than the minimum allowed run length (50 bp). Such windows will simply never accumulate peaks.
- BED windows supplied with `--per-window unique-positions` are merged prior to processing (see `into_flattened_reindexed`). That is expected because “unique” implies deduplicating overlaps, but it bears repeating: once merged, overlapping windows cannot be recovered in that mode.
- Global statistics. Maybe always add the global statistics?
- Allow extracting both peaks and statistics in the same run?
- COMPARE against Snyder! use their bins

## Depth generalization notes

- The pipeline expects coverage to be roughly normalized before thresholding. You can supply genomic scaling factors (`--scaling-factors`) which rescale raw coverage to an average near 1. This happens before smoothing and normalization, so the residuals become more depth-agnostic.
- After scaling, adjust `--min-peak-height` to match the dynamic range you see. The default 3.0 is a placeholder; the code contains a TODO to tune it empirically. In practice you should inspect peak distributions in representative samples (e.g., 1x cfDNA, 3x cfDNA, high-depth tissue) and set separate presets.
- Remaining depth effects stem from the rolling median window: low coverage produces noisier medians, and raising `--min-unmasked` beyond the available bases can blank entire regions. For ultra-low coverage consider reducing `--min-unmasked` so more centers survive (accepting that the median will be noisier because it is built from fewer bases).

## Differences relative to `peak_calling_logic.md`

| Topic               | Implemented behavior                                                                                                                                                                                                               | Notes                                                                                     |
| ------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------- |
| Extra halo          | The code adds a fixed 450 bp halo when fetching WPS and seeding `initial_segment_marker`. The logic doc mentions stretching halos when needed but does not specify this constant.                                                  | Document the constant and rationale (mirrors Snyder's 450 max run context).            |
| Stats stitching     | Distances across tiles are only counted when both peaks share a `segment_id`. The doc only states that masking should break runs but does not describe segment-aware stitching.                                                    | Doc should clarify that stats purposely ignore masked gaps even if windows overlap tiles. |
| Window merging      | Unique per-window outputs merge overlapping BED windows up front. Logic doc simply says “Peaks outside windows are ignored.”                                                                                                       | Document the pre-merge behavior so users know why duplicate positions disappear.          |
| Default threshold   | Implementation defaults to 3.0. Logic doc calls the Snyder cutoff “TODO: tune default empirically” but does not state the current value.                                                                                           | Keep doc synchronized.                                                                    |
| Experimental status | Config docstrings now mark `fragment-kmers`, `visualize-positions`, and `wps-peaks` as experimental features gated behind Cargo flags. `peak_calling_logic.md` does not mention that the subcommand can be disabled at build time. |

These differences do not block the current implementation but should be reconciled when the doc is next updated so that the design doc matches reality.
