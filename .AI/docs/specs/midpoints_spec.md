# midpoints Spec

`cfdna midpoints` builds grouped midpoint profiles around fixed-width genomic sites. The public output is dense, while tile internals are sparse.

## Input Contract

- `--intervals` is required and must be BED-like with columns `chromosome`, `start`, `end`, and `group_name`.
- The public CLI expects sorted intervals. The current loader sorts selected windows internally, but group index assignment still follows first observed group name in the input file.
- All selected intervals must have the same length.
- Sites with the same group name collapse into one group profile.
- Group indices are assigned by first observed group name during grouped BED loading.
- Optional interval strand tokens are `+`, `-`, or `.`. Files with six or more columns read strand only from column 6. Five-column files may read strand from column 5. A strand-like column 5 with a non-strand column 6 is rejected as ambiguous.
- Missing strand information is treated as unstranded.
- Chromosome filtering may remove rows, but the command must fail if no selected grouped windows remain.
- With blacklists, intervals whose output span plus `ceil(max_fragment_length / 2) + smoothing_flank` overlaps a blacklisted region are removed before counting unless `--keep-blacklisted-intervals` is set.

## Fragment And Midpoint Geometry

- Paired fragments span `[forward.pos, reverse.reference_end)`.
- Unpaired read-as-fragment mode spans `[read.pos, read.reference_end)`.
- Length bins are half-open and also define the accepted fragment length range.
- Even-length fragment midpoints use deterministic coordinate-derived random rounding. Duplicate fragments with the same coordinates choose the same midpoint base.
- Tile ownership is by midpoint position in the tile core, not by fragment start.
- Blacklist filtering is fragment-level and happens after midpoint tile ownership is checked so boundary fragments do not inflate blacklist statistics in neighboring tiles.

## Counting

- Each accepted fragment contributes to every grouped site window containing its midpoint.
- Forward and unstranded intervals use zero-based genomic offset from the site start.
- Reverse intervals mirror the site position, so the rightmost base in the half-open interval maps to profile position 0.
- The dense public array has shape `(group, length_bin, position)`.
- The flat storage order is group-major, then length-bin, then position.
- Fragment length lookup must use the shared `LengthAxis`, not a manual bin search.
- Weight per contribution is GC weight times optional scaling weight.
- Scaling is averaged over the full aligned fragment and applied to every selected midpoint window.
- Profiles are smoothed by default with order-3 Savitzky-Golay smoothing and a 165 bp window. `--smoothing none` disables smoothing. `--smoothing savgol=<odd_bp>` sets an explicit odd window size in base pairs.
- Smoothing flanks are derived from the selected smoothing mode: `--smoothing none` uses flank 0, and Savitzky-Golay profiles use `floor(window_bp / 2)`.
- Smoothing flanks are computation-only. When smoothing is enabled, counting uses flanked intervals, but output positions still refer to the original interval.
- Smoothing flanks must fit chromosome bounds. The command fails instead of edge-padding or clipping smoothed intervals.

## GC, Scaling, And Blacklists

- GC correction can come from `--gc-file` plus `--ref-2bit`, or from a two-byte BAM AUX `--gc-tag`.
- Invalid GC weights skip fragments by default. `--neutralize-invalid-gc` keeps them with weight 1.0.
- File-based GC correction reads reference sequence for the full tile fetch span.
- Scaling TSVs must fully cover every selected chromosome and pass GC-mode compatibility checks.
- Blacklist strategy uses the aligned fragment span and the shared blacklist strategies: `any`, `all`, `midpoint`, or `proportion=<threshold>`.
- `midpoint` blacklist strategy checks the single central base for odd fragments and either central base for even fragments. It does not use the randomized counted-midpoint tie break.
- `--keep-blacklisted-intervals` disables only interval-level blacklist prefiltering. Fragment-level blacklist filtering still applies.

## Tiling And Merge

- The command first selects grouped site windows overlapping a tile core, then narrows BAM fetch to the extreme selected site span plus maximum-fragment-length halo, clamped to tile fetch.
- Tile workers count into `SparseProfileGroupsCounts`.
- Non-empty tile partials are written as sparse `.npz` archives with arrays:
  - `idx.npy`: sorted flat dense indices as `u64`
  - `data.npy`: `f32` counts
  - `shape.npy`: `[group, length_bin, position]`
- Empty sparse tiles do not write partial files.
- Merge allocates one final dense output buffer and adds sparse partials in parallel into chunk locks.
- Sparse partials must validate shape, sorted indices, equal `idx`/`data` lengths, platform index fit, and destination bounds before merging.
- After merge, configured order-3 Savitzky-Golay smoothing runs on final grouped profiles along the position axis unless `--smoothing none` is used.
- Final binning averages adjacent positions after smoothing and flank trimming. `--bin-size 1` preserves base resolution.
- A shorter final position bin is divided by its actual width, not by `--bin-size`.

## Output Contract

- Main output is `<prefix>.midpoint_profiles.npy` with shape `(group, length_bin, position)`.
- The command writes one selected profile output. Users who want multiple transforms, such as unsmoothed and smoothed profiles, should run the command with separate output prefixes.
- Group index output is `<prefix>.group_index.tsv` with columns `group_idx`, `group_name`, and `eligible_intervals`.
- `eligible_intervals` is the number of intervals retained in each group after chromosome filtering and interval-level blacklist prefiltering. It is independent of fragment overlap, so an interval still counts when no fragment midpoint lands inside it.
- Settings output is `<prefix>.midpoint_profile_settings.json`. It records array axes, length-bin edges, output interval length, counted interval length, final position bin size, bin aggregation, last bin width, smoothing method, smoothing window, smoothing order, computation flank, correction flags, and interval blacklist prefilter state and margin.
- With the `plotters` feature, selected group indices can emit quick QC line plots and length-bin heatmaps.
- Run statistics include counted fragments, intervals after chromosome filtering, blacklist-prefiltered intervals, intervals retained for counting, blacklist exclusions, and GC failure summaries when relevant. When interval blacklist prefiltering is active, the reported run statistics also include the prefilter margin.

## Open Notes

! Warning: Projection from paired sites to midpoint-centered profiles remains future work. Keep projection plans outside finalized specs until implemented.

! Warning: File-based GC prefix construction still happens before no-window pruning for sparse targeted runs. This is a performance issue, not a correctness mismatch.
