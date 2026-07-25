# tiling Spec

Tiling provides parallelism and bounded memory. Command correctness depends on separating tile core ownership from fetch halos and reducer keys.

## Tile Model

- A `Tile` has a chromosome, BAM tid, zero-based tile index, core interval, and fetch interval.
- Core and fetch intervals are checked half-open coordinates.
- Fetch must fully contain core.
- Cores are non-overlapping and cover each selected chromosome from 0 to chromosome length.
- Fetch halos are clamped to chromosome bounds.
- Tile BAM tid must match the BAM reader tid. Negative tile tids are rejected before comparison.

## Building Tiles

- `build_tiles` creates tiles in selected chromosome order.
- `tile_bp` is the desired core width.
- `halo_bp` is added on both sides of each core for BAM/reference fetches.
- If `align_bp` is supplied and feasible, core width is rounded down to preserve fixed-window alignment.
- Alignment is guaranteed when `tile_bp` is already divisible by `align_bp`, or when at least ten aligned bins fit and rounding down is possible.
- Halos are not aligned. Only core starts and non-final core ends participate in alignment checks.

## Ownership Rules

Commands must state their tile ownership rule explicitly:

- Coverage and lengths-style fragment commands usually own a fragment when its aligned start lies in the tile core.
- Midpoint profiles own a fragment when its midpoint lies in the tile core.
- Positional coverage owns positions by tile core coordinates.

Do not infer ownership from fetch. Fetch is only the superset needed to reconstruct objects near core boundaries.

## Window Spans

- `TileWindowSpan` is a half-open index range into a start-sorted per-chromosome window slice.
- Empty spans have `first_idx == last_idx_exclusive`.
- `precompute_tile_window_spans` streams through tiles and windows with monotonic pointers.
- Left and right halos in span precomputation describe command-specific reach, not fetch padding.
- Cached spans are candidate sets. Some commands use them directly, while core-overlap commands filter true overlap with `overlapping_windows_for_tile`.

## Fetch Narrowing

Fetch narrowing must use the shared `BedFetchPolicy` vocabulary:

- `CoreOverlap`: BED windows are relevant if they overlap the tile core.
- `CandidateWindowExtent`: the caller has already selected candidate windows and provides a halo large enough to preserve needed aligned reads.
- `KeepTileFetch`: BED coordinates are not safe for aligned fetch narrowing, so keep the full tile fetch.

For fixed-size windows, fetch narrowing derives the first and last fixed bins touching the tile core and clamps that span with the requested halo.

## Temporary Directories

- Commands with temporary files create `tmp.<command-and-prefix>.work.<random>` under the directory selected with `--temp-dir` and `tmp.<command-and-prefix>.final.<random>` under the output directory.
- Both directories are visible, atomically created, unique per run, and owned by cleanup guards.
- Ctrl-C cleanup registration is shared across guarded temp dirs. The signal handler starts the hidden cleanup helper and exits the main process, so cleanup happens after in-process worker threads have stopped writing.
- Temporary files may be placed on node-local scratch with `--temp-dir`. Files awaiting their final rename remain under the selected output directory.

## Reducer Rules

- Reducers must consume explicit tile paths returned by workers.
- Reducers must not scan temp directories to discover work.
- Cross-tile metadata files should list only boundary-crossing rows or windows.
- Missing cross-index entries mean the row had one contribution.
- Fixed-size aggregate reducers key by full fixed-bin start before final chromosome clipping.
- BED aggregate reducers key by original BED index.
- Grouped reducers key by `group_idx`.
- Positional reducers merge by chromosome and tile index to preserve coordinate order.

## Fixed-Window Fast Path

- When tile cores align with fixed-size windows and no global post-scaling is needed, commands may concatenate per-tile final rows.
- Restore-mean length normalization and other global multipliers must force a reducer path.
- The final fixed-size bin is clipped to chromosome length in reducer output.

## Failure Policy

- Empty no-window tiles may skip worker output.
- A missing output for a chromosome with expected windows is an error.
- Duplicate reducer keys are errors unless the reducer explicitly defines summation for that key.
