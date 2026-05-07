# lengths Spec

`cfdna lengths` counts fragment lengths into a two-dimensional matrix. Rows are the selected output units, and columns are fragment length bins.

## Output Contract

- Main output is a `.npy` matrix with shape `(row, length_bin)`.
- Rows represent one of: global output, fixed/BED windows, or grouped BED groups.
- The length axis is shared by all rows and is described in the settings sidecar.
- Window metadata sidecars must preserve row interpretation for fixed/BED/grouped outputs.
- Settings JSON records length-axis intervals, aggregation level, window mode, indel mode, clip mode, assignment mode, GC settings, and whether scaling factors were used.

## Length Axis

- Length bins are half-open: `[start, end)`.
- CLI range syntax `start:end:step` creates contiguous bins ending before `end`.
- Explicit edges are allowed and must be strictly increasing.
- The first edge must be at least 10 bp, matching the shared minimum ACGT span needed for GC fraction logic.
- The final edge must stay within `MAX_SUPPORTED_FRAGMENT_LENGTH + 1`.
- The writer marks a bin definition as `stepped_range` when widths are regular, allowing a shorter final bin. Otherwise it records explicit edges.
- `LengthAxis` owns an O(1) lookup table from raw fragment length to bin index. Counting code should use this axis rather than searching bins manually.
- `LengthCounts` merge and stack operations must reject incompatible axes.

## Fragment Length Semantics

- Paired fragments use `reverse.reference_end - forward.pos`.
- Unpaired read-as-fragment mode uses `read.reference_end - read.pos`.
- Fragment length filtering is applied to the command's effective length after the selected indel and clip modes.
- Default CLI filtering is 30..1000 bp inclusive. Test fixtures must never use fragment lengths below the shared 10 bp minimum.

## Indels

- `ignore`: keep the aligned span length.
- `skip`: discard fragments affected by relevant indels.
- `adjust`: subtract deletion bases where both mates show a deletion at mate-overlap positions and add the shortest insertion per position.
- `adjust` cannot repair mate gaps. Gaps between mates remain outside the adjusted indel accounting.
- `max_deletion_bases` bounds per-fragment deletion accounting and defaults to 100, with a hard maximum of 256.

## Clipping

- `aligned`: use aligned reference span length.
- `skip`: discard fragments with relevant soft clipping.
- `adjust`: include terminal soft-clipped bases in the effective length.
- Hard-clipped fragments are not recovered as observed sequence.
- GC correction, scaling, and blacklist checks continue to use aligned reference geometry. For clipped-only `count-overlap` contributions, scaling uses the nearest aligned reference base.

## Window Assignment

- Global mode has one output row.
- Fixed-size and BED modes assign fragments by the configured shared `WindowAssigner`.
- `count-overlap` contributes a fractional count to each overlapped window.
- `any`, `all`, `midpoint`, and `proportion=<threshold>` select windows by fragment-position overlap and then count full fragment mass in selected rows.
- Even-length midpoint assignment uses deterministic coordinate-derived random rounding. Duplicate fragments with the same coordinates choose the same midpoint base.
- Grouped BED mode collapses windows with the same group name into one row.

## Blacklist, GC, And Scaling

- Blacklist strategy is fragment-level and uses aligned fragment geometry.
- File-based GC correction uses aligned reference span GC, independent of indel or clip adjustment.
- GC correction for lengths can average a correction package over a configured length range and weighting strategy (`equal`, `frequency`, or `max-frequency`).
- `--gc-length-trim-rare` trims rare length bins from the GC-length marginalization step.
- Scaling factors are averaged over aligned fragment bases and must fully cover every selected chromosome.

## Implementation Invariants

- Counting code should keep the matrix dense in final output and use sparse or tiled internals only as implementation detail.
- All per-window and per-group merges must preserve length-axis identity.
- Settings sidecars are part of the output contract. Add fields when they are necessary for reproducibility, but avoid duplicating values already obvious from filenames or matrix shape.
