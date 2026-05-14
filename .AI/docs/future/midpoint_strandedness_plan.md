# Midpoint Strandedness Plan

Status: future implementation plan, not current behavior.

Date: May 14th 2026.

## Goal

Add site-strand orientation to `midpoints` so grouped profiles can represent
directional interval sets such as TSS windows.

The output profile should stay in genomic coordinate order for unstranded and
forward-strand intervals. Reverse-strand intervals should be mirrored when a
midpoint is written into the group profile.

```text
+ or .: position = midpoint - window.start
-:     position = window.end - 1 - midpoint
```

The current helper in `src/commands/midpoints/strand.rs` already implements this
mapping. The remaining work is loading, carrying, and applying the strand values.

## Input Semantics

The midpoint interval file keeps the first four columns as they are today:

```text
chromosome  start  end  group_name
```

Optional strand values use BED/UCSC-style text tokens only:

```text
+  forward
-  reverse
.  unstranded
```

Do not support numeric strand tokens such as `1`, `-1`, or `0` in this command.
Those are ambiguous with score-like metadata and make column detection fragile.

## Strand Column Detection

Detection should inspect only the first 20 data rows. Data rows exclude blank
lines, comments, `track`, and `browser` directives.

The loader should choose one of three modes before parsing the full file:

- `Column6`:
    Use column 6 as strand.

- `Column5`:
    Use column 5 as strand for nonstandard grouped files that put strand
    directly after `group_name`.

- `Unstranded`:
    Ignore extra columns and treat every interval as unstranded.

Column 6 has priority because standard BED places strand in column 6:

```text
chrom  start  end  name  score  strand
```

Detection rules:

- If sampled rows with column 6 all have column 6 values in `+`, `-`, or `.`,
  use `Column6`.

- Else, if sampled rows with column 5 all have column 5 values in `+`, `-`, or
  `.`, and sampled rows do not look like standard BED6, use `Column5`.

- Else, use `Unstranded`.

The ambiguous case needs special handling:

- If sampled rows have 6 or more columns, column 6 is not a strand column, and
  column 5 looks like a strand column, log a warning and return an error. The
  message should say that strand-like values were found in column 5, but
  BED6-style files must put strand in column 6. This prevents silently ignoring
  strand information in a wide custom file.

For wide files where neither column 5 nor column 6 looks like strand, use
`Unstranded` and log at info or warn level that extra columns were ignored and no
strand column was detected. This keeps custom 20-column files usable.

## Full-File Parse Rules

After detection, parse the whole file in the chosen mode.

For `Column6`:

- Every retained row must have at least 6 columns.

- Column 6 must be `+`, `-`, or `.`.

- Any other value is an error.

For `Column5`:

- Every retained row must have at least 5 columns.

- Column 5 must be `+`, `-`, or `.`.

- Rows with 6 or more columns should error, because the chosen mode is the
  nonstandard 5-column strand format.

For `Unstranded`:

- Ignore all columns after `group_name`.

- Set every interval strand to `Strand::Unstranded` when strand parsing was
  requested.

The loader should log which mode was chosen and how many sampled rows informed
the decision.

## Shared Data Structure

Do not change `Interval` or `IndexedInterval`.

`IndexedInterval.idx()` already stores `group_idx` in grouped-window commands.
Adding strand there would mix two unrelated meanings into one field.

Instead, add optional strand metadata to `GroupedWindows`:

```rust
pub struct GroupedWindows {
    pub windows: Vec<IndexedInterval<u64>>,
    pub strands: Option<Vec<Strand>>,
    span: Span<i64>,
}
```

Use a shared strand type near the BED loader, not a midpoint-only type:

```rust
pub enum Strand {
    Forward,
    Reverse,
    Unstranded,
}
```

Invariant:

```text
if strands.is_some(), strands.as_ref().unwrap()[i] describes windows[i]
```

Existing commands that do not request strand metadata receive `strands: None`.
That keeps the current grouped BED behavior and avoids per-window overhead unless
a caller asks for strand parsing.

`GroupedWindows::new(...)` should receive the optional strand vector directly.
It does not need wrappers or a second constructor at call sites:

```rust
GroupedWindows::new(windows, None)
GroupedWindows::new(windows, Some(strands))
```

The unstranded form stores `strands: None`. The stranded form stores
`strands: Some(...)` and sorts windows and strands together. Existing unstranded
tests and call sites should pass `None` explicitly.

## Loader

Extend `load_grouped_windows_from_bed` with a direct boolean argument:

```rust
load_grouped_windows_from_bed(..., read_strands: bool)
```

Call sites should pass `false` when they want the current behavior and `true`
when they want strand metadata parsed into `GroupedWindows.strands`.

When `read_strands` is `false`:

- Ignore all columns after `group_name`, as today.

- Return `GroupedWindows { strands: None, ... }`.

When `read_strands` is `true`:

- Preserve the current group-name enumeration behavior.

- Filter chromosomes like the existing loader.

- Parse and validate the first four columns as today.

- Detect and parse strand according to the rules above.

- Sort windows and strands together by `(start, end)`.

This avoids copying the grouped BED loader while keeping strand handling explicit
at the call site.

## Window Preparation

Change `prepare_count_windows` to consume and return grouped windows with
optional strands:

```rust
FxHashMap<String, GroupedWindows>
```

It should keep the current in-place compaction pattern, but write the strand
value in lockstep when `strands` is present:

```rust
windows[write_idx] = expanded_window;
if let Some(strands) = strands.as_mut() {
    strands[write_idx] = strands[window_idx];
}
write_idx += 1;
```

Then truncate `windows` to `write_idx`. Also truncate `strands` when present.

The current checks for chromosome bounds, blacklist prefiltering, smoothing
flank expansion, and retained interval stats should remain unchanged.

## Tiling And Counting

Tile precomputation still uses only the window slice:

```rust
precompute_tile_window_spans(..., grouped_windows.as_slice(), ...)
```

Change tile narrowing to return tile-local grouped windows:

```rust
GroupedWindows
Interval<u64> fetch_span
```

The existing loop over `overlapping_windows_for_tile(...)` can push the window
and, when present, its strand by using the same source index.

Counting then becomes:

```rust
let window = core_overlapping_windows.windows[overlapped_window_idx];
let strand = core_overlapping_windows
    .strands
    .as_ref()
    .map(|strands| strands[overlapped_window_idx])
    .unwrap_or(Strand::Unstranded);
let window_position = stranded_window_position(window.interval, midpoint_u64, strand)?;
```

This should replace the current `midpoint - window.start()` calculation in both
the scaled and unscaled counting paths.

`find_overlapping_windows` can stay unchanged. Its `OverlappingWindow.idx` is
already the scan index into the supplied tile-local window slice.

## Tests

Keep the existing helper tests for `stranded_window_position`.

Add loader tests:

- Four-column grouped input becomes unstranded.

- Five-column input with `+`, `-`, `.` uses column 5.

- Standard BED6-style input uses column 6 and ignores score.

- Wide input with non-strand columns 5 and 6 is accepted as unstranded.

- Wide input with strand-like column 5 but non-strand column 6 errors with a
  message explaining that BED6-style strand belongs in column 6.

- Chosen strand mode rejects later rows with invalid strand values.

Add preparation tests:

- Blacklist compaction drops the matching strand entry with the dropped window.

- Smoothing expansion preserves the strand entry for retained windows.

Add counting-level unit coverage if practical:

- The same midpoint pattern in identical `+` and `-` windows writes mirrored
  group-profile arrays.

## Logging

Log the detected strand mode once per load:

```text
Detected midpoint interval strand column: column 6
Detected midpoint interval strand column: column 5
No midpoint interval strand column detected; treating intervals as unstranded
```

For wide custom files where columns 5 and 6 are not strand-like, warn or log at
info level that extra columns were ignored. Avoid warning on ordinary 4-column
or 5-column unstranded inputs.
