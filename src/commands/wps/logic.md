# Windowed Protection Score (WPS) Implementation Notes

### Windowed protection score semantics

The original Snyder et al. definition treats the two WPS components as disjoint groups. When `cfdna wps` evaluates a window centred at position `c`, a fragment:

- Contributes to the **protected** count when the window start lies at or to the right of the fragment start *and* the window end lies at or to the left of the fragment end. Exact edge alignment therefore still counts as "fully covered".
- Contributes to the **end-within** count only when an endpoint falls strictly inside the window (respecting the half-open convention we use for fragments). Endpoints that coincide with the window boundary are ignored.

This ensures a fragment can never be counted both as "fully spanning" and "ending inside" the same window, mirroring the reference implementation from the shendurelab `cfDNA` toolkit.

## Definition
- Let `window_size` be the user-configured span (default 120 bp).  
- For a genomic centre `c`, define the window as `[c - left_span, c + right_span)`, where:
  - `left_span = window_size / 2`
  - `right_span = window_size - left_span` (so even windows lean one base to the right).
- The per-position WPS value is
  ```
  WPS(c) = fragments_fully_covering_window(c) - fragment_ends_inside_window(c)
  ```
- A fragment fully covers the window when the window start is at or to the right of the fragment start *and* the window end is at or to the left of the fragment end.
- A fragment contributes to the "ends inside" term only when either endpoint sits strictly inside the window.

## Fragment Length Constraints
- Enforce `min_fragment_length >= window_size` to ensure every fragment can, in principle, span the full window.
- Enforce `window_size <= max_fragment_length` so fragments large enough to contribute exist.
- Keep `min_fragment_length` and `max_fragment_length` configurable, but default both bounds to include the canonical cfDNA band (e.g. 120-180 bp when `window_size = 120`).

## Tiling Strategy
- Process chromosomes in tiles (default 20 Mb) but extend each tile by a halo of `window_size` bases on both sides when gathering fragments.
- Accumulate WPS on the dilated span `[core_start - left_span, core_end + right_span)` and trim results back to the original tile core so each genomic position is written exactly once.

## Fragment Stream & Weights
- Pair reads into `Fragment` values using the minimal fragment iterator (indel-aware structures are unnecessary).
- Apply fragment-level filters already shared with `fcoverage` (orientation, mapping quality, length bounds, proper pairs optional).
- If genome scaling is enabled, multiply the final per-position WPS value by the scaling factor associated with that centre base. 
  This avoids calculating per-fragment weighting, which would only differ slightly around scaling bin edges.

## Difference-Array Accumulation
- For a tile covering the dilated reference coordinates `[tile_start, tile_end)`, allocate two `Vec<f32>` buffers with one entry per reference base:
  1. `overlap_diff`: tracks how many fragments fully span the window centred at each base.
  2. `end_diff`: tracks how many fragment ends fall inside the window centred at each base.
- Range updates use the standard "difference array" trick: to add weight `w` to every index `i` in `[a, b)`, increment `diff[a] += w` and decrement `diff[b] -= w`. A later prefix sum turns these markers into the per-base totals.
- For each fragment `[start, end)` and weight `w` (already clipped to the halo span):
- **Full-window span:** The window centred at `c` lies entirely inside the fragment whenever `start <= c - left_span` *and* `c + right_span <= end`. For integer centres this describes the inclusive range `[start + left_span, end - right_span]`, which we encode by adding `+w` at `start + left_span` and `-w` at `end - right_span + 1`.
- **Left endpoint contribution:** To count the left endpoint we require `c - left_span < start` (window start lies strictly before the endpoint) and `start < c + right_span` (endpoint sits strictly inside the upper bound). This is the open interval `(start - right_span, start + left_span)`, represented as `[start - right_span + 1, start + left_span)` in zero-based indices.
- **Right endpoint contribution:** The fragment interval is half-open, so the right endpoint sits at `end - 1`. We subtract when `end - right_span < c` and `c + left_span < end`, which maps to `(end - right_span, end + left_span - 1)` and becomes `[end - right_span, end + left_span - 1)` in the diff buffer.
- Illustration (window size 120 -> `left_span = right_span = 60`):

  ```text

  Fragment start = 100, fragment end = 250
  Window centred at c spans [c - 60, c + 60)

  Fully covered centres: c in [160, 191) (integers 160..190)
      overlap_diff[160] += w
      overlap_diff[191] -= w

  Left endpoint p = 100 (centres 41..160)
      end_diff[41]  += w
      end_diff[160] -= w

  Right endpoint p = 249 (centres 190..309)
      end_diff[190] += w
      end_diff[309] -= w

  ```

- Why `[41, 160)` instead of `[40, 160)`? With half-open windows `[c - left_span, c + right_span)`, the upper bound is excluded. At `c = 40` the window is `[-20, 100)`, so the base at 100 falls just outside. Starting at `c = 41` gives `[-19, 101)`, which genuinely contains 100. The interval `[41, 160)` therefore lists the 119 centres whose windows strictly contain the left endpoint; the centre at `c = 160` is excluded because the endpoint now sits on the boundary and is treated as fully covered.

- After all fragments have been processed for the tile, take a prefix sum over each diff buffer to recover the per-centre counts, then compute `wps = overlap_counts - end_counts`.

## Blacklist Handling
- Build a base-level mask from blacklist intervals via shared helpers.
- Dilate the mask by `max_fragment_length + max(left_span, right_span)` so any centre whose window touches a blacklisted base is excluded.
- Do **not** drop fragments outright; instead, mark WPS values as `NaN` wherever the dilated mask is true before writing aggregates. Positional writers skip these masked bases rather than emitting literal `NaN` rows.

## Output
- Round values per the `decimals` setting.
- Respect `keep_zero_runs`: when false, simply skip runs where `wps` remains zero so the bedGraph contains only non-zero intervals.
- Reuse `fcoverage` window writers for `per_window` modes (`unique-positions`, `indexed-positions`, `average`, `total`) so behaviour matches across tools.
