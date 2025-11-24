# Windowed Protection Score Overview

This note records how cfDNAlab computes the windowed protection score (WPS) and why the code
behaves the way it does. The goal is to keep the model identical to the Snyder et al. reference
implementation while explaining it in plain language.

## What WPS Measures

For each genomic center `c` we place a window of width `window_size` around that base:

- `left_span = window_size / 2`
- `right_span = window_size - left_span`
- Window coordinates: `[c - left_span, c + right_span)` (half-open).

The WPS at `c` is:

```
WPS(c) = protected_fragments(c) - endpoint_hits(c)
```

Fragments are never counted more than once for the same window:

1. **Protected fragment**  
   The window lies fully inside the fragment, including the case where the window edges match the
   fragment edges exactly.
2. **Endpoint hit**  
   One (or both) fragment endpoints falls strictly inside the window.  
   - If both fragment ends touch the window boundaries we keep the fragment in the protected term.  
   - If an endpoint lies on a boundary while the opposite side does not contain the window, the fragment
     is ignored for that center.
3. **Otherwise** the fragment does not affect that window.

These rules come directly from the Snyder Python script and guarantee that a fragment cannot be
both “protected” and an “endpoint hit” in the same window.

## Fragment Length Bounds

- Users configure `min_fragment_length` and `max_fragment_length`.
- We enforce `min_fragment_length >= window_size` so a valid fragment can always span the window.
- We also enforce `window_size <= max_fragment_length`.

## Tiling in Practice

Chromosomes are processed in tiles so we can stream large BAM files:

- Each tile has a *core* region that we will eventually write to disk.
- We extend that core with a halo of `max_fragment_length + window_size` on each side before
  [TODO: Language]
  reading fragments. This guarantees every window we later emit has seen all fragments that could
  influence it.
- After processing, we discard the halo and only write WPS values for the core.

## Counting with Difference Arrays

We maintain two difference arrays over the dilated tile:

1. `overlap_diff` – counts windows that are fully protected.
2. `end_diff` – counts windows that contain a fragment endpoint.

The helper `push_range(diff, start, end, weight)` adds `weight` to every center in `[start, end)`.
We then run a prefix sum to recover the actual per-center counts.

Given a fragment `[start, end)` (half-open, so the rightmost included base is `end - 1`):

- **Protected range**  
  A center is protected when `c - left_span >= start` and `c + right_span <= end`.  
  We encode this as the closed integer range `[start + left_span, end - right_span]`.  
  In diff-array form we call  
  `push_range(overlap_diff, start + left_span, end - right_span + 1, +1.0)`.

- **Left endpoint range**  
  The left endpoint is strictly inside when `c - left_span < start` and `start < c + right_span`.  
  Integer centers therefore satisfy `start - right_span + 1 <= c <= start + left_span - 1`.  
  We encode this as `push_range(end_diff, start - right_span + 1, start + left_span, +1.0)`.

- **Right endpoint range**  
  The right endpoint lives at `end - 1`. It is strictly inside when  
  `c - left_span < end - 1` and `end - 1 < c + right_span`.  
  Integer centers satisfy `end - right_span + 1 <= c <= end + left_span - 2`.  
  We encode this as `push_range(end_diff, end - right_span + 1, end + left_span - 1, +1.0)`.

The half-open intervals above ensure there is no overlap between the protected range and either
endpoint range, so the later subtraction `protected - endpoints` respects the intended logic.

## Blacklist and Valid Centers

We build a mask per tile that marks centers we must skip:

- Dilate each blacklist interval by `max_fragment_length + max(left_span, right_span)` on both sides
  before applying it. The dilation ensures that any window whose supporting fragments might cross a
  blacklisted base is marked invalid.
- Any center whose resulting window would touch a blacklisted base is masked out.
- Centers that would extend past the chromosome edges are masked out.

[TODO: Language]
Masked centres never contribute to the output. We suppress them rather than emitting `NaN`.

## Output Rules

- Values are rounded using the configured `decimals`.
- When `keep_zero_runs` is false we drop runs whose value is exactly zero so the bedGraph stays
  compact.
- Aggregated `--per-window` outputs reuse the shared writers from the fcoverage pipeline to keep the
  behaviours aligned.
