# Public API plan

This plan captures the initial public Rust API direction for the crate.

---

Before marking something `pub`, ask:

- Would we be comfortable supporting this name and behavior in a future release?
- Does this represent a real workflow, config type, or result type that another crate should use?
- Would exposing this make later cleanup harder?

---

## Decision summary

This top section is the actual decision layer. The detailed method inventory lower down is reference material.

### Public now or soon

- `GC correction`
  - Goal: let downstream Rust code load correction packages, build GC prefixes, and query multiplicative correction weights
- `Blacklisting`
  - Goal: let downstream Rust code load blacklist intervals, apply overlap strategies, and mask reference slices

### Public after cleanup

- `Intervals`
  - Goal: expose a small geometric core for half-open intervals, overlap math, and interval merging
  - Needs cleanup because current overlap code mixes general interval logic with command-window assignment logic
- `Genomic scaling`
  - Goal: expose validated scaling bins and scaling application helpers
  - Needs cleanup because the current function names are implementation-shaped and not yet a clear public vocabulary
- `Tiling`
  - Goal: expose tile planning and tile-window traversal helpers
  - Needs cleanup because the module currently mixes reusable tile logic with command orchestration and temp-file helpers

### Internal for alpha

- Command runners and command-level config wrappers
  - The CLI is the supported end-to-end interface for full analyses in the first release
- Temp-file orchestration helpers
  - These are tied to current command implementation details
- Fragment and iteration APIs beyond what is strictly needed internally
  - These likely have reusable value, but the shape is still too messy to publish confidently

### Structural direction

- A possible shared `Interval { start, end }` core is worth exploring
- That interval type should be a minimal field type for geometry only
- It should not replace domain structs like tiles, windows, peaks, or fragments

### Next implementation target

Pick one area and clean it up properly before widening visibility.

Recommended order:

1. `Intervals`
2. `GC correction`
3. `Blacklisting`
4. `Genomic scaling`
5. `Tiling`

## Current implementation mapping

This section lists the actual methods and types the current released commands depend on for GC correction, genomic smoothing, blacklisting, and tiling. These are the existing entrypoints to review when deciding what should become public.

### GC correction

Current commands load a correction package, build GC prefixes from a reference slice, and query multiplicative weights per fragment.

- `commands::gc_bias::correct::GCCorrector`
- `commands::gc_bias::correct::GCCorrector::from_package`
- `commands::gc_bias::correct::GCCorrector::correct_fragment`
- `commands::gc_bias::correct::GCCorrector::get_correction_weight`
- `commands::gc_bias::correct::LengthAgnosticGCCorrector`
- `commands::gc_bias::correct::LengthAgnosticGCCorrector::from_gc_corrector`
- `commands::gc_bias::correct::LengthAgnosticGCCorrector::correct_fragment`
- `commands::gc_bias::correct::LengthAgnosticGCCorrector::get_correction_weight`
- `commands::gc_bias::correct::MarginalizeLengthsWeightingScheme`
- `commands::gc_bias::correct::load_gc_corrector`
- `commands::gc_bias::correct::load_length_agnostic_gc_corrector`
- `commands::gc_bias::counting::GCPrefixes`
- `commands::gc_bias::counting::build_gc_prefixes`
- `commands::gc_bias::counting::get_gc_integer_percentage_for_window`

Likely command-level convenience only:

- `commands::cli_common::ApplyGCArgs`

### Genomic smoothing / scaling

Current commands either load a chromosome -> scaling-bin map, compute an average fragment/window scaling factor from overlapping bins, or apply per-base scaling in place to a positional signal.

- `shared::scale_genome::load_scaling_factors_tsv`
- `shared::scale_genome::apply_scaling_to_coverage_in_place`
- `shared::scale_genome::compute_per_window_scaling_over_overlap`
- `shared::scale_genome::compute_per_window_scaling_over_fragment`
- `shared::overlaps::OverlappingWindows`
- `shared::overlaps::OverlappingWindow`
- `shared::overlaps::find_overlapping_windows`
- `shared::overlaps::create_overlapping_bins_by_size`
- `shared::overlaps::half_open_intervals_overlap`
- `shared::overlaps::fraction_overlap_of_a`
- `shared::overlaps::overlap_len`

Likely command-level convenience only:

- `commands::cli_common::ScaleGenomeArgs`
- `commands::cli_common::load_scaling_map`

### Blacklisting

Current commands load merged blacklist intervals once per run, then either test fragment overlap during streaming or mask reference sequence slices before downstream counting.

- `shared::blacklist::BlacklistStrategy`
- `shared::blacklist::load_blacklists`
- `shared::blacklist::compute_blacklist_overlap`
- `shared::blacklist::is_blacklisted`
- `shared::blacklist::apply_blacklist_mask_to_seq`

Currently used internally but not re-exported from `shared::blacklist`:

- `shared::blacklist::load::merge_intervals`

Likely command-level convenience only:

- `commands::cli_common::load_blacklist_map`

### Tiling

Current tiled commands build tiles, cache overlapping window spans, narrow fetch ranges to the needed windows, and then iterate overlapping windows tile by tile.

Clearly reusable core:

- `shared::tiled_run::Tile`
- `shared::tiled_run::TileWindowSpan`
- `shared::tiled_run::TileWindowSpan::is_empty`
- `shared::tiled_run::build_tiles`
- `shared::tiled_run::precompute_tile_window_spans`
- `shared::tiled_run::overlapping_windows_for_tile`
- `shared::tiled_run::tile_window_min_max`
- `shared::tiled_run::clamp_fetch_to_window_span`
- `shared::tiled_run::windows_overlapping_core`

Currently used by command orchestration and temp-file merging, so these should be reviewed separately instead of exported by default:

- `shared::tiled_run::TileMode`
- `shared::tiled_run::parse_tile_index`
- `shared::tiled_run::make_temp_dir`

## Proposed public API shape

This section maps the current implementation methods to a cleaner public API structure.

The goal is not to freeze the current module names. The goal is to identify the concepts we want downstream users to work with.

### Possible shared interval core

One possible cleanup path is to introduce a small checked half-open `Interval { start, end }` type for the geometric part of existing domain structs.

This is not meant to replace domain structs like tiles, windows, peaks, or fragments.

Instead, it would replace repeated raw coordinate fields inside them, so the code can share:

- Half-open interval invariants
- Start/end validation in one place
- Overlap and merge helpers over one geometric type

Examples of the intended shape:

- `Tile { core: Interval<_>, fetch: Interval<_>, ... }`
- `IntermediateWindow { chrom: String, interval: Interval<_>, ... }`
- `Peak { interval: Interval<_>, ... }`
- `WindowOverlap { interval: Interval<_>, ... }`

This should stay minimal. It should not grow optional metadata fields like `idx`, since those belong to the domain structs that carry the interval.

### Intervals

This should be its own public category instead of being split across smoothing and blacklisting.

Proposed surface:

- Interval overlap predicates and overlap sizes
- Fixed-bin overlap enumeration
- Window-overlap collection for an interval
- Interval merging / coalescing

Maps from current code:

- Interval overlap predicates and sizes
  - `shared::overlaps::half_open_intervals_overlap`
  - `shared::overlaps::overlap_len`
  - `shared::overlaps::fraction_overlap_of_a`
- Fixed-bin overlap enumeration
  - `shared::overlaps::create_overlapping_bins_by_size`
- Window-overlap collection
  - `shared::overlaps::OverlappingWindow`
  - `shared::overlaps::OverlappingWindows`
  - `shared::overlaps::find_overlapping_windows`
- Interval merging / coalescing
  - `shared::blacklist::load::merge_intervals`

Notes:

- `merge_intervals` is conceptually general, even though it currently lives under blacklist loading
- If this becomes public, it likely belongs outside the blacklist-specific module path

### Genomic scaling

This should focus on scaling bins and scaling application, not on general interval logic.

Proposed surface:

- Load validated scaling bins
- Average scaling over an interval or overlap span
- Apply scaling to a positional signal in place

Maps from current code:

- Load validated scaling bins
  - `shared::scale_genome::load_scaling_factors_tsv`
- Apply scaling to a positional signal in place
  - `shared::scale_genome::apply_scaling_to_coverage_in_place`
- Average scaling over an interval or overlap span
  - `shared::scale_genome::compute_per_window_scaling_over_overlap`
  - `shared::scale_genome::compute_per_window_scaling_over_fragment`

Notes:

- The two `compute_window_scaling_*` functions likely need better public names before exposure
- `commands::cli_common::ScaleGenomeArgs` and `commands::cli_common::load_scaling_map` still look command-level, not utility-level

### GC correction

This should expose the pieces needed to turn a correction package and reference sequence into multiplicative fragment weights.

Proposed surface:

- Build GC prefixes from a reference sequence
- Load GC correction package helpers
- Query multiplicative GC correction weights
- Optional length-agnostic correction view

Maps from current code:

- Build GC prefixes from a reference sequence
  - `commands::gc_bias::counting::GCPrefixes`
  - `commands::gc_bias::counting::build_gc_prefixes`
  - `commands::gc_bias::counting::get_gc_integer_percentage_for_window`
- Load GC correction package helpers
  - `commands::gc_bias::correct::load_gc_corrector`
  - `commands::gc_bias::correct::load_length_agnostic_gc_corrector`
- Query multiplicative GC correction weights
  - `commands::gc_bias::correct::GCCorrector`
  - `commands::gc_bias::correct::GCCorrector::from_package`
  - `commands::gc_bias::correct::GCCorrector::correct_fragment`
  - `commands::gc_bias::correct::GCCorrector::get_correction_weight`
- Optional length-agnostic correction view
  - `commands::gc_bias::correct::LengthAgnosticGCCorrector`
  - `commands::gc_bias::correct::LengthAgnosticGCCorrector::from_gc_corrector`
  - `commands::gc_bias::correct::LengthAgnosticGCCorrector::correct_fragment`
  - `commands::gc_bias::correct::LengthAgnosticGCCorrector::get_correction_weight`
  - `commands::gc_bias::correct::MarginalizeLengthsWeightingScheme`

Notes:

- `commands::cli_common::ApplyGCArgs` still looks command-level, not utility-level

### Blacklisting

This should focus on loading, overlap decisions, and masking, while the general interval merging logic can move under `Intervals`.

Proposed surface:

- Load and merge blacklist intervals from BED-like files
- Query whether spans are blacklisted under a chosen strategy
- Compute blacklist overlap fractions
- Mask reference slices in place

Maps from current code:

- Load and merge blacklist intervals from BED-like files
  - `shared::blacklist::load_blacklists`
- Query whether spans are blacklisted under a chosen strategy
  - `shared::blacklist::BlacklistStrategy`
  - `shared::blacklist::is_blacklisted`
- Compute blacklist overlap fractions
  - `shared::blacklist::compute_blacklist_overlap`
- Mask reference slices in place
  - `shared::blacklist::apply_blacklist_mask_to_seq`

Notes:

- `commands::cli_common::load_blacklist_map` still looks like a command convenience wrapper

### Tiling

This should expose the reusable tile-planning and tile-window traversal pieces, while leaving temp-file orchestration separate.

Proposed surface:

- Tile descriptions and cached window spans
- Build tiles over contigs
- Precompute overlapping window spans per tile
- Iterate windows overlapping a tile core
- Tighten fetch ranges to the needed windows

Maps from current code:

- Tile descriptions and cached window spans
  - `shared::tiled_run::Tile`
  - `shared::tiled_run::TileWindowSpan`
  - `shared::tiled_run::TileWindowSpan::is_empty`
- Build tiles over contigs
  - `shared::tiled_run::build_tiles`
- Precompute overlapping window spans per tile
  - `shared::tiled_run::precompute_tile_window_spans`
- Iterate windows overlapping a tile core
  - `shared::tiled_run::overlapping_windows_for_tile`
  - `shared::tiled_run::windows_overlapping_core`
- Tighten fetch ranges to the needed windows
  - `shared::tiled_run::tile_window_min_max`
  - `shared::tiled_run::clamp_fetch_to_window_span`

Notes:

- `shared::tiled_run::TileMode`, `shared::tiled_run::parse_tile_index`, and `shared::tiled_run::make_temp_dir` still look tied to command orchestration and should be considered separately
 
