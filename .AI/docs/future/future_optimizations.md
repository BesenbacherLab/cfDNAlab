# Future Optimizations Spec

## Summary

This is the package-level future-optimizations list.

Keep it short.

Use it to track ideas worth revisiting later, with links to the relevant command-level specs when
those already exist.

## Future Ideas

### fcoverage

- Rolling local `f64` delta accumulation
  - Detail: `plans_and_specs/fcoverage_accumulation_numerics_spec_2026-04-16.md`
- Rolling GC prefix span loading in e.g. 1Mb sub tiles
- Consider `--by-grouped-bed` for aggregated outputs only (`average` and `total`) rather than
  positional per-window outputs in the first pass
  - Related detail: `plans_and_specs/GROUPED_BED_DISTRIBUTIONS_SPEC_2026-04-18.md`

### wps

- Consider `--by-grouped-bed` for aggregated outputs first, if grouped window semantics turn out
  to be useful there as well
  - Related detail: `plans_and_specs/GROUPED_BED_DISTRIBUTIONS_SPEC_2026-04-18.md`

### fragment-kmers

- Revisit grouped BED support later
  - Current guess: fragment-specific positional selection may still make grouped aggregation
    tractable, but it needs its own output-semantics design rather than being bundled into the
    first grouped-distribution work
  - Related detail: `plans_and_specs/GROUPED_BED_DISTRIBUTIONS_SPEC_2026-04-18.md`

### wps-peaks

- Do not assume grouped BED support makes sense here without a separate design pass
  - Indexed peak outputs and peak-stat outputs likely need different grouped semantics

### Package-wide

- Shorten buffer and temporary-state lifetimes where they are obviously longer than needed
- Narrow fetch or reference spans before loading heavy auxiliary state when that does not change semantics
