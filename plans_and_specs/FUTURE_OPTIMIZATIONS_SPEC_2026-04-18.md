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

### Package-wide

- Shorten buffer and temporary-state lifetimes where they are obviously longer than needed
- Narrow fetch or reference spans before loading heavy auxiliary state when that does not change semantics
