# Testing And API Guidelines

This repo is released as a CLI crate first. The public Rust API should stay small and intentional, even if some internal pieces may become reusable later.

The goal of this guide is simple: do not let test layout force a bad public API.

## Test placement

Use `tests/` for behavior that a downstream user can rely on:

- Public crate API
- CLI commands and output files
- End-to-end workflows
- Regression tests that exercise released behavior

Use module-local test files for internal logic that should stay private or `pub(crate)`:

- Parsing helpers
- Numerical kernels
- Windowing helpers
- Internal state machines
- Data reshaping helpers

The preferred pattern is a sibling `*_tests.rs` file included from the module that owns the private code.

Example:

```rust
#[cfg(test)]
mod tests {
    include!("interpolation_tests.rs");
}
```

This keeps tests in separate files without forcing internal items to become `pub`.

## Visibility rules

Choose visibility based on support intent, not on test convenience.

- Use private visibility for implementation details that are local to one module
- Use `pub(crate)` for helpers shared within the crate
- Use `pub` only for types and functions that are part of the intended external API

Before marking something `pub`, ask:

- Would we be comfortable supporting this name and behavior in a future release?
- Does this represent a real workflow, config type, or result type that another crate should use?
- Would exposing this make later cleanup harder?

If the answer is unclear, keep it private or `pub(crate)`.

## Top-level API policy

Do not export modules just because they exist.

The crate root should expose only a curated surface:

- High-level workflows that are genuinely reusable
- Stable config types that callers are expected to construct
- Stable result types that callers are expected to inspect

Do not expose:

- Command plumbing
- File naming helpers
- Progress reporting
- Intermediate structs that only reflect the current implementation
- Low-level helpers that are likely to change during cleanup

## Testing boundaries

Use the strongest boundary that still tests the right behavior.

- If a behavior is part of the public CLI or crate API, test it from `tests/`
- If a behavior is internal but important, test it in a module-local `*_tests.rs` file
- If an integration test needs broad internal access, the test is probably at the wrong layer

Do not make items `pub` just so a test in `tests/` can reach them.

## Migration direction

As the first release is CLI-focused, cleanup should move in this order:

1. Decide the small external API we actually want to support
2. Move internal tests behind module-local `include!` test files where needed
3. Tighten visibility from `pub` to `pub(crate)` or private
4. Update `tests/` to cover only public behavior
5. Re-export only the stable entry points from the crate root

This keeps the CLI shippable now and leaves room for a sane reusable Rust API later.
