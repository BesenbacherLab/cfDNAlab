# Public Testing Fixtures Plan

This plan covers a possible `cfdnalab::testing` API for reusable temporary-input builders and testing-only helper utilities.

## Problem

Current integration tests keep useful fixture builders in `tests/fixtures/mod.rs`. That works for public integration tests, but it does not work for module-local internal tests because `tests/fixtures` is outside the crate. Moving private-behavior tests into sibling `*_tests.rs` files therefore forces either local fixture duplication or visibility widening.

At the same time, some of these helpers would be useful to downstream users writing tests around cfDNAlab workflows. A public testing module can be reasonable, but only if it is treated as a real API surface rather than a shortcut for exposing private internals.

## Decision To Explore

Add a feature-gated public module:

```rust
#[cfg(feature = "testing")]
pub mod testing;
```

Downstream crates would opt in explicitly:

```toml
cfdnalab = { version = "...", features = ["testing"] }
```

The feature should not be enabled by default. It should expose documented fixture builders and testing utilities, not production command internals.

Integration tests that import `cfdnalab::testing` are external users of the
crate, so they must be gated with:

```rust
#![cfg(feature = "testing")]
```

Those tests run under `cargo test --features testing` or `cargo test
--all-features`, not under plain `cargo test`. This keeps the testing API
opt-in for downstream users while still allowing shared fixtures instead of
duplicating `tests/fixtures` code.

## Scope

The testing module may expose:

- Small synthetic BAM builders.
- Small two-bit reference builders.
- BED and scaling-factor file writers.
- Public-output artifact readers used in tests.
- Testing wrappers around package writers and loaders when the production API should stay narrower.

The testing module should not expose:

- Current command plumbing just because a test needs it.
- Private loader structs under their current implementation names.
- Helpers whose behavior depends on current temp-file orchestration.
- Broad command config constructors unless they are generalized into documented fixture scenarios.

## Raw Inventory

Current `tests/fixtures/mod.rs` helpers that contain reusable fixture logic:

- `BamFixture`
- `ReadSpec`
- `FragmentSpec`
- `paired_fragment`
- `bam_from_specs`
- `bam_from_specs_strict_identity`
- `bam_from_fragment_starts`
- `simple_inward_bam`
- `complex_bam_fixture`
- `long_fragment_bam`
- `single_read_bam_with_qualities`
- `TwoBitFixture`
- `twobit_from_sequences`
- `simple_reference_twobit`
- `complex_reference_twobit`
- `write_bed`
- `write_scaling_factors`
- `read_zst_to_string`
- `read_length_counts_text`
- `read_length_counts_tsv`
- `read_midpoint_zarr_counts`
- `read_midpoint_zarr_i32_1d`
- `read_midpoint_zarr_u32_1d`
- `touch_file`

Helpers that need heavier review before any public exposure:

- `build_real_neutral_gc_package`
- `build_real_neutral_gc_package_for_range`
- `build_real_non_neutral_gc_package`
- `write_constant_gc_package`
- `write_two_bin_gc_package`
- command-selection builders such as `single_position_selection` and `build_base_selection`

These names are not proposed public API names. They are an inventory of reusable code that must be generalized, renamed, documented, or rejected before becoming public.

Names such as `simple_inward_bam`, `complex_bam_fixture`, `simple_reference_twobit`, and `complex_reference_twobit` are test-local names, not good public names. They do not say enough about the fixture contract: contigs, fragment spans, read lengths, mapping qualities, CIGAR contents, pairing assumptions, or which behavior the fixture is meant to exercise.

## Generalization Work

### 1. Split Fixture Categories

Create separate submodules so the public testing surface is navigable:

```text
src/testing/
  mod.rs
  bam.rs
  reference.rs
  bed.rs
  scaling.rs
  output_readers.rs
  gc_packages.rs
```

Possible public paths:

```rust
cfdnalab::testing::bam::bam_from_fragments
cfdnalab::testing::reference::twobit_from_sequences
cfdnalab::testing::gc_packages::write_constant_gc_correction_package
```

Consider top-level re-exports only for the most commonly used fixtures.

### 2. Prefer General Builders Over Named Scenarios

The main public surface should be parameterized builders with documented defaults, not a collection of vague named fixtures.

Preferred shape:

```rust
let bam = TempBamBuilder::new()
    .contig("chr1", 200)
    .paired_fragment(PairedFragmentSpec {
        tid: 0,
        start: 20,
        fragment_length: 60,
        read_length: 20,
        mapq: 60,
    })
    .build()?;
```

Keep convenience constructors only when the scenario is both common and precisely named.

Acceptable examples:

- `single_contig_inward_pair_bam`
- `bam_with_stacked_fragments`
- `bam_with_indel_and_softclip_reads`
- `twobit_with_single_repeating_contig`
- `twobit_from_sequences`

Avoid names that describe complexity rather than contract:

- `simple_*`
- `complex_*`
- `default_*` unless every default is documented and intentionally stable
- command-specific names for general BAM/reference data

If a current named fixture is kept, rewrite it as a thin wrapper around a documented builder and document the exact coordinates and fragment spans it creates.

### 3. Document Ownership And Lifetime

Returned temporary-input structs must document that they own a temporary directory and keep generated paths valid for as long as the value is alive.

Example responsibilities:

- `TempBam` owns the BAM and BAI paths.
- `TempTwoBit` owns the two-bit path and stores the source sequences for assertions.
- Returned paths are invalid once the owning value is dropped.

### 4. Validate Fixture Inputs

Generalized fixture builders should fail clearly instead of silently creating invalid genomics fixtures.

Required checks include:

- Fragment length must be at least 10 bp.
- Read length must be positive.
- Paired fragment coordinates must keep the reverse read within the fragment span.
- BAM contig lengths must cover all read reference spans.
- Sequence strings should be uppercase-normalized and restricted to supported bases unless the helper explicitly documents ambiguous bases.
- Duplicate paired-read qnames should remain an error in the ordinary BAM builder.

### 5. Rename Implementation-Shaped Helpers

Public names should say what testing contract they provide.

Examples:

- `simple_inward_bam` -> builder call or `single_contig_inward_pair_bam`
- `complex_bam_fixture` -> split into narrower fixtures such as `bam_with_indel_and_softclip_reads` and `bam_with_cross_contig_mate_records`
- `simple_reference_twobit` -> `twobit_with_single_repeating_contig`
- `complex_reference_twobit` -> `twobit_with_repeating_contigs`
- `write_constant_gc_package` -> `write_constant_gc_correction_package` or `write_unit_gc_correction_package`
- `write_two_bin_gc_package` -> `write_two_bin_gc_correction_package`
- `build_real_neutral_gc_package` -> `build_command_produced_gc_correction_package_for_length`
- `build_real_neutral_gc_package_for_range` -> `build_command_produced_gc_correction_package_for_range`
- `build_real_non_neutral_gc_package` -> `build_command_produced_gc_correction_package_from_reference_windows`

The exact names can change during implementation, but the names should distinguish hand-authored fixture packages from artifacts created by running cfDNAlab command code.

### 6. Wrap Private Helpers Deliberately

If tests need access to useful non-public functionality, expose testing wrappers rather than production names.

For example, do not make the current `load_reference_gc_data` public just to support tests. Prefer a testing helper such as:

```rust
cfdnalab::testing::gc_packages::load_reference_gc_package_for_test(path)
```

That wrapper can return a documented test-facing summary type with only assertion-relevant fields. It should not leak internal structs like `ReferenceGCData` or `ReferenceGCMetadata` unless those types are intentionally redesigned as public package-inspection types.

### 7. Keep Feature Dependencies Explicit

Some fixture helpers require command features.

Examples:

- BAM and two-bit builders can likely live under `testing` alone.
- GC correction package writers require `cmd_gc_bias`.
- Command-produced GC correction package fixtures require both `cmd_gc_bias` and `cmd_ref_gc_bias`.
- Fragment-kmer fixtures require `cmd_fragment_kmers` only if they call command-specific APIs.

Use `#[cfg(...)]` at the helper or submodule level and document required feature combinations.

### 8. Migrate Existing Tests In Stages

Recommended migration order:

1. Add `testing` feature and introduce documented BAM, two-bit, BED, scaling, and output-reader builders under public names.
2. Reimplement current named fixtures as private integration-test wrappers over the new public builders.
3. Add `#![cfg(feature = "testing")]` to integration test files that import `cfdnalab::testing`.
4. Update integration tests to import stable builders from `cfdnalab::testing` where appropriate.
5. Update module-local tests that need general fixtures to use `crate::testing` under `#[cfg(test)]`; the crate-local module can be compiled with `#[cfg(any(feature = "testing", test))]`.
6. Generalize and document GC package fixture helpers.
7. Replace direct private loader usage in tests with module-local tests or explicit testing wrappers.
8. Remove or keep private the old vague named wrappers once no test needs them.

## Documentation Requirements

Each public fixture helper should document:

- What files or objects it creates.
- Whether the fixture is synthetic, hand-authored, or command-produced.
- What assumptions are built into the fixture, including contig names, fragment lengths, read lengths, mapping quality, orientation, and reference sequence.
- Which cargo features are required.
- Whether expected values can be derived directly from the helper's documented coordinates and fragment spans.

The module-level docs should say that `cfdnalab::testing` is for writing tests and examples. It is public and supported, but it is not the recommended production API for running analyses.

## Open Questions

- Should `testing` be included in the published crate package, or should this become a separate `cfdnalab-testing` crate later?
- Should command-produced fixture builders run command code directly, or should they stay in integration tests because they are heavier than ordinary fixture builders?
- Should public output readers in `testing` overlap with future Rust package-loader APIs, or stay explicitly test-oriented?
- Should reference GC package inspection become a production API, or only a testing wrapper?

## Acceptance Criteria

- `cargo check --features testing` succeeds.
- `cargo check --tests --features testing` succeeds after migrated tests are updated.
- Every integration test file that imports `cfdnalab::testing` has `#![cfg(feature = "testing")]`.
- Public testing helpers have module docs and per-helper docs for their fixture contract.
- No production-private type is made public solely to satisfy integration tests.
- Integration tests no longer depend on `tests/fixtures/mod.rs`.
- Module-local tests can reuse shared BAM, two-bit, BED, scaling, and package fixtures without duplicating large helper blocks.
