# Tile, window, and fetch execution workdoc

Date: 2026-04-01

## Purpose

This document constrains the next implementation steps so the work stays aligned with:

- [TILE_WINDOW_AND_FETCH_RULES_SPEC_2026-04-01.md](/Users/au547627/Documents/Development/rust/cfDNAlab/plans_and_specs/TILE_WINDOW_AND_FETCH_RULES_SPEC_2026-04-01.md)
- [TILE_WINDOW_AND_FETCH_TEST_SPEC_2026-04-01.md](/Users/au547627/Documents/Development/rust/cfDNAlab/plans_and_specs/TILE_WINDOW_AND_FETCH_TEST_SPEC_2026-04-01.md)

It exists because the recent drift pattern was:

- adding tests on helpers that are about to be decomposed
- silently treating migration guardrails as if they were target-contract tests
- replacing only part of the removed coverage at the correct boundary

This workdoc is the check against repeating that.

## Anti-drift rules

The next pass must obey all of the following:

1. Add no behavior changes.
2. Add only:
   - one working document
   - new helper stubs with `unimplemented!()`
   - tests required by the test spec
3. Do not add tests against old mixed-responsibility helpers just because they already exist.
4. If a rule survives decomposition, test it either:
   - on a stable new decomposed helper, or
   - at the command level
5. If a rule does not survive decomposition as an independent contract, do not test it directly.
6. Every new helper stub must have a docstring that states:
   - coordinate space
   - selection model
   - fragment ownership rule
   - counting or assignment interval assumption
   - whether aligned fetch narrowing is allowed
7. Each item in the test spec must map to at least one test before this pass is complete.
8. Nothing outside that mapping should be added in this pass.

## What is actually next

The next steps are exactly:

1. Define the new decomposed helper APIs with `unimplemented!()`.
2. Add helper-layer tests for the new APIs.
3. Add command-level tests for command contracts that remain valid after decomposition.
4. Compile with `cargo check --features cli,plotters --tests`.

No implementation of the new helpers belongs in this pass.

## Implementation pass status

The stub-and-test pass above has now been completed.

The current implementation pass is narrower than a general refactor. It is only:

1. Implement the new decomposed helpers already defined by the stub pass.
2. Wire `lengths`, `ends`, and other already-identified callers to the correct new helper/model.
3. Preserve the documented command classifications from the rules spec.
4. Change nothing outside those contracts.
5. Compile with `cargo check --features cli,plotters --tests`.

This implementation pass must still reject:

- new speculative behavior
- unrelated cleanup
- new tests that are not already required by the test spec
- helper changes that re-mix selection, aligned envelope derivation, and final clamp

## New helper layers to introduce now

The new helper stubs should make the decomposition explicit:

1. Candidate window selection by model
   - `CoreOverlap`
   - `ReachableFromTileOwnedFragments`
2. Aligned fetch-envelope derivation
   - only for models where an aligned-space derivation is valid
3. Final aligned clamp
   - already exists as `clamp_fetch_to_window_span(...)`

The new helper API should make it impossible to confuse:

- BED candidate selection
- aligned fetch-envelope derivation
- final fetch clamping

## Questions to ask before each added test

Before each test is added, the answer to all of these must be yes:

- Does this test target a contract that survives decomposition?
- Is this the strongest boundary that still tests that contract?
- If this is a command-level rule, is it tested at the command level rather than on an old helper?
- If this is a helper-level rule, is it attached to a new decomposed helper rather than a helper we
  intend to split apart?
- Is the ownership rule explicit in the fixture and the test name/comment?

If any answer is no, the test does not belong in this pass.

## Mapping discipline

For this pass, each major section of the test spec must be represented:

- Layer 1: new helper tests
- Layer 2: new helper tests
- Layer 3: existing clamp tests plus any missing additions
- Layer 4: command-level tests
- Layer 5: command-level BED vs fixed-size consistency tests

The pass is not complete until each required block in the test spec is covered by at least one
concrete test.

## Explicit execution checklist

This section must be updated before each new edit pass.

Each line below is a direct execution item:

- spec clause
- current status
- chosen boundary
- target file
- anti-drift justification

If an item is not listed here, it must not be added in the next pass.

### Checked items already covered

- Layer 1 candidate-window selection for `CoreOverlap`
  - status: covered
  - boundary: new helper-level tests
  - target: `src/shared/tiled_run_tests.rs`
  - justification: this survives decomposition as a stable selection-model contract

- Layer 1 candidate-window selection for aligned fragment reach
  - status: covered
  - boundary: new helper-level tests
  - target: `src/shared/tiled_run_tests.rs`
  - justification: this survives decomposition as a stable selection-model contract

- Layer 1 candidate-window selection for raw fragment reach
  - status: covered
  - boundary: new helper-level tests
  - target: `src/shared/tiled_run_tests.rs`
  - justification: this survives decomposition as a stable selection-model contract

- Layer 2 aligned fetch-envelope derivation for `CoreOverlap`
  - status: covered
  - boundary: new helper-level tests
  - target: `src/shared/window_fetch_tests.rs`
  - justification: this is the new decomposed helper contract, not an old mixed helper

- Layer 2 aligned fetch-envelope derivation for aligned fragment reach
  - status: covered
  - boundary: new helper-level tests
  - target: `src/shared/window_fetch_tests.rs`
  - justification: this is the new decomposed helper contract, not an old mixed helper

- Layer 2 raw `ends` BED fetch policy
  - status: covered
  - boundary: new helper-level tests plus raw command tests
  - target: `src/shared/window_fetch_tests.rs`, `tests/test_ends_command.rs`
  - justification: helper policy survives decomposition and must also hold at command level

- Layer 3 final clamp helper
  - status: covered
  - boundary: stable helper-level tests
  - target: `tests/test_tiling.rs`
  - justification: `clamp_fetch_to_window_span(...)` survives decomposition unchanged in role

- Layer 4 `lengths`
  - status: covered
  - boundary: command-level tests
  - target: `tests/test_lengths_command.rs`
  - justification: output semantics belong at the command boundary

- Layer 4 `ends`, raw mode
  - status: covered
  - boundary: command-level tests
  - target: `tests/test_ends_command.rs`
  - justification: output semantics belong at the command boundary

- Layer 4 `ends`, aligned and drop modes
  - status: covered
  - boundary: command-level tests
  - target: `tests/test_ends_command.rs`
  - justification: output semantics belong at the command boundary

- Layer 4 `midpoints`
  - status: covered
  - boundary: command-level tests
  - target: `tests/test_profile_groups_command.rs`
  - justification: the old midpoint fetch helper is not a target contract; final grouped output is
    the correct surviving boundary

- Layer 4 `fcoverage`
  - status: covered
  - boundary: command-level tests
  - target: `tests/test_fcoverage_command.rs`
  - justification: the old fcoverage fetch helper is not a target contract; final BED totals are
    the correct surviving boundary

- Layer 4 `fragment_kmers`
  - status: covered
  - boundary: command-specific helper tests
  - target: `src/commands/fragment_kmers/tiling_tests.rs`
  - justification: `determine_fetch_span(...)` is still the stable command-local contract for that
    command, not a shared mixed helper scheduled for decomposition

- Layer 4 `gc_bias`
  - status: covered
  - boundary: command-specific helper tests plus command-level cross-mode tests
  - target: `src/commands/gc_bias/windows_tests.rs`, `tests/test_gc_bias.rs`
  - justification: BED preparation semantics survive as a command-local contract and cross-mode
    agreement belongs at command level

- Layer 4 `ref_gc_bias`
  - status: covered
  - boundary: command-specific helper tests
  - target: `src/commands/ref_gc_bias/ref_gc_bias_tests.rs`
  - justification: tile-local BED clipping is a command-local contract inside `process_tile(...)`

- Layer 5 BED vs fixed-size consistency for `lengths`
  - status: covered
  - boundary: command-level tests
  - target: `tests/test_lengths_command.rs`
  - justification: cross-mode agreement belongs at the command boundary

- Layer 5 BED vs fixed-size consistency for `ends`, aligned and drop modes
  - status: covered
  - boundary: command-level tests
  - target: `tests/test_ends_command.rs`
  - justification: cross-mode agreement belongs at the command boundary

- Layer 5 BED vs fixed-size consistency for `gc_bias`
  - status: covered
  - boundary: command-level tests
  - target: `tests/test_gc_bias.rs`
  - justification: cross-mode agreement belongs at the command boundary

- Fixed-size explicit regression area for raw `ends`
  - status: covered
  - boundary: command-level tests
  - target: `tests/test_ends_command.rs`
  - justification: raw fixed-size endpoint behavior is a command contract, not a shared helper

### Remaining unchecked items

- none currently identified after the latest focused pass
  - status: pending re-audit before any next pass
  - boundary: n/a
  - target: n/a
  - justification: no new edits are allowed until a fresh spec-to-test audit identifies a real
    uncovered contract

### Current implementation pass status

- Implement decomposed candidate-window selection helpers
  - status: completed
  - boundary: shared helper implementation
  - target: `src/shared/tiled_run.rs`
  - justification: this is the stable Layer 1 contract from the rules and test specs

- Implement decomposed aligned BED support-envelope helpers
  - status: completed
  - boundary: shared helper implementation
  - target: `src/shared/window_fetch.rs`
  - justification: this is the stable Layer 2 contract from the rules and test specs

- Implement raw `ends` BED fetch policy
  - status: completed
  - boundary: shared helper implementation
  - target: `src/shared/window_fetch.rs`
  - justification: raw BED fetch policy is explicit in the rules spec and survives decomposition

- Wire `lengths` to fragment-reach aligned BED fetch logic
  - status: completed
  - boundary: command implementation
  - target: `src/commands/lengths/lengths.rs`
  - justification: `lengths` is classified as fragment-reach with aligned fetch narrowing

- Wire `ends` to aligned-vs-raw BED fetch logic
  - status: completed
  - boundary: command implementation
  - target: `src/commands/ends/ends.rs`
  - justification: `ends` must branch between aligned fragment-reach narrowing and raw full-tile
    fetch policy

- Keep `gc_bias` BED preparation on fragment-reach candidates instead of re-clamping to core overlap
  - status: completed
  - boundary: command-local implementation
  - target: `src/commands/gc_bias/windows.rs`
  - justification: `gc_bias` is classified as fragment-reach in the rules spec

- Compile with `cargo check --features cli,plotters --tests`
  - status: completed
  - boundary: repository compile check
  - target: whole workspace
  - justification: required repository verification step after code changes

## How the workdoc must be used

For every future pass in this area:

1. Re-read this checklist before editing.
2. Name the exact checklist item being worked on in the user update.
3. If a proposed test or edit does not match a checklist line, do not add it.
4. After the pass, update the status here before doing anything else.

If these steps are not followed, the pass is out of contract with this workdoc.
