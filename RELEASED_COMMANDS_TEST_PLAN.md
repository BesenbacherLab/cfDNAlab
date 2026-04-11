# Released Commands Test Plan

This plan is for the intended public command set from `CHANGELOG`, not the experimental feature-gated commands.

Scope:

- `cfdna fcoverage`
- `cfdna lengths`
- `cfdna midpoints`
- `cfdna gc-bias`
- `cfdna coverage-weights`
- `cfdna bam-to-bam`
- `cfdna bam-to-frag`
- `cfdna frag-to-bam`
- `cfdna ref-gc-bias`

This document was cleaned against the live test tree on 2026-03-25.

Only keep gaps that still look open after rereading the current tests.
If a gap is covered, delete it from this file instead of leaving stale fear behind.

## Goal

The goal is not "more tests".
The goal is release confidence that the commands preserve the scientific semantics they claim to preserve:

- Fragment semantics
- Window semantics
- Tiling invariance
- Weight composition
- Artifact contracts
- Cross-command interoperability
- Default behavior

## Current Coverage Snapshot

The repo is not thinly tested anymore.
Several claims in the old version of this plan are no longer true.

### Commands that already look strong

- `lengths`
  Command coverage is broad. It already covers global/by-size/by-bed paths, multi-chromosome cases, scaling, GC correction, GC length marginalization, indel modes, blacklist behavior, reducer logic, and key tiling helpers.

- `fcoverage`
  Command coverage is broad. It already covers positional and reduced outputs, tile-size invariance, aligned fast-path vs general reducer parity, blacklist behavior, `--ignore-gap`, GC file vs GC tag behavior, real `coverage-weights` TSV consumers, and real `ref-gc-bias -> gc-bias` workflows.

- `frag-to-bam`
  Coverage is strong for a format-boundary command. There are many workflow and roundtrip tests, including restoration of optional `GC` and `COV` tags through downstream consumers.

### Commands that now have substantial command-level coverage

- `gc-bias`
  The command is no longer helper-only. Current tests exercise default window behavior, explicit `--by-size`, explicit `--by-bed`, aligned vs misaligned tiles, empty middle tiles, `save_intermediates`, `min_window_acgt_pct`, real `ref-gc-bias -> gc-bias` workflows, outlier methods and scopes, loader failures, bin-merging knobs, and end-offset behavior.

- `ref-gc-bias`
  The command is no longer missing `run()` coverage. Current tests exercise package metadata, end-to-end count shapes, blacklist masking, BED flattening, end offsets, support-threshold plateaus, thread-count invariance with fixed seeds, and several producer-consumer contracts.

- `midpoints`
  Coverage is now substantial. Current tests exercise group and length-bin behavior, GC-from-file and GC-from-tag workflows, scaling TSV consumers, real `ref-gc-bias -> gc-bias` workflows, and the ndarray view/copy parity that used to be called out as a missing source-level test.

- `bam-to-frag`
  Coverage is now substantial. Current tests exercise chromosome ordering, default MAPQ, GC fallback behavior, schema/range failures, real `coverage-weights` TSV consumers, and real `ref-gc-bias -> gc-bias` parity against `bam-to-bam`.

### Commands that still look moderate rather than complete

- `bam-to-bam`
  The suite now covers filtering, blacklist behavior, window inclusion, chromosome ordering, COV tagging, GC fallback tagging, schema/range failures, and real GC/scaling producer workflows. It is no longer thin, but it still has a few obvious missing branches.

- `coverage-weights`
  The suite now covers hand-derived scaling values, chromosome coverage without gaps, valid/invalid stride checks, shared global mean, blacklist-aware normalization, output row ordering, and default MAPQ behavior. It still lacks a few edge cases with the highest numerical risk.

## Old Plan Items That Are No Longer Open

The old document had several claims that are now stale and should stay removed:

- It is no longer true that `ref_gc_bias::run()` lacks meaningful command-level coverage
- It is no longer true that `gc-bias` has little or no end-to-end `run()` coverage
- The midpoint ndarray stride/view parity test is now present
- The `ref-gc-bias` support-threshold step-function is now tested
- The `fcoverage` "BED window wider than 2x tile size" case is now tested
- The `bam-to-frag` GC fallback semantics are now tested

The plan should focus on what is still missing, not on historical gaps that have already been closed.

## Current Status After The 2026-03-25 Backfill

The concrete missing items from the cleaned 2026-03-25 audit were added:

- `ref-gc-bias`
  Added command-level smoothing and interpolation tests with hand-derived expectations.

- `ref-gc-bias`
  Added fixed-seed tile-size invariance.

- `gc-bias`
  Added real producer-consumer tests proving that reference-package smoothing and interpolation metadata are respected by the command.

- `bam-to-bam`
  Added combined `GC` + `COV` + `FLEN` tag tests, a combined BED/blacklist/scaling/GC test, and `skip_invalid_gc=true` coverage.

- `midpoints`
  Added even-length midpoint tie tests at window edges and an explicit blacklist-midpoint-vs-placement contract test.

- `lengths`
  Added an even-length midpoint assignment seam test.

- `fcoverage` and `lengths`
  Added a direct cross-command parity test on the same MAPQ-filtered fixture.

- `ref-gc-bias -> gc-bias -> coverage-weights -> fcoverage`
  Added one compact real artifact-chain workflow test as a release-confidence spine.

- Shared real-artifact consumer fixtures
  Added one shared real-artifact builder for a neutral `ref-gc-bias -> gc-bias` package plus a
  real `coverage-weights` TSV, then used it in isolated per-consumer contract tests for
  `bam-to-bam`, `lengths`, `midpoints`, and `fcoverage`.

## Remaining Work

No concrete high-priority command-specific test gaps are currently known from the latest reread.

If more test work is needed, it should now focus on optional structural consolidation rather than
plugging obvious missing branches:

- Continue deleting stale claims from this file whenever the tree changes

## Execution Order

1. Keep the release-spine workflow passing as the code evolves
2. Keep shared real-artifact fixtures readable and scoped to one consumer per test
3. Re-audit the plan only after rereading the current tests and source again

## What To Avoid

- Do not keep stale "under-tested" claims once tests exist
- Do not add gaps to this file without first rereading the tests for that command
- Do not inflate the plan with small CLI polish items while real scientific branches are still unpinned
- Do not treat historical audit notes as still open unless they are rechecked against the current tree

## First Concrete Deliverables

- The current command suites plus the compact release-spine workflow
