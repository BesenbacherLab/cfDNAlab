# `cfdna bam-to-frag` review

Date: 2026-05-04

Scope: `src/commands/bam_to_frag/*`, shared BAM fragment construction used by frag-file output, chromosome resolution, scaling-factor loading/application, GC correction, blacklist and BED-window filtering, and the existing `bam-to-frag` tests in `tests/test_bam_to_frag.rs` plus converter round-trip tests in `tests/test_frag_to_bam_command.rs`. I did not run tests.

Shared findings that affect this command:

- None active. The converter temporary path issue originally noted here as G-018 has since been implemented through shared temporary chromosome tokens.

## Release triage

Pre-release correctness/safety:

- None active.

Pre-release docs/API polish:

- None currently active.

Post-release performance:

- None currently active.

## Findings

No active command-specific correctness finding in this pass.

The main `bam-to-frag` data path is internally consistent with the command contract I re-read: it builds paired fragments as `forward.pos` to `reverse.reference_end`, writes `start`, `end`, minimum MAPQ, and read1 strand to the frag rows, and adds optional GC/scaling columns only when those inputs are configured ([config.rs](../../src/commands/bam_to_frag/config.rs#L8-L34), [bam_to_frag.rs](../../src/commands/bam_to_frag/bam_to_frag.rs#L415-L520)). The companion header is generated from the same option set as the row writer ([bam_to_frag.rs](../../src/commands/bam_to_frag/bam_to_frag.rs#L226-L251)).

## Existing coverage notes

The command has direct coverage for global and BED-window output, chromosome ordering, default MAPQ behavior, proper-pair/read filtering behavior through fixtures, GC correction including neutralization vs skipping invalid GC, coverage and count scaling columns, combined GC/scaling columns, scaling TSV chromosome coverage, scaling metadata compatibility, and real converter round-trips through `frag-to-bam`.

The previous path-safety gap for unusual chromosome names has been addressed through shared temporary chromosome tokens.
