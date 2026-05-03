# `cfdna bam-to-bam` review

Date: 2026-05-04

Scope: `src/commands/bam_to_bam/*`, CLI dispatch, shared BAM pairing and read-filter helpers, interval/overlap helpers used by the command, scaling-factor loading/application, blacklist filtering, and the existing `bam-to-bam` tests in `tests/test_bam_to_bam_command.rs` and cross-command converter tests. I did not run tests.

Shared findings that affect this command:

- G-017 in `00_shared_package_notes.md`: converter AUX tag names are longer than the BAM tag field.

## Release triage

Pre-release correctness/safety:

- B2B-001: default lexicographic chromosome processing can produce output that is not coordinate-sorted by BAM header target id.
- G-017: documented converter AUX tag names are not the actual serialized BAM tag keys.

Pre-release docs/API polish:

- None separate from the correctness findings above.

Post-release performance:

- None currently active.

## Findings

### B2B-001 - High - Default lexicographic chromosome order can produce BAMs that are not coordinate-sorted by header target id

The command help promises that `bam-to-bam` writes a coordinate-sorted BAM ([config.rs](../../src/commands/bam_to_bam/config.rs#L8-L13), [config.rs](../../src/commands/bam_to_bam/config.rs#L70-L81)). The same config says the command sorts selected chromosome names lexicographically by default, with `--skip-chromosome-sort` preserving the specified order instead ([config.rs](../../src/commands/bam_to_bam/config.rs#L96-L102)). The implementation does that sort directly after resolving chromosomes from the input BAM ([bam_to_bam.rs](../../src/commands/bam_to_bam/bam_to_bam.rs#L89-L93)).

The output header is copied from the input BAM ([bam_to_bam.rs](../../src/commands/bam_to_bam/bam_to_bam.rs#L168-L173)), and records are then written chromosome-by-chromosome in the potentially lexicographic or user-provided order ([bam_to_bam.rs](../../src/commands/bam_to_bam/bam_to_bam.rs#L175-L199)). That is not enough to guarantee coordinate-sorted BAM output, because BAM coordinate order is ordered by the header target ids, then position. If the input header order is `chr2, chr10, chr1`, lexicographic output order `chr1, chr10, chr2` produces a target-id sequence `2, 1, 0`, which is decreasing even though each chromosome is internally sorted.

Existing tests explicitly lock in this behavior for non-header chromosome order. One test constructs a header ordered `chr2, chr10, chr1` and expects default output records in `chr1, chr10, chr2` order ([test_bam_to_bam_command.rs](../../tests/test_bam_to_bam_command.rs#L315-L380)). Another checks that `--skip-chromosome-sort` follows the resolved header order only for `--chromosomes all` ([test_bam_to_bam_command.rs](../../tests/test_bam_to_bam_command.rs#L385-L430)), while user-provided chromosome lists can still request arbitrary order ([test_bam_to_bam_command.rs](../../tests/test_bam_to_bam_command.rs#L286-L309)).

Impact: output can be advertised as coordinate-sorted while violating the order expected by BAM indexers and fetch-based readers. This is release-relevant because a user may immediately index the converted BAM or feed it to another tool that assumes nondecreasing `(tid, pos)` order.

Recommended fix:

- For coordinate-sorted output, process selected chromosomes in input header target order, regardless of the textual chromosome-name order.
- If a custom chromosome order is still needed, either rebuild the output header and remap `tid`/`mtid` consistently, or make the mode explicit that the output is not coordinate-sorted and should not be indexed as coordinate-sorted.
- Add coverage that uses a non-lexicographic header and checks nondecreasing target ids in the output, ideally also checking that indexing/fetch expectations are not violated.

## Existing coverage notes

The command has direct coverage for MAPQ filtering, blacklist filtering, BED inclusion filtering, default MAPQ behavior, chromosome-order behavior, multi-chromosome output, scaling tags, paired/unpaired scaling equivalence, GC correction, neutralizing invalid GC, combined filters, scaling-file metadata mismatches, and scaling TSV chromosome coverage. The important gap from this review is that tag tests query overlong tag names and chromosome-order tests assert record name order without verifying BAM header target-id sortedness.
