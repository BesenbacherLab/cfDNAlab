# `cfdna frag-to-bam` review

Date: 2026-05-04

Scope: `src/commands/frag_to_bam/*`, chromosome-size loading, frag column/header detection, fragment parsing and filtering, temp staging, BAM header/record writing, optional AUX-tag transfer, and the existing `frag-to-bam` tests in `tests/test_frag_to_bam_command.rs`. I did not run tests.

Shared findings that affect this command:

- None active. The converter AUX-tag issue (G-017) and temporary path issue (G-018) originally noted here have since been implemented.

## Release triage

Pre-release correctness/safety:

- None active from shared findings.

Pre-release docs/API polish:

- F2B-001: headerless frag files with trailing columns ignore those columns even without `--ignore-extras`.

Post-release performance:

- None currently active.

## Findings

### F2B-001 - Low - Headerless frag files with trailing columns silently ignore those columns even without `--ignore-extras`

The CLI help says `frag-to-bam` accepts headerless 5-column frag files when no inline, explicit, or companion header is found, and says to use `--ignore-extras` when intentionally ignoring columns after the first five ([config.rs](../../src/commands/frag_to_bam/config.rs#L77-L91)). With no header source, the implementation resolves a fixed five-column layout and sets all recognized extra-column indices to `None` ([frag_to_bam.rs](../../src/commands/frag_to_bam/frag_to_bam.rs#L637-L691), [frag_to_bam.rs](../../src/commands/frag_to_bam/frag_to_bam.rs#L800-L811)). `parse_frag_line()` then parses required fields by index but does not reject additional trailing columns ([frag_to_bam.rs](../../src/commands/frag_to_bam/frag_to_bam.rs#L405-L465)).

There is an existing test that locks in the behavior: an `input.tsv` row with eight columns is accepted and the trailing values are ignored with no GC/cw/fl tags written ([test_frag_to_bam_command.rs](../../tests/test_frag_to_bam_command.rs#L821-L864)). That may be intentional, but it weakens the fail-fast story. A missing companion header or a mistyped filename can silently drop metadata that looks like `gc_weight`, scaling weights, or fragment length values.

Impact: not a BAM corruption bug, but a user can lose extra metadata without realizing it. This is especially easy when moving `bam-to-frag` output away from its companion `*.frag.header.tsv` file or renaming the frag file so companion-header auto-detection no longer applies.

Recommended fix:

- If no header is available and a data row has more than five columns, fail unless `--ignore-extras` is set.
- Keep accepting exactly five-column headerless files without extra flags.
- Add a regression for "headerless more-than-five columns without `--ignore-extras` fails", and keep a separate `--ignore-extras` success case for deliberate metadata discard.

## Existing coverage notes

The command has direct coverage for basic conversion, chromosome filtering, length and MAPQ filters, blacklist filtering, chromosome ordering checks, chromosome-size bounds, inline and explicit headers, companion header detection, supported and unsupported extra columns, `--ignore-extras`, `--allow-unknown-extras`, optional missing extra values, and cross-command `bam-to-frag` round-trips.

The serialized AUX-key and temporary path-safety gaps from this pass have since been implemented. The remaining important gap is fail-fast coverage for headerless files that contain trailing metadata columns.
