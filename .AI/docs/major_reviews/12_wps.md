# `cfdna wps` review

Date: 2026-05-09

Scope: not yet reviewed. This file exists to track known WPS issues discovered while working around failing WPS tests. A full command review still needs to read the config, run path, WPS aggregation, reducers/writers, masking semantics, and existing tests before drawing broader conclusions.

## Release triage

Pre-release docs/API polish:

- WPS-001: aggregate output column `blacklisted_positions` also counts non-blacklist masked positions.

## Findings

### WPS-001 - Low - Aggregate output column `blacklisted_positions` also counts non-blacklist masked positions

WPS builds one mask for both blacklist intervals and chromosome-edge centers whose WPS window would exceed chromosome bounds ([wps.rs](../../src/commands/wps/wps.rs#L1386-L1437)). Aggregate rows then compute `blacklisted_positions` as the row span minus the allowed-position count ([wps.rs](../../src/commands/wps/wps.rs#L1461-L1498)) and write that value through the shared aggregate output schema ([writers.rs](../../src/commands/fcoverage/writers.rs#L285-L305)).

The issue is exposed by edge-touching tests with no user blacklist: `by_size_total_handles_three_chromosomes`, `by_bed_total_handles_three_chromosomes`, and `by_bed_total_skips_chromosomes_without_windows_and_keeps_later_chromosomes` expect 0 `blacklisted_positions`, but current aggregation can report chromosome-edge masks in that column ([test_wps.rs](../../tests/test_wps.rs#L561-L738)). Technically, if the column is meant literally, those rows should report 0 blacklisted positions and track edge-masked/invalid-context positions separately.

Impact: in WPS aggregate outputs, `blacklisted_positions` is really "masked positions". A run with no user blacklist can still report nonzero `blacklisted_positions` for windows touching chromosome edges. That is defensible computationally, but the column name is misleading and may confuse downstream interpretation.

Recommended fix:

- Decide whether WPS aggregate outputs should keep the shared `blacklisted_positions` schema or expose a WPS-specific name such as `masked_positions`.
- If the shared schema stays, document that WPS counts both blacklist-derived and edge-context masks in this column.
- Keep regression coverage for edge windows without user blacklists so the eventual behavior remains explicit.
