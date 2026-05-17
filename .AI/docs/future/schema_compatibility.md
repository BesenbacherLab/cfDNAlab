# Schema Compatibility

This table records which cfDNAlab output schemas are written by which CLI
versions and which helper package versions can read them.

Use this for release planning and support questions when users have already
generated output files with an older cfDNAlab version.

## Policy

The schema name and schema version are the output contract:

```text
cfdnalab_schema
cfdnalab_schema_version
```

Rust, Python, and R package versions do not need to match. Helper packages must
validate schema names and versions and give useful errors when a file is too new
or too old.

When a schema changes, decide explicitly whether old output should remain
readable or whether users should rerun the command. Do not rely on implicit
package-version matching.

## Compatibility Table

| Schema | Schema version | Written by CLI versions | Python helper | R helper | Recommendation |
| --- | ---: | --- | --- | --- | --- |
| `midpoint_profiles` | 1 | `cfdnalab` 0.2.x | Python `cfdnalab` >= 0.1.0 | R `cfdnalab` >= 0.1.0 | Read normally. |
| `end_motif_counts` | 1 | `cfdnalab` 0.2.x | Python `cfdnalab` >= 0.1.0 | R `cfdnalab` >= 0.1.0 | Read normally. |
| `reference_gc_package` | 3 | `cfdnalab` 0.2.x | not applicable | not applicable | Use through `cfdna gc-bias --ref-gc-file`. |
| `gc_correction_package` | 3 | `cfdnalab` 0.2.x | not applicable | not applicable | Use through feature extraction `--gc-file`. |

## Recommendation Terms

- `Read normally`: current helper packages should read this schema.
- `Install older helper`: newer helpers intentionally dropped support; tell
  users the last compatible helper version.
- `Install newer helper`: output is valid, but the user's helper package is too
  old.
- `Rerun`: the schema was experimental, incomplete, or not worth supporting
  long-term.
- `Pending`: output exists or is planned, but helper package support is not
  ready.

## Release Update Checklist

When public output schemas change:

- update the writer schema version if the downstream interpretation changes
- update this table
- update Python and R helper validation/support notes
- update downstream compatibility tests
- update user-facing docs if the schema has been released
- decide whether existing files should be read with old helpers, new helpers, or
  rerun
