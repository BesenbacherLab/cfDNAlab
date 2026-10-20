# Output publication rollback: nice to have

Deferred by the owner on 2026-10-20. Discuss this as a package-wide policy before implementing
it. The outliers-specific implementation has been removed from the feature branch.

## Motivation

Commands already write complete temporary outputs before moving them to their final filenames.
The shared `FinalOutputFiles::move_into_place` helper moves them individually. A later rename
failure can leave a mixture of new, previous and missing outputs, with previous files replaced.

A future publication mechanism could retain previous versions until all moves succeed, remove
newly published outputs on failure, and restore previous versions. This would protect against
ordinary publication errors while the process remains running.

## Decisions before implementation

- Agree on a shared policy and its command scope, including regular files and directory outputs
  such as Zarr. Avoid command-specific methods, feature gates and backup-directory names.
- Decide handling of existing destinations, symlinks and concurrent runs using the same prefix.
- Define backup cleanup and manual recovery when restoration fails. Report both the original
  error and recovery failures, and keep recoverable originals outside disposable staging.
- State the limits clearly. Individual renames do not provide an atomic switch across filenames
  or protection against abrupt process termination.

If approved, derive deterministic tests for failure before publication, failure after earlier
moves, restoration failures, existing and new destinations, and directory outputs. Verify that
ordinary analysis failures retain the existing staging behavior.
