# Major command review plan

Date: 2026-04-24

Purpose: build one review document per released command, plus one shared document for package-level findings. The goal is to let shared fixes happen first, then command-specific fixes happen one command at a time without duplicated findings.

## Release boundary

The initial review scope is the command list in the README and the default Cargo feature set:

1. `cfdna fcoverage`
2. `cfdna midpoints`
3. `cfdna ends`
4. `cfdna lengths`
5. `cfdna ref-gc-bias`
6. `cfdna gc-bias`
7. `cfdna fragment-count-weights`
8. `cfdna coverage-weights`
9. `cfdna bam-to-bam`
10. `cfdna bam-to-frag`
11. `cfdna frag-to-bam`

Unreleased commands should be ignored except when released commands depend on their shared code.

## File layout

- `00_review_plan.md`: this workflow and command order.
- `00_shared_package_notes.md`: cross-command findings only.
- `01_fcoverage.md`, `02_midpoints.md`, etc.: one command-specific review per prompt.

## Review rules

- Do not repeat a finding in both shared and command-specific files. If a command is affected by a shared issue, the command file can point to the shared finding instead of restating it.
- Prefer evidence from current code over comments or old plans.
- Use source links with line references for every concrete finding.
- Separate bugs, scientific/semantic risks, performance risks, documentation inconsistencies, and test gaps.
- Keep active findings triaged into pre-release correctness/safety, pre-release docs/API polish, and post-release performance optimization. Do not leave completed findings in the active queue.
- Do not flag multiple outputs or auxiliary files as a problem by themselves. Only track output-related issues when they overwrite another output, obscure the documented output contract, make a requested primary artifact unavailable, or otherwise create a concrete correctness/safety risk.
- Treat backwards compatibility as out of scope unless explicitly requested.
- Do not run tests during review. Test suggestions should be derived by reading code and existing test coverage.

## Per-command checklist

For each command:

1. Read the config, CLI dispatch, run entrypoint, reducers/writers, and directly used shared helpers.
2. Read existing command tests enough to know what behavior is already covered.
3. Check input validation and fail-fast behavior.
4. Check output schema, filenames, sorting, compression, and sidecars.
5. Check fragment/window/blacklist/GC/scaling semantics where relevant.
6. Check performance bottlenecks that follow from the actual control flow.
7. Write only new findings into that command's markdown file.
