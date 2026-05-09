# Finalized Specs

This directory stores the current command and shared-behavior contracts that future development should preserve.

Rules:

- Filenames are lower-snake-case and do not include dates.
- Specs should describe current decided behavior, not planning history.
- Keep future ideas, temporary plans, and dated review notes outside this directory, usually under `.AI/docs/future/`.
- When a plan has been implemented, move only the lasting invariant into the relevant spec.
- When intended behavior and implementation disagree, add a short `! Warning:` note under the affected section or move the item to future review.
- Do not copy full CLI help here. Capture the invariants, output schemas, geometry choices, and edge cases that future developers are likely to break.
