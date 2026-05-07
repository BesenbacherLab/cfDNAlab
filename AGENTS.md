# AGENTS.md

This file is the authoritative entry point for repo-specific agent instructions.

## Always Apply

- Write for humans and prefer clear, direct language.
- Do not remove comments unless they are no longer true. If you remove one, replace it with an updated comment so documentation quality does not go down.
- Keep functions small and single-purpose when reasonable.
- Prefer explicit names over abbreviations. Do not use single-letter variable names.
- Do not fail silently. If something is wrong, the program should tell the user.
- Run `cargo check --features cli,plotters` after code changes and `cargo check --tests --features cli,plotters` after test code changes. After major refactors, run `cargo check --all-features`. When working on non-default commands (see cargo.toml), include their command features in these calls as well.
- Do NOT run tests. I run tests in this project. You use mental derivation only.
- If a file is changed, always read it before answering.
- If I ask for a new code review, never rely on memory to answer.
- Do not make conclusions about code you have not re-read.
- Base answers about existing functionality on actual code behavior, not comments etc. that might be outdated.
- For fragment code, preserve the project's domain semantics and vocabulary. Paired fragment spans are defined directionally as `forward.pos` to `reverse.reference_end`, and docs/comments should keep using `pos` / `end` / `reference_end` terminology even if the implementation stores checked intervals internally.
- Read the Interval and IndexedInterval API and default to using the helpers when working on interval-logic.
- Keep `.AI/docs/specs/` for concise current specs only. Do not put dated plan filenames there. Temporary plans, future ideas, review notes, and dated specs belong under `.AI/docs/future/` or another non-finalized docs folder. When implemented behavior becomes the current decision, distill only the lasting invariants into the relevant file under `.AI/docs/specs/`.
- The minimum allowed fragment length possible is 10bp. Do not use smaller values than that in test fixtures. Note that commands often set a minimum fragment length inclusion filter of 30bp, so check up on that. In general, check argument constraints before setting them in fixtures.
- Before doing anything due to "backwards compatibility", ask whether this is a concern first. In some cases it is, in some it's not and should not lead to clutter.

Always read [.AI/collaboration.md](/Users/au547627/Documents/Development/rust/cfDNAlab/.AI/collaboration.md) to avoid annoying sycophanting. I want truth not praise.

## Read These Files When Relevant

- For code style, comments, docstrings, CLI help, clap, and formatting-related rules, read [.AI/writing_style.md](/Users/au547627/Documents/Development/rust/cfDNAlab/.AI/writing_style.md)
- For testing rules and test philosophy, read [.AI/testing.md](/Users/au547627/Documents/Development/rust/cfDNAlab/.AI/testing.md)
- For test placement, visibility, and public API boundary rules, read [.AI/api_boundaries.md](/Users/au547627/Documents/Development/rust/cfDNAlab/.AI/api_boundaries.md)
- For reducers, aggregation, masking, tiling, normalization, and related scientific/counting code, read [.AI/scientific_code.md](/Users/au547627/Documents/Development/rust/cfDNAlab/.AI/scientific_code.md)
- For communication style, scope, backwards-compatibility assumptions, engineering choices, and non-code interaction style, read [.AI/collaboration.md](/Users/au547627/Documents/Development/rust/cfDNAlab/.AI/collaboration.md)
