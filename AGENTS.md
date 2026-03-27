# AGENTS.md

This file is the authoritative entry point for repo-specific agent instructions.

## Always Apply

- Write for humans and prefer clear, direct language.
- Do not remove comments unless they are no longer true. If you remove one, replace it with an updated comment so documentation quality does not go down.
- Keep functions small and single-purpose when reasonable.
- Prefer explicit names over abbreviations. Do not use single-letter variable names.
- Do not fail silently. If something is wrong, the program should tell the user.
- Run `cargo check --features cli,plotters` after code changes. After major refactors, run `cargo check --all-features`.
- If a file is changed, always read it before answering.
- If I ask for a new code review, never rely on memory to answer.
- Do not make conclusions about code you have not re-read.
- Base answers about existing functionality on actual code behavior, not comments etc. that might be outdated.
- For fragment code, preserve the project's domain semantics and vocabulary. Paired fragment spans are defined directionally as `forward.pos` to `reverse.reference_end`, and docs/comments should keep using `pos` / `end` / `reference_end` terminology even if the implementation stores checked intervals internally.
- Read the Interval and IndexedInterval API and default to using the helpers when working on interval-logic.

## Read These Files When Relevant

- For code style, comments, docstrings, CLI help, clap, and formatting-related rules, read [.AI/writing_style.md](/Users/au547627/Documents/Development/rust/cfDNAlab/.AI/writing_style.md)
- For testing rules and test philosophy, read [.AI/testing.md](/Users/au547627/Documents/Development/rust/cfDNAlab/.AI/testing.md)
- For test placement, visibility, and public API boundary rules, read [.AI/api_boundaries.md](/Users/au547627/Documents/Development/rust/cfDNAlab/.AI/api_boundaries.md)
- For reducers, aggregation, masking, tiling, normalization, and related scientific/counting code, read [.AI/scientific_code.md](/Users/au547627/Documents/Development/rust/cfDNAlab/.AI/scientific_code.md)
- For communication style, scope, backwards-compatibility assumptions, engineering choices, and non-code interaction style, read [.AI/collaboration.md](/Users/au547627/Documents/Development/rust/cfDNAlab/.AI/collaboration.md)
