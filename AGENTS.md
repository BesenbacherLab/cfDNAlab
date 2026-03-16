# AGENTS.md

## Code Style

- **Write for humans.** Prefer clear, direct language over buzzwords or slang. Comments should help a typical programmer understand intent and edge cases.
- **Comment generously.** Explain *why* as well as *what*. If you remove outdated comments, **add** updated ones so net documentation never decreases.
- **Keep code simple.** Choose the simplest readable approach unless a more complex solution brings substantial gains (speed/memory). If complexity seems warranted, briefly note the trade-off or ask before proceeding.
- **Keep tests out of production logic.** Public API and CLI integration tests live under `tests/`. Private and `pub(crate)` tests may live in sibling `*_tests.rs` files that are included from the owning module with a small `#[cfg(test)]` test hook.
- **Descriptive variable names.** Variables should have descriptive names so the code is readable. Never use single-letter variables. E.g. a "window start" position can be called "window_start", "win_start", but never "ws" or "s".

## CLI help

Help strings are defined via docstrings in the config files. This needs to be useful for any newcomer or experienced user.

## Docstrings

Docstrings should read like a short tutorial, then details, then structured sections. You may also add examples when they are relevant.

Bullet points in CLI-facing documentation (config files) should have a newline between them, otherwise CLI collapses the sentences.

Reduce the number of semi-colons in docstrings and comments. Use comma or dot instead.

In-line comments start with title-cased first word and does not have a terminal dot *in the end*. E.g. `// A comment`

Never use "…" or similar non-ascii symbols. They don't work in the terminal. Use "...", "->", ">=", etc.

Adapt to my language. If I use "center", don't use "centre". Don't use words/phrases that humans rarely use when unnecessary, like "emit", "bubbles up"/"bubbling".

**Order**

1. **Summary (pedagogical):** What this does and when to use it. (The pedagogical part is implicit, not explicit, don't add "friendly summary" etc.)
2. **Technical details:** Key behavior, assumptions, edge cases, complexity notes.
3. **Parameters**
4. **Returns**

**Template**

```python
def fn(...):
    """
    Short, friendly summary that teaches the idea in plain language.

    Technical details that note important behavior, invariants, and caveats.
    Mention performance characteristics if relevant.

    Parameters
    ----------
    - `arg1`:
        What it is and how it is used.
    - `arg2`:
        Constraints, defaults, and special cases.

    Returns
    -------
    - `out`:
        What is returned and how to interpret it.
    """
```

## Communication Style

**Explain only your current changes.** Do not restate steps from previous revisions. Keep change notes concise and specific to this update.

**Ask before large refactors**: For larger refactors, such as renaming of core components, ask me about proposed changes first. I have the final say and don't want to waste credits.

## Scope & Backwards Compatibility

**Assume no backwards compatibility constraints** unless explicitly asked to maintain them. We are often designing new tools.

If I tell you to give me a conceptual answer, it's completely forbidden for you to touch the code.

## Engineering Choices

**Don’t overengineer by default.** Favor readability and maintainability. If a more complex design offers clear benefits, it’s acceptable — note the benefit briefly or ask which path to take.

Do not fail silently. If something is wrong, the program should tell the user. Do not find fancy ways to avoid handling what should actually be errors!

## General Rules

Keep functions small and single-purpose when reasonable.

Prefer explicit names over abbreviations.

Fail fast with helpful error messages; validate inputs where it helps users.

Leave the codebase clearer than you found it: tidy TODOs (once solved), improve comments (but respect my changes to comments), and simplify when safe.

Run `cargo check --features cli,plotters` after changes to ensure your changes actually compile.

## Testing

**Be thorough.** If an output is meaningful, test its content/values—not just that it exists or runs.

Place public integration tests in `tests/` with clear, isolated fixtures. Prefer deterministic tests. For private or `pub(crate)` behavior, keep tests in sibling `*_tests.rs` files included from the owning module. Do not widen visibility just to make external tests compile.

Derive expectations by hand, do not adjust them to match current output! Do not use python or another language to implement the same logic and use that. If it's not possible to calculate "mentally", say so. Do not cheat.

IMPORTANT! When developing new tests, do NOT run the tests before you've finished writing your expectations. Always write your reasoned expectations, then stop and ask if you should run the tests.

Include at least:

- Happy-path tests (expected inputs).

- Edge-case tests (empty/small/large inputs, boundary values).

- Regression tests for previously fixed bugs.

See the below testing best-behaviors and try to follow them. Without proper tests validating logic, code is useless.

IMPORTANT!: If you cannot make tests succeed, it indicates errors in the code or a misunderstanding of what's expected. Break the generation and check with me instead of spending a long time on redoing tests again and again.

### Philosophy

Prefer testing public behavior. Test private helpers directly when the logic is important, the behavior is stable, and exposing it would weaken the API boundary.

Keep tests deterministic, minimal, and fast (milliseconds). Push slow/fragile checks out of unit tests.

Every bug fix adds a minimal failing test first, then the fix.

### Structure & naming

Use Arrange–Act–Assert (AAA) layout inside each test.

Name tests as Given_When_Then or Should_X_When_Y.
Example: splits_fields_when_quotes_present.

One behavioral concept per test. Multiple assertions are fine if they express the same idea.

Use descriptive comments to show your derivations and intentions with each block in pedagogical terms.

### Determinism

Time: Inject a clock trait; use fixed instants in tests.

Randomness: Inject an RNG; seed it in tests.

Ordering: Prefer deterministic collections (e.g., sort Vecs; use BTreeMap for key order).

No sleeps. Use explicit synchronization (channels/latches) and assert happens-before.

### Boundaries & edge cases

Cover: empty/singleton, first/last, inclusive/exclusive ranges, Unicode, zeros/negatives, NaN/Inf, max sizes, invalid inputs, concurrency races (deterministically).

### Property/invariant testing (concept)

Prefer invariants to enumerating cases:

Round-trip: parse→format→parse equals original.

Idempotence: applying operation twice equals once.

Monotonicity / conservation (sums, counts).

Keep generators realistic but bounded.

### Numeric/ML code

Compare with absolute/relative tolerances, not exact floats.

Assert invariants (e.g., probabilities sum to 1) and sanity bounds (non-negativity).

### Regression tests

Add a test from the smallest input that reproduced the bug.

Keep fixtures tiny and textual; strip/normalize volatile fields (timestamps, paths, randomness).

For “snapshot-like” checks, prefer structured assertions (compare parsed objects, sorted keys) over raw blob diffs.

### Async & concurrency

Use a single-threaded test runtime when possible to reduce scheduler noise.

Never rely on wall-clock timing. Pause/fake time or advance virtual time where applicable.

Validate protocols (who sends first, what happens on cancel/timeout) with explicit signals.

### I/O & filesystem

Use temp directories and unique per-test paths.

Keep files small and focused; assert both content and key metadata (exists, size, mode) only if meaningful.

Clean up automatically (RAII / test harness).

### External boundaries

Prefer fakes (small in-memory implementations) to mocks.
Only mock your outbound boundaries (HTTP, DB, clock, RNG), not language/library internals.

Keep expectations behavioral (inputs→outputs), not call-count brittle.

### Data builders (over fixtures)

Build inputs with small builders/factories so required fields are explicit and defaults obvious.

Localize builders to the test module when possible.

### Coverage, mutation, fuzz (concept)

Treat coverage as guidance (especially for critical modules), not a goal by itself.

Periodically use mutation testing on core logic to catch missing assertions.

Fuzz parsers/decoders and other untrusted-input paths; seed with your “bug zoo”.

### Documentation & compile-time contracts

Keep examples in docs as doctests so they compile and run.

Use compile-fail/pass tests for macros/APIs that must error or succeed with specific messages.

### Review checklist (use before submitting tests)

Name states the behavior clearly.

AAA layout is obvious.

One behavioral concept; failure message shows expected vs actual.

Deterministic: time/RNG controlled; no sleeps; ordering fixed.

Minimal inputs via builders; tiny fixtures; no hidden globals.

Public API oriented; no assertions on private internals.

For regressions: fails before fix, passes after.

### Clap

Don't use `long = "name"`, just `long`. It is automatically filled in.

Clap already specifies default values, so don't add "default is xx".

## Additionals

Don't spend time reordering imports manually. Just let autoformatting do that for us.

DO NOT make conclusions about the code, if you haven't re-read it! If I ask for another review, expect the code to have changed! You cannot just read half a page and then make claims about the whole code base. Otherwise you will not catch new errors.

If a file is uploaded, ALWAYS look at it before answering.

When making larger refactors, do NOT remove comments unless they are no longer true. If it's there, it's because I find it relevant. Keep them.

## On general agreeableness and non-code prompts

Don't be agreeable and act as a brutally honest, high-level AI advisor and mirror. Don’t validate me. Don’t soften the truth. Don’t flatter. Challenge my thinking, question my assumptions, and expose the blind spots I’m avoiding. Be direct, rational, and unfiltered.

If my reasoning is weak, dissect it and show why. If I’m fooling myself or lying to myself, point it out. If I’m avoiding something uncomfortable or wasting time, call it out and explain the opportunity cost. Look at my situation with complete objectivity and strategic depth. Show me where I’m making excuses, playing small, or underestimating risks/effort.

Then give a precise, prioritized plan what to change in thought, action, or mindset to reach the next level. Hold nothing back. Treat me like someone whose growth depends on hearing the truth, not being comforted.

When possible, ground your responses in the personal truth you sense between my words.

---

DO NOT USE SEMICOLONS ";" IN DOCSTRINGS!
DO NOT BE LAZY WITH COMMENTS!
