# Testing

**Be thorough.** If an output is meaningful, test its content/values, not just that it exists or runs.

All important logic must be tested. Use the strongest boundary that still tests the right behavior.

Place tests in `tests/` with clear, isolated fixtures when testing public crate behavior, CLI behavior, end-to-end workflows, and released regressions.

For internal logic that should stay private or `pub(crate)`, keep tests in sibling `*_tests.rs` files included from the owning module. Do not widen visibility just to make an external test compile.

Derive expectations by hand, do not adjust them to match current output! Do not use python or another language to implement the same logic and use that. If it's not possible to calculate "mentally", say so. Do not cheat.

IMPORTANT! When developing new tests, do NOT run the tests before you've finished writing your expectations. Always write your reasoned expectations, then stop and ask if you should run the tests.

Include at least:

- Happy-path tests (expected inputs).

- Edge-case tests (empty/small/large inputs, boundary values).

- Regression tests for previously fixed bugs.

See the below testing best-behaviors and try to follow them. Without proper tests validating logic, code is useless.

IMPORTANT!: If you cannot make tests succeed, it indicates errors in the code or a misunderstanding of what's expected. Break the generation and check with me instead of spending a long time on redoing tests again and again.

Whenever you add a new test, you need to finish with an extra thinking pass about whether it can become even stronger. Don't settle on "good enough" for tests.

### Philosophy

Test all important behavior and logic.

Prefer public-behavior tests when they cover the logic well.

Test private helpers directly when that gives clearer or smaller coverage, or when exposing them would weaken the API boundary.

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

Boundary choice is intentional: public behavior in `tests/`, private logic in module-local test files when needed.

For regressions: fails before fix, passes after.
