# What's missing for next release

Current state: close, but not first-release ready yet.

**Release blockers (non-feature)**
1. Public docs still contain TODO placeholders.
- [README.md:73](/Users/au547627/Documents/Development/rust/cfDNAlab/README.md:73), [README.md:340](/Users/au547627/Documents/Development/rust/cfDNAlab/README.md:340), [README.md:348](/Users/au547627/Documents/Development/rust/cfDNAlab/README.md:348), [README.md:437](/Users/au547627/Documents/Development/rust/cfDNAlab/README.md:437).

1. Validation coverage is uneven by command.
- Strong: `lengths`, `bam-to-bam`.
- Moderate: `fcoverage`, `coverage-weights`, `midpoints`, `bam-to-frag`.
- Weak: `gc-bias` and `ref-gc-bias` command-level E2E.
- Missing: `frag-to-bam` tests (no test target for this command in `tests/`).

1. Release/process hygiene is incomplete.
- No CI workflow directory (`.github/` absent).
- Cargo package metadata is minimal (no `license`, `description`, `repository`, etc. in [Cargo.toml](/Users/au547627/Documents/Development/rust/cfDNAlab/Cargo.toml)).

**Command-by-command validation plan before first release**
1. `cfdna gc-bias`
- Add full E2E test (tiny BAM + tiny 2bit + tiny ref-gc package).
- Validate output package schema and numeric sanity (finite, non-negative, expected matrix shape, compatible length/GC edges).
- Add failure tests for bad input combinations and malformed ref package.

2. `cfdna ref-gc-bias`
- Add E2E command test producing output dir artifacts and validating required files/metadata.
- Add failure tests for invalid smoothing settings and incompatible chromosome selections.

3. `cfdna coverage-weights`
- Expand current command tests to enforce contiguity and endpoint invariants of output bins.
- Add deterministic regression test with expected scaling values on a fixed fixture.

4. `cfdna fcoverage`
- Add tests for all `--per-window` modes (`average`, `total`, `unique-positions`, `indexed-positions`).
- Add explicit negative test for disallowed `--by-size` + positional-per-window modes (error path).

5. `cfdna lengths`
- Keep as near-release-grade baseline.
- Add one CLI-facing smoke test (help/options + output file naming contract).

6. `cfdna midpoints`
- Add tests for interval validation (unsorted, mixed widths, malformed group labels).
- Add tests for length-bin edge behavior and blacklist/scaling interactions.

7. `cfdna bam-to-bam`
- Add unpaired (`--reads-are-fragments`) coverage.
- Add GC-source validation tests (`--gc-file` handling, invalid/missing ref behavior).

8. `cfdna bam-to-frag`
- Expand beyond smoke: blacklist strategy behavior, by-bed filtering, optional extra columns with scaling/GC, stable sort order.

9. `cfdna frag-to-bam`
- Create full test module: parser errors, coordinate/order constraints, chromosome bounds, blacklist/mapq/length filters, empty-result behavior, BAM header correctness.

**Release gate checklist**
1. Docs gate
- Zero TODOs in README and user-facing help.
- All README commands copy-paste runnable.
- `cfdna <command> --help` manually reviewed for all 9 commands.

2. Validation gate
- `cargo check --features cli,plotters`
- Full command test suite green including new `frag-to-bam` and gc/ref-gc E2E tests.
- At least one end-to-end pipeline smoke test from README recipes.

3. Build/packaging gate
- `cargo build --release --features cli,plotters --bin cfdna`
- Add package metadata in `Cargo.toml` and a license file.
- Add CI workflow for check + tests + release build.

4. Scope gate
- Keep planned commands explicitly out of first release scope (`ends`, `fragment-kmers`, `wps-peaks` in [README.md:56](/Users/au547627/Documents/Development/rust/cfDNAlab/README.md:56)); ensure release notes only claim the 9 current commands.

If you want, I can turn this into an execution checklist file (`release_checklist.md`) and then start implementing Phase 1 (docs/help/install fixes) in one pass.