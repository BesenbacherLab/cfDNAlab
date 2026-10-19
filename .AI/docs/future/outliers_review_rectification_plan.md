# Outliers review checklist and remaining work

The first BAM-inspection version is implemented. Current behavior and its rationale are in the
[spec](../specs/outliers.md). This checklist tracks remaining work, not alternative specifications.

## Review inputs

- [Current spec](../specs/outliers.md): coverage, fitting, preservation, weighting, and outputs.
- [CLI help](../../../src/commands/outliers/config.rs): self-contained option reference.
- [Output README](../../../src/commands/outliers/readme.rs): files and column explanations.
- [Validation list](outliers_real_data_validation.md): later investigations and counterexamples.
- [Fitted ZIP/ZINB comparison](outliers_fragment_simulations/zinb_comparison_report.md): simulation
  evidence under its recorded settings, including the historical 10% preservation cutoff.

## Remaining implementation checks

- **ZIP error guidance (O-007):** `fit_two_stage_zip` still forwards the truncated-fit error
  without the chosen tail probability or advice to reduce it. Bring its insufficient-support
  errors in line with the ZINB path. Keep explicit failure rather than silently changing T1.
- **Command-level settings checks (O-004/O-014):** loader tests cover reference lengths and the
  distinction between analysed empty chromosomes and chromosomes not analysed. Confirm the
  corresponding `fcoverage` command coverage and detection-filter warning coverage in the next
  review. Do not claim these checks are absent solely from test names.

## What remains to decide

1. Inspect real BAM calls and matched raw/weighted coverage. Measure genomic extent separately
   from coverage reduction and the effect on relevant profiles or features. Decide whether the
   benefit warrants prioritizing further work before the application-note deadline.
2. Assess biological retention and technical sensitivity across the required 0.1x to 200x range,
   including narrow gains, moderate artifacts, and spikes inside broader gains. The range is a
   validation requirement, not a demonstrated guarantee.
3. Assess the initial 5% preservation rule and boundary behavior. Smarter evidence using regional
   coverage is a separate algorithm decision, justified by observed failures before implementation.
4. Assess dispersion, low-depth parameter uncertainty, convergence, and extreme-tail calibration.
   A good first BAM example does not establish behavior across libraries and depths. Historical
   simulation outcomes do not measure the current preservation rule.
5. Finalize CLI names and review documentation. Decide merge readiness after scientific and
   implementation checks. Inspection readiness alone does not establish merge readiness.

## Deferred work

- GC/scaling producer integration (O-002) and consumers beyond `fcoverage`, including midpoints.
- Moving the overlap resolver into shared code (O-011), when another consumer needs it.
- Smarter regional preservation evidence, alternative overlap weighting, and broader performance
  work, guided by the validation results.
- Weighting restore-mean statistics. This is a user-owned consideration. Do not investigate or
  implement it unless the user explicitly asks.

These items are not part of the first inspection milestone. No normalization changes or broad
refactors are authorized by this checklist.

## Verification

Project tests are authored and compiled. The user runs them. Do not infer correctness from
compilation alone. Derive new expectations independently, without duplicating existing shared
coverage tests for gap or read-as-fragment handling.

Use compile checks, including the command feature, after code changes:

```bash
cargo check --offline --locked
cargo check --offline --locked --tests --features "testing cmd_outliers"
cargo check --offline --locked --no-default-features --features cmd_outliers
```

Do not run project tests.
