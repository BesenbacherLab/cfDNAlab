# Historical outlier simulations: methods and selected outcomes

Recorded 2026-10-19 from the local historical report dated 2026-10-13 and its settings file.
The settings record base seed 20261014 and Git revision
`f5a752117861a4e53cb790b087fdd10557b5a728`. These dates are reproduced as recorded.
No simulations were rerun for this note. This tracked summary replaces links requiring the
ignored `outliers_fragment_simulations/` directory. The full research artifacts remain local.

The fitted-model follow-up contained 228 scenarios with three seeds each, or 684 paired
comparisons. Replicate seed = `20261014 + 100 * scenario_index + repetition`, with
zero-based scenario index and repetitions 0, 1 and 2. Both methods received the same simulated
fragments, masks and production-selected contexts. Ordinary positional fragment coverage ranged
from 0.1 to 200. This is not nominal sequenced-base depth. Underlying expected rates used gamma or lognormal fields, with rate CV
0, 0.05, 0.1, 0.2 or 0.3. CV is rate standard deviation divided by mean, before sampling.
Fields used 1 kb blocks, with selected 200 bp and 10 kb comparisons. Fragments were 167 bp.
Biological examples included requested 2.35-fold gains spanning 10 kb, 100 kb or 1 Mb.
Technical spikes added 500 or 5,000 fragments starting within 30 bp, including spikes inside gains.

ZIP used the repository fitter through a Rust harness. ZINB used a Python/SciPy 1.9.3
implementation with profiled zero inflation, three deterministic L-BFGS-B starts, an explicit
ZIP boundary and a converged Powell fallback when all gradient starts failed. It fitted both
stages with truncation normalization and received no known generating mean or CV. This is not
the current Rust Nelder-Mead optimizer. The recorded Python fitter SHA-256 is
`e97e11b462f82edd202839167de9ff4bc1b92dc9051c008ea58e15e29b97ba9e`.

Settings were alpha `1 / 2850000000`, active-component mean target, 100 kb preservation windows,
50 kb denominator floor, **10% cutoff**, 1,000 bp flanks and whole-fragment minimum weights.
Outcomes below are means over three replicates, reproduced from the report.

| Scenario | ZIP | ZINB |
| --- | --- | --- |
| Ordinary mean 200, lognormal rate CV 0.2 | 2.81% ordinary signal removed | 0% |
| 1 Mb gain, ordinary mean 17, rate CV 0.1 | 9.34% biological signal removed | 0% |
| 10 kb gain, ordinary mean 17, rate CV 0.1 | 36.06% biological signal removed | 18.05% |
| 10 kb gain, ordinary mean 200, rate CV 0.1 | 0% biological signal removed | 58.18% |
| 500-fragment spike in ordinary mean 200, CV 0.1 | 68.48% technical signal removed | 69.35% |
| 5,000-fragment spike inside 1 Mb gain, same background | 0% technical signal removed | 94.91% |

In the 58.18% loss case, the first replicate's window score was 10.063% for ZIP and 9.983%
for ZINB. Their opposite sides of the old cutoff explain the preservation difference.
In the large spike inside a gain, ZINB removed about 0.032% of biological gain signal.
Separating biological and technical mass was possible because simulated fragment origins were known.

The report records no ordinary-background calls under ZINB in its finite background replicates.
Separate fitted population-tail calculations still showed model dependence. At mean 200 and
rate CV 0.2, projected called ordinary positions over 2.85 billion positions were about 2.4
for gamma block interiors and 406 for lognormal block interiors. These are population calculations,
not measured false-call counts in an actual genome. An 8 Mb simulation cannot empirically
establish an extreme-tail probability near `1e-10`.

Of 32,836 optimizer attempts, 656 reported failure. Four fits required independently converged
Powell optimization. The record reports inclusive-threshold and retained-likelihood checks and
eight component-coverage replays. No project tests were run. The ignored research harness and
complete scenario/seed files are required for exact reproduction and are not provided by this note.

These outcomes motivate the [validation list](outliers_real_data_validation.md). They do not
validate the current Rust optimizer, current 5% preservation rule, or any biological release criterion.
