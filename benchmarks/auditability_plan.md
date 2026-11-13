# Benchmark auditability plan

## Goal and scope

A reviewer should be able to follow an experiment from its native command to its timing row and numerical comparison without understanding a general execution system. The main path should read as: prepare inputs, run warmups, measure repetitions, save outputs, compare.

Work only in `benchmarks/`. This is a refactor plan, not approval to change experimental settings. Follow the decisions in [native_usage_review.md](native_usage_review.md). Raise conflicting or unclear decisions in chat. Do not preserve unnecessary interfaces just for backwards compatibility, and do not add compatibility branches without discussion.

Some readers are necessary because the packages produce different formats. The aim is to make each conversion short, explicit and independently checkable, while removing unnecessary execution machinery.

## What must remain scientifically unchanged

- Keep the supplied BAM, indexes, references and site BED unchanged. Do not regenerate or select production sites during this refactor.
- Keep separate jobs for 1, 4, 12 and 32 cores. Pass N to documented native thread options. Keep explicitly selected serial tools in the comparison. Preserve the agreed warmup/repetition schedules, pilot resource policy and Griffin demo exception.
- Keep the approved native calls, version pins, filtering differences and timing boundaries. Picard timing is its standalone command. Its filtered histogram is produced outside timing for concordance. Griffin timing uses `griffin-gc`, retaining both raw and corrected tables.
- Keep base-resolution coverage MAE, weighted Pearson, residual histograms and Jaccard diagnostics. A histogram contains residual values and their base counts, not a residual row for every genomic position.
- Keep midpoint normalization over all 666 bins before comparison of the central 132 bins. Keep native smoothing order, 15 bp bins, MAPQ 20 and the distinction between cfDNAlab eligible intervals and Griffin listed sites.
- Keep the documented motif side approximation, cfDNAPro reference-contig restriction and native filtering differences. Do not make profiles agree by silently changing scope, shifting, smoothing or dropping positions.

## Main changes

| Current obstacle | Proposed change |
| --- | --- |
| Long `STEP` heredocs in four runners | Ordinary executable worker scripts, called directly by `srun` |
| Repeated CPU-placement and measurement code | A small shared implementation for these mechanical operations |
| `prefix`, numbered `fallback[...]` fields and temporary-file handshakes | One explicit record directory per invocation, with records written by field name |
| Generic helpers that contain tool-specific output rules | Output paths and preservation rules beside the relevant native call |
| Readers infer formats and comparison groups from name fragments | Explicit method-to-reader and quantity assignments within each feature |
| Tests extract and execute source fragments using `split()` | Invoke the actual worker with fake executables and check its arguments and files |
| Several copies of commands and settings in documentation | Executable command definitions are authoritative, with a concise comparison table linking to them |

Splitting the existing code across more files is insufficient. Remove the old indirection and duplicated rules as each replacement is completed.

## 1. Establish a short experiment contract

Before refactoring, reconcile the command report and specifications with the approved review. Make one concise table per feature containing:

- Method and measured quantity.
- The executable command definition, including multi-command workflows.
- Shared filters, native differences and scope.
- Work inside and outside the timer.
- Native output used for comparison, its units and coordinate convention.
- Any normalization or conversion applied afterward.

This is the review baseline. Do not treat the old implementation or its tests as proof that an option is scientifically justified. In particular, tests must not reinstate undocumented options such as the removed `bedcov` threading workaround.

Keep `experiment_commands.md` as the reviewer-facing explanation. Avoid maintaining a second full copy of executable recipes that can drift from the scripts. Link directly to the authoritative command definitions.

## 2. Simplify one complete feature first

Start with lengths. Finish its execution, saved outputs, comparison and tests before applying the structure to end motifs, midpoints and finally coverage.

- Keep the submission launcher responsible for resource requests. Keep each feature driver responsible for staging, selected methods and repetition loops.
- Replace each embedded step with a normal `<feature>/run_method.sh` executable. Give it explicit invocation arguments, including its run directory. Existing environment variables can remain the job settings.
- Put literal native calls in clearly named cases or functions. Keep a short multi-command workflow together. Retain the existing small R and Griffin helpers where they make those calls easier to inspect.
- Prepare a command once and execute that same command. Use Bash arrays for argument boundaries. Do not build a command generator, plugin registry or nested configuration language.
- Put each method's output filename, units and effective filter metadata beside its command. Generic timing code must not guess those details from a directory name.

For lengths, an auditor must immediately see that timed Picard reads the original BAM and concordance reads the separately filtered histogram.

## 3. Make measurement and saving explicit

Reduce `common.sh` to mechanical operations shared by the four drivers. It should not decide biological filters, select native output tables or assign quantities from method-name substrings.

- Create an invocation record before starting the tool. Record method, repetition, warmup status, settings and job identity explicitly.
- Use the passed run directory for timing, status and logs instead of shared `time.tsv`, `io.tsv`, `profile.path`, `stdout.path` and `step.id` handshakes.
- Assemble the final timing row from named fields in one small writer. Remove positional prefixes and numbered fallback arrays. Preserve useful output columns and their meanings without adding a migration system.
- Keep tool exit status, Slurm step status and missing measurements distinct. Handle expected command failure explicitly. File-saving errors must remain visible and fatal.
- Save each completed row immediately. Keep earlier measurements and outputs if a later tool fails. Never infer success from a missing status.
- Retain small stdout/stderr logs for successful and failed attempts. Do not duplicate stdout when it is the native data output, such as samtools stats.
- Preserve the exact worker/helper source and invocation arguments with results. Record expanded command arguments from the command actually executed, rather than from a separately maintained description.
- Keep staging, output preservation and comparison outside tool timing. Preserve the separate flush measurement and exact Slurm step IDs for memory accounting.

Keep one measured native profile per method as now. Delete redundant large outputs only after saving the required record and profile. Do not retain every genome-wide output just to simplify directory handling.

## 4. Separate native readers from metric calculations

Keep comparison after measurement, including automatic comparison where it currently applies. Each feature should have a readable sequence: load named outputs, check their declared meaning, align the agreed coordinates or labels, calculate metrics, write the small result tables.

- Use the cfDNAlab reader package and standard table/compression libraries where available. Do not reimplement those file formats.
- Give each native reader one explicit input format. Document the columns read, units, coordinate conversion and treatment of absent values beside that reader.
- Keep format conversion separate from normalization and metrics. For example, reading PanDepth coordinates should not also choose a comparison subset.
- Distinguish numeric zero, an empty field, undefined statistics and malformed data explicitly. Never use truthiness as a numeric missing-value test.
- Preserve legitimate smoothing undershoots. Keep documented NA/Inf results distinct from unexpected nonfinite input or arithmetic errors.
- Put formulas and their assumptions next to short metric functions. Retain useful shared numerical checks without introducing a configurable metrics framework.
- Make comparison groups explicit: positional read depth, fragment coverage, regional read depth and regional fragment counts are different quantities.

Coverage needs a streaming comparison because the genome is too large for dense whole-genome arrays. Keep the base-weighted sweep over constant-depth runs and the explicit reconciliation of contig order. The MAE must remain `sum(abs(left - right) * span_width) / compared_bases`, including declared zero positions.

### Residual storage decision

SQLite currently bounds memory even for fractional positional values. Replacing it with an ordinary Counter could simplify the current uncorrected integer-depth comparison, but support for fractional positional tracks and the storage guarantee remain an unresolved D7 decision.

Do not silently make that change. First verify the output types of the approved calls and discuss restricting positional support with the owner. Regional fractional fragment counts remain supported regardless. Until that decision, keep the current storage behavior isolated from the numerical calculation. Do not add interchangeable storage backends.

## 5. Replace brittle tests with small, meaningful checks

- Keep independent, manually derived numerical expectations. Do not regenerate expectations from the implementation being tested.
- Test real worker entry points with fake native executables. Check the actual arguments, timed work, untimed work and preserved files. Stop extracting code between comments or heredoc markers.
- Use tiny native-format fixtures with documented columns. Where feasible, retain a small real native output and its source version to check the reader assumptions.
- Keep explicit cases for numeric zero versus blank metadata, empty and constant profiles, coordinate boundaries, sparse zero filling, PanDepth's missing 1 bp tail, and full-window midpoint normalization before cropping.
- BAM fixtures must respect the real CLI constraints, including the 10 bp minimum possible fragment length and the experiment's inclusion filters.
- Keep focused failure tests showing that completed timing rows survive later failures. Do not recreate a general Slurm simulator in every feature's tests.
- Follow the repository's testing rules. The owner runs tests unless explicitly authorizing otherwise. Report syntax checks, compilation checks, executed tests and native validation separately.

Actual tool validation remains necessary. Use retained native outputs or the first scheduled jobs where appropriate. Do not require another expensive full Griffin run merely to validate the refactor, and do not assume a tiny synthetic BAM can support its GC estimator.

## Completion criteria

For each feature, demonstrate that a reviewer can answer these questions directly from the code:

1. What exact native command runs, and where did each nondefault option come from?
2. Which operations enter its reported time and Slurm memory measurement?
3. Which native output is compared, and how are coordinates and units converted?
4. Where are scope selection, normalization and metric formulas applied?
5. What survives a failed run, and how is failure distinguished from a valid zero?

Show a before/after trace for one experiment per feature. Reproduce at least one metric per feature with a short independent calculation from its saved native output, without calling the comparator being checked. Explain which bookkeeping and duplicated rules were removed, not just where they moved. Retire the old execution path as each feature is replaced.

Leave installer changes, summary-utility redesign, new metrics, new tool modes and unrelated data preparation outside this first refactor. Update the current specifications to describe the final structure. Keep existing review notes and unresolved plans available for the owner's review.
