# logging Spec

CLI logging separates normal run narration from warnings and errors. This keeps pipelines quiet when requested while preserving failure visibility.

## User Modes

`--log` accepts:

- `stdout`: normal run narrative goes to stdout.
- `quiet`: normal run narrative and progress bars are suppressed.
- `file`: normal run narrative goes to an auto-generated log file under the command output directory.
- `file=<path>`: normal run narrative goes to the exact path.

Invalid values fail during CLI parsing. `file=` with an empty path is invalid.

## Routing

- The primary sink carries banners, lifecycle info, and summary/statistics blocks.
- Warnings and errors always go to stderr.
- In file mode, warnings and errors are also written to the log file.
- In quiet mode, warnings and errors still go to stderr.
- Progress bars draw to stderr only when stderr is a terminal and logging is not quiet.
- Non-TTY progress output must be hidden to avoid corrupting redirected logs.

## Formatting

- Stdout mode omits tracing levels for normal info lines.
- File mode includes tracing levels for info, warnings, and errors.
- Stderr warning/error lines include level and target.
- File mode uses the plain banner signature instead of terminal-oriented formatting.
- Auto-generated log files use:

```text
<output_dir>/logs/<command>_<YYYY-MM-DD_HH-MM-SS>_<random>.log
```

## Command Boundary

- The top-level binary initializes logging once before command execution.
- Nested command calls, such as `coverage-weights` calling internal `fcoverage`, must not print nested top-level banners.
- Nested calls should log phase messages with their own targets, so users can see where work is happening.
- `prep-windows` and `visualize-positions` currently force stdout logging because they do not expose shared `LoggingArgs`.

## Implementation Invariants

- `PRIMARY_OUTPUT` is initialized once per process. Reinitialization is an error.
- `write_primary` and `write_primary_line` are the only APIs that should write command narration outside tracing.
- Top-level errors are rendered to stderr and duplicated to file mode after the command footer.
- Progress factories must call `logging::progress_enabled()` so quiet mode is honored globally.
