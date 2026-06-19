# Final Output Writes Spec

Final user-facing outputs are completion signals for workflow managers and ad hoc pipelines. Commands must not expose those outputs at their final paths until the command has finished writing every file in the output set.

## Final Output Set

- A final output is any user-facing file whose existence can be interpreted as command completion.
- Primary outputs, compressed tables, BED or bedGraph outputs, group-index files, settings files, headers, BAM files, package files, and other metadata needed to interpret a primary output are part of the final output set.
- Internal tile files and reducer scratch files are not final outputs. They stay in command temp directories and are not public completion signals.

## Write Contract

- Each final output must be written to a temporary path on the same filesystem as its final output directory.
- The temporary path must be under a command-owned temp directory, not beside the final file under a normal output-like name.
- Writers must be flushed and closed before the temp path is moved into place.
- Commands must record every temp-to-final path pair after the temp file has been fully written.
- Commands must move recorded files into place only after all final outputs in the set have been written successfully.
- If a final-output write fails before the move step, previously existing final paths should be left untouched when possible.
- Move failures are errors because they leave the command without the requested completed output set.

## Shared Helpers

Commands should use `shared::io::FinalOutputFiles` for final output placement:

- `FinalOutputFiles::new` creates the `final_outputs` subdirectory inside the supplied command temp directory.
- `temp_path_for` maps a final output path to the corresponding path in the final-output temp directory.
- `record` stores one temp-to-final path pair and rejects duplicate temp or final paths.
- `record_temp_files_with_same_names_in` records files created inside the final-output temp directory when a writer returns the paths it created.
- `move_into_place` renames all recorded temp files to their final paths.

The helper enforces path bookkeeping. Callers are still responsible for writing complete files, closing writers, and recording all required metadata files before calling `move_into_place`.

## Temp Directories and Cleanup

- Final-output temp directories live inside a command-owned temp directory.
- Command temp directories are unique per run and should be created under the selected output directory so temporary and final files stay on the expected filesystem.
- `TempDirGuard` owns cleanup on success, early return, and drop.
- Normal cleanup is best effort. Drop-time cleanup failures are warnings, not command failures after final outputs have been moved into place.
- Commands should call `TempDirGuard::remove` only when cleanup failure should become part of the command result for that specific path.

## Limits

- This contract is not a full transaction across multiple final files.
- A process kill during `move_into_place` can still leave only part of the output set moved.
- The contract does not change existing output filenames, schemas, or formats.
