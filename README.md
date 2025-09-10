# cfDNAlab

Toolkit for cfDNA analysis written in rust for *speed*.


## Commands
The following commands are currently available:

| Command         | Description                                                         | Minimal example                                                              |
| --------------- | ------------------------------------------------------------------- | ---------------------------------------------------------------------------- |
| `cfdna lengths` | Count fragment lengths (defined as `end(reverse) - start(forward)`) | `cfdna lengths --bam <> --output-dir <> --blacklist <> --assign-by midpoint` |
 

### Common options

 - **Windowing**: Perform the command in windows. Either a single global window (usually default), the windows in a BED given file, or via a window size. Assign fragments/reads/... to windows by how they overlap.
 - **Blacklist filtering**: Supply BED files with regions to blacklist. The implementation is specific to each tool (filtering of full fragments or just the overlapping positions).


## TODO

    - Bin chromosomes for higher parallelization where meaningful.