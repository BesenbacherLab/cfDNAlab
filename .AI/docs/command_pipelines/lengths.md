# `cfdna lengths`

Count fragment length distributions from a BAM file. The command turns alignments into fragments, assigns each kept fragment to output rows, and writes a length-count matrix that can be loaded directly in NumPy.

## Pipeline

```mermaid
flowchart TD
    input["Input BAM<br/>coordinate-sorted and indexed"]
    setup["Define length outputs<br/>length bins and optional windows"]
    tiles["Scan genome in tiles<br/>with fragment-length halos"]
    fragments["Create fragments<br/>pair mates or use each read as a fragment"]
    length_model["Resolve fragment length<br/>aligned, indel-adjusted, or clip-adjusted"]
    filters["Keep usable fragments<br/>MAPQ, pairing, length, blacklist, indel, and clip rules"]
    rows["Assign fragments to rows<br/>global, fixed bins, BED windows, or groups"]
    weights["Apply requested weights<br/>GC correction and scaling TSVs"]
    count["Count length bins<br/>one vector per output row"]
    merge["Merge tile counts<br/>including windows crossing tile edges"]
    matrix["Length count matrix<br/>length_counts.npy"]
    settings["Length settings<br/>fragment_length_settings.json"]
    metadata["Window metadata<br/>bins.tsv or group_index.tsv"]
    plot["QC plot<br/>overall length distribution"]
    stats["Run statistics<br/>fragments counted, skipped, and filtered"]

    region_inputs["Optional region inputs<br/>BED windows, grouped BED, fixed bins, blacklist"]
    correction_inputs["Optional correction inputs<br/>2bit reference, GC package, scaling TSVs"]

    input --> setup --> tiles --> fragments --> length_model --> filters --> rows --> weights --> count --> merge --> matrix --> stats
    setup --> settings
    rows --> metadata
    matrix --> plot

    region_inputs -.-> setup
    region_inputs -.-> filters
    region_inputs -.-> rows
    correction_inputs -.-> weights

    classDef core fill:#eef5ff,stroke:#3b73b9,color:#10233f;
    classDef optional fill:#f7f7f4,stroke:#777,color:#202020;
    classDef outputClass fill:#e9f8ef,stroke:#3e8f57,color:#102a17;

    class input,setup,tiles,fragments,length_model,filters,rows,weights,count,merge core;
    class region_inputs,correction_inputs,plot optional;
    class matrix,settings,metadata,stats outputClass;
```

## Length Model

For paired-end BAM input, `lengths` builds fragments from inward-facing mates and uses the fragment span from the forward read position to the reverse read reference end. In `--reads-are-fragments` mode, each accepted read is treated as its own fragment.

The counted length can be the aligned fragment length, an indel-adjusted length, or a clip-adjusted length. Length bins are half-open intervals, and each accepted fragment contributes to the bin containing its resolved length.

## Window Model

Without windowing, all counted fragments go into one global row. With fixed bins or BED windows, rows represent genomic intervals. With grouped BED input, rows represent group names, and windows with the same group are aggregated.

The default assignment counts a fragment in every overlapping row. Other assignment modes can require full containment, use the fragment midpoint, or weight rows by overlap proportion.

## Outputs

The main output is `<prefix>.length_counts.npy`, with one row per output bin or group and one column per length bin. The settings JSON records the length-bin and counting configuration. Windowed runs also write `bins.tsv` or `group_index.tsv` so matrix rows can be mapped back to genomic intervals or group names. When plotting is enabled, the command also writes an overall length-distribution PNG.
