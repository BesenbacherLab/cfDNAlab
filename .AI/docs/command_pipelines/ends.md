# `cfdna ends`

Count fragment end motifs from a BAM file. The command extracts bases just outside and inside each fragment end, orients the motif consistently, and writes motif counts per genome window or group.

## Pipeline

```mermaid
flowchart TD
    input["Input BAM<br/>coordinate-sorted and indexed"]
    setup["Define motif outputs<br/>inside/outside bases and windows"]
    tiles["Scan genome in tiles<br/>with fragment-length halos"]
    fragments["Create fragments<br/>pair mates or use each read as a fragment"]
    ends["Resolve fragment ends<br/>aligned, clipped, or raw-boundary mode"]
    filters["Keep usable end events<br/>MAPQ, length, indels, base quality, and blacklist rules"]
    windows["Assign motifs to rows<br/>global, fixed bins, BED windows, or groups"]
    weights["Apply requested weights<br/>GC correction and scaling TSVs"]
    count["Count end motifs<br/>outside_inside labels per output row"]
    merge["Merge tile counts"]
    output_choice{"Motif output"}
    sparse["Sparse observed motifs<br/>end_motifs.sparse.npz"]
    dense["Dense full motif set<br/>end_motifs.npy"]
    labels["Motif labels and settings<br/>end_motifs.txt and JSON"]
    metadata["Window metadata<br/>bins.tsv or group_index.tsv"]
    stats["Run statistics<br/>fragments and end motifs counted"]

    region_inputs["Optional region inputs<br/>BED windows, grouped BED, blacklists"]
    correction_inputs["Optional correction inputs<br/>2bit reference, GC package or tag, scaling TSVs"]

    input --> setup --> tiles --> fragments --> ends --> filters --> windows --> weights --> count --> merge --> output_choice
    output_choice -- "observed motifs" --> sparse
    output_choice -- "--all-motifs" --> dense
    sparse --> labels
    dense --> labels
    labels --> metadata --> stats

    region_inputs -.-> setup
    region_inputs -.-> filters
    region_inputs -.-> windows
    correction_inputs -.-> ends
    correction_inputs -.-> weights

    classDef core fill:#eef5ff,stroke:#3b73b9,color:#10233f;
    classDef optional fill:#f7f7f4,stroke:#777,color:#202020;
    classDef decision fill:#fff5d9,stroke:#b88700,color:#2f2500;
    classDef outputClass fill:#e9f8ef,stroke:#3e8f57,color:#102a17;

    class input,setup,tiles,fragments,ends,filters,windows,weights,count,merge core;
    class region_inputs,correction_inputs optional;
    class output_choice decision;
    class sparse,dense,labels,metadata,stats outputClass;
```

## Motif Model

For each kept fragment end, `ends` combines the requested outside-reference bases with the requested inside-fragment bases. Inside bases can come from the reads or from the reference. Right-end motifs are reverse-complemented so both ends use the same fragment-end-inward 5' to 3' orientation. Motif labels are written as `<outside>_<inside>`.

## Window Model

Without windowing, all counted motifs go into one global row. With fixed bins or BED windows, rows represent genomic intervals. With grouped BED input, rows represent group names, and windows with the same group are aggregated.

The default assignment counts each motif in the row containing its endpoint. Other modes can assign both end motifs by fragment overlap, midpoint, or overlap proportion.

## Outputs

By default, the command writes only observed motif columns in a sparse `.npz` matrix. With `--all-motifs`, it writes a dense `.npy` matrix containing every possible motif for the selected inside/outside lengths. The motif label file and settings JSON describe how to interpret the matrix, while `bins.tsv` or `group_index.tsv` describes the output rows when windowing is used.
