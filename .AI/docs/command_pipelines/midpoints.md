# `cfdna midpoints`

Build grouped midpoint profiles from a BAM file. The command counts fragment midpoints inside grouped windows and writes a profile tensor that can be loaded directly in NumPy.

## Pipeline

```mermaid
flowchart TD
    inputs["Input BAM and grouped BED<br/>coordinate-sorted BAM, fixed-size grouped windows"]
    setup["Prepare profile axes<br/>groups, positions, and fragment length bins"]
    interval_filter["Prefilter intervals<br/>chromosomes and blacklist safety margin"]
    tiles["Scan interval tiles<br/>with fragment-length halos"]
    fragments["Create fragments<br/>pair mates or use each read as a fragment"]
    midpoint["Place midpoint<br/>reproducible tie handling for even fragments"]
    filters["Keep usable midpoint events<br/>MAPQ, pairing, length, blacklist, and window hit"]
    weights["Apply requested weights<br/>GC tags or GC package, plus scaling TSVs"]
    count["Count into sparse tile profiles<br/>group x length bin x position"]
    merge["Merge sparse tile profiles"]
    postprocess["Post-process final profiles<br/>configured smoothing, trim flanks, final binning"]
    profile["Dense midpoint profile<br/>midpoint_profiles.npy"]
    group_index["Group index<br/>group_idx to group name"]
    settings["Settings<br/>axis, smoothing, and binning metadata"]
    plots["QC plots<br/>selected groups when plotting is enabled"]
    stats["Run statistics<br/>fragments counted, skipped, and filtered"]

    correction_inputs["Optional correction inputs<br/>GC package, 2bit reference, GC tag, scaling TSVs"]
    filter_inputs["Optional filter inputs<br/>blacklists and chromosome selection"]

    inputs --> setup --> interval_filter --> tiles --> fragments --> midpoint --> filters --> weights --> count --> merge --> postprocess --> profile --> stats
    setup --> group_index
    setup --> settings
    postprocess --> settings
    profile --> plots

    correction_inputs -.-> weights
    filter_inputs -.-> tiles
    filter_inputs -.-> interval_filter
    filter_inputs -.-> filters

    classDef core fill:#eef5ff,stroke:#3b73b9,color:#10233f;
    classDef optional fill:#f7f7f4,stroke:#777,color:#202020;
    classDef outputClass fill:#e9f8ef,stroke:#3e8f57,color:#102a17;

    class inputs,setup,interval_filter,tiles,fragments,midpoint,filters,weights,count,merge,postprocess core;
    class correction_inputs,filter_inputs,plots optional;
    class profile,group_index,settings,stats outputClass;
```

## Profile Model

Each grouped BED row contributes a fixed-width window. Windows with the same group name are collapsed into one profile, so the final array is shaped `(group, length_bin, position)`. Fragment length bins come from `--length-bins`, and each accepted fragment contributes to the bin containing its fragment length.

Savitzky-Golay smoothing is applied to final grouped profiles by default, not to each input interval. Use `--smoothing none` when you need unsmoothed profiles. Optional binning averages adjacent positions after smoothing and after computation flanks are trimmed.

## Midpoint Placement

For odd-length fragments, the midpoint is the center base. For even-length fragments, the command reproducibly assigns the midpoint to one of the two central bases so the profile does not always round ties in the same direction.

## Outputs

The main output is `<prefix>.midpoint_profiles.npy`. The companion `<prefix>.group_index.tsv` maps numeric group indices back to group names. The settings sidecar `<prefix>.midpoint_profile_settings.json` records the length axis, position binning, and smoothing settings needed to interpret the array. When plotting is enabled, selected groups also produce quick-look PNG profiles and length-bin heatmaps.

## Corrections And Filters

Optional blacklists remove fragments before counting. By default, blacklist input also removes intervals whose output span plus the fragment/smoothing safety margin overlaps a blacklisted region. Use `--keep-blacklisted-intervals` to keep those intervals while preserving fragment-level blacklist filtering.

Optional GC correction and genome scaling change each midpoint's count weight before it is added to the profile.
