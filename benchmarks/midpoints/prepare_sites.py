#!/usr/bin/env python3
"""Convert the supplied midpoint windows to cfDNAlab's and Griffin's input formats.

benchmark_midpoints.sh runs this once during setup, before any timed run, with the
setup Python. It only needs Python's standard library. Both tools receive the same
windows. Every site is kept, including repeated and overlapping rows. The native
tools check their own inputs, and this conversion adds no coordinate or strand
rules of its own. Only BED rows and the BAM header are read, no BAM reads.

Reads, in WORK_DIRECTORY:
- `input.bed`: the supplied grouped BED, with chromosome, start, end and group,
  optionally followed by a strand (BED5) or a score and a strand (BED6).
- `bam.chrom.sizes`: the BAM's contig names and lengths.

Writes, in WORK_DIRECTORY:
- `sites.bed`: BED6 for cfDNAlab, sorted by chromosome, start and end.
- `sites/<n>.tsv` and `sites.yaml`: Griffin's anchor table for each group, and the
  site-list YAML that names them.
- `chrom.sizes`: BAM lengths of the chromosomes with sites. Griffin uses it as its
  bigWig header.
- With `--griffin`, also the inputs for Griffin's untimed demo warmup
  (`warmup-sites.tsv`, `warmup-sites.yaml`, `warmup.chrom.sizes`) and
  `unit_gc_bias.tsv`.

Prints a tab-separated line with the window width, half width, number of sites,
number of groups and the comma-separated chromosomes with sites.
"""

import argparse
import csv
import json
from pathlib import Path


def read_bam_contigs(path):
    """Return the BAM's contig lengths from a two-column chrom.sizes file."""
    with path.open() as source:
        return {fields[0]: int(fields[1]) for fields in csv.reader(source, delimiter="\t")}


def read_sites(path):
    """
    Read the supplied BED rows as (chromosome, start, end, group, score, strand).

    BED4 has no strand, so its sites are unstranded ("."). Optional BED5/BED6
    strands follow cfDNAlab's convention: the fifth column of a BED5 row, or the
    sixth column of a BED6 row. Empty lines and "#" comment lines are skipped.
    """
    rows = []
    with path.open() as source:
        for fields in csv.reader(source, delimiter="\t"):
            if not fields or fields[0].startswith("#"):
                continue
            chromosome, start, end, group = fields[:4]
            strand = fields[5] if len(fields) >= 6 else fields[4] if len(fields) == 5 else "."
            rows.append((chromosome, int(start), int(end), group, 0, strand))
    if not rows:
        raise ValueError("No midpoint sites in the BED")
    return rows


def write_table(path, rows, header=None):
    """Write tab-separated rows, with an optional header row."""
    with path.open("w") as destination:
        writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
        if header:
            writer.writerow(header)
        writer.writerows(rows)


def write_griffin_sites(work, rows, half_width):
    """
    Write Griffin's anchor table for each group, and the YAML that lists them.

    Griffin describes a site by an anchor position and builds its window around it
    For odd widths Griffin uses [-floor(width/2), ceil(width/2)) around the anchor
    Its reverse-strand formula needs the other central boundary to recover the
    original BED bounds, so a "-" site is anchored at end - half_width. The native
    reversal then matches cfDNAlab's end - 1 - position indexing.

    Returns the number of groups.
    """
    groups = {}
    for chromosome, start, end, group, _score, strand in rows:
        anchor = end - half_width if strand == "-" else start + half_width
        groups.setdefault(group, []).append((chromosome, anchor, strand))
    (work / "sites").mkdir()
    site_paths = {}
    for group_index, (group, sites) in enumerate(groups.items()):
        path = work / "sites" / f"{group_index}.tsv"
        write_table(path, sites, header=("Chrom", "position", "Strand"))
        site_paths[group] = str(path)
    with (work / "sites.yaml").open("w") as destination:
        # JSON is valid YAML and avoids requiring PyYAML in the setup Python
        json.dump({"site_lists": site_paths}, destination)
    return len(groups)


def write_griffin_demo(work, bam_contigs):
    """
    Write the inputs for Griffin's untimed demo warmup.

    The release's fixed 1000-site CTCF demo gives Griffin a short warmup without
    changing the supplied production windows. It uses the same copied BAM and
    reference as the measured runs.
    """
    demo_path = work / "griffin-env/share/griffin/CTCF.hg38.1000.txt"
    demo_rows = []
    with demo_path.open() as source:
        for row in csv.DictReader(source, delimiter="\t"):
            demo_rows.append((row["Chrom"], int(row["position"]), "."))
    if len(demo_rows) != 1000:
        raise ValueError("Expected exactly 1000 sites in the pinned Griffin demo")
    write_table(work / "warmup-sites.tsv", demo_rows, header=("Chrom", "position", "Strand"))
    with (work / "warmup-sites.yaml").open("w") as destination:
        json.dump({"site_lists": {"CTCF-demo": str(work / "warmup-sites.tsv")}}, destination)
    demo_chromosomes = sorted(set(row[0] for row in demo_rows))
    write_table(
        work / "warmup.chrom.sizes",
        [(chromosome, bam_contigs[chromosome]) for chromosome in demo_chromosomes],
    )


def write_unit_gc_bias(work, min_length, max_length):
    """
    Write a GC bias table with weight 1 for every fragment length and GC count.

    Only the untimed Griffin demo uses these unit weights. Measured repetitions fit
    their own sample GC bias and keep both the raw and the corrected profile tables.
    """
    rows = [
        (fragment_length, gc_count, 1.0)
        for fragment_length in range(min_length, max_length + 1)
        for gc_count in range(fragment_length + 1)
    ]
    write_table(work / "unit_gc_bias.tsv", rows, header=("length", "num_GC", "smoothed_GC_bias"))


def main():
    """Prepare both tools' site inputs and print the site summary."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work_directory", type=Path)
    parser.add_argument("--griffin", action="store_true", help="Also prepare Griffin's demo warmup")
    parser.add_argument("--min-length", type=int, required=True)
    parser.add_argument("--max-length", type=int, required=True)
    arguments = parser.parse_args()
    work = arguments.work_directory

    bam_contigs = read_bam_contigs(work / "bam.chrom.sizes")
    rows = read_sites(work / "input.bed")
    # The supplied fixed-size windows define the profile width. Odd widths keep a
    # central base, e.g. 2001 positions span -1000 through +1000
    window_width = rows[0][2] - rows[0][1]
    half_width = window_width // 2
    rows.sort(key=lambda row: (row[0], row[1], row[2]))
    write_table(work / "sites.bed", rows)
    group_count = write_griffin_sites(work, rows, half_width)

    # Griffin sorts intervals lexicographically. Its bigWig header uses that same
    # contig order, so the exported entries stay correctly ordered
    chromosomes = sorted(set(row[0] for row in rows))
    write_table(work / "chrom.sizes", [(chromosome, bam_contigs[chromosome]) for chromosome in chromosomes])

    if arguments.griffin:
        write_griffin_demo(work, bam_contigs)
        write_unit_gc_bias(work, arguments.min_length, arguments.max_length)

    print(window_width, half_width, len(rows), group_count, ",".join(chromosomes), sep="\t")


if __name__ == "__main__":
    main()
