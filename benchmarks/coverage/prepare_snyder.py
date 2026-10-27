#!/usr/bin/env python3
"""Reconstruct short pairs and downsample fragments from Snyder's original BAM.

Mapped primary single reads whose names start with ``M_`` are merged consensus
reads in the authors' processing. Represent each as two fully overlapping mates,
preserving its aligned fragment span and assigning proper-pair status during
reconstruction. Sample by read name so both mates and any additional alignments
are kept or discarded together. Keep other retained records unchanged, including
existing pair flags. Do not filter chromosomes, MAPQ, flags or fragment lengths.

For BH01's reported 96-fold coverage, --keep-probability 0.15625 targets roughly
15-fold coverage before remapping. Doubling a merged read creates two mates for
the same fragment, so it does not double the fragment depth.

This reconstructs pairs for fragment analysis. It does not recover the original
read sequences or their separate quality scores. Preparation is outside timing.
See https://shendurelab.github.io/cfDNA/#primary-data-processing-of-sequencing-data
and https://github.com/shendurelab/cfDNA/blob/master/FilterUniqueBAM.py
"""

import argparse
import copy
import hashlib
from pathlib import Path
import shlex
import sys

import pysam


def reconstruct_pair(consensus):
    """Represent a merged fragment as two mates with the same alignment span.

    The consensus becomes read1, retaining its strand. Read2 has the opposite
    strand. BAM stores bases and qualities in alignment orientation, so both
    records retain the same stored sequence, qualities and CIGAR.
    """
    fragment_length = consensus.reference_end - consensus.reference_start
    if fragment_length <= 0:
        raise ValueError(f"Merged read has no reference span: {consensus.query_name}")

    # Retain the consensus alignment and quality information for both mates
    # Both mates represent the same fragment on the same contig and opposite strands
    # Assign proper-pair status as a documented reconstruction convention
    # This does not recover the original aligner's classification
    read1 = copy.copy(consensus)
    read1.is_paired = True
    read1.is_proper_pair = True
    read1.is_read1 = True
    read1.is_read2 = False
    read1.mate_is_unmapped = False
    read1.mate_is_reverse = not read1.is_reverse
    read1.next_reference_id = read1.reference_id
    read1.next_reference_start = read1.reference_start

    # For mates covering identical coordinates, SAM permits either TLEN sign
    # convention as long as the signs differ. Use positive for read1
    read1.template_length = fragment_length
    read1.set_tag("MC", read1.cigarstring)
    read1.set_tag("MQ", read1.mapping_quality)

    # Alternative-alignment tags describe the original read's orientation
    # Remove them from reconstructed records rather than leave misleading tags
    for tag in ("SA", "XA"):
        if read1.has_tag(tag):
            read1.set_tag(tag, None)

    read2 = copy.copy(read1)
    read2.is_read1 = False
    read2.is_read2 = True
    read2.is_reverse = not read1.is_reverse
    read2.mate_is_reverse = read1.is_reverse
    read2.template_length = -fragment_length
    return read1, read2


def main():
    """Sample fragments, reconstruct retained short pairs, and index the new BAM.

    Hashing the seed and read name gives the same sampling decision wherever a
    fragment's records occur in the coordinate-sorted input. This avoids storing
    names in memory and makes selection reproducible across Python processes.
    Sampling before reconstruction avoids copying discarded short fragments.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_bam", type=Path, help="Original coordinate-sorted Snyder BAM")
    parser.add_argument("output_bam", type=Path, help="New downsampled BAM with reconstructed pairs")
    parser.add_argument(
        "--keep-probability", type=float, required=True,
        help="Probability of retaining each fragment, greater than 0 and at most 1",
    )
    parser.add_argument("--seed", type=int, default=42, help="Seed for reproducible fragment sampling")
    parser.add_argument("--threads", type=int, default=4, help="BAM I/O threads per reader/writer")
    arguments = parser.parse_args()
    if not 0 < arguments.keep_probability <= 1:
        parser.error("--keep-probability must be greater than 0 and at most 1")
    if arguments.threads < 1:
        parser.error("--threads must be at least 1")
    if arguments.output_bam.exists() or Path(f"{arguments.output_bam}.bai").exists():
        parser.error("The output BAM or BAI already exists. Choose a new output path")

    input_records = 0
    retained_input_records = 0
    merged_fragments_seen = 0
    reconstructed_fragments = 0
    # A 64-bit hash spans [0, 2**64). Keep the requested fraction of that range
    sampling_threshold = int(arguments.keep_probability * 2**64)
    print(f"Preparing {arguments.input_bam} with pysam {pysam.__version__}", file=sys.stderr)
    print(
        f"Fragment retention probability: {arguments.keep_probability}, seed: {arguments.seed}",
        file=sys.stderr,
    )
    with pysam.AlignmentFile(arguments.input_bam, "rb", threads=arguments.threads) as source:
        # Record the transformation in the BAM itself. Preserve contig names,
        # reference order, read groups and existing program descriptions
        header = source.header.to_dict()
        program_id = "prepare_snyder"
        if any(program["ID"] == program_id for program in header.get("PG", [])):
            parser.error("This BAM already has a prepare_snyder program entry")
        header.setdefault("PG", []).append(
            {"ID": program_id, "PN": "prepare_snyder.py", "VN": "2", "CL": shlex.join(sys.argv)}
        )
        header.setdefault("CO", []).append(
            "Mapped primary unpaired M_ consensus reads were reconstructed as overlapping mates. "
            "Consensus bases and qualities were copied. Proper-pair flags were assigned during "
            "reconstruction. Existing pair flags were preserved."
        )
        header["CO"].append(
            f"Fragments were sampled by a 64-bit BLAKE2b hash of seed:read_name, "
            f"with retention probability {arguments.keep_probability} and seed {arguments.seed}. "
            "No chromosome, MAPQ, flag or fragment length filters were applied."
        )
        with pysam.AlignmentFile(
            arguments.output_bam, "wb", header=header, threads=arguments.threads
        ) as destination:
            # Sequential iteration needs no input index and keeps memory bounded
            # New mates have the original position, so coordinate order is retained
            for record in source:
                input_records += 1
                if input_records % 100_000_000 == 0:
                    print(f"Processed {input_records:,} input records", file=sys.stderr)
                if record.query_name is None:
                    parser.error(f"Input record {input_records:,} has no read name for sampling")

                # M_ distinguishes merged fragments from ordinary singletons
                # Do not reconstruct secondary, supplementary or unmapped records
                is_merged_fragment = (
                    not record.is_paired
                    and not record.is_unmapped
                    and not record.is_secondary
                    and not record.is_supplementary
                    and record.query_name.startswith("M_")
                )
                if is_merged_fragment:
                    merged_fragments_seen += 1

                # Use only the seed and read name, never position, flags or quality
                # Both existing mates get the same decision even when far apart
                # A retained consensus produces both reconstructed mates below
                fragment_hash = hashlib.blake2b(
                    f"{arguments.seed}:{record.query_name}".encode("utf-8"), digest_size=8
                ).digest()
                if int.from_bytes(fragment_hash, "big") >= sampling_threshold:
                    continue
                retained_input_records += 1

                if is_merged_fragment:
                    read1, read2 = reconstruct_pair(record)
                    destination.write(read1)
                    destination.write(read2)
                    reconstructed_fragments += 1
                else:
                    destination.write(record)

    if merged_fragments_seen == 0:
        parser.error("No merged M_ fragments were found. Check that this is an original Snyder BAM")
    if retained_input_records == 0:
        parser.error("Sampling retained no records. Increase --keep-probability")

    # Use BAI because FinaleToolkit's coverage reader requires it
    # Indexing also reports an error if the input did not preserve coordinate order
    print("Creating the BAI index", file=sys.stderr)
    pysam.index("-@", str(arguments.threads), str(arguments.output_bam))
    print(f"Input records: {input_records:,}", file=sys.stderr)
    print(f"Retained input records: {retained_input_records:,}", file=sys.stderr)
    print(f"Reconstructed fragments: {reconstructed_fragments:,}", file=sys.stderr)
    print(f"Output records: {retained_input_records + reconstructed_fragments:,}", file=sys.stderr)
    print(f"Ready: {arguments.output_bam} and {arguments.output_bam}.bai", file=sys.stderr)


if __name__ == "__main__":
    main()
