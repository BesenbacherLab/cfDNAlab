#!/usr/bin/env python3
"""Pair all unpaired records and downsample fragments from a Snyder BAM.

Represent every retained unpaired record as two mates, regardless of its name,
MAPQ, or secondary, supplementary, duplicate, QC-fail or unmapped status. Preserve
those filtering flags on both mates. Mapped reconstructed pairs cover the same
aligned reference span. Unmapped reconstructed pairs remain unmapped.

Sample by read name so both mates and any additional alignments are kept or
discarded together. Keep existing paired records unchanged. Do not filter
chromosomes, MAPQ, flags or fragment lengths.

For BH01's reported 96-fold coverage, --keep-probability 0.15625 targets roughly
15-fold coverage. Doubling a merged read creates two mates for
the same fragment, so it does not double the fragment depth.

This reconstructs pairs for fragment analysis. It does not recover the original
read sequences or their separate quality scores. For a singleton that is not a
complete merged fragment, the reconstructed span describes only the observed
alignment, not the missing original fragment boundaries. Preparation is outside
timing.
See https://shendurelab.github.io/cfDNA/#primary-data-processing-of-sequencing-data
and https://github.com/shendurelab/cfDNA/blob/master/FilterUniqueBAM.py
"""

import argparse
import copy
import hashlib
from pathlib import Path
import sys

import pysam


def reconstruct_pair(unpaired_record):
    """Represent an unpaired record as two mates, preserving its filtering flags.

    The original record becomes read1, retaining its strand. Read2 has the opposite
    strand. BAM stores bases and qualities in alignment orientation, so both
    records retain the same stored sequence, qualities and CIGAR.

    Unmapped records have no aligned reference span. Their reconstructed mates
    remain unmapped, with TLEN zero and proper-pair status unset.
    """
    fragment_length = 0
    if not unpaired_record.is_unmapped:
        reference_end = unpaired_record.reference_end
        if reference_end is None or reference_end <= unpaired_record.reference_start:
            raise ValueError(f"Mapped read has no reference span: {unpaired_record.query_name}")
        fragment_length = reference_end - unpaired_record.reference_start

    # Copy the record so secondary, supplementary, duplicate, QC-fail and unmapped
    # flags survive on both mates. Only pairing-related fields are changed
    # Assign proper-pair status to mapped pairs as a reconstruction convention
    # This does not recover the original aligner's classification
    read1 = copy.copy(unpaired_record)
    read1.is_paired = True
    read1.is_proper_pair = not unpaired_record.is_unmapped
    read1.is_read1 = True
    read1.is_read2 = False
    read1.mate_is_unmapped = unpaired_record.is_unmapped
    read1.mate_is_reverse = not read1.is_reverse
    read1.next_reference_id = read1.reference_id
    read1.next_reference_start = read1.reference_start

    # For mates covering identical coordinates, SAM permits either TLEN sign
    # convention as long as the signs differ. Use positive for read1
    # Unmapped pairs keep TLEN zero. A missing CIGAR removes any stale MC tag
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
    """Sample fragments, pair every retained unpaired record, and index the BAM.

    Hashing the seed and read name gives the same sampling decision wherever a
    fragment's records occur in the coordinate-sorted input. This avoids storing
    names in memory and makes selection reproducible across Python processes.
    Sampling before reconstruction avoids copying discarded records.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_bam", type=Path, help="Coordinate-sorted Snyder BAM")
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
    reconstructed_pairs = 0
    # A 64-bit hash spans [0, 2**64). Keep the requested fraction of that range
    sampling_threshold = int(arguments.keep_probability * 2**64)
    print(f"Preparing {arguments.input_bam} with pysam {pysam.__version__}", file=sys.stderr)
    print(
        f"Fragment retention probability: {arguments.keep_probability}, seed: {arguments.seed}",
        file=sys.stderr,
    )
    with pysam.AlignmentFile(arguments.input_bam, "rb", threads=arguments.threads) as source:
        # Preserve contig names, reference order, read groups and existing metadata
        # Do not add a processing-history entry or embed command-line file paths
        # Metadata already present in the input header is retained unchanged
        header = source.header.to_dict()
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

                # Use only the seed and read name, never position, flags or quality
                # Both existing mates get the same decision even when far apart
                # Every retained unpaired record produces two mates below
                fragment_hash = hashlib.blake2b(
                    f"{arguments.seed}:{record.query_name}".encode("utf-8"), digest_size=8
                ).digest()
                if int.from_bytes(fragment_hash, "big") >= sampling_threshold:
                    continue
                retained_input_records += 1

                # Pair every unpaired record, with no name or status restrictions
                # Downstream tools decide which alignment flags and MAPQ to accept
                if not record.is_paired:
                    read1, read2 = reconstruct_pair(record)
                    destination.write(read1)
                    destination.write(read2)
                    reconstructed_pairs += 1
                else:
                    destination.write(record)

    if retained_input_records == 0:
        parser.error("Sampling retained no records. Increase --keep-probability")

    # Use BAI because FinaleToolkit's coverage reader requires it
    # Indexing also reports an error if the input did not preserve coordinate order
    print("Creating the BAI index", file=sys.stderr)
    pysam.index("-@", str(arguments.threads), str(arguments.output_bam))
    print(f"Input records: {input_records:,}", file=sys.stderr)
    print(f"Retained input records: {retained_input_records:,}", file=sys.stderr)
    print(f"Reconstructed pairs: {reconstructed_pairs:,}", file=sys.stderr)
    print(f"Output records: {retained_input_records + reconstructed_pairs:,}", file=sys.stderr)
    print(f"Ready: {arguments.output_bam} and {arguments.output_bam}.bai", file=sys.stderr)


if __name__ == "__main__":
    main()
