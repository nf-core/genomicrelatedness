#!/usr/bin/env python3
"""
Split a BED file into multiple BED files with approximately equal genomic load.

- Preserves genomic order
- Splits long contigs if necessary
- Limits number of intervals per BED file
"""

import argparse
from pathlib import Path


def split_intervals(
    intervals: list[tuple[str, int, int]],
    target_number_files: int,
    max_intervals: int,
) -> list[list[tuple[str, int, int]]]:
    """
    Split intervals into groups with approximately equal total base pairs.

    Parameters
    ----------
    intervals
        List of (contig, start, end) intervals in genomic order.
    target_bp
        Target number of base pairs per BED file.
    max_intervals
        Maximum number of intervals per BED file.

    Returns
    -------
    List of BED files, each represented as a list of intervals.
    """
    beds: list[list[tuple[str, int, int]]] = []
    current_bed: list[tuple[str, int, int]] = []
    current_bp = 0

    target_bp = sum(end - start for _, start, end in intervals) // target_number_files

    def flush() -> None:
        nonlocal current_bed, current_bp
        if current_bed:
            beds.append(current_bed)
            current_bed = []
            current_bp = 0

    for contig, start, end in intervals:
        pos = start

        while pos < end:
            remaining = end - pos
            capacity = target_bp - current_bp

            if capacity <= 0 or len(current_bed) >= max_intervals:
                flush()
                continue

            chunk_len = min(remaining, capacity)
            chunk_end = pos + chunk_len

            current_bed.append((contig, pos, chunk_end))
            current_bp += chunk_len
            pos = chunk_end

    flush()
    return beds

def read_bed(path: Path) -> list[tuple[str, int, int]]:
    """
    Read a BED file and return a list of intervals.

    :param path: Path to the BED file.
    :type path: Path
    :return: List of intervals as (contig, start, end) tuples.
    :rtype: list[tuple[str, int, int]]
    """
    intervals = []
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            contig, start, end = line.rstrip().split()[:3]
            intervals.append((contig, int(start), int(end)))
    return intervals

def write_beds(
    beds: list[list[tuple[str, int, int]]],
    out_prefix: str,
    pad_width: int,
) -> None:
    """
    Write multiple BED files from a list of BED intervals.

    :param beds: List of BED files, each represented as a list of intervals.
    :type beds: list[list[tuple[str, int, int]]]
    :param out_prefix: Prefix for output BED file names.
    :type out_prefix: str
    :param pad_width: Zero-padding width for file numbering.
    :type pad_width: int
    """
    for i, bed in enumerate(beds, start=1):
        name = f"{out_prefix}_{i:0{pad_width}d}.bed"
        with open(name, "w") as fh:
            for contig, start, end in bed:
                fh.write(f"{contig}\t{start}\t{end}\n")

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Split BED intervals into balanced chunks"
    )
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--target-number-files", required=True, type=int)
    parser.add_argument("--out-prefix", default="interval")
    parser.add_argument("--max-intervals", type=int, default=100)

    args = parser.parse_args()

    intervals = read_bed(args.bed)
    beds = split_intervals(
        intervals,
        target_number_files=args.target_number_files,
        max_intervals=args.max_intervals,
    )

    pad_width = len(str(len(beds)))
    write_beds(beds, args.out_prefix, pad_width)

if __name__ == "__main__":
    main()
