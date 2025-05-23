#!/usr/bin/env python3

import sys
from Bio import SeqIO


def sort_fasta(input_file, output_file):
    # Read sequences
    records = list(SeqIO.parse(input_file, "fasta"))

    # Normalize headers
    # strip string polypolish if present
    for record in records:
        record.id = record.id.replace("polypolish", "")
        record.description = record.description.replace("polypolish", "")
        record.id = record.id.strip()
        record.description = record.description.strip()

    # Sort by length, reverse for largest to smallest
    records.sort(key=lambda x: len(x.seq), reverse=True)

    # append length to sequence description
    for record in records:
        record.description += f" length={len(record.seq)}"

    # Write sorted sequences
    SeqIO.write(records, output_file, "fasta")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python sort_fasta.py input.fasta output.fasta")
        sys.exit(1)

    sort_fasta(sys.argv[1], sys.argv[2])
