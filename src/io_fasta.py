#!/usr/bin/env python3

"""
Summary:
Script reading a FASTA file and returning a list of sequences as strings.

Args:
    -I / --input   Path to input FASTA file
    -O / --output  Path to output FASTA file

Input:
    FASTA file with nucleotide sequences

Output:
    FASTA file written from parsed sequences (for testing / passthrough)
"""

###############################
########### IMPORTS  ##########
###############################
from Bio import SeqIO
from pathlib import Path
import argparse


###############################
########## ARGUMENTS ##########
###############################
def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="io_fasta.py", description="Read fasta and return list of sequences"
    )
    parser.add_argument(
        "-I",
        "--input",
        type=Path,
        required=True,
        help="Path to your .fasta file with reads",
    )

    parser.add_argument(
        "-O",
        "--output",
        type=Path,
        required=True,
        help="Path to output. Produce single file",
    )

    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)

    # check if user provide proper extension, unless, make it properly
    if args.output.suffix.lower() != ".fasta":
        args.output = args.output.with_suffix(".fasta")

    if not args.input.exists():
        raise ValueError(f"Input FASTA file does not exist: {args.input}")

    if not args.input.is_file():
        raise ValueError(f"Input path is not a file: {args.input}")

    return args


###############################
########## FUNCTIONS ##########
###############################
def normalize_sequences(seq: str) -> str:
    """Make all nucleotides upper case and remove blank characters"""
    return seq.upper().strip()


def fasta_to_list(path: Path) -> list[str]:
    """
    Read fasta file, normalize it and return list of sequences
    """
    sequences: list[str] = []

    for seq_record in SeqIO.parse(path, "fasta"):
        seq = normalize_sequences(str(seq_record.seq))

        if not seq:
            continue

        sequences.append(seq)

    if not sequences:
        raise ValueError(f"No valid sequences found in FASTA file: {path}")

    return sequences


def write_fasta(seqs: list[str], path: Path, prefix: str = "seq") -> None:
    with open(path, "w") as handle:
        for i, seq in enumerate(seqs, start=1):
            header = f">{prefix}_{i}"
            handle.write(header + "\n")
            handle.write(seq + "\n")


###############################
######## MAIN LOGIC  ##########
###############################
def main(input: Path, output: Path) -> None:
    sequences = fasta_to_list(input)
    write_fasta(sequences, output)


if __name__ == "__main__":
    args = parse_arguments()
    main(args.input, args.output)
