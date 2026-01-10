#!/usr/bin/env python3

"""
Summary:
    Generate k-mers from reads and compute k-mer frequency statistics.

Input:
    - list of read sequences (strings, A/C/G/T only)
    - integer k

Output:
    - k-mer count histogram written to TSV file

    Interpretation
    count	n_kmers
    1	20340       # there are 20 340 distinct kmers that appear exactly once
    2	1653
    .
    .
    .
    10	18
    11	1           # there is just one unique kmer that appears eleven times

"""

###############################
########### IMPORTS  ##########
###############################
from pathlib import Path
from tqdm import tqdm
from collections import defaultdict
import argparse


###############################
########## ARGUMENTS ##########
###############################
def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="kmers.py",
        description="Script generates kmers, counts them and creates histogram (how many times every kmers occurs) -> .tsv",
    )
    parser.add_argument(
        "-I",
        "--input",
        required=True,
        type=Path,
        help="Path with input to your fasta files",
    )
    parser.add_argument(
        "-O",
        "--output",
        required=True,
        type=Path,
        help="Path where to save .tsv file with counted frequencies of kmers",
    )
    parser.add_argument(
        "-K", "--kmer-length", required=True, type=int, help="Kmer length"
    )

    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    # check if user provide proper extension, unless, make it properly
    if args.output.suffix.lower() != ".tsv":
        args.output = args.output.with_suffix(".tsv")

    if args.kmer_length <= 0:
        raise ValueError("K-mer length must be a positive integer")

    if not args.input.exists():
        raise ValueError(f"Input FASTA file does not exist: {args.input}")

    if not args.input.is_file():
        raise ValueError(f"Input path is not a file: {args.input}")

    return args


###############################
########## FUNCTIONS ##########
###############################
def iter_kmers(read: str, k: int):
    """
    Generate all kmers from a single read

    Args:
        read: nucleotide sequence
        k: kmer length
    Yields:
        succesive kmers as strings
    """
    read_len = len(read)
    if read_len < k:
        return

    for i in range(read_len - k + 1):
        yield read[i : i + k]


def count_kmers(reads: list[str], k: int) -> dict[str, int]:
    """
    Count k-mer frequencies across all reads

    args:
        reads: list of nucleotide seqeunces
        k: kmer length

    Returns: dictionary mapping kmer -> occurence count
    """

    # I do not understand this line. Explain me why irstly we paste variable name and then ':' ? And why convert kmer_counts into dict 'dict(kmer_counts) if kmer_counts was declated as dict earlier?
    kmer_counts: dict[str, int] = defaultdict(int)

    if reads == []:
        raise ValueError("Reads is empty list")

    # what is iter_kmers?
    for read in tqdm(reads, desc="Counting k-mers"):
        if k > len(read):
            continue

        for kmer in iter_kmers(read, k):
            kmer_counts[kmer] += 1
    return dict(kmer_counts)


def kmer_histogram(kmer_counts: dict[str, int]) -> dict[int, int]:
    """
    Build histogram of kmer counts

    args:
        kmer_counts: dict mapping kmer (str) -> count (int)

    returns:
        Dict mapping count -> number of kmers with that count
    """

    histogram: dict[int, int] = defaultdict(int)

    for count in kmer_counts.values():
        histogram[count] += 1

    return dict(histogram)


def write_kmer_histogram(histogram: dict[int, int], path: Path) -> None:
    """
    Save kmer count histogram to TSV file
    Format:
        count \t number_of_kmers
    """

    with open(path, "w") as handle:
        handle.write("count\tn_kmers\n")
        for count in sorted(histogram):
            handle.write(f"{count}\t{histogram[count]}\n")


###############################
######## MAIN LOGIC  ##########
###############################
def main():
    args = parse_arguments()

    from io_fasta import fasta_to_list

    reads = fasta_to_list(args.input)

    k = 31  # potem przesun to do osobnego pliku z konfiguracjami. I tak trzeba bedzie iterować i sprawdzac najlepsze opcje

    kmer_counts = count_kmers(reads, k)
    histogram = kmer_histogram(kmer_counts)

    write_kmer_histogram(histogram, args.output)


if __name__ == "__main__":
    main()
