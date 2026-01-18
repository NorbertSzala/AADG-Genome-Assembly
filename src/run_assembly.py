#!/usr/bin/env python3
"""
Summary:
Main genome assembly pipeline

Pipeline:
1) Read FASTA reads -> list[str]
2) Count k-mers + build k-mer histogram (TSV)
3) Build De Bruijn graph from filtered k-mers
4) Graph cleaning: islands, tips, optional simple bubbles
5) Traversal to contigs
6) Outputs: stats TSV, report TXT, params JSON

"""

###############################
########### IMPORTS  ##########
###############################
from __future__ import annotations

from pathlib import Path
import argparse

from assembly_core import run_assembly

###############################
########## ARGUMENTS ##########
###############################
def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="run_assembly.py",
        description=(
            "De novo assembly from single-end reads (same strand, single chromosome). "
            "Pipeline: FASTA -> k-mers -> DBG -> cleaning -> contigs."
        ),
    )
    parser.add_argument(
        "-I",
        "--input",
        required=True,
        type=Path,
        help="Input reads FASTA file (single-end).",
    )
    parser.add_argument(
        "-O",
        "--outdir",
        required=True,
        type=Path,
        help="Output directory for results (will be created if missing).",
    )
    parser.add_argument(
        "-K",
        "--kmer-length",
        type=int,
        default=31,
        help="K-mer length (positive integer). Default: 31",
    )
    parser.add_argument(
        "-M",
        "--min-kmer-count",
        type=int,
        default=2,
        help="Minimum k-mer count to keep when building DBG. Default: 2",
    )
    parser.add_argument(
        "--min-component-size",
        type=int,
        default=10,
        help="Minimum size (nodes) of connected component to keep. Default: 10",
    )
    parser.add_argument(
        "--tip-max-len",
        type=int,
        default=3,
        help="Maximum tip length (nodes) to remove. Default: 3",
    )
    parser.add_argument(
        "--pop-bubbles",
        action="store_true",
        help="Enable conservative bubble popping (simple bubbles only).",
    )
    parser.add_argument(
        "--max-bubble-len",
        type=int,
        default=5,
        help="Maximum bubble path length for simple bubble popping. Default: 5",
    )
    parser.add_argument(
        "--min-contig-len",
        type=int,
        default=300,
        help="Minimum contig length to report. Default: 300",
    )
    
    parser.add_argument(
        "--optimize",
        action="store_true",
        help="Run full parameter optimization (tests ~120 combinations, ~0.5 hours)",
    )

    args = parser.parse_args()

    # ----- validation -----
    if args.kmer_length <= 0:
        raise ValueError("K-mer length must be > 0")

    if args.min_kmer_count <= 0:
        raise ValueError("min-kmer-count must be > 0")

    if args.min_component_size <= 0:
        raise ValueError("min-component-size must be > 0")

    if args.tip_max_len < 0:
        raise ValueError("tip-max-len must be >= 0")

    if args.max_bubble_len <= 0:
        raise ValueError("max-bubble-len must be > 0")

    if args.min_contig_len <= 0:
        raise ValueError("min-contig-len must be > 0")

    if not args.input.exists():
        raise ValueError(f"Input FASTA file does not exist: {args.input}")

    if not args.input.is_file():
        raise ValueError(f"Input path is not a file: {args.input}")
    
    args.outdir.mkdir(parents=True, exist_ok=True)

    return args


###############################
######## MAIN LOGIC  ##########
###############################

def main() -> None:
    """Command-line interface wrapper for run_assembly()"""
    args = parse_arguments()

    if args.optimize:
        from optimize_parameters import optimize_parameters
        optimize_parameters(args.input, args.outdir)

    else:
        run_assembly(
            input=args.input,
            outdir=args.outdir,
            kmer_length=args.kmer_length,
            min_kmer_count=args.min_kmer_count,
            min_component_size=args.min_component_size,
            tip_max_len=args.tip_max_len,
            pop_bubbles=args.pop_bubbles,
            max_bubble_len=args.max_bubble_len,
            min_contig_len=args.min_contig_len,
            optimize=args.optimize
        )


if __name__ == "__main__":
    main()