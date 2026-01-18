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
import json
from datetime import datetime

from io_fasta import fasta_to_list, write_fasta
from kmers import count_kmers, kmer_histogram, write_kmer_histogram
from dbg import build_graph_from_kmers, graph_stats, write_graph_stats
from cleaning import remove_islands, remove_tips, pop_bubbles_simple
from traversal import extract_contigs, contig_stats

from velvet_optimizations import velvet_style_optimization
from read_correction import adaptive_correction
from optimize_parameters import optimize_parameters

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
        help="Run full parameter optimization (tests ~120 combinations, ~2 hours)",
    )
    parser.add_argument(
        "--optimize-quick",
        action="store_true",
        help="Run quick parameter optimization (tests ~15 combinations, ~15 minutes)",
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
########## HELPERS  ###########
###############################
def write_params_json(args: argparse.Namespace, path: Path) -> None:
    payload = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "input": str(args.input),
        "outdir": str(args.outdir),
        "kmer_length": args.kmer_length,
        "min_kmer_count": args.min_kmer_count,
        "min_component_size": args.min_component_size,
        "tip_max_len": args.tip_max_len,
        "pop_bubbles": args.pop_bubbles,
        "max_bubble_len": args.max_bubble_len,
        "min_contig_len": args.min_contig_len,
    }
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)


def write_report(path: Path, lines: list[str]) -> None:
    with open(path, "w") as f:
        for line in lines:
            f.write(line.rstrip() + "\n")


###############################
######## MAIN LOGIC  ##########
###############################
def main() -> None:
    args = parse_arguments()
    
    # ===============================================================
    # OPTIMIZATION MODE - Run parameter search instead of assembly
    # ===============================================================
    if args.optimize or args.optimize_quick:
        print("="*70)
        print("RUNNING IN OPTIMIZATION MODE")
        print("="*70)
        print(f"Input: {args.input}")
        print(f"Output: {args.outdir}")
        print()
        
        if args.optimize_quick:
            print("Mode: FULL optimization (~120 combinations, ~2 hours)")
            optimize_parameters(args.input, final_output_dir=str(args.outdir))
        
        print("\n✓ Optimization complete!")
        print(f"✓ Best assembly saved to: {args.outdir}/")
        return  # Exit - don't run normal assembly
    
    # ===============================================================
    # NORMAL ASSEMBLY MODE
    # ===============================================================
    
    outdir = args.outdir

    # Output files
    params_path = outdir / "params_used.json"
    report_path = outdir / "report.txt"
    hist_path = outdir / "kmer_histogram.tsv"
    stats_before_path = outdir / "graph_stats_before_cleaning.tsv"
    stats_after_path = outdir / "graph_stats_after_cleaning.tsv"
    contigs_path = outdir / "contigs.fasta"

    write_params_json(args, params_path)

    report_lines: list[str] = []
    report_lines.append("Genome assembly pipeline report")
    report_lines.append(f"Input: {args.input}")
    report_lines.append(f"Output dir: {outdir}")
    report_lines.append(f"k: {args.kmer_length}")
    report_lines.append(f"min_kmer_count: {args.min_kmer_count}")
    report_lines.append(f"min_component_size: {args.min_component_size}")
    report_lines.append(f"tip_max_len: {args.tip_max_len}")
    report_lines.append(f"pop_bubbles: {args.pop_bubbles}")
    if args.pop_bubbles:
        report_lines.append(f"max_bubble_len: {args.max_bubble_len}")
    report_lines.append("")

    # ---------- STEP 1: READ FASTA ----------
    print("[1/6] Reading reads...")
    reads = fasta_to_list(args.input)
    report_lines.append(f"Reads loaded: {len(reads)}")
    if reads:
        report_lines.append(f"Read length (example): {len(reads[0])}")
    report_lines.append("")

    reads = adaptive_correction(reads)

    # ---------- STEP 2: K-MERS ----------
    kmer_counts = count_kmers(reads, args.kmer_length)
    histogram = kmer_histogram(kmer_counts)
    write_kmer_histogram(histogram, hist_path)
    report_lines.append(f"Distinct k-mers (before filtering): {len(kmer_counts)}")
    report_lines.append(f"Histogram saved: {hist_path.name}")
    report_lines.append("")

    # ---------- STEP 3: BUILD DBG ----------
    graph = build_graph_from_kmers(
        kmer_counts=kmer_counts,
        k=args.kmer_length,
        min_kmer_count=args.min_kmer_count,
    )

    stats_before = graph_stats(graph)
    write_graph_stats(stats_before, stats_before_path)
    report_lines.append("Graph stats before cleaning:")
    report_lines.append(f"  nodes: {stats_before.get('nodes', 'NA')}")
    report_lines.append(f"  edges: {stats_before.get('edges', 'NA')}")
    report_lines.append(f"Stats saved: {stats_before_path.name}")
    report_lines.append("")


    graph = velvet_style_optimization(
        graph, 
        auto_cutoff=True,
        use_tourbus=True
    )

    # Debug info
    print(f"Before cleaning: {len(graph.nodes)} nodes")

    # ---------- STEP 4: CLEAN GRAPH ----------
    print("[4/6] Cleaning graph...")
    
    print("  Removing islands...")
    graph = remove_islands(graph, args.min_component_size)
    print(f"    After islands: {len(graph.nodes)} nodes")

    print("  Removing tips...")
    graph = remove_tips(graph, args.tip_max_len)
    print(f"    After tips: {len(graph.nodes)} nodes")
    
    if args.pop_bubbles:
        print("  Popping bubbles...")
        graph = pop_bubbles_simple(graph, args.max_bubble_len)
        print(f"    After bubbles: {len(graph.nodes)} nodes")

    stats_after = graph_stats(graph)
    write_graph_stats(stats_after, stats_after_path)
    report_lines.append("Graph stats after cleaning:")
    report_lines.append(f"  nodes: {stats_after.get('nodes', 'NA')}")
    report_lines.append(f"  edges: {stats_after.get('edges', 'NA')}")
    report_lines.append(f"Stats saved: {stats_after_path.name}")
    report_lines.append("")

    # ---------- STEP 5: TRAVERSAL TO CONTIGS ----------
    print("[5/6] Extracting contigs...")
    contigs = extract_contigs(
        graph, 
        k=args.kmer_length, 
        min_contig_len=args.min_contig_len
    )
    
    print(f"  Generated {len(contigs)} contigs >= {args.min_contig_len} bp")
    
    # ---------- STEP 6: WRITE OUTPUT ----------
    print("[6/6] Writing output...")
    write_fasta(contigs, contigs_path, prefix="contig")
    
    report_lines.append(f"Contigs written: {contigs_path.name}")
    report_lines.append(f"Contigs >= {args.min_contig_len} bp: {len(contigs)}")
    report_lines.append("")
    
    # Contig statistics
    contig_statistics = contig_stats(contigs)   
    for stat, value in contig_statistics.items():
        report_lines.append(f"{stat}: {value}")
    
    write_report(report_path, report_lines)
    
    print("\n" + "="*70)
    print("ASSEMBLY COMPLETE")
    print("="*70)
    print(f"Contigs: {len(contigs)}")
    if contigs:
        print(f"Longest: {contig_statistics['longest']} bp")
        print(f"Total length: {contig_statistics['total_length']} bp")
        print(f"N50: {contig_statistics['n50']} bp")
    print(f"\nResults in: {outdir}/")
    print("="*70)


if __name__ == "__main__":
    main()