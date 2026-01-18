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
from velvet_optimizations import velvet_style_optimization, suggest_optimal_k
from read_correction import adaptive_correction


from io_fasta import fasta_to_list, write_fasta
from kmers import count_kmers, kmer_histogram, write_kmer_histogram
from dbg import build_graph_from_kmers, graph_stats, write_graph_stats
from cleaning import remove_islands, remove_tips
from cleaning import pop_bubbles_simple

from traversal import extract_contigs, contig_stats


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

    outdir = args.outdir

    # Output files (portfolio/debug friendly)
    params_path = outdir / "params_used.json"
    report_path = outdir / "report.txt"
    hist_path = outdir / "kmer_histogram.tsv"
    stats_before_path = outdir / "graph_stats_before_cleaning.tsv"
    stats_after_path = outdir / "graph_stats_after_cleaning.tsv"
    contigs_path = (
        outdir / "contigs.fasta"
    )  # will be used after traversal is implemented

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
    reads = fasta_to_list(args.input)
    report_lines.append(f"Reads loaded: {len(reads)}")
    if reads:
        report_lines.append(f"Read length (example): {len(reads[0])}")
    report_lines.append("")

    reads = adaptive_correction(reads, verbose=True)

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

    print("\n[3b] Velvet-style optimization...")
    graph = velvet_style_optimization(
        graph, 
        auto_cutoff=True,  # Auto coverage cutoff
        use_tourbus=True,  # Resolve repeats
        verbose=True
    )
    print("przed czyszczeniem")
    print(f"Nodes with out_degree > 0: {sum(1 for n in graph.nodes if graph.out_degree.get(n,0) > 0)}")
    print(f"Nodes with in_degree > 0: {sum(1 for n in graph.nodes if graph.in_degree.get(n,0) > 0)}")

    # ---------- STEP 4: CLEAN GRAPH ----------
    graph = remove_islands(graph, args.min_component_size)
    print(f"After islands: {len(graph.nodes)} nodes")

    graph = remove_tips(graph, args.tip_max_len)
    print(f"After tips: {len(graph.nodes)} nodes")
    # if args.pop_bubbles:
    graph = pop_bubbles_simple(graph, args.max_bubble_len)

        # Po build_graph:
    print("po czyszczeniem")
    print(f"Nodes with out_degree > 0: {sum(1 for n in graph.nodes if graph.out_degree.get(n,0) > 0)}")
    print(f"Nodes with in_degree > 0: {sum(1 for n in graph.nodes if graph.in_degree.get(n,0) > 0)}")

    stats_after = graph_stats(graph)
    write_graph_stats(stats_after, stats_after_path)
    report_lines.append("Graph stats after cleaning:")
    report_lines.append(f"  nodes: {stats_after.get('nodes', 'NA')}")
    report_lines.append(f"  edges: {stats_after.get('edges', 'NA')}")
    report_lines.append(f"Stats saved: {stats_after_path.name}")
    report_lines.append("")

    # ---------- STEP 5: TRAVERSAL TO CONTIGS (TODO) ----------
    # Suggested interface for portfolio-quality structure:
    # from traversal import extract_contigs
    contigs = extract_contigs(graph, k=args.kmer_length, min_contig_len=args.min_contig_len)
    write_fasta(contigs, contigs_path, prefix="contig")
    
    report_lines.append(f"Contigs written: {contigs_path.name}")
    report_lines.append(f"Contigs >= {args.min_contig_len} bp: {len(contigs)}")
    
    # For now, we only record that traversal is not yet implemented.
    contig_statistics = contig_stats(contigs)   

    for stat, value in contig_statistics.items():
        report_lines.append(f"{stat}: {value}")
    report_lines.append(f"Planned contigs path: {contigs_path.name}")
    report_lines.append("")

    write_report(report_path, report_lines)


if __name__ == "__main__":
    main()
