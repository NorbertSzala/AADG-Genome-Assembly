from pathlib import Path
import json

from io_fasta import fasta_to_list, write_fasta
from kmers import count_kmers, kmer_histogram, write_kmer_histogram
from dbg import build_graph_from_kmers, graph_stats, write_graph_stats
from cleaning import remove_islands, remove_tips, pop_bubbles_simple
from traversal import extract_contigs, contig_stats

from velvet_optimizations import velvet_style_optimization
from read_correction import adaptive_correction


###############################
########## HELPERS  ###########
###############################
def write_params_json( kmer_length, outdir, min_kmer_count, min_component_size,
                       tip_max_len,pop_bubbles, max_bubble_len, min_contig_len, path: Path) -> None:
    payload = {
        "input": str(input),
        "outdir": str(outdir),
        "kmer_length": kmer_length,
        "min_kmer_count": min_kmer_count,
        "min_component_size": min_component_size,
        "tip_max_len": tip_max_len,
        "pop_bubbles": pop_bubbles,
        "max_bubble_len": max_bubble_len,
        "min_contig_len": min_contig_len,
    }
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)


def write_report(path: Path, lines: list[str]) -> None:
    with open(path, "w") as f:
        for line in lines:
            f.write(line.rstrip() + "\n")

###############################
########## FUNCTIONS ##########
###############################

def run_assembly(
    input: Path,
    outdir: Path,
    kmer_length: int = 31,
    min_kmer_count: int = 2,
    min_component_size: int = 10,
    tip_max_len: int = 3,
    pop_bubbles: bool = False,
    max_bubble_len: int = 5,
    min_contig_len: int = 300,
    optimize: bool = False
) -> dict:
    """
    Run genome assembly pipeline.
    
    Args:
        input_path: Path to input FASTA file
        outdir: Output directory
        kmer_length: K-mer length (default: 31)
        min_kmer_count: Minimum k-mer count (default: 2)
        min_component_size: Minimum component size (default: 10)
        tip_max_len: Maximum tip length (default: 3)
        pop_bubbles: Enable bubble popping (default: False)
        max_bubble_len: Maximum bubble length (default: 5)
        min_contig_len: Minimum contig length (default: 300)
        optimize: Run optimization mode (default: False)
    
    Returns:
        dict: Contig statistics including total_length, n50, longest, etc.
    """

    # ===============================================================
    # NORMAL ASSEMBLY MODE
    # ===============================================================
    
    outdir = outdir

    # Output files
    params_path = outdir / "params_used.json"
    report_path = outdir / "report.txt"
    hist_path = outdir / "kmer_histogram.tsv"
    stats_before_path = outdir / "graph_stats_before_cleaning.tsv"
    stats_after_path = outdir / "graph_stats_after_cleaning.tsv"
    contigs_path = outdir / "contigs.fasta"

    write_params_json(kmer_length, outdir, min_kmer_count, min_component_size,
                       tip_max_len,pop_bubbles, max_bubble_len, min_contig_len, params_path)

    report_lines: list[str] = []
    report_lines.append("Genome assembly pipeline report")
    report_lines.append(f"Input: {input}")
    report_lines.append(f"Output dir: {outdir}")
    report_lines.append(f"k: {kmer_length}")
    report_lines.append(f"min_kmer_count: {min_kmer_count}")
    report_lines.append(f"min_component_size: {min_component_size}")
    report_lines.append(f"tip_max_len: {tip_max_len}")
    report_lines.append(f"pop_bubbles: {pop_bubbles}")
    if pop_bubbles:
        report_lines.append(f"max_bubble_len: {max_bubble_len}")
    report_lines.append("")

    # ---------- STEP 1: READ FASTA ----------
    print("[1/6] Reading reads...")
    reads = fasta_to_list(input)
    report_lines.append(f"Reads loaded: {len(reads)}")
    if reads:
        report_lines.append(f"Read length (example): {len(reads[0])}")
    report_lines.append("")

    reads = adaptive_correction(reads)

    # ---------- STEP 2: K-MERS ----------
    kmer_counts = count_kmers(reads, kmer_length)
    histogram = kmer_histogram(kmer_counts)
    write_kmer_histogram(histogram, hist_path)
    report_lines.append(f"Distinct k-mers (before filtering): {len(kmer_counts)}")
    report_lines.append(f"Histogram saved: {hist_path.name}")
    report_lines.append("")

    # ---------- STEP 3: BUILD DBG ----------
    graph = build_graph_from_kmers(
        kmer_counts=kmer_counts,
        k=kmer_length,
        min_kmer_count=min_kmer_count,
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
    graph = remove_islands(graph, min_component_size)
    print(f"    After islands: {len(graph.nodes)} nodes")

    print("  Removing tips...")
    graph = remove_tips(graph, tip_max_len)
    print(f"    After tips: {len(graph.nodes)} nodes")
    
    if pop_bubbles:
        print("  Popping bubbles...")
        graph = pop_bubbles_simple(graph, max_bubble_len)
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
        k=kmer_length, 
        min_contig_len=min_contig_len
    )
    
    print(f"  Generated {len(contigs)} contigs >= {min_contig_len} bp")
    
    # ---------- STEP 6: WRITE OUTPUT ----------
    print("[6/6] Writing output...")
    write_fasta(contigs, contigs_path, prefix="contig")
    
    report_lines.append(f"Contigs written: {contigs_path.name}")
    report_lines.append(f"Contigs >= {min_contig_len} bp: {len(contigs)}")
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

    return contig_statistics
