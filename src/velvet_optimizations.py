#!/usr/bin/env python3
"""
Velvet-style optimizations for De Bruijn Graph assembly

Key improvements:
1. Coverage cutoff - automatic error removal based on k-mer coverage
2. Tour Bus algorithm - resolve simple repeats
3. Read threading - use read coherence information
"""

from collections import Counter, defaultdict
from dbg import DeBruijnGraph


###############################
#### COVERAGE CUTOFF ##########
###############################

def estimate_coverage_cutoff(graph: DeBruijnGraph, percentile: float = 0.2) -> int:
    """
    Estimate coverage cutoff automatically using k-mer coverage distribution.
    
    Theory:
    - Erroneous k-mers have low coverage (1-2×)
    - Real k-mers have coverage around genome coverage
    - Cutoff should be between these two peaks
    
    Args:
        graph: DeBruijnGraph with edge weights
        percentile: Use this percentile of coverage distribution (0.2 = 20th percentile)
    
    Returns:
        Suggested coverage cutoff value
    """
    # Collect all edge weights (k-mer counts)
    coverages = []
    for node, targets in graph.out_edges.items():
        for target, weight in targets.items():
            coverages.append(weight)
    
    if not coverages:
        return 2  # Default fallback
    
    # Sort and find percentile
    coverages.sort()
    cutoff_idx = int(len(coverages) * percentile)
    cutoff = coverages[cutoff_idx]
    
    # Don't go below 2
    return max(2, cutoff)

def apply_coverage_cutoff(graph: DeBruijnGraph, min_coverage: int) -> DeBruijnGraph:
    edges_to_remove = []
    
    # 1. Znajdź słabe krawędzie
    for u, targets in graph.out_edges.items():
        for v, weight in targets.items():
            if weight < min_coverage:
                edges_to_remove.append((u, v))
    
    print(f"  Removing {len(edges_to_remove)} low-coverage edges (< {min_coverage}×)")
    
    # 2. Usuń krawędzie i aktualizuj metadane
    for u, v in edges_to_remove:
        if u in graph.out_edges and v in graph.out_edges[u]:
            del graph.out_edges[u][v]
            
            # POPRAWKA: Bezpieczna dekrementacja
            if u in graph.out_degree:
                graph.out_degree[u] -= 1
            if v in graph.in_degree:
                graph.in_degree[v] -= 1
    
    # 3. CZYSZCZENIE: Usuń izolowane węzły
    nodes_to_check = set()
    for u, v in edges_to_remove:
        nodes_to_check.add(u)
        nodes_to_check.add(v)
    
    for node in nodes_to_check:
        in_d = graph.in_degree.get(node, 0)
        out_d = graph.out_degree.get(node, 0)
        
        if in_d == 0 and out_d == 0:
            graph.out_edges.pop(node, None)
            graph.in_degree.pop(node, None)
            graph.out_degree.pop(node, None)
            
            if hasattr(graph, 'nodes'):
                graph.nodes.discard(node)
    
    return graph




# ###############################
# #### TOUR BUS ALGORITHM #######
# ###############################
def _walk_with_coverage(graph, start_node, max_len):
    path = [start_node]
    total_cov = 0
    curr = start_node
    
    for step in range(max_len):
        # Sprawdzaj stopnie tylko ПОСЛЕ pierwszego węzła
        if step > 0:  # ← ZMIANA
            if graph.out_degree.get(curr, 0) != 1:
                break
            if graph.in_degree.get(curr, 0) > 1:
                break
        
        next_nodes = list(graph.out_edges.get(curr, {}).keys())
        if not next_nodes:
            break
        
        nxt = next_nodes[0]
        total_cov += graph.out_edges[curr][nxt]
        path.append(nxt)
        curr = nxt
    
    return path, total_cov


def tour_bus_optimized(graph, max_repeat_length=15, max_iterations=5):
    total_resolved = 0
    
    for iteration in range(max_iterations):
        resolved = 0
        
        for start in list(graph.out_edges.keys()):
            if start not in graph.out_edges:  # ← DODANE
                continue
            
            targets = list(graph.out_edges[start].keys())
            if len(targets) < 2:
                continue
            
            for i in range(len(targets)):
                for j in range(i + 1, len(targets)):
                    v1, v2 = targets[i], targets[j]
                    
                    # Check if still exists
                    if v1 not in graph.out_edges.get(start, {}) or \
                       v2 not in graph.out_edges.get(start, {}):
                        continue
                    
                    path1, cov1 = _walk_with_coverage(graph, v1, max_repeat_length)
                    path2, cov2 = _walk_with_coverage(graph, v2, max_repeat_length)
                    
                    if path1 and path2 and path1[-1] == path2[-1]:
                        to_remove_start = v1 if cov1 < cov2 else v2
                        
                        # Bezpieczne usunięcie
                        if to_remove_start in graph.out_edges[start]:
                            del graph.out_edges[start][to_remove_start]
                            
                            if start in graph.out_degree:
                                graph.out_degree[start] -= 1
                            if to_remove_start in graph.in_degree:
                                graph.in_degree[to_remove_start] -= 1
                            
                            resolved += 1
        
        total_resolved += resolved
        
        if resolved == 0:
            break
    
    if total_resolved > 0:
        print(f" [Tour Bus] Resolved {total_resolved} bubbles (in {iteration+1} iterations)")
    
    return graph
# ###############################
# #### AUTO K-MER SELECTION #####
# ###############################

def suggest_optimal_k(
    read_length: int,
    genome_size: int = 10000,
    coverage: float = None
) -> list[int]:
    """
    Suggest optimal k values to try, Velvet-style.
    
    Heuristics:
    - k should be odd (canonical k-mers)
    - k < read_length (otherwise no overlaps)
    - Higher coverage → can use larger k
    - Try multiple k values
    
    Args:
        read_length: Average read length
        genome_size: Estimated genome size
        coverage: Estimated coverage (optional)
    
    Returns:
        List of k values to try
    """
    max_k = min(read_length - 10, 41)  # Leave room for overlaps
    min_k = 15
    
    # Suggest 3-5 values
    if coverage and coverage > 15:
        # High coverage - try larger k
        k_values = [max_k - 4, max_k - 2, max_k]
    elif coverage and coverage < 8:
        # Low coverage - try smaller k
        k_values = [min_k, min_k + 4, min_k + 8]
    else:
        # Medium coverage - try range
        k_values = [min_k + 4, min_k + 10, max_k - 6]
    
    # Ensure odd values
    k_values = [k if k % 2 == 1 else k - 1 for k in k_values]
    
    return sorted(set(k_values))



# ###############################
# #### INTEGRATED PIPELINE ######
# ###############################

def velvet_style_optimization(
    graph: DeBruijnGraph,
    auto_cutoff: bool = True,
    use_tourbus: bool = True,
    verbose: bool = True
) -> DeBruijnGraph:
    """
    Apply Velvet-style optimizations to graph.
    
    Pipeline:
    1. Auto-detect coverage cutoff
    2. Remove low-coverage edges
    3. Tour Bus for short repeats
    
    Args:
        graph: Input graph
        auto_cutoff: Automatically estimate and apply coverage cutoff
        use_tourbus: Use Tour Bus algorithm
        verbose: Print progress
    
    Returns:
        Optimized graph
    """
    if verbose:
        print("\n=== Velvet-style Optimization ===")
    
    # 1. Coverage cutoff
    if auto_cutoff:
        if verbose:
            print("[1] Auto coverage cutoff...")
        
        cutoff = estimate_coverage_cutoff(graph, percentile=0.2)
        if verbose:
            print(f"  Estimated cutoff: {cutoff}×")
        
        graph = apply_coverage_cutoff(graph, min_coverage=cutoff)
    
    # 2. Tour Bus
    if use_tourbus:
        if verbose:
            print("[2] Tour Bus algorithm...")
        
        graph = tour_bus_optimized(graph, max_repeat_length=10)
    
    if verbose:
        print("=== Optimization complete ===\n")
    
    return graph