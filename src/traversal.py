#!/usr/bin/env python3
from dbg import DeBruijnGraph

"""
traversal.py

Summary:
Extract contigs from a cleaned De Bruijn graph via traversal.
"""

def find_start_nodes(graph: DeBruijnGraph) -> list[str]:
    """
    Find nodes that should be starting points for contig paths.
    
    A node is a start if: in_degree != 1 OR out_degree != 1
    """
    start_nodes = []
    
    for node in graph.nodes:
        in_deg = graph.in_degree.get(node, 0)
        out_deg = graph.out_degree.get(node, 0)
        
        # Not a simple path node (which would have in=1, out=1)
        if in_deg != 1 or out_deg != 1:
            start_nodes.append(node)
    
    return start_nodes

def walk_unitig(
    graph: DeBruijnGraph,
    start_node: str,
    used_edges: set[tuple[str, str]]):
    """
    Walk along an unambiguous path (unitig) starting from start_node.
    """
    path = [start_node]
    edges_used = set()
    current = start_node
    visited = {start_node}
    
    while True:
        # Check if current node has outgoing edges
        if current not in graph.out_edges:
            break
        
        # Find available (unused) successors
        available_successors = [
            (target, weight)
            for target, weight in graph.out_edges[current].items()
            if (current, target) not in used_edges
            and target not in visited
        ]
        
        # If not exactly one available successor, stop
        if len(available_successors) != 1:
            break
        
        next_node, weight = available_successors[0]
        
        if next_node in visited:  # ← DODAJ (double-check)
            break

        # Add to path
        path.append(next_node)
        edges_used.add((current, next_node))
        visited.add(next_node)
        
        # Check if next_node is a merge point (multiple predecessors)
        next_in_deg = graph.in_degree.get(next_node, 0)
        if next_in_deg > 1:
            break
        
        # Check if next_node is a branch point (multiple successors)
        next_out_deg = graph.out_degree.get(next_node, 0)
        if next_out_deg > 1:
            break
        
        # Continue walking
        current = next_node
    
    return path, edges_used


def path_to_sequence(path: list[str], node_len: int) -> str:
    """
    Konwertuje ścieżkę węzłów (k-merów) w gotową sekwencję DNA
    """
    if not path:
        return ""
    
    if len(path) == 1:
        return path[0]
    
    sequence = path[0]
    
    # dla kazdego kolejnego wezła w ściezce dodaj tylko jego ostatni znak
    for node in path[1:]:
        sequence += node[-1]
    
    return sequence


def extract_contigs(
    graph: DeBruijnGraph,
    k: int,
    min_contig_len: int = 100
) -> list[str]:
    """
    Ekstrakcja kontigów z wyczyszczonego grafu de Bruijna.
    Skupia się na wyznaczeniu ścieek jednoznacznych 
    """
    contigs = []
    used_edges = set()
    node_len = k - 1
    
    # identyfikacja węzłów startowych
    start_nodes = find_start_nodes(graph)
    print(f"Start nodes: {len(start_nodes)}")
    
    for start in start_nodes:
           
        # jeśli krawędź nie została jeszcze użyta w żadnym kontigu
        if start in graph.out_edges:
            for target in list(graph.out_edges[start].keys()):
                edge = (start, target)
                
                if edge not in used_edges:
                    # Idź wzdłuż unitigu
                    path, edges = walk_unitig(graph, start, used_edges)
                    
                    # Zaznaczenie krawędzi jako wykorzystanej
                    used_edges.update(edges)
                    
                    # Zamień na sekw nukleotydową
                    sequence = path_to_sequence(path, node_len)
                    
                    if len(sequence) >=  min_contig_len:
                        contigs.append(sequence)
    
    if hasattr(graph, 'nodes'):
        all_nodes = list(graph.nodes)
    else:
        all_nodes = list(set(graph.in_degree) | set(graph.out_degree))

    # Phase 2: Radzenie z cyklami i nieuzytymi nodes
    for node in all_nodes:
        if node in graph.out_edges:
            for target in list(graph.out_edges[node].keys()):
                edge = (node, target)
                
                if edge not in used_edges:

                    path, edges = walk_unitig(graph, node, used_edges)
                    used_edges.update(edges)                    
                    sequence = path_to_sequence(path, node_len)
                    
                    if len(sequence) >=  min_contig_len:
                        contigs.append(sequence)
    
    for cont in contigs:
        print(len(cont))

    unique_contigs = list(set(contigs))
    
    return unique_contigs


def contig_stats(contigs: list[str]) -> dict[str, any]:
    """Statystyki wyznaczonych kontigów"""
    if not contigs:
        return {
            "num_contigs": 0,
            "total_length": 0,
            "longest": 0,
            "shortest": 0,
            "mean_length": 0,
            "n50": 0
        }
    
    lengths = [len(c) for c in contigs]
    total = sum(lengths)
    
    # N50: sortuj malejąco i znajdź długość, przy której suma >= 50% total
    sorted_lengths = sorted(lengths, reverse=True)
    cumsum = 0
    n50 = 0
    for length in sorted_lengths:
        cumsum += length
        if cumsum >= total / 2:
            n50 = length
            break

    stats = {
        "num_contigs": len(contigs),
        "total_length": total,
        "longest": max(lengths),
        "shortest": min(lengths),
        "mean_length": total / len(contigs),
        "n50": n50
    }
    
    return stats
