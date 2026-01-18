#!/usr/bin/env python3
"""
Optymalizacje w stylu Velvet dla grafu de Bruijna.

Zawiera:
1. Coverage cutoff - automatyczne usuwanie błędów na podstawie pokrycia k-merów
2. Tour Bus - rozwiązywanie prostych powtórzeń
3. Sugestie optymalnego k
"""

from dbg import DeBruijnGraph


def estimate_coverage_cutoff(graph: DeBruijnGraph, percentile: float = 0.2) -> int:
    """
    Oszacuj automatycznie próg pokrycia na podstawie rozkładu wag krawędzi.
    
    Teoria: błędne k-mery mają niskie pokrycie (1-2×), prawdziwe ~pokrycie genomu.
    """
    coverages = []
    for node, targets in graph.out_edges.items():
        for target, weight in targets.items():
            coverages.append(weight)
    
    if not coverages:
        return 2
    
    coverages.sort()
    cutoff_idx = int(len(coverages) * percentile)
    cutoff = coverages[cutoff_idx]
    
    return max(2, cutoff)


def apply_coverage_cutoff(graph: DeBruijnGraph, min_coverage: int) -> DeBruijnGraph:
    """Usuń krawędzie o pokryciu poniżej progu (prawdopodobnie błędy)."""
    edges_to_remove = []
    
    # Znajdź słabe krawędzie
    for u, targets in graph.out_edges.items():
        for v, weight in targets.items():
            if weight < min_coverage:
                edges_to_remove.append((u, v))
    
    # Usuń krawędzie
    for u, v in edges_to_remove:
        if u in graph.out_edges and v in graph.out_edges[u]:
            del graph.out_edges[u][v]
            
            if u in graph.out_degree:
                graph.out_degree[u] -= 1
            if v in graph.in_degree:
                graph.in_degree[v] -= 1
    
    # Usuń izolowane węzły
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


def _walk_with_coverage(graph: DeBruijnGraph, start_node: str, max_len: int) -> tuple[list[str], int]:
    """
    Idź ścieżką od start_node, zwróć path i całkowite pokrycie.
    Zatrzymaj się na rozgałęzieniach lub po max_len krokach.
    """
    path = [start_node]
    total_cov = 0
    curr = start_node
    
    for step in range(max_len):
        # Sprawdź stopnie tylko PO pierwszym węźle (fix dla tour bus)
        if step > 0:
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


def tour_bus_optimized(graph: DeBruijnGraph, max_repeat_length: int = 15, max_iterations: int = 5) -> DeBruijnGraph:
    """
    Tour Bus - rozwiązuj krótkie powtórzenia używając pokrycia.
    
    Gdy dwie ścieżki rozchodzą się i łączą ponownie - wybierz ścieżkę o wyższym pokryciu.
    Wielokrotne iteracje dla trudnych przypadków.
    """
    total_resolved = 0
    
    for iteration in range(max_iterations):
        resolved = 0
        
        for start in list(graph.out_edges.keys()):
            if start not in graph.out_edges:
                continue
            
            targets = list(graph.out_edges[start].keys())
            if len(targets) < 2:
                continue
            
            # Sprawdź każdą parę gałęzi
            for i in range(len(targets)):
                for j in range(i + 1, len(targets)):
                    v1, v2 = targets[i], targets[j]
                    
                    # Sprawdź czy nadal istnieją
                    if v1 not in graph.out_edges.get(start, {}) or \
                       v2 not in graph.out_edges.get(start, {}):
                        continue
                    
                    # Idź obiema ścieżkami
                    path1, cov1 = _walk_with_coverage(graph, v1, max_repeat_length)
                    path2, cov2 = _walk_with_coverage(graph, v2, max_repeat_length)
                    
                    # Jeśli łączą się w tym samym węźle - usuń słabszą
                    if path1 and path2 and path1[-1] == path2[-1]:
                        to_remove_start = v1 if cov1 < cov2 else v2
                        
                        if to_remove_start in graph.out_edges[start]:
                            del graph.out_edges[start][to_remove_start]
                            
                            if start in graph.out_degree:
                                graph.out_degree[start] -= 1
                            if to_remove_start in graph.in_degree:
                                graph.in_degree[to_remove_start] -= 1
                            
                            resolved += 1
        
        total_resolved += resolved
        
        # Jeśli nic nie zmieniono - zakończ
        if resolved == 0:
            break
    
    return graph

def velvet_style_optimization(
    graph: DeBruijnGraph,
    auto_cutoff: bool = True,
    use_tourbus: bool = True) -> DeBruijnGraph:
    """
    Zastosuj optymalizacje w stylu Velvet.
    
    Pipeline:
    1. Auto-detekcja progu pokrycia
    2. Usunięcie krawędzi o niskim pokryciu
    3. Tour Bus dla krótkich powtórzeń
    
    Args:
        graph: Wejściowy graf
        auto_cutoff: Automatyczne oszacowanie i zastosowanie progu pokrycia
        use_tourbus: Użyj algorytmu Tour Bus
    
    Returns:
        Zoptymalizowany graf
    """
    # Coverage cutoff
    if auto_cutoff:
        cutoff = estimate_coverage_cutoff(graph, percentile=0.2)
        graph = apply_coverage_cutoff(graph, min_coverage=cutoff)
    
    # Tour Bus
    if use_tourbus:
        graph = tour_bus_optimized(graph, max_repeat_length=10)
    
    return graph