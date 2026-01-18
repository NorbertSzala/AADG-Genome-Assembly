#!/usr/bin/env python3
"""
velvet_optimizations.py - Optymalizacje w stylu Velvet dla grafu de Bruijna.

Zawiera:
1. Coverage cutoff - automatyczne usuwanie błędów (metoda valley detection)
2. Tour Bus - rozwiązywanie bubble/powtórzeń na podstawie coverage
3. Rock Band - łączenie ścieżek o podobnym coverage (opcjonalne)
"""

from collections import Counter
from typing import List, Tuple, Optional
from dbg import DeBruijnGraph

# Import helperów z cleaning.py
from cleaning import (
    get_predecessors,
    get_successors, 
    get_edge_weight,
    remove_node,
    remove_nodes,
    remove_edge,
    node_exists
)


###############################
####### COVERAGE ANALYSIS #####
###############################

def get_all_edge_weights(graph: DeBruijnGraph) -> List[int]:
    """Zbierz wagi wszystkich krawędzi w grafie."""
    weights = []
    for u, targets in graph.out_edges.items():
        for v, weight in targets.items():
            weights.append(weight)
    return weights


def build_coverage_histogram(graph: DeBruijnGraph) -> Counter:
    """Zbuduj histogram pokrycia krawędzi."""
    weights = get_all_edge_weights(graph)
    return Counter(weights)


def estimate_coverage_cutoff_valley(graph: DeBruijnGraph, max_search: int = 15) -> int:
    """
    Oszacuj próg pokrycia metodą "valley detection".
    
    Teoria: Histogram k-merów ma dwa piki:
    - Niski pik (coverage 1-3) = błędy sekwencjonowania
    - Wysoki pik (~średnie coverage) = prawdziwe k-mery
    - Między nimi jest "dolina" = optymalny cutoff
    
    Args:
        graph: Graf de Bruijna
        max_search: Maksymalny zakres szukania doliny
    
    Returns:
        Optymalny próg coverage
    """
    histogram = build_coverage_histogram(graph)
    
    if not histogram:
        return 2
    
    # Szukaj lokalnego minimum (doliny) w zakresie 1-max_search
    counts = sorted([c for c in histogram.keys() if 1 <= c <= max_search])
    
    if len(counts) < 3:
        return 2
    
    # Metoda 1: Znajdź dolinę (h[i-1] > h[i] < h[i+1])
    for i in range(1, len(counts) - 1):
        c = counts[i]
        prev_c, next_c = counts[i-1], counts[i+1]
        
        if histogram[c] < histogram[prev_c] and histogram[c] < histogram[next_c]:
            return c + 1  # Cutoff tuż za doliną
    
    # Metoda 2: Fallback - użyj stosunku singletonów
    total = sum(histogram.values())
    singletons = histogram.get(1, 0)
    singleton_ratio = singletons / total if total > 0 else 0
    
    if singleton_ratio > 0.6:
        return 4  # Dużo błędów
    elif singleton_ratio > 0.4:
        return 3
    else:
        return 2


###############################
####### COVERAGE CUTOFF #######
###############################

def apply_coverage_cutoff(graph: DeBruijnGraph, min_coverage: int) -> DeBruijnGraph:
    """
    Usuń krawędzie o pokryciu poniżej progu.
    
    Args:
        graph: Graf de Bruijna
        min_coverage: Minimalny próg coverage
    
    Returns:
        Zmodyfikowany graf
    """
    # Zbierz krawędzie do usunięcia
    edges_to_remove: List[Tuple[str, str]] = []
    
    for u, targets in list(graph.out_edges.items()):
        for v, weight in list(targets.items()):
            if weight < min_coverage:
                edges_to_remove.append((u, v))
    
    # Usuń krawędzie
    for u, v in edges_to_remove:
        remove_edge(graph, u, v)
    
    # Usuń izolowane węzły (in=0 i out=0)
    isolated = [
        node for node in list(graph.nodes)
        if graph.in_degree.get(node, 0) == 0 
        and graph.out_degree.get(node, 0) == 0
    ]
    
    remove_nodes(graph, set(isolated))
    
    return graph


###############################
######### TOUR BUS ############
###############################

def walk_path_with_coverage(
    graph: DeBruijnGraph, 
    start: str, 
    max_len: int,
    stop_at_branch: bool = True
) -> Tuple[List[str], int]:
    """
    Idź ścieżką od start, zbierając coverage.
    
    Args:
        graph: Graf
        start: Węzeł startowy
        max_len: Maksymalna długość ścieżki
        stop_at_branch: Zatrzymaj się na rozgałęzieniach
    
    Returns:
        (path, total_coverage)
    """
    path = [start]
    total_cov = 0
    current = start
    
    for step in range(max_len):
        # Sprawdź czy to prosta ścieżka (po pierwszym kroku)
        if step > 0 and stop_at_branch:
            out_deg = graph.out_degree.get(current, 0)
            in_deg = graph.in_degree.get(current, 0)
            
            if out_deg != 1 or in_deg > 1:
                break
        
        successors = get_successors(graph, current)
        if not successors:
            break
        
        # Wybierz następny węzeł (jeśli wiele - weź pierwszy)
        nxt = successors[0]
        
        # Dodaj coverage tej krawędzi
        edge_weight = get_edge_weight(graph, current, nxt)
        total_cov += edge_weight
        
        path.append(nxt)
        current = nxt
    
    return path, total_cov


def find_bubble_paths(
    graph: DeBruijnGraph,
    branch_node: str,
    max_bubble_len: int
) -> List[Tuple[List[str], int]]:
    """
    Znajdź wszystkie ścieżki wychodzące z węzła rozgałęzienia.
    
    Returns:
        Lista (path, coverage) dla każdej gałęzi
    """
    successors = get_successors(graph, branch_node)
    
    if len(successors) < 2:
        return []
    
    paths = []
    for succ in successors:
        path, cov = walk_path_with_coverage(graph, succ, max_bubble_len)
        # Dodaj coverage krawędzi branch_node -> succ
        cov += get_edge_weight(graph, branch_node, succ)
        paths.append((path, cov))
    
    return paths


def paths_merge_at_same_node(paths: List[Tuple[List[str], int]]) -> Optional[str]:
    """
    Sprawdź czy wszystkie ścieżki kończą się w tym samym węźle.
    
    Returns:
        Węzeł końcowy jeśli się łączą, None w przeciwnym razie
    """
    if len(paths) < 2:
        return None
    
    end_nodes = [path[0][-1] for path in paths if path[0]]
    
    if len(set(end_nodes)) == 1:
        return end_nodes[0]
    
    return None


def remove_weakest_bubble_path(
    graph: DeBruijnGraph,
    paths: List[Tuple[List[str], int]],
    branch_node: str
) -> bool:
    """
    Usuń najsłabszą ścieżkę bubble (najniższy coverage).
    
    Returns:
        True jeśli usunięto, False w przeciwnym razie
    """
    if len(paths) < 2:
        return False
    
    # Znajdź najsłabszą ścieżkę
    min_cov = float('inf')
    weakest_idx = 0
    
    for i, (path, cov) in enumerate(paths):
        if cov < min_cov:
            min_cov = cov
            weakest_idx = i
    
    weakest_path, _ = paths[weakest_idx]
    
    if not weakest_path:
        return False
    
    # Usuń krawędź branch_node -> pierwszy węzeł słabszej ścieżki
    first_node = weakest_path[0]
    remove_edge(graph, branch_node, first_node)
    
    # Usuń węzły słabszej ścieżki (oprócz ostatniego - punkt scalenia)
    nodes_to_remove = set(weakest_path[:-1])
    
    # Nie usuwaj węzłów które są częścią innych ścieżek
    for i, (path, _) in enumerate(paths):
        if i != weakest_idx:
            nodes_to_remove -= set(path)
    
    # Usuń tylko izolowane węzły ze słabszej ścieżki
    for node in nodes_to_remove:
        if node_exists(graph, node):
            in_deg = graph.in_degree.get(node, 0)
            out_deg = graph.out_degree.get(node, 0)
            
            # Usuń tylko jeśli stał się izolowany lub dead-end
            if in_deg == 0 or out_deg == 0:
                remove_node(graph, node)
    
    return True


def tour_bus(
    graph: DeBruijnGraph,
    max_bubble_len: int = 15,
    max_iterations: int = 10
) -> DeBruijnGraph:
    """
    Tour Bus algorithm - rozwiązuj bubble/powtórzenia na podstawie coverage.
    
    Gdy ścieżki rozchodzą się i łączą ponownie, usuń tę o niższym coverage.
    
    Args:
        graph: Graf de Bruijna
        max_bubble_len: Maksymalna długość ścieżki bubble
        max_iterations: Maksymalna liczba iteracji
    
    Returns:
        Zmodyfikowany graf
    """
    total_resolved = 0
    
    for iteration in range(max_iterations):
        resolved_this_round = 0
        
        # Znajdź wszystkie węzły rozgałęzienia
        branch_nodes = [
            node for node in list(graph.nodes)
            if graph.out_degree.get(node, 0) >= 2
        ]
        
        for branch in branch_nodes:
            if not node_exists(graph, branch):
                continue
            
            # Sprawdź czy nadal jest rozgałęzieniem
            if graph.out_degree.get(branch, 0) < 2:
                continue
            
            # Znajdź ścieżki
            paths = find_bubble_paths(graph, branch, max_bubble_len)
            
            if len(paths) < 2:
                continue
            
            # Sprawdź czy łączą się w tym samym węźle
            merge_node = paths_merge_at_same_node(paths)
            
            if merge_node:
                # To jest bubble - usuń słabszą ścieżkę
                if remove_weakest_bubble_path(graph, paths, branch):
                    resolved_this_round += 1
        
        total_resolved += resolved_this_round
        
        if resolved_this_round == 0:
            break
    
    return graph

###############################
###### MAIN OPTIMIZATION ######
###############################

def velvet_style_optimization(
    graph: DeBruijnGraph,
    auto_cutoff: bool = True,
    cutoff_value: int = None,
    use_tourbus: bool = True,
    tourbus_max_len: int = 15
) -> DeBruijnGraph:
    """
    Zastosuj optymalizacje w stylu Velvet.
    
    Pipeline:
    1. Auto-detekcja progu pokrycia (valley detection)
    2. Usunięcie krawędzi o niskim pokryciu
    3. Tour Bus dla rozwiązywania bubble
    
    Args:
        graph: Wejściowy graf
        auto_cutoff: Automatyczne oszacowanie progu pokrycia
        cutoff_value: Ręczny próg (jeśli auto_cutoff=False)
        use_tourbus: Użyj algorytmu Tour Bus
        tourbus_max_len: Maks długość bubble dla Tour Bus
        verbose: Wypisuj informacje diagnostyczne
    
    Returns:
        Zoptymalizowany graf
    """
    initial_nodes = len(graph.nodes)
    
    # 1. Coverage cutoff
    if auto_cutoff:
        cutoff = estimate_coverage_cutoff_valley(graph)
    elif cutoff_value is not None:
        cutoff = cutoff_value
    else:
        cutoff = 2
    
    graph = apply_coverage_cutoff(graph, min_coverage=cutoff)
    
    # 2. Tour Bus
    if use_tourbus:
        graph = tour_bus(graph, max_bubble_len=tourbus_max_len)
        
    return graph