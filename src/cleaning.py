#!/usr/bin/env python3
"""
cleaning.py - Graph cleaning utilities for De Bruijn graph.

Implements:
- Helper functions for node/edge operations
- Island removal (small disconnected components)
- Tip removal (short dead-end paths)
- Bubble popping (simple bubbles)
"""

from collections import deque
from typing import List, Set, Tuple
from dbg import DeBruijnGraph


###############################
######## HELPER FUNCTIONS #####
###############################

def get_predecessors(graph: DeBruijnGraph, node: str) -> List[str]:
    """Znajdź wszystkie węzły wchodzące do node."""
    return [u for u, targets in graph.out_edges.items() if node in targets]


def get_successors(graph: DeBruijnGraph, node: str) -> List[str]:
    """Znajdź wszystkie węzły wychodzące z node."""
    if node in graph.out_edges:
        return list(graph.out_edges[node].keys())
    return []


def get_neighbors_undirected(graph: DeBruijnGraph, node: str) -> Set[str]:
    """Pobierz wszystkich sąsiadów (ignorując kierunek krawędzi)."""
    neighbors = set()
    
    # Wychodzące
    if node in graph.out_edges:
        neighbors.update(graph.out_edges[node].keys())
    
    # Wchodzące
    neighbors.update(get_predecessors(graph, node))
    
    return neighbors


def get_edge_weight(graph: DeBruijnGraph, u: str, v: str) -> int:
    """Pobierz wagę krawędzi u -> v."""
    if u in graph.out_edges and v in graph.out_edges[u]:
        return graph.out_edges[u][v]
    return 0


def remove_edge(graph: DeBruijnGraph, u: str, v: str) -> None:
    """Usuń krawędź u -> v z grafu."""
    if u in graph.out_edges and v in graph.out_edges[u]:
        del graph.out_edges[u][v]
        graph.out_degree[u] -= 1
        graph.in_degree[v] -= 1
        
        # Wyczyść pusty słownik
        if not graph.out_edges[u]:
            del graph.out_edges[u]


def remove_node(graph: DeBruijnGraph, node: str) -> None:
    """Usuń węzeł i wszystkie jego krawędzie z grafu."""
    # Usuń wychodzące krawędzie
    if node in graph.out_edges:
        for v in list(graph.out_edges[node].keys()):
            graph.in_degree[v] -= 1
        del graph.out_edges[node]
    
    # Usuń wchodzące krawędzie
    for u in list(graph.out_edges.keys()):
        if node in graph.out_edges[u]:
            del graph.out_edges[u][node]
            graph.out_degree[u] -= 1
            
            # Wyczyść pusty słownik
            if not graph.out_edges[u]:
                del graph.out_edges[u]
    
    # Usuń z bookkeepingu
    graph.in_degree.pop(node, None)
    graph.out_degree.pop(node, None)
    graph.nodes.discard(node)


def remove_nodes(graph: DeBruijnGraph, nodes: Set[str]) -> None:
    """Usuń wiele węzłów naraz (bardziej efektywne)."""
    for node in nodes:
        remove_node(graph, node)


def node_exists(graph: DeBruijnGraph, node: str) -> bool:
    """Sprawdź czy węzeł istnieje w grafie."""
    return node in graph.nodes


###############################
######## ISLAND REMOVAL #######
###############################

def find_connected_component(graph: DeBruijnGraph, start: str, visited: Set[str]) -> Set[str]:
    """BFS aby znaleźć komponent połączony z start."""
    component = set()
    queue = deque([start])
    
    while queue:
        current = queue.popleft()
        if current in visited:
            continue
        
        visited.add(current)
        component.add(current)
        
        for neighbor in get_neighbors_undirected(graph, current):
            if neighbor not in visited:
                queue.append(neighbor)
    
    return component


def remove_islands(graph: DeBruijnGraph, min_component_size: int) -> DeBruijnGraph:
    """
    Usuń małe rozłączone komponenty (wyspy).
    
    Args:
        graph: Graf de Bruijna
        min_component_size: Minimalna liczba węzłów do zachowania komponentu
    
    Returns:
        Zmodyfikowany graf
    """
    visited: Set[str] = set()
    nodes_to_remove: Set[str] = set()
    
    for node in list(graph.nodes):
        if node in visited:
            continue
        
        component = find_connected_component(graph, node, visited)
        
        if len(component) < min_component_size:
            nodes_to_remove.update(component)
    
    remove_nodes(graph, nodes_to_remove)
    
    return graph


###############################
######### TIP REMOVAL #########
###############################

def trace_tip_backward(graph: DeBruijnGraph, start: str, max_len: int) -> List[str]:
    """
    Śledź tip wstecz od dead-end do rozgałęzienia.
    Zwraca ścieżkę lub None jeśli nie jest tipem.
    """
    path = [start]
    current = start
    
    while len(path) <= max_len:
        preds = get_predecessors(graph, current)
        
        if len(preds) != 1:
            break
        
        pred = preds[0]
        path.append(pred)
        
        # Dotarliśmy do rozgałęzienia - to jest tip
        if graph.out_degree.get(pred, 0) > 1:
            return path
        
        current = pred
    
    # Ścieżka za długa lub nie kończy się rozgałęzieniem
    return None


def trace_tip_forward(graph: DeBruijnGraph, start: str, max_len: int) -> List[str]:
    """
    Śledź tip do przodu od źródła do punktu scalenia.
    Zwraca ścieżkę lub None jeśli nie jest tipem.
    """
    path = [start]
    current = start
    
    while len(path) <= max_len:
        successors = get_successors(graph, current)
        
        if len(successors) != 1:
            break
        
        nxt = successors[0]
        path.append(nxt)
        
        # Dotarliśmy do punktu scalenia - to jest tip
        if graph.in_degree.get(nxt, 0) > 1:
            return path
        
        current = nxt
    
    return None


def calculate_path_coverage(graph: DeBruijnGraph, path: List[str]) -> int:
    """Oblicz sumę wag krawędzi na ścieżce."""
    total = 0
    for i in range(len(path) - 1):
        total += get_edge_weight(graph, path[i], path[i + 1])
    return total


def should_remove_tip(
    graph: DeBruijnGraph, 
    path: List[str], 
    min_coverage_ratio: float = 0.1
) -> bool:
    """
    Zdecyduj czy usunąć tip na podstawie coverage.
    Tip jest usuwany jeśli ma znacznie niższy coverage niż alternatywa.
    """
    if len(path) < 2:
        return True  # Bardzo krótkie tipy zawsze usuwamy
    
    tip_coverage = calculate_path_coverage(graph, path)
    
    # Znajdź węzeł rozgałęzienia (ostatni w ścieżce backward lub forward)
    branch_node = path[-1]
    
    # Znajdź coverage alternatywnej ścieżki
    alt_coverage = 0
    tip_neighbor = path[-2]
    
    for target in get_successors(graph, branch_node):
        if target != tip_neighbor:
            weight = get_edge_weight(graph, branch_node, target)
            alt_coverage = max(alt_coverage, weight)
    
    for pred in get_predecessors(graph, branch_node):
        if pred != tip_neighbor:
            weight = get_edge_weight(graph, pred, branch_node)
            alt_coverage = max(alt_coverage, weight)
    
    # Usuń jeśli coverage tipu jest znacznie niższy
    if alt_coverage > 0:
        return tip_coverage / alt_coverage < min_coverage_ratio
    
    return True  # Brak alternatywy - usuń tip


def remove_tips(
    graph: DeBruijnGraph,
    tip_max_len: int,
    min_coverage_ratio: float = 0.1,
    max_iterations: int = 20
) -> DeBruijnGraph:
    """
    Usuń krótkie dead-end ścieżki (tipy).
    
    Args:
        graph: Graf de Bruijna
        tip_max_len: Maksymalna długość tipu (w węzłach)
        min_coverage_ratio: Minimalny stosunek coverage tip/alternatywa
        max_iterations: Maksymalna liczba iteracji
    
    Returns:
        Zmodyfikowany graf
    """
    for iteration in range(max_iterations):
        nodes_to_remove: Set[str] = set()
        
        for node in list(graph.nodes):
            if node in nodes_to_remove:
                continue
            
            if not node_exists(graph, node):
                continue
            
            in_deg = graph.in_degree.get(node, 0)
            out_deg = graph.out_degree.get(node, 0)
            
            path = None
            
            # Tip na końcu (dead-end): out_degree == 0
            if out_deg == 0 and in_deg > 0:
                path = trace_tip_backward(graph, node, tip_max_len)
            
            # Tip na początku (źródło): in_degree == 0
            elif in_deg == 0 and out_deg > 0:
                path = trace_tip_forward(graph, node, tip_max_len)
            
            if path and should_remove_tip(graph, path, min_coverage_ratio):
                # Usuń wszystkie węzły oprócz węzła rozgałęzienia
                nodes_to_remove.update(path[:-1])
        
        if not nodes_to_remove:
            break  # Brak zmian - koniec
        
        remove_nodes(graph, nodes_to_remove)
    
    return graph


###############################
######## BUBBLE POPPING #######
###############################

def walk_simple_path(
    graph: DeBruijnGraph, 
    start: str, 
    max_len: int
) -> Tuple[str, List[str], int]:
    """
    Idź prostą ścieżką (in=1, out=1) od start.
    
    Returns:
        (end_node, path, total_weight)
    """
    path = [start]
    weight_sum = 0
    current = start
    
    for _ in range(max_len):
        out_deg = graph.out_degree.get(current, 0)
        in_deg = graph.in_degree.get(current, 0)
        
        # Zatrzymaj się jeśli nie jest to prosta ścieżka
        if out_deg != 1 or in_deg != 1:
            break
        
        successors = get_successors(graph, current)
        if not successors:
            break
        
        nxt = successors[0]
        weight_sum += get_edge_weight(graph, current, nxt)
        path.append(nxt)
        current = nxt
    
    return current, path, weight_sum


def pop_bubbles_simple(
    graph: DeBruijnGraph,
    max_bubble_len: int
) -> DeBruijnGraph:
    """
    Usuń proste bubble (rozgałęzienie -> dwie ścieżki -> scalenie).
    
    Strategia: Usuń słabszą ścieżkę (niższy coverage).
    
    Args:
        graph: Graf de Bruijna
        max_bubble_len: Maksymalna długość ścieżki bubble
    
    Returns:
        Zmodyfikowany graf
    """
    for u in list(graph.out_edges.keys()):
        if not node_exists(graph, u):
            continue
        
        targets = get_successors(graph, u)
        
        # Szukamy rozgałęzienia na dokładnie 2 ścieżki
        if len(targets) != 2:
            continue
        
        v1, v2 = targets
        
        # Idź obiema ścieżkami
        end1, path1, weight1 = walk_simple_path(graph, v1, max_bubble_len)
        end2, path2, weight2 = walk_simple_path(graph, v2, max_bubble_len)
        
        # Sprawdź czy łączą się w tym samym węźle
        if end1 != end2:
            continue
        
        # To jest bubble! Usuń słabszą ścieżkę
        if weight1 < weight2:
            weaker_path = path1
        else:
            weaker_path = path2
        
        # Usuń węzły słabszej ścieżki (oprócz końcowego)
        for node in weaker_path[:-1]:
            if node_exists(graph, node):
                remove_node(graph, node)
    
    return graph


###############################
######## COMBINED CLEANING ####
###############################

def clean_graph(
    graph: DeBruijnGraph,
    min_component_size: int = 10,
    tip_max_len: int = 3,
    pop_bubbles: bool = False,
    max_bubble_len: int = 5,
    min_coverage_ratio: float = 0.1
) -> DeBruijnGraph:
    """
    Pełne czyszczenie grafu w jednej funkcji.
    
    Args:
        graph: Graf de Bruijna
        min_component_size: Min węzłów w komponencie
        tip_max_len: Max długość tipu
        pop_bubbles: Czy usuwać bubble
        max_bubble_len: Max długość bubble
        min_coverage_ratio: Min stosunek coverage dla tipów
    
    Returns:
        Wyczyszczony graf
    """
    graph = remove_islands(graph, min_component_size)
    graph = remove_tips(graph, tip_max_len, min_coverage_ratio)
    
    if pop_bubbles:
        graph = pop_bubbles_simple(graph, max_bubble_len)
    
    # Drugie przejście tip removal po bubble popping
    if pop_bubbles:
        graph = remove_tips(graph, tip_max_len, min_coverage_ratio)
    
    return graph
