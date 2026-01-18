#!/usr/bin/env python3
"""
traversal.py - ulepszona wersja
"""
from dbg import DeBruijnGraph
from typing import List, Set, Tuple


def find_start_nodes(graph: DeBruijnGraph) -> List[str]:
    """Węzły startowe: nie są częścią prostej ścieżki."""
    starts = []
    for node in graph.nodes:
        in_deg = graph.in_degree.get(node, 0)
        out_deg = graph.out_degree.get(node, 0)
        
        # Start jeśli: źródło, rozgałęzienie, lub punkt scalenia
        if in_deg != 1 or out_deg != 1:
            starts.append(node)
    
    return starts


def get_predecessors(graph: DeBruijnGraph, node: str) -> List[str]:
    """Znajdź wszystkie węzły wchodzące do node."""
    preds = []
    for u, targets in graph.out_edges.items():
        if node in targets:
            preds.append(u)
    return preds


def walk_forward(
    graph: DeBruijnGraph,
    start: str,
    used_edges: Set[Tuple[str, str]]
) -> Tuple[List[str], Set[Tuple[str, str]]]:
    """Idź do przodu po jednoznacznej ścieżce."""
    path = [start]
    edges = set()
    current = start
    
    while True:
        if current not in graph.out_edges:
            break
        
        # Dostępne krawędzie (nieużyte)
        available = [
            (target, graph.out_edges[current][target])
            for target in graph.out_edges[current]
            if (current, target) not in used_edges
            and (current, target) not in edges
        ]
        
        if len(available) != 1:
            break
        
        next_node, weight = available[0]
        
        # Sprawdź czy next_node nie jest punktem scalenia
        if graph.in_degree.get(next_node, 0) > 1 and len(path) > 1:
            # Dodaj węzeł ale zakończ
            path.append(next_node)
            edges.add((current, next_node))
            break
        
        path.append(next_node)
        edges.add((current, next_node))
        
        # Sprawdź czy next_node nie jest rozgałęzieniem
        if graph.out_degree.get(next_node, 0) > 1:
            break
        
        current = next_node
    
    return path, edges


def walk_backward(
    graph: DeBruijnGraph,
    start: str,
    used_edges: Set[Tuple[str, str]]
) -> Tuple[List[str], Set[Tuple[str, str]]]:
    """Idź do tyłu po jednoznacznej ścieżce."""
    path = []
    edges = set()
    current = start
    
    while True:
        preds = get_predecessors(graph, current)
        
        # Filtruj nieużyte
        available = [
            p for p in preds
            if (p, current) not in used_edges
            and (p, current) not in edges
        ]
        
        if len(available) != 1:
            break
        
        pred = available[0]
        
        # Sprawdź czy pred nie jest rozgałęzieniem
        if graph.out_degree.get(pred, 0) > 1 and path:
            path.append(pred)
            edges.add((pred, current))
            break
        
        path.append(pred)
        edges.add((pred, current))
        
        # Sprawdź czy pred nie jest punktem scalenia
        if graph.in_degree.get(pred, 0) > 1:
            break
        
        current = pred
    
    return path, edges


def walk_unitig_bidirectional(
    graph: DeBruijnGraph,
    start: str,
    used_edges: Set[Tuple[str, str]]
) -> Tuple[List[str], Set[Tuple[str, str]]]:
    """
    Idź w OBU kierunkach od węzła startowego.
    Zwraca pełną ścieżkę unitigu.
    """
    # Do przodu
    forward_path, forward_edges = walk_forward(graph, start, used_edges)
    
    # Do tyłu
    backward_path, backward_edges = walk_backward(graph, start, used_edges)
    
    # Połącz: [backward odwrócony] + [forward bez pierwszego elementu]
    full_path = backward_path[::-1] + forward_path
    all_edges = forward_edges | backward_edges
    
    return full_path, all_edges


def path_to_sequence(path: List[str]) -> str:
    """Konwertuj ścieżkę węzłów na sekwencję DNA."""
    if not path:
        return ""
    if len(path) == 1:
        return path[0]
    
    sequence = path[0]
    for node in path[1:]:
        sequence += node[-1]
    
    return sequence


def extract_contigs(
    graph: DeBruijnGraph,
    k: int,
    min_contig_len: int = 100
) -> List[str]:
    """
    Ekstrakcja kontigów - ulepszona wersja.
    """
    contigs = []
    used_edges: Set[Tuple[str, str]] = set()
    
    # Faza 1: Zacznij od węzłów specjalnych (rozgałęzienia, źródła, ujścia)
    start_nodes = find_start_nodes(graph)
    
    for start in start_nodes:
        if start not in graph.out_edges:
            continue
        
        for target in graph.out_edges[start]:
            if (start, target) in used_edges:
                continue
            
            # Dwukierunkowy traversal
            path, edges = walk_unitig_bidirectional(graph, start, used_edges)
            used_edges.update(edges)
            
            seq = path_to_sequence(path)
            if len(seq) >= min_contig_len:
                contigs.append(seq)
    
    # Faza 2: Złap pozostałe komponenty (np. cykle)
    for node in graph.nodes:
        if node not in graph.out_edges:
            continue
        
        for target in graph.out_edges[node]:
            if (node, target) in used_edges:
                continue
            
            path, edges = walk_unitig_bidirectional(graph, node, used_edges)
            used_edges.update(edges)
            
            seq = path_to_sequence(path)
            if len(seq) >= min_contig_len:
                contigs.append(seq)
    
    # Usuń duplikaty
    return list(set(contigs))


def contig_stats(contigs: List[str]) -> dict:
    """Statystyki kontigów."""
    if not contigs:
        return {
            "num_contigs": 0, "total_length": 0, "longest": 0,
            "shortest": 0, "mean_length": 0, "n50": 0
        }
    
    lengths = sorted([len(c) for c in contigs], reverse=True)
    total = sum(lengths)
    
    # N50
    cumsum = 0
    n50 = 0
    for length in lengths:
        cumsum += length
        if cumsum >= total / 2:
            n50 = length
            break
    
    return {
        "num_contigs": len(contigs),
        "total_length": total,
        "longest": lengths[0],
        "shortest": lengths[-1],
        "mean_length": round(total / len(contigs), 1),
        "n50": n50
    }