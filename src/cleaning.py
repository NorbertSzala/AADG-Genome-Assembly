#!/usr/bin/env python3

"""
Summary:
Graph cleaning utilities for De Bruijn graph.

Implements removal of:
- islands (small disconnected components)
- tips (short dead-end paths)
- simple bubbles (optional, conservative)
"""


###############################
########### IMPORTS  ##########
###############################
from collections import deque
from dbg import DeBruijnGraph


###############################
########## FUNCTIONS ##########
###############################


# Helper function
def _neighbors_undirected(graph: DeBruijnGraph, node: str) -> set[str]:
    neighbors = set()

    if node in graph.out_edges:
        neighbors.update(graph.out_edges[node].keys())

    # incoming neighbors
    for u, targets in graph.out_edges.items():
        if node in targets:
            neighbors.add(u)
    return neighbors


def remove_islands(
    graph: DeBruijnGraph, min_component_size_nodes: int
) -> DeBruijnGraph:
    """
    Remove small disconected components (islands) from graph to clean up.

    Generally, real genome form some large conected graph but with many 'mistakes' (small disconected subraphs).
    We find all conected components (ignoring edge direction)

    Input:
    - DeBruijne Graph
    - min_component_size_nodes: minimal number of nodes a component must have to be kept

    OUTPUT:
    - the same graph object, modified in place
    """

    visited = set()  # Keep track of nodes that were already visited during BFS/DFS
    all_nodes = set(graph.in_degree) | set(
        graph.out_degree
    )  # All nodes in graoh (union of nodes with in/out degree)

    for node in list(all_nodes):
        if node in visited:
            continue

        # Start BFS from this node to find its connected component
        queue = deque([node])
        component = set([node])
        visited.add(node)

        # BFS loop
        while queue:
            current = queue.popleft()

            # get all neighbors of the current node, ignoring direction (incoming and outgoing)
            for neigh in _neighbors_undirected(graph, current):
                if neigh not in visited:
                    visited.add(neigh)
                    component.add(neigh)
                    queue.append(neigh)

        # now, 'component' contain all nodes connected to 'node'
        # Example:  component = {"ACGT", "CGTA", "GTAC"}

        # If the component is too small, consider it as an island
        if len(component) < min_component_size_nodes:
            # remove every node in this small component
            for n in component:
                # Remove outgoing edges from n
                if n in graph.out_edges:
                    for v in list(graph.out_edges[n]):
                        graph.in_degree[v] -= 1

                    del graph.out_edges[n]

                # remove incoming edges:
                for u in list(graph.out_edges):
                    if n in graph.out_edges[u]:
                        del graph.out_edges[u][n]
                        graph.out_degree[u] -= 1

                # remove degree bookkeeping for node n
                graph.in_degree.pop(n, None)
                graph.out_degree.pop(n, None)
                
                # Usuń z graph.nodes
                if hasattr(graph, 'nodes'):
                    graph.nodes.discard(n)
    
    return graph


def remove_tips(
    graph: DeBruijnGraph,
    tip_max_len: int,
) -> DeBruijnGraph:
    """
    Remove short dead-end paths (tips).

    As a 'tip' we understand a path that ends without outgoing edges or starts at node with in_degree==0 or out_degree==0

    Main idea is to find dead-end nodes, walk along the path as long as in_degree==1 and out_degree==1. Stop when we reach a branching node or the bath becomes too long.

    If path is short enough, delete it.

    Inpout:
        - graph: DeBruijnGraph
        - tip_max_len: maximum allowed length of a tip (in nodes)

    OUTPUT:
        - cleaned graph
    """

    # Iteruje wielokrotnie aż stabilizuje sie
    changed = True
    iterations = 0
    max_iterations = 30  
    
    while changed and iterations < max_iterations:
        changed = False
        iterations += 1
        nodes_to_remove = []
        nodes = set(graph.in_degree) | set(graph.out_degree)  # all nodes in the graph

        # check each node as a potential tip start
        for node in list(nodes):
            # Skip if already marked for removal
            if node in nodes_to_remove:
                continue
            
            # A dead-end node has no incoming OR no outgoing edges
            in_deg = graph.in_degree.get(node, 0)
            out_deg = graph.out_degree.get(node, 0)
            
            if out_deg == 0 or in_deg == 0:
                path = [node]
                current = node

                # FIX 6: Walk from dead end toward potential branch
                if out_deg == 0 and in_deg > 0:
                    # Tip at end - walk backwards via incoming edges
                    while len(path) < tip_max_len:
                        # Find predecessors
                        predecessors = [u for u in graph.out_edges if current in graph.out_edges[u]]
                        
                        if len(predecessors) != 1: break
                        
                        pred = predecessors[0]                        
                        # If predecessor is branch, we found a valid tip
                        if graph.out_degree.get(pred, 0) > 1:
                            path.append(pred)
                            break
                        
                        path.append(pred)
                        current = pred
                
                elif in_deg == 0 and out_deg > 0:
                    # Tip at start - walk forward
                    while len(path) <= tip_max_len:
                        next_nodes = list(graph.out_edges.get(current, {}).keys())
                        
                        if len(next_nodes) != 1:
                            break
                        
                        next_node = next_nodes[0]
                        
                        # If next is merge point, stop
                        if graph.in_degree.get(next_node, 0) > 1:
                            path.append(next_node)
                            break
                        
                        path.append(next_node)
                        current = next_node

                # If the path is short enough, mark for removal
                if 1 <= len(path) <= tip_max_len:
                    nodes_to_remove.extend(path)
                    changed = True

        # Remove collected nodes
        for n in nodes_to_remove:
            # Skip if already removed
            if hasattr(graph, 'nodes') and n not in graph.nodes:
                continue
            
            # remove outgoing
            if n in graph.out_edges:
                for v in list(graph.out_edges[n].keys()):
                    if v in graph.in_degree:
                        graph.in_degree[v] -= 1
                del graph.out_edges[n]

            # remove incoming
            for u in list(graph.out_edges.keys()):
                if n in graph.out_edges[u]:
                    del graph.out_edges[u][n]
                    if u in graph.out_degree:
                        graph.out_degree[u] -= 1

            # Remove degree bookkeeping
            graph.in_degree.pop(n, None)
            graph.out_degree.pop(n, None)
            
            if hasattr(graph, 'nodes'):
                graph.nodes.discard(n)

    return graph


def pop_bubbles_simple(
    graph: DeBruijnGraph,
    max_bubble_len: int,
) -> DeBruijnGraph:
    """
    Remove very simple bubbles from the graph.

    We consider a bubble as:
        - one node splits into two alternative paths, where these paths merge into one node later
        - They are often caused by sequencing errors or SNPs

    We only take into account:
        - exactly two outgoing edges (a simple split)
        - short paths
        - same merge node

    INPUT:
        - graph: DeBruijnGraph
        - max_bubble_len: maximum path length to consider a bubble

    OUTPUT:
        - graph with weaker bubble path removed
    """

    # Iterate over all nodes as potential split points
    for u in list(graph.out_edges):

        targets = list(graph.out_edges[u].keys())
        if len(targets) != 2:  # A split node must have exactly two outgoing edges
            continue  # not a split

        v1, v2 = targets

        # Helper function to walk a path forward
        def walk(start):
            path = [start]
            weight_sum = 0
            current = start

            for _ in range(max_bubble_len):  
                out_deg = graph.out_degree.get(current, 0)  
                in_deg = graph.in_degree.get(current, 0)    
                
                if out_deg != 1 or in_deg != 1:
                    break

                next_nodes = list(graph.out_edges.get(current, {}).keys())  
                if not next_nodes:  
                    break
                
                nxt = next_nodes[0]                
                weight_sum += graph.out_edges[current][nxt]
                path.append(nxt)
                current = nxt

            # Return:
            # - end node
            # - path nodes
            # - total edge weight
            return current, path, weight_sum

        # Walk both branches
        end1, path1, w1 = walk(v1)
        end2, path2, w2 = walk(v2)

        if end1 != end2:
            continue  # If they do not merge, this is not a bubble

        #  Decide which path is weaker (lower total weight)
        weaker = path1 if w1 < w2 else path2

        # Remove nodes belonging to the weaker path
        for n in weaker[:-1]:
            # Skip if already removed
            if hasattr(graph, 'nodes') and n not in graph.nodes:
                continue
            
            # Remove outgoing edges
            if n in graph.out_edges:
                for v in list(graph.out_edges[n]):
                    if v in graph.in_degree:
                        graph.in_degree[v] -= 1
                del graph.out_edges[n]

            # Remove incoming edges
            for x in list(graph.out_edges):
                if n in graph.out_edges[x]:
                    del graph.out_edges[x][n]
                    if x in graph.out_degree:
                        graph.out_degree[x] -= 1

            # Remove degree bookkeeping
            graph.in_degree.pop(n, None)
            graph.out_degree.pop(n, None)   

            # FIX 8: Usuń z graph.nodes
            if hasattr(graph, 'nodes'):
                graph.nodes.discard(n)

    return graph
