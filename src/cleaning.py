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
########## ARGUMENTS ##########
###############################


###############################
########## FUNCTIONS ##########
###############################


# Helper function
def _neighbors_undirected(graph: DeBruijnGraph, node: str) -> set[str]:
    neighbors = set()

    # outgoing neighbors
    if node is graph.out_degree:
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

    Generally, real genome formsone large coneccted graph but with many 'mistakes' (small disconected subraphs).
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

    nodes = set(graph.in_degree) | set(graph.out_degree)  # all nodes in the graph

    # check each node as a potential tip start
    for node in list(nodes):
        # A dead-end node has no incoming OR no outgoing edges
        if graph.out_degree[node] == 0 or graph.in_degree[node] == 0:
            path = [node]
            current = node

            while True:
                # Stop if this node is a branch (not a simple path)
                if graph.in_degree[current] != 1 or graph.out_degree[current] != 1:
                    break

                # follow unique edge
                next_nodes = list(graph.out_edges.get(current, {}).keys())
                if not next_nodes:
                    break

                next_node = next_nodes[0]
                path.append(next_node)
                current = next_node

                # stop if the path gets too long
                if len(path) > tip_max_len:
                    break

            # If the path is short, we consider it a tip
            if len(path) <= tip_max_len:
                # Remove all nodes in the tip path
                for n in path:
                    # remove outgoing
                    if n in graph.out_edges:
                        for v in list(graph.out_edges[n]):
                            graph.in_degree[v] -= 1
                        del graph.out_edges[n]

                    # remove incoming
                    for u in list(graph.out_edges):
                        if n in graph.out_edges[u]:
                            del graph.out_edges[u][n]
                            graph.out_degree[u] -= 1

                    # Remove degree bookkeeping
                    graph.in_degree.pop(n, None)
                    graph.out_degree.pop(n, None)

    return graph


def pop_bubbles_simple(
    graph: DeBruijnGraph,
    max_bubble_len: int,
) -> DeBruijnGraph:
    """
    Remove very simple bubbles from the graph.

    We consider a booble as:
        - one node splits into to alternative paths, where these paths merge into one node later
        - They are often caused by sequencing errors or SNPs

    We only tak into account:
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

            while True:
                # Stop if branch or dead-end
                if graph.out_degree[current] != 1 or graph.in_degree[current] != 1:
                    break
                nxt = list(graph.out_edges[current].keys())[0]
                weight_sum += graph.out_edges[current][nxt]
                path.append(nxt)
                current = nxt

                # Stop if path too long
                if len(path) > max_bubble_len:
                    break
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
        for n in weaker:
            # Remove outgoing edges
            if n in graph.out_edges:
                for v in list(graph.out_edges[n]):
                    graph.in_degree[v] -= 1
                del graph.out_edges[n]

            # Remove incoming edges
            for x in list(graph.out_edges):
                if n in graph.out_edges[x]:
                    del graph.out_edges[x][n]
                    graph.out_degree[x] -= 1

            # Remove degree bookkeeping
            graph.in_degree.pop(n, None)
            graph.out_degree.pop(n, None)

    return graph
