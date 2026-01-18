#!/usr/bin/env python3

"""
Summary:
    Build a De Bruijn graph from k-mer counts.

Description:
    This module implements data structures and functions required
    to construct a weighted directed De Bruijn graph where:
        - nodes are (k-1)-mers (strings)
        - edges represent k-mers
        - edge weight equals k-mer frequency

    This module contains NO file I/O and NO CLI.
    It is intended to be imported and used by run_assembly.py.
"""

###############################
########### IMPORTS  ##########
###############################
from pathlib import Path
from collections import defaultdict, deque


###############################
######### Graph Class #########
###############################

class DeBruijnGraph:
    """
    Class stores a directed weigthed graph where:
        - nodes = (k-1)mers - strings
        - edges = kmers
        - edge direction = prefix -> suffix
        - edge weigth = how many times this kmer occures

        # in example for k=5, kmer=ACGTA
        # prefix: ACGT (k-1)mer     (u)
        # sufix CGTA                (v)
        # edge: ACGT ---(weight w)---> CGTA
    """

    def __init__(self, k: int = None):
        self.node_len = k - 1 if k else None

        self.out_edges: dict[str, dict[str, int]] = defaultdict(
            dict
        )  # for each node 'u' store all outgoing edges
        # if we assume kmers: ACGTA (count 3) and ACGTC (count 2)
        # then: out_edges = {'ACGT" : {'CGTA': 3, 'CGTC': 2}}
        # from node ACGT three edges go to CGTA and two times to "CGTC"

        self.in_degree: dict[str, int] = defaultdict(
            int
        )  # number of incoming edges from str node
        # if these edges exist: ACGT -> CGTA, TTTA -> CGTA
        # then: in_degree = {
        # "CGTA": 2,
        # "ACGT": 0,
        # "TTTA": 0}

        self.out_degree: dict[str, int] = defaultdict(
            int
        )  # number of outgoing edges from str node
        # out_degree = {
        #     "ACGT": 2,
        #     "TTTA": 1,
        #     "CGTA": 0}

        self.nodes: set[str] = set()

    def add_edge(self, u: str, v: str, weight: int):
        self.nodes.add(u)  
        self.nodes.add(v) 

        if v in self.out_edges[u]:  # checks if an edge u-> v alreade exists
            # if exists, increase the edge weight
            self.out_edges[u][v] += weight
        else:
            # if not, create new directed edge and update other dicts
            self.out_edges[u][v] = weight
            self.out_degree[u] += 1
            self.in_degree[v] += 1

        # ensure nodes exist in degree dicts. In case if some nodes do not have outgoing or incoming edges
        self.in_degree.setdefault(u, self.in_degree[u])
        self.out_degree.setdefault(v, self.out_degree[v])


###############################
######### Functions  ##########
###############################
def build_graph_from_kmers(
    kmer_counts: dict[str, int], k: int, min_kmer_count: int
) -> DeBruijnGraph:
    """
    Build a De Bruijn graph from k-mer counts.

    Args:
        kmer_counts:
            dict mapping k-mer (string) -> occurrence count
        k:
            k-mer length
        min_kmer_count:
            minimum count threshold for k-mers

    Returns:
        DeBruijnGraph object
    """
    if not kmer_counts:
        raise ValueError("kmer_counts dictionary is empty")

    graph = DeBruijnGraph(k=k)

    for kmer, count in kmer_counts.items():
        if count < min_kmer_count:
            continue

        if len(kmer) != k:
            raise ValueError(f"Invalid k-mer length: expected {k}, got {len(kmer)}")

        prefix = kmer[:-1]
        suffix = kmer[1:]

        graph.add_edge(prefix, suffix, count)
    return graph


def graph_stats(graph: DeBruijnGraph) -> dict[str, int]:
    """
    Compute basic statistics of the De Bruijn graph.

    Returns:
        dict with keys:
            - nodes: number of nodes
            - edges: number of directed edges
    """

    n_nodes = len(graph.nodes)
    n_edges = sum(len(v) for v in graph.out_edges.values())

    return {"nodes": n_nodes, "edges": n_edges}


def write_graph_stats(stats: dict[str, int], path: Path) -> None:

    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as handle:
        for key, value in stats.items():
            handle.write(f"{key}\t{value}\n")
