#!/usr/bin/env python3
"""
Korekta błędów w readach na podstawie spektrum k-merów.

Strategia:
- K-mery o wysokiej liczności = prawidłowe
- K-mery o niskiej liczności = błędy
- Korekta przez substytucję baz
"""

from collections import Counter
from typing import List, Set


def count_all_kmers(reads: List[str], k: int) -> dict[str, int]:
    """Policz wszystkie k-mery we wszystkich readach."""
    kmer_counts = Counter()
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer_counts[read[i:i+k]] += 1
    return dict(kmer_counts)


def find_trusted_kmers(kmer_counts: dict[str, int], min_count: int = 3) -> Set[str]:
    """Znajdź zaufane k-mery (wysoką liczność = prawdopodobnie poprawne)."""
    return {kmer for kmer, count in kmer_counts.items() if count >= min_count}


def correct_kmer(kmer: str, trusted_kmers: Set[str]) -> str:
    """
    Spróbuj skorygować k-mer przez testowanie wszystkich substytucji 1-bazowych.
    Zwraca poprawiony k-mer jeśli znaleziony, inaczej oryginalny.
    """
    if kmer in trusted_kmers:
        return kmer
    
    bases = ['A', 'C', 'G', 'T']
    for i in range(len(kmer)):
        original = kmer[i]
        for new_base in bases:
            if new_base == original:
                continue
            candidate = kmer[:i] + new_base + kmer[i+1:]
            if candidate in trusted_kmers:
                return candidate
    
    return kmer


def correct_read(read: str, k: int, trusted_kmers: Set[str]) -> str:
    """Koryguj błędy w pojedynczym readzie używając spektrum k-merów."""
    read_list = list(read)
    
    for i in range(len(read) - k + 1):
        kmer = read[i:i+k]
        if kmer not in trusted_kmers:
            corrected = correct_kmer(kmer, trusted_kmers)
            if corrected != kmer:
                for j in range(k):
                    if kmer[j] != corrected[j]:
                        read_list[i+j] = corrected[j]
                        break
    
    return ''.join(read_list)


def estimate_error_rate(reads: List[str], k: int = 17) -> float:
    """
    Oszacuj procent błędów na podstawie spektrum k-merów.
    K-mery o liczności=1 to prawdopodobnie błędy.
    """
    kmer_counts = count_all_kmers(reads, k)
    singletons = sum(1 for count in kmer_counts.values() if count == 1)
    total_kmers = sum(kmer_counts.values())
    
    if total_kmers == 0:
        return 0.0
    
    error_rate = (singletons / total_kmers) / k
    return min(error_rate, 0.2)


def adaptive_correction(reads: List[str]) -> List[str]:
    """
    Automatyczna korekta błędów - wybiera parametry na podstawie danych.
    
    Proces:
    1. Oszacuj procent błędów
    2. Wybierz k, min_count na podstawie błędów
    3. Uruchom korekta (1-2 rundy)
    
    Args:
        reads: Lista sekwencji readów
    
    Returns:
        Lista skorygowanych readów
    """
    # Oszacuj błędy
    error_rate = estimate_error_rate(reads, k=17)
    
    # Wybierz parametry na podstawie błędów
    if error_rate < 0.02:  # <2% błędów
        k, min_count, rounds = 21, 2, 1
    elif error_rate < 0.04:  # 2-4% błędów
        k, min_count, rounds = 17, 3, 2
    else:  # >4% błędów
        k, min_count, rounds = 15, 4, 2
    
    # Korekta wielorundowa
    current_reads = reads
    for _ in range(rounds):
        kmer_counts = count_all_kmers(current_reads, k)
        trusted_kmers = find_trusted_kmers(kmer_counts, min_count)
        current_reads = [correct_read(read, k, trusted_kmers) for read in current_reads]
    
    return current_reads