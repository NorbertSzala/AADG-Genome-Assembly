#!/usr/bin/env python3
"""
read_correction.py - ulepszona wersja
"""
from collections import Counter
from typing import List, Set, Dict


def count_all_kmers(reads: List[str], k: int) -> Dict[str, int]:
    """Policz wszystkie k-mery."""
    kmer_counts = Counter()
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer_counts[read[i:i+k]] += 1
    return dict(kmer_counts)


def find_trusted_kmers(kmer_counts: Dict[str, int], min_count: int) -> Set[str]:
    """K-mery o wysokiej liczności = prawdopodobnie poprawne."""
    return {kmer for kmer, count in kmer_counts.items() if count >= min_count}


def estimate_error_rate(reads: List[str], k: int = 17) -> float:
    """Oszacuj błędy na podstawie singletonów."""
    kmer_counts = count_all_kmers(reads, k)
    if not kmer_counts:
        return 0.0
    
    singletons = sum(1 for c in kmer_counts.values() if c == 1)
    total = sum(kmer_counts.values())
    
    return (singletons / total) / k if total > 0 else 0.0


def correct_read_voting(read: str, k: int, trusted_kmers: Set[str]) -> str:
    """
    Korekta przez głosowanie - każda pozycja zbiera "głosy" na bazę.
    """
    if len(read) < k:
        return read
    
    n = len(read)
    # Dla każdej pozycji: {baza: liczba_głosów}
    votes = [Counter() for _ in range(n)]
    
    # Zbierz głosy z trusted k-merów
    for i in range(n - k + 1):
        kmer = read[i:i+k]
        
        if kmer in trusted_kmers:
            # Ten k-mer jest OK - głosuj za oryginalnymi bazami
            for j in range(k):
                votes[i + j][kmer[j]] += 1
        else:
            # Spróbuj naprawić przez substytucję
            for pos in range(k):
                for base in 'ACGT':
                    candidate = kmer[:pos] + base + kmer[pos+1:]
                    if candidate in trusted_kmers:
                        # Znaleziono korektę - głosuj
                        for j in range(k):
                            votes[i + j][candidate[j]] += 1
    
    # Wybierz bazę z największą liczbą głosów
    corrected = []
    for i, original_base in enumerate(read):
        if votes[i]:
            best_base = votes[i].most_common(1)[0][0]
            corrected.append(best_base)
        else:
            corrected.append(original_base)
    
    return ''.join(corrected)


def adaptive_correction(reads: List[str]) -> List[str]:
    """
    Automatyczna wielorundowa korekcja.
    """
    if not reads:
        return 
    
    current = reads
    
    # Oszacuj błędy
    error_rate = estimate_error_rate(current, k=15)
    
    # Parametry zależne od poziomu błędów
    if error_rate < 0.015:     
        params = [(21, 2)]       
    elif error_rate < 0.03:     
        params = [(17, 3), (21, 2)]  
    else:                       
        params = [(15, 4), (17, 3), (21, 2)]  
    
    for k, min_count in params:
        kmer_counts = count_all_kmers(current, k)
        trusted = find_trusted_kmers(kmer_counts, min_count)
        
        if len(trusted) < 100:  # Za mało trusted k-merów
            continue
            
        current = [correct_read_voting(r, k, trusted) for r in current]
    
    return current