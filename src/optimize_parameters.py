#!/usr/bin/env python3
"""
Automatic parameter optimization for assembly.

Tests different parameter combinations and finds the best one.
Metric: maximize total contig length + N50
"""

import json
from pathlib import Path
import sys
import shutil
import math

# Import run_assembly function directly
from assembly_core import run_assembly


# def score_result(metrics):
#     """
#     Calculate score for parameter combination.
    
#     Score = total_length * 0.6 + N50 * 0.3 + num_contigs_penalty
    
#     Higher is better.
#     """
#     if metrics is None or metrics.get('num_contigs', 0) == 0:
#         return 0
    
#     # We want: high total length, high N50, fewer contigs
#     total_score = metrics['total_length'] * 0.6
#     n50_score = metrics['n50'] * 0.3
    
#     # Penalty for too many contigs (fragmented assembly)
#     contig_penalty = -abs(metrics['num_contigs'] - 10) * 10
    
#     return total_score + n50_score + contig_penalty

def score_result(metrics):
    """
    Ocena bez znajomości długości genomu i referencji.
    """
    if metrics is None or metrics.get('num_contigs', 0) == 0:
        return 0
    
    # 1. N50 i Najdłuższy kontig
    continuity_score = math.log10(metrics['n50'] * metrics['longest'] + 1)
    
    # 2. Kara za fragmentację
    fragmentation_penalty = 1.0 / math.log(metrics['num_contigs'] + 4, 5)
    
    # 3. Bonus za średnią długość
    mean_bonus = min(1.0, metrics['mean_length'] / 1000)
    
    return continuity_score * fragmentation_penalty * (1 + mean_bonus)


def optimize_parameters(input_fasta, final_output_dir="best_assembly"):
    """
    Test multiple parameter combinations and save ONLY best.
    
    Args:
        input_fasta: Path to input FASTA file
        final_output_dir: Directory where best assembly will be saved
        
    Returns:
        dict: Best result with params, metrics, and score
    """
    input_fasta = Path(input_fasta)
    

    # Parameter grid
    # param_grid = {
    #     'kmer_length': [15, 17, 29, 31],
    #     'min_kmer_count': [2],
    #     'min_component_size': [3, 5, 10],
    #     'tip_max_len': [1, 2],
    #     'pop_bubbles': [True, False],
    #     'max_bubble_len': [5, 8],
    #     'min_contig_len': [300]
    # }
    # param_grid = {
    #     'kmer_length': [15, 17, 21, 25, 29, 31, 35],  # Szerszy zakres, niższe k dla wyższego error rate
    #     'min_kmer_count': [2, 3, 4, 5],  # Wyższe wartości filtrują błędne k-mery (error rate 5%)
    #     'min_component_size': [3, 5, 8, 10, 15],  # Więcej opcji dla filtrowania małych wysp
    #     'tip_max_len': [1, 2, 3, 5],  # Dłuższe tipy przy wyższym error rate
    #     'pop_bubbles': [True, False],
    #     'max_bubble_len': [3, 5, 8, 10],  # Dłuższe bubble paths przy error rate 5%
    #     'min_contig_len': [300]  # Różne progi długości
    # }
    param_grid = {
        'kmer_length': [15],  # Szerszy zakres, niższe k dla wyższego error rate
        'min_kmer_count': [2, 3],  # Wyższe wartości filtrują błędne k-mery (error rate 5%)
        'min_component_size': [3, 5],  # Więcej opcji dla filtrowania małych wysp
        'tip_max_len': [1, 2],  # Dłuższe tipy przy wyższym error rate
        'pop_bubbles': [True, False],
        'max_bubble_len': [3, 5, 8, 10],  # Dłuższe bubble paths przy error rate 5%
        'min_contig_len': [300]  # Różne progi długości
    }
    
    # Generate combinations
    combinations = []
    for k in param_grid['kmer_length']:
        for min_kmer in param_grid['min_kmer_count']:
            for min_comp in param_grid['min_component_size']:
                for tip_len in param_grid['tip_max_len']:
                    for pop_bub in param_grid['pop_bubbles']:
                        if pop_bub:
                            bubble_lens = param_grid['max_bubble_len']
                        else:
                            bubble_lens = [5]
                        for bubble_len in bubble_lens:
                            combinations.append({
                                'kmer_length': k,
                                'min_kmer_count': min_kmer,
                                'min_component_size': min_comp,
                                'tip_max_len': tip_len,
                                'pop_bubbles': pop_bub,
                                'max_bubble_len': bubble_len,
                                'min_contig_len': 300
                            })
    
    
    # Test all in temp directory
    results = []
    temp_dir = Path("temp_optimization")
    best_score = 0
    best_result = None

    for i, params in enumerate(combinations, 1):
        print(f"[{i}/{len(combinations)}] k={params['kmer_length']}, M={params['min_kmer_count']}")
        
        output_dir = temp_dir / f"run_{i:03d}"
        output_dir.mkdir(parents=True, exist_ok=True) 

        try:
            # Call run_assembly directly with parameters
            metrics = run_assembly(
                input=input_fasta,
                outdir=output_dir,
                **params
            )
            
            if metrics and metrics.get('num_contigs', 0) > 0:
                score = score_result(metrics)
                print(f"    → total={metrics['total_length']}, N50={metrics['n50']}, score={score:.0f}")
                results.append({
                    'params': params,
                    'metrics': metrics,
                    'score': score,
                    'output_dir': output_dir
                })

                if score > best_score:
                    best_score = score
                    # Usuń stary best
                    if best_result and best_result['output_dir'].exists():
                        shutil.rmtree(best_result['output_dir'])
                    
                    best_result = {
                        'params': params,
                        'metrics': metrics,
                        'score': score,
                        'output_dir': output_dir
                    }
                else:
                    # Usuń gorszy wynik od razu
                    shutil.rmtree(output_dir)
        
            else:
                print(f"FAILED (no contigs)")
                
        except Exception as e:
            print(f"FAILED ({e})")
    
    if not results:
        print("\nNo successful runs!")
        return None
    
    # Find best
    best = max(results, key=lambda x: x['score'])
    
    # Copy ONLY best result
    final_dir = Path(final_output_dir)
    if final_dir.exists():
        shutil.rmtree(final_dir)
    final_dir.mkdir(parents=True, exist_ok=True)
    
    for file in best['output_dir'].glob("*"):
        shutil.copy(file, final_dir / file.name)
    
    # Save best params
    with open(final_dir / "best_params.json", 'w') as f:
        json.dump(best['params'], f, indent=2)
    
    # Clean up ALL temp files
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
    
    print("BEST PARAMETERS")
    print("="*70)
    for key, value in best['params'].items():
        print(f"  {key}: {value}")
    print()
    print("Metrics:")
    for key, value in best['metrics'].items():
        print(f"  {key}: {value}")
    print()
    print(f"Saved to: {final_dir}/")
    
    return best


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 optimize_parameters.py <input.fasta>")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    optimize_parameters(input_fasta)


if __name__ == "__main__":
    main()