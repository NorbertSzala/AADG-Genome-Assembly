#!/usr/bin/env python3
"""
Automatic parameter optimization for assembly.

Tests different parameter combinations and finds the best one.
Metric: maximize total contig length + N50
"""

import subprocess
import json
from pathlib import Path
import sys


def run_assembly(input_fasta, output_dir, k, min_kmer_count, 
                min_component_size, tip_max_len, pop_bubbles, 
                max_bubble_len, min_contig_len):
    """
    Run assembly with given parameters.
    Returns dict with metrics or None if failed.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        "python3", "run_assembly.py",
        "-I", str(input_fasta),
        "-O", str(output_dir),
        "-K", str(k),
        "-M", str(min_kmer_count),
        "--min-component-size", str(min_component_size),
        "--tip-max-len", str(tip_max_len),
        "--min-contig-len", str(min_contig_len),
    ]
    
    if pop_bubbles:
        cmd.extend(["--pop-bubbles", "--max-bubble-len", str(max_bubble_len)])
    
   
    subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    
    # Read report to get metrics
    report_file = output_dir / "report.txt"
    if not report_file.exists():
        return None
    metrics = parse_report(report_file)

    metrics['contigs_file'] = output_dir / "contigs.fasta"

    return metrics
        

def parse_report(report_path):
    """Extract metrics from report.txt"""
    metrics = {
        'num_contigs': 0,
        'total_length': 0,
        'longest': 0,
        'n50': 0
    }
    
    with open(report_path, 'r') as f:
        for line in f:
            line = line.strip()
            if 'num_contigs:' in line:
                metrics['num_contigs'] = int(line.split(':')[1].strip())
            elif 'total_length:' in line:
                metrics['total_length'] = int(line.split(':')[1].strip())
            elif 'longest:' in line:
                metrics['longest'] = int(line.split(':')[1].strip())
            elif 'n50:' in line:
                metrics['n50'] = int(line.split(':')[1].strip())
    
    return metrics


def score_result(metrics):
    """
    Calculate score for parameter combination.
    
    Score = total_length * 0.6 + N50 * 0.3 + num_contigs_penalty
    
    Higher is better.
    """
    if metrics is None or metrics['num_contigs'] == 0:
        return 0
    
    # We want: high total length, high N50, fewer contigs
    total_score = metrics['total_length'] * 0.6
    n50_score = metrics['n50'] * 0.3
    
    # Penalty for too many contigs (fragmented assembly)
    contig_penalty = -abs(metrics['num_contigs'] - 10) * 10
    
    return total_score + n50_score + contig_penalty


def optimize_parameters(input_fasta, final_output_dir="best_assembly"):
    """
    Test multiple parameter combinations and save ONLY best.
    """
    print("="*70)
    print("PARAMETER OPTIMIZATION")
    print("="*70)
    print(f"Input: {input_fasta}")
    print()
    
    # Parameter grid (bez zmian)
    param_grid = {
        'k': [15, 17], #, 19, 21, 23],
        'min_kmer_count': [2], #3, 4],
        'min_component_size': [5, 10],
        'tip_max_len': [2], #3, 5],
        'pop_bubbles': [True, False],
        'max_bubble_len': [5, 8],
        'min_contig_len': [300]
    }
    
    # Generate combinations (bez zmian)
    combinations = []
    for k in param_grid['k']:
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
                                'k': k,
                                'min_kmer_count': min_kmer,
                                'min_component_size': min_comp,
                                'tip_max_len': tip_len,
                                'pop_bubbles': pop_bub,
                                'max_bubble_len': bubble_len,
                                'min_contig_len': 300
                            })
    
    print(f"Testing {len(combinations)} combinations...")
    print()
    
    # Test all in temp directory
    results = []
    temp_dir = Path("temp_optimization")
    
    for i, params in enumerate(combinations, 1):
        print(f"[{i}/{len(combinations)}] k={params['k']}, M={params['min_kmer_count']}")
        
        output_dir = temp_dir / f"run_{i:03d}"
        metrics = run_assembly(input_fasta, output_dir, **params)
        
        if metrics:
            score = score_result(metrics)
            print(f"    → total={metrics['total_length']}, N50={metrics['n50']}, score={score:.0f}")
            results.append({
                'params': params,
                'metrics': metrics,
                'score': score,
                'output_dir': output_dir
            })
        else:
            print(f"    → FAILED")
    
    if not results:
        print("\nNo successful runs!")
        return None
    
    # Find best
    best = max(results, key=lambda x: x['score'])
    
    # Copy ONLY best result
    import shutil
    
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
    shutil.rmtree(temp_dir)
    
    print("\n" + "="*70)
    print("BEST PARAMETERS")
    print("="*70)
    for key, value in best['params'].items():
        print(f"  {key}: {value}")
    print()
    print("Metrics:")
    for key, value in best['metrics'].items():
        if key != 'contigs_file':
            print(f"  {key}: {value}")
    print()
    print(f"✓ Saved to: {final_dir}/")
    print("="*70)
    
    return best
def main():
    if len(sys.argv) < 2:
        print("Usage: python3 optimize_parameters.py <input.fasta> [--quick]")
        print()
        print("Options:")
        print("  (default)  Full optimization (~2 hours, 120 tests)")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    

    optimize_parameters(input_fasta)


if __name__ == "__main__":
    main()