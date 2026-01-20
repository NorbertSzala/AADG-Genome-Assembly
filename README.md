# Genome Assembly Project – Script Overview

This project implements a simple **de novo genome assembler** based on a **De Bruijn graph**.  
Input: single-end reads from **one strand** of one chromosome.  
Output: assembled contigs in FASTA format.


---

## run_assembly.py
**Main entry point**

**What it does**
- Parses command-line arguments.
- Runs the full assembly pipeline.
- Optionally runs parameter optimization.

**How it works**
- Calls `run_assembly()` from `assembly_core.py`.
- Passes user parameters (k-mer length, tip removal, bubble popping, etc.).
- Validates input and creates output directories.

This is the main script you execute from the terminal.

---

## assembly_core.py
**Core assembly pipeline.**

**What it does**
- Connects all steps into one workflow:
  1. Read FASTA
  2. Correct reads
  3. Count k-mers
  4. Build De Bruijn graph
  5. Clean graph
  6. Extract contigs
  7. Write output and statistics

**How it works**
- Saves reports (`report.txt`), parameters (`params_used.json`), and statistics.
- Returns contig statistics such as total length and N50.

Main logic of the assembler.

---

## io_fasta.py
**FASTA input/output utilities.**

**What it does**
- Reads FASTA files into a list of sequences.
- Writes sequences back to a FASTA file.

**How it works**
- Uses `Bio.SeqIO` only for parsing.
- Normalizes sequences (uppercase, no whitespace).
- Keeps file handling separate from algorithmic code.

---

## kmers.py
**k-mer generation and statistics.**

**What it does**
- Generates k-mers from reads.
- Counts how many times each k-mer appears.
- Builds a histogram of k-mer frequencies.

**How it works**
- Slides a window of length `k` over each read.
- Uses dictionaries for fast counting.
- Writes the histogram to a TSV file.

This step helps identify sequencing noise.

---

## dbg.py
**De Bruijn graph implementation.**

**What it does**
- Defines the `DeBruijnGraph` class.
- Builds a weighted directed graph from k-mers.

**How it works**
- Nodes are `(k-1)`-mers.
- Edges represent k-mers.
- Edge weight equals k-mer frequency.
- Stores in-degree and out-degree for fast checks.

---

## cleaning.py
**Graph cleaning algorithms.**

**What it does**
- Removes common graph errors:
  - Small disconnected components (islands)
  - Short dead ends (tips)
  - Simple bubbles (optional)

**How it works**
- Use BFS to find connected components.
- Traces short paths to detect tips.
- Uses edge weights (coverage) to decide which paths to remove.
- Bubble popping removes the weaker path.

This improves contig quality.

---

## read_correction.py
**Read error correction.**

**What it does**
- Corrects sequencing errors before graph construction.

**How it works**
- Estimates error rate using rare k-mers.
- Finds “trusted” k-mers with high counts.
- Corrects reads by base voting from trusted k-mers.
- Runs multiple correction rounds depending on error level.

This reduces noise early in the pipeline.

---

## traversal.py
**Contig extraction from the graph.**

**What it does**
- Traverses the cleaned graph to build contigs.

**How it works**
- Finds start nodes (branches, sources, sinks).
- Walks uniquely determined paths in both directions.
- Converts node paths into DNA sequences.
- Filters contigs shorter than a minimum length.
- Computes contig statistics (N50, total length, etc.).

This is where final sequences are assembled.

---

## optimize_parameters.py
**Automatic parameter optimization.**

**What it does**
- Tests many parameter combinations.
- Selects the best assembly result.

**How it works**
- Runs the assembler multiple times with different parameters.
- Scores results using:
  - total contig length (coverage)
  - N50 (continuity)
- Keeps only the best run and removes the rest.

Useful for finding good parameters automatically.
