# Algorithms for Genomic Data Analysis  
## Assignment 3: Genome Assembly  
**Winter semester 2025/2026**

---

## Task

Design and implement an assembly algorithm that processes **single-end reads** originating from the **same strand of a single chromosome**.

In your program you:
- can use the codes from the classes,
- cannot use programs and libraries for read assembly, mapping, alignment, etc.,
- cannot use multiprocessing commands.

---

## Required Deliverables

The solution should include:
- a program file `assembly`, executable using the syntax:  
  `./assembly input_reads.fasta output_contigs.fasta`
- a `README` file with a short description of your approach,
- program source code (if the executable is binary).

---

## Input and Minimum Performance Requirements

Typical parameters of the input dataset:
- number of reads: **1000**
- read length: **80 bp**
- average percentage of mismatches: **≤ 5%**
- average coverage: **≥ 5×**

The typical input dataset should be processed:
- in **less than 1 hour** on a common laptop,
- using **up to 0.5 GB of memory**.

---

## Output and Evaluation

Solutions will be evaluated on **simulated data** (artificial reads generated from a reference sequence).

Output contigs will be **locally aligned** to the reference sequence. Resulting alignments will be processed as follows:
- ambiguous alignment fragments (sharing a reference sequence interval) will be trimmed away,
- alignments of length **< 300 bp** will be excluded.

---

## Scoring

Alignments passing the filtering criteria will be scored according to the formula:

```

S = (ref_cov · cont_cov · max(0.5, 10 · (ident − 0.9))) / log5(4 + n_alignments)

```

where:
- `ref_cov` is the proportion of the reference sequence covered by the alignments,
- `cont_cov` is the proportion of the contigs’ sequence covered by the alignments,
- `ident` is the identity proportion in the alignments,
- `n_alignments` is the number of alignments.

---

## Training Data

A training data package can be downloaded from Moodle. It consists of:
- directory `reference/` containing:
  - `reference.fasta`
  - Bowtie2 index files for the reference sequence,
- directory `reads/` containing simulated read files:
  - `reads1.fasta` – 1000 reads with ~1% mismatches,
  - `reads2.fasta` – 1000 reads with 1–3% mismatches,
  - `reads3.fasta` – 1000 reads with 3–5% mismatches,
- scripts to evaluate your assembly.

The evaluation scripts require:
- **Bowtie2**
- **Python module `pysam`**

Usage:

```
./evaluate.sh contigs.fasta
```


---

## Reference Results

Two reference algorithms were tested on the training datasets:

| Dataset        | Algorithm 1 | Algorithm 2 |
|---------------|-------------|-------------|
| reads1.fasta  | 0.06        | 0.59        |
| reads2.fasta  | 0.03        | 0.21        |
| reads3.fasta  | 0.03        | 0.08        |

---

## Terms and Conditions

The assignment can be completed:
- individually,
- in 2-person teams,
- or in 3-person teams.


---

## Assessment

Every solution that meets the minimum requirements receives **3 points**.

Additional points can be awarded for:
- **assembly quality**:
  - 8 points for scoring higher than reference algorithm 2,
  - 4 points for scoring higher than reference algorithm 1,
- **meeting deadlines and presentation quality**: up to 2 points,
- **team size**:
  - 2 points for 1 person,
  - 1 point for 2 persons.
