# Bacterial Phase Variation & Simple Sequence Repeat Analysis Pipeline

A single R script for identifying phase-variable genes in bacterial genomes by detecting homopolymeric repeat tracts and testing whether in silico expansion of those tracts restores or disrupts the reading frame.

---

## Background

Phase variation is the reversible, high-frequency switching of gene expression between ON and OFF states. In many bacteria, this is mediated by slipped-strand mispairing during replication of simple sequence repeats (SSRs) or short homopolymeric tracts of A, C, G, or T. Expansion or contraction of these tracts by a single nucleotide can shift the reading frame, creating a premature stop codon (OFF) or restoring a full-length protein (ON).

This pipeline takes a multi-FASTA file of bacterial gene sequences and:
1. Counts the longest consecutive run of each nucleotide in every sequence
2. Simulates SSR expansion by inserting +1, +2, and +3 copies of a chosen nucleotide into the longest tract
3. Identifies the longest open reading frame in each variant and expresses it as a percentage of the maximum possible protein length, using this to call ON/OFF expression state and flag phase-variable candidates

---

## Requirements

- **R** (version 4.0 or later recommended)
- **Bioconductor** and the **Biostrings** package

These are installed automatically the first time the script is run if not already present.

---

## Input

A **multi-FASTA file** of bacterial gene sequences (DNA). Sequences can be in single-line or multi-line (wrapped) FASTA format — both are handled natively.

Example:
```
>geneA
ATGAAACGCTTTATTTTTGTTTTCGCTTTTTTTTTTTTTTGCAAAAATAA
>geneB_allele1
GTGAAAGCACTCTTTCATTTATTTATTGTTTTATTTCTG
TTTTATATATCAGTTTATTGTTTTTATTCATATGTGTCA
```

All sequences should be submitted in the same reading frame (i.e. as annotated CDSs, not whole-gene loci with upstream/downstream sequence).

---

## Setup — Edit These Four Lines Before Running

Open `phase_variation_analysis.R` in RStudio and edit the configuration block near the top of the script. There are **four paths** and **one nucleotide** to set:

```r
INPUT_FASTA        <- "/path/to/your/sequences.fasta"
OUTPUT_STEP1       <- "/path/to/your/repeat_counts.csv"
OUTPUT_STEP2       <- "/path/to/your/expanded_sequences.fasta"
OUTPUT_STEP3       <- "/path/to/your/phase_variation_calls.csv"

EXPAND_NUCLEOTIDE  <- "G"    # Change to "A", "C", "G", or "T"
```

**`INPUT_FASTA`** — full path to your input multi-FASTA file.

**`OUTPUT_STEP1`** — where to save the repeat count table (CSV).

**`OUTPUT_STEP2`** — where to save the expanded sequence FASTA (used internally and useful for inspection).

**`OUTPUT_STEP3`** — where to save the final results table (CSV). This is the main output file.

**`EXPAND_NUCLEOTIDE`** — the nucleotide whose longest homopolymeric tract will be expanded. Must be exactly one of `"A"`, `"C"`, `"G"`, or `"T"`. Choose based on the repeat type associated with your gene(s) of interest (e.g. `"G"` for G-tracts in *rmpA2*, `"T"` for T-tracts in *rmpD*).

There is also one optional parameter:

```r
ON_THRESHOLD  <- 95    # percent — default 95%
```

This sets the minimum percentage of maximum possible protein length required to call a sequence ON. Sequences below this threshold are called OFF. The default of 95% is appropriate for most applications, but this can be lowered it only if your gene of interest has a known short C-terminal region that is not essential for function.

---

## Running the Script

Open `phase_variation_analysis.R` in RStudio, edit the configuration block as described above, then click **Source** (or press `Ctrl+Shift+Enter` / `Cmd+Shift+Enter` on Mac) to run the entire script.

Progress is printed to the R console as each sequence is processed.

---

## Output

### Step 1 output — `repeat_counts.csv`

One row per input sequence. Records the longest consecutive run of each nucleotide.

| Column | Description |
|---|---|
| `SequenceName` | FASTA header of the sequence |
| `Max_A` | Longest consecutive A-run (nt) |
| `Max_C` | Longest consecutive C-run (nt) |
| `Max_G` | Longest consecutive G-run (nt) |
| `Max_T` | Longest consecutive T-run (nt) |

### Step 2 output — `expanded_sequences.fasta`

A multi-FASTA file containing the original sequence plus three expanded variants for every input sequence, labelled `_G+0`, `_G+1`, `_G+2`, `_G+3` (or `_A+0` etc. depending on `EXPAND_NUCLEOTIDE`). Useful for manual inspection or downstream analysis.

### Step 3 output — `phase_variation_calls.csv`

One row per input gene, combining repeat counts and translation results across all four expansion states. **Phase-variable genes are sorted to the top.**

| Column | Description |
|---|---|
| `SequenceName` | Gene/allele name from FASTA header |
| `Max_A/C/G/T` | Longest run of each nucleotide (from Step 1) |
| `G0_PctOfMax` | % of maximum protein length at original sequence |
| `G0_State` | ON or OFF at original sequence |
| `G1_PctOfMax` | % of maximum protein length after +1 expansion |
| `G1_State` | ON or OFF after +1 expansion |
| `G2_PctOfMax` | % of maximum protein length after +2 expansion |
| `G2_State` | ON or OFF after +2 expansion |
| `G3_PctOfMax` | % of maximum protein length after +3 expansion |
| `G3_State` | ON or OFF after +3 expansion |
| `PhaseVariable` | **Yes** if ON/OFF status changes across any expansion; **No** if uniform |

*Column prefixes (G0, G1...) reflect the chosen `EXPAND_NUCLEOTIDE` — e.g. T0, T1... if `EXPAND_NUCLEOTIDE <- "T"`.*

---

## How Phase Variable Calls Are Made

A gene is called **Phase Variable = Yes** if the ON/OFF state is not identical across all four expansion states (original + three expansions). This captures:

- **OFF → ON**: the original sequence is frameshifted (truncated protein), and expansion by +1, +2, or +3 nt restores the reading frame — the classic signature of SSR-mediated phase variation
- **ON → OFF**: the original sequence is in-frame, but one or more expansions disrupt the reading frame — also indicative of SSR-mediated switching

A gene that is uniformly ON or uniformly OFF across all four states is **not** flagged as phase variable by this analysis.

### Translation and ORF detection

All six NCBI Genetic Code Table 11 alternative start codons are recognised (`ATG`, `GTG`, `TTG`, `CTG`, `ATT`, `ATC`, `ATA`), alongside all three stop codons (`TAA`, `TAG`, `TGA`). ORF detection works directly on DNA codons rather than the translated protein to ensure no start codons are missed.

---

## Important Caveats

- The pipeline expands only the **longest** tract of the chosen nucleotide in each sequence. If two tracts are tied in length, the first occurrence is expanded.
- Contraction events (SSR shortening) are not simulated — only expansion by +1, +2, +3.
- The tool is designed for annotated CDS sequences submitted in the correct reading frame. Submitting whole-gene loci or intergenic sequence will produce unreliable results.
- The `%_of_maximum` metric assumes the maximum possible protein runs the full length of the submitted sequence. Genes with naturally short C-termini may give values slightly below 100% even when fully in-frame — adjust `ON_THRESHOLD` accordingly if needed.

---

## Troubleshooting

| Problem | Likely cause | Fix |
|---|---|---|
| `Error reading FASTA file` | Wrong path or file not found | Double-check `INPUT_FASTA` path, including file extension |
| All sequences show 0% | Sequences not in reading frame, or no stop codon within submitted length | Check input sequences are annotated CDSs in frame |
| `Invalid EXPAND_NUCLEOTIDE` error | Value not in `"A","C","G","T"` | Check for typos and ensure quotes are included |
| No phase variable genes detected | Threshold too high, or wrong nucleotide selected | Try lowering `ON_THRESHOLD` or check which nucleotide the repeat tract consists of |
| Fewer sequences than expected | FASTA file has formatting issues | Ensure each sequence has a header line beginning with `>` |
