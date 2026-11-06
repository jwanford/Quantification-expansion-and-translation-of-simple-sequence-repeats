################################################################################
#                                                                              #
#           BACTERIAL PHASE VARIATION & SIMPLE SEQUENCE REPEAT ANALYSIS       #
#                                                                              #
#  This pipeline analyses repetitive DNA tracts in bacterial gene sequences   #
#  and their role in mediating phase variation via slipped-strand mispairing. #
#                                                                              #
#  WORKFLOW:                                                                   #
#    STEP 1 — Count maximum consecutive nucleotide runs (A, C, G, T)          #
#    STEP 2 — In silico G-tract expansion (+1, +2, +3 nt) to simulate         #
#             slipped-strand mispairing events                                 #
#    STEP 3 — Translate each sequence (original + expanded) and calculate     #
#             protein length as a % of maximum possible length to infer       #
#             ON/OFF expression state and flag phase-variable candidates       #
#                                                                              #
#  INPUT:  Multi-FASTA file of bacterial gene sequences (DNA)                 #
#  OUTPUT: Three output files (Steps 1 & 3 = CSV, Step 2 = FASTA)            #
#                                                                              #
#  Author:  [Joseph Wanford]                                                        #
#  Date:    [16.4.26]                                                             #
#  Notes:   Uses NCBI Genetic Code Table 11 (bacterial/archaeal)              #
#                                                                              #
################################################################################


# ==============================================================================
# 0. SETUP — Install and load required packages
# ==============================================================================

options(repos = c(CRAN = "https://cran.rstudio.com/"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}

library(Biostrings)

cat("=== Packages loaded successfully ===\n\n")


# ==============================================================================
# 1. USER CONFIGURATION — !! EDIT THIS SECTION BEFORE RUNNING !!
# ==============================================================================

# ------------------------------------------------------------------------------
# >>>  PATH 1 OF 4 — YOUR INPUT FILE                                        <<<
# >>>  Set this to the full path of your multi-FASTA file of DNA sequences  <<<
# ------------------------------------------------------------------------------
INPUT_FASTA  <- "/path/to/your/input_sequences.fasta"            # <-- CHANGE THIS

# ------------------------------------------------------------------------------
# >>>  PATH 2 OF 4 — STEP 1 OUTPUT                                          <<<
# >>>  CSV file: maximum consecutive nucleotide run counts per sequence     <<<
# ------------------------------------------------------------------------------
OUTPUT_STEP1 <- "/path/to/your/output_repeat_counts.csv"         # <-- CHANGE THIS

# ------------------------------------------------------------------------------
# >>>  PATH 3 OF 4 — STEP 2 OUTPUT                                          <<<
# >>>  FASTA file: original + G-expanded sequences (_G+0 to _G+3)          <<<
# ------------------------------------------------------------------------------
OUTPUT_STEP2 <- "/path/to/your/output_g_expanded_sequences.fasta" # <-- CHANGE THIS

# ------------------------------------------------------------------------------
# >>>  PATH 4 OF 4 — STEP 3 OUTPUT                                          <<<
# >>>  CSV file: translation results and phase variation calls              <<<
# ------------------------------------------------------------------------------
OUTPUT_STEP3 <- "/path/to/your/output_phase_variation_calls.csv" # <-- CHANGE THIS

# ------------------------------------------------------------------------------
# >>>  NUCLEOTIDE TO EXPAND — REQUIRED                                      <<<
# >>>  Specify which homopolymeric tract to expand in Step 2.               <<<
# >>>  Must be exactly one of: "A", "C", "G", or "T"                       <<<
# ------------------------------------------------------------------------------
EXPAND_NUCLEOTIDE <- "G"   # <-- CHANGE THIS to "A", "C", "G", or "T"

# ------------------------------------------------------------------------------
# Threshold: minimum % of maximum protein length to call a sequence as ON.
# Sequences below this threshold are called OFF (truncated / frameshifted).
# Default is 95% — adjust only if needed for your gene of interest.
# ------------------------------------------------------------------------------
ON_THRESHOLD <- 95   # percent

# --- Validate EXPAND_NUCLEOTIDE immediately so errors are caught early -------
if (!EXPAND_NUCLEOTIDE %in% c("A", "C", "G", "T")) {
  stop(paste0(
    "Invalid EXPAND_NUCLEOTIDE: \"", EXPAND_NUCLEOTIDE, "\"\n",
    "  Please set EXPAND_NUCLEOTIDE to exactly one of: \"A\", \"C\", \"G\", \"T\""
  ))
}
cat(paste0("  Nucleotide tract selected for expansion: ", EXPAND_NUCLEOTIDE, "\n\n"))


# ==============================================================================
# STEP 1 — COUNT MAXIMUM CONSECUTIVE NUCLEOTIDE RUNS
# ==============================================================================
# For each sequence in the FASTA file, this step finds the longest unbroken
# run of each nucleotide (A, C, G, T). These homopolymeric tracts are the
# substrates for slipped-strand mispairing during replication.
#
# readDNAStringSet() handles all standard FASTA formats natively, including
# sequences wrapped across multiple lines — no pre-processing of the FASTA
# file is needed or performed.
#
# Output: CSV with columns — SequenceName, Max_A, Max_C, Max_G, Max_T
# ==============================================================================

cat("==============================================================\n")
cat("STEP 1: Counting maximum consecutive nucleotide runs\n")
cat("==============================================================\n")


# --- Helper: find the longest run of a given nucleotide in a sequence ------

count_max_run <- function(sequence, char) {
  sequence_str  <- as.character(sequence)
  if (nchar(sequence_str) == 0) return(0)

  chars         <- unlist(strsplit(sequence_str, split = ""))
  max_count     <- 0
  current_count <- 0

  for (ch in chars) {
    if (ch == char) {
      current_count <- current_count + 1
      if (current_count > max_count) max_count <- current_count
    } else {
      current_count <- 0
    }
  }
  return(max_count)
}


# --- Main Step 1 function --------------------------------------------------

run_step1 <- function(fasta_file, output_file) {

  sequences <- readDNAStringSet(fasta_file)
  cat(paste("  Sequences detected:", length(sequences), "\n"))

  results <- data.frame(
    SequenceName = character(),
    Max_A        = integer(),
    Max_C        = integer(),
    Max_G        = integer(),
    Max_T        = integer(),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(sequences)) {
    seq_name <- names(sequences)[i]
    max_A    <- count_max_run(sequences[i], "A")
    max_C    <- count_max_run(sequences[i], "C")
    max_G    <- count_max_run(sequences[i], "G")
    max_T    <- count_max_run(sequences[i], "T")

    results <- rbind(results, data.frame(
      SequenceName = seq_name,
      Max_A = max_A,
      Max_C = max_C,
      Max_G = max_G,
      Max_T = max_T,
      stringsAsFactors = FALSE
    ))

    cat(paste0("  [", i, "/", length(sequences), "] ", seq_name,
               " | A:", max_A, " C:", max_C, " G:", max_G, " T:", max_T, "\n"))
  }

  write.csv(results, file = output_file, row.names = FALSE)
  cat(paste("\n  Step 1 complete. Results saved to:", output_file, "\n\n"))
  return(results)
}

# Run Step 1
step1_results <- run_step1(INPUT_FASTA, OUTPUT_STEP1)


# ==============================================================================
# STEP 2 — IN SILICO G-TRACT EXPANSION (+1, +2, +3 NUCLEOTIDES)
# ==============================================================================
# To simulate slipped-strand mispairing events that expand the G-tract by 1,
# 2, or 3 nucleotides, this step:
#   (a) Locates the longest G-run in each sequence
#   (b) Inserts 1, 2, or 3 additional G nucleotides immediately after that run
#   (c) Writes all expanded sequences (original + 3 variants) to a FASTA file
#
# Each output sequence is labelled with a suffix: _{NT}+0, _{NT}+1, _{NT}+2,
# _{NT}+3, where {NT} is the nucleotide chosen in EXPAND_NUCLEOTIDE.
# These labelled sequences are passed directly to Step 3 for translation.
#
# Note: if a sequence contains multiple runs of equal maximum length,
# the first occurrence is expanded (conservative approach).
#
# Output: FASTA file containing original + all expanded variants
# ==============================================================================

cat("==============================================================\n")
cat("STEP 2: In silico G-tract expansion (+1, +2, +3 nt)\n")
cat("==============================================================\n")


# --- Helper: remove any non-ACGTN characters from a sequence ---------------

clean_sequence <- function(sequence) {
  sequence_str     <- as.character(sequence)
  cleaned_sequence <- gsub("[^ACGTN]", "", sequence_str)
  return(DNAString(cleaned_sequence))
}


# --- Helper: insert n additional nucleotides after the longest homopolymer run
# nt        : the nucleotide to expand — "A", "C", "G", or "T"
# n_extra   : number of additional nucleotides to insert

expand_tract <- function(sequence, nt, n_extra) {
  sequence_str <- as.character(sequence)

  pattern     <- paste0(nt, "+")
  runs        <- gregexpr(pattern, sequence_str)[[1]]
  run_lengths <- attr(runs, "match.length")

  if (runs[1] == -1) {
    warning(paste0("  No ", nt, "-run found in sequence — returning original."))
    return(DNAString(sequence_str))
  }

  max_length   <- max(run_lengths)
  max_position <- which(run_lengths == max_length)[1]
  start_pos    <- runs[max_position]
  end_pos      <- start_pos + max_length - 1

  new_sequence <- paste0(
    substr(sequence_str, 1, end_pos),
    strrep(nt, n_extra),
    substr(sequence_str, end_pos + 1, nchar(sequence_str))
  )

  return(DNAString(new_sequence))
}


# --- Main Step 2 function --------------------------------------------------
# nt: the nucleotide whose longest homopolymeric run will be expanded.
#     Passed in from EXPAND_NUCLEOTIDE set in the config section above.

run_step2 <- function(fasta_file, output_file, nt) {

  sequences <- readDNAStringSet(fasta_file)
  cat(paste("  Sequences detected:", length(sequences), "\n"))
  cat(paste0("  Expanding longest ", nt, "-tract by +1, +2, +3 nt\n\n"))

  # Collect sequences and names as plain character vectors, then build the
  # DNAStringSet once at the end. Incremental c() on DNAStringSet objects
  # corrupts the internal GroupedIRanges index on large datasets.
  seq_char_vec <- character()
  name_vec     <- character()

  for (i in seq_along(sequences)) {
    seq_name <- names(sequences)[i]
    cleaned  <- clean_sequence(sequences[[i]])

    cat(paste0("  Processing: ", seq_name, "\n"))

    original   <- as.character(cleaned)
    expanded_1 <- as.character(expand_tract(cleaned, nt, 1))
    expanded_2 <- as.character(expand_tract(cleaned, nt, 2))
    expanded_3 <- as.character(expand_tract(cleaned, nt, 3))

    # Log original tract length and resulting sequence lengths
    run_orig <- max(attr(gregexpr(paste0(nt, "+"), original)[[1]], "match.length"))
    cat(paste0("    Original longest ", nt, "-run: ", run_orig, " nt\n"))
    cat(paste0("    ", nt, "+1 sequence length:    ", nchar(expanded_1), "\n"))
    cat(paste0("    ", nt, "+2 sequence length:    ", nchar(expanded_2), "\n"))
    cat(paste0("    ", nt, "+3 sequence length:    ", nchar(expanded_3), "\n"))

    seq_char_vec <- c(seq_char_vec, original, expanded_1, expanded_2, expanded_3)
    name_vec     <- c(name_vec,
                      paste0(seq_name, "_", nt, "+0"),
                      paste0(seq_name, "_", nt, "+1"),
                      paste0(seq_name, "_", nt, "+2"),
                      paste0(seq_name, "_", nt, "+3"))
  }

  names(seq_char_vec) <- name_vec
  all_sequences       <- DNAStringSet(seq_char_vec)

  writeXStringSet(all_sequences, filepath = output_file)

  cat(paste("\n  Step 2 complete.", length(all_sequences),
            "sequences written to:", output_file, "\n\n"))
  return(all_sequences)
}

# Run Step 2 — EXPAND_NUCLEOTIDE is passed in from the config section
step2_sequences <- run_step2(INPUT_FASTA, OUTPUT_STEP2, EXPAND_NUCLEOTIDE)


# ==============================================================================
# STEP 3 — TRANSLATION AND PHASE VARIATION CALLING
# ==============================================================================
# Each sequence (original _G+0, and expansions _G+1, _G+2, _G+3) is
# translated using NCBI Genetic Code Table 11 (standard bacterial code).
#
# For each translated protein, the script finds the longest ORF from any
# methionine (M) to the next stop codon (*), and expresses its length as a
# percentage of the maximum possible protein length (DNA length / 3).
#
#   >= ON_THRESHOLD%  -> likely IN-FRAME (full-length protein)  -> ON
#   <  ON_THRESHOLD%  -> likely FRAMESHIFTED (truncated protein) -> OFF
#
# PHASE VARIATION CALLING:
# A gene is flagged as PHASE VARIABLE if there is any change in ON/OFF status
# across the four expansion states (_G+0 through _G+3). This captures both:
#   - OFF -> ON transitions (expansion restores reading frame)
#   - ON  -> OFF transitions (expansion disrupts reading frame)
# A gene that is uniformly ON or uniformly OFF across all states is NOT
# considered phase variable by G-tract slippage.
#
# OUTPUT FORMAT: one row per gene, with expansion results as separate columns,
# repeat counts from Step 1 merged in, and a single PhaseVariable call.
# Columns: SequenceName, Max_A, Max_C, Max_G, Max_T,
#          G0_Pct, G0_State, G1_Pct, G1_State,
#          G2_Pct, G2_State, G3_Pct, G3_State, PhaseVariable
# ==============================================================================

cat("==============================================================\n")
cat("STEP 3: Translation and phase variation calling\n")
cat("==============================================================\n")


# --- Helper: find the longest in-frame ORF using DNA codons directly -------
#
# Works directly on the DNA codon sequence rather than the translated protein,
# so that all NCBI Genetic Code Table 11 start and stop codons are recognised
# explicitly without relying on translate() to render them correctly.
#
# NCBI Genetic Code 11 — start codons:
#   ATG (canonical), GTG, TTG, CTG, ATT, ATC, ATA
#
# NCBI Genetic Code 11 — stop codons:
#   TAA, TAG, TGA
#
# For every in-frame start codon, the function searches for the next
# downstream in-frame stop codon and records the ORF length (number of
# coding codons, excluding the stop). Returns the longest such ORF, or 0
# if no valid start-to-stop ORF is found.
#
# dna_sequence : character string of DNA (will be trimmed to codon boundary)

longest_orf_length <- function(dna_sequence) {

  # Bacterial start and stop codon sets — NCBI Genetic Code Table 11
  start_codons <- c("ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA")
  stop_codons  <- c("TAA", "TAG", "TGA")

  # Trim to nearest codon boundary and split into a vector of codons
  trimmed  <- substr(dna_sequence, 1, (nchar(dna_sequence) %/% 3) * 3)
  n_codons <- nchar(trimmed) %/% 3
  if (n_codons == 0) return(0)

  starts <- seq(1, n_codons * 3 - 2, by = 3)
  ends   <- seq(3, n_codons * 3,     by = 3)
  codons <- substring(trimmed, starts, ends)

  max_length <- 0

  for (i in seq_along(codons)) {
    if (codons[i] %in% start_codons) {
      # Search for the next in-frame stop codon from this start position
      for (j in i:length(codons)) {
        if (codons[j] %in% stop_codons) {
          orf_length <- j - i     # coding codons only, stop not counted
          if (orf_length > max_length) max_length <- orf_length
          break
        }
      }
    }
  }

  return(max_length)
}


# --- Main Step 3 function --------------------------------------------------
# Takes the expanded sequences from Step 2 and the repeat counts from Step 1,
# translates all variants, and returns one summary row per gene.

run_step3 <- function(expanded_sequences, repeat_counts, nt, on_threshold, output_file) {

  cat(paste("  Sequences to translate:", length(expanded_sequences), "\n"))
  cat(paste("  ON/OFF threshold: >=", on_threshold, "% of maximum protein length\n\n"))

  # --- Translate all sequences and collect raw per-variant results ----------
  raw <- data.frame(
    SequenceName     = character(),
    Expansion        = character(),
    PctOfMaximum     = numeric(),
    ExpressionState  = character(),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(expanded_sequences)) {
    full_name    <- names(expanded_sequences)[i]
    dna_sequence <- as.character(expanded_sequences[[i]])
    dna_length   <- nchar(dna_sequence)

    # Parse gene name and expansion label from the sequence name suffix.
    # The suffix is _{NT}+0 ... _{NT}+3 where {NT} is EXPAND_NUCLEOTIDE.
    suffix_pattern  <- paste0("_", nt, "[+][0-9]+")
    expansion_label <- sub(paste0(".*(", suffix_pattern, ")$"), "\\1", full_name)
    seq_name        <- sub(paste0(suffix_pattern, "$"),         "",      full_name)

    # longest_orf_length works directly on DNA codons — no pre-translation needed
    max_protein_len <- longest_orf_length(dna_sequence)
    max_possible    <- dna_length %/% 3

    pct_of_max <- if (max_possible > 0) {
      round((max_protein_len / max_possible) * 100, 2)
    } else { 0 }

    expression_state <- if (pct_of_max >= on_threshold) "ON" else "OFF"

    raw <- rbind(raw, data.frame(
      SequenceName    = seq_name,
      Expansion       = expansion_label,
      PctOfMaximum    = pct_of_max,
      ExpressionState = expression_state,
      stringsAsFactors = FALSE
    ))

    cat(paste0("  ", full_name, " | ORF: ", max_protein_len,
               " aa | ", pct_of_max, "% | ", expression_state, "\n"))
  }

  # --- Reshape to one row per gene with expansion columns -------------------
  # Each gene gets four columns per expansion: G0_Pct, G0_State, etc.

  unique_genes <- unique(raw$SequenceName)
  expansions   <- paste0("_", nt, "+", 0:3)
  col_prefixes <- paste0(nt,  0:3)

  wide <- data.frame(SequenceName = unique_genes, stringsAsFactors = FALSE)

  for (k in seq_along(expansions)) {
    exp_label  <- expansions[k]
    col_pct    <- paste0(col_prefixes[k], "_PctOfMax")
    col_state  <- paste0(col_prefixes[k], "_State")

    exp_rows   <- raw[raw$Expansion == exp_label, ]

    wide[[col_pct]]   <- exp_rows$PctOfMaximum[   match(wide$SequenceName, exp_rows$SequenceName)]
    wide[[col_state]] <- exp_rows$ExpressionState[ match(wide$SequenceName, exp_rows$SequenceName)]
  }

  # --- Merge Step 1 repeat counts -------------------------------------------
  # Brings Max_A, Max_C, Max_G, Max_T into the final table.
  # The SequenceName in repeat_counts may include extra annotation after a
  # space (common in FASTA headers); match on the portion before the first
  # space to be safe.

  repeat_counts$MatchName <- sub(" .*", "", repeat_counts$SequenceName)
  wide$MatchName           <- sub(" .*", "", wide$SequenceName)

  wide <- merge(wide,
                repeat_counts[, c("MatchName", "Max_A", "Max_C", "Max_G", "Max_T")],
                by = "MatchName", all.x = TRUE)
  wide$MatchName <- NULL   # drop the temporary matching column

  # Reorder columns for readability: name, repeats, then expansion results
  pct_cols   <- paste0(nt, 0:3, "_PctOfMax")
  state_cols_ordered <- paste0(nt, 0:3, "_State")
  expand_cols <- as.vector(rbind(pct_cols, state_cols_ordered))
  wide <- wide[, c("SequenceName", "Max_A", "Max_C", "Max_G", "Max_T", expand_cols)]

  # --- Phase variation calling ----------------------------------------------
  # A gene is PHASE VARIABLE if ON/OFF status is not identical across all
  # four expansion states. This flags both:
  #   OFF->ON (expansion restores frame — classic PV by slippage)
  #   ON->OFF (expansion disrupts frame — also indicative of SSR-mediated PV)
  # Genes uniformly ON or uniformly OFF across all expansions are NOT flagged.

  state_cols <- paste0(nt, 0:3, "_State")

  wide$PhaseVariable <- apply(wide[, state_cols], 1, function(states) {
    states <- states[!is.na(states)]
    if (length(unique(states)) > 1) "Yes" else "No"
  })

  # Sort: phase variable genes first, then alphabetically by name
  wide <- wide[order(wide$PhaseVariable == "No", wide$SequenceName), ]

  write.csv(wide, file = output_file, row.names = FALSE)

  n_pv <- sum(wide$PhaseVariable == "Yes")
  cat(paste0("\n  Step 3 complete.\n"))
  cat(paste0("  Phase variable genes detected: ", n_pv, " of ", nrow(wide), "\n"))
  cat(paste0("  Results saved to: ", output_file, "\n\n"))

  return(wide)
}

# Run Step 3 — pass step1_results so repeat counts are merged into the output
step3_results <- run_step3(step2_sequences, step1_results, EXPAND_NUCLEOTIDE, ON_THRESHOLD, OUTPUT_STEP3)


# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("==============================================================\n")
cat("PIPELINE COMPLETE\n")
cat("==============================================================\n")
cat(paste("  Step 1 output (repeat counts):    ", OUTPUT_STEP1, "\n"))
cat(paste("  Step 2 output (expanded FASTA):   ", OUTPUT_STEP2, "\n"))
cat(paste("  Step 3 output (phase var. calls): ", OUTPUT_STEP3, "\n"))
cat("\nPhase variable gene summary:\n")

phase_var_genes <- step3_results[step3_results$PhaseVariable == "Yes", ]

if (nrow(phase_var_genes) > 0) {
  for (j in seq_len(nrow(phase_var_genes))) {
    g <- phase_var_genes[j, ]
    cat(paste0("  * ", g$SequenceName, "\n"))
    max_run_col <- paste0("Max_", EXPAND_NUCLEOTIDE)
    cat(paste0("      Max ", EXPAND_NUCLEOTIDE, "-run: ", g[[max_run_col]], " nt\n"))
    for (n in 0:3) {
      pct_col   <- paste0(EXPAND_NUCLEOTIDE, n, "_PctOfMax")
      state_col <- paste0(EXPAND_NUCLEOTIDE, n, "_State")
      cat(paste0("      ", EXPAND_NUCLEOTIDE, "+", n, ": ",
                 g[[pct_col]], "% -> ", g[[state_col]], "\n"))
    }
  }
} else {
  cat("  No phase variable genes identified under current threshold.\n")
}

cat("\nDone.\n")
