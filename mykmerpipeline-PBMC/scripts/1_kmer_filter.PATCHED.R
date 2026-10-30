# ============================================================
# Script: 1_kmer_filter.R
# Author: Megan Huang
#
# Description:
# [OPTIONAL] Loads a k-mer count matrix and removes highly correlated,
# overlapping k-mers based on suffix-prefix matching.
# Correlated k-mers above a specified threshold are filtered
# to reduce redundancy in downstream analyses.
#
# Inputs:
#   1. kmer_matrix.mtx     - Sparse matrix of k-mer counts (Matrix Market format)
#   2. cell_barcodes.txt   - List of cell/sample barcodes (columns)
#   3. kmer_ids.txt        - List of k-mer sequences (rows)
#
# Parameters:
#   - cor_threshold        - Correlation cutoff for filtering [default: 0.95]
#   - num_cores            - Number of CPU cores for parallel processing
#   - OUTDIR               - Directory containing inputs + outputs [MUST EDIT]
#
# Outputs:
#   - filtered_kmer_ids.txt        : Filtered list of k-mers
#   - filtered_kmer_matrix.mtx     : Filtered sparse k-mer matrix
#
# TODO:
#   - permanently incorporate into nextflow script (0_create_kmer_counts.nf) 
#     once confirming whether this consistently helps on more sample datasets
#     
# ============================================================

library(Matrix)
library(stringr)
library(data.table)
library(parallel)

# ============================
# PARAMETERS
# ============================

cor_threshold <- 0.95  # correlation threshold
num_cores <- detectCores() - 1

# Centralized directory (EDIT THIS ONLY)
OUTDIR <- "/u/scratch/m/mhuang/scripts/output_files_default"  # [EDIT THIS BASED ON NF WORK FOLDER]

# ============================
# Load input
# ============================

kmer_matrix <- readMM(file.path(OUTDIR, "kmer_matrix.mtx"))
barcodes   <- readLines(file.path(OUTDIR, "cell_barcodes.txt"))
kmer_ids   <- readLines(file.path(OUTDIR, "kmer_ids.txt"))

rownames(kmer_matrix) <- kmer_ids
colnames(kmer_matrix) <- barcodes

# Convert to dense if manageable
dense_mat <- as.matrix(kmer_matrix)

# ============================
# Generate suffix-prefix match table
# ============================

k <- nchar(kmer_ids[1])
kmer_dt <- data.table(
  id = seq_along(kmer_ids),
  kmer = kmer_ids,
  prefix = substr(kmer_ids, 1, k - 1),
  suffix = substr(kmer_ids, 2, k)
)

# Get overlapping k-mer pairs (i, j)
overlap_pairs <- kmer_dt[kmer_dt, on = .(suffix = prefix), allow.cartesian=TRUE, nomatch=0]
setnames(overlap_pairs, c("query_id", "query_kmer", "query_prefix", "query_suffix",
                          "target_id", "target_kmer", "target_suffix"))

# Remove self-matches or backward comparisons
overlap_pairs <- overlap_pairs[query_id < target_id]

# ============================
# Parallel correlation filtering
# ============================

# BUGFIX: setnames() above renamed `id` to `query_id`, so overlap_pairs$id was NULL
# and seq_along(NULL) gave a zero-length grouping factor -- split() returned no
# chunks, mclapply ran over nothing, and to_remove stayed all-FALSE (no filtering).
chunks <- split(seq_len(nrow(overlap_pairs)),
                ceiling(seq_len(nrow(overlap_pairs)) / 1000))

to_remove <- rep(FALSE, length(kmer_ids))

removal_list <- mclapply(chunks, function(idx_range) {
  local_remove <- rep(FALSE, length(kmer_ids))
  for (idx in idx_range) {
    i <- overlap_pairs$query_id[idx]
    j <- overlap_pairs$target_id[idx]
    if (local_remove[j] || local_remove[i]) next
    cor_val <- suppressWarnings(cor(dense_mat[i, ], dense_mat[j, ]))
    if (!is.na(cor_val) && cor_val > cor_threshold) {
      local_remove[j] <- TRUE
    }
  }
  return(local_remove)
}, mc.cores = num_cores)

# Combine results
for (r in removal_list) {
  to_remove <- to_remove | r
}

# ============================
# Save outputs
# ============================

filtered_mat   <- kmer_matrix[!to_remove, ]
filtered_kmers <- kmer_ids[!to_remove]

writeLines(filtered_kmers, file.path(OUTDIR, "filtered_kmer_ids.txt"))
writeMM(as(filtered_mat, "dgTMatrix"),
        file.path(OUTDIR, "filtered_kmer_matrix.mtx"))