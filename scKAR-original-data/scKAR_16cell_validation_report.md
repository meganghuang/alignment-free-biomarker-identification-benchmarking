# scKAR run on GSE73121 — 16-cell validation

Run on Hoffman2, project root `/u/scratch/m/mhuang/claude/scKAR`.
Completed 2026-10-22 (driver jobs 15242968 + 15246051).

## Dataset

GSE73121 / SRP063840 — Kim et al., single-cell transcriptomes of metastatic
renal cell carcinoma patient-derived cells. This is the source of the example
accessions in scKAR's own README (GSM1887215-18). Illumina HiSeq 2500,
paired-end, one library per cell — the layout scKAR expects.

16 cells, balanced across the two PDX groups:

| group | condition | cells |
|---|---|---|
| PDX_mRCC | A | SC_04, 21, 34, 52, 60, 63, 90, 92 |
| PDX_pRCC | B | SC_13, 24, 34, 35, 37, 58, 70, 73 |

Inputs assembled at `data/`:
- `reads/<cell>.fastq.gz` — R1 and R2 concatenated per cell, from ENA (9.3 GB)
- `expression_matrix/expression.csv` — 57,820 genes x 16 cells, built from the
  per-cell RSEM TPM files in GEO's GSE73121_RAW.tar
- `genome_assembly/genome.2bit` — hg38 from UCSC (797 MB)
- `metadata/clustering.csv`, `metadata/bipartitions.csv` — group labels

## Configuration

`MODE="CUSTOM_METADATA"`, `FILTER_REFERENCE_KMERS=false`, k=31,
`PSEUDOBULK_SIZE=500000`, `MIN_ROW_THREDSHOLD=20`, `MIN_ROW_COUNT=10`,
`MIN_COL_COUNT=10`, `LOG2FC=1.0`, `PVAL=0.05`, `PVAL_COMBINATION_METHOD=hmp`,
`NUMBER_OF_THREADS=8`.

`MIN_ROW_THREDSHOLD=20` rather than the default 40 follows the authors' own
comment at the foot of `adj_to_mat.cpp` ("renal cell pre_filter = 20").

## What the run produced

| stage | output |
|---|---|
| jellyfish 31-mer count + filter | 16 per-cell k-mer tables, 23 GB; 4.2-5.0e8 retained k-mer instances per cell |
| adjacency list | 8 chunks, 18 GB |
| conditions | 1 bipartition, 8 A vs 8 B |
| F-test (one-way ANOVA + BH) | 80 matrix chunks, 17,518 k-mers retained |
| DESeq2 (LRT, glmGamPoi) | 79 result tables |
| merged DE k-mers | **7,233** at padj <= 0.05 and abs(log2FC) >= 1 |
| ABySS unitigs (k=27) | 11 + 13 contigs |
| BLAT vs hg38 | 57 + 50 blast9 hit lines |
| contig p-values (harmonic mean) | 11 + 13 contigs scored |

Contig lengths: 57-217 bp (median 69) for the `A_contigs.fasta` set,
27-545 bp (median 52) for `B_contigs.fasta`.

Results live at `data/final_results/bipartition_0/`:
`filtered_kmers.tsv`, `A_kmers.fasta`/`A_kmers.tsv`, `B_kmers.fasta`/`B_kmers.tsv`,
and `abyss/` with the contigs, their BLAT alignments and HMP p-values.

## Direction of the A/B split — read the labels with care

`deseq2_glmgamPoi.R` extracts `results(dds, name="conditionB")` with levels
`c("A","B")`, so a NEGATIVE log2FoldChange means higher in condition A.
`merge_deseq_results.py` then writes `log2FoldChange > 0` into `A_kmers.fasta`
and `< 0` into `B_kmers.fasta` — i.e. the file names are inverted with respect
to the condition labels.

Spot check on the three most significant negative-log2FC k-mers, using their
per-cell TPM rows from the F-test matrices (columns 1-8 = PDX_mRCC = A,
9-16 = PDX_pRCC = B):

| k-mer (first 20 nt) | mean TPM mRCC | mean TPM pRCC |
|---|---|---|
| TCCACGGGTCATAGGTAAAA... | 0.9 | 0.1 |
| TCCTCAAGATACAATAACTT... | 0.9 | 0.1 |
| TCCATCAACTTCTATAAAAG... | 1.0 | 0.1 |

All three are enriched in mRCC (condition A) yet are written to
`B_kmers.fasta`. So: the 4,829 k-mers in `B_kmers.fasta` are the mRCC-enriched
set and the 2,404 in `A_kmers.fasta` are the pRCC-enriched set. This rests on
the DESeq2 coefficient convention plus the three-k-mer spot check above; it was
not recomputed across all 7,233 k-mers.

## Fixes required to get the pipeline through

1. **Missing wait in `mother_script.sh`** (correctness). The per-cell filtering
   loop runs 5-wide in the background and the script proceeds to
   `create_adj_list` without waiting for the final batch. With 16 cells the
   adjacency was built from **15**, silently dropping one. Added `check_pids`.
2. **`tpm_sum.csv` is unsorted** (correctness). Parallel workers append in
   arbitrary order, and that row order defines the adjacency column order,
   while DESeq2 labels those columns from the name-sorted bipartition file.
   Unsorted, sample labels and per-cell TPM totals attach to the wrong columns.
   Now sorted (LC_ALL=C) before the adjacency step.
3. **Resume guards test the wrong filename.** `jellyfish_count.sh` and
   `mother_script.sh` look for `<base>_1_filtered.csv`, but `filter_kmers.py`
   writes `<base>.fa.csv_1_filtered.csv`, so no stage was ever skippable and
   every rerun redid the 19-minute counting step. Both fixed.
4. **C++17.** The sources use `std::filesystem`; the repo's bare `g++ x.cpp -o y`
   calls now carry `-O2 -std=c++17 -pthread -static-libstdc++ -static-libgcc
   ... -lstdc++fs`, and `gcc/9.5.0` must be loaded.
5. **`BLAT_PATH` undefined** in `config.env` though `run.sh` passes it on.
6. **`blat.sh` is not executable** in the repo and `create_contigs.sh` invokes it
   as `./blat.sh` — the BLAT step failed with "Permission denied" and the whole
   run still reported success. `chmod +x` fixes it.
7. **`sc.tl.tsne` has no small-n guard**: sklearn requires perplexity < n_samples,
   so with 16 cells `graph-clustering.py` aborted *after* clustering but *before*
   writing `cluster.csv`. Patched to `perplexity=min(30, n_obs-1)`.

## Why CUSTOM_METADATA rather than the default GENE_EXPRESSION mode

In `GENE_EXPRESSION` mode, bipartitions come from internal edges of the cluster
dendrogram. With 16 cells leiden gave one cluster at the default resolution
0.6 (dendrogram impossible) and two clusters at 1.0 — and a two-leaf dendrogram
has **no internal edges**, so `graph-clustering.py` emitted zero bipartitions
and DESeq2 had nothing to test. The dendrogram route needs >= 3 clusters.
Those two leiden clusters also split 12/4 across the biology rather than
separating mRCC from pRCC, so they would not have been the contrast of interest.
`CUSTOM_METADATA` with the published group labels gives the balanced 8-vs-8
comparison directly.

## Caveats

- 16 cells is a smoke test, not a finding. Read the 7,233 k-mers as evidence
  the pipeline runs end to end on real data, not as a biological result.
- Reference k-mer filtering was off, so these k-mers are not restricted to
  non-reference sequence — the thing scKAR is actually designed to isolate.
- Reads were not trimmed; the authors' `misc/trimmomatic_*_qc.sh` step was skipped.
- `jellyfish dump` was given `-L 2` and `-s` raised to 1G. The `-L 2` exactly
  matches the `count > 1` filter `filter_kmers.py` applies anyway, but without
  it the per-cell dumps reach ~8 GB and the pandas filter step will not fit.
- One of the 80 F-test chunks produced no DESeq2 table (79 of 80).
