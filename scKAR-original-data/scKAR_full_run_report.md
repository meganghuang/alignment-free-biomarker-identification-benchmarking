# scKAR full run on GSE73121 — all 118 single cells

Hoffman2, `/u/scratch/m/mhuang/claude/scKAR`. Pipeline job 15261046, one attempt,
Sat 2026-10-24 02:41:30 -> 18:23:32 PDT (**15 h 42 min**).

## Inputs

- **118 single cells**, the complete GSE73121 / SRP063840 single-cell set
  (47 PDX_pRCC, 36 PDX_mRCC, 35 Pt_mRCC; the 3 pooled libraries excluded).
  74 GiB of FASTQ, every file byte-exact against ENA's reported sizes.
- Gene x cell TPM matrix, 57,820 x 118, from GEO's per-cell RSEM files.
- **Reference filtering ON**: GENCODE v50 transcripts -> 6.88 GB sorted 31-mer
  list (~215 M k-mers), byte-order verified. All 63 adjacency chunks were passed
  through `ref_filter`, so the tested k-mers are non-reference sequence.
- hg38 `genome.2bit` for BLAT.
- Defaults otherwise: k=31, leiden resolution 0.6, N_PCS=40, N_NEIGHBORS=10,
  PSEUDOBULK_SIZE=500000, MIN_ROW_THREDSHOLD=20, LOG2FC=1.0, PVAL=0.05, hmp.

## Clustering and the contrast tested

Leiden at the **default resolution 0.6** gave 3 clusters, and they track the
study design closely:

| cluster | n | composition |
|---|---|---|
| 0 | 45 | PDX_pRCC (45 of 47) |
| 1 | 44 | Pt_mRCC (35) + PDX_mRCC (9) |
| 2 | 29 | PDX_mRCC (27) + PDX_pRCC (2) |

Three dendrogram leaves give exactly one internal edge, so there is a single
bipartition: **A = clusters 1+2 (73 cells, mRCC-dominant) vs B = cluster 0
(45 cells, PDX_pRCC)**. That is the primary disease contrast of the study,
recovered without supervision.

## Results

| stage | output | wall time |
|---|---|---|
| jellyfish + filter + adjacency + reference filtering | 118 k-mer tables (187 GB); 63 adjacency chunks (115 GB), all reference-filtered | 11 h |
| clustering | 3 clusters, 1 bipartition | 81 s |
| F-test (ANOVA + BH) | 69 matrix chunks, 2.1 GB | 2 h |
| DESeq2 (LRT, glmGamPoi) | 68 result tables, 434 MB | ~2.7 h incl. assembly |
| merged DE k-mers | **2,571,204** at padj <= 0.05 and abs(log2FC) >= 1 | |
| ABySS unitigs (k=27) | 8,180 + 23,562 = **31,742** contigs, 1.56 Mb assembled | |
| BLAT vs hg38 | 153,523 + 141,717 blast9 hit lines | |
| contig p-values (harmonic mean) | all 31,742 contigs scored | |

Contig lengths: set A n=8,180, median 50 bp, max 1,606 bp, 572,803 bp total;
set B n=23,562, median 43 bp, max 1,956 bp, 986,081 bp total.

Disk footprint ~305 GB, dominated by the per-cell k-mer tables (187 GB) and the
adjacency (115 GB).

## Direction of the split

The same label inversion found in the 16-cell run applies: `results(dds,
name="conditionB")` makes negative log2FoldChange mean "higher in A", while
`merge_deseq_results.py` writes positive log2FC to `A_kmers.fasta`. So:

- `A_kmers.fasta` — 1,375,438 k-mers — **enriched in cluster 0 (PDX_pRCC)**
- `B_kmers.fasta` — 1,195,766 k-mers — **enriched in clusters 1+2 (mRCC)**

This follows from the coefficient convention plus the per-cell TPM check done on
the 16-cell run; it was not re-verified against this run's matrices.

## Statistical caveats — read 2.57 M with care

1. **Multiple testing is corrected per chunk, not globally.** Each of the 69
   `f_test_result_*.csv` files is a separate DESeq2 run, so `padj` is a BH
   correction within that chunk (~78 k rows each after filtering), and
   `merge_deseq_results.py` then concatenates the already-adjusted values. There
   is no correction across the ~5.4 M k-mers tested overall. The same applies to
   the BH step inside the F-test, which runs per pseudobulk batch.
2. **k-mers are not independent.** Adjacent 31-mers from one transcript rise and
   fall together, so 2.57 M significant k-mers represent far fewer independent
   signals. The contig-level harmonic-mean p-values exist precisely to collapse
   this — 31,742 contigs is the more meaningful unit, and 1.56 Mb of assembled
   sequence the more meaningful scale.
3. The A/B groups are **clusters, not randomised conditions**, and they were
   derived from gene expression of the same cells. Testing k-mer abundance
   between clusters defined on those cells is circular to a degree; the contrast
   happens to align with the study's own PDX-pRCC / mRCC distinction, which is
   what makes it interpretable.

## Known issues carried into this run

All patched in the clone, documented in the 16-cell report: the missing `wait`
in `mother_script.sh`, unsorted `tpm_sum.csv`, broken resume guards, the C++17
flags, undefined `BLAT_PATH`, non-executable `blat.sh`, and `sc.tl.tsne`'s
small-n crash. One more was needed for this run: **`ref_filter.cpp` read
`ref_data[ref_index]` before its bounds check**, which is undefined behaviour at
the tail of every adjacency chunk once the reference list is exhausted. Patched
to test the index first; with that fix all 63 chunks filtered cleanly.

One of the 69 F-test chunks produced no DESeq2 table (68 of 69).
