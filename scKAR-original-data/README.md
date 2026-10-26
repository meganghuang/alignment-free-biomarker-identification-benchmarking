# scKAR on its original dataset (GSE73121)

Output of running [scKAR](https://github.com/RedwanulKarim612/scKAR) end to end on
**GSE73121 / SRP063840** — Kim et al., single-cell transcriptomes of metastatic renal
cell carcinoma patient-derived cells. This is the dataset the scKAR README's own example
accessions (GSM1887215-18) come from, so it serves as the reference point for
benchmarking other alignment-free biomarker methods.

Run on UCLA Hoffman2, 2026-10-24, 15 h 42 min wall, single pass.

## What was run

- **All 118 single cells** (47 PDX_pRCC, 36 PDX_mRCC, 35 Pt_mRCC; the 3 pooled
  libraries excluded). 74 GiB of FASTQ from ENA, every file byte-exact against the
  ENA `fastq_bytes` field (see `metadata/ena_expected_bytes.tsv`).
- `MODE="GENE_EXPRESSION"` — conditions derived unsupervised by leiden clustering of
  the gene x cell TPM matrix, at scKAR's **default resolution 0.6**.
- `FILTER_REFERENCE_KMERS=true` against GENCODE v50 transcripts (6.88 GB sorted
  31-mer list, ~215 M k-mers), so the tested k-mers are non-reference sequence.
- k=31, `PSEUDOBULK_SIZE=500000`, `MIN_ROW_THREDSHOLD=20`, `LOG2FC=1.0`, `PVAL=0.05`,
  harmonic-mean p-value combination. Full parameters in `metadata/config.env`.

## The contrast

Leiden gave 3 clusters that track the study design:

| cluster | n | composition |
|---|---|---|
| 0 | 45 | PDX_pRCC (45 of 47) |
| 1 | 44 | Pt_mRCC (35) + PDX_mRCC (9) |
| 2 | 29 | PDX_mRCC (27) + PDX_pRCC (2) |

Three dendrogram leaves give exactly one internal edge, hence a single bipartition:
**A = clusters 1+2 (73 mRCC-dominant cells) vs B = cluster 0 (45 PDX_pRCC)** — the
study's own disease contrast, recovered without supervision.

## Results

| quantity | value |
|---|---|
| differential k-mers (padj <= 0.05, abs(log2FC) >= 1) | 2,571,204 |
| ABySS unitigs (k=27) | 31,742 (8,180 + 23,562) |
| assembled sequence | 1.56 Mb |
| contig lengths | A: median 50 bp, max 1,606; B: median 43 bp, max 1,956 |
| BLAT vs hg38 | 153,523 + 141,717 blast9 hit lines |

**Label direction.** scKAR's A/B file names are inverted relative to its condition
labels: `deseq2_glmgamPoi.R` reads `results(dds, name="conditionB")` with levels
`c("A","B")`, so negative log2FoldChange means higher in condition A, but
`merge_deseq_results.py` writes positive log2FC into `A_kmers.fasta`. Therefore:

- `A_kmers.fasta.gz` (1,375,438 k-mers) and `contigs/A_*` are **enriched in cluster 0
  (PDX_pRCC)**
- `B_kmers.fasta.gz` (1,195,766 k-mers) and `contigs/B_*` are **enriched in clusters
  1+2 (mRCC)**

Verified empirically on the 16-cell validation run (the three strongest negative-log2FC
k-mers averaged 0.9-1.0 TPM in mRCC vs ~0.1 in pRCC, yet were written to
`B_kmers.fasta`); not re-verified against the full run's matrices.

## Contents

```
results/contigs/   assembled unitigs, their harmonic-mean p-values, and BLAT hg38 alignments
results/kmers/     the differential 31-mers (gzipped FASTA)
metadata/          cluster assignments, the bipartition, per-cell k-mer TPM totals, config
logs/              pipeline, reference-k-mer and download-verification logs
patches/           local source changes vs upstream, plus the four job scripts used
```

`patches/scKAR-local-changes.diff` is the diff against upstream
`088c3092` ("fix environment setup", 2026-03-08) and is required to reproduce this run.

## Not included (too large for git)

| file | size | location |
|---|---|---|
| `filtered_kmers.tsv` | 327 MB (138 MB gzipped) | `data/final_results/bipartition_0/` |
| `A_kmers.tsv` / `B_kmers.tsv` | 174 / 153 MB | same |
| per-cell k-mer tables | 187 GB | `data/jellyfish/` |
| k-mer x cell adjacency | 115 GB | `data/adj/` |

All on Hoffman2 under `/u/scratch/m/mhuang/claude/scKAR/` — note that is **purge-eligible
scratch**. `filtered_kmers.tsv` is the per-k-mer DESeq2 table (baseMean,
log2FoldChange, lfcSE, stat, pvalue, padj) and is the one worth preserving elsewhere,
via Git LFS or an archive.

## Caveats for benchmarking

1. **Multiple testing is corrected per chunk, not globally.** Each of the 69
   `f_test_result_*.csv` files is a separate DESeq2 run, so `padj` is a BH correction
   within ~78 k rows; `merge_deseq_results.py` concatenates already-adjusted values.
   There is no correction across the ~5.4 M k-mers tested. The F-test's own BH step is
   likewise per pseudobulk batch. Any cross-method comparison should account for this.
2. **k-mers are not independent** — adjacent 31-mers from one transcript co-vary, so
   2.57 M significant k-mers represent far fewer independent signals. The 31,742
   contigs (1.56 Mb) are the more defensible unit.
3. **The groups are clusters, not randomised conditions**, and were derived from gene
   expression of the same cells, so the test is circular to a degree. It is
   interpretable here only because the split happens to align with the study's own
   PDX-pRCC / mRCC distinction.
4. Reads were **not trimmed** (the authors' `misc/trimmomatic_*_qc.sh` step was skipped).
5. `jellyfish dump` was run with `-L 2` and `-s 1G`. The `-L 2` exactly matches the
   `count > 1` filter `filter_kmers.py` applies anyway, but without it the per-cell
   dumps reach ~8 GB each.

See `scKAR_full_run_report.md` for stage timings and `scKAR_16cell_validation_report.md`
for the smaller shakedown run and the full list of upstream bugs that had to be fixed.
