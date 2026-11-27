# Environment — scKAR PBMC v4 (k-mer clustering)

Identical to `../v3-5580bc/ENVIRONMENT.md` — Rocky Linux 9.8, Slurm 24.05.0,
gcc 11.5.0, jellyfish 2.3.1, ABySS 2.2.3, Python 3.12.9 (numpy 2.3.2,
scanpy 1.12), R 4.1.3 (DESeq2 1.34.0, scran 1.22.1, glmGamPoi 1.6.0),
BLAT v39x1, GENCODE v41.

## Job shapes specific to v4

| stage | resources | wall time |
|---|---|---|
| k-mer matrix build + Leiden clustering | 16 cores, 16 GB/core | ~3 h |
| F-test | same | ~4 min |
| DESeq2 (13 bipartitions x 6 chunks) | same | 6 h 24 m |
| ABySS assembly | same | ~1 min |
| contig annotation (BLAT + GENCODE) | 1 core, 48 GB | ~1 min |

Stage 1 (jellyfish counting, reference subtraction, adjacency) was **not**
re-run: those products are independent of clustering mode and were copied from
`../v3-5580bc/`. They were copied rather than symlinked because `run_stages.sh`
runs `rm -f $D/adj/adj_*.csv` whenever the stage sentinel is absent, which
through a symlink would have destroyed the v3 originals.
