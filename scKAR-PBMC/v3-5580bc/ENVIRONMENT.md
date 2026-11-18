# Environment — scKAR PBMC v3 (5,580 cells)

Run on UCLA Hoffman2 (OARC Slurm cluster), compute nodes.

| component | version |
|---|---|
| OS | Rocky Linux 9.8 (Blue Onyx) |
| scheduler | Slurm 24.05.0 |
| gcc | 11.5.0 20240719 (Red Hat 11.5.0-14) |
| jellyfish | 2.3.1 |
| ABySS | 2.2.3 (`abyss-pe`, run with `np=1`) |
| Python | 3.12.9 — numpy 2.3.2, scanpy 1.12 |
| R | 4.1.3 (2022-03-10) |
| DESeq2 | 1.34.0 |
| scran | 1.22.1 |
| glmGamPoi | 1.6.0 |

Conda environments (`/u/scratch/m/mhuang/claude/scKAR/envs/`): `sckar_py`,
`sckar_r`, `sckar_tools`. No `module load` is needed on compute nodes; system
gcc is new enough to build scKAR's C++ once the static-link flags are removed.

Reference for k-mer subtraction: GENCODE v50 transcriptome 31-mers,
215,142,544 sorted entries (`data/ref_31mers.txt`, 6.9 GB — not in this repo).

## Job shapes

| stage | resources | wall time |
|---|---|---|
| demultiplex + UMI collapse | 16 cores, 8 GB/core | ~2 h |
| jellyfish + reference subtraction + adjacency | 16 cores, 12 GB/core | 7 h 43 m |
| clustering + F-test | 16 cores, 12 GB/core | ~3 min |
| DESeq2 (10 bipartitions × 6 chunks) | 16 cores, 12 GB/core | 2 h 53 m |
| ABySS assembly | 16 cores, 12 GB/core | ~1 min |

Slurm jobs route by account only (`--account=pi_mcdb`, no `--partition`), which
lands on `pi_shared` with `qos=pi_shared24` and a 24-hour walltime ceiling.
