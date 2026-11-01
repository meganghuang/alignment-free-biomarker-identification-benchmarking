# Environment

Host: UCLA Hoffman2, CentOS 7.9 compute node (job 15370678, 1 slot, h_data=4G).

R stack: conda env at `/u/scratch/m/mhuang/claude/envs/seurat5_r44`

- R 4.4.3
- Seurat 5.5.1, SeuratObject 5.4.0
- Matrix 1.7.6
- hdf5r 1.3.16 (required by `Read10X_h5`)
- ggplot2, patchwork, plyr

Run with `R_LIBS_USER=""` so the personal library does not shadow the env's
SeuratObject. Full `sessionInfo()` is at the end of `logs/run_log.txt`.

The original analysis was authored against Seurat 4; this rerun used Seurat 5.
Because cell-type labels are joined per barcode rather than re-derived from the
new clustering, the labels are unaffected by the version change, and the
cluster sizes reproduce the original run exactly.

Input data (not in this repo):
`/u/scratch/m/mhuang/std_scripts/counts/5k_Human_Donor4_PBMC_3p_gem-x_*_filtered_feature_bc_matrix.h5`
