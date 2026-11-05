# Environment

Host: UCLA Hoffman2 (job 15370899, 4 slots, h_data=12G, `pod_smp.q@n1215`).
Started 2026-11-02 12:36:19 PST, completed 17:22:47 PST — 4 h 46 min, single attempt.

Tooling (micromamba envs under `/u/scratch/m/mhuang/claude/scKAR/envs`, reused):

- `sckar_py` — Python 3.12.9, scanpy 1.12, numpy 2.3.2
- `sckar_r` — R 4.1.3, DESeq2 1.34.0, scran 1.22.1, glmGamPoi 1.6.0, harmonicmeanp 3.0
- `sckar_tools` — jellyfish 2.3.1 (bioconda `kmer-jellyfish`), ABySS 2.2.3
- blat at `~/userApps/bin/blat`; gcc/9.5.0 for the C++ sources

scKAR source: github.com/RedwanulKarim612/scKAR @ 088c3092, with the local
patches carried over from the RCC run (C++17 flags, `ref_filter.cpp`
out-of-bounds guard, `tpm_sum.csv` byte-order sort, `blat.sh` exec bit,
`mother_script.sh` wait-for-background-batch). The patched clone used here is a
separate copy at `compare3/sckar_pbmc/repo` so the RCC results stay reproducible.

Stage wall times: jellyfish counting of 52 pseudo-cells 1,781 s; filtering,
adjacency (16 chunks) and GENCODE reference subtraction to ~15:00; F-test 51
chunks; DESeq2 50 tables; assembly + BLAT + harmonic-mean p-values to 17:22.
