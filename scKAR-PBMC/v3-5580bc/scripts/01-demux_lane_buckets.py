#!/usr/bin/env python3
"""01-demux_lane_buckets.py — pass 1 of the full-scale demultiplexer.

At 5,580 cells we cannot hold one output handle per cell (the open-file limit is
typically 1024), so this pass routes each read into one of N_BUCKETS files keyed
by a hash of the cell barcode. Barcode and UMI are written into the read name;
pass 2 splits a bucket into its per-cell FASTQs.

Original 10x lane FASTQs are opened read-only and never modified.
10x 3' GEM-X geometry: R1 = 16 nt cell barcode + 12 nt UMI.
"""
import sys, os, gzip, zlib

lane  = sys.argv[1]
ROOT  = sys.argv[2]
FQ    = "/u/scratch/m/mhuang/scripts/fastq"
BC_LEN, UMI_LEN, N_BUCKETS = 16, 12, 128

want = set()
with open(os.path.join(ROOT, "meta", "cells_full.tsv")) as fh:
    next(fh)
    for line in fh:
        want.add(line.split("\t")[0].encode())

out_dir = os.path.join(ROOT, "data", "buckets")
os.makedirs(out_dir, exist_ok=True)
bh = [gzip.open(f"{out_dir}/{lane}_b{i:03d}.fastq.gz", "wb", compresslevel=1)
      for i in range(N_BUCKETS)]

r1 = open(f"{FQ}/5k_Human_Donor4_PBMC_3p_gem-x_GEX_S4_{lane}_R1_001.fastq", "rb")
r2 = open(f"{FQ}/5k_Human_Donor4_PBMC_3p_gem-x_GEX_S4_{lane}_R2_001.fastq", "rb")

n = matched = 0
while True:
    h1 = r1.readline()
    if not h1:
        break
    s1_ = r1.readline(); r1.readline(); r1.readline()
    h2 = r2.readline(); s2 = r2.readline(); p2 = r2.readline(); q2 = r2.readline()
    n += 1
    bc = s1_[:BC_LEN]
    if bc in want:
        umi = s1_[BC_LEN:BC_LEN + UMI_LEN]
        f = bh[zlib.crc32(bc) % N_BUCKETS]
        f.write(b"@" + bc + b":" + umi + b"\n")
        f.write(s2); f.write(p2); f.write(q2)
        matched += 1
for f in bh:
    f.close()
r1.close(); r2.close()
print(f"{lane}\tpairs={n}\tmatched={matched}\tfrac={matched/max(n,1):.4f}", flush=True)
