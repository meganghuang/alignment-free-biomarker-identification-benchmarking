#!/usr/bin/env python3
"""01-demux_lane_umi.py — route 10x reads to per-cell FASTQ, carrying the UMI.

Reads ONE lane of the original 10x FASTQs (read-only; never modified) and writes
one gzipped FASTQ per wanted barcode containing the R2 (cDNA) records. The 12 nt
UMI from R1 is appended to each read name as "UMI:<umi>" so the next step can
collapse PCR duplicates without a second pass over the lanes.

10x 3' GEM-X geometry: R1 = 16 nt cell barcode + 12 nt UMI.
"""
import sys, os, gzip

lane = sys.argv[1]
FQ   = "/u/scratch/m/mhuang/scripts/fastq"
ROOT = "/u/scratch/m/mhuang/claude/compare3/sckar_percell"
BC_LEN, UMI_LEN = 16, 12

want = {}
with open(os.path.join(ROOT, "meta", "cells.tsv")) as fh:
    next(fh)
    for line in fh:
        cell, bc, _, _ = line.rstrip("\n").split("\t")
        want[cell.encode()] = cell

out = os.path.join(ROOT, "data", "lanes_umi", lane)
os.makedirs(out, exist_ok=True)
fh_out = {k: gzip.open(f"{out}/{v}.fastq.gz", "wb", compresslevel=1) for k, v in want.items()}

r1 = open(f"{FQ}/5k_Human_Donor4_PBMC_3p_gem-x_GEX_S4_{lane}_R1_001.fastq", "rb")
r2 = open(f"{FQ}/5k_Human_Donor4_PBMC_3p_gem-x_GEX_S4_{lane}_R2_001.fastq", "rb")

n = matched = 0
while True:
    h1 = r1.readline()
    if not h1:
        break
    s1 = r1.readline(); r1.readline(); r1.readline()
    h2 = r2.readline(); s2 = r2.readline(); p2 = r2.readline(); q2 = r2.readline()
    n += 1
    f = fh_out.get(s1[:BC_LEN])
    if f is not None:
        umi = s1[BC_LEN:BC_LEN + UMI_LEN]
        f.write(h2.rstrip(b"\n") + b" UMI:" + umi + b"\n")
        f.write(s2); f.write(p2); f.write(q2)
        matched += 1
for f in fh_out.values():
    f.close()
r1.close(); r2.close()
print(f"{lane}\tpairs={n}\tmatched={matched}\tfrac={matched/max(n,1):.4f}", flush=True)
