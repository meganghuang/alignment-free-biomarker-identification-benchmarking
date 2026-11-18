#!/usr/bin/env python3
"""02-split_dedup_bucket.py — pass 2: per-cell FASTQs with UMI collapsing.

Reads one bucket across all four lanes, groups records by cell barcode, keeps a
single read per distinct UMI (first occurrence wins), and writes
data/reads/<barcode>.fastq.gz plus a per-cell before/after count table.

scKAR as published has no UMI model; collapsing is our addition.
"""
import sys, os, gzip
from collections import defaultdict

bucket = int(sys.argv[1])
ROOT   = sys.argv[2]
LANES  = ["L001", "L002", "L003", "L004"]
src    = os.path.join(ROOT, "data", "buckets")
dst    = os.path.join(ROOT, "data", "reads")
os.makedirs(dst, exist_ok=True)

seen  = defaultdict(set)   # barcode -> set(umi)
recs  = defaultdict(list)
total = defaultdict(int)

for lane in LANES:
    p = f"{src}/{lane}_b{bucket:03d}.fastq.gz"
    if not os.path.exists(p):
        continue
    with gzip.open(p, "rb") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline(); pl = fh.readline(); q = fh.readline()
            bc, umi = h[1:].rstrip(b"\n").split(b":", 1)
            total[bc] += 1
            if umi in seen[bc]:
                continue
            seen[bc].add(umi)
            recs[bc].append((h, s, pl, q))

rep = open(os.path.join(ROOT, f"dedup_counts.b{bucket:03d}.tsv"), "w")
for bc, rs in recs.items():
    name = bc.decode()
    with gzip.open(os.path.join(dst, f"{name}.fastq.gz"), "wb", compresslevel=6) as o:
        for h, s, pl, q in rs:
            o.write(h); o.write(s); o.write(pl); o.write(q)
    rep.write(f"{name}\t{total[bc]}\t{len(rs)}\t{1 - len(rs)/total[bc]:.4f}\n")
rep.close()
print(f"bucket {bucket} cells={len(recs)}", flush=True)
