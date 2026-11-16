#!/usr/bin/env python3
"""02-dedup_umi.py — collapse PCR duplicates per cell and write scKAR's input reads.

For each wanted barcode, merge that cell's records from all four lanes and keep
ONE read per distinct UMI (first occurrence wins). Writes data/reads/<cell>.fastq.gz,
the per-cell FASTQ that scKAR's Methods specify as input, and a TSV of
before/after counts so the duplicate rate is a reported number.

scKAR as published has no UMI model; this collapsing step is our addition and is
documented as a deviation in the arm's README.
"""
import os, gzip, sys

ROOT  = "/u/scratch/m/mhuang/claude/compare3/sckar_percell"
LANES = ["L001", "L002", "L003", "L004"]
src   = os.path.join(ROOT, "data", "lanes_umi")
dst   = os.path.join(ROOT, "data", "reads")
os.makedirs(dst, exist_ok=True)

cells = []
with open(os.path.join(ROOT, "meta", "cells.tsv")) as fh:
    next(fh)
    for line in fh:
        cells.append(line.split("\t")[0])

shard, nshard = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (0, 1)
rep = open(os.path.join(ROOT, f"dedup_counts.{shard}.tsv"), "w")
rep.write("cell\treads_in\treads_out\tdup_frac\n")

for i, cell in enumerate(cells):
    if i % nshard != shard:
        continue
    seen, kept, total = set(), [], 0
    for lane in LANES:
        p = os.path.join(src, lane, f"{cell}.fastq.gz")
        if not os.path.exists(p):
            continue
        with gzip.open(p, "rb") as fh:
            while True:
                h = fh.readline()
                if not h:
                    break
                s = fh.readline(); pl = fh.readline(); q = fh.readline()
                total += 1
                umi = h.rsplit(b"UMI:", 1)[-1].strip() if b"UMI:" in h else b""
                if umi in seen:
                    continue
                seen.add(umi)
                kept.append((h, s, pl, q))
    with gzip.open(os.path.join(dst, f"{cell}.fastq.gz"), "wb", compresslevel=6) as o:
        for h, s, pl, q in kept:
            o.write(h); o.write(s); o.write(pl); o.write(q)
    dup = 1 - (len(kept) / total) if total else 0.0
    rep.write(f"{cell}\t{total}\t{len(kept)}\t{dup:.4f}\n")
    rep.flush()
rep.close()
print(f"shard {shard}/{nshard} done", flush=True)
