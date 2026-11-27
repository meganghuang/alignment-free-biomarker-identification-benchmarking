#!/usr/bin/env python3
"""05-collect_contigs.py — merge every bipartition's contigs into one FASTA.

BLAT reloads the whole 3.2 Gb genome on each invocation, so all contigs are
concatenated into a single query file and mapped in one pass. Each sequence is
renamed <bipartition>|<arm>|<original_id> so the hits can be split apart again.
"""
import os, sys, glob

F   = sys.argv[1]
out = sys.argv[2]
n = 0
with open(out, "w") as o:
    for bip in range(13):
        for arm in ("A", "B"):
            p = f"{F}/data/final_results/bipartition_{bip}/abyss/{arm}_contigs.fasta"
            if not os.path.exists(p) or os.path.getsize(p) == 0:
                continue
            with open(p) as fh:
                for line in fh:
                    if line.startswith(">"):
                        n += 1
                        o.write(f">bip{bip}|{arm}|{line[1:].strip().split()[0]}\n")
                    else:
                        o.write(line)
print(f"contigs_written={n}", flush=True)
