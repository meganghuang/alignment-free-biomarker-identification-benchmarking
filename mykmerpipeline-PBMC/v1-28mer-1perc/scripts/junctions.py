#!/usr/bin/env python
"""Phase 3a: junctions from spliced contig placements, classified against GENCODE v41."""
import sys, os, csv, collections, subprocess, json

BASE = "/u/scratch/m/mhuang/claude"
GDIR = "/u/scratch/m/mhuang/genomes/GRCh38"
OUT  = f"{BASE}/junctions"
os.makedirs(OUT, exist_ok=True)
MINBLK = 8     # ignore micro-blocks that are likely alignment noise

# ---- contig-defined junctions from the PSL -------------------------------
best = {}
with open(f"{BASE}/placement/all_contigs.psl") as f:
    for line in f:
        a = line.rstrip("\n").split("\t")
        if len(a) < 21 or not a[0].isdigit():
            continue
        q, score = a[9], int(a[0]) - int(a[1])
        if q not in best or score > best[q][0]:
            best[q] = (score, a)

junc = []
n_spliced = 0
for q, (score, a) in best.items():
    nb = int(a[17])
    if nb < 2:
        continue
    sizes = [int(x) for x in a[18].rstrip(",").split(",") if x]
    starts = [int(x) for x in a[20].rstrip(",").split(",") if x]
    chrom, strand = a[13], a[8][0]
    spliced_here = False
    for i in range(nb - 1):
        if sizes[i] < MINBLK or sizes[i+1] < MINBLK:
            continue
        don = starts[i] + sizes[i]          # 0-based end of left block
        acc = starts[i+1]                   # 0-based start of right block
        if acc - don < 20:                  # gap too small to be an intron
            continue
        junc.append((chrom, don, acc, q, strand, sizes[i], sizes[i+1]))
        spliced_here = True
    n_spliced += spliced_here
print(f"contigs with a spliced best-placement: {n_spliced:,}", flush=True)
print(f"candidate junctions (both flanks >={MINBLK} bp, gap >=20 bp): {len(junc):,}", flush=True)
if not junc:
    print("NO_JUNCTIONS -- nothing further to classify", flush=True)

# ---- annotated junction set from GENCODE ---------------------------------
print("building the GENCODE v41 junction set ...", flush=True)
tx = collections.defaultdict(list)
with open(f"{GDIR}/gencode.v41.annotation.gtf") as f:
    for line in f:
        if line.startswith("#"):
            continue
        a = line.split("\t")
        if a[2] != "exon":
            continue
        tid = a[8].split('transcript_id "')[1].split('"')[0]
        tx[(a[0], tid, a[6])].append((int(a[3]) - 1, int(a[4])))
annotated = set()
for (chrom, tid, strand), exons in tx.items():
    exons.sort()
    for (s1, e1), (s2, e2) in zip(exons, exons[1:]):
        if s2 > e1:
            annotated.add((chrom, e1, s2))
print(f"annotated junctions in GENCODE v41: {len(annotated):,}", flush=True)

# ---- flanking dinucleotides ----------------------------------------------
fa = f"{GDIR}/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
SAM = "/u/home/m/mhuang/local/bin/samtools"
def dinucs(batch):
    regions = []
    for chrom, don, acc in batch:
        regions += [f"{chrom}:{don+1}-{don+2}", f"{chrom}:{acc-1}-{acc}"]
    if not regions:
        return {}
    p = subprocess.run([SAM, "faidx", fa] + regions, capture_output=True, text=True)
    seqs, cur = [], None
    for line in p.stdout.splitlines():
        if line.startswith(">"):
            if cur is not None: seqs.append(cur)
            cur = ""
        else:
            cur += line.strip()
    if cur is not None: seqs.append(cur)
    out = {}
    for i, key in enumerate(batch):
        if 2*i + 1 < len(seqs):
            out[key] = (seqs[2*i].upper(), seqs[2*i+1].upper())
    return out

uniq = sorted({(c, d, a) for c, d, a, *_ in junc})
dn = {}
for i in range(0, len(uniq), 200):
    dn.update(dinucs(uniq[i:i+200]))

with open(f"{OUT}/contig_junctions.tsv", "w") as fh:
    fh.write("chrom\tdonor_end\tacceptor_start\tcontig_id\tcluster\tstrand\t"
             "left_block\tright_block\tintron_len\tstatus\tdonor_dinuc\tacceptor_dinuc\tmotif\n")
    stat = collections.Counter()
    for chrom, don, acc, q, strand, lb, rb in sorted(junc):
        status = "annotated" if (chrom, don, acc) in annotated else "novel"
        d, a2 = dn.get((chrom, don, acc), ("NA", "NA"))
        motif = f"{d}-{a2}"
        canon = motif in ("GT-AG", "CT-AC", "GC-AG", "CT-GC", "AT-AC", "GT-AT")
        stat[(q.split("_")[0], status, canon)] += 1
        fh.write(f"{chrom}\t{don}\t{acc}\t{q}\t{q.split('_')[0].lstrip('c')}\t{strand}\t"
                 f"{lb}\t{rb}\t{acc-don}\t{status}\t{d}\t{a2}\t{motif}\n")
print("\n=== junction classification (cluster, status, canonical motif) ===", flush=True)
for k in sorted(stat):
    print(f"  {k[0]:>4}  {k[1]:>9}  canonical={k[2]}  n={stat[k]}", flush=True)
json.dump({"|".join(map(str, k)): v for k, v in stat.items()},
          open(f"{OUT}/junction_counts.json", "w"), indent=1)
print("JUNCTIONS_DONE", flush=True)
