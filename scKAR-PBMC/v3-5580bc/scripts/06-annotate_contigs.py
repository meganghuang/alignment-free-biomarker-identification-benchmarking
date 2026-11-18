#!/usr/bin/env python3
"""06-annotate_contigs.py — assign genes to BLAT hits.

Reads the combined PSL, keeps the best-scoring alignment per contig (and counts
how many loci each contig hit), then assigns overlapping GENCODE genes by
interval overlap. gffutils is installed, but building a database over the 1.5 GB
GTF is far slower than indexing its `gene` features directly, and the answer is
the same.

Outputs:
  contig_annotation.tsv  one row per contig: best locus, multi-mapping count, genes
  annotation_summary.tsv per-bipartition/arm placement and annotation counts
"""
import sys, bisect
from collections import defaultdict

PSL, GTF, OUT_ANN, OUT_SUM = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

# ---- PSL -> best hit per contig -------------------------------------------
best, nloci = {}, defaultdict(int)
with open(PSL) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 21 or not f[0].isdigit():
            continue                      # header / blank
        match, mismatch = int(f[0]), int(f[1])
        qname, strand   = f[9], f[8]
        tname, tstart, tend = f[13], int(f[15]), int(f[16])
        score = match - mismatch
        nloci[qname] += 1
        if qname not in best or score > best[qname][0]:
            best[qname] = (score, tname, tstart, tend, strand, int(f[10]))
print(f"contigs_with_hits={len(best)}", flush=True)

# ---- GENCODE gene intervals ------------------------------------------------
genes = defaultdict(list)
with open(GTF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t", 9)
        if f[2] != "gene":
            continue
        attr = f[8]
        def grab(key):
            i = attr.find(key + ' "')
            if i < 0:
                return ""
            j = attr.find('"', i + len(key) + 2)
            return attr[i + len(key) + 2:j]
        genes[f[0]].append((int(f[3]) - 1, int(f[4]), grab("gene_name"), grab("gene_type")))
for c in genes:
    genes[c].sort()
starts  = {c: [g[0] for g in v] for c, v in genes.items()}
# longest gene on each chromosome bounds how far back an overlap can start
maxlen  = {c: max(ge - gs for gs, ge, _, _ in v) for c, v in genes.items()}
print(f"gtf_chroms={len(genes)} genes={sum(len(v) for v in genes.values())} "
      f"longest_gene={max(maxlen.values())}", flush=True)

def overlap(chrom, s, e):
    """Genes overlapping [s, e). Exact: scan back only as far as the longest
    gene on this chromosome, so no long gene starting upstream is missed."""
    v = genes.get(chrom)
    if not v:
        return []
    hi  = bisect.bisect_right(starts[chrom], e)          # genes starting before e
    lo  = bisect.bisect_left(starts[chrom], s - maxlen[chrom])
    return [(name, gtype) for gs, ge, name, gtype in v[lo:hi] if ge > s]

# ---- write -----------------------------------------------------------------
summary = defaultdict(lambda: {"contigs": 0, "placed": 0, "unique_locus": 0, "annotated": 0})
with open(OUT_ANN, "w") as o:
    o.write("contig\tbipartition\tarm\tscore\tchrom\tstart\tend\tstrand\tn_loci\tgenes\tgene_types\n")
    for qname, (score, tname, ts, te, strand, qsize) in sorted(best.items()):
        bip, arm, cid = qname.split("|", 2)
        hits = overlap(tname, ts, te)
        names = ",".join(sorted({h[0] for h in hits if h[0]})) or "."
        types = ",".join(sorted({h[1] for h in hits if h[1]})) or "."
        o.write(f"{cid}\t{bip}\t{arm}\t{score}\t{tname}\t{ts}\t{te}\t{strand}\t{nloci[qname]}\t{names}\t{types}\n")
        k = (bip, arm)
        summary[k]["placed"] += 1
        if nloci[qname] == 1:
            summary[k]["unique_locus"] += 1
        if names != ".":
            summary[k]["annotated"] += 1

with open(OUT_SUM, "w") as o:
    o.write("bipartition\tarm\tplaced\tunique_locus\tannotated\n")
    for (bip, arm), d in sorted(summary.items()):
        o.write(f"{bip}\t{arm}\t{d['placed']}\t{d['unique_locus']}\t{d['annotated']}\n")
print("done", flush=True)
