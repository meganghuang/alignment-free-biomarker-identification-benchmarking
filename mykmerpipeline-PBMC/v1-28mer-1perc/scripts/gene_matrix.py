#!/usr/bin/env python
"""Full gene x cluster table of DISTINCT marker k-mers (all, and well-placed <=10 loci)."""
import csv, glob, os, collections

gdir = "/u/scratch/m/mhuang/claude/runs/15231823_seurat/plots/geneIDs"
UNIQ_MAX = 10
all_counts, uniq_counts, clusters = {}, {}, []

for path in sorted(glob.glob(os.path.join(gdir, "*_cluster_geneID.csv"))):
    cl = os.path.basename(path).split("_")[0]
    clusters.append(cl)
    kmer_n = collections.Counter()
    kmer_genes = collections.defaultdict(set)
    with open(path, newline="") as f:
        r = csv.reader(f); next(r, None)
        for row in r:
            if len(row) < 6: continue
            kmer_n[row[0]] += 1
            if row[5]: kmer_genes[row[0]].update(row[5].split(","))
    a, u = collections.Counter(), collections.Counter()
    for kmer, genes in kmer_genes.items():
        wp = kmer_n[kmer] <= UNIQ_MAX
        for g in genes:
            a[g] += 1
            if wp: u[g] += 1
    all_counts[cl], uniq_counts[cl] = a, u
    print(f"cluster {cl}: {len(a)} genes", flush=True)

genes = sorted(set().union(*[set(c) for c in all_counts.values()]))
for name, data in (("gene_by_cluster_all.csv", all_counts), ("gene_by_cluster_wellplaced.csv", uniq_counts)):
    with open(os.path.join(gdir, name), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["gene"] + [f"c{c}" for c in clusters])
        for g in genes:
            w.writerow([g] + [data[c].get(g, 0) for c in clusters])
    print("wrote", name, flush=True)
print("genes total:", len(genes), flush=True)
