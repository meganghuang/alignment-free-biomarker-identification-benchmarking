#!/usr/bin/env python
"""Summarize HTSeq gene assignments per cluster.

Counts DISTINCT marker k-mers per gene (not alignments: one k-mer can contribute
thousands of alignments when it multi-maps, which would dominate any naive count).
Reports the full ranking plus the subset of k-mers with few alignments, where a
gene assignment is actually trustworthy.
"""
import csv, glob, os, sys, collections, json

gdir = "/u/scratch/m/mhuang/claude/runs/15231823_seurat/plots/geneIDs"
out_json = os.path.join(gdir, "gene_summary.json")
UNIQ_MAX = 10          # k-mers aligning to <=10 loci count as "well-placed"
summary = {}

for path in sorted(glob.glob(os.path.join(gdir, "*_cluster_geneID.csv"))):
    cl = os.path.basename(path).split("_")[0]
    n_aln = 0
    kmer_nalign = collections.Counter()
    kmer_genes = collections.defaultdict(set)
    with open(path, newline="") as f:
        r = csv.reader(f)
        next(r, None)
        for row in r:
            if len(row) < 6:
                continue
            kmer, genes = row[0], row[5]
            n_aln += 1
            kmer_nalign[kmer] += 1
            if genes:
                kmer_genes[kmer].update(genes.split(","))

    gene_kmers = collections.Counter()
    gene_kmers_uniq = collections.Counter()
    for kmer, genes in kmer_genes.items():
        well_placed = kmer_nalign[kmer] <= UNIQ_MAX
        for g in genes:
            gene_kmers[g] += 1
            if well_placed:
                gene_kmers_uniq[g] += 1

    n_kmers = len(kmer_nalign)
    n_well = sum(1 for k, v in kmer_nalign.items() if v <= UNIQ_MAX)
    summary[cl] = dict(
        alignments=n_aln,
        kmers_aligned=n_kmers,
        kmers_well_placed=n_well,
        median_alignments_per_kmer=sorted(kmer_nalign.values())[n_kmers // 2] if n_kmers else 0,
        genes_hit=len(gene_kmers),
        top_genes=gene_kmers.most_common(15),
        top_genes_well_placed=gene_kmers_uniq.most_common(15),
    )
    s = summary[cl]
    print(f"\n=== cluster {cl} ===", flush=True)
    print(f"  alignments {s['alignments']:,} | k-mers {s['kmers_aligned']:,} | "
          f"well-placed (<= {UNIQ_MAX} loci) {s['kmers_well_placed']:,} "
          f"({100*s['kmers_well_placed']/max(1,s['kmers_aligned']):.1f}%) | "
          f"median loci/k-mer {s['median_alignments_per_kmer']} | genes {s['genes_hit']:,}")
    print("  top genes (all k-mers)   :", ", ".join(f"{g}({n})" for g, n in s["top_genes"][:10]))
    print("  top genes (well-placed)  :", ", ".join(f"{g}({n})" for g, n in s["top_genes_well_placed"][:10]))

with open(out_json, "w") as f:
    json.dump(summary, f, indent=1)
print("\nwrote", out_json, flush=True)
