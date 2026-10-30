#!/usr/bin/env python
"""Multi-BAM variant of pos_to_gene_htseq_PARALLEL.py.

Identical output format and filenames; the only change is that the 1.5 GB GTF is
parsed ONCE and reused across all cluster BAMs instead of reloaded per array task,
and paths are derived from the BAM location instead of hardcoded.
"""
import HTSeq, csv, os, sys

gtf_file = "/u/scratch/m/mhuang/genomes/GRCh38/gencode.v41.annotation.gtf"
bams = sys.argv[1:]
if not bams:
    sys.exit("usage: pos_to_gene_multi.py <bam> [<bam> ...]")

out_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(bams[0]))), "geneIDs")
os.makedirs(out_dir, exist_ok=True)
print("output dir:", out_dir, flush=True)

print("loading GTF ...", flush=True)
gene_intervals = HTSeq.GenomicArrayOfSets("auto", stranded=False)
n_genes = 0
for feature in HTSeq.GFF_Reader(gtf_file):
    if feature.type == "gene":
        gene_id = feature.attr.get("gene_name")
        if gene_id:
            gene_intervals[feature.iv] += gene_id
            n_genes += 1
print(f"GTF loaded: {n_genes} gene features", flush=True)

for bam_file in bams:
    cluster = os.path.basename(bam_file).replace("_cluster_markers.sorted.bam", "")
    out_csv = os.path.join(out_dir, f"{cluster}_cluster_geneID.csv")
    n_aln = n_with_gene = 0
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["read_id", "chrom", "start", "end", "strand", "gene_id"])
        for aln in HTSeq.BAM_Reader(bam_file):
            if not aln.aligned:
                continue
            n_aln += 1
            iv = aln.iv
            genes = set()
            for _, val in gene_intervals[iv].steps():
                genes.update(val)
            if genes:
                n_with_gene += 1
            w.writerow([aln.read.name, iv.chrom, iv.start, iv.end, iv.strand, ",".join(sorted(genes))])
    print(f"cluster {cluster}: {n_aln} alignments, {n_with_gene} overlapping a gene -> {out_csv}", flush=True)
print("DONE", flush=True)
