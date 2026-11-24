#!/usr/bin/env python
"""Parse BLAT PSL -> best placement per contig -> exonic/intronic/intergenic.

Precedence (stated): a contig is EXONIC if any aligned block overlaps a GENCODE exon;
else INTRONIC if any block overlaps an annotated intron (or falls inside a gene body);
else INTERGENIC. Multi-mapping is reported, never hidden: a contig whose second-best
alignment scores within 5% of the best is flagged ambiguous.
"""
import sys, os, csv, collections, subprocess, json

PSL  = sys.argv[1]
OUT  = sys.argv[2]
GDIR = "/u/scratch/m/mhuang/genomes/GRCh38"
BT   = "/u/home/m/mhuang/bedtools2/bin/bedtools"

rows = []
with open(PSL) as f:
    for line in f:
        a = line.rstrip("\n").split("\t")
        if len(a) < 21 or not a[0].isdigit():
            continue
        rows.append(dict(matches=int(a[0]), mismatches=int(a[1]), strand=a[8],
                         qname=a[9], qsize=int(a[10]), qstart=int(a[11]), qend=int(a[12]),
                         tname=a[13], tstart=int(a[15]), tend=int(a[16]),
                         blocks=int(a[17]),
                         blockSizes=[int(x) for x in a[18].rstrip(",").split(",") if x],
                         tStarts=[int(x) for x in a[20].rstrip(",").split(",") if x]))
print(f"PSL alignments parsed: {len(rows):,}", flush=True)

by_q = collections.defaultdict(list)
for r in rows:
    by_q[r["qname"]].append(r)
print(f"contigs with >=1 alignment: {len(by_q):,}", flush=True)

best, ambiguous = {}, set()
for q, al in by_q.items():
    al.sort(key=lambda r: -(r["matches"] - r["mismatches"]))
    b = al[0]
    cov = (b["qend"] - b["qstart"]) / b["qsize"]
    if cov < 0.90:
        continue
    best[q] = b
    s0 = b["matches"] - b["mismatches"]
    n_near = sum(1 for r in al if (r["matches"] - r["mismatches"]) >= 0.95 * s0)
    b["n_alignments"] = len(al)
    b["n_near_best"] = n_near
    if n_near > 1:
        ambiguous.add(q)
print(f"contigs placed (>=90% query covered): {len(best):,}; ambiguous (multi-mapping): {len(ambiguous):,}",
      flush=True)
print(f"spliced placements (>1 block): {sum(1 for b in best.values() if b['blocks'] > 1):,}", flush=True)

os.makedirs(OUT, exist_ok=True)
blocks_bed = f"{OUT}/contig_blocks.bed"
with open(blocks_bed, "w") as fh:
    for q, b in sorted(best.items()):
        for sz, st in zip(b["blockSizes"], b["tStarts"]):
            fh.write(f"{b['tname']}\t{st}\t{st+sz}\t{q}\t{sz}\t{b['strand'][0]}\n")
subprocess.run(f"sort -k1,1 -k2,2n {blocks_bed} -o {blocks_bed}", shell=True, check=True)

def overlap_set(ref):
    out = collections.defaultdict(set)
    p = subprocess.run([BT, "intersect", "-a", blocks_bed, "-b", f"{GDIR}/{ref}", "-wa", "-wb"],
                       capture_output=True, text=True)
    for line in p.stdout.splitlines():
        a = line.split("\t")
        out[a[3]].add(a[9] if len(a) > 9 else "")
    return out

ex  = overlap_set("exons.bed")
intr = overlap_set("introns.bed")
gen = overlap_set("genes.bed")
print(f"block overlaps -> exonic {len(ex):,} | intronic {len(intr):,} | within-gene {len(gen):,}", flush=True)

cat_counts = collections.Counter()
with open(f"{OUT}/contig_annotation.tsv", "w") as fh:
    fh.write("contig_id\tcluster\tlength\tn_kmers\tchrom\tstart\tend\tstrand\tblocks\t"
             "n_alignments\tn_near_best\tambiguous\tcategory\tgenes\n")
    for q, b in sorted(best.items()):
        cl = q.split("_")[0].lstrip("c")
        nk = q.split("_n")[1].split("_")[0]
        if q in ex:      cat = "exonic"
        elif q in intr:  cat = "intronic"
        elif q in gen:   cat = "intronic_genebody"
        else:            cat = "intergenic"
        genes = sorted((ex.get(q) or intr.get(q) or gen.get(q) or set()) - {""})
        cat_counts[(cl, cat)] += 1
        fh.write(f"{q}\t{cl}\t{b['qsize']}\t{nk}\t{b['tname']}\t{b['tstart']}\t{b['tend']}\t"
                 f"{b['strand'][0]}\t{b['blocks']}\t{b['n_alignments']}\t{b['n_near_best']}\t"
                 f"{int(q in ambiguous)}\t{cat}\t{','.join(genes[:8])}\n")

# unplaced contigs matter as much as placed ones
all_ids = []
for cl in range(8):
    fa = f"/u/scratch/m/mhuang/claude/contigs/cluster{cl}_contigs.fasta"
    if os.path.exists(fa):
        all_ids += [l[1:].strip() for l in open(fa) if l.startswith(">")]
unplaced = [q for q in all_ids if q not in best]
with open(f"{OUT}/contigs_unplaced.txt", "w") as fh:
    fh.write("\n".join(unplaced) + "\n")
print(f"contigs total {len(all_ids):,}; placed {len(best):,}; unplaced {len(unplaced):,}", flush=True)

print("\n=== category by cluster ===", flush=True)
cats = ["exonic", "intronic", "intronic_genebody", "intergenic"]
print("cluster  " + "  ".join(f"{c:>18}" for c in cats))
for cl in [str(i) for i in range(8)]:
    tot = sum(cat_counts[(cl, c)] for c in cats)
    if tot:
        print(f"{cl:>7}  " + "  ".join(f"{cat_counts[(cl,c)]:>18}" for c in cats) + f"   (n={tot})")
json.dump({f"{k[0]}|{k[1]}": v for k, v in cat_counts.items()},
          open(f"{OUT}/category_counts.json", "w"), indent=1)
print("ANNOTATION_DONE", flush=True)
