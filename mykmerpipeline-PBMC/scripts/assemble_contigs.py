#!/usr/bin/env python
"""Phase 1: assemble cluster-specific marker k-mers into contigs.

Merge rule (from the aim): join k-mers that share a k-1 = 27 bp suffix-prefix overlap
AND whose per-cell count vectors are similar. Unitig walk breaks at branch points, so a
contig is a maximal unambiguous path -- no consensus/greedy extension is attempted.
"""
import csv, sys, os, json, collections
import numpy as np
from scipy import sparse

BASE   = "/u/scratch/m/mhuang/claude"
COUNTS = "/u/scratch/m/mhuang/scripts/bkc_counts"
MARKERS = f"{BASE}/runs/15231823_seurat/plots/cluster_markers.csv"
OUT    = os.environ.get("ASM_OUT", f"{BASE}/contigs")
K      = 28
PADJ   = 0.05
THRESHOLDS = [0.5, 0.7, 0.8, 0.9]
PRIMARY = float(os.environ.get("ASM_R", "0.8"))
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- markers
print("reading markers ...", flush=True)
by_cluster = collections.defaultdict(set)
n_all = n_sig = 0
with open(MARKERS, newline="") as f:
    for row in csv.DictReader(f):
        n_all += 1
        if float(row["p_val_adj"]) < PADJ:
            n_sig += 1
            by_cluster[row["cluster"]].add(row["gene"])
clusters = sorted(by_cluster, key=int)
wanted = set().union(*by_cluster.values())
print(f"marker rows {n_all:,}; significant (padj<{PADJ}) {n_sig:,}; "
      f"unique significant k-mers {len(wanted):,}", flush=True)
for c in clusters:
    print(f"  cluster {c}: {len(by_cluster[c]):,} marker k-mers", flush=True)

# ---------------------------------------------------------------- row indices
print("mapping k-mers to matrix rows ...", flush=True)
row_of = {}
with open(f"{COUNTS}/kmer_ids.txt") as f:
    for i, line in enumerate(f):          # .mtx rows are 1-based in file order
        s = line.strip()
        if s in wanted:
            row_of[s] = i
n_cells = sum(1 for _ in open(f"{COUNTS}/cell_barcodes.txt"))
print(f"located {len(row_of):,} / {len(wanted):,} marker k-mers; {n_cells:,} cells", flush=True)
missing = len(wanted) - len(row_of)
if missing:
    print(f"WARNING: {missing} marker k-mers absent from kmer_ids.txt", flush=True)

# ---------------------------------------------------------------- stream the mtx
print("streaming kmer_matrix.mtx (keeping marker rows only) ...", flush=True)
keep = {v: n for n, v in enumerate(sorted(row_of.values()))}   # matrix row -> compact row
rows, cols, vals = [], [], []
with open(f"{COUNTS}/kmer_matrix.mtx") as f:
    header = f.readline()
    dims = f.readline().split()
    nrow, ncol, nnz = int(dims[0]), int(dims[1]), int(dims[2])
    print(f"  mtx dims {nrow:,} x {ncol:,}, {nnz:,} non-zeros", flush=True)
    for n, line in enumerate(f, 1):
        a = line.split()
        r = int(a[0]) - 1
        c2 = keep.get(r)
        if c2 is not None:
            rows.append(c2); cols.append(int(a[1]) - 1); vals.append(float(a[2]))
        if n % 50_000_000 == 0:
            print(f"  ... {n:,} / {nnz:,} entries, kept {len(vals):,}", flush=True)
M = sparse.csr_matrix((vals, (rows, cols)), shape=(len(keep), n_cells), dtype=np.float32)
del rows, cols, vals
print(f"marker submatrix: {M.shape[0]:,} x {M.shape[1]:,}, {M.nnz:,} non-zeros "
      f"({M.nnz/(M.shape[0]*M.shape[1]):.1%} dense)", flush=True)
sparse.save_npz(f"{OUT}/marker_matrix.npz", M)

# per-row moments for Pearson r
n = float(n_cells)
s1 = np.asarray(M.sum(axis=1)).ravel()
s2 = np.asarray(M.multiply(M).sum(axis=1)).ravel()
mean = s1 / n
var = s2 / n - mean**2
sd = np.sqrt(np.maximum(var, 0))
compact_of = {kmer: keep[row_of[kmer]] for kmer in row_of}

def pearson(i, j):
    if sd[i] <= 0 or sd[j] <= 0:
        return 0.0
    dot = M[i].multiply(M[j]).sum()
    return float((dot / n - mean[i] * mean[j]) / (sd[i] * sd[j]))

# ---------------------------------------------------------------- assembly
def assemble(kmers, r_min):
    """unitig walk over the k-1 overlap graph, gated on count-vector correlation"""
    kmers = [k for k in kmers if k in compact_of]
    by_prefix = collections.defaultdict(list)
    for k in kmers:
        by_prefix[k[:K-1]].append(k)
    succ, pred = {}, collections.defaultdict(list)
    n_cand = n_pass = 0
    for a in kmers:
        outs = []
        for b in by_prefix.get(a[1:], ()):
            if b == a:
                continue
            n_cand += 1
            if pearson(compact_of[a], compact_of[b]) >= r_min:
                n_pass += 1
                outs.append(b)
        if len(outs) == 1:
            succ[a] = outs[0]; pred[outs[0]].append(a)
        elif len(outs) > 1:
            for b in outs:
                pred[b].append(a)           # branch: record, do not extend
    # unitig starts: no unique predecessor, or predecessor branches
    contigs = []
    used = set()
    def is_start(k):
        ps = [p for p in pred.get(k, []) if succ.get(p) == k]
        return len(ps) != 1
    for k in kmers:
        if k in used or not is_start(k):
            continue
        path = [k]; used.add(k)
        cur = k
        while True:
            nxt = succ.get(cur)
            if nxt is None or nxt in used:
                break
            ps = [p for p in pred.get(nxt, []) if succ.get(p) == nxt]
            if len(ps) != 1:
                break
            path.append(nxt); used.add(nxt); cur = nxt
        contigs.append(path)
    for k in kmers:                          # isolated / cycle remnants
        if k not in used:
            contigs.append([k]); used.add(k)
    return contigs, n_cand, n_pass

def seq_of(path):
    s = path[0]
    for nxt in path[1:]:
        s += nxt[-1]
    return s

summary = {}
sweep = {}
for c in clusters:
    kmers = sorted(by_cluster[c])
    sweep[c] = {}
    for r_min in THRESHOLDS:
        cg, nc, npass = assemble(kmers, r_min)
        lens = [len(seq_of(p)) for p in cg]
        sweep[c][r_min] = dict(n_contigs=len(cg), n_multi=sum(1 for p in cg if len(p) > 1),
                               max_len=max(lens), mean_len=round(float(np.mean(lens)), 1),
                               cand_edges=nc, passed=npass)
        print(f"cluster {c} r>={r_min}: {len(cg):,} contigs "
              f"({sum(1 for p in cg if len(p)>1):,} multi-kmer), max {max(lens)} bp, "
              f"edges {npass:,}/{nc:,} passed", flush=True)

    contigs, nc, npass = assemble(kmers, PRIMARY)
    contigs.sort(key=lambda p: -len(p))
    fa = f"{OUT}/cluster{c}_contigs.fasta"
    tb = f"{OUT}/cluster{c}_contigs.tsv"
    with open(fa, "w") as fh, open(tb, "w") as th:
        th.write("contig_id\tlength\tn_kmers\tmean_total_count\tsequence\tmember_kmers\n")
        for i, path in enumerate(contigs):
            s = seq_of(path)
            idx = [compact_of[k] for k in path]
            mt = float(np.mean(s1[idx]))
            cid = f"c{c}_contig{i}_n{len(path)}_L{len(s)}"
            fh.write(f">{cid}\n{s}\n")
            th.write(f"{cid}\t{len(s)}\t{len(path)}\t{mt:.1f}\t{s}\t{','.join(path)}\n")
    lens = [len(seq_of(p)) for p in contigs]
    summary[c] = dict(n_marker_kmers=len(kmers), n_contigs=len(contigs),
                      n_multi_kmer=sum(1 for p in contigs if len(p) > 1),
                      max_len=max(lens), median_len=float(np.median(lens)),
                      n_ge_56bp=sum(1 for L in lens if L >= 2*K),
                      kmers_collapsed=len(kmers) - len(contigs))
    print(f"  -> wrote {fa}", flush=True)

json.dump(dict(primary_r=PRIMARY, padj=PADJ, k=K, summary=summary, sweep=sweep),
          open(f"{OUT}/assembly_summary.json", "w"), indent=1)
print("\n=== ASSEMBLY SUMMARY (r >= %.2f) ===" % PRIMARY, flush=True)
for c in clusters:
    s = summary[c]
    print(f"cluster {c}: {s['n_marker_kmers']:,} k-mers -> {s['n_contigs']:,} contigs "
          f"({s['n_multi_kmer']:,} multi-k-mer, {s['n_ge_56bp']:,} >=56bp), "
          f"max {s['max_len']} bp, collapsed {s['kmers_collapsed']:,}", flush=True)
print("ASSEMBLY_DONE", flush=True)
