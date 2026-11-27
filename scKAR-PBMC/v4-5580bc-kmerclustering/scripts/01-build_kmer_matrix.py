#!/usr/bin/env python3
"""01-build_kmer_matrix.py — prevalence-filtered k-mer x cell matrix for clustering.

scKAR's own KMER_ABUNDANCE mode calls merge_adj_for_clustering, which holds a
DENSE vector<double> of cell_count per k-mer in memory and writes a dense CSV.
At 21,235,492 adjacency k-mers x 5,580 cells that is ~948 GB of RAM and a
multi-terabyte file, so the mode cannot run on droplet-scale data as written.

This builds the same matrix restricted to k-mers present in at least
MIN_CELLS cells. The filter is label-free — it uses only per-k-mer prevalence,
never any clustering or cell annotation — so the clustering that follows remains
independent of the gene-expression pipeline. MIN_CELLS matches the
MIN_ROW_THREDSHOLD the F-test applies downstream (56 = 1% of 5,580).

Output layout matches expression.csv so scKAR's graph-clustering.py reads it
unchanged: header of cell barcodes, one row per k-mer.
"""
import sys, os, glob

F         = sys.argv[1]
MIN_CELLS = int(sys.argv[2])
out       = sys.argv[3]

# cell index -> barcode, in the order the adjacency files index them
# (mother_script.sh sorts tpm_sum.csv, so line order IS index order)
barcodes = []
with open(os.path.join(F, "data", "tpm_sum.csv")) as fh:
    for line in fh:
        barcodes.append(line.split("\t")[0])
n_cells = len(barcodes)
print(f"cells={n_cells}", flush=True)

kept = scanned = 0
with open(out, "w") as o:
    o.write("," + ",".join(barcodes) + "\n")
    for path in sorted(glob.glob(os.path.join(F, "data", "adj", "adj_*.filtered.csv"))):
        with open(path) as fh:
            for line in fh:
                scanned += 1
                parts = line.rstrip("\n").split(",")
                n = (len(parts) - 1) // 2
                if n < MIN_CELLS:
                    continue
                row = ["0"] * n_cells
                it = iter(parts[1:])
                for idx, val in zip(it, it):
                    row[int(idx)] = val
                o.write(parts[0] + "," + ",".join(row) + "\n")
                kept += 1
                if kept % 20000 == 0:
                    print(f"  kept={kept} scanned={scanned}", flush=True)
print(f"scanned={scanned} kept={kept} min_cells={MIN_CELLS}", flush=True)
