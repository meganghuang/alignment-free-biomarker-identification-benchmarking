
import sys, gzip, os
lane, mapfile, outdir, fqdir = sys.argv[1:5]
bc2pc = {}
with open(mapfile,'rb') as fh:
    for line in fh:
        bc, pc = line.rstrip().split(b'\t')
        bc2pc[bc] = pc.decode()
pcs = sorted(set(bc2pc.values()))
os.makedirs(outdir, exist_ok=True)
H = {pc: gzip.open(os.path.join(outdir, pc + '.fastq.gz'), 'wb', compresslevel=1) for pc in pcs}
base = os.path.join(fqdir, "5k_Human_Donor4_PBMC_3p_gem-x_GEX_S4_%s" % lane)
r1 = open(base + "_R1_001.fastq", 'rb', buffering=1<<23)
r2 = open(base + "_R2_001.fastq", 'rb', buffering=1<<23)
tot = matched = 0
get = bc2pc.get
rl1, rl2 = r1.readline, r2.readline
while True:
    h1 = rl1()
    if not h1: break
    s1 = rl1(); rl1(); rl1()
    h2 = rl2(); s2 = rl2(); p2 = rl2(); q2 = rl2()
    tot += 1
    pc = get(s1[:16])
    if pc is not None:
        w = H[pc].write
        w(h2); w(s2); w(p2); w(q2)
        matched += 1
for h in H.values(): h.close()
with open(os.path.join(outdir, 'stats.txt'), 'w') as f:
    f.write("%s\ttotal\t%d\tmatched\t%d\n" % (lane, tot, matched))
print(lane, "total", tot, "matched", matched, flush=True)
