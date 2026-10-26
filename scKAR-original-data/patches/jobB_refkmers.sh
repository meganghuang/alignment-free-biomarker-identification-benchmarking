#!/bin/bash
#$ -cwd
#$ -o logs/full_refkmers.log
#$ -j y
#$ -l h_rt=4:00:00,h_data=8G
#$ -pe shared 8
set -u
R=/u/scratch/m/mhuang/claude/scKAR; D=$R/data
export TMPDIR=$R/.tmp; mkdir -p "$TMPDIR" $D/ref
source $R/env.sh
echo "=== refkmers start $(hostname) $(date) ==="
FA=$D/ref/gencode.v50.transcripts.fa.gz
[ -s $FA ] || curl -sSf --retry 4 -o $FA https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v50.transcripts.fa.gz
ls -lh $FA
G=$R/repo/preprocessing_script/sorted_31mers_generator_RT.sh
sed -i 's/-s 100M/-s 2G/' $G
sed -i "s@awk '!/^>/' \$1.fa | sort > \$2@awk '!/^>/' \$1.fa | LC_ALL=C sort -S 20G -T \$TMPDIR > \$2@" $G
cat $G
if [ ! -s $D/ref_31mers.txt ]; then
  bash $G $FA $D/ref_31mers.txt
  echo "generator exit $?"
fi
ls -lh $D/ref_31mers.txt; echo "ref 31-mers: $(wc -l < $D/ref_31mers.txt)"
head -2 $D/ref_31mers.txt; LC_ALL=C sort -c $D/ref_31mers.txt && echo "sorted order OK (C collation)"
rm -f $FA.fa $FA.jf
echo "=== refkmers done $(date) ==="
