#!/bin/bash
#$ -cwd
#$ -o logs/full_download.log
#$ -j y
#$ -l h_rt=4:00:00,h_data=4G
#$ -pe shared 8
set -u
R=/u/scratch/m/mhuang/claude/scKAR; D=$R/data
mkdir -p $D/reads $R/stage/fq
echo "=== download start $(hostname) $(date) ==="
tr -d '\r' < $R/stage/manifest_full.tsv > $R/stage/manifest_full.clean.tsv
dl(){ cell=$1; u1=$2; u2=$3; out=$D/reads/$cell.fastq.gz
  if [ -s "$out" ]; then return; fi
  curl -sSf --retry 4 --retry-delay 5 -o $R/stage/fq/${cell}_1.gz "https://$u1" &&
  curl -sSf --retry 4 --retry-delay 5 -o $R/stage/fq/${cell}_2.gz "https://$u2" &&
  cat $R/stage/fq/${cell}_1.gz $R/stage/fq/${cell}_2.gz > $out &&
  rm -f $R/stage/fq/${cell}_1.gz $R/stage/fq/${cell}_2.gz || { echo "FAIL $cell"; rm -f $out; }
}
n=0
while IFS=$'\t' read -r cell srr u1 u2; do
  dl "$cell" "$u1" "$u2" & n=$((n+1)); [ $((n % 8)) -eq 0 ] && wait
done < $R/stage/manifest_full.clean.tsv
wait
echo "reads present: $(ls -1 $D/reads/*.fastq.gz | wc -l) of 118"; du -sh $D/reads
missing=$(comm -13 <(ls -1 $D/reads | sed 's/\.fastq\.gz$//' | LC_ALL=C sort) <(cut -f1 $R/stage/manifest_full.clean.tsv | LC_ALL=C sort))
if [ -n "$missing" ]; then echo "MISSING: $missing"; else touch $D/reads/.all118; echo "all 118 present"; fi
echo "=== download done $(date) ==="
