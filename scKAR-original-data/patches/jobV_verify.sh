#!/bin/bash
#$ -cwd
#$ -o logs/verify_reads.log
#$ -j y
#$ -l h_rt=4:00:00,h_data=4G
#$ -pe shared 8
set -u
R=/u/scratch/m/mhuang/claude/scKAR; D=$R/data
echo "=== verify start $(hostname) $(date) ==="
rm -f $D/reads/.all118
curl -sSf "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP063840&result=read_run&fields=run_accession,sample_title,fastq_ftp,fastq_bytes&format=tsv" -o $R/stage/ena_full.tsv
awk -F'\t' 'NR>1 && $2 ~ /_SC_/ {n=split($4,b,";"); t=0; for(i=1;i<=n;i++) t+=b[i]; print $2"\t"t}' $R/stage/ena_full.tsv | LC_ALL=C sort > $R/stage/expected_bytes.tsv
echo "expected entries: $(wc -l < $R/stage/expected_bytes.tsv)"
BAD=$R/stage/bad_cells.txt; : > $BAD
while IFS=$'\t' read -r cell exp; do
  f=$D/reads/$cell.fastq.gz
  if [ ! -f "$f" ]; then echo "$cell MISSING"; echo "$cell" >> $BAD; continue; fi
  act=$(stat -c %s "$f")
  if [ "$act" -ne "$exp" ]; then echo "$cell SIZE $act != $exp"; echo "$cell" >> $BAD; fi
done < $R/stage/expected_bytes.tsv
echo "size-mismatched or missing: $(wc -l < $BAD)"
# re-download the bad ones
if [ -s $BAD ]; then
  tr -d '\r' < $R/stage/manifest_full.tsv > $R/stage/manifest_full.clean.tsv
  n=0
  while read -r cell; do
    line=$(grep -P "^$cell\t" $R/stage/manifest_full.clean.tsv); u1=$(echo "$line"|cut -f3); u2=$(echo "$line"|cut -f4)
    ( rm -f $D/reads/$cell.fastq.gz
      curl -sSf --retry 6 --retry-delay 10 -o $R/stage/fq/${cell}_1.gz "https://$u1" &&
      curl -sSf --retry 6 --retry-delay 10 -o $R/stage/fq/${cell}_2.gz "https://$u2" &&
      cat $R/stage/fq/${cell}_1.gz $R/stage/fq/${cell}_2.gz > $D/reads/$cell.fastq.gz &&
      rm -f $R/stage/fq/${cell}_1.gz $R/stage/fq/${cell}_2.gz && echo "refetched $cell" || echo "REFETCH FAIL $cell" ) &
    n=$((n+1)); [ $((n % 6)) -eq 0 ] && wait
  done < $BAD
  wait
fi
# final pass: every cell present and byte-exact
FAIL=0
while IFS=$'\t' read -r cell exp; do
  f=$D/reads/$cell.fastq.gz
  if [ ! -f "$f" ] || [ "$(stat -c %s "$f")" -ne "$exp" ]; then echo "STILL BAD: $cell"; FAIL=1; fi
done < $R/stage/expected_bytes.tsv
echo "reads: $(ls -1 $D/reads/*.fastq.gz | wc -l)/118  $(du -sh $D/reads | cut -f1)"
if [ $FAIL -eq 0 ]; then echo verified > $D/reads/.all118; echo "ALL 118 BYTE-EXACT"; else echo "VERIFY INCOMPLETE"; fi
echo "=== verify done $(date) ==="
