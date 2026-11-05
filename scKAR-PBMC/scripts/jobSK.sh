#!/bin/bash
#$ -cwd
#$ -o /u/scratch/m/mhuang/claude/compare3/logs/sckar.$JOB_ID.out
#$ -j y
#$ -l h_rt=23:59:00,h_data=12G
#$ -pe shared 4
set -u
R=/u/scratch/m/mhuang/claude/compare3/sckar_pbmc; D=$R/data
cd $R; source $R/env.sh
export TMPDIR=$R/.tmp; mkdir -p "$TMPDIR"
A=$(cat $R/.attempt 2>/dev/null || echo 0); A=$((A+1)); echo $A > $R/.attempt
echo "=== scKAR PBMC attempt $A on $(hostname) slots $NSLOTS $(date) ==="

NR=$(ls -1 $D/reads/*.fastq.gz 2>/dev/null | wc -l)
if [ ! -f $R/.demuxdone ] || [ "$NR" -lt 52 ]; then
  echo "reads incomplete ($NR/52) -> requeue"
  if [ "$A" -lt 10 ]; then cd $R; qsub jobSK.sh; fi
  exit 0
fi
[ -s $D/ref_31mers.txt ] || { echo "ABORT ref_31mers.txt missing"; exit 1; }
[ -s $D/expression_matrix/expression.csv ] || { echo "ABORT expression.csv missing"; exit 1; }
. $R/repo/src/config.env
cat $R/repo/src/config.env
echo "reads: $NR  $(du -sh $D/reads | cut -f1)"

# stage 1: jellyfish count + filter + adjacency + reference subtraction
if [ -f $D/adj/.done ]; then echo "[1] skip"; else
  echo "[1] mother_script $(date)"
  cd $R/repo/src/kmer-filtering
  bash ./mother_script.sh "$INPUT_DIR" "$FILTER_REFERENCE_KMERS" "$REFERENCE_TRANSCRIPTOME_SORTED_31MERS_PATH" "$NUMBER_OF_THREADS"
  s=$?; echo "[1] exit $s $(date)"; if [ $s -eq 0 ]; then touch $D/adj/.done; fi
fi
echo "adj chunks: $(ls -1 $D/adj/*.csv 2>/dev/null | wc -l)  $(du -sh $D/adj 2>/dev/null | cut -f1)"
echo "jellyfish: $(du -sh $D/jellyfish 2>/dev/null | cut -f1)"
wc -l < $D/tpm_sum.csv

# stage 2: leiden clustering of pseudo-cells -> bipartitions
if [ -s $D/clustering_results/bipartitions/bipartition_0.csv ]; then echo "[2] skip"; else
  echo "[2] clustering $(date)"; cd $R/repo/src/clustering
  for RES in $RESOLUTION 1.0 1.5 2.0; do
    python3 graph-clustering.py "$INPUT_DIR"/expression_matrix/expression.csv $CLUSTERING_ALGO $MIN_GENES $MIN_CELLS $N_NEIGHBORS $N_PCS $RES 2>&1 | tail -3
    if [ -s $D/clustering_results/bipartitions/bipartition_0.csv ]; then echo "[2] ok at resolution=$RES"; break; fi
  done
fi
echo "clusters:"; tail -n +2 $D/clustering_results/cluster.csv 2>/dev/null | cut -d, -f2 | sort | uniq -c
echo "bipartitions: $(ls -1 $D/clustering_results/bipartitions 2>/dev/null | wc -l)"

# stage 3: F-test
if [ -f $D/f_test_results/.done ]; then echo "[3] skip"; else
  echo "[3] F-test $(date)"; cd $R/repo/src/f-test
  python3 preprocess_clustering_results.py "$INPUT_DIR"/tpm_sum.csv "$INPUT_DIR"/clustering_results/cluster.csv "$INPUT_DIR"/cluster_tpm.csv
  rm -f $D/f_test_results/f_test_result_*.csv
  ./f_test.out "$INPUT_DIR" $PSEUDOBULK_SIZE $MIN_ROW_THREDSHOLD
  s=$?; echo "[3] exit $s $(date)"; if [ $s -eq 0 ]; then touch $D/f_test_results/.done; fi
fi
echo "f_test chunks: $(ls -1 $D/f_test_results/*.csv 2>/dev/null | wc -l)  $(du -sh $D/f_test_results 2>/dev/null | cut -f1)"

# stage 4: DESeq2
echo "[4] DESeq2 $(date)"; cd $R/repo/src/de-test
bash ./deseq_runner.sh "$INPUT_DIR" $MIN_ROW_COUNT $MIN_COL_COUNT $NUMBER_OF_THREADS $LOG2FC $PVAL $BASE_MEAN_THRESHOLD
echo "[4] exit $? $(date); tables: $(ls -1 $D/deseq_results 2>/dev/null | wc -l)"

# stage 5: ABySS assembly + BLAT + contig p-values
echo "[5] assembly $(date)"; cd $R/repo/src/merge_kmers
for f in $(find $D/final_results -mindepth 1 -maxdepth 1 -type d 2>/dev/null); do
  echo "--- $f"
  bash ./create_contigs.sh $f $BLAT_PATH "$INPUT_DIR"/genome_assembly/genome.2bit
  bash ./combine_p_values.sh $f $PVAL_COMBINATION_METHOD
done
echo "[5] exit $? $(date)"
find $D/final_results -type f 2>/dev/null | head -40

NB=$(ls -1 $D/clustering_results/bipartitions 2>/dev/null | wc -l)
NP=$(find $D/final_results -name '*contigs_pvalues.tsv' 2>/dev/null | wc -l)
echo "progress: bipartitions=$NB contig_pvalue_files=$NP attempt=$A"
if [ "$NB" -gt 0 ] && [ "$NP" -ge $((NB*2)) ]; then
  echo "=== PIPELINE COMPLETE ==="; echo DONE > $R/.sckardone
elif [ "$A" -lt 10 ]; then
  echo "incomplete -> resubmitting"; cd $R; qsub jobSK.sh
fi
echo "=== attempt $A done $(date) ==="
