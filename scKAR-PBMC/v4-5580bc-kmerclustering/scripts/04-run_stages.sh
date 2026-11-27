#!/bin/bash
: "${NSLOTS:=${SLURM_CPUS_PER_TASK:-8}}"
export NSLOTS
#$ -cwd
#$ -o logs/stages.log
#$ -j y
#$ -l h_rt=4:00:00,h_data=8G
#$ -pe shared 8
set -u
ROOT=/u/scratch/m/mhuang/claude/compare3/sckar_kmerclust
D=$ROOT/data
cd $ROOT
echo "=== host $(hostname) slots $NSLOTS start $(date) ==="
source $ROOT/env.sh
export TMPDIR=$ROOT/.tmp; mkdir -p "$TMPDIR"
. $ROOT/repo/src/config.env

# ---- stage 1: k-mer counting + filtering + adjacency ----
if [ -f $D/adj/.done ]; then echo "[1] adjacency already done, skipping"; else
  echo "[1] mother_script (jellyfish/filter cached, then adjacency) $(date)"
  rm -f $D/adj/adj_*.csv
  cd $ROOT/repo/src/kmer-filtering
  bash ./mother_script.sh "$INPUT_DIR" "$FILTER_REFERENCE_KMERS" "$REFERENCE_TRANSCRIPTOME_SORTED_31MERS_PATH" "$NUMBER_OF_THREADS"
  s=$?; echo "[1] exit $s $(date)"; [ $s -eq 0 ] && touch $D/adj/.done
fi
ls -l $D/adj | tail -3

# ---- stage 2: clustering / bipartitions ----
if [ -s $D/clustering_results/bipartitions/bipartition_0.csv ]; then echo "[2] conditions already done, skipping"; else
  echo "[2] condition generation (MODE=$MODE) $(date)"
  cd $ROOT/repo/src/clustering
  if [ "$MODE" = "CUSTOM_METADATA" ]; then
    python3 generate_bipartitions_from_metadata.py "$INPUT_DIR"
  elif [ "$MODE" = "KMER_ABUNDANCE" ]; then
    if [ ! -s "$INPUT_DIR"/expression_matrix/kmer_abundance.csv ]; then
      python3 $ROOT/scripts/01-build_kmer_matrix.py $ROOT $MIN_ROW_THREDSHOLD "$INPUT_DIR"/expression_matrix/kmer_abundance.csv
    fi
    ls -lh "$INPUT_DIR"/expression_matrix/kmer_abundance.csv
    for RES in $RESOLUTION 1.0 1.5 2.0 3.0; do
      echo "[2] trying resolution=$RES"
      python3 graph-clustering.py "$INPUT_DIR"/expression_matrix/kmer_abundance.csv $CLUSTERING_ALGO $MIN_GENES $MIN_CELLS $N_NEIGHBORS $N_PCS $RES 2>&1 | tail -3
      if [ -s $D/clustering_results/bipartitions/bipartition_0.csv ]; then echo "[2] succeeded at resolution=$RES"; break; fi
    done
  else
    for RES in $RESOLUTION 1.0 1.5 2.0 3.0; do
      echo "[2] trying resolution=$RES"
      python3 graph-clustering.py "$INPUT_DIR"/expression_matrix/expression.csv $CLUSTERING_ALGO $MIN_GENES $MIN_CELLS $N_NEIGHBORS $N_PCS $RES 2>&1 | tail -3
      if [ -s $D/clustering_results/bipartitions/bipartition_0.csv ]; then echo "[2] succeeded at resolution=$RES"; break; fi
    done
  fi
  echo "[2] exit $? $(date)"
fi
echo "clusters:"; cut -d, -f2 $D/clustering_results/cluster.csv 2>/dev/null | sort | uniq -c
echo "bipartitions:"; ls -1 $D/clustering_results/bipartitions 2>/dev/null | wc -l

# ---- stage 3: F-test ----
if [ -f $D/f_test_results/.done ]; then echo "[3] F-test already done, skipping"; else
  echo "[3] F-test $(date)"
  cd $ROOT/repo/src/f-test
  python3 preprocess_clustering_results.py "$INPUT_DIR"/tpm_sum.csv "$INPUT_DIR"/clustering_results/cluster.csv "$INPUT_DIR"/cluster_tpm.csv
  head -3 "$INPUT_DIR"/cluster_tpm.csv; wc -l "$INPUT_DIR"/cluster_tpm.csv
  rm -f $D/f_test_results/f_test_result_*.csv
  ./f_test.out "$INPUT_DIR" $PSEUDOBULK_SIZE $MIN_ROW_THREDSHOLD
  s=$?; echo "[3] exit $s $(date)"; [ $s -eq 0 ] && touch $D/f_test_results/.done
fi
echo "f_test chunks:"; ls -1 $D/f_test_results/*.csv 2>/dev/null | wc -l; du -sh $D/f_test_results 2>/dev/null

# ---- stage 4: DESeq2 (resumable per file) ----
echo "[4] DESeq2 $(date)"
cd $ROOT/repo/src/de-test
bash ./deseq_runner.sh "$INPUT_DIR" $MIN_ROW_COUNT $MIN_COL_COUNT $NUMBER_OF_THREADS $LOG2FC $PVAL $BASE_MEAN_THRESHOLD
echo "[4] exit $? $(date)"
ls -1 $D/deseq_results 2>/dev/null | wc -l

# ---- stage 5: assembly + contig p-values ----
echo "[5] assembly $(date)"
cd $ROOT/repo/src/merge_kmers
for folder in $(find $D/final_results -mindepth 1 -maxdepth 1 -type d); do
  bash ./create_contigs.sh $folder $BLAT_PATH "$INPUT_DIR"/genome_assembly/genome.2bit
  bash ./combine_p_values.sh $folder $PVAL_COMBINATION_METHOD
done
echo "[5] exit $? $(date)"
find $D/final_results -type f | head -30
echo "=== done $(date) ==="
