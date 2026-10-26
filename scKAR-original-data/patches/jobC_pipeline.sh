#!/bin/bash
#$ -cwd
#$ -o logs/full_pipeline.log
#$ -j y
#$ -l h_rt=23:59:00,h_data=8G
#$ -pe shared 8
set -u
R=/u/scratch/m/mhuang/claude/scKAR; D=$R/data
cd $R; source $R/env.sh
export TMPDIR=$R/.tmp; mkdir -p "$TMPDIR"
A=$(cat $R/.attempt 2>/dev/null || echo 0); A=$((A+1)); echo $A > $R/.attempt
echo "=== pipeline attempt $A on $(hostname) $(date) ==="
if [ ! -f $D/reads/.all118 ]; then
  echo "reads incomplete ($(ls -1 $D/reads/*.fastq.gz 2>/dev/null | wc -l)/118) -> requeueing pipeline"
  [ "$A" -lt 10 ] && qsub jobC_pipeline.sh
  exit 0
fi
[ -s $D/ref_31mers.txt ] || { echo "ABORT: ref_31mers.txt missing"; exit 1; }

# one-time setup on first attempt
if [ ! -f $R/.full_setup_done ]; then
  echo "--- setup: 118-cell matrix, config, ref_filter bounds fix ---"
  cp $R/stage/expression_full.csv $D/expression_matrix/expression.csv
  awk -F, 'NR==1{print "expression.csv cols: "NF-1} NR>1{n++} END{print "genes: "n}' $D/expression_matrix/expression.csv
  sed -i -e 's|^FILTER_REFERENCE_KMERS=.*|FILTER_REFERENCE_KMERS=true|' \
         -e 's|^MODE=.*|MODE="GENE_EXPRESSION"|' -e 's|^N_PCS=.*|N_PCS=40|' \
         -e 's|^N_NEIGHBORS=.*|N_NEIGHBORS=10|' -e 's|^RESOLUTION=.*|RESOLUTION=0.6|' \
         -e "s|^REFERENCE_TRANSCRIPTOME_SORTED_31MERS_PATH=.*|REFERENCE_TRANSCRIPTOME_SORTED_31MERS_PATH=\"$D/ref_31mers.txt\"|" $R/repo/src/config.env
  F=$R/repo/src/kmer-filtering/ref_filter.cpp
  sed -i 's|^        string ref_kmer = ref_data\[ref_index\];|        string ref_kmer = (ref_index < ref_data_size) ? ref_data[ref_index] : string();|' $F
  sed -i 's|^            ref_kmer = ref_data\[ref_index\];|            ref_kmer = (ref_index < ref_data_size) ? ref_data[ref_index] : string();|' $F
  grep -n 'ref_index < ref_data_size' $F
  rm -rf $D/adj $D/f_test_results $D/deseq_results $D/final_results $D/clustering_results $D/cluster_tpm.csv
  touch $R/.full_setup_done
fi
cat $R/repo/src/config.env
. $R/repo/src/config.env

# stage 1: k-mers + adjacency + reference filtering
if [ -f $D/adj/.done ]; then echo "[1] skip"; else
  echo "[1] mother_script $(date)"; rm -f $D/adj/adj_*.csv
  cd $R/repo/src/kmer-filtering
  bash ./mother_script.sh "$INPUT_DIR" "$FILTER_REFERENCE_KMERS" "$REFERENCE_TRANSCRIPTOME_SORTED_31MERS_PATH" "$NUMBER_OF_THREADS"
  s=$?; echo "[1] exit $s $(date)"; [ $s -eq 0 ] && touch $D/adj/.done
fi
ls -1 $D/adj/*.csv 2>/dev/null | wc -l; du -sh $D/adj 2>/dev/null; wc -l < $D/tpm_sum.csv

# stage 2: conditions
if [ -s $D/clustering_results/bipartitions/bipartition_0.csv ]; then echo "[2] skip"; else
  echo "[2] clustering $(date)"; cd $R/repo/src/clustering
  for RES in $RESOLUTION 1.0 1.5 2.0; do
    python3 graph-clustering.py "$INPUT_DIR"/expression_matrix/expression.csv $CLUSTERING_ALGO $MIN_GENES $MIN_CELLS $N_NEIGHBORS $N_PCS $RES 2>&1 | tail -2
    [ -s $D/clustering_results/bipartitions/bipartition_0.csv ] && { echo "[2] ok at resolution=$RES"; break; }
  done
  echo "[2] exit $? $(date)"
fi
echo "clusters:"; cut -d, -f2 $D/clustering_results/cluster.csv 2>/dev/null | tail -n +2 | sort | uniq -c
echo "bipartitions: $(ls -1 $D/clustering_results/bipartitions 2>/dev/null | wc -l)"

# stage 3: F-test
if [ -f $D/f_test_results/.done ]; then echo "[3] skip"; else
  echo "[3] F-test $(date)"; cd $R/repo/src/f-test
  python3 preprocess_clustering_results.py "$INPUT_DIR"/tpm_sum.csv "$INPUT_DIR"/clustering_results/cluster.csv "$INPUT_DIR"/cluster_tpm.csv
  wc -l < "$INPUT_DIR"/cluster_tpm.csv
  rm -f $D/f_test_results/f_test_result_*.csv
  ./f_test.out "$INPUT_DIR" $PSEUDOBULK_SIZE $MIN_ROW_THREDSHOLD
  s=$?; echo "[3] exit $s $(date)"; [ $s -eq 0 ] && touch $D/f_test_results/.done
fi
echo "f_test chunks: $(ls -1 $D/f_test_results/*.csv 2>/dev/null | wc -l) rows: $(cat $D/f_test_results/*.csv 2>/dev/null | wc -l)"

# stage 4: DESeq2 (resumable per output file)
echo "[4] DESeq2 $(date)"; cd $R/repo/src/de-test
bash ./deseq_runner.sh "$INPUT_DIR" $MIN_ROW_COUNT $MIN_COL_COUNT $NUMBER_OF_THREADS $LOG2FC $PVAL $BASE_MEAN_THRESHOLD
echo "[4] exit $? $(date); tables: $(ls -1 $D/deseq_results 2>/dev/null | wc -l)"

# stage 5: assembly
echo "[5] assembly $(date)"; cd $R/repo/src/merge_kmers
for f in $(find $D/final_results -mindepth 1 -maxdepth 1 -type d 2>/dev/null); do
  bash ./create_contigs.sh $f $BLAT_PATH "$INPUT_DIR"/genome_assembly/genome.2bit
  bash ./combine_p_values.sh $f $PVAL_COMBINATION_METHOD
done
echo "[5] exit $? $(date)"
find $D/final_results -name '*_pvalues.tsv' 2>/dev/null | wc -l

# resubmit if not finished
NB=$(ls -1 $D/clustering_results/bipartitions 2>/dev/null | wc -l)
NP=$(find $D/final_results -name '*_contigs_pvalues.tsv' 2>/dev/null | wc -l)
echo "progress: bipartitions=$NB contig_pvalue_files=$NP attempt=$A"
if [ "$NP" -lt $((NB*2)) ] && [ "$A" -lt 10 ]; then
  echo "incomplete -> resubmitting"; cd $R; qsub jobC_pipeline.sh
else
  echo "=== PIPELINE COMPLETE ==="
fi
echo "=== attempt $A done $(date) ==="
