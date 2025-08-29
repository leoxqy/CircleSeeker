#!/usr/bin/env bash
set -euo pipefail

echo "[cleanup] Removing caches and temporary work..."
find . -type d -name "__pycache__" -prune -exec rm -rf {} + || true
rm -rf .pytest_cache || true

echo "[cleanup] Removing pipeline temp workdirs (.tmp_work)..."
find . -type d -name ".tmp_work" -prune -exec rm -rf {} + || true

echo "[cleanup] Removing checkpoints and logs..."
find . -type f -name "*.checkpoint" -delete || true
find . -type f -name "*.checkpoint.bak" -delete || true
find . -type f -name "*.log" -delete || true

echo "[cleanup] Removing large test artifacts under tests/ (fa,fasta,mmi,bam,bai,bed,tsv,csv,clstr)..."
find tests -type f \( \
  -name "*.mmi" -o -name "*.fa" -o -name "*.fasta" -o -name "*.fai" -o \
  -name "*.fq" -o -name "*.fastq" -o -name "*.bam" -o -name "*.bai" -o \
  -name "*.sam" -o -name "*.bed" -o -name "*.bedpe" -o -name "*.tsv" -o \
  -name "*.csv" -o -name "*.clstr" \) -delete 2>/dev/null || true

echo "[cleanup] Done. Stage with: git add -A"

