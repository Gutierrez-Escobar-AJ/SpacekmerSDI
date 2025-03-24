#!/bin/bash

# ================================
# Test Run for SpacekmerSDI_v3.py
# ================================

# Activate conda environment (uncomment if needed)
# conda activate spacekmersdi

# Define input and output paths
INPUT_DIR="./data/genomes"
MASK="11100111"
MAX_CLUSTERS=5
OUTPUT="results/test_spaced_kmer_SDI_results.csv"
MATRIX="results/test_spaced_kmer_SDI_matrix.csv"
DISTANCE="wasserstein"
CLUSTERING="agglomerative"
LINKAGE="ward"

# Create output directory if it doesn't exist
mkdir -p results

# Run the tool
python3 SpacekmerSDI_v3.py \
    --input "$INPUT_DIR" \
    --mask "$MASK" \
    --max_clusters "$MAX_CLUSTERS" \
    --output "$OUTPUT" \
    --output_matrix "$MATRIX" \
    --distance_metric "$DISTANCE" \
    --clustering "$CLUSTERING" \
    --linkage "$LINKAGE" \
    --verbose
