# 📘 SpacekmerSDI User Manual

**Tool Name:** SpacekmerSDI  
**Version:** 3.0  
**Main Script:** SpacekmerSDI_v3.py  
**Author:** [Your Name]  
**Repository:** https://github.com/yourusername/SpacekmerSDI

---

## 1. 🔎 Introduction

**SpacekmerSDI** is a command-line Python tool designed for analyzing genomic structural variation using **spaced k-mers**, **entropy**, and **clustering**. It is particularly useful for researchers working in microbial genomics, epigenomics, and comparative genomics. The tool takes in multiple FASTA files, extracts spaced k-mers based on a binary mask, calculates entropy, clusters genomes, and computes the Structural Differentiation Index (SDI) between the clusters.

---

## 2. 📥 Input Requirements

- **FASTA Directory:** A folder containing `.fasta` files (assembled genomes, contigs, or methylation-masked sequences).
- **Spaced k-mer Mask:** A binary string (`e.g., 11100111`) where each `1` indicates positions included in the k-mer.

---

## 3. ⚙️ How It Works

### 3.1 Spaced k-mer Extraction
The tool extracts k-mers only from positions specified as `1` in the binary mask and skips k-mers containing non-standard nucleotides (non-ATGC).

### 3.2 Feature Calculation
For each genome:
- **SpacedKmer_Count:** Number of valid spaced k-mers.
- **Entropy:** Shannon entropy based on k-mer frequency.
- **Sequence_Length:** Total number of nucleotides.
- **NonStandard_Count:** Number of non-ATGC nucleotides.
- **Skipped_Kmers:** Number of invalid or skipped k-mers.

### 3.3 Clustering
The genomes are clustered using:
- **Agglomerative Clustering** or **KMeans**
- **Silhouette Score** is used to automatically select the optimal number of clusters.

### 3.4 SDI Matrix Computation
- **Wasserstein Distance** or **Euclidean Distance** is used to measure dissimilarities between clusters.
- Result is a log-transformed, symmetric SDI matrix.

---

## 4. 🚀 How To Run

### Basic Command:

```bash
python3 SpacekmerSDI_v3.py \
  --input ./genomes \
  --mask 11100111 \
  --max_clusters 10 \
  --output spaced_kmer_SDI_results.csv \
  --output_matrix spaced_kmer_SDI_matrix.csv \
  --distance_metric wasserstein \
  --clustering agglomerative \
  --linkage ward \
  --verbose
```

---

## 5. 🧾 Output Files

| File Name                        | Description                                           |
|----------------------------------|-------------------------------------------------------|
| `spaced_kmer_SDI_results.csv`    | Table with features + cluster ID + silhouette score   |
| `spaced_kmer_SDI_matrix.csv`     | Pairwise SDI matrix between genome clusters           |
| `*_stats.txt`                    | Non-standard base and skipped k-mer stats per genome |
| Log messages                     | Info/debug messages if `--verbose` is enabled         |

---

## 6. 🧠 Internals / Methods

### 6.1 Functions

- `spaced_kmer()`: Extracts valid spaced k-mers using mask.
- `calculate_spaced_kmer_frequencies()`: Computes k-mer entropy.
- `parse_fasta_file()`: Processes each FASTA, validates nucleotides.
- `process_genomes()`: Multiprocessing parser for all genomes.
- `perform_clustering()`: Uses `AgglomerativeClustering` or `KMeans`.
- `compute_sdi_matrix()`: Builds SDI matrix between clusters.

### 6.2 Distance Metrics

- **Wasserstein Distance:** Robust against distributional variation.
- **Euclidean Distance:** Fast and interpretable mean distance.

---

## 7. 🧪 Performance Tips

- Use pre-filtered FASTA sequences (remove Ns).
- Use binary masks of length 8–12 for spaced k-mers.
- Set `--max_clusters` to slightly more than expected biological groups.
- Enable `--verbose` for debugging and tracing performance.

---

## 8. 🔧 Troubleshooting

| Issue                               | Solution                                                |
|------------------------------------|----------------------------------------------------------|
| `No FASTA files found`             | Make sure input directory contains `.fasta` files        |
| `Mask must contain only '0' or '1'`| Check mask string format (`e.g., 11100111`)              |
| `No valid genome data was processed`| Ensure FASTA sequences have only valid ATGC bases        |
| Clustering too coarse/fine         | Adjust `--max_clusters` and/or try `--clustering kmeans` |

---

## 9. 📄 License

MIT License – see `LICENSE` file for details.

---

## 10. 📬 Contact

For bugs, issues, or feature requests, please open an [Issue](https://github.com/yourusername/SpacekmerSDI/issues) on GitHub.

---

