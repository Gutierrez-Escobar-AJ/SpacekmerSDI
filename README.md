# 🧬 SpacekmerSDI

**Spaced k-mer Based Structural Differentiation Index Analysis**  
**Version:** 3.0  
**Author:** Gutierrez-Escobar-AJ  
**Repository:** https://github.com/yourusername/SpacekmerSDI  

---

## 🔬 Overview

**SpacekmerSDI** is a lightweight yet powerful command-line tool for genome-wide structural differentiation analysis using **spaced k-mer frequencies** and **unsupervised clustering**. It enables rapid comparative analysis of genomes by extracting spaced k-mer signatures using a binary mask, computing entropy-based features, and clustering genomes based on their feature space. The tool then calculates **Structural Differentiation Index (SDI)** matrices using **Wasserstein** or **Euclidean** distances to quantify pairwise dissimilarities.

This approach is highly suitable for:

- **Comparative genomics**
- **Microbial population structure analysis**
- **Phylogeographic clustering**
- **Methylation or epigenomic profile stratification (with pre-masked input)**

---

## ✨ Key Features

- ✅ **Spaced k-mer extraction** with user-defined binary masks (`e.g., 11100111`)
- 📈 **Entropy-based** feature calculation for measuring sequence complexity
- 🤖 **Unsupervised clustering** (Agglomerative or K-Means)
- 📊 **Silhouette Score Optimization** for automatic cluster number selection
- 🔁 **Pairwise SDI matrix computation** with Wasserstein or Euclidean distances
- 📦 **Multi-core parallel FASTA processing**
- 📁 Generates log files and non-standard nucleotide summaries
- 🧪 Supports genome, contig, or methylation-masked FASTA input

---

## ⚙️ Installation

```bash
git clone https://github.com/yourusername/SpacekmerSDI.git
cd SpacekmerSDI
conda env create -f environment.yml
conda activate spacekmersdi
```

> ✅ Requires Conda and Python 3.7+.

---

## ▶️ Usage

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

## 📂 Input Requirements

- **Input Directory:** Folder containing `.fasta` files (e.g., assembled genomes, contigs, or methylation-masked sequences).
- **Mask Format:** Binary string where `1` denotes positions to include in the k-mer (e.g., `11100111`).

---

## 📤 Output Files

| File | Description |
|------|-------------|
| `spaced_kmer_SDI_results.csv` | Feature table with entropy, k-mer counts, and assigned cluster |
| `spaced_kmer_SDI_matrix.csv` | SDI matrix representing inter-cluster dissimilarities |
| `*_stats.txt` | Summary of non-standard nucleotides and skipped k-mers per genome |
| Logs | Detailed processing logs (if `--verbose` is used) |

---

## 📊 Example Use Cases

- Explore structural variation across *H. pylori* strains.
- Differentiate methylated vs. unmethylated sequences using masked input.
- Detect population clusters in metagenomic or pan-genomic datasets.
- Identify outlier genomes or contaminants via entropy/skipped-kmer profiles.

---

## 🧠 Under the Hood

- **Spaced k-mer Frequency & Entropy Calculation:**  
  Extracts k-mers from positions marked `1` in a binary mask and computes frequency distribution entropy.

- **Unsupervised Clustering with Optimization:**  
  Uses **silhouette score** to select the optimal number of clusters before applying KMeans or Agglomerative clustering.

- **SDI Matrix Calculation:**  
  Computes Wasserstein or Euclidean distances between feature distributions of clusters, then applies `log1p()` transformation.

---

## 🧾 License

This project is licensed under the [MIT License](LICENSE).

---

## 📣 Citation

If you use **SpacekmerSDI** in your research, please cite:

```
Title: "SpacekmerSDI: Spaced k-mer Based Structural Differentiation Index Analysis"
version: 3.0.0
doi: 10.5281/zenodo.15079737
date-released: 2024-03-24
url: https://github.com/Gutierrez-Escobar-AJ/SpacekmerSDI
repository-code: https://github.com/Gutierrez-Escobar-AJ/SpacekmerSDI
license: MIT
}
```

---

## 💡 Logo

![SpacekmerSDI Logo](./logo.png)

---

## 🤝 Contributions

Feel free to fork this repo and submit pull requests. Feature suggestions and bug reports are welcome via GitHub Issues.

---

## 🙌 Acknowledgments

Developed with inspiration from spaced k-mer approaches and their applications in comparative genomics and epigenomics.

---
