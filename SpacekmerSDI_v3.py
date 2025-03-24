#!/usr/bin/env python3
import os
import sys
import logging
import argparse
from collections import Counter
from itertools import combinations
from typing import List, Tuple

import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.stats import wasserstein_distance
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
from time import time

# Optional: Try to import tqdm for progress bars; if not available, define a dummy function.
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **kwargs):
        return iterable

# Configure logging with more detailed format.
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# ===========================
# Robustness: Input validation
# ===========================
def validate_mask(mask: str) -> None:
    """Validate that the mask is non-empty and only contains '0' and '1' characters."""
    if not mask:
        raise ValueError("Mask must be non-empty.")
    if any(c not in '01' for c in mask):
        raise ValueError("Mask must contain only '0' and '1' characters.")

# ===========================
# Spaced k-mer Extraction and Frequency Calculation
# ===========================
def spaced_kmer(sequence: str, mask: str) -> List[str]:
    """
    Extract spaced k-mers from a sequence using a binary mask.
    Skips k-mers containing non-standard nucleotides (non-ATGC) in masked positions.
    """
    k = len(mask)
    valid_nucleotides = {'A', 'T', 'C', 'G'}
    valid_kmers = []
    
    for i in range(len(sequence) - k + 1):
        valid = True
        # Check all masked positions for valid nucleotides
        for j, m in enumerate(mask):
            if m == '1':
                try:
                    if sequence[i+j].upper() not in valid_nucleotides:
                        valid = False
                        break
                except IndexError:
                    valid = False
                    break
        if valid:
            kmer = ''.join([sequence[i + j] for j, m in enumerate(mask) if m == '1'])
            valid_kmers.append(kmer)
            
    return valid_kmers

def calculate_spaced_kmer_frequencies(sequence: str, mask: str) -> Tuple[int, float, int, int]:
    """
    Compute spaced k-mer counts, entropy, sequence length, and skipped k-mers.
    Mod: Added validation for non-standard nucleotides in masked positions.
    """
    k = len(mask)
    if len(sequence) < k:
        return 0, 0.0, len(sequence), 0
    
    spaced_kmers = spaced_kmer(sequence, mask)
    total = len(spaced_kmers)
    total_possible = len(sequence) - k + 1
    skipped_kmers = total_possible - total
    
    if total == 0:
        return 0, 0.0, len(sequence), skipped_kmers
    
    counts = Counter(spaced_kmers)
    freqs = np.array(list(counts.values())) / total
    entropy = -np.sum(freqs * np.log2(freqs + 1e-10))
    return total, entropy, len(sequence), skipped_kmers

# ===========================
# FASTA File Processing
# ===========================
def parse_fasta_file(file_path: str, mask: str) -> List[List]:
    """
    Process a FASTA file to extract spaced k-mer features and track non-standard nucleotides.
    Mod: Added non-standard nucleotide tracking and k-mer validation.
    """
    try:
        records = list(SeqIO.parse(file_path, "fasta"))
        if not records:
            logging.warning(f"No records found in {file_path}.")
            return []
            
        # Concatenate all sequences and count non-standard nucleotides
        full_seq = ''.join(str(r.seq) for r in records)
        seq_len = len(full_seq)
        non_standard = sum(1 for c in full_seq if c.upper() not in {'A', 'T', 'C', 'G'})
        
        # Calculate spaced k-mer features
        total, entropy, _, skipped = calculate_spaced_kmer_frequencies(full_seq, mask)
        
        return [[os.path.basename(file_path), total, entropy, seq_len, non_standard, skipped]]
    except Exception as e:
        logging.error(f"Error processing {file_path}: {str(e)}")
        return []

def process_genomes(folder_path: str, mask: str) -> pd.DataFrame:
    """
    Process all FASTA files in a directory with parallel execution.
    Mod: Added collection of non-standard nucleotide statistics.
    """
    if not os.path.isdir(folder_path):
        raise ValueError(f"Input folder '{folder_path}' does not exist or is not a directory.")

    files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith(".fasta")]
    if not files:
        raise ValueError(f"No FASTA files found in {folder_path}")

    logging.info(f"Found {len(files)} FASTA files. Processing...")
    from multiprocessing import Pool, cpu_count
    with Pool(min(cpu_count(), len(files))) as pool:
        results = list(tqdm(pool.starmap(parse_fasta_file, [(f, mask) for f in files]),
                            total=len(files), desc="Processing files"))
    
    # Flatten the list of results
    data = [item for sublist in results for item in sublist]
    if not data:
        raise ValueError("No valid genome data was processed.")
    return pd.DataFrame(data, columns=['Genome', 'SpacedKmer_Count', 'Entropy', 
                                      'Sequence_Length', 'NonStandard_Count', 'Skipped_Kmers'])

# ===========================
# Clustering Functions
# ===========================
def determine_optimal_clusters(scaled_data: np.ndarray, max_clusters: int,
                               clustering_method: str, clustering_linkage: str) -> Tuple[int, float]:
    """
    Determine the optimal number of clusters based on silhouette score.
    Mod: Now supports both 'agglomerative' and 'kmeans' methods.
    Returns a tuple of (optimal_clusters, best_silhouette_score)
    """
    best_score = -1
    optimal_clusters = 2
    
    for n in range(2, max_clusters + 1):
        if clustering_method == "agglomerative":
            model = AgglomerativeClustering(n_clusters=n, linkage=clustering_linkage)
            labels = model.fit_predict(scaled_data)
        elif clustering_method == "kmeans":
            model = KMeans(n_clusters=n, random_state=42, n_init=10)
            labels = model.fit_predict(scaled_data)
        else:
            raise ValueError(f"Unsupported clustering method: {clustering_method}")
        
        score = silhouette_score(scaled_data, labels)
        if score > best_score:
            best_score = score
            optimal_clusters = n
    
    logging.info(f"Optimal clusters determined: {optimal_clusters} (silhouette: {best_score:.3f})")
    return optimal_clusters, best_score

def perform_clustering(df: pd.DataFrame, max_clusters: int, clustering_method: str,
                       clustering_linkage: str) -> Tuple[pd.Series, int, float]:
    """
    Perform clustering on the feature data and determine the optimal number of clusters.
    Mod: Supports configurable clustering methods.
    """
    features = ['SpacedKmer_Count', 'Entropy']
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(df[features])
    
    # Determine optimal clusters and record silhouette score once.
    optimal_clusters, best_silhouette = determine_optimal_clusters(scaled_data, max_clusters,
                                                                   clustering_method, clustering_linkage)
    
    # Perform final clustering using the optimal number of clusters.
    if clustering_method == "agglomerative":
        model = AgglomerativeClustering(n_clusters=optimal_clusters, linkage=clustering_linkage)
        clusters = model.fit_predict(scaled_data)
    elif clustering_method == "kmeans":
        model = KMeans(n_clusters=optimal_clusters, random_state=42, n_init=10)
        clusters = model.fit_predict(scaled_data)
    else:
        raise ValueError(f"Unsupported clustering method: {clustering_method}")
    
    logging.info(f"Clustering completed using method '{clustering_method}'.")
    return pd.Series(clusters, index=df.index, name='Cluster'), optimal_clusters, best_silhouette

# ===========================
# SDI Matrix Computation with Configurable Distance Metric
# ===========================
def compute_sdi_matrix(df: pd.DataFrame, distance_metric: str) -> pd.DataFrame:
    """
    Compute the SDI (Spaced k-mer Distance Index) matrix between clusters.
    Mod: Added support for 'euclidean' distance metric.
         Optimized the nested loop using combinations to avoid redundant computations.
    """
    cluster_ids = sorted(df['Cluster'].unique())
    sdi_matrix = pd.DataFrame(index=cluster_ids, columns=cluster_ids, dtype=float)

    # Pre-compute cluster feature subsets.
    cluster_features = {
        cid: df[df['Cluster'] == cid][['SpacedKmer_Count', 'Entropy']]
        for cid in cluster_ids
    }
    
    # Compute pairwise distances only once per unique pair.
    for i, j in combinations(cluster_ids, 2):
        cluster_i = cluster_features[i]
        cluster_j = cluster_features[j]
        
        if distance_metric == "wasserstein":
            wd_kmer = wasserstein_distance(cluster_i['SpacedKmer_Count'], cluster_j['SpacedKmer_Count'])
            wd_entropy = wasserstein_distance(cluster_i['Entropy'], cluster_j['Entropy'])
            d = wd_kmer + wd_entropy
        elif distance_metric == "euclidean":
            # Compute Euclidean distance between the means of the clusters.
            mean_i = cluster_i.mean()
            mean_j = cluster_j.mean()
            d = np.linalg.norm(mean_i - mean_j)
        else:
            raise ValueError(f"Unsupported distance metric: {distance_metric}")
        
        # Log-transform and ensure a minimum value.
        sdi = np.log1p(d)
        sdi = max(sdi, 0.01)
        
        sdi_matrix.loc[i, j] = sdi
        sdi_matrix.loc[j, i] = sdi  # Symmetric matrix

    # Fill diagonal with zeros.
    for cid in cluster_ids:
        sdi_matrix.loc[cid, cid] = 0.0

    logging.info(f"SDI matrix computed using '{distance_metric}' distance.")
    return sdi_matrix

# ===========================
# Statistics File Generation
# ===========================
def generate_stats_file(df: pd.DataFrame, output_path: str) -> None:
    """Generate statistics file about non-standard nucleotides and skipped k-mers."""
    try:
        with open(output_path, 'w') as f:
            # Write header
            f.write("Genome\tNonStandard_Count\tNonStandard_Percent(%)\tSkipped_Kmers\n")
            
            # Write data rows
            for _, row in df.iterrows():
                genome = row['Genome']
                non_std = row['NonStandard_Count']
                seq_len = row['Sequence_Length']
                skipped = row['Skipped_Kmers']
                
                percent = (non_std / seq_len * 100) if seq_len > 0 else 0.0
                f.write(f"{genome}\t{non_std}\t{percent:.2f}\t{skipped}\n")
                
        logging.info(f"Non-standard nucleotide statistics saved to {output_path}")
    except Exception as e:
        logging.error(f"Failed to generate statistics file: {str(e)}")

# ===========================
# Main Execution Flow
# ===========================
def main():
    parser = argparse.ArgumentParser(description='Enhanced Spaced k-mer SDI Analysis')
    parser.add_argument('--input', required=True, help='Input directory with FASTA files')
    parser.add_argument('--mask', required=True, help='Spaced k-mer mask (e.g., 11100111)')
    parser.add_argument('--max_clusters', type=int, default=10, 
                        help='Maximum clusters for automatic determination')
    parser.add_argument('--output', default='improved_spaced_kmerSDI_results.csv', 
                        help='Results output file')
    parser.add_argument('--output_matrix', default='improved_spaced_kmerSDI_matrix.csv', 
                        help='SDI matrix output file')
    parser.add_argument('--distance_metric', choices=['wasserstein', 'euclidean'], default='wasserstein',
                        help="Distance metric to compute SDI matrix ('wasserstein' or 'euclidean')")
    parser.add_argument('--clustering', choices=['agglomerative', 'kmeans'], default='agglomerative',
                        help="Clustering method to use ('agglomerative' or 'kmeans')")
    parser.add_argument('--linkage', choices=['ward', 'complete', 'average', 'single'], default='ward',
                        help="Linkage criterion for Agglomerative Clustering (ignored for kmeans)")
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        validate_mask(args.mask)
        start_time = time()
        
        logging.info("Starting genome data processing...")
        df = process_genomes(args.input, args.mask)
        logging.info(f"Processed {len(df)} genome files.")

        logging.info("Generating statistics file...")
        stats_path = os.path.splitext(args.output)[0] + '_stats.txt'
        generate_stats_file(df, stats_path)

        logging.info("Performing optimized clustering...")
        clusters, n_clusters, best_silhouette = perform_clustering(df, args.max_clusters, args.clustering, args.linkage)
        df = df.join(clusters)

        logging.info("Computing SDI matrix...")
        sdi_matrix = compute_sdi_matrix(df, args.distance_metric)
        
        # Add overall silhouette score
        features = ['SpacedKmer_Count', 'Entropy']
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(df[features])
        overall_silhouette = silhouette_score(scaled_data, df['Cluster'])
        df['Silhouette_Score'] = overall_silhouette
        
        # Save outputs
        df.to_csv(args.output, index=False)
        sdi_matrix.to_csv(args.output_matrix)
        
        elapsed = time() - start_time
        logging.info(f"Analysis complete. Total processing time: {elapsed:.2f} seconds")
        
    except Exception as e:
        logging.error(f"Critical error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()  
 
# python3 SpacekmerSDI_v3.py --input ./mycoplasma --mask 11100111 --max_clusters 10 --output results.csv --output_matrix sdi_matrix.csv --distance_metric wasserstein --clustering agglomerative --linkage ward --verbose   
