import matplotlib.pyplot as plt
from collections import Counter
import os

# --- Settings ---
FASTA_FILE_TO_ANALYZE = "MN908947.3.fasta"
KMER_LENGTH = 6
TOP_N_TO_PLOT = 20
# -----------------

def parse_fasta(filename):
    if not os.path.exists(filename):
        print(f"Error: File not found at {filename}")
        return None, 0

    sequence = ""
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            sequence += line.strip()
    
    print(f"Parsed {filename}. Genome length: {len(sequence)} bp")
    return sequence, len(sequence)

def find_frequent_kmers(sequence, k, top_n):
    kmer_counts = Counter()
    
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i : i+k]
        kmer_counts[kmer] += 1
        
    most_common_list = kmer_counts.most_common(top_n)
    
    return dict(most_common_list)

def plot_kmer_frequencies(kmer_counts, filename, k, top_n, genome_len):
    labels = list(kmer_counts.keys())
    counts = list(kmer_counts.values())
    
    labels.reverse()
    counts.reverse()
    
    plt.figure(figsize=(10, 8))
    plt.barh(labels, counts, color='steelblue')
    
    title = f"Top {top_n} Most Frequent {k}-mers"
    subtitle = f"in {filename} (Genome Length: {genome_len} bp)"
    plt.title(f"{title}\n{subtitle}")
    plt.xlabel("Frequency (Count)")
    plt.ylabel(f"{k}-mer Sequence")
    
    plt.tight_layout()
    output_filename = f"{os.path.splitext(filename)[0]}_top{top_n}_{k}mers.png"
    plt.savefig(output_filename)
    
    print(f"Chart saved as {output_filename}")

if __name__ == "__main__":
    genome_sequence, genome_length = parse_fasta(FASTA_FILE_TO_ANALYZE)
    
    if genome_sequence:
        top_kmers = find_frequent_kmers(genome_sequence, KMER_LENGTH, TOP_N_TO_PLOT)
        
        plot_kmer_frequencies(top_kmers, FASTA_FILE_TO_ANALYZE, KMER_LENGTH, TOP_N_TO_PLOT, genome_length)
        
        print("\nTop k-mers found:")
        # Sort the final dictionary for a clean printout
        for kmer, count in sorted(top_kmers.items(), key=lambda item: item[1], reverse=True):
            print(f"{kmer}: {count}")