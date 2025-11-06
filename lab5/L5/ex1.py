import random
import time
import sys
import os
from Bio import SeqIO
sys.setrecursionlimit(3000)
READ_COUNT = 2000         
READ_MIN_LEN = 100
READ_MAX_LEN = 150
MIN_OVERLAP = 10
def load_sequence_from_file(filename="lab5.fasta"):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Error: The file '{filename}' was not found. Please ensure it is in the same directory.")

    print(f"Loading original DNA sequence from {filename}...")
    
    try:
        with open(filename, "r") as handle:
            seq_record = next(SeqIO.parse(handle, "fasta"))
        
        sequence = str(seq_record.seq).upper().replace('N', '').replace('-', '')
        
        seq_id = seq_record.id
        description = seq_record.description

        if not (1000 <= len(sequence) <= 3000):
            print(f"Warning: Sequence length ({len(sequence)}) is outside the recommended 1000-3000 nt range.")
            
        return sequence, seq_id, description
    
    except Exception as e:
        raise IOError(f"Failed to parse FASTA file '{filename}'. Check the file format. Error: {e}")


def generate_reads(sequence, count=READ_COUNT, min_len=READ_MIN_LEN, max_len=READ_MAX_LEN):
    reads = []
    seq_len = len(sequence)
    
    if seq_len < max_len:
         raise ValueError("Original sequence is too short for the specified read length.")

    for _ in range(count):
        read_len = random.randint(min_len, max_len)
        start = random.randint(0, seq_len - read_len)
        
        sample = sequence[start : start + read_len]
        reads.append(sample)
    
    print(f"Generated {len(reads):,} short reads (samples) from the original sequence.")
    return reads
def find_best_overlap(s1, s2, min_overlap=MIN_OVERLAP):
    best_overlap_length = 0
    merged_contig = None
    max_k = min(len(s1), len(s2))
    
    for k in range(min_overlap, max_k + 1):
        if s1[-k:] == s2[:k]:
            if k > best_overlap_length:
                best_overlap_length = k
                merged_contig = s1 + s2[k:]
                
    return best_overlap_length, merged_contig

def greedy_assembly(reads):
    contigs = list(reads)
    merge_operations = 0
    
    print(f"Starting assembly with {len(contigs)} initial reads.")
    while True:
        best_i, best_j = -1, -1
        max_overlap = MIN_OVERLAP - 1
        best_merged = None
        for i in range(len(contigs)):
            for j in range(len(contigs)):
                if i == j:
                    continue 
                
                c1 = contigs[i]
                c2 = contigs[j]
                overlap, merged = find_best_overlap(c1, c2)
                
                if overlap > max_overlap:
                    max_overlap = overlap
                    best_i, best_j = i, j
                    best_merged = merged
        if max_overlap < MIN_OVERLAP:
            break
        indices_to_pop = sorted([best_i, best_j], reverse=True)
        for index in indices_to_pop:
            contigs.pop(index)
            
        contigs.append(best_merged)
        merge_operations += 1
        if merge_operations % 100 == 0 or len(contigs) == 1:
             print(f"Merge #{merge_operations}. Remaining contigs: {len(contigs)}. Max overlap: {max_overlap} bp.")
        
        if len(contigs) == 1:
            break

    if not contigs:
        return ""
    return max(contigs, key=len)


def compare_sequences(original, assembled):
    if not assembled:
        print("Assembly failed: No sequence was reconstructed.")
        return

    def longest_common_substring(s1, s2):
        m = [[0] * (1 + len(s2)) for _ in range(1 + len(s1))]
        longest = 0
        x_longest = 0
        for x in range(1, 1 + len(s1)):
            for y in range(1, 1 + len(s2)):
                if s1[x - 1] == s2[y - 1]:
                    m[x][y] = m[x - 1][y - 1] + 1
                    if m[x][y] > longest:
                        longest = m[x][y]
                        x_longest = x
                else:
                    m[x][y] = 0
        return longest, s1[x_longest - longest : x_longest]

    lcs_length, _ = longest_common_substring(original, assembled)
    
    if original in assembled:
        print("\nAssembly Accuracy: PERFECT MATCH (100% Contiguity)")
    else:
        print("\nAssembly Accuracy: IMPERFECT MATCH")
        print(f"Longest Contiguous Match with Original: {lcs_length:,} nt")
        print(f"Contiguity Percentage: {(lcs_length / len(original)) * 100:.2f}%")

if __name__ == "__main__":
    try:
        original_sequence, seq_id, description = load_sequence_from_file("lab5.fasta")
    except (FileNotFoundError, IOError, ValueError) as e:
        print(e)
        sys.exit(1)

    print("=" * 70)
    print(f"Original Sequence ID: {seq_id}")
    print(f"Description: {description}")
    print(f"Length of Original Sequence: {len(original_sequence):,} nt")
    print("=" * 70)
    reads = generate_reads(original_sequence)
    print("\n--- Starting Greedy Assembly ")
    start_time = time.time()
    reconstructed_seq = greedy_assembly(reads)
    end_time = time.time()
    print("\n--- Final Assembly Results ---")
    print(f"Assembly Time: {end_time - start_time:.2f} seconds")
    print(f"Final Assembled Contig Length: {len(reconstructed_seq):,} nt")
    
    compare_sequences(original_sequence, reconstructed_seq)

    print("The simple greedy algorithm always picks the longest overlap, and if a sequence has long, identical repeats, reads from those different sections look the same. This forces the assembler to merge them wrong, resulting in a scrambled sequence or the process breaking down into fragments. Real assembly needs better tricks, like using graph theory, to map the repeats correctly.")