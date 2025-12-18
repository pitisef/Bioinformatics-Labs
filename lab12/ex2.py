import math
import matplotlib.pyplot as plt

motifs = [
    "CAGGTTGAG", "CATGTTGAG", "CAGGTAGAG", "CAGGTTGAT", "CAGGTTAAG",
    "CAGGTTCAG", "CAGGTGGAG", "CAGGTTGGG", "CAGGTTGCG"
]

k = len(motifs[0])
n = len(motifs)
pseudocount = 1
background = 0.25

counts = {'A': [0]*k, 'C': [0]*k, 'G': [0]*k, 'T': [0]*k}
for seq in motifs:
    for i in range(k):
        counts[seq[i]][i] += 1

log_matrix = {'A': [], 'C': [], 'G': [], 'T': []}
total_col = n + (4 * pseudocount)

for base in ['A', 'C', 'G', 'T']:
    for i in range(k):
        weight = (counts[base][i] + pseudocount) / total_col
        rel_freq = weight / background
        log_matrix[base].append(math.log2(rel_freq))

TARGET_FILES = [
    "NC_026433",
    "NC_026434",
    "NC_026435",
    "NC_026436",
    "NC_026437",
    "NC_026438",
    "NC_026432",
    "NC_026431",
    "NC_007373",
    "NC_007366"
]

def read_fasta(filename):
    seq = []
    # Try opening file directly or with .fasta extension
    possible_names = [filename, filename + ".fasta", filename + ".fa"]
    
    file_found = False
    for name in possible_names:
        try:
            with open(name, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line.startswith(">"):
                        seq.append(line.upper())
            file_found = True
            break
        except FileNotFoundError:
            continue
            
    if not file_found:
        print(f"File not found: {filename}")
        return ""
        
    return "".join(seq)

def scan_genome(sequence, matrix, k_len):
    scores = []
    positions = []
    
    limit = len(sequence) - k_len + 1
    
    for i in range(limit):
        window = sequence[i : i + k_len]
        score = 0
        valid_window = True
        
        for j in range(k_len):
            nuc = window[j]
            if nuc in matrix:
                score += matrix[nuc][j]
            else:
                valid_window = False
                break
        
        if valid_window:
            scores.append(score)
            positions.append(i)
            
    return positions, scores

print(f"Scanning {len(TARGET_FILES)} genomes for motif matches...")

for filename in TARGET_FILES:
    dna_seq = read_fasta(filename)
    
    if not dna_seq:
        continue
        
    print(f"Processing {filename} (Length: {len(dna_seq)} bp)...")
    
    dna_seq = "".join([b for b in dna_seq if b in "ACGT"])
    
    pos, val = scan_genome(dna_seq, log_matrix, k)
    
    if not pos:
        print(f"  No valid windows found in {filename}.")
        continue

    threshold = 5.0
    hits_x = []
    hits_y = []
    for p, s in zip(pos, val):
        if s > threshold:
            hits_x.append(p)
            hits_y.append(s)

    plt.figure(figsize=(12, 4))
    plt.plot(pos, val, color='gray', linewidth=0.5, label='Score Signal')
    
    if hits_x:
        plt.scatter(hits_x, hits_y, color='red', s=10, label=f'Likely Motifs (>{threshold})')
    
    plt.title(f"Motif Scan: {filename}")
    plt.xlabel("Position (bp)")
    plt.ylabel("Log-Likelihood Score")
    plt.axhline(y=0, color='blue', linestyle='--', linewidth=0.8, label="Random Baseline")
    plt.legend(loc='upper right')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()