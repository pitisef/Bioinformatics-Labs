import math

motifs = [
    "CAGGTTGAG",
    "CATGTTGAG",
    "CAGGTAGAG",
    "CAGGTTGAT",
    "CAGGTTAAG",
    "CAGGTTCAG",
    "CAGGTGGAG",
    "CAGGTTGGG",
    "CAGGTTGCG"
]

S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"

k = len(motifs[0]) 
n = len(motifs)    
background_prob = 0.25 

print(f"Analyzing {n} motifs of length {k}...")
print("-" * 50)
#count matrix
counts = {
    'A': [0] * k,
    'C': [0] * k,
    'G': [0] * k,
    'T': [0] * k
}

for seq in motifs:
    for i in range(k):
        nucleotide = seq[i]
        counts[nucleotide][i] += 1

print("1. COUNT MATRIX")
for base in ['A', 'C', 'G', 'T']:
    print(f"{base}: {counts[base]}")
print("-" * 50)

#weight matrix
pseudocount = 1
total_col = n + (4 * pseudocount)

weights = {'A': [], 'C': [], 'G': [], 'T': []}

for base in ['A', 'C', 'G', 'T']:
    for i in range(k):
        val = (counts[base][i] + pseudocount) / total_col
        weights[base].append(round(val, 4))

print("2. WEIGHT MATRIX (Probabilities)")
for base in ['A', 'C', 'G', 'T']:
    print(f"{base}: {weights[base]}")
print("-" * 50)

#relative frequencies matrix
rel_freqs = {'A': [], 'C': [], 'G': [], 'T': []}

for base in ['A', 'C', 'G', 'T']:
    for i in range(k):
        val = weights[base][i] / background_prob
        rel_freqs[base].append(round(val, 4))

print("3. RELATIVE FREQUENCIES MATRIX")
for base in ['A', 'C', 'G', 'T']:
    print(f"{base}: {rel_freqs[base]}")
print("-" * 50)

#Log-likelihoods matrix
log_matrix = {'A': [], 'C': [], 'G': [], 'T': []}

for base in ['A', 'C', 'G', 'T']:
    for i in range(k):
        val = rel_freqs[base][i]
        log_val = math.log2(val)
        log_matrix[base].append(round(log_val, 4))

print("4. LOG-LIKELIHOODS MATRIX (PWM)")
for base in ['A', 'C', 'G', 'T']:
    print(f"{base}: {log_matrix[base]}")
print("-" * 50)

#analyze sequence s
print(f"5. ANALYZING SEQUENCE S")
print(f"Sequence: {S}")
print(f"Length:   {len(S)}")
print("-" * 50)

limit = len(S) - k + 1
threshold = 0
found_signal = False
best_score = -999
best_pos = -1

print(f"{'Pos':<5} {'Window':<15} {'Score':<10}")

for i in range(limit):
    window = S[i : i + k]
    score = 0
    for j in range(k):
        nuc = window[j]
        if nuc in log_matrix:
            score += log_matrix[nuc][j]
    
    print(f"{i:<5} {window:<15} {score:.4f}")
    
    if score > best_score:
        best_score = score
        best_pos = i
        
    if score > 5.0:
        found_signal = True

print("-" * 50)
print(f"Best Score: {best_score:.4f} at Position {best_pos}")

if found_signal:
    print("\nCONCLUSION: YES, signals indicate an exon-intron border.")
    print(f"it matches the motif profile strongly.")
else:
    print("\nCONCLUSION: NO strong signals found.")