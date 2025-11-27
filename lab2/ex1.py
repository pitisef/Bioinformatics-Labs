from itertools import product

def calculate_percentage(S, length):
    freq = {}
    n = len(S) - length + 1
    
    for i in range(n):
        pattern = S[i:i+length]
        freq[pattern] = freq.get(pattern, 0) + 1
        
    # Convert to percentage
    for pattern in freq:
        freq[pattern] = (freq[pattern] / n) * 100
        
    return freq

S = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"
bases = ['A', 'C', 'G', 'T']

dinucleotides = [''.join(p) for p in product(bases, repeat=2)]
trinucleotides = [''.join(p) for p in product(bases, repeat=3)]

print("All Dinucleotides:", dinucleotides)
print("All Trinucleotides:", trinucleotides)

di_stats = calculate_percentage(S, 2)
print("\nDinucleotide Percentages:")
for k in sorted(di_stats):
    print(f"{k}: {di_stats[k]:.2f}%")

tri_stats = calculate_percentage(S, 3)
print("\nTrinucleotide Percentages:")
for k in sorted(tri_stats):
    print(f"{k}: {tri_stats[k]:.2f}%")