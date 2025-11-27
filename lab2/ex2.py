from itertools import product

nucleotides = ['A', 'C', 'G', 'T']
dinucleotides = [''.join(p) for p in product(nucleotides, repeat=2)]
trinucleotides = [''.join(p) for p in product(nucleotides, repeat=3)]

print("All Dinucleotides:", dinucleotides)
print("All Trinucleotides:", trinucleotides)
print("-" * 20)

def percentage(S, length):
    freq = {}
    limit = len(S) - length + 1

    for i in range(limit):
        sub = S[i : i + length]
        if sub in freq:
            freq[sub] += 1
        else:
            freq[sub] = 1

    for sub in freq:
        freq[sub] = round((freq[sub] / limit) * 100, 2)

    return freq

S = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"

di_stats = percentage(S, 2)
print("Percentages for ALL dinucleotides:")
for comb in dinucleotides:
    val = 0
    if comb in di_stats:
        val = di_stats[comb]
    print(f"{comb}: {val}%")

print("-" * 20)

tri_stats = percentage(S, 3)
print("Percentages for ALL trinucleotides:")
for comb in trinucleotides:
    val = 0
    if comb in tri_stats:
        val = tri_stats[comb]
    print(f"{comb}: {val}%")