import random
import matplotlib.pyplot as plt

def nucleotide_percentage(S):
    freq = {}
    n = len(S)
    for char in S:
        freq[char] = freq.get(char, 0) + 1
    for char in freq:
        freq[char] = (freq[char]/n) * 100
    return freq

def cg_percentage(S):
    count = 0
    n = len(S)
    for i in range(n-1):
        if S[i:i+2] == "CG":
            count += 1
    return (count/n) * 100

S = "".join(random.choice("ACGT") for _ in range(200))
print('Sequence:', S)

print('Nucleotide %:', nucleotide_percentage(S))
cg = cg_percentage(S)
print('CG %:', cg)

plt.bar(['CG'], [cg])
plt.show()