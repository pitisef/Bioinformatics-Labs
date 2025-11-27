import random

def generate_dna(length):
    bases = "ACGT"
    seq = ""
    for i in range(length):
        seq += random.choice(bases)
    return seq

S = generate_dna(200)
print('The generated sequence is:', S)