#use a random dna sequence of about 50 letters. use this sequence to compute the transition probability between letters. your output should be the transition matrix, stored as a JSON file.
import random
import json

def generate_sequence(length):
    bases = ['A', 'C', 'G', 'T']
    seq = ""
    for i in range(length):
        seq += random.choice(bases)
    return seq

def compute_transitions(seq):
    counts = {'A': {'A':0, 'C':0, 'G':0, 'T':0},
              'C': {'A':0, 'C':0, 'G':0, 'T':0},
              'G': {'A':0, 'C':0, 'G':0, 'T':0},
              'T': {'A':0, 'C':0, 'G':0, 'T':0}}
    
    total = {'A':0, 'C':0, 'G':0, 'T':0}
    
    for i in range(len(seq) - 1):
        current = seq[i]
        next_base = seq[i+1]
        counts[current][next_base] += 1
        total[current] += 1
        
    matrix = {'A': {}, 'C': {}, 'G': {}, 'T': {}}
    
    for base1 in ['A', 'C', 'G', 'T']:
        for base2 in ['A', 'C', 'G', 'T']:
            if total[base1] > 0:
                prob = counts[base1][base2] / total[base1]
                matrix[base1][base2] = round(prob, 2)
            else:
                matrix[base1][base2] = 0.0
                
    return matrix

random_len = random.randint(50, 100)
sequence = generate_sequence(random_len)

transition_matrix = compute_transitions(sequence)
print("transition matrix:")
print(transition_matrix)

with open('transition_matrix.json', 'w') as f:
    json.dump(transition_matrix, f, indent=4)