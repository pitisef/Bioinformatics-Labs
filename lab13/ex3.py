#2. use a random english text of about 300 letters(that implies spaces, punctuation) and compute the transition probabilities between words. store the transition matrix as a JSON file. For ease of implementation, you can represent each new word by using a symbol of your choice (ASCII)
import random
import json
import string

def generate_text(length):
    words = [
        "the", "quick", "brown", "fox", "jumps", "over", "lazy", "dog", 
        "hello", "world", "AI", "code", "python", "is", "great", 
        "run", "fast", "eat", "sleep", "repeat", "sun", "moon", "stars", 
        "sky", "blue", "red", "green", "yellow", "cat", "mouse", 
        "computer", "keyboard", "screen", "algorithm", "data", "science",
        "learning", "neural", "network", "matrix", "vector", "analysis",
        "bioinformatics", "gene", "protein", "cell", "life", "future",
        "coding", "programming", "developer", "system", "web", "internet"
    ]
    punctuation = [".", ",", "!", "?"]
    
    text = ""
    while len(text) < length:
        w = random.choice(words)
        if random.random() < 0.2:
            w += random.choice(punctuation)
        text += w + " "
        
    return text[:length]

def get_symbols(text):
    tokens = text.split()
    unique = sorted(list(set(tokens)))
    
    mapping = {}
    reverse_mapping = {}
    start = 33 
    
    for i, w in enumerate(unique):
        code = start + (i % 94) 
        mapping[w] = chr(code)
        reverse_mapping[chr(code)] = w
        
    return tokens, mapping, reverse_mapping

def compute_transitions(tokens, mapping):
    unique = sorted(list(set(tokens)))
    
    counts = {w: {w2: 0 for w2 in unique} for w in unique}
    totals = {w: 0 for w in unique}
    
    for i in range(len(tokens) - 1):
        curr = tokens[i]
        nxt = tokens[i+1]
        counts[curr][nxt] += 1
        totals[curr] += 1
        
    matrix = {}
    
    for w1 in unique:
        sym1 = mapping[w1]
        matrix[sym1] = {}
        for w2 in unique:
            sym2 = mapping[w2]
            if totals[w1] > 0:
                prob = counts[w1][w2] / totals[w1]
                matrix[sym1][sym2] = round(prob, 2)
            else:
                matrix[sym1][sym2] = 0.0
                
    return matrix

def synthesize_text(matrix, reverse_map, length=20):
    if not matrix:
        return ""
    
    current_sym = random.choice(list(matrix.keys()))
    result = [reverse_map[current_sym]]
    
    for _ in range(length - 1):
        probs = matrix[current_sym]
        
        candidates = []
        weights = []
        
        for sym, p in probs.items():
            if p > 0:
                candidates.append(sym)
                weights.append(p)
        
        if not candidates:
            current_sym = random.choice(list(matrix.keys()))
        else:
            current_sym = random.choices(candidates, weights=weights, k=1)[0]
            
        result.append(reverse_map[current_sym])
        
    return " ".join(result)

text = generate_text(300)
print(f"Original Text ({len(text)} chars): {text}")

tokens, maps, rev_maps = get_symbols(text)
print("\nSymbol Mapping:")
print(maps)

result = compute_transitions(tokens, maps)
print("\nTransition Matrix:")
print(result)

with open('word_transition_matrix.json', 'w') as f:
    json.dump(result, f, indent=4)

print("\nSynthesized Sequence:")
print(synthesize_text(result, rev_maps, 20))