import sys

def get_rev_comp(seq):
    mapping = str.maketrans("ACGT", "TGCA")
    return seq.translate(mapping)[::-1]

def read_fasta(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    seq = ""
    for line in lines:
        if not line.startswith(">"):
            seq += line.strip()
    return seq

def detect_unknown_tes(dna_sequence):
    found_elements = []
    n = len(dna_sequence)
    
    min_ir = 4
    max_ir = 6
    
    min_dist = 50 
    max_dist = 2000

    for i in range(n - min_dist):
        for k in range(min_ir, max_ir + 1):
            if i + k > n:
                break
            
            start_ir = dna_sequence[i : i+k]
            target_end_ir = get_rev_comp(start_ir)
            
            search_start = i + k + min_dist
            search_end = min(i + k + max_dist, n)
            
            window = dna_sequence[search_start : search_end]
            
            offset = 0
            while True:
                match_pos = window.find(target_end_ir, offset)
                if match_pos == -1:
                    break
                
                total_end_index = search_start + match_pos + k
                
                found_elements.append((i, total_end_index, total_end_index - i))
                
                offset = match_pos + 1

    unique_elements = sorted(list(set(found_elements)))
    return unique_elements

filenames = ["genome1.fasta", "genome2.fasta", "genome3.fasta"]

for fname in filenames:
    try:
        print(f"Processing {fname}...")
        sequence = read_fasta(fname)
        print(f"Genome Length: {len(sequence)} bp")
        
        results = detect_unknown_tes(sequence)
        
        print(f"Found {len(results)} potential candidates.")
        print(f"{'Start':<10} {'End':<10} {'Length':<10}")
        
        limit = 20 
        for idx, (s, e, l) in enumerate(results):
            if idx >= limit: 
                print("... (more results hidden) ...")
                break
            print(f"{s:<10} {e:<10} {l:<10}")
        print("-" * 30)
        
    except FileNotFoundError:
        print(f"Error: Could not find file {fname}")