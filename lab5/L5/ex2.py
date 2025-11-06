import random
import time
import sys
from collections import Counter

sys.setrecursionlimit(3000)

READ_COUNT = 200
READ_MIN_LEN = 100
READ_MAX_LEN = 150
MIN_OVERLAP = 10

def generate_mock_sequence(length, gc_percent, num_repeats=0, repeat_unit="ATGCATGCATGCATGC"):
    gc_count = int(length * (gc_percent / 100.0))
    at_count = length - gc_count
    
    g_count = c_count = gc_count // 2
    a_count = t_count = at_count // 2
    
    bases = (['G'] * g_count + ['C'] * c_count + ['A'] * a_count + ['T'] * t_count)
    random.shuffle(bases)
    base_sequence = "".join(bases)
    
    if num_repeats > 0:
        for _ in range(num_repeats):
            pos = random.randint(0, len(base_sequence) - len(repeat_unit))
            base_sequence = base_sequence[:pos] + repeat_unit + base_sequence[pos:]
            
        base_sequence = base_sequence[:length]

    return base_sequence

VIRAL_GENOMES = [
    {
        "id": "NC_003977", "name": "Hepatitis B Virus (HBV)", "length": 3215,
        "gc_content": 48.0, "repeats": 5, "sequence": generate_mock_sequence(3215, 48.0, 5)
    },
    {
        "id": "NC_001401", "name": "Adeno-associated Virus (AAV)", "length": 4680,
        "gc_content": 40.0, "repeats": 2, "sequence": generate_mock_sequence(4680, 40.0, 2)
    },
    {
        "id": "NC_001416", "name": "Lambda Phage (Fragment)", "length": 4800,
        "gc_content": 50.0, "repeats": 10, "sequence": generate_mock_sequence(4800, 50.0, 10, "GATTACA")
    },
    {
        "id": "K03455", "name": "HIV-1 (complete)", "length": 9719,
        "gc_content": 42.0, "repeats": 3, "sequence": generate_mock_sequence(9719, 42.0, 3)
    },
    {
        "id": "NC_039401", "name": "Zika Virus (ZIKV)", "length": 10860,
        "gc_content": 48.5, "repeats": 1, "sequence": generate_mock_sequence(10860, 48.5, 1)
    },
    {
        "id": "K02718", "name": "Human Papillomavirus 16 (HPV)", "length": 7904,
        "gc_content": 49.0, "repeats": 4, "sequence": generate_mock_sequence(7904, 49.0, 4)
    },
    {
        "id": "V01149", "name": "Poliovirus", "length": 7441,
        "gc_content": 46.0, "repeats": 2, "sequence": generate_mock_sequence(7441, 46.0, 2)
    },
    {
        "id": "M38515", "name": "Canine Parvovirus (CPV)", "length": 5200,
        "gc_content": 40.5, "repeats": 6, "sequence": generate_mock_sequence(5200, 40.5, 6, "TAGTAGTAGTAG")
    },
    {
        "id": "J02400", "name": "SV40 (complete)", "length": 5243,
        "gc_content": 41.0, "repeats": 1, "sequence": generate_mock_sequence(5243, 41.0, 1)
    },
    {
        "id": "NC_001422", "name": "Phi-X174 Phage", "length": 5386,
        "gc_content": 44.0, "repeats": 8, "sequence": generate_mock_sequence(5386, 44.0, 8, "GGGGCCCCCC")
    }
]

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
        
        if len(contigs) == 1:
            break

    if not contigs:
        return ""
    return max(contigs, key=len)

def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    total_bases = len(sequence)
    if total_bases == 0:
        return 0.0
    return (gc_count / total_bases) * 100

def run_assembly_analysis(genome_data):
    
    original_sequence = genome_data['sequence']
    
    gc_percent = calculate_gc_content(original_sequence)
    
    reads = generate_reads(original_sequence)
    
    start_time = time.perf_counter()
    reconstructed_seq = greedy_assembly(reads)
    end_time = time.perf_counter()
    
    assembly_time_ms = (end_time - start_time) * 1000
    
    if original_sequence in reconstructed_seq:
        contiguity_status = "PERFECT"
    else:
        contiguity_status = "FRAGMENTED" 

    return {
        "id": genome_data['id'],
        "name": genome_data['name'],
        "gc_percent": round(gc_percent, 2),
        "assembly_time_ms": round(assembly_time_ms, 2),
        "assembled_length": len(reconstructed_seq),
        "original_length": len(original_sequence),
        "contiguity": contiguity_status
    }

if __name__ == "__main__":
    
    results = []

    for virus in VIRAL_GENOMES:
        results.append(run_assembly_analysis(virus))

    print(f"{'ID':<10} | {'Virus Name':<30} | {'GC % (X)':>8} | {'Time (ms) (Y)':>12} | {'Original Length':>15} | {'Assembly Status':>17}")
    print("-" * 70)
    
    for r in results:
        print(f"{r['id']:<10} | {r['name']:<30} | {r['gc_percent']:>8.2f} | {r['assembly_time_ms']:>12.2f} | {r['original_length']:>15,} | {r['contiguity']:>17}")