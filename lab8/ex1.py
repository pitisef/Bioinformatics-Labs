import random

def get_rev_comp(seq):
    mapping = str.maketrans("ACGT", "TGCA")
    return seq.translate(mapping)[::-1]

def generate_dna(length):
    return ''.join(random.choice('ACGT') for _ in range(length))

def insert_te(dna, transposon, ir):
    seq_list = list(dna)
    pos = random.randint(0, len(seq_list))
    
    ir_rev = get_rev_comp(ir)
    element = ir + transposon + ir_rev
    
    seq_list.insert(pos, element)
    return ''.join(seq_list), pos, pos + len(element)

def scan_for_tes(dna, ir_seq):
    found = []
    n = len(dna)
    ir_len = len(ir_seq)
    ir_rev = get_rev_comp(ir_seq)
    
    i = 0
    while i < n - ir_len:
        if dna[i:i+ir_len] == ir_seq:
            start = i
            for j in range(i + ir_len, n - ir_len + 1):
                if dna[j:j+ir_len] == ir_rev:
                    end = j + ir_len
                    found.append((start, end, dna[start:end]))
                    i = end - 1 
                    break
        i += 1
    return found

base_len = 300
ir_marker = "GTCA"
te_variants = ["ATGCGTTT", "CCGGTTAA", "TGCAACGT"]

dna_seq = generate_dna(base_len)
print(f"Original Sequence Length: {len(dna_seq)}")

print("\n--- Inserting Elements ---")
ground_truth = []
for val in te_variants:
    dna_seq, s, e = insert_te(dna_seq, val, ir_marker)
    ground_truth.append((s, e))
    print(f"Inserted at: {s} - {e}")

print(f"Final Sequence Length: {len(dna_seq)}")

print("\n--- Running Detection ---")
detected = scan_for_tes(dna_seq, ir_marker)

for idx, (s, e, seq) in enumerate(detected):
    print(f"TE {idx+1} detected: {s} - {e} | Len: {len(seq)}")

if len(detected) == len(te_variants):
    print("\nAll elements recovered successfully.")
else:
    print(f"\nMismatch: Expected {len(te_variants)}, found {len(detected)}")