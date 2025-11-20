import time
import sys
import os

sys.setrecursionlimit(3000)

MIN_UNIT_LENGTH = 3
MAX_UNIT_LENGTH = 6
MIN_REPETITIONS = 2

def load_sequence_from_file(filename="escherichia_coli.fasta"):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File '{filename}' not found.")

    sequence_lines = []
    
    with open(filename, 'r') as f:
        lines = f.readlines()
        
        if not lines or not lines[0].startswith('>'):
            raise IOError("Invalid FASTA format.")
        
        header = lines[0].strip()
        sequence_lines = [line.strip() for line in lines[1:] if line.strip()]
            
    sequence = "".join(sequence_lines).upper().replace('N', '').replace('-', '')
        
    seq_id = header.split()[0].replace('>', '')
    description = ' '.join(header.split()[1:])
    if not description:
        description = "Unknown fragment"

    return sequence, seq_id, description

def find_tandem_repeats(sequence):
    results = {}
    seq_len = len(sequence)

    for unit_len in range(MIN_UNIT_LENGTH, MAX_UNIT_LENGTH + 1):
        for i in range(seq_len - unit_len):
            unit = sequence[i : i + unit_len]
            current_pos = i + unit_len
            repeat_count = 1
            
            while True:
                next_unit = sequence[current_pos : current_pos + unit_len]
                
                if next_unit == unit:
                    repeat_count += 1
                    current_pos += unit_len
                else:
                    break
            
            if repeat_count >= MIN_REPETITIONS:
                total_length = repeat_count * unit_len
                
                key = (i, total_length)
                
                results[key] = {
                    "start": i + 1,
                    "unit": unit,
                    "length": unit_len,
                    "count": repeat_count,
                    "end": i + total_length
                }

    return sorted(results.values(), key=lambda x: x['start'])

if __name__ == "__main__":
    
    try:
        sequence, seq_id, description = load_sequence_from_file("escherichia_coli.fasta")
        
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except IOError as e:
        print(f"Error parsing file: {e}")
        sys.exit(1)

    print(f"Target: {description} ({seq_id})")
    print(f"Sequence Length: {len(sequence):,} bases")
    print("-" * 60)
    
    repeats = find_tandem_repeats(sequence)
    
    print(f"Found {len(repeats)} repetition regions:")
    
    if not repeats:
        print("No repeats found matching the unit size (3-6 bases, N >= 2).")
    else:
        print(f"{'Start':<8}{'End':<8}{'Unit Size':<10}{'Count':<8}{'Pattern'}")
        print("-" * 60)
        for r in repeats:
            total_pattern = r['unit'] * r['count']
            print(f"{r['start']:<8}{r['end']:<8}{r['length']:<10}{r['count']:<8}{total_pattern}")
        print("-" * 60)