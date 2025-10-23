from collections import Counter
import matplotlib.pyplot as plt
import numpy as np

CODON_TABLE = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': 'STOP', 'UAG': 'STOP',
    'UGU': 'C', 'UGC': 'C', 'UGA': 'STOP', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Amino Acid names
AMINO_ACID_NAMES = {
    'F': 'Phenylalanine', 'L': 'Leucine', 'I': 'Isoleucine', 'M': 'Methionine',
    'V': 'Valine', 'S': 'Serine', 'P': 'Proline', 'T': 'Threonine', 'A': 'Alanine',
    'Y': 'Tyrosine', 'STOP': 'STOP', 'H': 'Histidine', 'Q': 'Glutamine', 'N': 'Asparagine',
    'K': 'Lysine', 'D': 'Aspartic acid', 'E': 'Glutamic acid', 'C': 'Cysteine',
    'W': 'Tryptophan', 'R': 'Arginine', 'G': 'Glycine',
}

def parse_fasta(filename):
    sequence = ""
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                sequence += line.strip()
    return sequence

def analyze_genome(sequence):
    rna_sequence = sequence.upper().replace('T', 'U')
    valid_bases = {'A', 'C', 'G', 'U'}
    cleaned_sequence = ''.join([base for base in rna_sequence if base in valid_bases])
    
    codon_counts = Counter()
    for i in range(0, len(cleaned_sequence) - 2, 3):
        codon = cleaned_sequence[i:i+3]
        if len(codon) == 3:
            codon_counts[codon] += 1
            
    amino_acid_counts = Counter()
    for codon, count in codon_counts.items():
        if codon in CODON_TABLE:
            amino_acid = CODON_TABLE[codon]
            amino_acid_counts[amino_acid] += count
            
    return codon_counts, amino_acid_counts

def plot_top_10_codons(codon_counts, title):
    top_10 = codon_counts.most_common(10)
    codons, counts = zip(*top_10)
    
    labels = []
    for codon in codons:
        aa_short = CODON_TABLE.get(codon, '?')
        aa_long = AMINO_ACID_NAMES.get(aa_short, '?')
        labels.append(f"{codon}\n({aa_long})")

    plt.figure(figsize=(12, 7))
    plt.bar(labels, counts, color='skyblue')
    plt.title(title, fontsize=16)
    plt.ylabel("Frequency", fontsize=12)
    plt.xlabel("Codon", fontsize=12)
    plt.tight_layout()
    plt.show()

def plot_aa_comparison(covid_aa, flu_aa, title):
    covid_aa = Counter({aa: c for aa, c in covid_aa.items() if aa != 'STOP'})
    flu_aa = Counter({aa: c for aa, c in flu_aa.items() if aa != 'STOP'})

    covid_top_3 = covid_aa.most_common(3)
    flu_top_3 = flu_aa.most_common(3)
    
    covid_top_3_labels = [aa[0] for aa in covid_top_3]
    flu_top_3_labels = [aa[0] for aa in flu_top_3]
    
    combined_labels = sorted(list(set(covid_top_3_labels + flu_top_3_labels)))
    
    covid_counts = [covid_aa.get(aa, 0) for aa in combined_labels]
    flu_counts = [flu_aa.get(aa, 0) for aa in combined_labels]
    x_labels = [f"{AMINO_ACID_NAMES.get(aa)}\n({aa})" for aa in combined_labels]

    x = np.arange(len(x_labels))
    width = 0.35 

    fig, ax = plt.subplots(figsize=(12, 7))
    rects1 = ax.bar(x - width/2, covid_counts, width, label='COVID-19', color='coral')
    rects2 = ax.bar(x + width/2, flu_counts, width, label='Influenza A', color='deepskyblue')

    ax.set_ylabel('Total Frequency')
    ax.set_title(title, fontsize=16)
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels)
    ax.legend()

    ax.bar_label(rects1, padding=3)
    ax.bar_label(rects2, padding=3)

    fig.tight_layout()
    plt.show()
    
    return [AMINO_ACID_NAMES.get(aa) for aa in covid_top_3_labels]

def print_final_answers(top_3_list, covid_top_3_data, flu_top_3_data):
    
    print("\n" + "="*50)
    print("(D) Top 3 Most Frequent Amino Acids (Console Output)")
    print("="*50)
    
    print("COVID-19 (SARS-CoV-2):")
    for i, (aa, count) in enumerate(covid_top_3_data):
        print(f"  {i+1}. {AMINO_ACID_NAMES[aa]:<16} ({aa}): {count:,}")
        
    print("\nInfluenza A:")
    for i, (aa, count) in enumerate(flu_top_3_data):
        print(f"  {i+1}. {AMINO_ACID_NAMES[aa]:<16} ({aa}): {count:,}")

    aa1, aa2, aa3 = top_3_list[0], top_3_list[1], top_3_list[2]
    
    print("\n" + "="*50)
    print(f"(E) Foods Lacking: {aa1}, {aa2}, and {aa3}")
    print("="*50)
    print(f"""

1.  Fats and Oils:
    * Olive oil, vegetable oil, coconut oil, butter

2.  Sugars and Starches:
    * Table sugar, cornstarch, high-fructose corn syrup

3.  Most Fruits and Non-Starchy Vegetables:
    * Apples, grapes, oranges, melons
    * Cucumbers, lettuce, bell peppers, carrots
""")

if __name__ == "__main__":
    
    COVID_FILE = "covid.fasta"
    INFLUENZA_FILE = "influenza.fasta"
    
    try:
        print(f"Reading {COVID_FILE}...")
        covid_seq = parse_fasta(COVID_FILE)
        covid_codons, covid_aa = analyze_genome(covid_seq)
        
        print(f"Reading {INFLUENZA_FILE}...")
        flu_seq = parse_fasta(INFLUENZA_FILE)
        flu_codons, flu_aa = analyze_genome(flu_seq)

        print("Generating plots...")

        # Plot (A)
        plot_top_10_codons(covid_codons, "(A) COVID-19: Top 10 Codons")
        
        # Plot (B)
        plot_top_10_codons(flu_codons, "(B) Influenza A: Top 10 Codons")
        
        # Plot (C)
        top_3_names_for_e = plot_aa_comparison(
            covid_aa, 
            flu_aa, 
            "(C) Comparison of Top Amino Acid Frequencies"
        )

        covid_aa_non_stop = Counter({aa: c for aa, c in covid_aa.items() if aa != 'STOP'})
        flu_aa_non_stop = Counter({aa: c for aa, c in flu_aa.items() if aa != 'STOP'})

        print_final_answers(
            top_3_names_for_e,
            covid_aa_non_stop.most_common(3),
            flu_aa_non_stop.most_common(3)
        )

    except FileNotFoundError as e:
        print(f"\n--- ERROR ---")
        print(f"File not found: '{e.filename}'")
        print(f"Make sure '{COVID_FILE}' and '{INFLUENZA_FILE}' are in the same folder as this script.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")