GENETIC_CODE = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
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

def transcribe_dna_to_rna(dna_sequence):
    """Converts a DNA string into an RNA string by replacing T with U."""
    return dna_sequence.upper().replace('T', 'U')

def translate_rna_to_protein(rna_sequence):
    """
    Translates an RNA sequence into a protein chain, starting from the
    first 'AUG' (START) codon and ending at a STOP codon.
    """
    
    start_index = rna_sequence.find('AUG')
    
    if start_index == -1:
        return "Error: No START codon (AUG) found in sequence."
        
    protein_chain = []
    current_position = start_index
    
    while current_position + 3 <= len(rna_sequence):
        codon = rna_sequence[current_position : current_position + 3]
        
        amino_acid = GENETIC_CODE.get(codon, 'X')
        

        if amino_acid == '*':
            break
            
        protein_chain.append(amino_acid)
        
        current_position += 3
        
    return "".join(protein_chain)

def main():
    """
    Main function to run the transcription and translation.
    """
    dna_strand = "CGTATTATGGCTTCGACGAGTAATAGCAT"
    
    rna_strand = transcribe_dna_to_rna(dna_strand)
    
    protein_sequence = translate_rna_to_protein(rna_strand)
    
    print(f"Original DNA: {dna_strand}")
    print(f"Transcribed RNA: {rna_strand}")
    print(f"Translated Protein: {protein_sequence}")
    
    print("-" * 20)
    
    dna_strand_2 = "AATGGCTTCGACG"
    rna_strand_2 = transcribe_dna_to_rna(dna_strand_2)
    protein_sequence_2 = translate_rna_to_protein(rna_strand_2)
    print(f"Original DNA: {dna_strand_2}")
    print(f"Transcribed RNA: {rna_strand_2}")
    print(f"Translated Protein: {protein_sequence_2}")


if __name__ == "__main__":
    main()
