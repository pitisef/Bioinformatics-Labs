#DOWNLOAD 10 INFLUENZA virus variants from NCBI and analyze their genome by using the application from assignment 1 (previous assignment)
#make an electrophoresys gel for each genome
#eliminate all lines that are in common between the gel simul such that the differences will be shown
#merge all electrophoresys gel simulations (that show only the differences) in one general electrophoresys gel
import tkinter as tk
from tkinter import messagebox
import matplotlib.pyplot as plt
import os

enzymes = {
    "EcoRI": ("GAATTC", 1),
    "BamHI": ("GGATCC", 1),
    "HindIII": ("AAGCTT", 1),
    "TaqI": ("TCGA", 1),
    "HaeIII": ("GGCC", 2)
}

# Hardcoded list of files to import automatically
TARGET_FILES = [
    "NC_026433.1.fasta",
    "escherichia_coli.fasta",
    "lab9.fasta",
    "MN908947.3.fasta",
    "NC_001477.1.fasta",
    "NC_001498.1.fasta",
    "NC_001542.1.fasta",
    "NC_002058.3.fasta",
    "NC_002200.1.fasta",
    "NC_002549.1.fasta"
]

def digest_sequence(sequence):
    results = {}
    seq_len = len(sequence)
    
    for name, (pattern, offset) in enzymes.items():
        cuts = [0]
        i = sequence.find(pattern)
        while i != -1:
            cut_pos = i + offset
            cuts.append(cut_pos)
            i = sequence.find(pattern, i + 1)
        cuts.append(seq_len)
        cuts.sort()
        
        fragments = []
        for k in range(len(cuts) - 1):
            length = cuts[k+1] - cuts[k]
            fragments.append(length)
            
        results[name] = fragments
    return results

def load_specific_files():
    sequences = {}
    
    print("Attempting to load files...")
    
    for filename in TARGET_FILES:
        try:
            with open(filename, "r") as f:
                header = None
                seq_parts = []
                for line in f:
                    line = line.strip()
                    if not line: continue
                    if line.startswith(">"):
                        # Use the filename or the first word of header as key
                        header = filename 
                        break # We only need the header line to start
                
                # Reset file pointer to read full content properly if needed, 
                # but standard FASTA parsing:
                f.seek(0)
                current_seq = []
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        continue
                    current_seq.append(line.upper())
                
                full_seq = "".join(current_seq)
                clean_seq = "".join([c for c in full_seq if c in "ACGT"])
                
                if clean_seq:
                    sequences[filename] = clean_seq
                    print(f"Successfully loaded: {filename} ({len(clean_seq)} bp)")
                else:
                    print(f"Warning: {filename} was empty or invalid.")

        except FileNotFoundError:
            print(f"Error: Could not find {filename}")
    
    return sequences

def analyze():
    # Load the specific files defined in TARGET_FILES
    genomes = load_specific_files()
    
    if not genomes:
        messagebox.showerror("Error", "No valid sequences loaded.\nMake sure the FASTA files are in the same folder.")
        return

    print(f"\nTotal genomes loaded: {len(genomes)}")

    all_data = {}
    for name, seq in genomes.items():
        all_data[name] = digest_sequence(seq)

    common_map = {}
    
    # 1. Identify common fragments
    first_genome = list(all_data.keys())[0]
    
    for enz in enzymes:
        common_frags = set(all_data[first_genome][enz])
        
        for name in genomes:
            frags = set(all_data[name][enz])
            # Intersection: keep only fragments present in ALL genomes
            common_frags = common_frags.intersection(frags)
        
        common_map[enz] = common_frags
        
    # 2. Filter data
    filtered_data = {}
    for name in genomes:
        filtered_data[name] = {}
        for enz in enzymes:
            original = all_data[name][enz]
            commons = common_map[enz]
            unique = [f for f in original if f not in commons]
            filtered_data[name][enz] = unique

    # 3. Plot
    plt.figure(figsize=(16, 9))
    
    variant_names = list(genomes.keys())
    num_vars = len(variant_names)
    gap = 3
    current_x = 0
    
    # Colors for different files to distinguish them better
    colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'purple', 'orange', 'grey']
    
    for enz in enzymes:
        # Label the enzyme group
        group_center = current_x + (num_vars / 2)
        plt.text(group_center, -100, enz, ha='center', fontsize=12, fontweight='bold')
        
        for i, v_name in enumerate(variant_names):
            frags = filtered_data[v_name][enz]
            lane_pos = current_x + i
            
            # Pick color cyclically
            col = colors[i % len(colors)]
            
            for f in frags:
                plt.hlines(y=f, xmin=lane_pos-0.4, xmax=lane_pos+0.4, linewidth=1.5, color=col)
        
        current_x += num_vars + gap

    # Create a custom legend for the files
    patches = [plt.Line2D([0], [0], color=colors[i % len(colors)], lw=4, label=v_name) for i, v_name in enumerate(variant_names)]
    plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc='upper left', title="Genomes")

    plt.title("Comparative Restriction Analysis (Unique Fragments Only)")
    plt.ylabel("Fragment Size (bp)")
    
    plt.gca().invert_yaxis()
    
    # Remove X axis ticks
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    
    plt.grid(axis='y', linestyle=':', alpha=0.3)
    plt.tight_layout()
    plt.show()

root = tk.Tk()
root.title(" viral Genome Analysis")
root.geometry("350x150")

label = tk.Label(root, text="Analysis of 10 Pre-loaded Genomes", font=("Arial", 10))
label.pack(pady=10)

btn = tk.Button(root, text="Run Analysis", command=analyze, font=("Arial", 12), bg="#dddddd")
btn.pack(expand=True, pady=10)

root.mainloop()