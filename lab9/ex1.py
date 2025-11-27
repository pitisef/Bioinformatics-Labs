import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt

enzymes = {
    "EcoRI": ("GAATTC", 1),
    "BamHI": ("GGATCC", 1),
    "HindIII": ("AAGCTT", 1),
    "TaqI": ("TCGA", 1),
    "HaeIII": ("GGCC", 2)
}

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
            
        results[name] = {
            "cuts": len(cuts) - 2,
            "positions": cuts[1:-1],
            "fragments": fragments
        }
    return results

def analyze():
    filepath = filedialog.askopenfilename(initialfile="lab9.fasta", filetypes=[("FASTA Files", "*.fasta *.fa *.txt")])
    if not filepath:
        return

    with open(filepath, "r") as f:
        lines = f.readlines()
    
    sequence = ""
    for line in lines:
        line = line.strip()
        if not line.startswith(">"):
            sequence += line.upper()

    final_seq = ""
    for char in sequence:
        if char in "ACGT":
            final_seq += char
            
    if len(final_seq) < 100:
        messagebox.showwarning("Warning", "Sequence might be too short for digestion.")

    data = digest_sequence(final_seq)
    
    print("-" * 40)
    print(f"Sequence Length: {len(final_seq)} bp")
    print("-" * 40)
    
    plt.figure(figsize=(10, 6))
    
    x_positions = range(len(data))
    names = list(data.keys())
    
    for idx, name in enumerate(names):
        info = data[name]
        print(f"Enzyme: {name}")
        print(f"Cuts: {info['cuts']}")
        print(f"Positions: {info['positions']}")
        print(f"Fragments: {info['fragments']}")
        print("-" * 20)
        
        for frag_len in info['fragments']:
            plt.hlines(y=frag_len, xmin=idx-0.3, xmax=idx+0.3, linewidth=3, color='black')

    plt.xticks(x_positions, names)
    plt.ylabel("Fragment Size (bp)")
    plt.title("Restriction Enzyme Digestion (Gel Simulation)")
    plt.gca().invert_yaxis()
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.show()

root = tk.Tk()
root.title("Restriction Enzyme Digest")
root.geometry("300x150")

btn = tk.Button(root, text="Load FASTA & Simulate Gel", command=analyze, font=("Arial", 12))
btn.pack(expand=True)

root.mainloop()