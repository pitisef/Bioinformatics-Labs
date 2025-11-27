# Use the AI to design an application that uses the sliding window methodology
# in order to scan a DNA sequence from a FASTA file, and display the melting 
# temperature along the seq., by using a chart. The chart should have 2 signals,
# one for each formula.
# Note: The sliding window should have 9 positions.

import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt
import math

def analyze_and_plot():
    filepath = filedialog.askopenfilename(initialfile="ex3.fasta", filetypes=[("FASTA Files", "*.fasta *.fa *.txt")])
    
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
    
    if len(final_seq) < 9:
        messagebox.showerror("Error", "Sequence is too short (min 9)")
        return
    window_size = 9
    limit = len(final_seq) - window_size + 1
    
    na_conc = 0.001
    
    tm1_list = []
    tm2_list = []
    positions = []

    for i in range(limit):
        window = final_seq[i : i + window_size]
        
        A = 0
        T = 0
        C = 0
        G = 0
        
        for base in window:
            if base == 'A': A += 1
            if base == 'T': T += 1
            if base == 'C': C += 1
            if base == 'G': G += 1
        val1 = 4 * (G + C) + 2 * (A + T)
        gc_percent = ((G + C) / window_size) * 100
        val2 = 81.5 + 16.6 * math.log10(na_conc) + 0.41 * gc_percent - (600 / window_size)
        
        tm1_list.append(val1)
        tm2_list.append(val2)
        positions.append(i)

    plt.figure(figsize=(10, 6))
    plt.plot(positions, tm1_list, label="Basic Tm", color="blue")
    plt.plot(positions, tm2_list, label="Salt Adjusted Tm", color="red")
    
    plt.title("Melting Temperature (Window 9)")
    plt.xlabel("Position")
    plt.ylabel("Temperature (C)")
    plt.legend()
    plt.show()

root = tk.Tk()
root.title("Tm Analysis")
root.geometry("300x150")

label = tk.Label(root, text="Analyze Melting Temp", font=("Arial", 14))
label.pack(pady=20)

btn = tk.Button(root, text="Select FASTA File", command=analyze_and_plot)
btn.pack()

root.mainloop()