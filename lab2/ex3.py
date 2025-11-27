#Design an application using the AI which contains GUI which allows the user 
#select a FASTA file. The content of the file should be analyzed by using a 
#sliding window of 30 positions. The content of each sliding window should be 
#used in order to extract the relative frequencies of the symbols found in the 
#alphabet of the sequence. Thus your input should be the DNA seq from the 
#FASTA file, and the output should be the values of relative frequencies of 
#each symbol in the alphabet of the sequence. Translate in lines on a chart, 
#thus your chart in the case of DNA should have 4 lines which reflect the 
#values found over the sequence

import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt

def analyze_and_plot():
    initialfile="ex3.fasta",
    filepath = filedialog.askopenfilename(filetypes=[("FASTA Files", "*.fasta *.fa *.txt")])
    if not filepath:
        return

    try:
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
        
        if len(final_seq) < 30:
            messagebox.showerror("Error", "Sequence is too short (min 30).")
            return
        window_size = 30
        limit = len(final_seq) - window_size + 1
        
        positions = []
        freq_a = []
        freq_c = []
        freq_g = []
        freq_t = []

        for i in range(limit):
            segment = final_seq[i : i + window_size]
            a = 0
            c = 0
            g = 0
            t = 0
            
            for base in segment:
                if base == 'A': a += 1
                elif base == 'C': c += 1
                elif base == 'G': g += 1
                elif base == 'T': t += 1
            
            freq_a.append(a / window_size)
            freq_c.append(c / window_size)
            freq_g.append(g / window_size)
            freq_t.append(t / window_size)
            positions.append(i)

        plt.figure(figsize=(10, 6))
        plt.plot(positions, freq_a, label="A", color="blue", linewidth=1)
        plt.plot(positions, freq_c, label="C", color="green", linewidth=1)
        plt.plot(positions, freq_g, label="G", color="orange", linewidth=1)
        plt.plot(positions, freq_t, label="T", color="red", linewidth=1)

        plt.title("Nucleotide Frequencies (Window 30)")
        plt.xlabel("Position")
        plt.ylabel("Frequency")
        plt.legend()
        plt.show()

    except Exception as e:
        messagebox.showerror("Error", f"Something went wrong: {e}")

root = tk.Tk()
root.title("DNA Sliding Window")
root.geometry("300x150")

label = tk.Label(root, text="Analyze DNA Frequencies", font=("Arial", 14))
label.pack(pady=20)

btn = tk.Button(root, text="Select FASTA File", command=analyze_and_plot)
btn.pack()

root.mainloop()