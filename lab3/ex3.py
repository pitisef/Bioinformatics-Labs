# Show the minimum & maximum values of 2 signals. Also, allow the user to set
# the trashold(like a filter) that is able to take it into consideration only
# the values above the trashold. These values above the trashold should be 
# shown to the user on a second chart as horizontal bars. Thus the charts of the
# signal that are above the trashold signal are shown as a horizontal line over the sequence
# Wherever the signal is bellow the trashold the chart should show empty space
import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt
import math

def analyze():
    filepath = filedialog.askopenfilename(initialfile="lab3.fasta", filetypes=[("FASTA Files", "*.fasta *.fa *.txt")])
    if not filepath:
        return

    try:
        threshold = float(entry_thresh.get())
    except ValueError:
        messagebox.showerror("Error", "Please enter a valid number for threshold")
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
        messagebox.showerror("Error", "Sequence too short")
        return

    window = 9
    limit = len(final_seq) - window + 1
    na = 0.001
    
    tm1 = []
    tm2 = []
    positions = []

    for i in range(limit):
        sub = final_seq[i : i + window]
        
        A = 0
        T = 0
        C = 0
        G = 0
        for base in sub:
            if base == 'A': A += 1
            elif base == 'T': T += 1
            elif base == 'C': C += 1
            elif base == 'G': G += 1
            
        val1 = 4 * (G + C) + 2 * (A + T)
        
        gc_perc = ((G + C) / window) * 100
        val2 = 81.5 + 16.6 * math.log10(na) + 0.41 * gc_perc - (600 / window)
        
        tm1.append(val1)
        tm2.append(val2)
        positions.append(i)

    print(f"Tm1 Min: {min(tm1)}, Max: {max(tm1)}")
    print(f"Tm2 Min: {min(tm2)}, Max: {max(tm2)}")

    tm1_filtered = []
    tm2_filtered = []
    
    for val in tm1:
        if val > threshold:
            tm1_filtered.append(val)
        else:
            tm1_filtered.append(None)
            
    for val in tm2:
        if val > threshold:
            tm2_filtered.append(val)
        else:
            tm2_filtered.append(None)

    plt.figure(figsize=(10, 8))
    
    plt.subplot(2, 1, 1)
    plt.plot(positions, tm1, label="Basic Tm", color="blue")
    plt.plot(positions, tm2, label="Salt Tm", color="red")
    plt.axhline(y=threshold, color='black', linestyle='--', label="Threshold")
    plt.title("Full Signals")
    plt.ylabel("Tm")
    plt.legend()

    plt.subplot(2, 1, 2)
    plt.plot(positions, tm1_filtered, label="Basic > Threshold", color="blue")
    plt.plot(positions, tm2_filtered, label="Salt > Threshold", color="red")
    plt.title(f"Signals Above Threshold ({threshold})")
    plt.ylabel("Tm")
    plt.xlabel("Position")
    plt.legend()
    
    plt.tight_layout()
    plt.show()

root = tk.Tk()
root.title("Tm Threshold Analysis")
root.geometry("300x200")

label = tk.Label(root, text="Tm Threshold Filter", font=("Arial", 14))
label.pack(pady=10)

frame = tk.Frame(root)
frame.pack(pady=5)
lbl_t = tk.Label(frame, text="Threshold:")
lbl_t.pack(side=tk.LEFT)
entry_thresh = tk.Entry(frame)
entry_thresh.pack(side=tk.LEFT)

btn = tk.Button(root, text="Select File & Analyze", command=analyze)
btn.pack(pady=20)

root.mainloop()