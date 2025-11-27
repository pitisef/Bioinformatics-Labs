# The melting temperature (Tm) is the temp at which 1/2 of a particular
# DNA seq. will dissociate and become a single strand of DNA. Primer length & seq.
# are of critical importance in the design of the params. over successful amplification
# (PCR). The melting temp. of a nucleic adic duplex increases both wit its length
# A simple formula used for the calculation of tm is 
# Tm = 4 * (G + C) + 2 * (A + T)
# The actual Tm is influenced by the concentration of Mg2+ , K+ , and cosolvents. An alternative formula is:
# Tm = 81.5 + 16.6(log10([Na+])) + .41*(%GC) – 600/length
# where Na+ is the concentration of the solution and has a value of 0.001 
# Input 6-12 letters
import math

seq = input("Enter a DNA sequence (6-12 bases): ")
seq = seq.upper()

if len(seq) < 6 or len(seq) > 12:
    print("Invalid length. Please enter 6-12 letters.")
else:
    A = 0
    T = 0
    C = 0
    G = 0
    
    for base in seq:
        if base == 'A': A += 1
        if base == 'T': T += 1
        if base == 'C': C += 1
        if base == 'G': G += 1
        
    tm1 = 4 * (G + C) + 2 * (A + T)
    
    length = len(seq)
    gc_percent = ((G + C) / length) * 100
    na = 0.001
    
    tm2 = 81.5 + 16.6 * math.log10(na) + 0.41 * gc_percent - (600 / length)
    
    print("Tm Formula 1:", tm1)
    print("Tm Formula 2:", round(tm2, 2))