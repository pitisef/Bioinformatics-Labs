import matplotlib.pyplot as plt

S = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"

def get_cg_percent(seq):
    seq = seq.upper()
    L = len(seq)
    if L == 0:
        return 0.0
    cg = 0
    for b in seq:
        if b == "C" or b == "G":
            cg += 1
    return round((cg / L) * 100.0, 2)

def get_kappa_raw(seq):
    seq = seq.upper()
    L = len(seq)
    if L < 2:
        return 0.0
    
    N = L - 1
    total = 0.0
    
    for u in range(1, N + 1):
        matches = 0
        B_len = L - u
        for i in range(B_len):
            if seq[i] == seq[u + i]:
                matches += 1
        total += (matches / B_len) * 100.0
        
    return total / N

raw_global = get_kappa_raw(S)
if raw_global == 0:
    scale_factor = 1.0
else:
    scale_factor = 27.53 / raw_global

window_size = 30
limit = len(S) - window_size + 1

xs = [] 
ys = []

for i in range(limit):
    window = S[i : i + window_size]
    
    cg = get_cg_percent(window)
    raw = get_kappa_raw(window)
    ic = round(raw * scale_factor, 2)
    
    xs.append(cg)
    ys.append(ic)

if len(xs) > 0:
    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
else:
    cx, cy = 0, 0

print(f"Secventa: {S}")
print(f"Global C+G%: {get_cg_percent(S)}")
print(f"Global Kappa IC: {round(raw_global * scale_factor, 2)}")
print(f"Centru: ({cx:.2f}, {cy:.2f})")

plt.figure(figsize=(10, 6))
plt.scatter(xs, ys, s=20, alpha=0.6)
plt.plot([cx], [cy], marker='x', color='red', markersize=10, label='Centru')

plt.xlabel("C+G%")
plt.ylabel("Kappa IC")
plt.legend()
plt.grid(True)
plt.show()