import math
import matplotlib.pyplot as plt

def get_counts(seq):
    chars = "ACGT"
    counts = {c: {c2: 0 for c2 in chars} for c in chars}
    sums = {c: 0 for c in chars}
    
    for i in range(len(seq) - 1):
        curr = seq[i]
        next_c = seq[i+1]
        counts[curr][next_c] += 1
        sums[curr] += 1
        
    probs = {c: {c2: 0.0 for c2 in chars} for c in chars}
    for c in chars:
        for c2 in chars:
            if sums[c] > 0:
                probs[c][c2] = counts[c][c2] / sums[c]
            else:
                probs[c][c2] = 0.0
    return probs

def get_llm(pos, neg):
    chars = "ACGT"
    matrix = {c: {c2: 0.0 for c2 in chars} for c in chars}
    eps = 1e-10
    
    for c in chars:
        for c2 in chars:
            p = pos[c][c2]
            n = neg[c][c2]
            matrix[c][c2] = math.log2((p + eps) / (n + eps))
            
    return matrix

def score_seq(seq, matrix):
    s = 0.0
    for i in range(len(seq) - 1):
        c = seq[i]
        n = seq[i+1]
        s += matrix[c][n]
    return s

def plot_heatmap(data, title, ax):
    chars = "ACGT"
    grid = [[data[r][c] for c in chars] for r in chars]
    
    im = ax.imshow(grid, cmap='viridis')
    ax.set_title(title)
    ax.set_xticks(range(4))
    ax.set_yticks(range(4))
    ax.set_xticklabels(list(chars))
    ax.set_yticklabels(list(chars))
    
    for i in range(4):
        for j in range(4):
            val = grid[i][j]
            color = "black"
            if abs(val) > 0.5: 
                 color = "white"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", color=color)

S1 = "ATCGATTCGATATCATACACGTAT"
S2 = "CTCGACTAGTATGAAGTCCACGCTTG"
target = "CAGGTTGGAAACGTAA"

pos_probs = get_counts(S1)
neg_probs = get_counts(S2)
llm = get_llm(pos_probs, neg_probs)
final_score = score_seq(target, llm)

print("Positive Matrix (+):")
for r in "ACGT":
    print(r, [round(pos_probs[r][c], 2) for c in "ACGT"])

print("\nNegative Matrix (-):")
for r in "ACGT":
    print(r, [round(neg_probs[r][c], 2) for c in "ACGT"])

print("\nLog-Likelihood Matrix:")
for r in "ACGT":
    print(r, [round(llm[r][c], 2) for c in "ACGT"])

print("\nSequence:", target)
print("Score:", round(final_score, 4))

if final_score > 0:
    print("Result: Likely CpG Island")
else:
    print("Result: Not CpG Island")

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 5))
plot_heatmap(pos_probs, "Observed (+)", ax1)
plot_heatmap(neg_probs, "Expected (-)", ax2)
plot_heatmap(llm, "Log-Likelihood", ax3)
plt.tight_layout()
plt.show()