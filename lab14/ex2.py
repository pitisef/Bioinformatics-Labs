# We are inside a courtroom, Mihai is accused of plagiarism.
# The lawyer must proove his innocence, and asks you for your help
# as a specialist in infractions. You take 2 poetries: one of M. Eminescu
# and the other fo N. Stanescu. You use these 2 texts to create transition
# matrices, which is able to capture  the probability of transition
# between word. You combine the matrices into a Log LikelyHood Matrix.
# We use the LLM to scan, by using a sliding window, the text of Mihai, the 
# guy accused of plagiarism.
# negative scores -> one model
# positive score -> other model 
# zero values -> neither
# simulate the test of mihai using AI and ask it to ake it a combination between 2 poetries of Eminescu  & Stanescu
#at the end of this process you should be able to tell which parts of the text 
# in written by Eminescu or by Stanescu or by neither

import math
import matplotlib.pyplot as plt

text_eminescu = """
A fost odata ca-n povesti a fost ca niciodata
Din rude mari imparatesti o prea frumoasa fata
Si era una la parinti si mandra-n toate cele
Cum e Fecioara intre sfinti si luna intre stele
Din umbra falnicelor bolti ea pasul si-l indreapta
Langa fereastra unde-n colt Luceafarul asteapta
Privind in zari cum rasareau a marii stele-n-cheaguri
Carari neulite cadeau pe miscatoarele catarge
Cobori in jos luceafar bland alunecand pe o raza
Patrunde in casa si in gand si viata imi lumineaza
El asculta tremurator se aprindea mai tare
Si se arunca fulgerator se cufunda in mare
"""

text_stanescu = """
Leoaica tanara iubirea mi-a sarit in fata
Ma pandise in incordare mai demult
Colti albi mi-a infipt in fata m-a muscat leoaica azi de fata
Si in jurul unui clopot m-a facut sa ma invart
Si mai zise celui care l-a lovit cu barda-n frunte
Chiar si el a fost lovit de-o leoaica mult mai crunta
Privirea-n sus tasni curcubeu taiat in doua
Si auzul o-ntalni tocmai langa ciocarlie
Mi-am dus mana la spranceana la tampla si la barbie
Dar mana nu mi-o mai stiu
Aluneca in nestire pe un desert in stralucire
Peste care trece-alene o leoaica aramie
Cu miscarile viclene inca-o vreme
Si-inca-o vreme
"""

text_mihai = """
cobori in jos luceafar bland leoaica tanara iubirea 
mi-a sarit in fata alunecand pe o raza, 
yeye, sunt un fraier yeye
"""

def get_tokens(text):
    text = text.lower()
    for char in [".", ",", "!", "?", "-", "\n"]:
        text = text.replace(char, " ")
    tokens = []
    for w in text.split(" "):
        if w.strip() != "":
            tokens.append(w.strip())
    return tokens

tokens_e = get_tokens(text_eminescu)
tokens_s = get_tokens(text_stanescu)
tokens_m = get_tokens(text_mihai)

vocab = sorted(list(set(tokens_e + tokens_s + tokens_m)))

def calculate_probs(tokens, vocabulary):
    counts = {}
    sums = {}
    
    for w in vocabulary:
        counts[w] = {}
        sums[w] = 0
        for w2 in vocabulary:
            counts[w][w2] = 0
            
    for i in range(len(tokens) - 1):
        curr = tokens[i]
        nxt = tokens[i+1]
        counts[curr][nxt] += 1
        sums[curr] += 1
        
    probs = {}
    alpha = 0.1
    v_len = len(vocabulary)
    
    for w in vocabulary:
        probs[w] = {}
        for w2 in vocabulary:
            num = counts[w][w2] + alpha
            den = sums[w] + (alpha * v_len)
            probs[w][w2] = num / den
            
    return probs

probs_e = calculate_probs(tokens_e, vocab)
probs_s = calculate_probs(tokens_s, vocab)

llm = {}
for w in vocab:
    llm[w] = {}
    for w2 in vocab:
        pe = probs_e[w][w2]
        ps = probs_s[w][w2]
        llm[w][w2] = math.log2(pe / ps)

window_size = 3
scores = []
indices = []

for i in range(len(tokens_m) - window_size + 1):
    win = tokens_m[i : i + window_size]
    score = 0
    for j in range(len(win) - 1):
        w1 = win[j]
        w2 = win[j+1]
        score += llm[w1][w2]
    scores.append(score)
    indices.append(i)

print("Pos   Window                       Score    Verdict")
print("-" * 55)

for i, s in zip(indices, scores):
    snippet = " ".join(tokens_m[i:i+window_size])
    verdict = "Neutral"
    if s > 0.5:
        verdict = "Eminescu"
    elif s < -0.5:
        verdict = "Stanescu"
        
    print(f"{i:<5} {snippet:<28} {s:.2f}     {verdict}")

plt.figure(figsize=(10, 6))
plt.plot(indices, scores, marker='o', color='purple', label='Score')
plt.axhline(0, color='black', linestyle='--')

plt.fill_between(indices, scores, 0, where=[s > 0 for s in scores], facecolor='blue', alpha=0.3, label='Eminescu')
plt.fill_between(indices, scores, 0, where=[s < 0 for s in scores], facecolor='red', alpha=0.3, label='Stanescu')

plt.title("Plagiarism Detection Analysis")
plt.xlabel("Window Position")
plt.ylabel("Log-Likelihood Score")
plt.legend()
plt.grid(True, alpha=0.3)

matrix_data = []
for w1 in vocab:
    row = []
    for w2 in vocab:
        row.append(llm[w1][w2])
    matrix_data.append(row)

plt.figure(figsize=(12, 10))
plt.imshow(matrix_data, cmap='bwr_r', aspect='auto')
plt.colorbar(label='Log-Likelihood (Blue=Eminescu, Red=Stanescu)')
plt.title("Log-Likelihood Matrix Heatmap")
plt.xlabel("Next Word Index")
plt.ylabel("Current Word Index")

plt.show()