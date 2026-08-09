import random

Dna = ["TTACCTTAAC",
       "GATGTCTGTC",
       "ACGGCGTTAG",
       "CCCTAACGAG",
       "CGTCAGAGGT"]
k = 3
t = 5

def RandomMotifs(Dna, k, t):
    """Return a list of t random k-mers, one from each string in Dna."""
    motifs = []
    n = len(Dna[0])                                    # length of each DNA string
    for i in range(t):                                 # iterate over each string in Dna
        start = random.randint(0, n - k)               # random start position for the k-mer
        motifs.append(Dna[i][start:start + k])         # extract k-mer at random position
    return motifs

print(RandomMotifs(Dna, k, t))
