text = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
k = 5
profile = {"A": [0.2, 0.2, 0.3, 0.2, 0.3],
           "C": [0.4, 0.3, 0.1, 0.5, 0.1],
           "G": [0.3, 0.3, 0.5, 0.2, 0.4],
           "T": [0.1, 0.2, 0.1, 0.1, 0.2]}

def Pr(text, profile):
    """Return the probability that Profile generates Text."""
    p = 1                                  # initialize probability
    for i in range (len(text)):            # iterate over each position in Text
       p = p * profile[text[i]][i]         # multiply probability for nucleotide at position i
    return p

def ProfileMostProbableKmer(text, k, profile):
    """Return the most probable k-mer in text given a profile matrix."""
    best_kmer = ""                        # the best k-mer (start with nothing)
    best_prob = -1                        # the highest probability (start with nothing)
    for i in range(len(text)-k+1):        # sliding window
        kmer = text[i:i+k]                # actual k-mer (first: "ACCTG")
        prob = Pr(kmer, profile)          # probability of this k-mer
        if prob > best_prob:              # probability is higher than actual probability?
            best_prob = prob              # if yes -> becomes the new highest probability
            best_kmer = kmer              # if yes -> becomes the new best k-mer
    return best_kmer

print(ProfileMostProbableKmer(text, k, profile))
