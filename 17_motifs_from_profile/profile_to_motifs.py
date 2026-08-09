Profile = {
    'A': [0.8, 0.0, 0.0, 0.2],
    'C': [0.0, 0.6, 0.2, 0.0],
    'G': [0.2, 0.2, 0.8, 0.0],
    'T': [0.0, 0.2, 0.0, 0.8]
}

Dna = ["TTACCTTAAC",
       "GATGTCTGTC",
       "ACGGCGTTAG",
       "CCCTAACGAG",
       "CGTCAGAGGT"]


def Motifs(Profile, Dna):
    """Return the profile-most probable k-mer from each string in Dna."""
    motifs = []
    k = len(Profile["A"])                                             # k-mer length from profile width
    for i in range(len(Dna)):                                         # iterate over each string in Dna
        motifs.append(ProfileMostProbableKmer(Dna[i], k, Profile))    # add best k-mer for this string
    return motifs


def ProfileMostProbableKmer(text, k, profile):
    """Return the most probable k-mer in text given a profile matrix."""
    best_kmer = ""                       # the best k-mer (start with nothing)
    best_prob = -1                       # the highest probability (start with nothing)
    for i in range(len(text)-k+1):       # sliding window
        kmer = text[i:i+k]               # extract k-mer at position i
        prob = Pr(kmer, profile)         # probability of this k-mer
        if prob > best_prob:             # probability is higher than actual probability?
            best_prob = prob             # if yes -> becomes the new highest probability
            best_kmer = kmer             # if yes -> becomes the new best k-mer
    return best_kmer


def Pr(text, profile):
    """Return the probability that profile generates text."""
    p = 1                                                            # initialize probability
    for i in range(len(text)):                                       # iterate over each position in text
        p = p * profile[text[i]][i]                                  # multiply probability for nucleotide at position i
    return p


print(Motifs(Profile, Dna))
