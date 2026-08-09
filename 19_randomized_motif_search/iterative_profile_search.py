import random

Dna = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
       "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
       "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
       "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
       "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]
k = 8
t = 5

def RandomizedMotifSearch(Dna, k, t):
    """Return the best motifs found by iteratively improving a random starting collection."""
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs


def RandomMotifs(Dna, k, t):
    """Return a list of t random k-mers, one from each string in Dna."""
    motifs = []
    n = len(Dna[0])                                    # length of each DNA string
    for i in range(t):                                 # iterate over each string in Dna
        start = random.randint(0, n - k)               # random start position for the k-mer
        motifs.append(Dna[i][start:start + k])         # extract k-mer at random position
    return motifs


def CountWithPseudocounts(Motifs):
    """Count nucleotide occurrences per column of Motifs with pseudocounts added."""
    count = {}
    t = len(Motifs)                       # number of motifs (rows)
    k = len(Motifs[0])                    # length of each motif (number of columns)
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)       # initialize with pseudocount instead of zero
    for i in range(t):                    # fill count matrix
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1         # increment count for nucleotide at column j
    return count


def ProfileWithPseudocounts(Motifs):
    """Return the profile matrix of Motifs as nucleotide frequencies, based on pseudocounts."""
    k = len(Motifs[0])                                           # length of each motif (number of columns)
    profile = {}
    count = CountWithPseudocounts(Motifs)                        # get nucleotide count matrix (with pseudocounts)
    t = sum(count[symbol][0] for symbol in "ACGT")               # column sum after pseudocounts
    for symbol in "ACGT":
        profile[symbol] = []                                     # initialize empty list for each nucleotide
        for j in range(k):
            profile[symbol].append(count[symbol][j] / t)         # divide by column total incl. pseudocounts
    return profile


def Motifs(Profile, Dna):
    """Return the profile-most probable k-mer from each string in Dna."""
    motifs = []
    k = len(Profile["A"])                              # k-mer length from profile width
    for i in range(len(Dna)):                           # iterate over each string in Dna
        motifs.append(ProfileMostProbableKmer(Dna[i], k, Profile))    # add best k-mer for this string
    return motifs


def ProfileMostProbableKmer(text, k, profile):
    """Return the most probable k-mer in text given a profile matrix."""
    best_kmer = ""                       # the best k-mer (start with nothing)
    best_prob = -1                       # the highest probability (start with nothing)
    for i in range(len(text) - k + 1):   # sliding window
        kmer = text[i:i + k]             # extract k-mer at position i
        prob = Pr(kmer, profile)         # probability of this k-mer
        if prob > best_prob:             # higher than current best?
            best_prob = prob             # if yes -> becomes the new highest probability
            best_kmer = kmer             # if yes -> becomes the new best k-mer
    return best_kmer


def Pr(text, profile):
    """Return the probability that profile generates text."""
    p = 1                                                            # initialize probability
    for i in range(len(text)):                                       # iterate over each position in text
        p = p * profile[text[i]][i]                                  # multiply probability for nucleotide at position i
    return p


def Score(Motifs):
    """Return the score of Motifs as the total number of mismatches with the consensus string."""
    count = 0
    consensus = Consensus(Motifs)               # get consensus string
    for i in range(len(Motifs)):                # iterate over all motifs (rows)
        for j in range(len(Motifs[0])):         # iterate over all columns
            if Motifs[i][j] != consensus[j]:    # mismatch with consensus
                count += 1                      # increment score
    return count


def Consensus(Motifs):
    """Return the consensus string by finding the most frequent nucleotide in each column."""
    k = len(Motifs[0])                          # length of each motif (number of columns)
    count = CountWithPseudocounts(Motifs)       # get nucleotide count matrix (with pseudocounts)
    consensus = ""                              # initialize empty consensus string
    for j in range(k):
        m = 0                                   # track maximum count
        frequentSymbol = ""                     # most frequent nucleotide in column j
        for symbol in "ACGT":
            if count[symbol][j] > m:            # new maximum found
                m = count[symbol][j]            # update if new maximum found
                frequentSymbol = symbol
        consensus += frequentSymbol             # append most frequent nucleotide to consensus
    return consensus

print(RandomizedMotifSearch(Dna, k, t))
