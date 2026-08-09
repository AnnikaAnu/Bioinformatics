import random

Dna = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
       "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
       "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
       "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
       "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]
k = 8
t = 5
N = 100

def GibbsSampler(Dna, k, t, N):
    """Return the best motifs found using iterative single-motif Gibbs sampling."""
    Motifs = RandomMotifs(Dna, k, t)
    BestMotifs = Motifs
    for j in range(N):
        i = random.randint(0, t - 1)                     # randomly pick a string to remove
        reduced_Motifs = Motifs[:i] + Motifs[i+1:]         # all motifs except the i-th
        Profile = ProfileWithPseudocounts(reduced_Motifs)
        Motifs[i] = ProfileGeneratedString(Dna[i], Profile, k)   # replace with weighted-random k-mer
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
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


def Pr(text, profile):
    """Return the probability that profile generates text."""
    p = 1                                                            # initialize probability
    for i in range(len(text)):                                       # iterate over each position in text
        p = p * profile[text[i]][i]                                  # multiply probability for nucleotide at position i
    return p


def Normalize(Probabilities):
    """Rescale a dictionary of probabilities so that all values sum to 1."""
    total = sum(Probabilities.values())                 # sum of all probabilities
    for key in Probabilities:                           # iterate over each k-mer
        Probabilities[key] = Probabilities[key] / total # rescale probability
    return Probabilities


def WeightedDie(Probabilities):
    """Return a k-mer chosen randomly, weighted by the values in Probabilities."""
    kmer = ''
    p = random.uniform(0, 1)                            # random number between 0 and 1
    running_total = 0                                   # cumulative probability so far
    for key in Probabilities:                           # iterate over each k-mer
        running_total += Probabilities[key]              # add this k-mer's probability
        if p < running_total:                            # p falls into this k-mer's range
            kmer = key
            break
    return kmer


def ProfileGeneratedString(Text, profile, k):
    """Return a randomly generated k-mer from Text, chosen with probabilities from profile."""
    n = len(Text)
    probabilities = {}
    for i in range(0, n - k + 1):                       # iterate over all k-mers in Text
        probabilities[Text[i:i + k]] = Pr(Text[i:i + k], profile)   # probability of each k-mer
    probabilities = Normalize(probabilities)             # rescale so probabilities sum to 1
    return WeightedDie(probabilities)                    # roll weighted die to pick a k-mer


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

print(GibbsSampler(Dna, k, t, N))
