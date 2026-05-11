Dna = ["GGCGTTCAGGCA",
       "AAGAATCAGTCA",
       "CAAGGAGTTCGC",
       "CACGTCAATCAC",
       "CAATAATATTCG"]
k = 3
t = 5


def Score(Motifs):
    """Return the score of Motifs as the total number of mismatches with the consensus string."""
    count = 0
    consensus = Consensus(Motifs)               # call the most frequent nucleotide
    for i in range(len(Motifs)):                # iterate over all motifs (rows)  
        for j in range(len(Motifs[0])):         # iterate over all columns
            if Motifs[i][j] != consensus[j]:    # mismatch with consensus
                count += 1                      # increment score
    return count


def Count(Motifs):
    """Count nucleotide occurrences in each column of Motifs and return a count matrix."""
    count = {} 

    k = len(Motifs[0])                         # length of each motif (number of columns)
    for symbol in "ACGT":                      # initialize count matrix with zeros
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)

    t = len(Motifs)                           # number of motifs (rows)
    for i in range(t):                        # fill count matrix
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] +=1              # increment count for nucleotide at column j

    return count


def Profile(Motifs):
    """Return the profile matrix of Motifs as nucleotide frequencies per column."""
    t = len(Motifs)                          # number of motifs (rows)
    k = len(Motifs[0])                       # length of each motif (number of columns)
    profile = {}
    count = Count(Motifs)                    # get nucleotide count matrix
    for symbol in "ACGT":
        profile[symbol] = []                 # initialize empty list for each nucleotide
        for j in range(k):
            profile[symbol].append(count[symbol][j] / t)     # divide count by total rows to get frequency
    return profile


def Consensus(Motifs):
    """Return the consensus string by finding the most frequent nucleotide in each column."""
    k = len(Motifs[0])                      # length of each motif (number of columns)
    count = Count(Motifs)                   # get nucleotide count matrix
    consensus = ""                          # initialize empty consensus string
    for j in range(k):
        m = 0                               # track maximum count
        frequentSymbol = ""                 # most frequent nucleotide in column j
        for symbol in "ACGT":
            if count[symbol][j] > m:        # new maximum found
                m = count[symbol][j]        # update if new maximum found
                frequentSymbol = symbol
        consensus += frequentSymbol         # append most frequent nucleotide to consensus
    return consensus


def ProfileMostProbableKmer(text, k, profile):
    """Return the most probable k-mer in text given a profile matrix."""
    best_kmer = ""                                                    # the best k-mer (start with nothing)
    best_prob = -1                                                    # the highest probability (start with nothing)
    for i in range(len(text)-k+1):                                    # sliding window
        kmer = text[i:i+k]                                            # actual k-mer (first: "ACCTG")
        prob = Pr(kmer, profile)                                      # probability of this k-mer
        if prob > best_prob:                                          # probability is higher than actual probability?
            best_prob = prob                                          # if yes -> becomes the new highest probability
            best_kmer = kmer                                          # if yes -> becomes the new best k-mer
    return best_kmer


def Pr(text, profile):
    """Return the probability that Profile generates Text."""
    p = 1                                                             # initialize probability
    for i in range (len(text)):                                       # iterate over each position in Text
       p = p * profile[text[i]][i]                                    # multiply probability for nucleotide at position i
    return p


def GreedyMotifSearch(Dna, k, t):
    """Return the best scoring motifs from Dna using a greedy approach."""
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])                                # initialize BestMotifs with first k-mer of each string
    n = len(Dna[0])                                                   # length of each DNA string
    for i in range(n-k+1):                                            # try every k-mer in Dna[0] as starting motif
        Motifs = []
        Motifs.append(Dna[0][i:i+k])                                  # start with i-th k-mer from Dna[0]
        for j in range(1, t):                                         # extend motifs to remaining strings
            P = Profile(Motifs[0:j])                                  # build profile from current motifs
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))      # find most probable k-mer in Dna[j]
        if Score(Motifs) < Score(BestMotifs):                         # update if better motifs found
            BestMotifs = Motifs
    return BestMotifs


print(GreedyMotifSearch(Dna, k, t))

