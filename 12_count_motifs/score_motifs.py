def Consensus(Motifs):
    """Return the consensus string by finding the most frequent nucleotide in each column."""
    k = len(Motifs[0])                        # length of each motif (number of columns)
    count = Count(Motifs)                     # get nucleotide count matrix
    consensus = ""                            # initialize empty consensus string
    for j in range(k):
        m = 0                                 # track maximum count
        frequentSymbol = ""                   # most frequent nucleotide in column j
        for symbol in "ACGT":
            if count[symbol][j] > m:          # new maximum found
                m = count[symbol][j]          # update if new maximum found
                frequentSymbol = symbol
        consensus += frequentSymbol           # append most frequent nucleotide to consensus
    return consensus


def Count(Motifs):
    """Count nucleotide occurrences in each column of Motifs and return a count matrix."""
    count = {} 

    k = len(Motifs[0])                        # length of each motif (number of columns)
    for symbol in "ACGT":                     # initialize count matrix with zeros
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)

    t = len(Motifs)                            # number of motifs (rows)
    for i in range(t):                         # fill count matrix
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] +=1               # increment count for nucleotide at column j

    return count


Motifs = ["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"]


def Score(Motifs):
    """Return the score of Motifs as the total number of mismatches with the consensus string."""
    count = 0
    consensus = Consensus(Motifs)               # call the most frequent nucleotide
    for i in range(len(Motifs)):                # iterate over all motifs (rows)  
        for j in range(len(Motifs[0])):         # iterate over all columns
            if Motifs[i][j] != consensus[j]:    # mismatch with consensus
                count += 1                      # increment score
    return count
            

print(Score(Motifs))
