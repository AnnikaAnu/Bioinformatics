Motifs = ["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"]

def ProfileWithPseudocounts(Motifs):
    """Return the profile matrix of Motifs as nucleotide frequencies, based on pseudocounts."""
    k = len(Motifs[0])
    profile = {}
    count = CountWithPseudocounts(Motifs)                        # get nucleotide count matrix (with pseudocounts)
    t = sum(count[symbol][0] for symbol in "ACGT")               # column sum after pseudocounts
    
    for symbol in "ACGT":
        profile[symbol] = []                                     # initialize empty list for each nucleotide
        for j in range(k):
            profile[symbol].append(count[symbol][j] / t)         # divide by column total incl. pseudocounts
    return profile

def CountWithPseudocounts(Motifs):
    """Count nucleotide occurrences per column of Motifs with pseudocounts added."""
    
    t = len(Motifs)                       # number of motifs (rows)
    k = len(Motifs[0])                    # length of each motif (number of columns)
    count = {}
    
    for symbol in "ACGT":                 # initialize with pseudocount instead of zero
        count[symbol] = []
        for j in range(k):                
            count[symbol].append(1)      

    for i in range(t):                    # fill count matrix
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1         # increment count for nucleotide at column j

    return count

print(ProfileWithPseudocounts(Motifs))
