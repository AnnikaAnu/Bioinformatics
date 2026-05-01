Motifs = ["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"]

def Count(Motifs):
    """Count nucleotide occurrences in each column of Motifs and return a count matrix."""
    count = {} 

    k = len(Motifs[0])                   # length of each motif (number of columns)
    for symbol in "ACGT":                # initialize count matrix with zeros
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)

    t = len(Motifs)                       # number of motifs (rows)
    for i in range(t):                    # fill count matrix
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] +=1          # increment count for nucleotide at column j

    return count

print(Count(Motifs))
