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


Motifs = ["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"]


def Profile(Motifs):
    """Return the profile matrix of Motifs as nucleotide frequencies per column."""
    t = len(Motifs)                      # number of motifs (rows)
    k = len(Motifs[0])                   # length of each motif (number of columns)
    profile = {}
    count = Count(Motifs)                # get nucleotide count matrix
    for symbol in "ACGT":
        profile[symbol] = []             # initialize empty list for each nucleotide
        for j in range(k):
            profile[symbol].append(count[symbol][j] / t)     # divide count by total rows to get frequency
    return profile


print(Profile(Motifs))

