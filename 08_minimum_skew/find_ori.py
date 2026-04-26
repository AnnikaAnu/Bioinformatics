Genome = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"

def MinimumSkew(Genome):
    """Return all positions where the skew array reaches its minimum value."""
    positions = []             
    skew = SkewArray(Genome)
    min_skew =  min(skew)          # find minimum skew value
    for i in range(len(skew)):     
        if skew[i] == min_skew:    # position achieves minimum skew
            positions.append(i)
    return positions

def SkewArray(Genome):
    """Compute the skew array as a list of running G-C differences across the genome."""
    skew = [0]                     # initialize list with Skew[0] = 0                   
    n = len(Genome)
    for i in range(n):
        # update skew: +1 for G, -1 for C, 0 for A or T
        if Genome [i] == "A":
            skew.append(skew[-1])
        elif Genome [i] == "T":
             skew.append(skew[-1])
        elif Genome [i] == "G":
             skew.append(skew[-1] +1)
        elif Genome [i] == "C":
             skew.append(skew[-1] -1)
    return skew

print (MinimumSkew(Genome))
