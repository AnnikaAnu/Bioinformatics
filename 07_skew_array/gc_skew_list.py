Genome = "CATGGGCATCGGCCATACGCC"

def SkewArray(Genome):
    """Compute the skew array as a list of running G-C differences across the genome."""
    skew = [0]                                             # initialize list with Skew[0] = 0                   
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

print(SkewArray(Genome))
