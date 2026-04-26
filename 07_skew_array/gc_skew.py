Genome = "GAGCCACCGCGATA"

def skew_array(Genome):
    """Compute the skew array as the running difference of G and C counts across the genome."""
    skew = {0: 0}                     # initialize skew array with Skew[0] = 0
    n = len(Genome)
    for i in range(n):
        # update skew: +1 for G, -1 for C, 0 for A or T
        if Genome [i] == "A":
            skew[i+1] = skew[i]
        elif Genome [i] == "T":
            skew[i+1] = skew[i]
        elif Genome [i] == "G":
            skew[i+1] = skew[i] +1
        elif Genome [i] == "C":
            skew[i+1] = skew[i] -1
    return skew

print(skew_array(Genome))
