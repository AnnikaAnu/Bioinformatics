import random

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


Probabilities = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}

print(WeightedDie(Probabilities))
