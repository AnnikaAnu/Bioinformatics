def ProfileGeneratedString(Text, profile, k):
    """Return a randomly generated k-mer from Text, chosen with probabilities from profile."""
    n = len(Text)
    probabilities = {}
    for i in range(0, n - k + 1):                       # iterate over all k-mers in Text
        probabilities[Text[i:i + k]] = Pr(Text[i:i + k], profile)   # probability of each k-mer
    probabilities = Normalize(probabilities)             # rescale so probabilities sum to 1
    return WeightedDie(probabilities)                    # roll weighted die to pick a k-mer


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


Text = "AAACCCAAACCC"
profile = {'A': [0.5, 0.1], 'C': [0.3, 0.2], 'G': [0.2, 0.4], 'T': [0.0, 0.3]}
k = 2

print(ProfileGeneratedString(Text, profile, k))
