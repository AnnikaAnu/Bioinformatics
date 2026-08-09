def Normalize(Probabilities):
    """Rescale a dictionary of probabilities so that all values sum to 1."""
    total = sum(Probabilities.values())                 # sum of all probabilities
    for key in Probabilities:                           # iterate over each k-mer
        Probabilities[key] = Probabilities[key] / total # rescale probability
    return Probabilities


Probabilities = {'A': 0.1, 'C': 0.1, 'G': 0.1, 'T': 0.1}

print(Normalize(Probabilities))
