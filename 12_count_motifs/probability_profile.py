Text = "ACGGGGATTACC"
Profile = {"A": [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
           "C": [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
           "G": [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
           "T": [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]}

def Pr(Text, Profile):
    """Return the probability that Profile generates Text."""
    p = 1                                   # initialize probability
    for i in range (len(Text)):            # iterate over each position in Text
       p = p * Profile[Text[i]][i]         # multiply probability for nucleotide at position i
    return p

print(Pr(Text, Profile))
