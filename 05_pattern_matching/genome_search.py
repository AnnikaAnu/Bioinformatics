# Find all occurrences of a pattern in a DNA sequence (genome)
def PatternMatching(Pattern, Genome):                    
    positions = []                                    # initialize empty list
    for i in range(len(Genome)-len(Pattern)+1):       # slide window across sequence
        if Genome[i:i+len(Pattern)] == Pattern:       # k-mer match found
            positions.append(i)                       # append position to list
    return positions                                  # return all matching positions
  
print(PatternMatching("ATAT","GATATATGCATATACTT"))    # positions: 1, 3, 9
