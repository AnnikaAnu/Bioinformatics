Genome = "AAAAGGGG"
symbol = "A"

def FasterSymbolArray(Genome, symbol):
    """Build symbol array efficiently by tracking changes between adjacent windows."""
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]             # extend genome to handle circular DNA

    array[0] = PatternCount(symbol, Genome[0:n//2])      # count symbol in first window

    for i in range(1, n):                                # start at 1 since array[0] is already computed
        array[i] = array[i-1]                            # start with previous window count
        if ExtendedGenome[i-1] == symbol:                # if symbol leaves window on the left
            array[i] = array[i]-1                        # decrease count
        if ExtendedGenome[i+(n//2)-1] == symbol:         # if symbol enters window on the right
            array[i] = array[i]+1                        # increase count
    return array

def PatternCount(Pattern, Text):
    """Count the occurrences of a k-mer (pattern) in a  DNA sequence (text)."""
    count = 0                                        
    for i in range(len(Text)-len(Pattern)+1):        # slide window across sequence
        if Text[i:i+len(Pattern)] == Pattern:        # k-mer match found
            count +=1                                
    return count                                    

import sys
lines = sys.stdin.read().splitlines()
print(FasterSymbolArray(lines[0],lines[1]))



