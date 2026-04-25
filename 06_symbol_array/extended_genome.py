Genome = "AAAAGGGG"
symbol = "A"

def SymbolArray(Genome, symbol):
    """Count occurrences of symbol in each half-genome window."""
    array = {}                                      
    n = len(Genome)                                 
    ExtendedGenome = Genome + Genome[0:n//2]                          # extend genome to handle circular DNA
    for i in range(n):                                                # iterate over each position in genome
        array[i] = PatternCount(ExtendedGenome[i:i+(n//2)], symbol)   # count symbol occurrences in window starting at position i
    return array

def PatternCount(Text, Pattern):
    """Count the occurrences of a k-mer (pattern) in a  DNA sequence (text)."""
    count = 0                                       
    for i in range(len(Text)-len(Pattern)+1):                        # slide window across sequence
        if Text[i:i+len(Pattern)] == Pattern:                        # k-mer match found
            count +=1                                                # increment counter
    return count                                     

print(SymbolArray(Genome, symbol))
