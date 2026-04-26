Pattern = "ATTCTGGA"
Text = "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT"
d = 3

def ApproximatePatternMatching(Text, Pattern, d):
    """Return all positions where Pattern appears in Text with at most d mismatches."""
    positions = [] 
    for i in range(len(Text)-len(Pattern)+1):                          # slide window across sequence
        if HammingDistance(Pattern, Text[i:i+len(Pattern)]) <= d:      # pattern matches with at most d mismatches
            positions.append(i)                             
    return positions

def HammingDistance(p, q):
    """Return the Hamming distance between two equal-length strings."""
    count = 0
    for i in range (len(p)):
        if p[i] != q [i]:       #mismatch found
            count += 1
    return count

print(ApproximatePatternMatching(Text, Pattern, d))
