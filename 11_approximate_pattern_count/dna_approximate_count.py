Pattern = "GAGG"
Text = "TTTAGAGCCTTCAGAGG"
d = 2

def ApproximatePatternCount(Pattern, Text, d):
    """Count occurrences of Pattern in Text with at most d mismatches."""
    count = 0
    for i in range(len(Text)-len(Pattern)+1):                          # slide window across sequence
        if HammingDistance(Pattern, Text[i:i+len(Pattern)]) <= d:      # pattern matches with at most d mismatches
            count += 1                             
    return count

def HammingDistance(p, q):
    """Return the Hamming distance between two equal-length strings."""
    count = 0
    for i in range (len(p)):
        if p[i] != q [i]:       #mismatch found
            count += 1
    return count

print(ApproximatePatternCount(Pattern, Text, d))
