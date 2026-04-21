# Find the most frequent k-mers in a DNA sequence
def FrequentWords(Text, k):
    words = []                                          # initialize empty list
    freq = FrequencyMap(Text, k)                        # all k-mers with frequency
    m = max(freq.values())                              # find maximum frequency
    for key in freq:                                    # loop over all k-mers
        if freq[key] == m and key not in words:         # the most frequent k-mer? # no duplicates
            words.append(key)                           # yes -> append to list
    return words                                        # return all most frequent k-mers

# Generate a frequency map of all k-mers in a DNA sequence
def FrequencyMap(Text, k):
    freq = {}                        # initialize dictionary
    n = len(Text)                    # calculate the length of the DNA sequence
    for i in range(n-k+1):           # slide window across sequence
        Pattern = Text[i:i+k]        # extract k-mer
        if Pattern not in freq:      # first occurrence
            freq[Pattern] = 1        # initialize counter
        else:
            freq[Pattern] += 1       #  increment counter
                   
    return freq                      # return frequency dictionary

print(FrequentWords("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4))    # most frequent 4-mers: CATG, GCAT

