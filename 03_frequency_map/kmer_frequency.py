# Generate a frequency map of all k-mers in a DNA sequence

def FrequencyMap(Text, k):
    freq = {}                        # initialize dictionary
    n = len(Text)                    # calculate the length of the DNA sequence
    for i in range(n-k+1):           # slide window across sequence
        Pattern = Text[i:i+k]        # extract k-mer
        freq[Pattern] = 0            # set counter at 0
    for i in range(n-k+1):           # slide window across sequence
        Pattern = Text[i:i+k]        # extract k-mer
        freq[Pattern] += 1           # increment counter
    return freq                      # return frequency dictionary

print(FrequencyMap("CGATATATCCATAG", 3))     # most frequent k-mer: ATA (3 times)
