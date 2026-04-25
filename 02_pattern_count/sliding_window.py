# Count the occurrences of a k-mer (pattern) in a  DNA sequence (text)

def PatternCount(Text, Pattern):                    
    count = 0                                        # initialize counter
    for i in range(len(Text)-len(Pattern)+1):        # slide window across sequence
        if Text[i:i+len(Pattern)] == Pattern:        # k-mer match found
            count +=1                                # increment counter
    return count                                     # return total match count
  
print(PatternCount("GCGCG", "GCG"))                 # sample: GCG appears twice
