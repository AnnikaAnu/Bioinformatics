# Count the occurrences of a k-mer (pattern) in a  DNA sequence (text)

def pattern_count(text, pattern):                    
    count = 0                                        # initialize counter
    for i in range(len(text)-len(pattern)+1):        # slide window across sequence
        if text[i:i+len(pattern)] == pattern:        # k-mer match found
            count +=1                                # increment counter
    return count                                     # return total match count

# Vibrio cholerae oriC — searching for potential DnaA box (TGATCA)
text = "ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC"
pattern = "TGATCA"

print(pattern_count(text, pattern))                  # TGATCA appears 8 times 
