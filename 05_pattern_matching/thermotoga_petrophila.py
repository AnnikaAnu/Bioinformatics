# Count occurrences of a k-mer (pattern) in a DNA sequence (text)
def PatternCount(Text, Pattern):                    
    count = 0                                        # initialize counter
    for i in range(len(Text)-len(Pattern)+1):        # slide window across sequence
        if Text[i:i+len(Pattern)] == Pattern:        # k-mer match found
            count +=1                                # increment counter
    return count                                     # return total match count

# Thermotoga petrophila oriC region
Text = "AACTCTATACCTCCTTTTTGTCGAATTTGTGTGATTTATAGAGAAAATCTTATTAACTGAAACTAAAATGGTAGGTTTGGTGGTAGGTTTTGTGTACATTTTGTAGTATCTGATTTTTAATTACATACCGTATATTGTATTAAATTGACGAACAATTGCATGGAATTGAATATATGCAAAACAAACCTACCACCAAACTCTGTATTGACCATTTTAGGACAACTTCAGGGTGGTAGGTTTCTGAAGCTCTCATCAATAGACTATTTTAGTCTTTACAAACAATATTACCGTTCAGATTCAAGATTCTACAACGCTGTTTTAATGGGCGTTGCAGAAAACTTACCACCTAAAATCCAGTATCCAAGCCGATTTCAGAGAAACCTACCACTTACCTACCACTTACCTACCACCCGGGTGGTAAGTTGCAGACATTATTAAAAACCTCATCAGAAGCTTGTTCAAAAATTTCAATACTCGAAACCTACCACCTGCGTCCCCTATTATTTACTACTACTAATAATAGCAGTATAATTGATCTGA"

count_1 = PatternCount(Text, "ATGATCAAG")            # count DnaA box forward strand 
count_2 = PatternCount(Text, "CTTGATCAT")            # count DnaA box reverse strand
total = count_1 + count_2                            # total DnaA box occurrences
print(total)                                         # combined occurrences in Thermotoga petrophila oriC (total: 0)

