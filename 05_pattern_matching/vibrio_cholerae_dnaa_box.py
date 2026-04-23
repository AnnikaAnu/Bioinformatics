# Locate all DnaA box (CTTGATCAT) positions in the Vibrio cholerae genome
def PatternMatching(Pattern, Genome):                    
    Positions = []                                    # initialize empty list
    for i in range(len(Genome)-len(Pattern)+1):       # slide window across sequence
        if Genome[i:i+len(Pattern)] == Pattern:       # k-mer match found
            Positions.append(i)                       # append position to list
    return Positions                                  # return all matching positions

# Vibrio cholerae genome source:
# https://bioinformaticsalgorithms.com/data/realdatasets/Replication/Vibrio_cholerae.txt
Pattern = "CTTGATCAT"                                 # DnaA box pattern
Positions = PatternMatching(Pattern, Genome)          # search in Vibrio cholerae genome
print(Positions)                                      # print all matching positions, 16 DnaA box positions found       
