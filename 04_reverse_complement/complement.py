# Complement a DNA string
def Complement (pattern):
    comp = ""                    # initialize empty string
    for char in pattern:         # iterate over each nucleotide
        # replace each nucleotide with its Watson-Crick complement
        if char == "A":        
            comp += "T"
        elif char == "T":
            comp += "A"
        elif char == "C":
            comp += "G"
        elif char == "G":
            comp += "C"
    return comp                   # return complementary DNA string
        
print(Complement("AAAACCCGGT"))   # TTTTGGGCCA

