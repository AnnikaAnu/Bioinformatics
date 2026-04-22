# Compute the reverse complement of a DNA string
def ReverseComplement(Pattern):         
    Pattern = Reverse(Pattern)            # reverse DNA string
    Pattern = Complement(Pattern)         # complement DNA string
    return Pattern                        # return reversed and complemented DNA string
  
# Reverse a DNA string
def Reverse(Pattern):
    rev = ""                              # initialize empty string
    for char in Pattern:                  # iterate over each nucleotide
        rev = char + rev                  # add nucleotide to front of string
    return rev                            # return reversed DNA string

# Complement a DNA string
def Complement (Pattern):
    comp = ""                             # initialize empty string
    for char in Pattern:                  # iterate over each nucleotide
        # replace each nucleotide with its Watson-Crick complement
        if char == "A":        
            comp += "T"
        elif char == "T":
            comp += "A"
        elif char == "C":
            comp += "G"
        elif char == "G":
            comp += "C"
    return comp                          # return complementary DNA string
        
print(ReverseComplement("AAAACCCGGT"))   # ACCGGGTTTT
