# Reverse a DNA string
def Reverse(Pattern):
    rev = ""                   # initialize empty string
    for char in Pattern:       # iterate over each nucleotide
        rev = char + rev       # add nucleotide to front of string
    return rev                 # return reversed DNA string
        
print(Reverse("AAAACCCGGT"))   # TGGCCCAAAA
