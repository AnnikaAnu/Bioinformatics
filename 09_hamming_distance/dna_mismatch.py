p = "GGGCCGTTGGT"
q = "GGACCGTTGAC"

def HammingDistance(p, q):
    """Return the Hamming distance between two equal-length strings."""
    count = 0
    for i in range (len(p)):
        if p[i] != q [i]:       #mismatch found
            count += 1
    return count

print(HammingDistance(p, q))
