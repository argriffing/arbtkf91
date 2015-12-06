import math
import random
from fractions import gcd

def decompose_pairwise_coprime(A):
    B = set()
    C = set(A)
    count = 0
    while C:
        count += 1
        #print B
        #print C
        c = C.pop()
        #print 'removed', c, 'from C'
        found_divisor = False
        for b in B:
            d = gcd(b, c)
            if d != 1:
                found_divisor = True
                break
        if found_divisor:
            B.remove(b)
            #print 'removed', b, 'from B'
            for value in c // d, b // d, d:
                if value != 1:
                    C.add(value)
                    #print 'added', value, 'to C'
        else:
            B.add(c)
            #print 'added', c, 'to B'
        if 1 in B.union(C):
            #print c, b, d
            break
    #print B
    #print C
    return B, count


def main():
    #A = [random.randrange(2, 1000) for i in range(5)]
    #A.append(A[0] * A[0])
    #A.append(A[1] * A[3])
    A = [2*3, 2*5]
    B, count = decompose_pairwise_coprime(A)
    A = sorted(A)
    B = sorted(B)
    print sorted(A)
    print sorted(B)
    print count, 'iterations'

    # Represent each element of A in terms of powers of elements of B.
    M = []
    for a in A:
        v = [0] * len(B)
        for i, b in enumerate(B):
            count = 0
            while a != 1 and a % b == 0:
                a //= b
                v[i] += 1
        M.append(v)

    for v in M:
        print v

    print "{",
    for v in M:
        print "{" + ','.join(str(x) for x in v) + "},",
    print "}",

main()
