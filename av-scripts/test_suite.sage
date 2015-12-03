polRing.<x> = PolynomialRing(ZZ)

t = True
CF = ComplexField(53)
for p in primes(25):
    print(p)
    l = []
    for i in range(1,5):
        ans, count = roots_on_unit_circle(1 + (p*x^2)^i)
        l.append(ans)
        for i in ans:
            if not (i(-x) in ans):
                t = False
                print "Error: ", i
        for i in ans:
            for j in i.roots(CF):
                if (abs(p*abs(j[0])^2 - 1) > 0.01):
                    t = False
                    print "Error: ", i
    for (i,j) in [(1,1), (1,2), (1,3), (2,2)]:
        for t1 in l[i-1]:
            for t2 in l[j-1]:
                if not (t1*t2 in l[i+j-1]):
                    t = False
                    print "Error: ", t1, t2
    #Todo: compute zeta functions of random hyperelliptic curves
    #Todo: compare against an exhaustive search in a box
                    
if t: print("All tests passed")
