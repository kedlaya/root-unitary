polRing.<x> = PolynomialRing(ZZ)

t = True
CF = ComplexField(53)
for p in primes(10):
    print(p)
    l = []
    for i in range(1,5):
    	ans = list(WeilPolynomials(2*i, p, 1))
	l.append(ans)
        for j in ans:
            if not (j(-x) in ans):
                t = False
                print "Error: ", j
        for j in ans:
            for k in j.roots(CF):
                if (abs(abs(k[0])^2 - p) > 0.01):
                    t = False
                    print "Error: ", j
    for (j,k) in [(1,1), (1,2), (1,3), (2,2)]:
        for t1 in l[j-1]:
            for t2 in l[k-1]:
                if not (t1*t2 in l[j+k-1]):
                    t = False
                    print "Error: ", t1, t2
    #Todo: compute zeta functions of random hyperelliptic curves
    #Todo: compare against an exhaustive search in a box
                    
if t: print("All tests passed")
