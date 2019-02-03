# Depends on weil_polynomials.pyx, trac 23948

import time
polRing.<x> = PolynomialRing(ZZ)
ans = [polRing(3)]
t = time.time()
c = 0
for i in range(1, 11):
    wp = WeilPolynomials(2*i, 1, sign=1, lead=3, parallel=True, squarefree=True)
    l = [j for j in wp if not j.has_cyclotomic_factor()]
    ans += l
    c += wp.node_count()
    print len(l), "polynomials added"
    print c, "nodes enumerated"
    print "time so far: ", time.time() - t, " seconds"

f = open("k3f3-lines.txt", "wb")
for i in ans:
    f.write(str(i.list()) + "\n")
f.close()
