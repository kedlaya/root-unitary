# Depends on weil_polynomials.pyx, trac 23948

import time
polRing.<x> = PolynomialRing(ZZ)
ans = [polRing(2)]
t = time.time()
c = 0
for i in range(1, 12):
    wp = WeilPolynomials(2*i, 1, sign=1, lead=2, num_threads=512)
    l = [j for j in wp if not j.has_cyclotomic_factor()]
    ans += l
    c += wp.count
    print len(l), "polynomials added"
    print c, "nodes enumerated"
    print "time so far: ", time.time() - t, " seconds"

f = open("cubic4foldf2-lines.txt", "wb")
for i in ans:
    f.write(str(i.list()) + "\n")
f.close()
