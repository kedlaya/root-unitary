from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials

import time
polRing.<x> = PolynomialRing(ZZ)
ans = [polRing(2)]
t = time.time()
c = 0
for i in range(1, 11):
    wp = WeilPolynomials(2*i, 1, sign=1, lead=2, squarefree=True)
    wp.num_threads = 512
    l = [j for j in wp if not j.has_cyclotomic_factor()]
    ans += l
    c += wp.node_count()
    print(len(l), "polynomials added")
    print(c, "nodes enumerated")
    print("time so far: ", time.time() - t, " seconds")

with open("k3f2-lines.txt", "w") as f:
    for i in ans:
        f.write(str(i.list()) + "\n")
