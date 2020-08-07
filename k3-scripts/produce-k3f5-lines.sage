from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials

polRing.<x> = PolynomialRing(ZZ)
f = open("k3f5-lines.txt", "w")
f.write("[5]\n")
c = 0
for i in range(1, 11):
    wp = WeilPolynomials(2*i, 1, sign=1, lead=5, squarefree=True, parallel=False)
    for i in wp:
        if not i.has_cyclotomic_factor():
            f.write(str(list(i)) + "\n")
            c += 1
    print(c)
f.close()
