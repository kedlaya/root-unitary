# Depends on weil_polynomials.pyx, trac 23948

polRing.<x> = PolynomialRing(ZZ)
f = open("k3f5-lines.txt", "wb")
f.write("[5]\n")
for i in range(1, 11):
    wp = WeilPolynomials(2*i, 1, sign=1, lead=5, num_threads=512)
    for i in wp:
        if not i.has_cyclotomic_factor():
	    f.write(str(list(i)) + "\n")
            c += 1
    print c
f.close()
