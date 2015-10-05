load("prescribed_roots.sage")
polRing.<x> = PolynomialRing(ZZ)
f = open("k3f5-lines.txt", "wb")
f.write("[5]\n")
for i in range(1, 11):
    temp, count = roots_on_unit_circle(5*(x^(2*i)+1), 
                                       filter=no_roots_of_unity,
                                       num_threads=768)
    print len(temp)
    for i in temp:
        f.write(str(list(i)))
    f.write("\n")
f.close()
