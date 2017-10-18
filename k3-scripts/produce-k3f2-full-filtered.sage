# Depends on ej-test.sage

polRing.<x> = PolynomialRing(ZZ)
powRing.<y> = PowerSeriesRing(QQ, 5)
l1 = []
for i in range(11):
    l1.append([])
f = open("k3f1-lines.txt", "rb")
for i in f:
    j = eval(i)
    k = (len(j)-2)/2
    l1[k].append(polRing(j))
f.close()
print "Loaded", sum(len(i) for i in l1), "1-polynomials"
l2 = []
for i in range(11):
    l2.append([])
f = open("k3f2-lines.txt", "rb")
for i in f:
    j = eval(i)
    k = (len(j)-1)/2
    l2[k].append(polRing(j))
f.close()
print "Loaded", sum(len(i) for i in l2), "2-polynomials"
l3 = []
for i in range(11):
    for j in l1[i]:
        for k in l2[10-i]:
	    l3.append(j*k)
print "Full set of products:", len(l3), "polynomials"
l4 = [i for i in l3 if ej_test(i)]
print "Satisfying Artin-Tate condition:", len(l4), "polynomials"
l5 = [i for i in l4 if -i[1]+7 >= 0 and i[1]^2 - 4*i[2] + 21 >= -i[1] + 7 and
      -i[1]^3 + 6*i[1]*i[2] - 12*i[3] + 73 >= -i[1] + 7 and
      i[1]^4 - 8*i[1]^2*i[2] + 8*i[2]^2 + 16*i[1]*i[3] - 32*i[4] + 273 >=
      i[1]^2 - 4*i[2] + 21]
print "Satisfying nonnegativity:", len(l5), "polynomials"
f = open("k3f2-full-filtered.txt", "wb")
for i in l5:
    f.write(str(i.list()) + "\n")
f.close()



