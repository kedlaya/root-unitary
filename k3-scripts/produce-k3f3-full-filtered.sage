# Depends on ej-test.sage

polRing.<x> = PolynomialRing(ZZ)
powRing.<y> = PowerSeriesRing(QQ, 8)
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
f = open("k3f3-lines.txt", "rb")
for i in f:
    j = eval(i)
    k = (len(j)-1)/2
    l2[k].append(polRing(j))
f.close()
print "Loaded", sum(len(i) for i in l2), "3-polynomials"
l3 = []
for i in range(11):
    for j in l1[i]:
        for k in l2[10-i]:
	    l3.append(j*k)
print "Full set of products:", len(l3), "polynomials"
l4 = [i for i in l3 if ej_test(i)]
print "Satisfying Artin-Tate condition:", len(l4), "polynomials"
del l3
l5 = [i for i in l4 if -i[1] +13 >= 0 and i[1]^2-6*i[2]+91 >= -i[2]+13]
print "Satisfying nonnegativity:", len(l5), "polynomials"
f = open("k3f3-full-filtered.txt", "wb")
for i in l5:
    f.write(str(list(i)) + "\n")
f.close()


