load("../prescribed_roots.sage")
polRing.<x> = PolynomialRing(ZZ)
l = []
l.append([x+1, x-1])
for i in range(1, 11):
    l.append([])
    ans, count = roots_on_unit_circle(x^(2*i)+1)
    for j in ans:
        l[-1].append(j * (x+1))
        l[-1].append(j * (x-1))
f = open("k3f1-list.txt", "wb")
f.write(str(l))
f.close()
