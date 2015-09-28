load("../prescribed_roots.sage")
polRing.<x> = PolynomialRing(ZZ)
ans = []
ans.append([polRing(2)])
for i in range(1, 11):
    temp, count = roots_on_unit_circle(2*(x^(2*i)+1), filter=no_roots_of_unity)
    print len(temp)
    ans.append(temp)
f = open("k3f2-list.txt", "wb")
f.write(str(ans))
f.close()
