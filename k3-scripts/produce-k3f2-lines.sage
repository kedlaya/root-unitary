load("prescribed_roots.sage")
polRing.<x> = PolynomialRing(ZZ)
ans = [polRing(2)]
for i in range(1, 11):
    temp, count = roots_on_unit_circle(2*(x^(2*i)+1), filter=no_roots_of_unity)
    print len(temp)
    ans += temp
f = open("k3-scripts/k3f2-list.txt", "wb")
for i in range(11):
    for j in l[i]:
        f.write(str(j.list()))
        f.write("\n")
f.close()
