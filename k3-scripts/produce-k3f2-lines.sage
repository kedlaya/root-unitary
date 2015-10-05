load("prescribed_roots.sage")
polRing.<x> = PolynomialRing(ZZ)
ans = [polRing(2)]
for i in range(1, 11):
    temp, count = roots_on_unit_circle(2*(x^(2*i)+1), 
                                       filter=no_roots_of_unity,
                                       num_threads=512)
    print len(temp)
    ans += temp
f = open("k3-scripts/k3f2-lines.txt", "wb")
for i in ans:
    f.write(str(i.list()))
    f.write("\n")
f.close()
