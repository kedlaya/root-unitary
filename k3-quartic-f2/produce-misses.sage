polRing.<x> = PolynomialRing(ZZ)
powRing.<y> = PowerSeriesRing(QQ)
l1 = {}
f = open("k3points.txt", "rb")
for i in f:
    j = i.split()
    k = tuple([int(j[t]) for t in range(1,13)])
    if k in l1:
        l1[k].append(j[0])
    else:
        l1[k] = [j[0]]
f.close()
print "Loaded", len(l1.keys()), "distinct point count lists"
l2 = []
for i in l1.keys():
    k1 = powRing([0] + [QQ(i[j-1])/j for j in range(1,13)])
    k2 = (exp(k1)*(1-y)*(1-2*y)*(1-4*y)).inverse()
    m = [2^(1-j)*k2[j] for j in range(13)]
    if m[11] == 0 and m[12] == 0:
        k = 13
        while m[21-k] == 0: k += 1
        for j in l1[i]:
            l2.append([j, k])
print len(l2), " ambiguous cases"
f = open("k3misses.txt", "wb")
for i in l2:
    f.write(i[0])
    f.write(" ")
    f.write(str(i[1]))
    f.write("\n")
f.close()
