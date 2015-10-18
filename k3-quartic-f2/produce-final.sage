polRing.<x> = PolynomialRing(ZZ)
powRing.<y> = PowerSeriesRing(QQ, 30)
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
for h in range(13, 17):
    s = "k3miss" + str(h) + "pts.txt"
    f = open(s, "rb")
    for i in f:
        j = i.split()
        k = tuple([int(j[t]) for t in range(1,h+1)])
        k1 = tuple([int(j[t]) for t in range(1,13)])
        l1[k1].remove(j[0])
        if len(l1[k1]) == 0:
            del l1[k1]
        if k in l1:
            l1[k].append(j[0])
        else:
            l1[k] = [j[0]]
    f.close()
print "Loaded", len(l1.keys()), "distinct point count lists"
l2 = []
l3 = []
for i in l1.keys():
    k1 = powRing([0] + [QQ(i[j-1])/j for j in range(1,len(i)+1)])
    k2 = (exp(k1)*(1-y)*(1-2*y)*(1-4*y)).inverse()
    m = [2^(1-j)*k2[j] for j in range(len(i)+1)]
    j = 11
    while j < len(m) and m[j]==0: 
        j += 1
    if j >= len(m):
        l3.append([l1[i], i])
        continue
    s = sign(m[21-j]) // sign(m[j])
    m1 = m[:]
    for k in range(len(i)+1, 22):
        m1.append(s*m[21-k])
    l2.append([polRing(m1), [int(j) for j in l1[i]]])
print len(l2), " zeta functions uniquely reconstructed"
f = open("k3final.txt", "wb")
for i in l2:
    f.write(str(i))
    f.write("\n")
f.close()
print sum(len(i[0]) for i in l3), " ambiguous cases remain"
