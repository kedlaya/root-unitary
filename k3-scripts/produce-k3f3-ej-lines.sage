polRing.<x> = PolynomialRing(ZZ)
f = open("k3f1-list.txt", "rb")
l = sage_eval(f.read(), locals={'x':x})
f.close()
l2 = []
f = open("k3f3-lines.txt", "rb")
for i in f:
    l2.append(eval(i))
f.close()
l3 = []
for i in range(11):
    l3.append([])
for i in l2:
    j = (len(i)-1)/2
    l3[j].append(polRing(i))
l3[0].append(polRing(3))
l4 = []
for i in range(11):
    for j in l[i]:
        for k in l3[10-i]:
            t = j*k
            u = t
            while not u(1):
                u //= x-1
            v = u(-1)
            if u.degree()%2:
                v *= 3
            if v.is_square():
                l4.append(t)
f = open("k3f3-ej-lines.txt", "wb")
for i in l4:
    f.write(str(list(i)))
    f.write("\n")
f.close()


