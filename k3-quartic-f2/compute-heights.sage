polRing.<x> = PolynomialRing(Integers())
f = open("k3final.txt", "rb")
l = [sage_eval(i, locals={'x':x}) for i in f]
f.close()
ht1 = [0,] * 11
ht2 = [0,] * 11
for i in l:
    jj = 0
    for j in factor(i[0]):
        if abs(j[0][0]) == 2:
            jj = j[0].change_ring(GF(2))
    if jj != 0:
        k = jj.ord(x)-1
        ht1[k] += 1
        ht2[k] += len(i[1])
    else:
        ht1[10] += 1
        ht2[10] += len(i[1])
print ht1, ht2
