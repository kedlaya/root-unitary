import time
load("prescribed_roots.sage")

polRing.<x> = PolynomialRing(ZZ)
l = []
l.append(x+1)
l.append(x-1)
t = time.time()
for i in range(1, 11):
    ans, count = roots_on_unit_circle(x^(2*i)+1, num_threads=512)
    for j in ans:
        l.append(j * (x+1))
        l.append(j * (x-1))
    print len(ans)*2, "polynomials added"
    print "time so far: ", time.time() - t, " seconds"

f = open("k3-scripts/k3f1-lines.txt", "wb")
for i in l:
    f.write(str(i.list()))
    f.write("\n")
f.close()
