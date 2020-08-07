import time
polRing.<x> = PolynomialRing(ZZ)
l = []
l.append(1+x)
l.append(1-x)
t = time.time()
for i in range(1, 11):
    wp = WeilPolynomials(2*i, 1, 1, parallel=True)
    ans = list(wp)
    for j in ans:
        l.append(j * (1+x))
        l.append(j * (1-x))
    print(len(ans)*2, "polynomials added")
    print("time so far: ", time.time() - t, " seconds")

with open("k3f1-lines.txt", "w") as f:
    for i in l:
        f.write(str(i.list()))
        f.write("\n")

