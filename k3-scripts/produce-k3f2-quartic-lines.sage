load("prescribed_roots.sage")
serRing.<x> = PowerSeriesRing(QQ, 50)
f = open("k3points.txt", "rb")
ans = []
for i in f:
    l = i.split()
    t = sum(int(l[j])*x^j/j for j in range(1,11))
    t = exp(t)
    t *= (1 - x)*(1-2*x)*(1-4*x)
    t = t.inverse()
    m = [2^(1-j)*t[j] for j in range(11)]
    
[2, 7, 14, 23, 29, 28, 23, 15, 6, 2, 1, -1, -2, -6, -15, -23, -28, -29, -23, -14, -7, -2]
