polRing.<x> = PolynomialRing(ZZ)
f = open("k3f1-list.txt", "rb")
l = sage_eval(f.read(), locals={'x':x})
f.close()

l2 = []
c = [cyclotomic_polynomial(i, var=x) for i in range(1,200) 
     if euler_phi(i) <= 21]
w = WeightedIntegerVectors([i.degree() for i in c])
for i in range(11):
    l2.append([])
    w1 = w.graded_component(2*i+1)
    for j in w1:
        l2[i].append(prod(c[k]^j[k] for k in range(len(j))))
    print set(l[i]) == set(l2[i])


