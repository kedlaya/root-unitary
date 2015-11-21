polRing.<x> = PolynomialRing(ZZ)
f = open("k3f2-full-filtered.txt", "rb")
l = [tuple(eval(i)) for i in f]
f.close()
print "Loaded", len(l), "candidate polynomials"
ls = set(l)
f = open("k3final.txt", "rb")
l2 = [tuple(sage_eval(i, locals={'x':x})[0].list()) for i in f]
f.close()
print "Loaded", len(l2), "realized polynomials"
l3 = [i for i in l2 if not i in ls]
print "Missing", len(l3), "candidates"

