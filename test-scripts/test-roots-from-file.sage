def test_roots_from_file(filename, q):
    t = True
    threshold = 0.0001 * q
    polRing.<x> = PolynomialRing(ZZ)
    f = open(filename, "rb")
    for i in f:
        j = polRing(eval(i))
        l = j.roots(ComplexField())
        for k in l:
            if abs(abs(k[0])^2 - q) > threshold:
                t = False
                break    
    f.close()
    return(t)
