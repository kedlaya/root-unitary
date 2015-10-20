#%pushd ..
#load("prescribed_roots.sage")
#%popd
polRing.<x> = PolynomialRing(ZZ)
ans, count = roots_on_unit_circle(2*(x^20+1), 
                                       filter=no_roots_of_unity,
                                       num_threads=512)
f = open("f2-candidates.txt", "wb")
for i in ans:
    f.write(str(i.list()))
f.close()
