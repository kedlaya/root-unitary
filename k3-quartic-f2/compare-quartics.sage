polRing.<x> = PolynomialRing(ZZ)
n = 19
print n
f = open("k3f2-lines.txt", "rb")
la0 = [polRing(eval(i)) for i in f]
f.close()
la1 = [i for i in la0 if i.degree() == 20]
la2 = [i for i in la1 if (i*(1+x)).list()[11:n+1] != [0]*(n-10)]
la3 = [i for i in la2 if i(1) == 2]
#i(1) > 0 and i(1)%2==0 and (i(1)//2).is_square()]
la4a = [i for i in la3 if i(-1) != 2]
la4 = [i for i in la3 if i(-1) > 2 and (i(-1)//2).is_squarefree()]
la = [i for i in la4 if i[1]%2==1]
print len(la1), len(la2), len(la3), len(la4a), len(la4), len(la)
f = open("k3final.txt", "rb")
lb0 = [polRing(sage_eval(i, locals={'x':x})[0]) for i in f]
f.close()
lb1 = [i//(1+x) for i in lb0 if i(-1)==0 and (i//(1+x)).is_irreducible()]
lb2 = [i for i in lb1 if (i*(1+x)).list()[11:n+1] != [0]*(n-10)]
lb3 = [i for i in lb2 if i(1) == 2]
#i(1) > 0 and i(1)%2==0 and (i(1)//2).is_square()]
lb4a = [i for i in lb3 if i(-1) != 2]
lb4 = [i for i in lb3 if i(-1) > 2 and (i(-1)//2).is_squarefree()]
lb = [i for i in lb4 if i[1]%2==1]
print len(lb1), len(lb2), len(lb3), len(lb4a), len(lb4), len(lb)

