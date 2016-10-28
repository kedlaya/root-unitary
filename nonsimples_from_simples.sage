#non-simples from simples
load("table_functions.sage")
load("simple_maker.sage")
from itertools import product, combinations_with_replacement, imap

def multiplicity_dict(L):
    L = list(L)
    return {a: L.count(a) for a in set(L)}

def nonsimples_from_simples(g,q,rootdir=None):
    filename = "weil-all-g%s-q%s"
    p,r = q.is_prime_power(get_data=True)
    polyRingF.<t> = PolynomialRing(Qp(p))
    D = load_previous_polys(g, rootdir)
    for split in Partitions(g):
        split_mD = multiplicity_dict(split)
        it = product(*(combinations_with_replacement(D[a,q],c) for a,c in sorted(split_mD.items())))
        for factors in imap(flatten, it):
            Lpoly = product(poly for lbl,poly in factors)
            Lpoly_mD = multiplicity_dict([factor[0] for factor in factors])
            seen = set()
            labels_uniq = [lbl for lbl,poly in factors if not (lbl in seen or seen.add(lbl))] # keeps order
            decomposition = str([[quote_me(lbl), Lpoly_mD[lbl]] for lbl in labels_uniq]).replace(" ","")
            label = make_label(g, q, Lpoly)
            polydata = PolyData(g, q, p, r, label, decomposition, polyRingF)
            
            #need to compute a method for computing the invs and newton poly
            #return a single factor
            line = '['
            #label
            my_label = make_label(g,q,Lpoly)
            line = line + quote_me(my_label) + ','
            #dim
            line = line + str(g) + ','
            #q
            line = line + str(q) + ','
            #coeffs
            line = line + str(coeffs) + ','
            #angle numbers
            line = line + str(sorted(angles(Lpoly))) + ','
            #angle ranks
            line = line + '"",'
            #number_field label
            #line = line + WebNumberField.from_polynomial(Ppoly).label + ','
            #FIXME THE FOLLOWING IS STUPID
            slopes,p_rank = newton_and_prank(p,r,Ppolyt)
            #line = line + ("%s" % slopes) + ','
            #p-rank
            line = line + str(p_rank) + ','
            #slopes
            s = [quote_me(x) for x in str(sorted(slopes))]
            line = line + str(s) + ','
            #a-counts
            a_counts = abelian_counts(g,p,r,Lpoly)
            line = line + str(a_counts) + ','
            #c-counts
            c_counts = curve_counts(g,q,Lpoly)
            line = line + str(c_counts) + ','
            #known jacobian
            line = line + '0,'
            #principally polarizable
            line = line + '0,'
            #decomposition
            factor_labels = [quote_me(factor[0]) for factor in factors]
            line = line + '[' + quote_me(factor_labels) + '],'
            #brauers
            line = line + '" "' + ','
            #primitive models
            line = line + '""'
            target.write(line + ']' + '\n')
            #The new L-function is the product of the L_i where the L_i is in the list associated to r_i,q we need to extract the label.   
