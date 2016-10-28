#non-simples from simples
load("table_functions.sage")
load("simple_maker.sage")
from itertools import product, combinations_with_replacement, imap
import shutil

def multiplicity_dict(L):
    L = list(L)
    return {a: L.count(a) for a in set(L)}

def nonsimples_from_simples(g,q,rootdir=None):
    simplefilename = "weil-simple-g%s-q%s.txt"%(g, q)
    filename = "weil-all-g%s-q%s.txt"%(g, q)
    shutil.copy(simplefilename, filename)
    p,r = q.is_prime_power(get_data=True)
    polyRingF.<t> = PolynomialRing(Qp(p))
    D = load_previous_polys(g, rootdir)
    with open(filename, 'a') as target:
        for split in Partitions(g):
            if len(split) == 1:
                continue
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
                target.write(create_line(Lpoly, polydata, False))
