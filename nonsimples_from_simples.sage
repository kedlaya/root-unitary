#non-simples from simples
load("table_functions.sage")
load("simple_maker.sage")
from itertools import combinations_with_replacement, imap
import shutil

def multiplicity_dict(L):
    L = list(L)
    return {a: L.count(a) for a in set(L)}

class CombosWithReplacement(object):
    def __init__(self, L, m):
        self.L = L
        self.m = m
    def __len__(self):
        return binomial(len(self.L)+self.m-1,self.m)
    def __iter__(self):
        return combinations_with_replacement(self.L, self.m)

def nonsimples_from_simples(g,q,D=None,rootdir=None):
    if rootdir is None:
        rootdir = os.path.abspath(os.curdir)
    print "Generating from simples g=%s, q=%s"%(g, q)
    simplefilename = os.path.join(rootdir, "weil-simple-g%s-q%s.txt"%(g, q))
    filename = os.path.join(rootdir, "weil-all-g%s-q%s.txt"%(g, q))
    shutil.copy(simplefilename, filename)
    p,r = q.is_prime_power(get_data=True)
    polyRingF.<t> = PolynomialRing(Qp(p))
    if D is None:
        D = {}
        for a in range(1,g):
            D.update(load_previous_polys(q, a, rootdir))
    with open(filename, 'a') as target:
        for split in Partitions(g):
            if len(split) == 1:
                continue
            split_mD = multiplicity_dict(split)
            #for a,c in sorted(split_mD.iteritems()):
            #    print a, c, len(list(combinations_with_replacement(D[a,q],c)))
            it = cartesian_product_iterator([CombosWithReplacement(D[a,q],c) for a,c in sorted(split_mD.items())])
            for factors in it:
                if len(factors) == 1:
                    factors = factors[0]
                else:
                    factors = sum(factors,())
                Lpoly = prod(poly for lbl,poly in factors)
                if Lpoly.degree() != 2*g:
                    print factors
                    raise RuntimeError
                Lpoly_mD = multiplicity_dict([factor[0] for factor in factors])
                seen = set()
                labels_uniq = [lbl for lbl,poly in factors if not (lbl in seen or seen.add(lbl))] # keeps order
                decomposition = [[lbl, Lpoly_mD[lbl]] for lbl in labels_uniq]
                label = make_label(g, q, Lpoly)
                polydata = PolyData(g, q, p, r, label, decomposition, polyRingF)
                target.write(create_line(Lpoly, polydata, False))

def all_nonsimples_from_simples(qstart=None, gstart=None, rootdir=None):
    if rootdir is None:
        rootdir = os.path.abspath(os.curdir)
    gs = defaultdict(list)
    for filename in os.listdir(rootdir):
        match = gmatcher.match(filename)
        if match:
            gf, qf = map(Integer, match.groups())
            gs[qf].append(gf)
    for q in gs:
        if qstart is not None and q < qstart:
            continue
        gs[q].sort()
        if gs[q][-1] > 1:
            D = {}
            for a in range(1,gs[q][-1]):
                D.update(load_previous_polys(q, a, rootdir))
        for g in gs[q]:
            if gstart is not None:
                if g < gstart:
                    continue
                gstart = None
            if g == 1:
                simplefilename = "weil-simple-g%s-q%s.txt"%(g, q)
                filename = "weil-all-g%s-q%s.txt"%(g, q)
                shutil.copy(simplefilename, filename)
            else:
                nonsimples_from_simples(g, q, D, rootdir)
