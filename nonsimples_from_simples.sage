#####################################################
#### Process for computing and uploading new data ###
#####################################################
## In the following, I assume that the root-unitary
## and lmfdb folders are contained in a common root.
# `cd /PATH/TO/lmfdb`
# `./warwick.sh &`
# `sage`
## From the sage prompt:
# run `import lmfdb`
# `cd ../root-unitary`
# run `load('nonsimples_from_simples.sage')`
# run `make_simples(g, q)` for new (g, q=p^m)
## this will create new weil-simple-g%-q%.txt
## Make sure that you still have the files
## weil-simple-g%-q%.txt for smaller g (same q) and
## weil-all-g%-q%.txt for divisors of q (same g)
# run `nonsimples_from_simples(g, q)`
## this will create new weil-all-g%-q%.txt
## Alternatively, you can run
## all_nonsimples_from_simples with a range of (g, q)
# create dictionary qs with
# `qs[g] = [p^r for r in range(1,m+1)]`
# create dictionary models with
# `models = _make_models_dict(qs)`
# run `_fill_primitive_models(models, qs)`
## this will fill in the primitive models field
## Add more g and q to qs to fill more models
# run `_fill_angle_ranks(g, q)`
## this will fill in the angle_ranks field
# run `load("../lmfdb/lmfdb/import_scripts/abvar/fq_isog.py")`
# run `do_import_all('weil-all-g%s-q%s.txt'%(g, q))`

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
                Lpoly = prod(factor.poly for factor in factors)
                if Lpoly.degree() != 2*g:
                    print factors
                    raise RuntimeError
                Lpoly_mD = multiplicity_dict([factor.label for factor in factors])
                seen = set()
                labels_uniq = [factor.label for factor in factors if not (factor.label in seen or seen.add(factor.label))] # keeps order
                decomposition = [[lbl, Lpoly_mD[lbl]] for lbl in labels_uniq]
                label = make_label(g, q, Lpoly)
                angle_numbers = sorted(sum(factor.angle_numbers for factor in factors, []))
                p_rank = sum(factor.p_rank for factor in factors)
                slopes = [str(a) for a in sorted(sum([QQ(s) for s in factor.slopes] for factor in factors, []))]
                invs = sum(factor.invs, [])
                places = sum(factor.places, [])
                jac = know_nonsimple_jac(g, q, p, r, [factor.poly for factor in factors], p_rank)
                polydata = NonsimplePolyData(g, q, p, r, label, angle_numbers, p_rank, slopes, jac, decomposition, invs, places, polyRingF)
                target.write(create_line(Lpoly, polydata, False))

def all_nonsimples_from_simples(qstart=None, gstart=None, rootdir=None):
    if rootdir is None:
        rootdir = os.path.abspath(os.curdir)
    gs = defaultdict(list)
    for filename in os.listdir(rootdir):
        match = simplematcher.match(filename)
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

def _fill_angle_rank(g, q, rootdir=None):
    rootdir = rootdir or os.path.abspath(os.curdir)
    R = PolynomialRing(ZZ,'x')
    filename = 'weil-all-g%s-q%s.txt'%(g, q)
    print filename
    s = ""
    early_break = False
    with open(filename) as F:
        for line in F.readlines():
            data = json.loads(line.strip())
            Ppoly = R(list(reversed(data[3])))
            if data[5] != "":
                early_break = True
                break
            data[5] = int(num_angle_rank(Ppoly))
            s += json.dumps(data) + "\n"
    if not early_break:
        with open(filename, 'w') as F:
            F.write(s)

def fill_angle_ranks(rootdir=None):
    rootdir = rootdir or os.path.abspath(os.curdir)
    gs = defaultdict(list)
    for filename in os.listdir(rootdir):
        match = allmatcher.match(filename)
        if match:
            gf, qf = map(int, match.groups())
            gs[qf].append(gf)
    for q in sorted(gs.keys()):
        gs[q].sort()
        for g in gs[q]:
            _fill_angle_rank(g, q, rootdir)

def _fill_primitive_models(models, qs, rootdir=None):
    rootdir = rootdir or os.path.abspath(os.curdir)
    R = PolynomialRing(QQ,'x')
    for g in qs:
        for q in qs[g]:
            filename = os.path.join(rootdir, "weil-all-g%s-q%s.txt"%(g, q))
            s = ""
            print filename
            with open(filename) as F:
                for line in F.readlines():
                    data = json.loads(line.strip())
                    Lpoly = R(data[3])
                    if Lpoly in models:
                        data[14] = models[Lpoly]
                    else:
                        data[14] = [data[0]]
                    s += json.dumps(data) + "\n"
            try:
                with open(filename, 'w') as F:
                    F.write(s)
            except KeyboardInterrupt:
                with open(filename, 'w') as F:
                    F.write(s)
                raise

def _make_models_dict(qs, rootdir=None):
    rootdir = rootdir or os.path.abspath(os.curdir)
    models = defaultdict(list)
    D = {}
    for g in qs:
        qs[g].sort()
        for q in qs[g]:
            if q^2 in qs[g]:
                D.update(load_previous_polys(q,g,rootdir,True))
                r = 2
                while q^r in qs[g]:
                    for factor in D[g, q]:
                        if factor.poly in models: # already a base change
                            continue
                        models[base_change(factor.poly, r)].append(factor.label)
                    r += 1
    return models


def fill_primitive_models(rootdir=None):
    # Assumes that, for fixed g, if q^r exists then all smaller powers do as well.
    rootdir = rootdir or os.path.abspath(os.curdir)
    qs = defaultdict(list)
    for filename in os.listdir(rootdir):
        match = allmatcher.match(filename)
        if match:
            gf, qf = map(int, match.groups())
            qs[gf].append(qf)
    models = _make_models_dict(qs, rootdir)
    _fill_primitive_models(models, qs, rootdir)
