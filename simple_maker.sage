load("table_functions.sage")
load("prescribed_roots.sage")
from sage.databases.cremona import cremona_letter_code
from sage.databases.cremona import class_to_int
from collections import namedtuple
import json
PolyData = namedtuple("PolyData","g q p r label decomposition polyRingF")
#from lmfdb.WebNumberField import WebNumberField

def create_line(Lpoly, polydata, check_simple = True):
    Ppoly = Lpoly.reverse()
    coeffs = Lpoly.coefficients(sparse=False)
    Ppolyt = polydata.polyRingF(Ppoly) #p-adic (automatically changes variable from x to t)
    if check_simple:
        factors = Lpoly.factor()
        if len(factors) != 1:
            return ""
        factor, power = factors[0]
        invs, newton_slopes = find_invs_and_slopes(polydata.p, polydata.r, factor.reverse())
        e = lcm([a.denominator() for a in invs])
        if e != power:
            return ""
        invs = [str(inv) for inv in invs]
    else:
        invs = ""
    my_label = make_label(polydata.g, polydata.q, Lpoly)
    slopes, p_rank = newton_and_prank(polydata.p, polydata.r, Ppolyt)
    line = [my_label, # label
            int(polydata.g), # dim
            int(polydata.q), # q
            map(int, coeffs), # polynomial coefficients
            map(float, sorted(angles(Lpoly))), # angle numbers
            "", # angle ranks
            p_rank, # p-rank
            [str(s) for s in sorted(slopes)], # slopes
            map(int, abelian_counts(polydata.g, polydata.p, polydata.r, Lpoly)), # A_counts
            map(int, curve_counts(polydata.g, polydata.q, Lpoly)), # C_counts
            0r, # known jacobian
            0r, # principally polarizable
            polydata.decomposition, # decomposition
            invs, # brauer invariants
            ""] # primitive models
    return json.dumps(line) + "\n"

@cached_function
def symfunc(i, r):
    Sym = SymmetricFunctions(QQ)
    p = Sym.powersum()
    if i == 0:
        return p.one()
    e = Sym.elementary()
    return e(p(e[i]).map_support(lambda A: Partition([r*c for c in list(A)])))

@cached_function
def basechange_transform(g, r, q):
    f = [symfunc(i, r) for i in range(g+1)]
    coeffs = [b.coefficients() for b in f]
    exps = [[{a: list(elem).count(a) for a in set(elem)} for elem in sorted(b.support()) if list(elem) and max(elem) <= 2*g] for b in f]
    def bc(Lpoly):
        # Assume that Lpoly has constant coefficient 1.
        R = Lpoly.parent()
        signed_coeffs = [(-1)^j * c for j, c in enumerate(Lpoly)]
        bc_coeffs = [1]
        for i in range(1, g+1):
            bc_coeffs.append((-1)^i*sum(c*prod(signed_coeffs[j]^e for j,e in D.iteritems()) for c, D in zip(coeffs[i], exps[i])))
        for i in range(1,g+1):
            # a_{g+i} = q^(ri) * a_{g-i}
            bc_coeffs.append(q^(r*i) * bc_coeffs[g-i])
        return R(bc_coeffs)
    return bc

def base_change(Lpoly, r, algorithm='sym', g = None, q = None, prec=53):
    if g is None:
        g = Lpoly.degree()
        assert g % 2 == 0
        g = g // 2
    if q is None:
        q = Lpoly.leading_coefficient().nth_root(g)
    if algorithm == 'approx':
        C = ComplexField(prec)
        R = RealField(prec)
        LC = Lpoly.change_ring(C)
        x = LC.parent().gen()
        approx = prod((1 - x/alpha^r)^e for alpha, e in LC.roots())
        approx_coeffs = approx.list()
        acceptable_error = R(2)^-(prec//2)
        exact_coeffs = [c.real().round() for c in approx_coeffs]
        if max(abs(ap - ex) for ap, ex in zip(approx_coeffs, exact_coeffs)) > acceptable_error:
            raise RuntimeError
        return Lpoly.parent()(exact_coeffs)
    else:
        return basechange_transform(g, r, q)(Lpoly)

def make_simples(g, q):
    print "starting g=%s, q=%s..." % (g,q)
    p,r = q.is_prime_power(get_data=True)

    #open a file for writing
    filename = "weil-simple-g%s-q%s.txt"%(g, q)

    polyRing.<x> = PolynomialRing(ZZ)
    polyRingF.<t> = PolynomialRing(Qp(p)) #p-adic version

    Lpolys,_ = roots_on_unit_circle( 1+ (q*x^2)^g )
    npolys = len(Lpolys)
    with open(filename, 'w') as target:
        for polycount,Lpoly in enumerate(Lpolys):
            if polycount % 500 == 0 and polycount != 0:
                print "%s / %s polynomials complete" % (polycount, npolys)
            label = make_label(g, q, Lpoly)
            polydata = PolyData(g, q, p, r, label, [[label, 1]], polyRingF)
            target.write(create_line(Lpoly, polydata))

def revise_file(g, q):
    print "revising %s, %s"%(g, q)
    p,r = ZZ(q).is_prime_power(get_data=True)
    oldfilename = "weil-%s-%s.txt"%(g, q)
    newfilename = "weil-simple-g%s-q%s.txt"%(g, q)
    s = ""
    polyRing.<x> = PolynomialRing(ZZ)
    polyRingF.<t> = PolynomialRing(Qp(p)) #p-adic version
    with open(oldfilename) as F:
        for line in F.readlines():
            try:
                data = json.loads(line.strip())
                label = data[0]
                data[12] = [[label, 1r]] # decomposition
                s += json.dumps(data) + "\n"
            except ValueError:
                Lpoly = polyRing(json.loads(line[line.find('[',1):line.find(']')+1]))
                label = make_label(g, q, Lpoly)
                polydata = PolyData(g, q, p, r, label, [[label, 1]], polyRingF)
                s += create_line(Lpoly, polydata)
    with open(newfilename,'w') as F:
        F.write(s)

def revise_files(rootdir=None):
    if rootdir is None:
        rootdir = os.path.abspath(os.curdir)
    for filename in os.listdir(rootdir):
        match = oldmatcher.match(filename)
        if match:
            gf, qf = map(int,match.groups())
            revise_file(gf, qf)

def fill_primitive_models(rootdir=None):
    # Assumes that, for fixed g, if q^r exists then all smaller powers do as well.
    if rootdir is None:
        rootdir = os.path.abspath(os.curdir)
    models = defaultdict(list)
    qs = defaultdict(list)
    for filename in os.listdir(rootdir):
        match = allmatcher.match(filename)
        if match:
            gf, qf = map(int, match.groups())
            qs[gf].append(qf)
    D = {}
    for g in qs:
        qs[g].sort()
        for q in qs[g]:
            if q^2 in qs[g]:
                D.update(load_previous_polys(q,g,rootdir,True))
                r = 2
                while q^r in qs[g]:
                    for label, Lpoly in D[g, q]:
                        if Lpoly in models: # already a base change
                            continue
                        models[base_change(Lpoly, r)].append(label)
                    r += 1
    R = PolynomialRing(QQ,'x')
    for g in qs:
        for q in qs[g]:
            filename = os.path.join(rootdir, "weil-all-g%s-q%s.txt"%(g, q))
            s = ""
            if g == 5 and q == 3:
                continue
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
