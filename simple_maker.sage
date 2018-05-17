##############################################################
##################### Data specification #####################
##############################################################
# 0  label                   string
# 1  g                       int
# 2  q                       int
# 3  polynomial              list of ints
# 4  angle_numbers           list of floats
# 5  angle_ranks             int ("" if unfilled)
# 6  p_rank                  int
# 7  slopes                  list of strings
# 8  A_counts                list of ints
# 9  C_counts                list of ints
# 10 known_jacobian          int
# 11 principally_polarizable int
# 12 decomposition           list of pairs [string, int]
# 13 brauer_invariants       list of strings
# 14 places                  list of lists of lists of strings
# 15 primitive_models        list of strings ("" if unfilled)
# 16 number_field            string ("" if unknown/non-simple)
# 17 galois_n                int ("" if non-simple)
# 18 galois_t                int ("" if unknown/non-simple)

load("table_functions.sage")
load("prescribed_roots.sage")
from collections import namedtuple
import json
NonsimplePolyData = namedtuple("NonsimplePolyData","g q p r label angle_numbers p_rank slopes jac decomposition invs places polyRingF")
SimplePolyData = namedtuple("SimplePolyData","g q p r label polyRingF")
from lmfdb.WebNumberField import WebNumberField

def check_dataspec(filename):
    bad = {}
    with open(filename) as F:
        for line in F.readlines():
            if not line.strip():
                bad[100] = 'empty lines'
                continue
            data = json.loads(line.strip())
            if len(data) != 19:
                print "%s: %s has length %s"%(filename, data[0], len(data))
                break
            if not (isinstance(data[14], list) and len(data[14]) > 0 and
                    isinstance(data[14][0], list) and len(data[14][0]) > 0 and
                    isinstance(data[14][0][0], basestring)):
                bad[14] = data[0]
        else:
            if bad:
                for key in sorted(bad.keys()):
                    print "%s: %s"%(filename, bad[key])
            else:
                print filename, "ok"

def check_all(rootdir=None):
    rootdir = rootdir or os.path.abspath(os.curdir)
    for filename in os.listdir(rootdir):
        if simplematcher.match(filename):
            check_dataspec(filename)

def know_simple_pp(g, q, coeffs):
    if g == 1:
        return 1r
    if g == 2:
        # P-poly: x^4 = ax^3 + bx^2 + aqx + q^2
        # Howe, Maisner, Nart, Ritzenthaler
        # "Principally polarizable isogeny classes of abelian surfaces over finite fields"
        a = ZZ(coeffs[1])
        b = ZZ(coeffs[2])
        if (a^2 - b == q and b < 0 and all(p % 3 == 1 for p in b.prime_divisors())):
            return -1r
        else:
            return 1r
    return 0r

def know_nonsimple_pp(g, q, coeffs):
    if g == 2:
        return know_simple_pp(g, q, coeffs)
    return 0r

def know_simple_jac(g, q, p, r, coeffs, p_rank):
    if g == 1:
        return 1r
    if g == 2:
        # P-poly: x^4 + ax^3 + bx^2 + aqx + q^2
        # Howe, Nart, Ritzenthaler
        # "Jacobians in isogeny classes of abelian surfaces over finite fields"
        a = ZZ(coeffs[1])
        b = ZZ(coeffs[2])
        if (a^2 - b == q and b < 0 and all(p % 3 == 1 for p in b.prime_divisors())
            or p_rank == 2 and a == 0 and (b == 1-2*q or p > 2 and b == 2-2*q)
            or p_rank == 0 and a == 0 and b == -q and (p % 12 == 11 and r % 2 == 0 or
                                                       p == 3 and r % 2 == 0 or
                                                       p == 2 and r % 2 == 1)
            or p_rank == 0 and a == 0 and b == -2*q and (q == 2 or q == 3)):
            return -1r
        else:
            return 1r
    return 0r

def know_nonsimple_jac(g, q, p, r, factors, p_rank):
    if g == 2:
        # P-poly: (x^2 - sx + q)(x^2 - tx + q) with |s| >= |t|
        # Howe, Nart, Ritzenthaler
        # "Jacobians in isogeny classes of abelian surfaces over finite fields"
        s = -ZZ(factors[0][1])
        t = -ZZ(factors[1][1])
        if abs(t) > abs(s):
            s, t = t, s
        if (abs(s - t) == 1
            or p_rank == 2 and (s == t and t^2 - 4*q in [-3, -4, -7] or
                                q == 2 and abs(s) == abs(t) == 1 and s != t)
            or p_rank == 1 and r % 2 == 0 and s^2 == 4*q and (s-t).is_squarefree()
            or p_rank == 0 and (p > 3 and abs(s) != abs(t) or
                                p == 3 and r % 2 == 1 and s^2 == t^2 == 3*q or
                                p == 3 and r % 2 == 0 and (s - t) % (3*p^(r//2)) != 0 or
                                p == 2 and (s^2 - t^2) % (2*q) != 0 or
                                q in [2,3] and s == t or
                                q in [4,9] and s^2 == t^2 == 4*q)):
            return -1r
        else:
            return 1r
    return 0r

def create_line(Lpoly, polydata, simple = True):
    Ppoly = Lpoly.reverse()
    coeffs = Lpoly.coefficients(sparse=False)
    Ppolyt = polydata.polyRingF(Ppoly) #p-adic (automatically changes variable from x to t)
    if simple:
        factors = Ppoly.factor()
        if len(factors) != 1:
            return ""
        factor, power = factors[0]
        invs, places, newton_slopes = find_invs_places_and_slopes(polydata.p, polydata.r, factor)
        e = lcm([a.denominator() for a in invs])
        # When q is not a square, the case 1-qx^2 must be handled separately.
        if (not is_square(polydata.q)) and factor.degree() == 2 and factor[1] == 0 and factor[0] < 0:
            e = 2
        if e != power:
            return ""
        angle_numbers = map(float, sorted(angles(Lpoly)))
        slopes, p_rank = newton_and_prank(polydata.p, polydata.r, Ppolyt)
        slopes = [str(s) for s in sorted(slopes)]
        jac = know_simple_jac(polydata.g, polydata.q, polydata.p, polydata.r, coeffs, p_rank)
        pp = know_simple_pp(polydata.g, polydata.q, coeffs)
        decomposition = [[polydata.label, 1r]]
        invs = [str(inv) for inv in invs]
        places = [places]
        galois_n = int(factor.degree())
        nf = WebNumberField.from_polynomial(factor)
        if nf.label == 'a':
            nf_label = galois_t = ""
        else:
            nf_label = nf.label
            galois_t = nf.galois_t()
    else:
        angle_numbers = polydata.angle_numbers
        p_rank = polydata.p_rank
        slopes = polydata.slopes
        jac = polydata.jac
        pp = know_nonsimple_pp(polydata.g, polydata.q, coeffs)
        decomposition = polydata.decomposition
        invs = polydata.invs
        places = polydata.places
        nf_label = galois_t = galois_n = ""
    C_counts = map(int, curve_counts(polydata.g, polydata.q, Lpoly))
    if any(c < 0r for c in C_counts):
        if jac == 1r:
            raise RuntimeError("Incorrect value for Known Jacobian")
        jac = -1r
    line = [polydata.label, # label
            int(polydata.g), # dim
            int(polydata.q), # q
            map(int, coeffs), # polynomial coefficients
            angle_numbers, # angle numbers
            "", # angle ranks
            int(p_rank), # p-rank
            slopes, # slopes
            map(int, abelian_counts(polydata.g, polydata.p, polydata.r, Lpoly)), # A_counts
            C_counts, # C_counts
            jac, # known jacobian
            pp, # principally polarizable
            decomposition, # decomposition
            invs, # brauer invariants
            places, # corresponding ideals
            "", # primitive models
            nf_label, # number field
            galois_n, galois_t] # galois group
    return json.dumps(line) + "\n"

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
            polydata = SimplePolyData(g, q, p, r, label, polyRingF)
            target.write(create_line(Lpoly, polydata))

def regen(g, q, rootdir):
    print "Regenerating g=%s q=%s"%(g, q)
    rootdir = rootdir or os.path.abspath(os.curdir)
    filename = os.path.join(rootdir, "weil-simple-g%s-q%s.txt"%(g, q))
    p,r = Integer(q).is_prime_power(get_data=True)
    polyRing.<x> = PolynomialRing(ZZ)
    polyRingF.<t> = PolynomialRing(Qp(p)) #p-adic version
    s = ""
    counter = 0
    with open(filename) as F:
        for line in F.readlines():
            data = json.loads(line.strip())
            label = data[0]
            if counter % 100 == 0 and counter != 0:
                print "%s (%s)"%(label, counter)
            counter += 1
            Lpoly = polyRing(data[3])
            polydata = SimplePolyData(g, q, p, r, label, polyRingF)
            s += create_line(Lpoly, polydata)
    with open(filename, 'w') as F:
        F.write(s)

def regen_qrange(qmin, qmax, rootdir=None):
    rootdir = rootdir or os.path.abspath(os.curdir)
    for filename in os.listdir(rootdir):
        match = simplematcher.match(filename)
        if match:
            gf, qf = map(int, match.groups())
            if qf >= qmin and qf <= qmax:
                regen(gf, qf, rootdir)

def regen_all(rootdir=None):
    rootdir = rootdir or os.path.abspath(os.curdir)
    for filename in os.listdir(rootdir):
        match = simplematcher.match(filename)
        if match:
            gf, qf = map(int, match.groups())
            if gf == 5 and qf == 3:
                continue
            regen(gf, qf, rootdir)

# def revise_file(g, q):
#     print "revising %s, %s"%(g, q)
#     p,r = ZZ(q).is_prime_power(get_data=True)
#     oldfilename = "weil-%s-%s.txt"%(g, q)
#     newfilename = "weil-simple-g%s-q%s.txt"%(g, q)
#     s = ""
#     polyRing.<x> = PolynomialRing(ZZ)
#     polyRingF.<t> = PolynomialRing(Qp(p)) #p-adic version
#     with open(oldfilename) as F:
#         for line in F.readlines():
#             try:
#                 data = json.loads(line.strip())
#                 label = data[0]
#                 data[12] = [[label, 1r]] # decomposition
#                 s += json.dumps(data) + "\n"
#             except ValueError:
#                 Lpoly = polyRing(json.loads(line[line.find('[',1):line.find(']')+1]))
#                 label = make_label(g, q, Lpoly)
#                 polydata = SimplePolyData(g, q, p, r, label, polyRingF)
#                 s += create_line(Lpoly, polydata)
#     with open(newfilename,'w') as F:
#         F.write(s)

# def revise_files(rootdir=None):
#     rootdir = rootdir or os.path.abspath(os.curdir)
#     for filename in os.listdir(rootdir):
#         match = oldmatcher.match(filename)
#         if match:
#             gf, qf = map(int,match.groups())
#             revise_file(gf, qf)

# def fill_nf_data(rootdir=None):
#     rootdir = rootdir or os.path.abspath(os.curdir)
#     R = PolynomialRing(QQ,'x')
#     def insert_data(data, newdata):
#         data[10] = newdata[0]
#         data[13:14] = newdata[1:3]
#         data.extend(newdata[3:])
#     for filename in os.listdir(rootdir):
#         match = simplematcher.match(filename)
#         if match:
#             g, q = map(int,match.groups())
#             if g == 6 and q == 2:
#                 continue
#             print filename
#             p,r = Integer(q).is_prime_power(get_data=True)
#             D = {}
#             s = ""
#             done_already = False
#             with open(filename) as F:
#                 for line in F.readlines():
#                     data = json.loads(line.strip())
#                     if len(data) == 19:
#                         done_already = True
#                         print "Skipping"
#                         break
#                     label = data[0]
#                     print label
#                     Ppoly = R(data[3]).reverse()
#                     C_counts = data[9]
#                     if any(c < 0 for c in C_counts):
#                         jac = -1r
#                     else:
#                         jac = 0r
#                     factor, power = Ppoly.factor()[0]
#                     galois_n = int(factor.degree())
#                     nf = WebNumberField.from_polynomial(factor)
#                     if nf.label == 'a':
#                         nf_label = galois_t = ""
#                     else:
#                         nf_label = nf.label
#                         galois_t = int(nf.galois_t())
#                     invs, places, slopes = find_invs_places_and_slopes(p, r, factor)
#                     newdata = (jac, [str(inv) for inv in invs], [places], nf_label, galois_n, galois_t)
#                     D[label] = newdata
#                     insert_data(data, newdata)
#                     s += json.dumps(data) + '\n'
#             if done_already:
#                 continue
#             with open(filename, 'w') as F:
#                 F.write(s)
#             s = ""
#             filename = "weil-all-g%s-q%s.txt"%(g, q)
#             print filename
#             with open(filename) as F:
#                 for line in F.readlines():
#                     data = json.loads(line.strip())
#                     label = data[0]
#                     print label
#                     if label in D:
#                         insert_data(data, D[label])
#                     else:
#                         C_counts = data[9]
#                         if any(c < 0 for c in C_counts):
#                             jac = -1r
#                         else:
#                             jac = 0r
#                         insert_data(data, (jac, "", "", "", "", ""))
#                     s += json.dumps(data) + '\n'
#             with open(filename, 'w') as F:
#                 F.write(s)

def clean(qmin, qmax, rootdir=None):
    rootdir = rootdir or os.path.abspath(os.curdir)
    qs = set()
    for q in srange(qmin, qmax+1):
        p, r = q.is_prime_power(get_data=True)
        if r != 0: # prime power
            revise_file(1,q)
            for a in range(1,r+1):
                qs.add(p^a)
    qs = {1r: sorted(list(qs))}
    models = _make_models_dict(qs, rootdir)
    _fill_primitive_models(models, qs, rootdir)

def fix_places_and_jac_pp(rootdir=None):
    rootdir = rootdir or os.path.abspath(os.curdir)
    for filename in os.listdir(rootdir):
        match = simplematcher.match(filename)
        if match:
            print filename
            g, q = map(Integer, match.groups())
            p, r = q.is_prime_power(get_data=True)
            # if g = 2
            # q power of 3 or q=r^2 where r=2 (mod 3)
            # then pp is True
            s = ""
            with open(filename) as F:
                for line in F.readlines():
                    data = json.loads(line.strip())
                    assert len(data) == 19
                    coeffs = data[3]
                    p_rank = data[6]
                    data[10] = know_simple_jac(g, q, p, r, coeffs, p_rank)
                    data[11] = know_simple_pp(g, q, coeffs)
                    data[14] = [data[14]]
                    s += json.dumps(data) + '\n'
            with open(filename, 'w') as F:
                F.write(s)
