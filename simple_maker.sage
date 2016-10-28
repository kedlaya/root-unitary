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
    #if check_simple:
    factors = Lpoly.factor()
    if len(factors) != 1:
        return ""
    factor, power = factors[0]
    invs, newton_slopes = find_invs_and_slopes(polydata.p, polydata.r, factor.reverse())
    e = lcm([a.denominator() for a in invs])
    if e != power:
        return ""

    #return a single factor
    line = '['
    #label
    my_label = make_label(polydata.g, polydata.q, Lpoly)
    line = line + quote_me(my_label) + ','
    #dim
    line = line + str(polydata.g) + ','
    #q
    line = line + str(polydata.q) + ','
    #polynomial coefficients
    line = line + str(coeffs) + ','
    #angle numbers
    line = line + str(sorted(angles(Lpoly))) + ','
    #angle ranks
    line = line + '"",'
    #number_field label
    #line = line + WebNumberField.from_polynomial(Ppoly).label + ','
    slopes,p_rank = newton_and_prank(polydata.p, polydata.r, Ppolyt)
    #p-rank
    line = line + str(p_rank) + ','
    #slopes
    printed_slopes = '['
    for s in sorted(slopes):
        if printed_slopes != '[':
            printed_slopes += ','
        printed_slopes += quote_me(s)
    #s = [quote_me(x) for x in sorted(slopes)]
    line = line + printed_slopes + '],'
    #A_counts
    a_counts = abelian_counts(polydata.g, polydata.p, polydata.r, Lpoly)
    line = line + str(a_counts) + ','
    #C_counts
    c_counts = curve_counts(polydata.g, polydata.q, Lpoly)
    line = line + str(c_counts) + ','
    #known jacobian
    line = line + '0,'
    #principally polarizable
    line = line + '0,'
    #decomposition
    line = line + '[[' + quote_me(my_label) + ',1]],'
    #brauer invariants
    line += '['
    printed_invs = '['
    for inv in invs:
        if printed_invs != '[':
            printed_invs += ','
        printed_invs += quote_me(inv)
    #invs_str = [quote_me(inv) for inv in invs]
    line = line + printed_invs + '],'
    #primitive models
    line = line + '""'
    return line + ']\n'

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
            polydata = PolyData(g, q, p, r, label, '[[%s,1]]'%(quote_me(label)), polyRingF)
            target.write(create_line(Lpoly, polydata))

def revise_file(g, q):
    print "revsing %s, %s"%(g, q)
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
                data[12] = [[quote_me(label), 1r]]
                s += json.dumps(data) + "\n"
            except ValueError:
                Lpoly = polyRing(json.loads(line[line.find('[',1):line.find(']')+1]))
                label = make_label(g, q, Lpoly)
                polydata = PolyData(g, q, p, r, label, '[[%s,1]]'%(quote_me(label)), polyRingF)
                s += create_line(Lpoly, polydata, False)
    with open(newfilename,'w') as F:
        F.write(s)

def revise_files(rootdir=None):
    if rootdir is None:
        rootdir = os.path.abspath(os.curdir)
    for dirpath, dirnames, filenames in os.walk(rootdir):
        for filename in filenames:
            match = oldmatcher.match(filename)
            if match:
                gf, qf = map(int,match.groups())
                revise_file(gf, qf)
