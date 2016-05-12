
"""
Functions for generating the LMFDB data on isogeny classes of abelian varieties over FF_q

AUTHORS:
  -- (2016-05-11) Taylor Dupuy, Kiran S. Kedlaya, David Roe, Christelle Vincent 


Fields we want to populate with an example

label: "2.9.12.20"
polynomial: ["1","-1","3","-9","81"]
angle_numbers: (doubles): [0.23756..., 0.69210...]
number_field: "4.0.213413.1"
p-rank: 1
slopes: ["0","1/2","1/2","1"]
A_counts: ["75", "7125"]
C_counts: ["9", "87"]
known_jacobian (0,1,-1): 1
decomposition: ["9.2.-1._3"]
pricipally_polarizable (0,1,-1): 1
Brauer Invariants: inv_v( End(A_{FFbar_q})_{QQ} )=(v(\pi)/v(q))*[QQ(pi)_{v}: QQ(pi): v\vert p place of QQ(\pi)], these are stored as elements of QQ.
Primative models: 
"""

######################################################################################################

#for the make_label function
from sage.databases.cremona import cremona_letter_code
from sage.databases.cremona import class_to_int

load("prescribed_roots.sage")
#this command should be replaced in the big list
#polyx =  x^4-x^3+3*x^2-9*x+81

#create p-adic and symbolic versions
#coeffs = poly.coefficients(sparse=False)
#polyt = polyx(t) #p-adic 
#polyz = 0 #symbolic
#for i in range(2*g+1):
#    polyz = polyz + coeffs[i]*z^i

######################################################################################################

def newton_and_prank(p,r,poly_in_t):
    """
    Newton Polygons
    Needs:
    F = Qp(p)
    polyRingF.<t> = PolynomialRing(F)
    """
    d = poly_in_t.degree()
    mynp=poly_in_t.newton_polygon()
    slopes = []
    for i in range(d):
        rise = (-mynp(i+1) + mynp(i))/r
        slopes.append(rise)
    p_rank = slopes.count(0)
    return slopes,p_rank


def angles(poly_in_z):
    #returns the angles (as a multiple of pi) in the upperhalf plane
    roots = poly_in_z.roots(CC)
    return sum([[root.arg()/pi]*e for (root, e) in roots if root.imag() >= 0],[])


#def make_label(g,q,poly_in_x):
#    """
#    Needs: poly.<x> = PolynomialRing(ZZ)
#
#    we are going to replace the a_k label with our new labels 
#        l_k = a_k - c_k/((-1)^k k) + q^k 2g/k/
#    What are these crazy numbers?
#    We are using Newton's identifies: 
#        ka_k = \sum_{i=0}^k (-1)^{i-1} a_{k-1} s_i 
#    we are which that a_k in an interval around some c_k of size the Lang-Weil bound (LW/k = q^{k/2
#    2g/k ) on the coefficients (also the trivial bound on the power sums s_k)
#        c_k = sum_{i=0}^{k-1} (-1)^{k-1} a_{k-i}s_i.
#    We then take all the integers in the interval |x - c_k| < LW and relabel them with element
#    [0,2*LW]. 
#    This is done by shifting a_k - c_k by LW/k.
#    """
#    a=poly_in_x.coefficients(sparse=False)
#    s=[2*g]+[0 for i in range(2*g)]
#    c =[0 for i in range(2*g)]
#    #computes the power sums
#    for i in range(1,2*g):
#        cilist = [ (-1)^(j+i)*a[i-j]*s[j] for j in range(1,i)]
#        c[i] = sum( cilist )
#        s[i] = ((-1)^(i)*i*a[i]-c[i])/a[0]
#    #computes the new labels
#    #TEST ME (in our test case we get 12,22)
#    l = [0 for i in range(0,g+1)]
#    for i in range(1,g+1):
#        l[i] = ceil(a[i] - (-1)^i*c[i]/i + q^(i/2+g)*(2*g)/i)
#    label = [l[i] for i in range(1,g+1)]
#    return a,s,c,label 


def abelian_counts(g,p,r,L):
    return [L.resultant(x^i-1) for i in range(1,g+1)]

def curve_counts(g,q,L):
    S = PowerSeriesRing(QQ, 'x', g+1)
    g = S(L)/((1-x)*(1-q*x))
    return g.log().derivative().coefficients()

def alternating(pol, m):
    """
    This appears to take forever but there is a precomputed version elsewhere. 
    """
    d = pol.degree()
    pl = pol.list()
    e = SymmetricFunctions(QQ).e()
    dm = binomial(d, m)
    l = [(-1)^i*e[i](e[m]).restrict_parts(d) for i in range(dm+1)]
    P = pol.parent()
    R = P.base_ring()
    ans = []
    for i in range(dm,-1,-1):
        s = R.zero()
        u = tuple(l[i])
        for j, c in u:
            s += R(c) * prod(pl[d-k] for k in j)
        ans.append(s)
    return P(ans)
    
def find_invs_and_slopes(p,r,L):
    poly = L.change_ring(QQ)
    K.<a> = NumberField(poly)
    l = K.primes_above(p)
    invs = []
    slopes = []
    for v in l:
        vslope = a.valuation(v)/K(p^r).valuation(v)
        slope.append(vslope)
        vdeg = v.residue_class_degree()*v.ramification_index()
        invs.append(vslope*vdeg)
    return invs,slopes
        
def make_label(coeffs):
    


#######################################################################################################
#######################################################################################################

def make_table(g,q):
    """
    For a dimension g and a prime power q, generate a file "weilgp.txt" which contains the possible
    weil polynomials.
    FIXM: We need to pass an intelligent answer_count=1000
    """
    p,r = q.is_prime_power(get_data=True)
    currentfilename = "weil" + "-" + ( "%s" % g ) + "-"+ ("%s" % q) + ".txt"
    target = open(currentfilename, 'w')
    
    polyRing.<x> = PolynomialRing(ZZ)
    F = Qp(p) #p-adic version
    polyRingF.<t> = PolynomialRing(F)
    
    weil_polys,some_number_i_dont_understand = roots_on_unit_circle(1+(q*x^2)^g)
    for Lpoly in weil_polys:
        #create the two types of polynomial objects needed for algorithms 
        Ppoly = Lpoly.reverse()
        coeffs = Lpoly.coefficients(sparse=False)
        Ppolyt = polyRingF(Ppoly) #p-adic (automatically changes variable from x to t) 
        
        line_for_file = ""
        
        #label
        #FOR NEWTON LABELING CODE:
        #a,s,c,label = make_label(g,q,Ppoly)
        line_for_file = line_for_file + ("%s" % g) + "." + ("%s" % q)
        
        #for i in range(g):
        #    line_for_file = line_for_file + (".%s" % label[i]) 
        #line_for_file = line_for_file + ","
        
        
        #FOR NAIVE LABELING: 
        #convert to base 26 and use a leading 'a' for negative signs. 
        for i in range(1,g+1):
            c = coeffs[i]
            if sign(c) == -1:
                line_for_file = line_for_file + ".a" + cremona_letter_code((-1)*c)
            else
                line_for_file = line_for_file + "." + cremona_letter_code(c)
        line_for_file = line_for_file + ","
        
        #polynomial_coeffs
        line_for_file = line_for_file + ("%s" % coeffs) + ","
                
        #FIXME: number_field
        line_for_file = line_for_file + ","
        
        #slopes and p-rank
        slopes,p_rank = newton_and_prank(p,r,Ppolyt)
        line_for_file = line_for_file + ("%s" % slopes) + "," + ("%s" % p_rank)
        
        #A-counts
        a_counts = abelian_counts(g,p,r,Lpoly)
        line_for_file = line_for_file + ("%s" % a_counts) + "," 
        
        #C-counts
        c_counts = curve_counts(g,q,Lpoly)
        line_for_file = line_for_file + ("%s" % c_counts) + "," 
        
        #Known Jacobian
        line_for_file = line_for_file + "0" + ","
        
        #Has principal polarization
        
        #Check_Invs
        
        target.write(line_for_file + "\n")
####################################################################################################################################################

