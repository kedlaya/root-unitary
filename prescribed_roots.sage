"""
Finding polynomials with roots in prescribed regions

AUTHOR:
  -- Kiran S. Kedlaya (2007-05-28): initial version
  -- Kiran S. Kedlaya (2015-08-29): updated version; switch from NTL to FLINT
        
EXAMPLES:
    sage: from prescribed_roots import roots_on_unit_circle
    sage: polRing.<x> = PolynomialRing(Rationals())
    sage: P0 = 3*x^21 + 5*x^20 + 6*x^19 + 7*x^18 + 5*x^17 + 4*x^16 + 2*x^15 - x^14 - 3*x^13 - 5*x^12 - 5*x^11 - 5*x^10 - 5*x^9 - 3*x^8 - x^7 + 2*x^6 + 4*x^5 + 5*x^4 + 7*x^3 + 6*x^2 + 5*x + 3
    sage: ans, count = roots_on_unit_circle(P0, 3^2, 1)
    sage: print "Number of terminal nodes:", count
    Number of terminal nodes: 1157
    sage: print ans
    [3*x^21 + 5*x^20 + 6*x^19 + 7*x^18 + 5*x^17 + 4*x^16 + 2*x^15 - x^14
    - 3*x^13 - 5*x^12 - 5*x^11 - 5*x^10 - 5*x^9 - 3*x^8 - x^7 + 2*x^6
    + 4*x^5 + 5*x^4 + 7*x^3 + 6*x^2 + 5*x + 3]

NOTES:

"""

attach("prescribed_roots_pyx.spyx")

## Auxiliary function for detecting roots of unity
def no_roots_of_unity(pol):
    """
    Return True if the given polynomial is irreducible and not cyclotomic,
    and False if it is divisible by a cyclotomic polynomial. Other inputs
    give undetermined results.

    INPUT:
        pol -- polynomial with rational coefficients

    OUTPUT:
        True if pol is irreducible and not cyclotomic.
        False if pol has a cyclotomic factor.

    EXAMPLES:
        sage: pol.<x> = PolynomialRing(Rationals())
        sage: no_roots_of_unity(x^5-1)
        False
        sage: no_roots_of_unity(x^5-2)
        True
    """
    polRing = pol.parent()
    x = polRing.gen()

    pol1 = pol
    pol2 = pol(-x)
    while pol2 == pol1: # Force pol1 not to be even
        l = pol1.list()
        l1 = [l[i] for i in range(0, len(l), 2)]
        pol1 = polRing(l1)
        pol2 = pol1(-x)
    if not pol1.gcd(pol2).is_constant(): return(False) # zeta_{*}, v_2(*) >= 2
    l = (pol1*pol2).list()
    l1 = [l[i] for i in range(0, len(l), 2)]
    pol3 = polRing(l1)
    if not pol1.gcd(pol3).is_constant(): return(False) # zeta_{*}, v_2(*) = 0
    pol4 = pol3(-x)
    if not pol1.gcd(pol4).is_constant(): return(False) # zeta_{*}, v_2(*) = 1
    return(True)

def eh_test(pol, q):
    polRing = pol.parent()
    x = polRing.gen()

    pol1 = pol
    while pol1(1) == 0:
        pol1 //= x-1
    d = pol1.degree()
    v = pol1(-1)
    if d%2:
        v *= q
    return(v.is_square())

def asymmetrize(P):
    """
    Convert a self-inversive polynomial to asymmetric form.

    INPUT:
        P -- polynomial with rational coefficients, which must be
        self-inversive

    OUTPUT:
        Q, R -- two polynomials such that R is a monic divisor of (x+1)(x-1),
        and P(x) = Q(x + 1/x) x^{\deg(Q)} R(x).

    EXAMPLES:
        sage: pol.<x> = PolynomialRing(Rationals())
        sage: asymmetrize(x^5-1)
        (x^2 + x - 1, x - 1)
        
    """
    polRing = P.parent()
    x = polRing.gen()
    d = P.degree()
    P0 = P.reverse()
    if (P == P0):
        sign = +1
    elif (P == -P0):
        sign = -1
    else:
        raise RuntimeError, "Polynomial not self-inversive"
    if P.degree() % 2 == 1:
        if sign == 1:
            cofactor = x+1
        else:
            cofactor = x-1
    elif sign == 1:
        cofactor = polRing(1)
    else:
        cofactor = x^2 - 1
    Q = P // cofactor
    coeffs = []
    m = Q.degree() // 2
    for i in reversed(range(m+1)):
        coeffs.insert(0, Q.constant_coefficient())
        Q = (Q % (x^2+1)^i) // x
    return polRing(coeffs), cofactor

def symmetrize(Q, R=1):
    """
    Convert a self-inversive polynomial from asymmetric form.

    INPUT:
        Q, R -- polynomials with rational coefficients

    OUTPUT:
        P -- the polynomial P(x) = Q(x + 1/x) x^{\deg(Q)} R(x).

    EXAMPLES: 
        sage: pol.<x> = PolynomialRing(Rationals())
        sage: symmetrize(x^2+x-1, x-1)
        x^5 - 1
    
    """
    polRing = Q.parent()
    x = polRing.gen()
    return polRing(x^(Q.degree()) * Q(x + 1/x)) * R

def roots_on_unit_circle(P0, modulus=1, n=1,
                         answer_count=None,
                         verbosity=None, node_count=None, filter=None):
    """
    Find polynomials with roots on the unit circle under extra restrictions.

    INPUT:
        P0 -- polynomial with rational coefficients, which must be
           self-inversive
        m -- positive integer or list of positive integers
        n -- positive integer
        answer_count -- positive integer or None
        node_count -- positive integer or None; if not None, an exception will
            be raised if this many nodes of the tree are encountered.
	filter -- function or None; if not None, only polynomials for which 
            this function evaluates to True will be returned.

    OUTPUT:
        list -- a list of all polynomials P with roots on the unit circle
            such that P is congruent to P0 modulo m and shares its highest
            n coefficients with P0. If answer_count is not None, return at
            most answer_count polynomials, otherwise return all of them.
	integer -- the number of terminal nodes in the tree enumerated
            in order to compute the list.

    EXAMPLES:
        sage: pol.<x> = PolynomialRing(Rationals())
        sage: roots_on_unit_circle(x^5 - 1, 2, 1)
        ([x^5 - 1, x^5 - 2*x^4 + 2*x^3 - 2*x^2 + 2*x - 1], 4)
        sage: roots_on_unit_circle(x^5 - 1, 4, 1)
        ([x^5 - 1], 2)
        
    """
    polRing = P0.parent()
    x = polRing.gen()

    Q0, cofactor = asymmetrize(P0)
    sign = cmp(Q0.leading_coefficient(), 0)
    Q0 *= sign
    d = Q0.degree()
    lead = Q0.leading_coefficient()

    count = 0

    try:
        modlist = list(modulus)
    except TypeError:
        modlist = [modulus]
    modlist = [0]*n + modlist
    modlist += [modlist[-1]] * (d+1 - len(modlist)) + [1]

    process = process_queue(d, n, lead, modlist, node_count,
                      verbosity, Q0)
    ans = []
    anslen = 0
    try:
        while True:
            t = process.exhaust_next_answer()
            if t>0:
                Q2 = polRing(process.Q1_array.tolist())
                Q3 = sign * symmetrize(Q2, cofactor)
                if filter == None or filter(Q3):
                    ans.append(Q3)
                    anslen += 1
                    if answer_count != None and anslen >= answer_count:
                        break
            elif t==0:
                break
            else:
                raise RuntimeError, "Node count (" + str(self.node_count) + ") exceeded"
    finally:
        process.clear()
    return(ans, process.count)
