"""
Polynomial filter functions
"""
def ej_test(pol): # Elsenhans-Jahnel condition based on Artin-Tate formula
    """
    Test a polynomial for the Elsenhans-Jahnel condition.

    EXAMPLES::

    """
    polRing = pol.parent()
    x = polRing.gen()

    pol1 = pol[0].sign() * pol // (1-x)^pol.ord(1-x)
    return(pol1(-1).is_square())

def has_cyclotomic_factor(pol, assume_irreducible=False):
    """
    Return True if the given polynomial has a nontrivial cyclotomic factor.

    If the polynomial is specified to be irreducible, a more efficient test is used.

    EXAMPLES::
        sage: pol.<x> = PolynomialRing(Rationals())
        sage: has_cyclotomic_factor(x^5-1)
        True
        sage: has_cyclotomic_factor(x^5-2)
        False

    .. SEEALSO::
        :meth:`cyclotomic_part`
    """
    polRing = pol.parent()
    x = polRing.gen()

    if assume_irreducible:
        pol1 = pol
        pol2 = pol(-x)
        while pol2 == pol1: # Force pol1 not to be even
            pol1 = polRing(list(pol1)[::2])
            pol2 = pol1(-x)
        if not pol1.gcd(pol2).is_constant(): return(True) # zeta_{*}, v_2(*) >= 2
        pol3 = polRing(list(pol1*pol2)[::2])
        if not pol1.gcd(pol3).is_constant(): return(True) # zeta_{*}, v_2(*) = 0
        pol4 = pol3(-x)
        if not pol1.gcd(pol4).is_constant(): return(True) # zeta_{*}, v_2(*) = 1
        return(False)
    pol1 = pol
    pol2 = pol1.gcd(pol1(-x))
    while not pol2.is_constant():
        pol1 = (pol1 // pol2) * polRing(list(pol2)[::2])
        pol2 = pol1.gcd(pol1(-x))
    pol1 = polRing(list(pol1*pol1(-x))[::2])
    while True:
        if pol1.is_constant(): return(False)
        pol2 = pol1.gcd(polRing(list(pol1*pol1(-x))[::2]))
        if pol1.degree() == pol2.degree(): return(True)
        pol1 = pol2

def ej_test(pol): # Elsenhans-Jahnel condition based on Artin-Tate formula
    """
    Test a polynomial for the Elsenhans-Jahnel condition.

    EXAMPLES::

    """
    polRing = pol.parent()
    x = polRing.gen()

    pol1 = pol[0].sign() * pol // (1-x)^pol.ord(1-x)
    return(pol1(-1).is_square())

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
    if P[0] == 0:
        raise ValueError, "Polynomial " + str(P) + " not self-inversive"
    d = P.degree()
    sg = (P[d]/P[0]).sign()
    q = abs(P[d]/P[0])^(2/d)
    if not q in Rationals():
        raise ValueError, "Polynomial " + str(P) + " not self-inversive"
    for i in range(d+1):
        if P[i] != sg*P[d-i]/q^(d/2-i):
            raise ValueError, "Polynomial " + str(P) + " not self-inversive"
    cofactor = polRing(1)
    Q = P
    if sg == -1:
        cofactor *= 1-q*x
        Q //= 1-q*x
    if Q.degree() %2 == 1:
        cofactor *= 1+q*x
        Q //= 1+q*x
    coeffs = []
    m = Q.degree() // 2
    for i in reversed(range(m+1)):
        coeffs.insert(0, Q.constant_coefficient())
        Q = (Q % (1 + q*x^2)^i) // x
    return polRing(coeffs), cofactor, q

def symmetrize(Q, R=1, q=1):
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
    return polRing(x^(Q.degree()) * Q(q*x + 1/x)) * R


