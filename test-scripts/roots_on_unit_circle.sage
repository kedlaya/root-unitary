# This requires the WeilPolynomials package. It also depends on
# trac ticket 23947 for the trace polynomial.

def roots_on_unit_circle(P0, modulus=1, n=1,
                         answer_count=None, node_limit=None, filter=None,
                         num_threads=1, return_nodes=False):
    """
    Find polynomials with roots on the unit circle under extra restrictions.

    INPUT:
        P0 -- polynomial with rational coefficients, which must be
           self-inversive
        m -- positive integer or list of positive integers
        n -- positive integer
        answer_count -- positive integer or None
        node_limit -- positive integer or None; if not None, an exception will
            be raised if this many nodes of the tree are encountered.
	filter -- function or None; if not None, only polynomials for which 
            this function evaluates to True will be returned.
        return_nodes -- boolean (default False)

    OUTPUT:
        list -- a list of all polynomials P with roots on the unit circle
            such that P is congruent to P0 modulo m and shares its highest
            n coefficients with P0. If answer_count is not None, return at
            most answer_count polynomials, otherwise return all of them.

        If return_nodes is True, there is a second return value giving
	the number of terminal nodes in the tree enumerated
            in order to compute the list.

    EXAMPLES::
        sage: polRing.<x> = PolynomialRing(Rationals())
        sage: P0 = polRing([3,5,6,7,5,4,2,-1,-3,-5,-5,-5,-5,-3,-1,2,4,5,7,6,5,3])
        sage: roots_on_unit_circle(P0, 3^2, 1)
        [3*x^21 + 5*x^20 + 6*x^19 + 7*x^18 + 5*x^17 + 4*x^16 + 2*x^15 - x^14
        - 3*x^13 - 5*x^12 - 5*x^11 - 5*x^10 - 5*x^9 - 3*x^8 - x^7 + 2*x^6
        + 4*x^5 + 5*x^4 + 7*x^3 + 6*x^2 + 5*x + 3]

        sage: pol.<x> = PolynomialRing(Rationals())
        sage: roots_on_unit_circle(x^5 - 1, 2, 1, return_nodes=True)
        ([x^5 - 1, x^5 - 2*x^4 + 2*x^3 - 2*x^2 + 2*x - 1], 4)
        sage: roots_on_unit_circle(x^5 - 1, 4, 1)
        [x^5 - 1]
        
    """
    pol = P0.parent()
    x = pol.gen()
    Q0, cofactor, q = P0.trace_polynomial()
    d = Q0.degree()
    coefflist = list(Q0)
    coefflist.reverse()
    modlist = [0 for _ in range(n)] + [modulus for _ in range(d+1-n)]
    lead = [(coefflist[i], modlist[i]) for i in range(d+1)]

    temp = WeilPolynomials(P0.degree(), q, cofactor.constant_coefficient().sign(),
                           lead, node_limit, num_threads)
    ans = []
    anslen = 0
    for i in temp:
        if filter == None or filter(i):
            ans.append(i)
            anslen += 1
            if answer_count != None and anslen >= answer_count:
                break    
    if not return_nodes:
        return ans
    return(ans, temp.node_count())
