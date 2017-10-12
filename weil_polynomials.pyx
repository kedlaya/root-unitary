r"""
Iterator for Weil polynomials

AUTHOR:
  -- Kiran S. Kedlaya (2007-05-28): initial version
  -- (2015-08-29): switch from NTL to FLINT
  -- (2017-10-03): consolidate Sage layer into .pyx file
                   define WeilPolynomials iterator
                   reverse convention for polynomials
                   pass multiprecision integers to/from C
        
NOTES:
    This depends on `trac`:23947: for the trace polynomial construction.

TODO:
    - Implement real root isolation (e-antic)
"""

#clang c
#cinclude $SAGE_LOCAL/include/flint/
#clib gomp
#cargs -fopenmp
#cfile all_roots_in_interval.c
#cfile power_sums.c

from cpython cimport array
from cython.parallel import prange
from libc.stdlib cimport malloc, free
cimport cython

from sage.rings.integer cimport Integer
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.libs.gmp.types cimport *
from sage.libs.gmp.mpz cimport *
from sage.libs.flint.types cimport *
from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_vec cimport *

cdef extern from "power_sums.h":
    ctypedef struct ps_static_data_t:
        pass

    ctypedef struct ps_dynamic_data_t:
        int flag
        long node_count
        fmpz *sympol # Holds a return value

    ps_static_data_t *ps_static_init(int d, int q, int coeffsign, fmpz_t lead,
    		     		     int cofactor, 
                                     fmpz *modlist,
                                     long node_limit)
    ps_dynamic_data_t *ps_dynamic_init(int d, fmpz *coefflist)
    ps_dynamic_data_t *ps_dynamic_split(ps_dynamic_data_t *dy_data)
    void ps_static_clear(ps_static_data_t *st_data)
    void ps_dynamic_clear(ps_dynamic_data_t *dy_data)
    void next_pol(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data) nogil

# Data structure to manage parallel depth-first search.
cdef class dfs_manager:
    cdef public long count   
    cdef int d
    cdef int num_threads
    cdef long node_limit
    cdef ps_static_data_t *ps_st_data
    cdef ps_dynamic_data_t **dy_data_buf

    def __cinit__(self, int d, q, coefflist, modlist, int sign, int cofactor,
                        int node_limit, int num_threads):
        cdef int i
        cdef fmpz_t temp_lead
        cdef fmpz *temp_array
        self.count = 0
        self.d = d
        self.num_threads = num_threads
        self.dy_data_buf = <ps_dynamic_data_t **>malloc(num_threads*cython.sizeof(cython.pointer(ps_dynamic_data_t)))
        for i in range(self.num_threads):
            self.dy_data_buf[i] = NULL
        self.node_limit = node_limit
        fmpz_init(temp_lead)
        fmpz_set_mpz(temp_lead, Integer(coefflist[-1]).value)
        temp_array = _fmpz_vec_init(d+1)
        for i in range(d+1):
            fmpz_set_mpz(temp_array+i, Integer(modlist[i]).value)
        self.ps_st_data = ps_static_init(d, q, sign, temp_lead, cofactor,
                                    temp_array, node_limit)
        fmpz_clear(temp_lead)
        for i in range(d+1):
            fmpz_set_mpz(temp_array+i, Integer(coefflist[i]).value)
        self.dy_data_buf[0] = ps_dynamic_init(d, temp_array)
        _fmpz_vec_clear(temp_array, d+1)
        for i in range(1, self.num_threads):
            self.dy_data_buf[i] = NULL

    def __init__(self, d, q, coefflist, modlist, sign, cofactor,
                        node_limit, num_threads):
        pass

    def __dealloc__(self):
        ps_static_clear(self.ps_st_data)
        self.ps_st_data = NULL
        if self.dy_data_buf != NULL:
            for i in range(self.num_threads):
                ps_dynamic_clear(self.dy_data_buf[i])
            free(self.dy_data_buf)
            self.dy_data_buf = NULL

    cpdef object advance_exhaust(self):
        """
        Advance the tree exhaustion.
        """
        cdef int i, j, k, d = self.d, t=1, np = self.num_threads
        cdef long ans_count = 100*np
        cdef mpz_t z
        cdef Integer temp
        ans = []

        k=0
        while (t>0 and len(ans) < ans_count):
            t = 0
            if np>1: k = k%(np-1) + 1
            with nogil: # Drop GIL for this parallel loop
                for i in prange(np, schedule='dynamic', num_threads=np):
                    if self.dy_data_buf[i] != NULL:
                        next_pol(self.ps_st_data, self.dy_data_buf[i])
            for i in range(np):
                if self.dy_data_buf[i] != NULL:
                    if self.dy_data_buf[i].flag > 0:
                        t += 1
                        l = []
                        # Convert a vector of fmpz's into mpz's and then Integers.
                        for j in range(2*d+3):
                            flint_mpz_init_set_readonly(z, &self.dy_data_buf[i].sympol[j])
                            temp = Integer(0)
                            mpz_set(temp.value, z)
                            l.append(temp)
                            flint_mpz_clear_readonly(z)
                        ans.append(l)
                    elif self.dy_data_buf[i].flag < 0:
                        raise RuntimeError("Node limit (" + str(self.node_limit) + ") exceeded")
                    else:
                        self.count += self.dy_data_buf[i].node_count
                        ps_dynamic_clear(self.dy_data_buf[i])
                        self.dy_data_buf[i] = NULL
                if self.dy_data_buf[i] == NULL: # Steal work
                    j = (i-k) % np
                    self.dy_data_buf[i] = ps_dynamic_split(self.dy_data_buf[j])
        return ans

class WeilPolynomials():
    r"""
    Iterator for Weil polynomials, i.e., integer polynomials with all complex roots
    having a particular absolute value.

    If num_threads is specified, a parallel search (using the specified
    number of threads) is carried out. In this case, the order of values is not
    specified.

    EXAMPLES:

    By Kronecker's theorem, a monic integer polynomial has all roots of absolute
    value 1 if and only if it is a product of cyclotomic polynomials. For such a
    product to have positive sign of the functional equation, the factors `x-1`
    and `x+1` must each occur with even multiplicity. This code confirms Kronecker's
    theorem for polynomials of degree 6::
        sage: P.<x> = PolynomialRing(ZZ)
        sage: d = 6
        sage: ans1 = list(WeilPolynomials(d, 1, 1))
        sage: ans1.sort()
        sage: l = [(x-1)^2, (x+1)^2] + [cyclotomic_polynomial(n,x) 
        ....:     for n in range(3, 2*d*d) if euler_phi(n) <= d]
        sage: w = WeightedIntegerVectors(d, [i.degree() for i in l])
        sage: ans2 = [prod(l[i]^v[i] for i in range(len(l))) for v in w]
        sage: ans2.sort()
        sage: print(ans1 == ans2)
        True
    """
    def __init__(self, d, q, sign=1, lead=1, node_limit=None, num_threads=1):
        r"""
        Initialize an iterator for Weil polynomials.

        INPUT:
        
            -- ``d`` - degree of the Weil polynomials
            -- ``q`` - square absolute value of the roots
            -- ``sign`` - sign of the functional equation (default: 1)
            -- ``lead`` - one or more leading coefficients; see below (default: 1)
            -- ``node_limit`` -- if specified, maximum number of nodes to allow in
                   a single tree traversal without raising a RuntimeError
            -- ``num_threads`` -- number of threads to use for parallel computation
                   (default: 1)

        If ``lead`` cannot be parsed as a list, it is treated as a single integer
        which will be the leading coefficient of all returned polynomials. Otherwise,
        each entry is parsed either as a single integer, which is a prescribed coefficient,
        or a pair `(i,j)` in which `i` is a coefficient which is prescribed modulo `j`.
        Unless `d` is even and `sign` is 1, the values of `j` are required to be monotone
        under reverse divisibility (that is, the first value must be a multiple of the second
        and so on, with omitted values taken to be 0).

        """
        self.pol = PolynomialRing(QQ, name='x')
        x = self.pol.gen()
        d = Integer(d)
        if sign != 1 and sign != -1:
            return ValueError("Invalid sign")
        if d%2==0:
            if sign==1:
                d2 = d//2
                num_cofactor = 0
            else:
                d2 = d//2 - 1
                num_cofactor = 3
        else:
            d2 = d//2
            if sign==1: num_cofactor = 1
            else: num_cofactor = 2
        try:
            leadlist = list(lead)
        except TypeError:
            leadlist = [(lead, 0)]
        coefflist = []
        modlist = []
        for i in leadlist:
            try:
                (j,k) = i
            except TypeError:
                (j,k) = (i,0)
            j = Integer(j)
            k = Integer(k)
            if len(modlist) == 0 and k != 0:
                raise ValueError("Leading coefficient must be specified exactly")
            if num_cofactor != 0 and len(modlist) > 0 and modlist[-1]%k != 0:
                raise ValueError("Invalid moduli")
            coefflist.append(j)
            modlist.append(k)
        for _ in range(d2+1-len(coefflist)):
            coefflist.append(0)
            modlist.append(1)
        coeffsign = cmp(coefflist[0], 0)
        coefflist = map(lambda x: x*coeffsign, coefflist)
        coefflist.reverse()
        self.process = dfs_manager(d2, q, coefflist, modlist, coeffsign, num_cofactor, 
                                     (-1 if node_limit is None else node_limit), num_threads)
        self.ans = []

    def __iter__(self):
        return(self)

    def next(self):
        if len(self.ans) == 0:
            if self.process is None:
                raise StopIteration
            self.ans = self.process.advance_exhaust()
        if len(self.ans) == 0:
            self.count = self.process.count
            self.process = None
            raise StopIteration
        return self.pol(self.ans.pop())

    def node_count(self):
        if self.process is None:
            return self.count
        return self.process.count

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
#    num_cofactor = [1, x+q.sqrt(), x-q.sqrt(), x**2-q].index(cofactor)
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
    return(ans, temp.count)
