r"""
Iterator for Weil polynomials

AUTHOR:
  -- Kiran S. Kedlaya (2007-05-28): initial version
  -- (2015-08-29): switch from NTL to FLINT
  -- (2017-10-03): consolidate Sage layer into .pyx file
                   define WeilPolynomials iterator
                   reverse convention for polynomials
        
EXAMPLES::
    sage: polRing.<x> = PolynomialRing(Rationals())
    sage: P0 = polRing([3,5,6,7,5,4,2,-1,-3,-5,-5,-5,-5,-3,-1,2,4,5,7,6,5,3])
    sage: roots_on_unit_circle(P0, 3^2, 1)
    [3*x^21 + 5*x^20 + 6*x^19 + 7*x^18 + 5*x^17 + 4*x^16 + 2*x^15 - x^14
    - 3*x^13 - 5*x^12 - 5*x^11 - 5*x^10 - 5*x^9 - 3*x^8 - x^7 + 2*x^6
    + 4*x^5 + 5*x^4 + 7*x^3 + 6*x^2 + 5*x + 3]

NOTES:
    This depends on `trac`:23947: for the trace polynomial construction.

TODO:
    - Pass coefficients as FLINT integers rather than C integers
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

from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

cdef extern from "power_sums.h":
    ctypedef struct ps_static_data_t:
        pass
    ctypedef struct ps_dynamic_data_t:
        long count

    ps_static_data_t *ps_static_init(int d, int lead, int sign, int q,
    		     		     int cofactor, 
                                     int *modlist,
                                     long node_limit)
    ps_dynamic_data_t *ps_dynamic_init(int d, int *Q0)
    ps_dynamic_data_t *ps_dynamic_clone(ps_dynamic_data_t *dy_data)
    ps_dynamic_data_t *ps_dynamic_split(ps_dynamic_data_t *dy_data)
    void extract_symmetrized_pol(int *Q, ps_dynamic_data_t *dy_data)
    void ps_static_clear(ps_static_data_t *st_data)
    void ps_dynamic_clear(ps_dynamic_data_t *dy_data)
    int next_pol(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data) nogil

# Data structure to manage parallel depth-first search.
cdef class dfs_manager:
    cdef public int d
    cdef public long count
    cdef public int k
    cdef public array.array Qsym_array

    cdef array.array Q0_array
    cdef int[:] Q0
    cdef int[:] Qsym
    cdef int[:] modlist
    cdef array.array modlist_array
    cdef int sign
    cdef int cofactor
    cdef int num_threads
    cdef long node_limit
    cdef ps_static_data_t *ps_st_data
    cdef ps_dynamic_data_t **dy_data_buf

    def __cinit__(self, int d, int n, int lead, int sign, int q, int cofactor,
                  int num_threads, modlist, node_limit, Q):
        cdef int i
        self.d = d
        self.k = d
        self.sign = sign
        self.cofactor = cofactor
        self.num_threads = num_threads
        self.Q0_array = array.array('i', [0,] * (d+1))
        self.Q0 = self.Q0_array
        self.Qsym_array = array.array('i', [0,] * (2*d+3))
        self.Qsym = self.Qsym_array
        self.modlist_array = array.array('i', [0,] * (d+1))
        self.modlist = self.modlist_array
        for i in range(d+1):
            self.modlist[i] = modlist[d-i]
            self.Q0[i] = Q[i]
        if node_limit == None:
            self.node_limit = -1
        else:
            self.node_limit = node_limit
        self.count = 0
        self.ps_st_data = ps_static_init(d, lead, sign, q, cofactor,
                                    self.modlist_array.data.as_ints,
                                         self.node_limit)
        self.dy_data_buf = <ps_dynamic_data_t **>malloc(num_threads*cython.sizeof(cython.pointer(ps_dynamic_data_t)))
        self.dy_data_buf[0] = ps_dynamic_init(d, self.Q0_array.data.as_ints)
        for i in range(1, self.num_threads):
            self.dy_data_buf[i] = NULL

    def __init__(self, int d, int n, int lead, int sign, int q, int cofactor,
                 int num_threads, modlist, node_limit, Q):
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
        cdef int i, j, k, m, d = self.d, t=1, np = self.num_threads
        cdef long ans_count = 100*np
        ans = []
        cdef int *res = <int *>malloc(np*sizeof(int)) # Must be Python-free to avoid GIL

        k=0
        while (t>0 and len(ans) < ans_count):
            t = 0
            if np>1: k = k%(np-1) + 1
            with nogil: # Drop GIL for this parallel loop
                for i in prange(np, schedule='dynamic', num_threads=np):
                    if self.dy_data_buf[i] != NULL:
                        res[i] = next_pol(self.ps_st_data, self.dy_data_buf[i])
            for i in range(np):
                if self.dy_data_buf[i] != NULL:
                    if res[i] > 0:
                        t += 1
                        extract_symmetrized_pol(self.Qsym_array.data.as_ints,
                                                self.dy_data_buf[i])
                        ans.append(list(self.Qsym_array))
                    elif res[i] < 0:
                        free(res)
                        raise RuntimeError("Node limit (" + str(self.node_limit) + ") exceeded")
                    else: 
                        self.count += self.dy_data_buf[i].count
                        ps_dynamic_clear(self.dy_data_buf[i])
                        self.dy_data_buf[i] = NULL
                if self.dy_data_buf[i] == NULL: # Steal work
                    j = (i-k) % np
                    self.dy_data_buf[i] = ps_dynamic_split(self.dy_data_buf[j])
        free(res)
        return ans

class WeilPolynomials():
    r"""
    Iterator for Weil polynomials, i.e., integer polynomials with all complex roots
    having a particular absolute value.

    By default, Weil polynomials are computed using a serial depth-first search.
    Alternatively, if num_threads is specified, a parallel search (using the specified
    number of threads) is carried out at creation time.

    EXAMPLES:

    """
    def __init__(self, d, q, sign=1, lead=1, P0=None, modlist=None, node_limit=None,
                 num_threads=1):
        r"""
        Initialize an iterator for Weil polynomials.
        """
        if P0 is None:
            self.pol = PolynomialRing(QQ, name='x')
            x = self.pol.gen()
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
            n = 1
            coeffsign = 1
            Q0 = x**d2
            if modlist is None:
                modlist = []
            modlist = [0] + [1 for _ in range(d2+1-len(modlist))]
        else:
            self.pol = P0.parent()
            x = self.pol.gen()
            Q0, cofactor, q = P0.trace_polynomial()
            num_cofactor = [1, x+q.sqrt(), x-q.sqrt(), x**2-q].index(cofactor)
            coeffsign = cmp(Q0.leading_coefficient(), 0)
            Q0 *= coeffsign
            lead = Q0.leading_coefficient()
            d2 = Q0.degree()
            modlist += [modlist[-1] for _ in range(d2+1 - len(modlist))]
        self.process = dfs_manager(d2, n, lead, coeffsign, q, num_cofactor, 
                                     num_threads, modlist, node_limit, Q0)
        self.ans = []

    def __iter__(self):
        return(self)

    def next(self):
        if len(self.ans) == 0:
            if self.process is None:
                raise StopIteration
            self.ans = self.process.advance_exhaust()
        if len(self.ans) == 0:
            self.process = None
            raise StopIteration
        return self.pol(self.ans.pop())

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
        sage: pol.<x> = PolynomialRing(Rationals())
        sage: roots_on_unit_circle(x^5 - 1, 2, 1, return_nodes=True)
        ([x^5 - 1, x^5 - 2*x^4 + 2*x^3 - 2*x^2 + 2*x - 1], 4)
        sage: roots_on_unit_circle(x^5 - 1, 4, 1)
        [x^5 - 1]
        
    """
    try:
        modlist = list(modulus)
    except TypeError:
        modlist = [modulus]
    modlist = [0 for _ in range(n)] + modlist

    temp = WeilPolynomials(None, None, None, None, P0, modlist, node_limit, num_threads)
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
