"""
Finding polynomials with roots in prescribed regions

AUTHOR:
  -- Kiran S. Kedlaya (2007-05-28): initial version
  -- Kiran S. Kedlaya (2015-08-29): updated version; switch from NTL to FLINT
        
EXAMPLES:
    sage: polRing.<x> = PolynomialRing(Rationals())
    sage: P0 = 3*x^21 + 5*x^20 + 6*x^19 + 7*x^18 + 5*x^17 + 4*x^16 + 2*x^15 - 
    ....: x^14 - 3*x^13 - 5*x^12 - 5*x^11 - 5*x^10 - 5*x^9 - 3*x^8 - x^7 + 
    ....: 2*x^6 + 4*x^5 + 5*x^4 + 7*x^3 + 6*x^2 + 5*x + 3
    sage: ans, count = roots_on_unit_circle(P0, 3^2, 1)
    sage: print "Number of terminal nodes:", count
    Number of terminal nodes: 1404
    sage: print ans
    [3*x^21 + 5*x^20 + 6*x^19 + 7*x^18 + 5*x^17 + 4*x^16 + 2*x^15 - x^14
    - 3*x^13 - 5*x^12 - 5*x^11 - 5*x^10 - 5*x^9 - 3*x^8 - x^7 + 2*x^6
    + 4*x^5 + 5*x^4 + 7*x^3 + 6*x^2 + 5*x + 3]

NOTES:

"""

#clang c
#cinclude $SAGE_LOCAL/include/flint/
#clib gomp
#cargs -fopenmp
#cfile all_roots_in_interval.c
#cfile power_sums.c

from cpython cimport array
import array
from cython.parallel import prange
from libc.stdlib cimport malloc, free, rand
cimport cython

cdef extern from "power_sums.h":
    ctypedef struct ps_static_data_t:
        pass
    ctypedef struct ps_dynamic_data_t:
        long count

    ps_static_data_t *ps_static_init(int d, int lead, int sign, int q,
    		     		     int cofactor, 
                                     int *modlist,
                                     int verbosity, long node_count)
    ps_dynamic_data_t *ps_dynamic_init(int d, int *Q0)
    ps_dynamic_data_t *ps_dynamic_clone(ps_dynamic_data_t *dy_data)
    ps_dynamic_data_t *ps_dynamic_split(ps_dynamic_data_t *dy_data)
    void extract_pol(int *Q, ps_dynamic_data_t *dy_data)
    void extract_symmetrized_pol(int *Q, ps_dynamic_data_t *dy_data)
    void ps_static_clear(ps_static_data_t *st_data)
    void ps_dynamic_clear(ps_dynamic_data_t *dy_data)
    int next_pol(ps_static_data_t *st_data, ps_dynamic_data_t *dy_data) nogil

cdef class process_queue:
    cdef int d, verbosity
    cdef long node_count
    cdef public long count
    cdef public int k
    cdef public array.array Q0_array
    cdef int[:] Q0
    cdef public array.array Qsym_array
    cdef int[:] Qsym
    cdef int[:] modlist
    cdef public array.array modlist_array
    cdef int sign
    cdef int cofactor
    cdef ps_static_data_t *ps_st_data
    cdef ps_dynamic_data_t *ps_dy_data

    def __init__(self, int d, int n, int lead, int sign, int q, int cofactor,
                 modlist, node_count, verbosity, Q):
        cdef int i
        self.d = d
        self.k = d
        self.sign = sign
        self.cofactor = cofactor
        self.Q0_array = array.array('i', [0,] * (d+1))
        self.Q0 = self.Q0_array
        self.Qsym_array = array.array('i', [0,] * (2*d+3))
        self.Qsym = self.Qsym_array
        self.modlist_array = array.array('i', [0,] * (d+1))
        self.modlist = self.modlist_array
        for i in range(d+1):
            self.modlist[i] = modlist[d-i]
            self.Q0[i] = Q[i]
        if verbosity == None:
            self.verbosity = -1
        else:
            self.verbosity = verbosity
        if node_count == None:
            self.node_count = -1
        else:
            self.node_count = node_count
        self.count = 0
        self.ps_st_data = ps_static_init(d, lead, sign, q, cofactor,
                                    self.modlist_array.data.as_ints,
                                         self.verbosity, self.node_count)
        self.ps_dy_data = ps_dynamic_init(d, self.Q0_array.data.as_ints)

    def clear(self):
        ps_static_clear(self.ps_st_data)
        ps_dynamic_clear(self.ps_dy_data)

    cpdef int exhaust_next_answer(self):
        cdef int t
        t = next_pol(self.ps_st_data, self.ps_dy_data)
        extract_symmetrized_pol(self.Qsym_array.data.as_ints, self.ps_dy_data)
        self.count = self.ps_dy_data.count
        return(t)

    cpdef object parallel_exhaust(process_queue self, int num_processes, f=None):
        cdef ps_dynamic_data_t **dy_data_buf
        cdef int i, j, k, m, d = self.d, t=1, np = num_processes
        ans = []
        dy_data_buf = <ps_dynamic_data_t **>malloc(np*cython.sizeof(cython.pointer(ps_dynamic_data_t)))
        dy_data_buf[0] = ps_dynamic_clone(self.ps_dy_data)
        cdef int *res = <int *>malloc(np*sizeof(int))

        for i in range(1, np):
            dy_data_buf[i] = NULL
        k=0
        while (t>0):
            t = 0
            k = k%(np-1) + 1
            with nogil: # Drop GIL for this parallel loop
                for i in prange(np, schedule='dynamic', num_threads=np):
                    if dy_data_buf[i] != NULL:
                        res[i] = next_pol(self.ps_st_data, dy_data_buf[i])
            for i in range(np):
                if dy_data_buf[i] != NULL:
                    if res[i] > 0:
                        t += 1
                        extract_symmetrized_pol(self.Qsym_array.data.as_ints,
                                                dy_data_buf[i])
                        if (f != None):
                            f.write(str(list(self.Qsym_array)))
                            f.write("\n")
                        else: ans.append(list(self.Qsym_array))
                    else:
                        self.count += dy_data_buf[i].count
                        ps_dynamic_clear(dy_data_buf[i])
                        dy_data_buf[i] = NULL
                if dy_data_buf[i] == NULL:
                    j = (i-k) % np
                    dy_data_buf[i] = ps_dynamic_split(dy_data_buf[j])
        free(dy_data_buf)
        free(res)
        if (f != None): return None
        else: return(ans)

def roots_on_unit_circle(P0, modulus=1, n=1,
                         answer_count=None,
                         verbosity=None, node_count=None, filter=None,
                         num_threads=None, output=None):
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

    Q0, cofactor, q = P0.inverse_reciprocal_transform()
    num_cofactor = [1, x+sqrt(q), x-sqrt(q), x^2-q].index(cofactor)
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
    if len(modlist) < d+1:
        modlist += [modlist[-1]] * (d+1 - len(modlist))

    process = process_queue(d, n, lead, sign, q, num_cofactor, 
                            modlist, node_count, verbosity, Q0)
    ans = []
    anslen = 0
    if (num_threads): # parallel version
        ans1 = process.parallel_exhaust(num_threads, output)
        if output != None:
            return process.count
        for i in ans1:
            Q2 = polRing(i)
            if filter == None or filter(Q2):
                ans.append(Q2)
                anslen += 1
                if answer_count != None and anslen >= answer_count:
                    break
        process.clear()
        return(ans, process.count)

    try:
        while True:
            t = process.exhaust_next_answer()
            if t>0:
                Q2 = polRing(process.Qsym_array.tolist())
                if filter == None or filter(Q2):
                    if output != None: output.write(str(list(Q2)))
                    else: ans.append(Q2)
                    anslen += 1
                    if answer_count != None and anslen >= answer_count:
                        break
            elif t==0:
                break
            else:
                raise RuntimeError, "Node count (" + str(process.node_count) + ") exceeded"
    finally:
        process.clear()
    if output != None: return(process.count)
    return(ans, process.count)
