# root-unitary
Search code for root-unitary (Weil) polynomials

This package implements a search strategy for Weil polynomials
(integer polynomials with all roots on a fixed circle) described in:

> K.S. Kedlaya, Search techniques for root-unitary polynomials, in 
>    Computational Arithmetic Geometry, Contemporary Mathematics 463, 
>    American Mathematical Society, 2008, 71-82;

plus some additional improvements described in:

> K.S. Kedlaya and A.V. Sutherland, A census of zeta functions of
>    quartic K3 surfaces over F_2, LMS Journal of Computing 19
>    special issue A (2016), 1-11.

Many cosmetic changes were made in preparation for migrating this package into Sage.
For the version circa 2015, see the branch archive-2015.

The source files for the main package are:
* weil_polynomials.pyx: Cython source file for Sage, providing the
   class WeilPolynomials to iterate over Weil polynomials
* power_sums.c: C code (using FLINT) to enumerate the tree based on
    Rolle's theorem, Sturm's theorem, and bounds computed from power sums
* power_sums.h: associated header file

From a Sage prompt, type
```
  sage: load("weil_polynomials.pyx")
```
and everything should compile automatically.

There are also some test scripts. See the README files in the following directories:
* test-scripts: Miscellaneous tests
* av-scripts: Build tables associated to abelian varieties
* k3-scripts: Build tables associated to K3 surfaces
* k3-quartic-f2: Scripts and data associated to smooth quartic surfaces over F2

POSSIBLE TODO LIST: 
* Improve parallel computation in the current Cython model.
* Add a port from Sage (based on Python) to Nemo (based on Julia).
* Use real root isolation instead of (or in addition to) Sturm's theorem. One
   possible implementation, again using FLINT, is available in e-antic.

