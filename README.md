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

For the version circa 2015, see the branch archive-2015.

As of December 2019, the code is slated for inclusion into Sage: try
```
sage: R.<x> = ZZ[]
sage: R.weil_polynomials?
```
to learn more. I plan to maintain this repository with bugfixes to the
underlying C code, as well as auxiliary code that is not intended for 
inclusion in Sage. That includes scripts from the aforementioned papers 
as well as code for abelian varieties in the L-Functions and Modular 
Forms Database (https://www.lmfdb.org).

The source files for the main package are:
* weil_polynomials.pyx: Cython source file for Sage, providing the
   class WeilPolynomials to iterate over Weil polynomials
* power_sums.c: C code (using FLINT) to enumerate the tree based on
    Rolle's theorem, Sturm's theorem, power sum bounds, etc.
* power_sums.h: associated header file

To use the Cython file, from a Sage prompt type
```
  sage: load("weil_polynomials.pyx")
```
and everything should compile automatically.

If the OpenMP library is available systemwide, you can also do some computations
in parallel. To enable this, you will need to uncomment a couple of lines in
weil_polynomials.pyx, and possibly change the value of the environment variable
OMP_NUM_THREADS (thanks to Edgar Costa for this tip).

There are also some test scripts. See the README files in the following directories:
* test-scripts: Miscellaneous tests
* av-scripts: Build tables associated to abelian varieties
* k3-scripts: Build tables associated to K3 surfaces
* k3-quartic-f2: Scripts and data associated to smooth quartic surfaces over F2

POSSIBLE TODO LIST: 
* Improve scheduling of parallel computation in the current Cython model.
* Add a Julia wrapper for use in Nemo.
* Use real root isolation instead of (or in addition to) Sturm's theorem. One
   possible implementation, again using FLINT, is available in e-antic.
* Use the mts library for budgeted reverse search:
  http://cgm.cs.mcgill.ca/~avis/doc/tutorial.html

