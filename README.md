![test workflow](https://github.com/mrpnk/COLDAEpp/actions/workflows/test.yml/badge.svg)

# Coldae++
Coldae++ is a C++ port of NetLib's COLDAE.

COLDAE uses a collocation (COL) method to solve semi-explicit differential-algebraic equations (DAE) of the form

<p align="middle">
  <img width=450 src=/doc/equation_1.jpg>
</p>





The original Fortran 77 code by Uri Ascher and Ray Spiteri is available on [Netlib](http://www.netlib.org/ode/coldae.f).



### How to use it

- Initialize the {fmt} submodule :
  `git submodule update --init --recursive`
- `#include "ColDAEpp.hpp"`
- Derive a class from `coldae::system` and provide the subroutines. (see `test\systems.hpp`)
- Create a `coldae::solver` instance and call `COLDAE()` on it. The parameters  have the same meaning as in the original. (see `test\test.hpp`)

 

### Benchmark

The C++ version takes about 20% less time than the original version but much more memory.

