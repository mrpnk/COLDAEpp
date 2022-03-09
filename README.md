# Coldae++
Coldae++ is a C++ port of NetLib's COLDAE.

COLDAE uses a collocation (COL) method to solve semi-explicit differential-algebraic equations (DAE).

It supports boundary values of type





The original Fortran 77 code by Uri Ascher and Ray Spiteri is available on [Netlib](http://www.netlib.org/ode/coldae.f).



Initialize the {fmt} submodule :

`git submodule update --init --recursive`

### Benchmark

The C++ version takes (117+-2)% of time and about 10 times the memory of the original version (without fast-math).

60% with fast-math

