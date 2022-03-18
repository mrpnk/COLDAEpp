![test workflow](https://github.com/mrpnk/COLDAEpp/actions/workflows/test.yml/badge.svg)

# Coldae++
Coldae++ is a C++ port of NetLib's COLDAE.

COLDAE uses a collocation (COL) method to solve semi-explicit differential-algebraic equations (DAE) of the form

![equation_1](D:\Simon\Dokumente\Entwicklung\coldae++\equation_1.png)







The original Fortran 77 code by Uri Ascher and Ray Spiteri is available on [Netlib](http://www.netlib.org/ode/coldae.f).



Initialize the {fmt} submodule :

`git submodule update --init --recursive`



### Benchmark

The C++ version takes about 20% less time than the original version but much more memory.

