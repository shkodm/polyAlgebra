polyAlg - Algebra of multivariative polynomials for R
=========

[![Build Status](https://travis-ci.org/llaniewski/polyAlgebra.svg?branch=master)](https://travis-ci.org/llaniewski/polyAlgebra)

Algebra of functions which looks like:

3 * x^2 * y^1 + 5 * z

And lot of other neat stuff.

Example:
```r
x = PV("x")
y = PV("y")
z = x^3+y+x*y
C(z) # prints z in C style
```
