/*************************************************************************
Cephes Math Library Release 2.8:  June, 2000
Copyright by Stephen L. Moshier

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from C to
      pseudocode.

See subroutines comments for additional copyrights.

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************/

#ifndef _bessel_h
#define _bessel_h

#include "ap.h"
#include "ialglib.h"

/*************************************************************************
Bessel function of order zero

Returns Bessel function of order zero of the argument.

The domain is divided into the intervals [0, 5] and
(5, infinity). In the first interval the following rational
approximation is used:


       2         2
(w - r  ) (w - r  ) P (w) / Q (w)
      1         2    3       8

           2
where w = x  and the two r's are zeros of the function.

In the second interval, the Hankel asymptotic expansion
is employed with two rational functions of degree 6/6
and 7/7.

ACCURACY:

                     Absolute error:
arithmetic   domain     # trials      peak         rms
   IEEE      0, 30       60000       4.2e-16     1.1e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*************************************************************************/
double besselj0(double x);


/*************************************************************************
Bessel function of order one

Returns Bessel function of order one of the argument.

The domain is divided into the intervals [0, 8] and
(8, infinity). In the first interval a 24 term Chebyshev
expansion is used. In the second, the asymptotic
trigonometric representation is employed using two
rational functions of degree 5/5.

ACCURACY:

                     Absolute error:
arithmetic   domain      # trials      peak         rms
   IEEE      0, 30       30000       2.6e-16     1.1e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*************************************************************************/
double besselj1(double x);


/*************************************************************************
Bessel function of integer order

Returns Bessel function of order n, where n is a
(possibly negative) integer.

The ratio of jn(x) to j0(x) is computed by backward
recurrence.  First the ratio jn/jn-1 is found by a
continued fraction expansion.  Then the recurrence
relating successive orders is applied until j0 or j1 is
reached.

If n = 0 or 1 the routine for j0 or j1 is called
directly.

ACCURACY:

                     Absolute error:
arithmetic   range      # trials      peak         rms
   IEEE      0, 30        5000       4.4e-16     7.9e-17


Not suitable for large n or x. Use jv() (fractional order) instead.

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double besseljn(int n, double x);


/*************************************************************************
Bessel function of the second kind, order zero

Returns Bessel function of the second kind, of order
zero, of the argument.

The domain is divided into the intervals [0, 5] and
(5, infinity). In the first interval a rational approximation
R(x) is employed to compute
  y0(x)  = R(x)  +   2 * log(x) * j0(x) / PI.
Thus a call to j0() is required.

In the second interval, the Hankel asymptotic expansion
is employed with two rational functions of degree 6/6
and 7/7.



ACCURACY:

 Absolute error, when y0(x) < 1; else relative error:

arithmetic   domain     # trials      peak         rms
   IEEE      0, 30       30000       1.3e-15     1.6e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*************************************************************************/
double bessely0(double x);


/*************************************************************************
Bessel function of second kind of order one

Returns Bessel function of the second kind of order one
of the argument.

The domain is divided into the intervals [0, 8] and
(8, infinity). In the first interval a 25 term Chebyshev
expansion is used, and a call to j1() is required.
In the second, the asymptotic trigonometric representation
is employed using two rational functions of degree 5/5.

ACCURACY:

                     Absolute error:
arithmetic   domain      # trials      peak         rms
   IEEE      0, 30       30000       1.0e-15     1.3e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*************************************************************************/
double bessely1(double x);


/*************************************************************************
Bessel function of second kind of integer order

Returns Bessel function of order n, where n is a
(possibly negative) integer.

The function is evaluated by forward recurrence on
n, starting with values computed by the routines
y0() and y1().

If n = 0 or 1 the routine for y0 or y1 is called
directly.

ACCURACY:
                     Absolute error, except relative
                     when y > 1:
arithmetic   domain     # trials      peak         rms
   IEEE      0, 30       30000       3.4e-15     4.3e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double besselyn(int n, double x);


/*************************************************************************
Modified Bessel function of order zero

Returns modified Bessel function of order zero of the
argument.

The function is defined as i0(x) = j0( ix ).

The range is partitioned into the two intervals [0,8] and
(8, infinity).  Chebyshev polynomial expansions are employed
in each interval.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,30        30000       5.8e-16     1.4e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double besseli0(double x);


/*************************************************************************
Modified Bessel function of order one

Returns modified Bessel function of order one of the
argument.

The function is defined as i1(x) = -i j1( ix ).

The range is partitioned into the two intervals [0,8] and
(8, infinity).  Chebyshev polynomial expansions are employed
in each interval.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0, 30       30000       1.9e-15     2.1e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double besseli1(double x);


/*************************************************************************
Modified Bessel function, second kind, order zero

Returns modified Bessel function of the second kind
of order zero of the argument.

The range is partitioned into the two intervals [0,8] and
(8, infinity).  Chebyshev polynomial expansions are employed
in each interval.

ACCURACY:

Tested at 2000 random points between 0 and 8.  Peak absolute
error (relative when K0 > 1) was 1.46e-14; rms, 4.26e-15.
                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0, 30       30000       1.2e-15     1.6e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double besselk0(double x);


/*************************************************************************
Modified Bessel function, second kind, order one

Computes the modified Bessel function of the second kind
of order one of the argument.

The range is partitioned into the two intervals [0,2] and
(2, infinity).  Chebyshev polynomial expansions are employed
in each interval.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0, 30       30000       1.2e-15     1.6e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double besselk1(double x);


/*************************************************************************
Modified Bessel function, second kind, integer order

Returns modified Bessel function of the second kind
of order n of the argument.

The range is partitioned into the two intervals [0,9.55] and
(9.55, infinity).  An ascending power series is used in the
low range, and an asymptotic expansion in the high range.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,30        90000       1.8e-8      3.0e-10

Error is high only near the crossover point x = 9.55
between the two expansions used.

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1988, 2000 by Stephen L. Moshier
*************************************************************************/
double besselkn(int nn, double x);


#endif

