/*************************************************************************
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

#ifndef _chebyshev_h
#define _chebyshev_h

#include "ap.h"
#include "ialglib.h"

/*************************************************************************
Calculation of the value of the Chebyshev polynomials of the
first and second kinds.

Parameters:
    r   -   polynomial kind, either 1 or 2.
    n   -   degree, n>=0
    x   -   argument, -1 <= x <= 1

Result:
    the value of the Chebyshev polynomial at x
*************************************************************************/
double chebyshevcalculate(const int& r, const int& n, const double& x);


/*************************************************************************
Summation of Chebyshev polynomials using Clenshaw’s recurrence formula.

This routine calculates
    c[0]*T0(x) + c[1]*T1(x) + ... + c[N]*TN(x)
or
    c[0]*U0(x) + c[1]*U1(x) + ... + c[N]*UN(x)
depending on the R.

Parameters:
    r   -   polynomial kind, either 1 or 2.
    n   -   degree, n>=0
    x   -   argument

Result:
    the value of the Chebyshev polynomial at x
*************************************************************************/
double chebyshevsum(const ap::real_1d_array& c,
     const int& r,
     const int& n,
     const double& x);


/*************************************************************************
Representation of Tn as C[0] + C[1]*X + ... + C[N]*X^N

Input parameters:
    N   -   polynomial degree, n>=0

Output parameters:
    C   -   coefficients
*************************************************************************/
void chebyshevcoefficients(const int& n, ap::real_1d_array& c);


/*************************************************************************
Conversion of a series of Chebyshev polynomials to a power series.

Represents A[0]*T0(x) + A[1]*T1(x) + ... + A[N]*Tn(x) as
B[0] + B[1]*X + ... + B[N]*X^N.

Input parameters:
    A   -   Chebyshev series coefficients
    N   -   degree, N>=0
    
Output parameters
    B   -   power series coefficients
*************************************************************************/
void fromchebyshev(const ap::real_1d_array& a,
     const int& n,
     ap::real_1d_array& b);


#endif

