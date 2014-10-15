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

#include <stdafx.h>
#include "laguerre.h"

/*************************************************************************
Calculation of the value of the Laguerre polynomial.

Parameters:
    n   -   degree, n>=0
    x   -   argument

Result:
    the value of the Laguerre polynomial Ln at x
*************************************************************************/
double laguerrecalculate(const int& n, const double& x)
{
    double result;
    double a;
    double b;
    double i;

    result = 1;
    a = 1;
    b = 1-x;
    if( n==1 )
    {
        result = b;
    }
    i = 2;
    while(ap::fp_less_eq(i,n))
    {
        result = ((2*i-1-x)*b-(i-1)*a)/i;
        a = b;
        b = result;
        i = i+1;
    }
    return result;
}


/*************************************************************************
Summation of Laguerre polynomials using Clenshaw’s recurrence formula.

This routine calculates c[0]*L0(x) + c[1]*L1(x) + ... + c[N]*LN(x)

Parameters:
    n   -   degree, n>=0
    x   -   argument

Result:
    the value of the Laguerre polynomial at x
*************************************************************************/
double laguerresum(const ap::real_1d_array& c, const int& n, const double& x)
{
    double result;
    double b1;
    double b2;
    int i;

    b1 = 0;
    b2 = 0;
    for(i = n; i >= 0; i--)
    {
        result = (2*i+1-x)*b1/(i+1)-(i+1)*b2/(i+2)+c(i);
        b2 = b1;
        b1 = result;
    }
    return result;
}


/*************************************************************************
Representation of Ln as C[0] + C[1]*X + ... + C[N]*X^N

Input parameters:
    N   -   polynomial degree, n>=0

Output parameters:
    C   -   coefficients
*************************************************************************/
void laguerrecoefficients(const int& n, ap::real_1d_array& c)
{
    int i;

    c.setbounds(0, n);
    c(0) = 1;
    for(i = 0; i <= n-1; i++)
    {
        c(i+1) = -c(i)*(n-i)/(i+1)/(i+1);
    }
}




