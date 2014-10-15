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
#include "legendre.h"

/*************************************************************************
Calculation of the value of the Legendre polynomial Pn.

Parameters:
    n   -   degree, n>=0
    x   -   argument

Result:
    the value of the Legendre polynomial Pn at x
*************************************************************************/
double legendrecalculate(const int& n, const double& x)
{
    double result;
    double a;
    double b;
    int i;

    result = 1;
    a = 1;
    b = x;
    if( n==0 )
    {
        result = a;
        return result;
    }
    if( n==1 )
    {
        result = b;
        return result;
    }
    for(i = 2; i <= n; i++)
    {
        result = ((2*i-1)*x*b-(i-1)*a)/i;
        a = b;
        b = result;
    }
    return result;
}


/*************************************************************************
Summation of Legendre polynomials using Clenshaw’s recurrence formula.

This routine calculates
    c[0]*P0(x) + c[1]*P1(x) + ... + c[N]*PN(x)

Parameters:
    n   -   degree, n>=0
    x   -   argument

Result:
    the value of the Legendre polynomial at x
*************************************************************************/
double legendresum(const ap::real_1d_array& c, const int& n, const double& x)
{
    double result;
    double b1;
    double b2;
    int i;

    b1 = 0;
    b2 = 0;
    for(i = n; i >= 0; i--)
    {
        result = (2*i+1)*x*b1/(i+1)-(i+1)*b2/(i+2)+c(i);
        b2 = b1;
        b1 = result;
    }
    return result;
}


/*************************************************************************
Representation of Pn as C[0] + C[1]*X + ... + C[N]*X^N

Input parameters:
    N   -   polynomial degree, n>=0

Output parameters:
    C   -   coefficients
*************************************************************************/
void legendrecoefficients(const int& n, ap::real_1d_array& c)
{
    int i;

    c.setbounds(0, n);
    for(i = 0; i <= n; i++)
    {
        c(i) = 0;
    }
    c(n) = 1;
    for(i = 1; i <= n; i++)
    {
        c(n) = c(n)*(n+i)/2/i;
    }
    for(i = 0; i <= n/2-1; i++)
    {
        c(n-2*(i+1)) = -c(n-2*i)*(n-2*i)*(n-2*i-1)/2/(i+1)/(2*(n-i)-1);
    }
}




