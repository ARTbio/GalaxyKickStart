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
#include "hermite.h"

/*************************************************************************
Calculation of the value of the Hermite polynomial.

Parameters:
    n   -   degree, n>=0
    x   -   argument

Result:
    the value of the Hermite polynomial Hn at x
*************************************************************************/
double hermitecalculate(const int& n, const double& x)
{
    double result;
    int i;
    double a;
    double b;

    
    //
    // Prepare A and B
    //
    a = 1;
    b = 2*x;
    
    //
    // Special cases: N=0 or N=1
    //
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
    
    //
    // General case: N>=2
    //
    for(i = 2; i <= n; i++)
    {
        result = 2*x*b-2*(i-1)*a;
        a = b;
        b = result;
    }
    return result;
}


/*************************************************************************
Summation of Hermite polynomials using Clenshaw’s recurrence formula.

This routine calculates
    c[0]*H0(x) + c[1]*H1(x) + ... + c[N]*HN(x)

Parameters:
    n   -   degree, n>=0
    x   -   argument

Result:
    the value of the Hermite polynomial at x
*************************************************************************/
double hermitesum(const ap::real_1d_array& c, const int& n, const double& x)
{
    double result;
    double b1;
    double b2;
    int i;

    b1 = 0;
    b2 = 0;
    for(i = n; i >= 0; i--)
    {
        result = 2*(x*b1-(i+1)*b2)+c(i);
        b2 = b1;
        b1 = result;
    }
    return result;
}


/*************************************************************************
Representation of Hn as C[0] + C[1]*X + ... + C[N]*X^N

Input parameters:
    N   -   polynomial degree, n>=0

Output parameters:
    C   -   coefficients
*************************************************************************/
void hermitecoefficients(const int& n, ap::real_1d_array& c)
{
    int i;

    c.setbounds(0, n);
    for(i = 0; i <= n; i++)
    {
        c(i) = 0;
    }
    c(n) = exp(n*log(double(2)));
    for(i = 0; i <= n/2-1; i++)
    {
        c(n-2*(i+1)) = -c(n-2*i)*(n-2*i)*(n-2*i-1)/4/(i+1);
    }
}




