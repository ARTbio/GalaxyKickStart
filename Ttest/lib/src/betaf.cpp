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

#include <stdafx.h>
#include "betaf.h"

/*************************************************************************
Beta function


                  -     -
                 | (a) | (b)
beta( a, b )  =  -----------.
                    -
                   | (a+b)

For large arguments the logarithm of the function is
evaluated using lgam(), then exponentiated.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE       0,30       30000       8.1e-14     1.1e-14

Cephes Math Library Release 2.0:  April, 1987
Copyright 1984, 1987 by Stephen L. Moshier
*************************************************************************/
double beta(double a, double b)
{
    double result;
    double y;
    double sg;
    double s;

    sg = 1;
    ap::ap_error::make_assertion(ap::fp_greater(a,0)||ap::fp_neq(a,ap::ifloor(a)), "Overflow in Beta");
    ap::ap_error::make_assertion(ap::fp_greater(b,0)||ap::fp_neq(b,ap::ifloor(b)), "Overflow in Beta");
    y = a+b;
    if( ap::fp_greater(fabs(y),171.624376956302725) )
    {
        y = lngamma(y, s);
        sg = sg*s;
        y = lngamma(b, s)-y;
        sg = sg*s;
        y = lngamma(a, s)+y;
        sg = sg*s;
        ap::ap_error::make_assertion(ap::fp_less_eq(y,log(ap::maxrealnumber)), "Overflow in Beta");
        result = sg*exp(y);
        return result;
    }
    y = gamma(y);
    ap::ap_error::make_assertion(ap::fp_neq(y,0), "Overflow in Beta");
    if( ap::fp_greater(a,b) )
    {
        y = gamma(a)/y;
        y = y*gamma(b);
    }
    else
    {
        y = gamma(b)/y;
        y = y*gamma(a);
    }
    result = y;
    return result;
}




