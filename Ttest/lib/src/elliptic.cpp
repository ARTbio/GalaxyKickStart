/*************************************************************************
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier

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
#include "elliptic.h"

/*************************************************************************
Complete elliptic integral of the first kind

Approximates the integral



           pi/2
            -
           | |
           |           dt
K(m)  =    |    ------------------
           |                   2
         | |    sqrt( 1 - m sin t )
          -
           0

using the approximation

    P(x)  -  log x Q(x).

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE       0,1        30000       2.5e-16     6.8e-17

Cephes Math Library, Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double ellipticintegralk(double m)
{
    double result;

    result = ellipticintegralkhighprecision(1.0-m);
    return result;
}


/*************************************************************************
Complete elliptic integral of the first kind

Approximates the integral



           pi/2
            -
           | |
           |           dt
K(m)  =    |    ------------------
           |                   2
         | |    sqrt( 1 - m sin t )
          -
           0

where m = 1 - m1, using the approximation

    P(x)  -  log x Q(x).

The argument m1 is used rather than m so that the logarithmic
singularity at m = 1 will be shifted to the origin; this
preserves maximum accuracy.

K(0) = pi/2.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE       0,1        30000       2.5e-16     6.8e-17

Алгоритм взят из библиотеки Cephes
*************************************************************************/
double ellipticintegralkhighprecision(double m1)
{
    double result;
    double p;
    double q;

    if( ap::fp_less_eq(m1,ap::machineepsilon) )
    {
        result = 1.3862943611198906188E0-0.5*log(m1);
    }
    else
    {
        p = 1.37982864606273237150E-4;
        p = p*m1+2.28025724005875567385E-3;
        p = p*m1+7.97404013220415179367E-3;
        p = p*m1+9.85821379021226008714E-3;
        p = p*m1+6.87489687449949877925E-3;
        p = p*m1+6.18901033637687613229E-3;
        p = p*m1+8.79078273952743772254E-3;
        p = p*m1+1.49380448916805252718E-2;
        p = p*m1+3.08851465246711995998E-2;
        p = p*m1+9.65735902811690126535E-2;
        p = p*m1+1.38629436111989062502E0;
        q = 2.94078955048598507511E-5;
        q = q*m1+9.14184723865917226571E-4;
        q = q*m1+5.94058303753167793257E-3;
        q = q*m1+1.54850516649762399335E-2;
        q = q*m1+2.39089602715924892727E-2;
        q = q*m1+3.01204715227604046988E-2;
        q = q*m1+3.73774314173823228969E-2;
        q = q*m1+4.88280347570998239232E-2;
        q = q*m1+7.03124996963957469739E-2;
        q = q*m1+1.24999999999870820058E-1;
        q = q*m1+4.99999999999999999821E-1;
        result = p-q*log(m1);
    }
    return result;
}


/*************************************************************************
Incomplete elliptic integral of the first kind F(phi|m)

Approximates the integral



               phi
                -
               | |
               |           dt
F(phi_\m)  =    |    ------------------
               |                   2
             | |    sqrt( 1 - m sin t )
              -
               0

of amplitude phi and modulus m, using the arithmetic -
geometric mean algorithm.




ACCURACY:

Tested at random points with m in [0, 1] and phi as indicated.

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE     -10,10       200000      7.4e-16     1.0e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double incompleteellipticintegralk(double phi, double m)
{
    double result;
    double a;
    double b;
    double c;
    double e;
    double temp;
    double pio2;
    double t;
    double k;
    int d;
    int md;
    int s;
    int npio2;

    pio2 = 1.57079632679489661923;
    if( ap::fp_eq(m,0) )
    {
        result = phi;
        return result;
    }
    a = 1-m;
    if( ap::fp_eq(a,0) )
    {
        result = log(tan(0.5*(pio2+phi)));
        return result;
    }
    npio2 = ap::ifloor(phi/pio2);
    if( npio2%2!=0 )
    {
        npio2 = npio2+1;
    }
    if( npio2!=0 )
    {
        k = ellipticintegralk(1-a);
        phi = phi-npio2*pio2;
    }
    else
    {
        k = 0;
    }
    if( ap::fp_less(phi,0) )
    {
        phi = -phi;
        s = -1;
    }
    else
    {
        s = 0;
    }
    b = sqrt(a);
    t = tan(phi);
    if( ap::fp_greater(fabs(t),10) )
    {
        e = 1.0/(b*t);
        if( ap::fp_less(fabs(e),10) )
        {
            e = atan(e);
            if( npio2==0 )
            {
                k = ellipticintegralk(1-a);
            }
            temp = k-incompleteellipticintegralk(e, m);
            if( s<0 )
            {
                temp = -temp;
            }
            result = temp+npio2*k;
            return result;
        }
    }
    a = 1.0;
    c = sqrt(m);
    d = 1;
    md = 0;
    while(ap::fp_greater(fabs(c/a),ap::machineepsilon))
    {
        temp = b/a;
        phi = phi+atan(t*temp)+md*ap::pi();
        md = ap::trunc((phi+pio2)/ap::pi());
        t = t*(1.0+temp)/(1.0-temp*t*t);
        c = 0.5*(a-b);
        temp = sqrt(a*b);
        a = 0.5*(a+b);
        b = temp;
        d = d+d;
    }
    temp = (atan(t)+md*ap::pi())/(d*a);
    if( s<0 )
    {
        temp = -temp;
    }
    result = temp+npio2*k;
    return result;
}


/*************************************************************************
Complete elliptic integral of the second kind

Approximates the integral


           pi/2
            -
           | |                 2
E(m)  =    |    sqrt( 1 - m sin t ) dt
         | |
          -
           0

using the approximation

     P(x)  -  x log x Q(x).

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE       0, 1       10000       2.1e-16     7.3e-17

Cephes Math Library, Release 2.8: June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*************************************************************************/
double ellipticintegrale(double m)
{
    double result;
    double p;
    double q;

    ap::ap_error::make_assertion(ap::fp_greater_eq(m,0)&&ap::fp_less_eq(m,1), "Domain error in EllipticIntegralE: m<0 or m>1");
    m = 1-m;
    if( ap::fp_eq(m,0) )
    {
        result = 1;
        return result;
    }
    p = 1.53552577301013293365E-4;
    p = p*m+2.50888492163602060990E-3;
    p = p*m+8.68786816565889628429E-3;
    p = p*m+1.07350949056076193403E-2;
    p = p*m+7.77395492516787092951E-3;
    p = p*m+7.58395289413514708519E-3;
    p = p*m+1.15688436810574127319E-2;
    p = p*m+2.18317996015557253103E-2;
    p = p*m+5.68051945617860553470E-2;
    p = p*m+4.43147180560990850618E-1;
    p = p*m+1.00000000000000000299E0;
    q = 3.27954898576485872656E-5;
    q = q*m+1.00962792679356715133E-3;
    q = q*m+6.50609489976927491433E-3;
    q = q*m+1.68862163993311317300E-2;
    q = q*m+2.61769742454493659583E-2;
    q = q*m+3.34833904888224918614E-2;
    q = q*m+4.27180926518931511717E-2;
    q = q*m+5.85936634471101055642E-2;
    q = q*m+9.37499997197644278445E-2;
    q = q*m+2.49999999999888314361E-1;
    result = p-q*m*log(m);
    return result;
}


/*************************************************************************
Incomplete elliptic integral of the second kind

Approximates the integral


               phi
                -
               | |
               |                   2
E(phi_\m)  =    |    sqrt( 1 - m sin t ) dt
               |
             | |
              -
               0

of amplitude phi and modulus m, using the arithmetic -
geometric mean algorithm.

ACCURACY:

Tested at random arguments with phi in [-10, 10] and m in
[0, 1].
                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE     -10,10      150000       3.3e-15     1.4e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1993, 2000 by Stephen L. Moshier
*************************************************************************/
double incompleteellipticintegrale(double phi, double m)
{
    double result;
    double pio2;
    double a;
    double b;
    double c;
    double e;
    double temp;
    double lphi;
    double t;
    double ebig;
    int d;
    int md;
    int npio2;
    int s;

    pio2 = 1.57079632679489661923;
    if( ap::fp_eq(m,0) )
    {
        result = phi;
        return result;
    }
    lphi = phi;
    npio2 = ap::ifloor(lphi/pio2);
    if( npio2%2!=0 )
    {
        npio2 = npio2+1;
    }
    lphi = lphi-npio2*pio2;
    if( ap::fp_less(lphi,0) )
    {
        lphi = -lphi;
        s = -1;
    }
    else
    {
        s = 1;
    }
    a = 1.0-m;
    ebig = ellipticintegrale(m);
    if( ap::fp_eq(a,0) )
    {
        temp = sin(lphi);
        if( s<0 )
        {
            temp = -temp;
        }
        result = temp+npio2*ebig;
        return result;
    }
    t = tan(lphi);
    b = sqrt(a);
    
    //
    // Thanks to Brian Fitzgerald <fitzgb@mml0.meche.rpi.edu>
    // for pointing out an instability near odd multiples of pi/2
    //
    if( ap::fp_greater(fabs(t),10) )
    {
        
        //
        // Transform the amplitude
        //
        e = 1.0/(b*t);
        
        //
        // ... but avoid multiple recursions.
        //
        if( ap::fp_less(fabs(e),10) )
        {
            e = atan(e);
            temp = ebig+m*sin(lphi)*sin(e)-incompleteellipticintegrale(e, m);
            if( s<0 )
            {
                temp = -temp;
            }
            result = temp+npio2*ebig;
            return result;
        }
    }
    c = sqrt(m);
    a = 1.0;
    d = 1;
    e = 0.0;
    md = 0;
    while(ap::fp_greater(fabs(c/a),ap::machineepsilon))
    {
        temp = b/a;
        lphi = lphi+atan(t*temp)+md*ap::pi();
        md = ap::trunc((lphi+pio2)/ap::pi());
        t = t*(1.0+temp)/(1.0-temp*t*t);
        c = 0.5*(a-b);
        temp = sqrt(a*b);
        a = 0.5*(a+b);
        b = temp;
        d = d+d;
        e = e+c*sin(lphi);
    }
    temp = ebig/ellipticintegralk(m);
    temp = temp*((atan(t)+md*ap::pi())/(d*a));
    temp = temp+e;
    if( s<0 )
    {
        temp = -temp;
    }
    result = temp+npio2*ebig;
    return result;
    return result;
}




