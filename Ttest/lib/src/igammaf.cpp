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
#include "igammaf.h"

/*************************************************************************
Incomplete gamma integral

The function is defined by

                          x
                           -
                  1       | |  -t  a-1
 igam(a,x)  =   -----     |   e   t   dt.
                 -      | |
                | (a)    -
                          0


In this implementation both arguments must be positive.
The integral is evaluated by either a power series or
continued fraction expansion, depending on the relative
values of a and x.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,30       200000       3.6e-14     2.9e-15
   IEEE      0,100      300000       9.9e-14     1.5e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double incompletegamma(double a, double x)
{
    double result;
    double igammaepsilon;
    double ans;
    double ax;
    double c;
    double r;
    double tmp;

    igammaepsilon = 0.000000000000001;
    if( ap::fp_less_eq(x,0)||ap::fp_less_eq(a,0) )
    {
        result = 0;
        return result;
    }
    if( ap::fp_greater(x,1)&&ap::fp_greater(x,a) )
    {
        result = 1-incompletegammac(a, x);
        return result;
    }
    ax = a*log(x)-x-lngamma(a, tmp);
    if( ap::fp_less(ax,-709.78271289338399) )
    {
        result = 0;
        return result;
    }
    ax = exp(ax);
    r = a;
    c = 1;
    ans = 1;
    do
    {
        r = r+1;
        c = c*x/r;
        ans = ans+c;
    }
    while(ap::fp_greater(c/ans,igammaepsilon));
    result = ans*ax/a;
    return result;
}


/*************************************************************************
Complemented incomplete gamma integral

The function is defined by


 igamc(a,x)   =   1 - igam(a,x)

                           inf.
                             -
                    1       | |  -t  a-1
              =   -----     |   e   t   dt.
                   -      | |
                  | (a)    -
                            x


In this implementation both arguments must be positive.
The integral is evaluated by either a power series or
continued fraction expansion, depending on the relative
values of a and x.

ACCURACY:

Tested at random a, x.
               a         x                      Relative error:
arithmetic   domain   domain     # trials      peak         rms
   IEEE     0.5,100   0,100      200000       1.9e-14     1.7e-15
   IEEE     0.01,0.5  0,100      200000       1.4e-13     1.6e-15

Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double incompletegammac(double a, double x)
{
    double result;
    double igammaepsilon;
    double igammabignumber;
    double igammabignumberinv;
    double ans;
    double ax;
    double c;
    double yc;
    double r;
    double t;
    double y;
    double z;
    double pk;
    double pkm1;
    double pkm2;
    double qk;
    double qkm1;
    double qkm2;
    double tmp;

    igammaepsilon = 0.000000000000001;
    igammabignumber = 4503599627370496.0;
    igammabignumberinv = 2.22044604925031308085*0.0000000000000001;
    if( ap::fp_less_eq(x,0)||ap::fp_less_eq(a,0) )
    {
        result = 1;
        return result;
    }
    if( ap::fp_less(x,1)||ap::fp_less(x,a) )
    {
        result = 1-incompletegamma(a, x);
        return result;
    }
    ax = a*log(x)-x-lngamma(a, tmp);
    if( ap::fp_less(ax,-709.78271289338399) )
    {
        result = 0;
        return result;
    }
    ax = exp(ax);
    y = 1-a;
    z = x+y+1;
    c = 0;
    pkm2 = 1;
    qkm2 = x;
    pkm1 = x+1;
    qkm1 = z*x;
    ans = pkm1/qkm1;
    do
    {
        c = c+1;
        y = y+1;
        z = z+2;
        yc = y*c;
        pk = pkm1*z-pkm2*yc;
        qk = qkm1*z-qkm2*yc;
        if( ap::fp_neq(qk,0) )
        {
            r = pk/qk;
            t = fabs((ans-r)/r);
            ans = r;
        }
        else
        {
            t = 1;
        }
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        if( ap::fp_greater(fabs(pk),igammabignumber) )
        {
            pkm2 = pkm2*igammabignumberinv;
            pkm1 = pkm1*igammabignumberinv;
            qkm2 = qkm2*igammabignumberinv;
            qkm1 = qkm1*igammabignumberinv;
        }
    }
    while(ap::fp_greater(t,igammaepsilon));
    result = ans*ax;
    return result;
}


/*************************************************************************
Inverse of complemented imcomplete gamma integral

Given p, the function finds x such that

 igamc( a, x ) = p.

Starting with the approximate value

        3
 x = a t

 where

 t = 1 - d - ndtri(p) sqrt(d)

and

 d = 1/9a,

the routine performs up to 10 Newton iterations to find the
root of igamc(a,x) - p = 0.

ACCURACY:

Tested at random a, p in the intervals indicated.

               a        p                      Relative error:
arithmetic   domain   domain     # trials      peak         rms
   IEEE     0.5,100   0,0.5       100000       1.0e-14     1.7e-15
   IEEE     0.01,0.5  0,0.5       100000       9.0e-14     3.4e-15
   IEEE    0.5,10000  0,0.5        20000       2.3e-13     3.8e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double invincompletegammac(double a, double y0)
{
    double result;
    double igammaepsilon;
    double iinvgammabignumber;
    double x0;
    double x1;
    double x;
    double yl;
    double yh;
    double y;
    double d;
    double lgm;
    double dithresh;
    int i;
    int dir;
    double tmp;

    igammaepsilon = 0.000000000000001;
    iinvgammabignumber = 4503599627370496.0;
    x0 = iinvgammabignumber;
    yl = 0;
    x1 = 0;
    yh = 1;
    dithresh = 5*igammaepsilon;
    d = 1/(9*a);
    y = 1-d-invnormaldistribution(y0)*sqrt(d);
    x = a*y*y*y;
    lgm = lngamma(a, tmp);
    i = 0;
    while(i<10)
    {
        if( ap::fp_greater(x,x0)||ap::fp_less(x,x1) )
        {
            d = 0.0625;
            break;
        }
        y = incompletegammac(a, x);
        if( ap::fp_less(y,yl)||ap::fp_greater(y,yh) )
        {
            d = 0.0625;
            break;
        }
        if( ap::fp_less(y,y0) )
        {
            x0 = x;
            yl = y;
        }
        else
        {
            x1 = x;
            yh = y;
        }
        d = (a-1)*log(x)-x-lgm;
        if( ap::fp_less(d,-709.78271289338399) )
        {
            d = 0.0625;
            break;
        }
        d = -exp(d);
        d = (y-y0)/d;
        if( ap::fp_less(fabs(d/x),igammaepsilon) )
        {
            result = x;
            return result;
        }
        x = x-d;
        i = i+1;
    }
    if( ap::fp_eq(x0,iinvgammabignumber) )
    {
        if( ap::fp_less_eq(x,0) )
        {
            x = 1;
        }
        while(ap::fp_eq(x0,iinvgammabignumber))
        {
            x = (1+d)*x;
            y = incompletegammac(a, x);
            if( ap::fp_less(y,y0) )
            {
                x0 = x;
                yl = y;
                break;
            }
            d = d+d;
        }
    }
    d = 0.5;
    dir = 0;
    i = 0;
    while(i<400)
    {
        x = x1+d*(x0-x1);
        y = incompletegammac(a, x);
        lgm = (x0-x1)/(x1+x0);
        if( ap::fp_less(fabs(lgm),dithresh) )
        {
            break;
        }
        lgm = (y-y0)/y0;
        if( ap::fp_less(fabs(lgm),dithresh) )
        {
            break;
        }
        if( ap::fp_less_eq(x,0.0) )
        {
            break;
        }
        if( ap::fp_greater_eq(y,y0) )
        {
            x1 = x;
            yh = y;
            if( dir<0 )
            {
                dir = 0;
                d = 0.5;
            }
            else
            {
                if( dir>1 )
                {
                    d = 0.5*d+0.5;
                }
                else
                {
                    d = (y0-yl)/(yh-yl);
                }
            }
            dir = dir+1;
        }
        else
        {
            x0 = x;
            yl = y;
            if( dir>0 )
            {
                dir = 0;
                d = 0.5;
            }
            else
            {
                if( dir<-1 )
                {
                    d = 0.5*d;
                }
                else
                {
                    d = (y0-yl)/(yh-yl);
                }
            }
            dir = dir-1;
        }
        i = i+1;
    }
    result = x;
    return result;
}




