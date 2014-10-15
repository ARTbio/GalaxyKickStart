/*************************************************************************
Copyright (c) 2007-2009, Sergey Bochkanov (ALGLIB project).

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
#include "ratint.h"

static const int brcvnum = 10;

static void barycentricnormalize(barycentricinterpolant& b);
static void barycentriccalcbasis(const barycentricinterpolant& b,
     double t,
     ap::real_1d_array& y);
static void barycentricfitwcfixedd(ap::real_1d_array x,
     ap::real_1d_array y,
     const ap::real_1d_array& w,
     int n,
     ap::real_1d_array xc,
     ap::real_1d_array yc,
     const ap::integer_1d_array& dc,
     int k,
     int m,
     int d,
     int& info,
     barycentricinterpolant& b,
     barycentricfitreport& rep);

/*************************************************************************
Rational interpolation using barycentric formula

F(t) = SUM(i=0,n-1,w[i]*f[i]/(t-x[i])) / SUM(i=0,n-1,w[i]/(t-x[i]))

Input parameters:
    B   -   barycentric interpolant built with one of model building
            subroutines.
    T   -   interpolation point

Result:
    barycentric interpolant F(t)

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
double barycentriccalc(const barycentricinterpolant& b, double t)
{
    double result;
    double s1;
    double s2;
    double s;
    double v;
    int i;

    
    //
    // special case: N=1
    //
    if( b.n==1 )
    {
        result = b.sy*b.y(0);
        return result;
    }
    
    //
    // Here we assume that task is normalized, i.e.:
    // 1. abs(Y[i])<=1
    // 2. abs(W[i])<=1
    // 3. X[] is ordered
    //
    s = fabs(t-b.x(0));
    for(i = 0; i <= b.n-1; i++)
    {
        v = b.x(i);
        if( ap::fp_eq(v,t) )
        {
            result = b.sy*b.y(i);
            return result;
        }
        v = fabs(t-v);
        if( ap::fp_less(v,s) )
        {
            s = v;
        }
    }
    s1 = 0;
    s2 = 0;
    for(i = 0; i <= b.n-1; i++)
    {
        v = s/(t-b.x(i));
        v = v*b.w(i);
        s1 = s1+v*b.y(i);
        s2 = s2+v;
    }
    result = b.sy*s1/s2;
    return result;
}


/*************************************************************************
Differentiation of barycentric interpolant: first derivative.

Algorithm used in this subroutine is very robust and should not fail until
provided with values too close to MaxRealNumber  (usually  MaxRealNumber/N
or greater will overflow).

INPUT PARAMETERS:
    B   -   barycentric interpolant built with one of model building
            subroutines.
    T   -   interpolation point

OUTPUT PARAMETERS:
    F   -   barycentric interpolant at T
    DF  -   first derivative
    
NOTE


  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void barycentricdiff1(const barycentricinterpolant& b,
     double t,
     double& f,
     double& df)
{
    double v;
    double vv;
    int i;
    int k;
    double n0;
    double n1;
    double d0;
    double d1;
    double s0;
    double s1;
    double xk;
    double xi;
    double xmin;
    double xmax;
    double xscale1;
    double xoffs1;
    double xscale2;
    double xoffs2;
    double xprev;

    
    //
    // special case: N=1
    //
    if( b.n==1 )
    {
        f = b.sy*b.y(0);
        df = 0;
        return;
    }
    if( ap::fp_eq(b.sy,0) )
    {
        f = 0;
        df = 0;
        return;
    }
    ap::ap_error::make_assertion(ap::fp_greater(b.sy,0), "BarycentricDiff1: internal error");
    
    //
    // We assume than N>1 and B.SY>0. Find:
    // 1. pivot point (X[i] closest to T)
    // 2. width of interval containing X[i]
    //
    v = fabs(b.x(0)-t);
    k = 0;
    xmin = b.x(0);
    xmax = b.x(0);
    for(i = 1; i <= b.n-1; i++)
    {
        vv = b.x(i);
        if( ap::fp_less(fabs(vv-t),v) )
        {
            v = fabs(vv-t);
            k = i;
        }
        xmin = ap::minreal(xmin, vv);
        xmax = ap::maxreal(xmax, vv);
    }
    
    //
    // pivot point found, calculate dNumerator and dDenominator
    //
    xscale1 = 1/(xmax-xmin);
    xoffs1 = -xmin/(xmax-xmin)+1;
    xscale2 = 2;
    xoffs2 = -3;
    t = t*xscale1+xoffs1;
    t = t*xscale2+xoffs2;
    xk = b.x(k);
    xk = xk*xscale1+xoffs1;
    xk = xk*xscale2+xoffs2;
    v = t-xk;
    n0 = 0;
    n1 = 0;
    d0 = 0;
    d1 = 0;
    xprev = -2;
    for(i = 0; i <= b.n-1; i++)
    {
        xi = b.x(i);
        xi = xi*xscale1+xoffs1;
        xi = xi*xscale2+xoffs2;
        ap::ap_error::make_assertion(ap::fp_greater(xi,xprev), "BarycentricDiff1: points are too close!");
        xprev = xi;
        if( i!=k )
        {
            vv = ap::sqr(t-xi);
            s0 = (t-xk)/(t-xi);
            s1 = (xk-xi)/vv;
        }
        else
        {
            s0 = 1;
            s1 = 0;
        }
        vv = b.w(i)*b.y(i);
        n0 = n0+s0*vv;
        n1 = n1+s1*vv;
        vv = b.w(i);
        d0 = d0+s0*vv;
        d1 = d1+s1*vv;
    }
    f = b.sy*n0/d0;
    df = (n1*d0-n0*d1)/ap::sqr(d0);
    if( ap::fp_neq(df,0) )
    {
        df = ap::sign(df)*exp(log(fabs(df))+log(b.sy)+log(xscale1)+log(xscale2));
    }
}


/*************************************************************************
Differentiation of barycentric interpolant: first/second derivatives.

INPUT PARAMETERS:
    B   -   barycentric interpolant built with one of model building
            subroutines.
    T   -   interpolation point

OUTPUT PARAMETERS:
    F   -   barycentric interpolant at T
    DF  -   first derivative
    D2F -   second derivative

NOTE: this algorithm may fail due to overflow/underflor if  used  on  data
whose values are close to MaxRealNumber or MinRealNumber.  Use more robust
BarycentricDiff1() subroutine in such cases.


  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void barycentricdiff2(const barycentricinterpolant& b,
     double t,
     double& f,
     double& df,
     double& d2f)
{
    double v;
    double vv;
    int i;
    int k;
    double n0;
    double n1;
    double n2;
    double d0;
    double d1;
    double d2;
    double s0;
    double s1;
    double s2;
    double xk;
    double xi;

    f = 0;
    df = 0;
    d2f = 0;
    
    //
    // special case: N=1
    //
    if( b.n==1 )
    {
        f = b.sy*b.y(0);
        df = 0;
        d2f = 0;
        return;
    }
    if( ap::fp_eq(b.sy,0) )
    {
        f = 0;
        df = 0;
        d2f = 0;
        return;
    }
    ap::ap_error::make_assertion(ap::fp_greater(b.sy,0), "BarycentricDiff: internal error");
    
    //
    // We assume than N>1 and B.SY>0. Find:
    // 1. pivot point (X[i] closest to T)
    // 2. width of interval containing X[i]
    //
    v = fabs(b.x(0)-t);
    k = 0;
    for(i = 1; i <= b.n-1; i++)
    {
        vv = b.x(i);
        if( ap::fp_less(fabs(vv-t),v) )
        {
            v = fabs(vv-t);
            k = i;
        }
    }
    
    //
    // pivot point found, calculate dNumerator and dDenominator
    //
    xk = b.x(k);
    v = t-xk;
    n0 = 0;
    n1 = 0;
    n2 = 0;
    d0 = 0;
    d1 = 0;
    d2 = 0;
    for(i = 0; i <= b.n-1; i++)
    {
        if( i!=k )
        {
            xi = b.x(i);
            vv = ap::sqr(t-xi);
            s0 = (t-xk)/(t-xi);
            s1 = (xk-xi)/vv;
            s2 = -2*(xk-xi)/(vv*(t-xi));
        }
        else
        {
            s0 = 1;
            s1 = 0;
            s2 = 0;
        }
        vv = b.w(i)*b.y(i);
        n0 = n0+s0*vv;
        n1 = n1+s1*vv;
        n2 = n2+s2*vv;
        vv = b.w(i);
        d0 = d0+s0*vv;
        d1 = d1+s1*vv;
        d2 = d2+s2*vv;
    }
    f = b.sy*n0/d0;
    df = b.sy*(n1*d0-n0*d1)/ap::sqr(d0);
    d2f = b.sy*((n2*d0-n0*d2)*ap::sqr(d0)-(n1*d0-n0*d1)*2*d0*d1)/ap::sqr(ap::sqr(d0));
}


/*************************************************************************
This subroutine performs linear transformation of the argument.

INPUT PARAMETERS:
    B       -   rational interpolant in barycentric form
    CA, CB  -   transformation coefficients: x = CA*t + CB

OUTPUT PARAMETERS:
    B       -   transformed interpolant with X replaced by T

  -- ALGLIB PROJECT --
     Copyright 19.08.2009 by Bochkanov Sergey
*************************************************************************/
void barycentriclintransx(barycentricinterpolant& b, double ca, double cb)
{
    int i;
    int j;
    double v;

    
    //
    // special case, replace by constant F(CB)
    //
    if( ap::fp_eq(ca,0) )
    {
        b.sy = barycentriccalc(b, cb);
        v = 1;
        for(i = 0; i <= b.n-1; i++)
        {
            b.y(i) = 1;
            b.w(i) = v;
            v = -v;
        }
        return;
    }
    
    //
    // general case: CA<>0
    //
    for(i = 0; i <= b.n-1; i++)
    {
        b.x(i) = (b.x(i)-cb)/ca;
    }
    if( ap::fp_less(ca,0) )
    {
        for(i = 0; i <= b.n-1; i++)
        {
            if( i<b.n-1-i )
            {
                j = b.n-1-i;
                v = b.x(i);
                b.x(i) = b.x(j);
                b.x(j) = v;
                v = b.y(i);
                b.y(i) = b.y(j);
                b.y(j) = v;
                v = b.w(i);
                b.w(i) = b.w(j);
                b.w(j) = v;
            }
            else
            {
                break;
            }
        }
    }
}


/*************************************************************************
This  subroutine   performs   linear  transformation  of  the  barycentric
interpolant.

INPUT PARAMETERS:
    B       -   rational interpolant in barycentric form
    CA, CB  -   transformation coefficients: B2(x) = CA*B(x) + CB

OUTPUT PARAMETERS:
    B       -   transformed interpolant

  -- ALGLIB PROJECT --
     Copyright 19.08.2009 by Bochkanov Sergey
*************************************************************************/
void barycentriclintransy(barycentricinterpolant& b, double ca, double cb)
{
    int i;
    double v;

    for(i = 0; i <= b.n-1; i++)
    {
        b.y(i) = ca*b.sy*b.y(i)+cb;
    }
    b.sy = 0;
    for(i = 0; i <= b.n-1; i++)
    {
        b.sy = ap::maxreal(b.sy, fabs(b.y(i)));
    }
    if( ap::fp_greater(b.sy,0) )
    {
        v = 1/b.sy;
        ap::vmul(&b.y(0), 1, ap::vlen(0,b.n-1), v);
    }
}


/*************************************************************************
Extracts X/Y/W arrays from rational interpolant

INPUT PARAMETERS:
    B   -   barycentric interpolant

OUTPUT PARAMETERS:
    N   -   nodes count, N>0
    X   -   interpolation nodes, array[0..N-1]
    F   -   function values, array[0..N-1]
    W   -   barycentric weights, array[0..N-1]

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void barycentricunpack(const barycentricinterpolant& b,
     int& n,
     ap::real_1d_array& x,
     ap::real_1d_array& y,
     ap::real_1d_array& w)
{
    double v;

    n = b.n;
    x.setlength(n);
    y.setlength(n);
    w.setlength(n);
    v = b.sy;
    ap::vmove(&x(0), 1, &b.x(0), 1, ap::vlen(0,n-1));
    ap::vmove(&y(0), 1, &b.y(0), 1, ap::vlen(0,n-1), v);
    ap::vmove(&w(0), 1, &b.w(0), 1, ap::vlen(0,n-1));
}


/*************************************************************************
Serialization of the barycentric interpolant

INPUT PARAMETERS:
    B   -   barycentric interpolant

OUTPUT PARAMETERS:
    RA      -   array of real numbers which contains interpolant,
                array[0..RLen-1]
    RLen    -   RA lenght

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void barycentricserialize(const barycentricinterpolant& b,
     ap::real_1d_array& ra,
     int& ralen)
{

    ralen = 2+2+3*b.n;
    ra.setlength(ralen);
    ra(0) = ralen;
    ra(1) = brcvnum;
    ra(2) = b.n;
    ra(3) = b.sy;
    ap::vmove(&ra(4), 1, &b.x(0), 1, ap::vlen(4,4+b.n-1));
    ap::vmove(&ra(4+b.n), 1, &b.y(0), 1, ap::vlen(4+b.n,4+2*b.n-1));
    ap::vmove(&ra(4+2*b.n), 1, &b.w(0), 1, ap::vlen(4+2*b.n,4+3*b.n-1));
}


/*************************************************************************
Unserialization of the barycentric interpolant

INPUT PARAMETERS:
    RA  -   array of real numbers which contains interpolant,

OUTPUT PARAMETERS:
    B   -   barycentric interpolant

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void barycentricunserialize(const ap::real_1d_array& ra,
     barycentricinterpolant& b)
{

    ap::ap_error::make_assertion(ap::round(ra(1))==brcvnum, "BarycentricUnserialize: corrupted array!");
    b.n = ap::round(ra(2));
    b.sy = ra(3);
    b.x.setlength(b.n);
    b.y.setlength(b.n);
    b.w.setlength(b.n);
    ap::vmove(&b.x(0), 1, &ra(4), 1, ap::vlen(0,b.n-1));
    ap::vmove(&b.y(0), 1, &ra(4+b.n), 1, ap::vlen(0,b.n-1));
    ap::vmove(&b.w(0), 1, &ra(4+2*b.n), 1, ap::vlen(0,b.n-1));
}


/*************************************************************************
Copying of the barycentric interpolant

INPUT PARAMETERS:
    B   -   barycentric interpolant

OUTPUT PARAMETERS:
    B2  -   copy(B1)

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void barycentriccopy(const barycentricinterpolant& b,
     barycentricinterpolant& b2)
{

    b2.n = b.n;
    b2.sy = b.sy;
    b2.x.setlength(b2.n);
    b2.y.setlength(b2.n);
    b2.w.setlength(b2.n);
    ap::vmove(&b2.x(0), 1, &b.x(0), 1, ap::vlen(0,b2.n-1));
    ap::vmove(&b2.y(0), 1, &b.y(0), 1, ap::vlen(0,b2.n-1));
    ap::vmove(&b2.w(0), 1, &b.w(0), 1, ap::vlen(0,b2.n-1));
}


/*************************************************************************
Rational interpolant from X/Y/W arrays

F(t) = SUM(i=0,n-1,w[i]*f[i]/(t-x[i])) / SUM(i=0,n-1,w[i]/(t-x[i]))

INPUT PARAMETERS:
    X   -   interpolation nodes, array[0..N-1]
    F   -   function values, array[0..N-1]
    W   -   barycentric weights, array[0..N-1]
    N   -   nodes count, N>0

OUTPUT PARAMETERS:
    B   -   barycentric interpolant built from (X, Y, W)

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void barycentricbuildxyw(const ap::real_1d_array& x,
     const ap::real_1d_array& y,
     const ap::real_1d_array& w,
     int n,
     barycentricinterpolant& b)
{

    ap::ap_error::make_assertion(n>0, "BarycentricBuildXYW: incorrect N!");
    
    //
    // fill X/Y/W
    //
    b.x.setlength(n);
    b.y.setlength(n);
    b.w.setlength(n);
    ap::vmove(&b.x(0), 1, &x(0), 1, ap::vlen(0,n-1));
    ap::vmove(&b.y(0), 1, &y(0), 1, ap::vlen(0,n-1));
    ap::vmove(&b.w(0), 1, &w(0), 1, ap::vlen(0,n-1));
    b.n = n;
    
    //
    // Normalize
    //
    barycentricnormalize(b);
}


/*************************************************************************
Rational interpolant without poles

The subroutine constructs the rational interpolating function without real
poles  (see  'Barycentric rational interpolation with no  poles  and  high
rates of approximation', Michael S. Floater. and  Kai  Hormann,  for  more
information on this subject).

Input parameters:
    X   -   interpolation nodes, array[0..N-1].
    Y   -   function values, array[0..N-1].
    N   -   number of nodes, N>0.
    D   -   order of the interpolation scheme, 0 <= D <= N-1.
            D<0 will cause an error.
            D>=N it will be replaced with D=N-1.
            if you don't know what D to choose, use small value about 3-5.

Output parameters:
    B   -   barycentric interpolant.

Note:
    this algorithm always succeeds and calculates the weights  with  close
    to machine precision.

  -- ALGLIB PROJECT --
     Copyright 17.06.2007 by Bochkanov Sergey
*************************************************************************/
void barycentricbuildfloaterhormann(const ap::real_1d_array& x,
     const ap::real_1d_array& y,
     int n,
     int d,
     barycentricinterpolant& b)
{
    double s0;
    double s;
    double v;
    int i;
    int j;
    int k;
    ap::integer_1d_array perm;
    ap::real_1d_array wtemp;

    ap::ap_error::make_assertion(n>0, "BarycentricFloaterHormann: N<=0!");
    ap::ap_error::make_assertion(d>=0, "BarycentricFloaterHormann: incorrect D!");
    
    //
    // Prepare
    //
    if( d>n-1 )
    {
        d = n-1;
    }
    b.n = n;
    
    //
    // special case: N=1
    //
    if( n==1 )
    {
        b.x.setlength(n);
        b.y.setlength(n);
        b.w.setlength(n);
        b.x(0) = x(0);
        b.y(0) = y(0);
        b.w(0) = 1;
        barycentricnormalize(b);
        return;
    }
    
    //
    // Fill X/Y
    //
    b.x.setlength(n);
    b.y.setlength(n);
    ap::vmove(&b.x(0), 1, &x(0), 1, ap::vlen(0,n-1));
    ap::vmove(&b.y(0), 1, &y(0), 1, ap::vlen(0,n-1));
    tagsortfastr(b.x, b.y, n);
    
    //
    // Calculate Wk
    //
    b.w.setlength(n);
    s0 = 1;
    for(k = 1; k <= d; k++)
    {
        s0 = -s0;
    }
    for(k = 0; k <= n-1; k++)
    {
        
        //
        // Wk
        //
        s = 0;
        for(i = ap::maxint(k-d, 0); i <= ap::minint(k, n-1-d); i++)
        {
            v = 1;
            for(j = i; j <= i+d; j++)
            {
                if( j!=k )
                {
                    v = v/fabs(b.x(k)-b.x(j));
                }
            }
            s = s+v;
        }
        b.w(k) = s0*s;
        
        //
        // Next S0
        //
        s0 = -s0;
    }
    
    //
    // Normalize
    //
    barycentricnormalize(b);
}


/*************************************************************************
Weghted rational least  squares  fitting  using  Floater-Hormann  rational
functions  with  optimal  D  chosen  from  [0,9],  with  constraints   and
individual weights.

Equidistant  grid  with M node on [min(x),max(x)]  is  used to build basis
functions. Different values of D are tried, optimal D (least WEIGHTED root
mean square error) is chosen.  Task  is  linear,  so  linear least squares
solver  is  used.  Complexity  of  this  computational  scheme is O(N*M^2)
(mostly dominated by the least squares solver).

SEE ALSO
* BarycentricFitFloaterHormann(), "lightweight" fitting without invididual
  weights and constraints.

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    W   -   weights, array[0..N-1]
            Each summand in square  sum  of  approximation deviations from
            given  values  is  multiplied  by  the square of corresponding
            weight. Fill it by 1's if you don't  want  to  solve  weighted
            task.
    N   -   number of points, N>0.
    XC  -   points where function values/derivatives are constrained,
            array[0..K-1].
    YC  -   values of constraints, array[0..K-1]
    DC  -   array[0..K-1], types of constraints:
            * DC[i]=0   means that S(XC[i])=YC[i]
            * DC[i]=1   means that S'(XC[i])=YC[i]
            SEE BELOW FOR IMPORTANT INFORMATION ON CONSTRAINTS
    K   -   number of constraints, 0<=K<M.
            K=0 means no constraints (XC/YC/DC are not used in such cases)
    M   -   number of basis functions ( = number_of_nodes), M>=2.

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearWC() subroutine.
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
                        -3 means inconsistent constraints
                        -1 means another errors in parameters passed
                           (N<=0, for example)
    B   -   barycentric interpolant.
    Rep -   report, same format as in LSFitLinearWC() subroutine.
            Following fields are set:
            * DBest         best value of the D parameter
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroitine doesn't calculate task's condition number for K<>0.

SETTING CONSTRAINTS - DANGERS AND OPPORTUNITIES:

Setting constraints can lead  to undesired  results,  like ill-conditioned
behavior, or inconsistency being detected. From the other side,  it allows
us to improve quality of the fit. Here we summarize  our  experience  with
constrained barycentric interpolants:
* excessive  constraints  can  be  inconsistent.   Floater-Hormann   basis
  functions aren't as flexible as splines (although they are very smooth).
* the more evenly constraints are spread across [min(x),max(x)],  the more
  chances that they will be consistent
* the  greater  is  M (given  fixed  constraints),  the  more chances that
  constraints will be consistent
* in the general case, consistency of constraints IS NOT GUARANTEED.
* in the several special cases, however, we CAN guarantee consistency.
* one of this cases is constraints on the function  VALUES at the interval
  boundaries. Note that consustency of the  constraints  on  the  function
  DERIVATIVES is NOT guaranteed (you can use in such cases  cubic  splines
  which are more flexible).
* another  special  case  is ONE constraint on the function value (OR, but
  not AND, derivative) anywhere in the interval

Our final recommendation is to use constraints  WHEN  AND  ONLY  WHEN  you
can't solve your task without them. Anything beyond  special  cases  given
above is not guaranteed and may result in inconsistency.

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************/
void barycentricfitfloaterhormannwc(const ap::real_1d_array& x,
     const ap::real_1d_array& y,
     const ap::real_1d_array& w,
     int n,
     const ap::real_1d_array& xc,
     const ap::real_1d_array& yc,
     const ap::integer_1d_array& dc,
     int k,
     int m,
     int& info,
     barycentricinterpolant& b,
     barycentricfitreport& rep)
{
    int d;
    int i;
    double wrmscur;
    double wrmsbest;
    barycentricinterpolant locb;
    barycentricfitreport locrep;
    int locinfo;

    if( n<1||m<2||k<0||k>=m )
    {
        info = -1;
        return;
    }
    
    //
    // Find optimal D
    //
    // Info is -3 by default (degenerate constraints).
    // If LocInfo will always be equal to -3, Info will remain equal to -3.
    // If at least once LocInfo will be -4, Info will be -4.
    //
    wrmsbest = ap::maxrealnumber;
    rep.dbest = -1;
    info = -3;
    for(d = 0; d <= ap::minint(9, n-1); d++)
    {
        barycentricfitwcfixedd(x, y, w, n, xc, yc, dc, k, m, d, locinfo, locb, locrep);
        ap::ap_error::make_assertion(locinfo==-4||locinfo==-3||locinfo>0, "BarycentricFitFloaterHormannWC: unexpected result from BarycentricFitWCFixedD!");
        if( locinfo>0 )
        {
            
            //
            // Calculate weghted RMS
            //
            wrmscur = 0;
            for(i = 0; i <= n-1; i++)
            {
                wrmscur = wrmscur+ap::sqr(w(i)*(y(i)-barycentriccalc(locb, x(i))));
            }
            wrmscur = sqrt(wrmscur/n);
            if( ap::fp_less(wrmscur,wrmsbest)||rep.dbest<0 )
            {
                barycentriccopy(locb, b);
                rep.dbest = d;
                info = 1;
                rep.rmserror = locrep.rmserror;
                rep.avgerror = locrep.avgerror;
                rep.avgrelerror = locrep.avgrelerror;
                rep.maxerror = locrep.maxerror;
                rep.taskrcond = locrep.taskrcond;
                wrmsbest = wrmscur;
            }
        }
        else
        {
            if( locinfo!=-3&&info<0 )
            {
                info = locinfo;
            }
        }
    }
}


/*************************************************************************
Rational least squares fitting, without weights and constraints.

See BarycentricFitFloaterHormannWC() for more information.

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************/
void barycentricfitfloaterhormann(const ap::real_1d_array& x,
     const ap::real_1d_array& y,
     int n,
     int m,
     int& info,
     barycentricinterpolant& b,
     barycentricfitreport& rep)
{
    ap::real_1d_array w;
    ap::real_1d_array xc;
    ap::real_1d_array yc;
    ap::integer_1d_array dc;
    int i;

    if( n<1 )
    {
        info = -1;
        return;
    }
    w.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        w(i) = 1;
    }
    barycentricfitfloaterhormannwc(x, y, w, n, xc, yc, dc, 0, m, info, b, rep);
}


/*************************************************************************
Normalization of barycentric interpolant:
* B.N, B.X, B.Y and B.W are initialized
* B.SY is NOT initialized
* Y[] is normalized, scaling coefficient is stored in B.SY
* W[] is normalized, no scaling coefficient is stored
* X[] is sorted

Internal subroutine.
*************************************************************************/
static void barycentricnormalize(barycentricinterpolant& b)
{
    ap::integer_1d_array p1;
    ap::integer_1d_array p2;
    int i;
    int j;
    int j2;
    double v;

    
    //
    // Normalize task: |Y|<=1, |W|<=1, sort X[]
    //
    b.sy = 0;
    for(i = 0; i <= b.n-1; i++)
    {
        b.sy = ap::maxreal(b.sy, fabs(b.y(i)));
    }
    if( ap::fp_greater(b.sy,0)&&ap::fp_greater(fabs(b.sy-1),10*ap::machineepsilon) )
    {
        v = 1/b.sy;
        ap::vmul(&b.y(0), 1, ap::vlen(0,b.n-1), v);
    }
    v = 0;
    for(i = 0; i <= b.n-1; i++)
    {
        v = ap::maxreal(v, fabs(b.w(i)));
    }
    if( ap::fp_greater(v,0)&&ap::fp_greater(fabs(v-1),10*ap::machineepsilon) )
    {
        v = 1/v;
        ap::vmul(&b.w(0), 1, ap::vlen(0,b.n-1), v);
    }
    for(i = 0; i <= b.n-2; i++)
    {
        if( ap::fp_less(b.x(i+1),b.x(i)) )
        {
            tagsort(b.x, b.n, p1, p2);
            for(j = 0; j <= b.n-1; j++)
            {
                j2 = p2(j);
                v = b.y(j);
                b.y(j) = b.y(j2);
                b.y(j2) = v;
                v = b.w(j);
                b.w(j) = b.w(j2);
                b.w(j2) = v;
            }
            break;
        }
    }
}


/*************************************************************************
Internal subroutine, calculates barycentric basis functions.
Used for efficient simultaneous calculation of N basis functions.

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
static void barycentriccalcbasis(const barycentricinterpolant& b,
     double t,
     ap::real_1d_array& y)
{
    double s2;
    double s;
    double v;
    int i;
    int j;

    
    //
    // special case: N=1
    //
    if( b.n==1 )
    {
        y(0) = 1;
        return;
    }
    
    //
    // Here we assume that task is normalized, i.e.:
    // 1. abs(Y[i])<=1
    // 2. abs(W[i])<=1
    // 3. X[] is ordered
    //
    // First, we decide: should we use "safe" formula (guarded
    // against overflow) or fast one?
    //
    s = fabs(t-b.x(0));
    for(i = 0; i <= b.n-1; i++)
    {
        v = b.x(i);
        if( ap::fp_eq(v,t) )
        {
            for(j = 0; j <= b.n-1; j++)
            {
                y(j) = 0;
            }
            y(i) = 1;
            return;
        }
        v = fabs(t-v);
        if( ap::fp_less(v,s) )
        {
            s = v;
        }
    }
    s2 = 0;
    for(i = 0; i <= b.n-1; i++)
    {
        v = s/(t-b.x(i));
        v = v*b.w(i);
        y(i) = v;
        s2 = s2+v;
    }
    v = 1/s2;
    ap::vmul(&y(0), 1, ap::vlen(0,b.n-1), v);
}


/*************************************************************************
Internal Floater-Hormann fitting subroutine for fixed D
*************************************************************************/
static void barycentricfitwcfixedd(ap::real_1d_array x,
     ap::real_1d_array y,
     const ap::real_1d_array& w,
     int n,
     ap::real_1d_array xc,
     ap::real_1d_array yc,
     const ap::integer_1d_array& dc,
     int k,
     int m,
     int d,
     int& info,
     barycentricinterpolant& b,
     barycentricfitreport& rep)
{
    ap::real_2d_array fmatrix;
    ap::real_2d_array cmatrix;
    ap::real_1d_array y2;
    ap::real_1d_array w2;
    ap::real_1d_array sx;
    ap::real_1d_array sy;
    ap::real_1d_array sbf;
    ap::real_1d_array xoriginal;
    ap::real_1d_array yoriginal;
    ap::real_1d_array tmp;
    lsfitreport lrep;
    double v0;
    double v1;
    double mx;
    barycentricinterpolant b2;
    int i;
    int j;
    int relcnt;
    double xa;
    double xb;
    double sa;
    double sb;
    double decay;

    if( n<1||m<2||k<0||k>=m )
    {
        info = -1;
        return;
    }
    for(i = 0; i <= k-1; i++)
    {
        info = 0;
        if( dc(i)<0 )
        {
            info = -1;
        }
        if( dc(i)>1 )
        {
            info = -1;
        }
        if( info<0 )
        {
            return;
        }
    }
    
    //
    // weight decay for correct handling of task which becomes
    // degenerate after constraints are applied
    //
    decay = 10000*ap::machineepsilon;
    
    //
    // Scale X, Y, XC, YC
    //
    lsfitscalexy(x, y, n, xc, yc, dc, k, xa, xb, sa, sb, xoriginal, yoriginal);
    
    //
    // allocate space, initialize:
    // * FMatrix-   values of basis functions at X[]
    // * CMatrix-   values (derivatives) of basis functions at XC[]
    //
    y2.setlength(n+m);
    w2.setlength(n+m);
    fmatrix.setlength(n+m, m);
    if( k>0 )
    {
        cmatrix.setlength(k, m+1);
    }
    y2.setlength(n+m);
    w2.setlength(n+m);
    
    //
    // Prepare design and constraints matrices:
    // * fill constraints matrix
    // * fill first N rows of design matrix with values
    // * fill next M rows of design matrix with regularizing term
    // * append M zeros to Y
    // * append M elements, mean(abs(W)) each, to W
    //
    sx.setlength(m);
    sy.setlength(m);
    sbf.setlength(m);
    for(j = 0; j <= m-1; j++)
    {
        sx(j) = double(2*j)/double(m-1)-1;
    }
    for(i = 0; i <= m-1; i++)
    {
        sy(i) = 1;
    }
    barycentricbuildfloaterhormann(sx, sy, m, d, b2);
    mx = 0;
    for(i = 0; i <= n-1; i++)
    {
        barycentriccalcbasis(b2, x(i), sbf);
        ap::vmove(&fmatrix(i, 0), 1, &sbf(0), 1, ap::vlen(0,m-1));
        y2(i) = y(i);
        w2(i) = w(i);
        mx = mx+fabs(w(i))/n;
    }
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= m-1; j++)
        {
            if( i==j )
            {
                fmatrix(n+i,j) = decay;
            }
            else
            {
                fmatrix(n+i,j) = 0;
            }
        }
        y2(n+i) = 0;
        w2(n+i) = mx;
    }
    if( k>0 )
    {
        for(j = 0; j <= m-1; j++)
        {
            for(i = 0; i <= m-1; i++)
            {
                sy(i) = 0;
            }
            sy(j) = 1;
            barycentricbuildfloaterhormann(sx, sy, m, d, b2);
            for(i = 0; i <= k-1; i++)
            {
                ap::ap_error::make_assertion(dc(i)>=0&&dc(i)<=1, "BarycentricFit: internal error!");
                barycentricdiff1(b2, xc(i), v0, v1);
                if( dc(i)==0 )
                {
                    cmatrix(i,j) = v0;
                }
                if( dc(i)==1 )
                {
                    cmatrix(i,j) = v1;
                }
            }
        }
        for(i = 0; i <= k-1; i++)
        {
            cmatrix(i,m) = yc(i);
        }
    }
    
    //
    // Solve constrained task
    //
    if( k>0 )
    {
        
        //
        // solve using regularization
        //
        lsfitlinearwc(y2, w2, fmatrix, cmatrix, n+m, m, k, info, tmp, lrep);
    }
    else
    {
        
        //
        // no constraints, no regularization needed
        //
        lsfitlinearwc(y, w, fmatrix, cmatrix, n, m, k, info, tmp, lrep);
    }
    if( info<0 )
    {
        return;
    }
    
    //
    // Generate interpolant and scale it
    //
    ap::vmove(&sy(0), 1, &tmp(0), 1, ap::vlen(0,m-1));
    barycentricbuildfloaterhormann(sx, sy, m, d, b);
    barycentriclintransx(b, 2/(xb-xa), -(xa+xb)/(xb-xa));
    barycentriclintransy(b, sb-sa, sa);
    
    //
    // Scale absolute errors obtained from LSFitLinearW.
    // Relative error should be calculated separately
    // (because of shifting/scaling of the task)
    //
    rep.taskrcond = lrep.taskrcond;
    rep.rmserror = lrep.rmserror*(sb-sa);
    rep.avgerror = lrep.avgerror*(sb-sa);
    rep.maxerror = lrep.maxerror*(sb-sa);
    rep.avgrelerror = 0;
    relcnt = 0;
    for(i = 0; i <= n-1; i++)
    {
        if( ap::fp_neq(yoriginal(i),0) )
        {
            rep.avgrelerror = rep.avgrelerror+fabs(barycentriccalc(b, xoriginal(i))-yoriginal(i))/fabs(yoriginal(i));
            relcnt = relcnt+1;
        }
    }
    if( relcnt!=0 )
    {
        rep.avgrelerror = rep.avgrelerror/relcnt;
    }
}




