/*************************************************************************
Copyright (c) 2006-2009, Sergey Bochkanov (ALGLIB project).

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
#include "spline1d.h"

static const int spline1dvnum = 11;

static void spline1dfitinternal(int st,
     ap::real_1d_array x,
     ap::real_1d_array y,
     const ap::real_1d_array& w,
     int n,
     ap::real_1d_array xc,
     ap::real_1d_array yc,
     const ap::integer_1d_array& dc,
     int k,
     int m,
     int& info,
     spline1dinterpolant& s,
     spline1dfitreport& rep);
static void heapsortpoints(ap::real_1d_array& x, ap::real_1d_array& y, int n);
static void heapsortdpoints(ap::real_1d_array& x,
     ap::real_1d_array& y,
     ap::real_1d_array& d,
     int n);
static void solvetridiagonal(ap::real_1d_array a,
     ap::real_1d_array b,
     ap::real_1d_array c,
     ap::real_1d_array d,
     int n,
     ap::real_1d_array& x);
static void solvecyclictridiagonal(const ap::real_1d_array& a,
     ap::real_1d_array b,
     const ap::real_1d_array& c,
     const ap::real_1d_array& d,
     int n,
     ap::real_1d_array& x);
static double diffthreepoint(double t,
     double x0,
     double f0,
     double x1,
     double f1,
     double x2,
     double f2);

/*************************************************************************
This subroutine builds linear spline interpolant

INPUT PARAMETERS:
    X   -   spline nodes, array[0..N-1]
    Y   -   function values, array[0..N-1]
    N   -   points count, N>=2
    
OUTPUT PARAMETERS:
    C   -   spline interpolant


ORDER OF POINTS

Subroutine automatically sorts points, so caller may pass unsorted array.

  -- ALGLIB PROJECT --
     Copyright 24.06.2007 by Bochkanov Sergey
*************************************************************************/
void spline1dbuildlinear(ap::real_1d_array x,
     ap::real_1d_array y,
     int n,
     spline1dinterpolant& c)
{
    int i;

    ap::ap_error::make_assertion(n>1, "Spline1DBuildLinear: N<2!");
    
    //
    // Sort points
    //
    heapsortpoints(x, y, n);
    
    //
    // Build
    //
    c.periodic = false;
    c.n = n;
    c.k = 3;
    c.x.setlength(n);
    c.c.setlength(4*(n-1));
    for(i = 0; i <= n-1; i++)
    {
        c.x(i) = x(i);
    }
    for(i = 0; i <= n-2; i++)
    {
        c.c(4*i+0) = y(i);
        c.c(4*i+1) = (y(i+1)-y(i))/(x(i+1)-x(i));
        c.c(4*i+2) = 0;
        c.c(4*i+3) = 0;
    }
}


/*************************************************************************
This subroutine builds cubic spline interpolant.

INPUT PARAMETERS:
    X           -   spline nodes, array[0..N-1].
    Y           -   function values, array[0..N-1].
    N           -   points count, N>=2
    BoundLType  -   boundary condition type for the left boundary
    BoundL      -   left boundary condition (first or second derivative,
                    depending on the BoundLType)
    BoundRType  -   boundary condition type for the right boundary
    BoundR      -   right boundary condition (first or second derivative,
                    depending on the BoundRType)

OUTPUT PARAMETERS:
    C           -   spline interpolant


ORDER OF POINTS

Subroutine automatically sorts points, so caller may pass unsorted array.

SETTING BOUNDARY VALUES:

The BoundLType/BoundRType parameters can have the following values:
    * -1, which corresonds to the periodic (cyclic) boundary conditions.
          In this case:
          * both BoundLType and BoundRType must be equal to -1.
          * BoundL/BoundR are ignored
          * Y[last] is ignored (it is assumed to be equal to Y[first]).
    *  0, which  corresponds  to  the  parabolically   terminated  spline
          (BoundL and/or BoundR are ignored).
    *  1, which corresponds to the first derivative boundary condition
    *  2, which corresponds to the second derivative boundary condition

PROBLEMS WITH PERIODIC BOUNDARY CONDITIONS:

Problems with periodic boundary conditions have Y[first_point]=Y[last_point].
However, this subroutine doesn't require you to specify equal  values  for
the first and last points - it automatically forces them to be equal.

  -- ALGLIB PROJECT --
     Copyright 23.06.2007 by Bochkanov Sergey
*************************************************************************/
void spline1dbuildcubic(ap::real_1d_array x,
     ap::real_1d_array y,
     int n,
     int boundltype,
     double boundl,
     int boundrtype,
     double boundr,
     spline1dinterpolant& c)
{
    ap::real_1d_array a1;
    ap::real_1d_array a2;
    ap::real_1d_array a3;
    ap::real_1d_array b;
    ap::real_1d_array d;
    ap::real_1d_array dt;
    int i;
    double v;

    ap::ap_error::make_assertion(n>=2, "Spline1DBuildCubic: N<2!");
    ap::ap_error::make_assertion(boundltype==-1||boundltype==0||boundltype==1||boundltype==2, "Spline1DBuildCubic: incorrect BoundLType!");
    ap::ap_error::make_assertion(boundrtype==-1||boundrtype==0||boundrtype==1||boundrtype==2, "Spline1DBuildCubic: incorrect BoundRType!");
    ap::ap_error::make_assertion(boundrtype==-1&&boundltype==-1||boundrtype!=-1&&boundltype!=-1, "Spline1DBuildCubic: incorrect BoundLType/BoundRType!");
    
    //
    // Special cases:
    // * N=2, parabolic terminated boundary condition on both ends
    // * N=2, periodic boundary condition
    //
    if( n==2&&boundltype==0&&boundrtype==0 )
    {
        
        //
        // Change task type
        //
        boundltype = 2;
        boundl = 0;
        boundrtype = 2;
        boundr = 0;
    }
    if( n==2&&boundltype==-1&&boundrtype==-1 )
    {
        
        //
        // Change task type
        //
        boundltype = 1;
        boundl = 0;
        boundrtype = 1;
        boundr = 0;
        y(1) = y(0);
    }
    
    //
    // Periodic and non-periodic boundary conditions are
    // two separate classes
    //
    if( boundrtype==-1&&boundltype==-1 )
    {
        
        //
        // Periodic boundary conditions
        //
        a1.setlength(n-1);
        a2.setlength(n-1);
        a3.setlength(n-1);
        b.setlength(n-1);
        
        //
        // Sort points.
        //
        heapsortpoints(x, y, n);
        y(n-1) = y(0);
        
        //
        // Boundary conditions at N-1 points
        // (one point less because last point is the same as first point).
        //
        a1(0) = x(1)-x(0);
        a2(0) = 2*(x(1)-x(0)+x(n-1)-x(n-2));
        a3(0) = x(n-1)-x(n-2);
        b(0) = 3*(y(n-1)-y(n-2))/(x(n-1)-x(n-2))*(x(1)-x(0))+3*(y(1)-y(0))/(x(1)-x(0))*(x(n-1)-x(n-2));
        for(i = 1; i <= n-2; i++)
        {
            
            //
            // Altough last point is [N-2], we use X[N-1] and Y[N-1]
            // (because of periodicity)
            //
            a1(i) = x(i+1)-x(i);
            a2(i) = 2*(x(i+1)-x(i-1));
            a3(i) = x(i)-x(i-1);
            b(i) = 3*(y(i)-y(i-1))/(x(i)-x(i-1))*(x(i+1)-x(i))+3*(y(i+1)-y(i))/(x(i+1)-x(i))*(x(i)-x(i-1));
        }
        
        //
        // Solve, add last point (with index N-1)
        //
        solvecyclictridiagonal(a1, a2, a3, b, n-1, dt);
        d.setlength(n);
        ap::vmove(&d(0), 1, &dt(0), 1, ap::vlen(0,n-2));
        d(n-1) = d(0);
        
        //
        // Now problem is reduced to the cubic Hermite spline
        //
        spline1dbuildhermite(x, y, d, n, c);
        c.periodic = true;
    }
    else
    {
        
        //
        // Non-periodic boundary condition
        //
        a1.setlength(n);
        a2.setlength(n);
        a3.setlength(n);
        b.setlength(n);
        
        //
        // Sort points.
        //
        heapsortpoints(x, y, n);
        
        //
        // Left boundary conditions
        //
        if( boundltype==0 )
        {
            a1(0) = 0;
            a2(0) = 1;
            a3(0) = 1;
            b(0) = 2*(y(1)-y(0))/(x(1)-x(0));
        }
        if( boundltype==1 )
        {
            a1(0) = 0;
            a2(0) = 1;
            a3(0) = 0;
            b(0) = boundl;
        }
        if( boundltype==2 )
        {
            a1(0) = 0;
            a2(0) = 2;
            a3(0) = 1;
            b(0) = 3*(y(1)-y(0))/(x(1)-x(0))-0.5*boundl*(x(1)-x(0));
        }
        
        //
        // Central conditions
        //
        for(i = 1; i <= n-2; i++)
        {
            a1(i) = x(i+1)-x(i);
            a2(i) = 2*(x(i+1)-x(i-1));
            a3(i) = x(i)-x(i-1);
            b(i) = 3*(y(i)-y(i-1))/(x(i)-x(i-1))*(x(i+1)-x(i))+3*(y(i+1)-y(i))/(x(i+1)-x(i))*(x(i)-x(i-1));
        }
        
        //
        // Right boundary conditions
        //
        if( boundrtype==0 )
        {
            a1(n-1) = 1;
            a2(n-1) = 1;
            a3(n-1) = 0;
            b(n-1) = 2*(y(n-1)-y(n-2))/(x(n-1)-x(n-2));
        }
        if( boundrtype==1 )
        {
            a1(n-1) = 0;
            a2(n-1) = 1;
            a3(n-1) = 0;
            b(n-1) = boundr;
        }
        if( boundrtype==2 )
        {
            a1(n-1) = 1;
            a2(n-1) = 2;
            a3(n-1) = 0;
            b(n-1) = 3*(y(n-1)-y(n-2))/(x(n-1)-x(n-2))+0.5*boundr*(x(n-1)-x(n-2));
        }
        
        //
        // Solve
        //
        solvetridiagonal(a1, a2, a3, b, n, d);
        
        //
        // Now problem is reduced to the cubic Hermite spline
        //
        spline1dbuildhermite(x, y, d, n, c);
    }
}


/*************************************************************************
This subroutine builds Catmull-Rom spline interpolant.

INPUT PARAMETERS:
    X           -   spline nodes, array[0..N-1].
    Y           -   function values, array[0..N-1].
    N           -   points count, N>=2
    BoundType   -   boundary condition type:
                    * -1 for periodic boundary condition
                    *  0 for parabolically terminated spline
    Tension     -   tension parameter:
                    * tension=0   corresponds to classic Catmull-Rom spline
                    * 0<tension<1 corresponds to more general form - cardinal spline

OUTPUT PARAMETERS:
    C           -   spline interpolant


ORDER OF POINTS

Subroutine automatically sorts points, so caller may pass unsorted array.

PROBLEMS WITH PERIODIC BOUNDARY CONDITIONS:

Problems with periodic boundary conditions have Y[first_point]=Y[last_point].
However, this subroutine doesn't require you to specify equal  values  for
the first and last points - it automatically forces them to be equal.

  -- ALGLIB PROJECT --
     Copyright 23.06.2007 by Bochkanov Sergey
*************************************************************************/
void spline1dbuildcatmullrom(ap::real_1d_array x,
     ap::real_1d_array y,
     int n,
     int boundtype,
     double tension,
     spline1dinterpolant& c)
{
    ap::real_1d_array a1;
    ap::real_1d_array a2;
    ap::real_1d_array a3;
    ap::real_1d_array b;
    ap::real_1d_array d;
    ap::real_1d_array dt;
    int i;
    double v;

    ap::ap_error::make_assertion(n>=2, "Spline1DBuildCatmullRom: N<2!");
    ap::ap_error::make_assertion(boundtype==-1||boundtype==0, "Spline1DBuildCatmullRom: incorrect BoundType!");
    
    //
    // Special cases:
    // * N=2, parabolic terminated boundary condition on both ends
    // * N=2, periodic boundary condition
    //
    if( n==2&&boundtype==0 )
    {
        
        //
        // Just linear spline
        //
        spline1dbuildlinear(x, y, n, c);
        return;
    }
    if( n==2&&boundtype==-1 )
    {
        
        //
        // Same as cubic spline with periodic conditions
        //
        spline1dbuildcubic(x, y, n, -1, 0.0, -1, 0.0, c);
        return;
    }
    
    //
    // Periodic or non-periodic boundary conditions
    //
    if( boundtype==-1 )
    {
        
        //
        // Sort points.
        //
        heapsortpoints(x, y, n);
        y(n-1) = y(0);
        
        //
        // Periodic boundary conditions
        //
        d.setlength(n);
        d(0) = (y(1)-y(n-2))/(2*(x(1)-x(0)+x(n-1)-x(n-2)));
        for(i = 1; i <= n-2; i++)
        {
            d(i) = (1-tension)*(y(i+1)-y(i-1))/(x(i+1)-x(i-1));
        }
        d(n-1) = d(0);
        
        //
        // Now problem is reduced to the cubic Hermite spline
        //
        spline1dbuildhermite(x, y, d, n, c);
        c.periodic = true;
    }
    else
    {
        
        //
        // Sort points.
        //
        heapsortpoints(x, y, n);
        
        //
        // Non-periodic boundary conditions
        //
        d.setlength(n);
        for(i = 1; i <= n-2; i++)
        {
            d(i) = (1-tension)*(y(i+1)-y(i-1))/(x(i+1)-x(i-1));
        }
        d(0) = 2*(y(1)-y(0))/(x(1)-x(0))-d(1);
        d(n-1) = 2*(y(n-1)-y(n-2))/(x(n-1)-x(n-2))-d(n-2);
        
        //
        // Now problem is reduced to the cubic Hermite spline
        //
        spline1dbuildhermite(x, y, d, n, c);
    }
}


/*************************************************************************
This subroutine builds Hermite spline interpolant.

INPUT PARAMETERS:
    X           -   spline nodes, array[0..N-1]
    Y           -   function values, array[0..N-1]
    D           -   derivatives, array[0..N-1]
    N           -   points count, N>=2

OUTPUT PARAMETERS:
    C           -   spline interpolant.


ORDER OF POINTS

Subroutine automatically sorts points, so caller may pass unsorted array.

  -- ALGLIB PROJECT --
     Copyright 23.06.2007 by Bochkanov Sergey
*************************************************************************/
void spline1dbuildhermite(ap::real_1d_array x,
     ap::real_1d_array y,
     ap::real_1d_array d,
     int n,
     spline1dinterpolant& c)
{
    int i;
    double delta;
    double delta2;
    double delta3;

    ap::ap_error::make_assertion(n>=2, "BuildHermiteSpline: N<2!");
    
    //
    // Sort points
    //
    heapsortdpoints(x, y, d, n);
    
    //
    // Build
    //
    c.x.setlength(n);
    c.c.setlength(4*(n-1));
    c.periodic = false;
    c.k = 3;
    c.n = n;
    for(i = 0; i <= n-1; i++)
    {
        c.x(i) = x(i);
    }
    for(i = 0; i <= n-2; i++)
    {
        delta = x(i+1)-x(i);
        delta2 = ap::sqr(delta);
        delta3 = delta*delta2;
        c.c(4*i+0) = y(i);
        c.c(4*i+1) = d(i);
        c.c(4*i+2) = (3*(y(i+1)-y(i))-2*d(i)*delta-d(i+1)*delta)/delta2;
        c.c(4*i+3) = (2*(y(i)-y(i+1))+d(i)*delta+d(i+1)*delta)/delta3;
    }
}


/*************************************************************************
This subroutine builds Akima spline interpolant

INPUT PARAMETERS:
    X           -   spline nodes, array[0..N-1]
    Y           -   function values, array[0..N-1]
    N           -   points count, N>=5

OUTPUT PARAMETERS:
    C           -   spline interpolant


ORDER OF POINTS

Subroutine automatically sorts points, so caller may pass unsorted array.

  -- ALGLIB PROJECT --
     Copyright 24.06.2007 by Bochkanov Sergey
*************************************************************************/
void spline1dbuildakima(ap::real_1d_array x,
     ap::real_1d_array y,
     int n,
     spline1dinterpolant& c)
{
    int i;
    ap::real_1d_array d;
    ap::real_1d_array w;
    ap::real_1d_array diff;

    ap::ap_error::make_assertion(n>=5, "BuildAkimaSpline: N<5!");
    
    //
    // Sort points
    //
    heapsortpoints(x, y, n);
    
    //
    // Prepare W (weights), Diff (divided differences)
    //
    w.setlength(n-1);
    diff.setlength(n-1);
    for(i = 0; i <= n-2; i++)
    {
        diff(i) = (y(i+1)-y(i))/(x(i+1)-x(i));
    }
    for(i = 1; i <= n-2; i++)
    {
        w(i) = fabs(diff(i)-diff(i-1));
    }
    
    //
    // Prepare Hermite interpolation scheme
    //
    d.setlength(n);
    for(i = 2; i <= n-3; i++)
    {
        if( ap::fp_neq(fabs(w(i-1))+fabs(w(i+1)),0) )
        {
            d(i) = (w(i+1)*diff(i-1)+w(i-1)*diff(i))/(w(i+1)+w(i-1));
        }
        else
        {
            d(i) = ((x(i+1)-x(i))*diff(i-1)+(x(i)-x(i-1))*diff(i))/(x(i+1)-x(i-1));
        }
    }
    d(0) = diffthreepoint(x(0), x(0), y(0), x(1), y(1), x(2), y(2));
    d(1) = diffthreepoint(x(1), x(0), y(0), x(1), y(1), x(2), y(2));
    d(n-2) = diffthreepoint(x(n-2), x(n-3), y(n-3), x(n-2), y(n-2), x(n-1), y(n-1));
    d(n-1) = diffthreepoint(x(n-1), x(n-3), y(n-3), x(n-2), y(n-2), x(n-1), y(n-1));
    
    //
    // Build Akima spline using Hermite interpolation scheme
    //
    spline1dbuildhermite(x, y, d, n, c);
}


/*************************************************************************
Weighted fitting by cubic  spline,  with constraints on function values or
derivatives.

Equidistant grid with M-2 nodes on [min(x,xc),max(x,xc)] is  used to build
basis functions. Basis functions are cubic splines with continuous  second
derivatives  and  non-fixed first  derivatives  at  interval  ends.  Small
regularizing term is used  when  solving  constrained  tasks  (to  improve
stability).

Task is linear, so linear least squares solver is used. Complexity of this
computational scheme is O(N*M^2), mostly dominated by least squares solver

SEE ALSO
    Spline1DFitHermiteWC()  -   fitting by Hermite splines (more flexible,
                                less smooth)
    Spline1DFitCubic()      -   "lightweight" fitting  by  cubic  splines,
                                without invididual weights and constraints

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    W   -   weights, array[0..N-1]
            Each summand in square  sum  of  approximation deviations from
            given  values  is  multiplied  by  the square of corresponding
            weight. Fill it by 1's if you don't  want  to  solve  weighted
            task.
    N   -   number of points, N>0.
    XC  -   points where spline values/derivatives are constrained,
            array[0..K-1].
    YC  -   values of constraints, array[0..K-1]
    DC  -   array[0..K-1], types of constraints:
            * DC[i]=0   means that S(XC[i])=YC[i]
            * DC[i]=1   means that S'(XC[i])=YC[i]
            SEE BELOW FOR IMPORTANT INFORMATION ON CONSTRAINTS
    K   -   number of constraints, 0<=K<M.
            K=0 means no constraints (XC/YC/DC are not used in such cases)
    M   -   number of basis functions ( = number_of_nodes+2), M>=4.

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearWC() subroutine.
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
                        -3 means inconsistent constraints
                        -1 means another errors in parameters passed
                           (N<=0, for example)
    S   -   spline interpolant.
    Rep -   report, same format as in LSFitLinearWC() subroutine.
            Following fields are set:
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroitine doesn't calculate task's condition number for K<>0.


ORDER OF POINTS

Subroutine automatically sorts points, so caller may pass unsorted array.

SETTING CONSTRAINTS - DANGERS AND OPPORTUNITIES:

Setting constraints can lead  to undesired  results,  like ill-conditioned
behavior, or inconsistency being detected. From the other side,  it allows
us to improve quality of the fit. Here we summarize  our  experience  with
constrained regression splines:
* excessive constraints can be inconsistent. Splines are  piecewise  cubic
  functions, and it is easy to create an example, where  large  number  of
  constraints  concentrated  in  small  area will result in inconsistency.
  Just because spline is not flexible enough to satisfy all of  them.  And
  same constraints spread across the  [min(x),max(x)]  will  be  perfectly
  consistent.
* the more evenly constraints are spread across [min(x),max(x)],  the more
  chances that they will be consistent
* the  greater  is  M (given  fixed  constraints),  the  more chances that
  constraints will be consistent
* in the general case, consistency of constraints IS NOT GUARANTEED.
* in the several special cases, however, we CAN guarantee consistency.
* one of this cases is constraints  on  the  function  values  AND/OR  its
  derivatives at the interval boundaries.
* another  special  case  is ONE constraint on the function value (OR, but
  not AND, derivative) anywhere in the interval

Our final recommendation is to use constraints  WHEN  AND  ONLY  WHEN  you
can't solve your task without them. Anything beyond  special  cases  given
above is not guaranteed and may result in inconsistency.


  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************/
void spline1dfitcubicwc(const ap::real_1d_array& x,
     const ap::real_1d_array& y,
     const ap::real_1d_array& w,
     int n,
     const ap::real_1d_array& xc,
     const ap::real_1d_array& yc,
     const ap::integer_1d_array& dc,
     int k,
     int m,
     int& info,
     spline1dinterpolant& s,
     spline1dfitreport& rep)
{

    spline1dfitinternal(0, x, y, w, n, xc, yc, dc, k, m, info, s, rep);
}


/*************************************************************************
Weighted  fitting  by Hermite spline,  with constraints on function values
or first derivatives.

Equidistant grid with M nodes on [min(x,xc),max(x,xc)] is  used  to  build
basis functions. Basis functions are Hermite splines.  Small  regularizing
term is used when solving constrained tasks (to improve stability).

Task is linear, so linear least squares solver is used. Complexity of this
computational scheme is O(N*M^2), mostly dominated by least squares solver

SEE ALSO
    Spline1DFitCubicWC()    -   fitting by Cubic splines (less flexible,
                                more smooth)
    Spline1DFitHermite()    -   "lightweight" Hermite fitting, without
                                invididual weights and constraints

INPUT PARAMETERS:
    X   -   points, array[0..N-1].
    Y   -   function values, array[0..N-1].
    W   -   weights, array[0..N-1]
            Each summand in square  sum  of  approximation deviations from
            given  values  is  multiplied  by  the square of corresponding
            weight. Fill it by 1's if you don't  want  to  solve  weighted
            task.
    N   -   number of points, N>0.
    XC  -   points where spline values/derivatives are constrained,
            array[0..K-1].
    YC  -   values of constraints, array[0..K-1]
    DC  -   array[0..K-1], types of constraints:
            * DC[i]=0   means that S(XC[i])=YC[i]
            * DC[i]=1   means that S'(XC[i])=YC[i]
            SEE BELOW FOR IMPORTANT INFORMATION ON CONSTRAINTS
    K   -   number of constraints, 0<=K<M.
            K=0 means no constraints (XC/YC/DC are not used in such cases)
    M   -   number of basis functions (= 2 * number of nodes),
            M>=4,
            M IS EVEN!

OUTPUT PARAMETERS:
    Info-   same format as in LSFitLinearW() subroutine:
            * Info>0    task is solved
            * Info<=0   an error occured:
                        -4 means inconvergence of internal SVD
                        -3 means inconsistent constraints
                        -2 means odd M was passed (which is not supported)
                        -1 means another errors in parameters passed
                           (N<=0, for example)
    S   -   spline interpolant.
    Rep -   report, same format as in LSFitLinearW() subroutine.
            Following fields are set:
            * RMSError      rms error on the (X,Y).
            * AvgError      average error on the (X,Y).
            * AvgRelError   average relative error on the non-zero Y
            * MaxError      maximum error
                            NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroitine doesn't calculate task's condition number for K<>0.

IMPORTANT:
    this subroitine supports only even M's


ORDER OF POINTS

Subroutine automatically sorts points, so caller may pass unsorted array.

SETTING CONSTRAINTS - DANGERS AND OPPORTUNITIES:

Setting constraints can lead  to undesired  results,  like ill-conditioned
behavior, or inconsistency being detected. From the other side,  it allows
us to improve quality of the fit. Here we summarize  our  experience  with
constrained regression splines:
* excessive constraints can be inconsistent. Splines are  piecewise  cubic
  functions, and it is easy to create an example, where  large  number  of
  constraints  concentrated  in  small  area will result in inconsistency.
  Just because spline is not flexible enough to satisfy all of  them.  And
  same constraints spread across the  [min(x),max(x)]  will  be  perfectly
  consistent.
* the more evenly constraints are spread across [min(x),max(x)],  the more
  chances that they will be consistent
* the  greater  is  M (given  fixed  constraints),  the  more chances that
  constraints will be consistent
* in the general case, consistency of constraints is NOT GUARANTEED.
* in the several special cases, however, we can guarantee consistency.
* one of this cases is  M>=4  and   constraints  on   the  function  value
  (AND/OR its derivative) at the interval boundaries.
* another special case is M>=4  and  ONE  constraint on the function value
  (OR, BUT NOT AND, derivative) anywhere in [min(x),max(x)]

Our final recommendation is to use constraints  WHEN  AND  ONLY  when  you
can't solve your task without them. Anything beyond  special  cases  given
above is not guaranteed and may result in inconsistency.

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************/
void spline1dfithermitewc(const ap::real_1d_array& x,
     const ap::real_1d_array& y,
     const ap::real_1d_array& w,
     int n,
     const ap::real_1d_array& xc,
     const ap::real_1d_array& yc,
     const ap::integer_1d_array& dc,
     int k,
     int m,
     int& info,
     spline1dinterpolant& s,
     spline1dfitreport& rep)
{

    spline1dfitinternal(1, x, y, w, n, xc, yc, dc, k, m, info, s, rep);
}


/*************************************************************************
Least squares fitting by cubic spline.

This subroutine is "lightweight" alternative for more complex and feature-
rich Spline1DFitCubicWC().  See  Spline1DFitCubicWC() for more information
about subroutine parameters (we don't duplicate it here because of length)

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************/
void spline1dfitcubic(const ap::real_1d_array& x,
     const ap::real_1d_array& y,
     int n,
     int m,
     int& info,
     spline1dinterpolant& s,
     spline1dfitreport& rep)
{
    int i;
    ap::real_1d_array w;
    ap::real_1d_array xc;
    ap::real_1d_array yc;
    ap::integer_1d_array dc;

    if( n>0 )
    {
        w.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            w(i) = 1;
        }
    }
    spline1dfitcubicwc(x, y, w, n, xc, yc, dc, 0, m, info, s, rep);
}


/*************************************************************************
Least squares fitting by Hermite spline.

This subroutine is "lightweight" alternative for more complex and feature-
rich Spline1DFitHermiteWC().  See Spline1DFitHermiteWC()  description  for
more information about subroutine parameters (we don't duplicate  it  here
because of length).

  -- ALGLIB PROJECT --
     Copyright 18.08.2009 by Bochkanov Sergey
*************************************************************************/
void spline1dfithermite(const ap::real_1d_array& x,
     const ap::real_1d_array& y,
     int n,
     int m,
     int& info,
     spline1dinterpolant& s,
     spline1dfitreport& rep)
{
    int i;
    ap::real_1d_array w;
    ap::real_1d_array xc;
    ap::real_1d_array yc;
    ap::integer_1d_array dc;

    if( n>0 )
    {
        w.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            w(i) = 1;
        }
    }
    spline1dfithermitewc(x, y, w, n, xc, yc, dc, 0, m, info, s, rep);
}


/*************************************************************************
This subroutine calculates the value of the spline at the given point X.

INPUT PARAMETERS:
    C   -   spline interpolant
    X   -   point

Result:
    S(x)

  -- ALGLIB PROJECT --
     Copyright 23.06.2007 by Bochkanov Sergey
*************************************************************************/
double spline1dcalc(const spline1dinterpolant& c, double x)
{
    double result;
    int l;
    int r;
    int m;
    double t;

    ap::ap_error::make_assertion(c.k==3, "Spline1DCalc: internal error");
    
    //
    // correct if periodic
    //
    if( c.periodic )
    {
        apperiodicmap(x, c.x(0), c.x(c.n-1), t);
    }
    
    //
    // Binary search in the [ x[0], ..., x[n-2] ] (x[n-1] is not included)
    //
    l = 0;
    r = c.n-2+1;
    while(l!=r-1)
    {
        m = (l+r)/2;
        if( ap::fp_greater_eq(c.x(m),x) )
        {
            r = m;
        }
        else
        {
            l = m;
        }
    }
    
    //
    // Interpolation
    //
    x = x-c.x(l);
    m = 4*l;
    result = c.c(m)+x*(c.c(m+1)+x*(c.c(m+2)+x*c.c(m+3)));
    return result;
}


/*************************************************************************
This subroutine differentiates the spline.

INPUT PARAMETERS:
    C   -   spline interpolant.
    X   -   point

Result:
    S   -   S(x)
    DS  -   S'(x)
    D2S -   S''(x)

  -- ALGLIB PROJECT --
     Copyright 24.06.2007 by Bochkanov Sergey
*************************************************************************/
void spline1ddiff(const spline1dinterpolant& c,
     double x,
     double& s,
     double& ds,
     double& d2s)
{
    int l;
    int r;
    int m;
    double t;

    ap::ap_error::make_assertion(c.k==3, "Spline1DCalc: internal error");
    
    //
    // correct if periodic
    //
    if( c.periodic )
    {
        apperiodicmap(x, c.x(0), c.x(c.n-1), t);
    }
    
    //
    // Binary search
    //
    l = 0;
    r = c.n-2+1;
    while(l!=r-1)
    {
        m = (l+r)/2;
        if( ap::fp_greater_eq(c.x(m),x) )
        {
            r = m;
        }
        else
        {
            l = m;
        }
    }
    
    //
    // Differentiation
    //
    x = x-c.x(l);
    m = 4*l;
    s = c.c(m)+x*(c.c(m+1)+x*(c.c(m+2)+x*c.c(m+3)));
    ds = c.c(m+1)+2*x*c.c(m+2)+3*ap::sqr(x)*c.c(m+3);
    d2s = 2*c.c(m+2)+6*x*c.c(m+3);
}


/*************************************************************************
This subroutine makes the copy of the spline.

INPUT PARAMETERS:
    C   -   spline interpolant.

Result:
    CC  -   spline copy

  -- ALGLIB PROJECT --
     Copyright 29.06.2007 by Bochkanov Sergey
*************************************************************************/
void spline1dcopy(const spline1dinterpolant& c, spline1dinterpolant& cc)
{

    cc.periodic = c.periodic;
    cc.n = c.n;
    cc.k = c.k;
    cc.x.setlength(cc.n);
    ap::vmove(&cc.x(0), 1, &c.x(0), 1, ap::vlen(0,cc.n-1));
    cc.c.setlength((cc.k+1)*(cc.n-1));
    ap::vmove(&cc.c(0), 1, &c.c(0), 1, ap::vlen(0,(cc.k+1)*(cc.n-1)-1));
}


/*************************************************************************
This subroutine unpacks the spline into the coefficients table.

INPUT PARAMETERS:
    C   -   spline interpolant.
    X   -   point

Result:
    Tbl -   coefficients table, unpacked format, array[0..N-2, 0..5].
            For I = 0...N-2:
                Tbl[I,0] = X[i]
                Tbl[I,1] = X[i+1]
                Tbl[I,2] = C0
                Tbl[I,3] = C1
                Tbl[I,4] = C2
                Tbl[I,5] = C3
            On [x[i], x[i+1]] spline is equals to:
                S(x) = C0 + C1*t + C2*t^2 + C3*t^3
                t = x-x[i]

  -- ALGLIB PROJECT --
     Copyright 29.06.2007 by Bochkanov Sergey
*************************************************************************/
void spline1dunpack(const spline1dinterpolant& c,
     int& n,
     ap::real_2d_array& tbl)
{
    int i;
    int j;

    tbl.setbounds(0, c.n-2, 0, 2+c.k);
    n = c.n;
    
    //
    // Fill
    //
    for(i = 0; i <= n-2; i++)
    {
        tbl(i,0) = c.x(i);
        tbl(i,1) = c.x(i+1);
        for(j = 0; j <= c.k; j++)
        {
            tbl(i,2+j) = c.c((c.k+1)*i+j);
        }
    }
}


/*************************************************************************
This subroutine performs linear transformation of the spline argument.

INPUT PARAMETERS:
    C   -   spline interpolant.
    A, B-   transformation coefficients: x = A*t + B
Result:
    C   -   transformed spline

  -- ALGLIB PROJECT --
     Copyright 30.06.2007 by Bochkanov Sergey
*************************************************************************/
void spline1dlintransx(spline1dinterpolant& c, double a, double b)
{
    int i;
    int j;
    int n;
    double v;
    double dv;
    double d2v;
    ap::real_1d_array x;
    ap::real_1d_array y;
    ap::real_1d_array d;

    n = c.n;
    
    //
    // Special case: A=0
    //
    if( ap::fp_eq(a,0) )
    {
        v = spline1dcalc(c, b);
        for(i = 0; i <= n-2; i++)
        {
            c.c((c.k+1)*i) = v;
            for(j = 1; j <= c.k; j++)
            {
                c.c((c.k+1)*i+j) = 0;
            }
        }
        return;
    }
    
    //
    // General case: A<>0.
    // Unpack, X, Y, dY/dX.
    // Scale and pack again.
    //
    ap::ap_error::make_assertion(c.k==3, "Spline1DLinTransX: internal error");
    x.setbounds(0, n-1);
    y.setbounds(0, n-1);
    d.setbounds(0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        x(i) = c.x(i);
        spline1ddiff(c, x(i), v, dv, d2v);
        x(i) = (x(i)-b)/a;
        y(i) = v;
        d(i) = a*dv;
    }
    spline1dbuildhermite(x, y, d, n, c);
}


/*************************************************************************
This subroutine performs linear transformation of the spline.

INPUT PARAMETERS:
    C   -   spline interpolant.
    A, B-   transformation coefficients: S2(x) = A*S(x) + B
Result:
    C   -   transformed spline

  -- ALGLIB PROJECT --
     Copyright 30.06.2007 by Bochkanov Sergey
*************************************************************************/
void spline1dlintransy(spline1dinterpolant& c, double a, double b)
{
    int i;
    int j;
    int n;

    n = c.n;
    for(i = 0; i <= n-2; i++)
    {
        c.c((c.k+1)*i) = a*c.c((c.k+1)*i)+b;
        for(j = 1; j <= c.k; j++)
        {
            c.c((c.k+1)*i+j) = a*c.c((c.k+1)*i+j);
        }
    }
}


/*************************************************************************
This subroutine integrates the spline.

INPUT PARAMETERS:
    C   -   spline interpolant.
    X   -   right bound of the integration interval [a, x],
            here 'a' denotes min(x[])
Result:
    integral(S(t)dt,a,x)

  -- ALGLIB PROJECT --
     Copyright 23.06.2007 by Bochkanov Sergey
*************************************************************************/
double spline1dintegrate(const spline1dinterpolant& c, double x)
{
    double result;
    int n;
    int i;
    int j;
    int l;
    int r;
    int m;
    double w;
    double v;
    double t;
    double intab;
    double additionalterm;

    n = c.n;
    
    //
    // Periodic splines require special treatment. We make
    // following transformation:
    //
    //     integral(S(t)dt,A,X) = integral(S(t)dt,A,Z)+AdditionalTerm
    //
    // here X may lie outside of [A,B], Z lies strictly in [A,B],
    // AdditionalTerm is equals to integral(S(t)dt,A,B) times some
    // integer number (may be zero).
    //
    if( c.periodic&&(ap::fp_less(x,c.x(0))||ap::fp_greater(x,c.x(c.n-1))) )
    {
        
        //
        // compute integral(S(x)dx,A,B)
        //
        intab = 0;
        for(i = 0; i <= c.n-2; i++)
        {
            w = c.x(i+1)-c.x(i);
            m = (c.k+1)*i;
            intab = intab+c.c(m)*w;
            v = w;
            for(j = 1; j <= c.k; j++)
            {
                v = v*w;
                intab = intab+c.c(m+j)*v/(j+1);
            }
        }
        
        //
        // map X into [A,B]
        //
        apperiodicmap(x, c.x(0), c.x(c.n-1), t);
        additionalterm = t*intab;
    }
    else
    {
        additionalterm = 0;
    }
    
    //
    // Binary search in the [ x[0], ..., x[n-2] ] (x[n-1] is not included)
    //
    l = 0;
    r = n-2+1;
    while(l!=r-1)
    {
        m = (l+r)/2;
        if( ap::fp_greater_eq(c.x(m),x) )
        {
            r = m;
        }
        else
        {
            l = m;
        }
    }
    
    //
    // Integration
    //
    result = 0;
    for(i = 0; i <= l-1; i++)
    {
        w = c.x(i+1)-c.x(i);
        m = (c.k+1)*i;
        result = result+c.c(m)*w;
        v = w;
        for(j = 1; j <= c.k; j++)
        {
            v = v*w;
            result = result+c.c(m+j)*v/(j+1);
        }
    }
    w = x-c.x(l);
    m = (c.k+1)*l;
    v = w;
    result = result+c.c(m)*w;
    for(j = 1; j <= c.k; j++)
    {
        v = v*w;
        result = result+c.c(m+j)*v/(j+1);
    }
    result = result+additionalterm;
    return result;
}


/*************************************************************************
Internal spline fitting subroutine

  -- ALGLIB PROJECT --
     Copyright 08.09.2009 by Bochkanov Sergey
*************************************************************************/
static void spline1dfitinternal(int st,
     ap::real_1d_array x,
     ap::real_1d_array y,
     const ap::real_1d_array& w,
     int n,
     ap::real_1d_array xc,
     ap::real_1d_array yc,
     const ap::integer_1d_array& dc,
     int k,
     int m,
     int& info,
     spline1dinterpolant& s,
     spline1dfitreport& rep)
{
    ap::real_2d_array fmatrix;
    ap::real_2d_array cmatrix;
    ap::real_1d_array y2;
    ap::real_1d_array w2;
    ap::real_1d_array sx;
    ap::real_1d_array sy;
    ap::real_1d_array sd;
    ap::real_1d_array tmp;
    ap::real_1d_array xoriginal;
    ap::real_1d_array yoriginal;
    lsfitreport lrep;
    double v0;
    double v1;
    double v2;
    double mx;
    spline1dinterpolant s2;
    int i;
    int j;
    int relcnt;
    double xa;
    double xb;
    double sa;
    double sb;
    double bl;
    double br;
    double decay;

    ap::ap_error::make_assertion(st==0||st==1, "Spline1DFit: internal error!");
    if( st==0&&m<4 )
    {
        info = -1;
        return;
    }
    if( st==1&&m<4 )
    {
        info = -1;
        return;
    }
    if( n<1||k<0||k>=m )
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
    if( st==1&&m%2!=0 )
    {
        
        //
        // Hermite fitter must have even number of basis functions
        //
        info = -2;
        return;
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
    // * SX     -   grid for basis functions
    // * SY     -   values of basis functions at grid points
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
    if( st==0 )
    {
        
        //
        // allocate space for cubic spline
        //
        sx.setlength(m-2);
        sy.setlength(m-2);
        for(j = 0; j <= m-2-1; j++)
        {
            sx(j) = double(2*j)/double(m-2-1)-1;
        }
    }
    if( st==1 )
    {
        
        //
        // allocate space for Hermite spline
        //
        sx.setlength(m/2);
        sy.setlength(m/2);
        sd.setlength(m/2);
        for(j = 0; j <= m/2-1; j++)
        {
            sx(j) = double(2*j)/double(m/2-1)-1;
        }
    }
    
    //
    // Prepare design and constraints matrices:
    // * fill constraints matrix
    // * fill first N rows of design matrix with values
    // * fill next M rows of design matrix with regularizing term
    // * append M zeros to Y
    // * append M elements, mean(abs(W)) each, to W
    //
    for(j = 0; j <= m-1; j++)
    {
        
        //
        // prepare Jth basis function
        //
        if( st==0 )
        {
            
            //
            // cubic spline basis
            //
            for(i = 0; i <= m-2-1; i++)
            {
                sy(i) = 0;
            }
            bl = 0;
            br = 0;
            if( j<m-2 )
            {
                sy(j) = 1;
            }
            if( j==m-2 )
            {
                bl = 1;
            }
            if( j==m-1 )
            {
                br = 1;
            }
            spline1dbuildcubic(sx, sy, m-2, 1, bl, 1, br, s2);
        }
        if( st==1 )
        {
            
            //
            // Hermite basis
            //
            for(i = 0; i <= m/2-1; i++)
            {
                sy(i) = 0;
                sd(i) = 0;
            }
            if( j%2==0 )
            {
                sy(j/2) = 1;
            }
            else
            {
                sd(j/2) = 1;
            }
            spline1dbuildhermite(sx, sy, sd, m/2, s2);
        }
        
        //
        // values at X[], XC[]
        //
        for(i = 0; i <= n-1; i++)
        {
            fmatrix(i,j) = spline1dcalc(s2, x(i));
        }
        for(i = 0; i <= k-1; i++)
        {
            ap::ap_error::make_assertion(dc(i)>=0&&dc(i)<=2, "Spline1DFit: internal error!");
            spline1ddiff(s2, xc(i), v0, v1, v2);
            if( dc(i)==0 )
            {
                cmatrix(i,j) = v0;
            }
            if( dc(i)==1 )
            {
                cmatrix(i,j) = v1;
            }
            if( dc(i)==2 )
            {
                cmatrix(i,j) = v2;
            }
        }
    }
    for(i = 0; i <= k-1; i++)
    {
        cmatrix(i,m) = yc(i);
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
    }
    y2.setlength(n+m);
    w2.setlength(n+m);
    ap::vmove(&y2(0), 1, &y(0), 1, ap::vlen(0,n-1));
    ap::vmove(&w2(0), 1, &w(0), 1, ap::vlen(0,n-1));
    mx = 0;
    for(i = 0; i <= n-1; i++)
    {
        mx = mx+fabs(w(i));
    }
    mx = mx/n;
    for(i = 0; i <= m-1; i++)
    {
        y2(n+i) = 0;
        w2(n+i) = mx;
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
    // Generate spline and scale it
    //
    if( st==0 )
    {
        
        //
        // cubic spline basis
        //
        ap::vmove(&sy(0), 1, &tmp(0), 1, ap::vlen(0,m-2-1));
        spline1dbuildcubic(sx, sy, m-2, 1, tmp(m-2), 1, tmp(m-1), s);
    }
    if( st==1 )
    {
        
        //
        // Hermite basis
        //
        for(i = 0; i <= m/2-1; i++)
        {
            sy(i) = tmp(2*i);
            sd(i) = tmp(2*i+1);
        }
        spline1dbuildhermite(sx, sy, sd, m/2, s);
    }
    spline1dlintransx(s, 2/(xb-xa), -(xa+xb)/(xb-xa));
    spline1dlintransy(s, sb-sa, sa);
    
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
            rep.avgrelerror = rep.avgrelerror+fabs(spline1dcalc(s, xoriginal(i))-yoriginal(i))/fabs(yoriginal(i));
            relcnt = relcnt+1;
        }
    }
    if( relcnt!=0 )
    {
        rep.avgrelerror = rep.avgrelerror/relcnt;
    }
}


/*************************************************************************
Internal subroutine. Heap sort.
*************************************************************************/
static void heapsortpoints(ap::real_1d_array& x, ap::real_1d_array& y, int n)
{
    int i;
    int j;
    int k;
    int t;
    double tmp;
    bool isascending;
    bool isdescending;

    
    //
    // Test for already sorted set
    //
    isascending = true;
    isdescending = true;
    for(i = 1; i <= n-1; i++)
    {
        isascending = isascending&&ap::fp_greater(x(i),x(i-1));
        isdescending = isdescending&&ap::fp_less(x(i),x(i-1));
    }
    if( isascending )
    {
        return;
    }
    if( isdescending )
    {
        for(i = 0; i <= n-1; i++)
        {
            j = n-1-i;
            if( j<=i )
            {
                break;
            }
            tmp = x(i);
            x(i) = x(j);
            x(j) = tmp;
            tmp = y(i);
            y(i) = y(j);
            y(j) = tmp;
        }
        return;
    }
    
    //
    // Special case: N=1
    //
    if( n==1 )
    {
        return;
    }
    
    //
    // General case
    //
    i = 2;
    do
    {
        t = i;
        while(t!=1)
        {
            k = t/2;
            if( ap::fp_greater_eq(x(k-1),x(t-1)) )
            {
                t = 1;
            }
            else
            {
                tmp = x(k-1);
                x(k-1) = x(t-1);
                x(t-1) = tmp;
                tmp = y(k-1);
                y(k-1) = y(t-1);
                y(t-1) = tmp;
                t = k;
            }
        }
        i = i+1;
    }
    while(i<=n);
    i = n-1;
    do
    {
        tmp = x(i);
        x(i) = x(0);
        x(0) = tmp;
        tmp = y(i);
        y(i) = y(0);
        y(0) = tmp;
        t = 1;
        while(t!=0)
        {
            k = 2*t;
            if( k>i )
            {
                t = 0;
            }
            else
            {
                if( k<i )
                {
                    if( ap::fp_greater(x(k),x(k-1)) )
                    {
                        k = k+1;
                    }
                }
                if( ap::fp_greater_eq(x(t-1),x(k-1)) )
                {
                    t = 0;
                }
                else
                {
                    tmp = x(k-1);
                    x(k-1) = x(t-1);
                    x(t-1) = tmp;
                    tmp = y(k-1);
                    y(k-1) = y(t-1);
                    y(t-1) = tmp;
                    t = k;
                }
            }
        }
        i = i-1;
    }
    while(i>=1);
}


/*************************************************************************
Internal subroutine. Heap sort.
*************************************************************************/
static void heapsortdpoints(ap::real_1d_array& x,
     ap::real_1d_array& y,
     ap::real_1d_array& d,
     int n)
{
    int i;
    int j;
    int k;
    int t;
    double tmp;
    bool isascending;
    bool isdescending;

    
    //
    // Test for already sorted set
    //
    isascending = true;
    isdescending = true;
    for(i = 1; i <= n-1; i++)
    {
        isascending = isascending&&ap::fp_greater(x(i),x(i-1));
        isdescending = isdescending&&ap::fp_less(x(i),x(i-1));
    }
    if( isascending )
    {
        return;
    }
    if( isdescending )
    {
        for(i = 0; i <= n-1; i++)
        {
            j = n-1-i;
            if( j<=i )
            {
                break;
            }
            tmp = x(i);
            x(i) = x(j);
            x(j) = tmp;
            tmp = y(i);
            y(i) = y(j);
            y(j) = tmp;
            tmp = d(i);
            d(i) = d(j);
            d(j) = tmp;
        }
        return;
    }
    
    //
    // Special case: N=1
    //
    if( n==1 )
    {
        return;
    }
    
    //
    // General case
    //
    i = 2;
    do
    {
        t = i;
        while(t!=1)
        {
            k = t/2;
            if( ap::fp_greater_eq(x(k-1),x(t-1)) )
            {
                t = 1;
            }
            else
            {
                tmp = x(k-1);
                x(k-1) = x(t-1);
                x(t-1) = tmp;
                tmp = y(k-1);
                y(k-1) = y(t-1);
                y(t-1) = tmp;
                tmp = d(k-1);
                d(k-1) = d(t-1);
                d(t-1) = tmp;
                t = k;
            }
        }
        i = i+1;
    }
    while(i<=n);
    i = n-1;
    do
    {
        tmp = x(i);
        x(i) = x(0);
        x(0) = tmp;
        tmp = y(i);
        y(i) = y(0);
        y(0) = tmp;
        tmp = d(i);
        d(i) = d(0);
        d(0) = tmp;
        t = 1;
        while(t!=0)
        {
            k = 2*t;
            if( k>i )
            {
                t = 0;
            }
            else
            {
                if( k<i )
                {
                    if( ap::fp_greater(x(k),x(k-1)) )
                    {
                        k = k+1;
                    }
                }
                if( ap::fp_greater_eq(x(t-1),x(k-1)) )
                {
                    t = 0;
                }
                else
                {
                    tmp = x(k-1);
                    x(k-1) = x(t-1);
                    x(t-1) = tmp;
                    tmp = y(k-1);
                    y(k-1) = y(t-1);
                    y(t-1) = tmp;
                    tmp = d(k-1);
                    d(k-1) = d(t-1);
                    d(t-1) = tmp;
                    t = k;
                }
            }
        }
        i = i-1;
    }
    while(i>=1);
}


/*************************************************************************
Internal subroutine. Tridiagonal solver. Solves

( B[0] C[0]                      )
( A[1] B[1] C[1]                 )
(      A[2] B[2] C[2]            )
(            ..........          ) * X = D
(            ..........          )
(           A[N-2] B[N-2] C[N-2] )
(                  A[N-1] B[N-1] )

*************************************************************************/
static void solvetridiagonal(ap::real_1d_array a,
     ap::real_1d_array b,
     ap::real_1d_array c,
     ap::real_1d_array d,
     int n,
     ap::real_1d_array& x)
{
    int k;
    double t;

    x.setbounds(0, n-1);
    a(0) = 0;
    c(n-1) = 0;
    for(k = 1; k <= n-1; k++)
    {
        t = a(k)/b(k-1);
        b(k) = b(k)-t*c(k-1);
        d(k) = d(k)-t*d(k-1);
    }
    x(n-1) = d(n-1)/b(n-1);
    for(k = n-2; k >= 0; k--)
    {
        x(k) = (d(k)-c(k)*x(k+1))/b(k);
    }
}


/*************************************************************************
Internal subroutine. Cyclic tridiagonal solver. Solves

( B[0] C[0]                 A[0] )
( A[1] B[1] C[1]                 )
(      A[2] B[2] C[2]            )
(            ..........          ) * X = D
(            ..........          )
(           A[N-2] B[N-2] C[N-2] )
( C[N-1]           A[N-1] B[N-1] )
*************************************************************************/
static void solvecyclictridiagonal(const ap::real_1d_array& a,
     ap::real_1d_array b,
     const ap::real_1d_array& c,
     const ap::real_1d_array& d,
     int n,
     ap::real_1d_array& x)
{
    int k;
    double t;
    double alpha;
    double beta;
    double gamma;
    ap::real_1d_array y;
    ap::real_1d_array z;
    ap::real_1d_array u;

    beta = a(0);
    alpha = c(n-1);
    gamma = -b(0);
    b(0) = 2*b(0);
    b(n-1) = b(n-1)-alpha*beta/gamma;
    u.setlength(n);
    for(k = 0; k <= n-1; k++)
    {
        u(k) = 0;
    }
    u(0) = gamma;
    u(n-1) = alpha;
    solvetridiagonal(a, b, c, d, n, y);
    solvetridiagonal(a, b, c, u, n, z);
    x.setlength(n);
    for(k = 0; k <= n-1; k++)
    {
        x(k) = y(k)-(y(0)+beta/gamma*y(n-1))/(1+z(0)+beta/gamma*z(n-1))*z(k);
    }
}


/*************************************************************************
Internal subroutine. Three-point differentiation
*************************************************************************/
static double diffthreepoint(double t,
     double x0,
     double f0,
     double x1,
     double f1,
     double x2,
     double f2)
{
    double result;
    double a;
    double b;

    t = t-x0;
    x1 = x1-x0;
    x2 = x2-x0;
    a = (f2-f0-x2/x1*(f1-f0))/(ap::sqr(x2)-x1*x2);
    b = (f1-f0-a*ap::sqr(x1))/x1;
    result = 2*a*t+b;
    return result;
}




