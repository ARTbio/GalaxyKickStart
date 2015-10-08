/*************************************************************************
Copyright (c) 2006-2010, Sergey Bochkanov (ALGLIB project).

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
#include "pspline.h"

static void pspline2par(const ap::real_2d_array& xy,
     int n,
     int pt,
     ap::real_1d_array& p);
static void pspline3par(const ap::real_2d_array& xy,
     int n,
     int pt,
     ap::real_1d_array& p);

/*************************************************************************
This function  builds  non-periodic 2-dimensional parametric spline  which
starts at (X[0],Y[0]) and ends at (X[N-1],Y[N-1]).

INPUT PARAMETERS:
    XY  -   points, array[0..N-1,0..1].
            XY[I,0:1] corresponds to the Ith point.
            Order of points is important!
    N   -   points count, N>=5 for Akima splines, N>=2 for other types  of
            splines.
    ST  -   spline type:
            * 0     Akima spline
            * 1     parabolically terminated Catmull-Rom spline (Tension=0)
            * 2     parabolically terminated cubic spline
    PT  -   parameterization type:
            * 0     uniform
            * 1     chord length
            * 2     centripetal

OUTPUT PARAMETERS:
    P   -   parametric spline interpolant


NOTES:
* this function  assumes  that  there all consequent points  are distinct.
  I.e. (x0,y0)<>(x1,y1),  (x1,y1)<>(x2,y2),  (x2,y2)<>(x3,y3)  and  so on.
  However, non-consequent points may coincide, i.e. we can  have  (x0,y0)=
  =(x2,y2).

  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline2build(ap::real_2d_array xy,
     int n,
     int st,
     int pt,
     pspline2interpolant& p)
{
    ap::real_1d_array tmp;
    double v;
    int i;

    ap::ap_error::make_assertion(st>=0&&st<=2, "PSpline2Build: incorrect spline type!");
    ap::ap_error::make_assertion(pt>=0&&pt<=2, "PSpline2Build: incorrect parameterization type!");
    if( st==0 )
    {
        ap::ap_error::make_assertion(n>=5, "PSpline2Build: N<5 (minimum value for Akima splines)!");
    }
    else
    {
        ap::ap_error::make_assertion(n>=2, "PSpline2Build: N<2!");
    }
    
    //
    // Prepare
    //
    p.n = n;
    p.periodic = false;
    tmp.setlength(n);
    
    //
    // Build parameterization, check that all parameters are distinct
    //
    pspline2par(xy, n, pt, p.p);
    ap::ap_error::make_assertion(apservaredistinct(p.p, n), "PSpline2Build: consequent points are too close!");
    
    //
    // Build splines
    //
    if( st==0 )
    {
        ap::vmove(&tmp(0), 1, &xy(0, 0), xy.getstride(), ap::vlen(0,n-1));
        spline1dbuildakima(p.p, tmp, n, p.x);
        ap::vmove(&tmp(0), 1, &xy(0, 1), xy.getstride(), ap::vlen(0,n-1));
        spline1dbuildakima(p.p, tmp, n, p.y);
    }
    if( st==1 )
    {
        ap::vmove(&tmp(0), 1, &xy(0, 0), xy.getstride(), ap::vlen(0,n-1));
        spline1dbuildcatmullrom(p.p, tmp, n, 0, 0.0, p.x);
        ap::vmove(&tmp(0), 1, &xy(0, 1), xy.getstride(), ap::vlen(0,n-1));
        spline1dbuildcatmullrom(p.p, tmp, n, 0, 0.0, p.y);
    }
    if( st==2 )
    {
        ap::vmove(&tmp(0), 1, &xy(0, 0), xy.getstride(), ap::vlen(0,n-1));
        spline1dbuildcubic(p.p, tmp, n, 0, 0.0, 0, 0.0, p.x);
        ap::vmove(&tmp(0), 1, &xy(0, 1), xy.getstride(), ap::vlen(0,n-1));
        spline1dbuildcubic(p.p, tmp, n, 0, 0.0, 0, 0.0, p.y);
    }
}


/*************************************************************************
This function  builds  non-periodic 3-dimensional parametric spline  which
starts at (X[0],Y[0],Z[0]) and ends at (X[N-1],Y[N-1],Z[N-1]).

Same as PSpline2Build() function, but for 3D, so we  won't  duplicate  its
description here.

  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline3build(ap::real_2d_array xy,
     int n,
     int st,
     int pt,
     pspline3interpolant& p)
{
    ap::real_1d_array tmp;
    double v;
    int i;

    ap::ap_error::make_assertion(st>=0&&st<=2, "PSpline3Build: incorrect spline type!");
    ap::ap_error::make_assertion(pt>=0&&pt<=2, "PSpline3Build: incorrect parameterization type!");
    if( st==0 )
    {
        ap::ap_error::make_assertion(n>=5, "PSpline3Build: N<5 (minimum value for Akima splines)!");
    }
    else
    {
        ap::ap_error::make_assertion(n>=2, "PSpline3Build: N<2!");
    }
    
    //
    // Prepare
    //
    p.n = n;
    p.periodic = false;
    tmp.setlength(n);
    
    //
    // Build parameterization, check that all parameters are distinct
    //
    pspline3par(xy, n, pt, p.p);
    ap::ap_error::make_assertion(apservaredistinct(p.p, n), "PSpline3Build: consequent points are too close!");
    
    //
    // Build splines
    //
    if( st==0 )
    {
        ap::vmove(&tmp(0), 1, &xy(0, 0), xy.getstride(), ap::vlen(0,n-1));
        spline1dbuildakima(p.p, tmp, n, p.x);
        ap::vmove(&tmp(0), 1, &xy(0, 1), xy.getstride(), ap::vlen(0,n-1));
        spline1dbuildakima(p.p, tmp, n, p.y);
        ap::vmove(&tmp(0), 1, &xy(0, 2), xy.getstride(), ap::vlen(0,n-1));
        spline1dbuildakima(p.p, tmp, n, p.z);
    }
    if( st==1 )
    {
        ap::vmove(&tmp(0), 1, &xy(0, 0), xy.getstride(), ap::vlen(0,n-1));
        spline1dbuildcatmullrom(p.p, tmp, n, 0, 0.0, p.x);
        ap::vmove(&tmp(0), 1, &xy(0, 1), xy.getstride(), ap::vlen(0,n-1));
        spline1dbuildcatmullrom(p.p, tmp, n, 0, 0.0, p.y);
        ap::vmove(&tmp(0), 1, &xy(0, 2), xy.getstride(), ap::vlen(0,n-1));
        spline1dbuildcatmullrom(p.p, tmp, n, 0, 0.0, p.z);
    }
    if( st==2 )
    {
        ap::vmove(&tmp(0), 1, &xy(0, 0), xy.getstride(), ap::vlen(0,n-1));
        spline1dbuildcubic(p.p, tmp, n, 0, 0.0, 0, 0.0, p.x);
        ap::vmove(&tmp(0), 1, &xy(0, 1), xy.getstride(), ap::vlen(0,n-1));
        spline1dbuildcubic(p.p, tmp, n, 0, 0.0, 0, 0.0, p.y);
        ap::vmove(&tmp(0), 1, &xy(0, 2), xy.getstride(), ap::vlen(0,n-1));
        spline1dbuildcubic(p.p, tmp, n, 0, 0.0, 0, 0.0, p.z);
    }
}


/*************************************************************************
This  function  builds  periodic  2-dimensional  parametric  spline  which
starts at (X[0],Y[0]), goes through all points to (X[N-1],Y[N-1]) and then
back to (X[0],Y[0]).

INPUT PARAMETERS:
    XY  -   points, array[0..N-1,0..1].
            XY[I,0:1] corresponds to the Ith point.
            XY[N-1,0:1] must be different from XY[0,0:1].
            Order of points is important!
    N   -   points count, N>=3 for other types of splines.
    ST  -   spline type:
            * 1     Catmull-Rom spline (Tension=0) with cyclic boundary conditions
            * 2     cubic spline with cyclic boundary conditions
    PT  -   parameterization type:
            * 0     uniform
            * 1     chord length
            * 2     centripetal

OUTPUT PARAMETERS:
    P   -   parametric spline interpolant


NOTES:
* this function  assumes  that there all consequent points  are  distinct.
  I.e. (x0,y0)<>(x1,y1), (x1,y1)<>(x2,y2),  (x2,y2)<>(x3,y3)  and  so  on.
  However, non-consequent points may coincide, i.e. we can  have  (x0,y0)=
  =(x2,y2).
* last point of sequence is NOT equal to the first  point.  You  shouldn't
  make curve "explicitly periodic" by making them equal.

  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline2buildperiodic(ap::real_2d_array xy,
     int n,
     int st,
     int pt,
     pspline2interpolant& p)
{
    ap::real_2d_array xyp;
    ap::real_1d_array tmp;
    double v;
    int i;

    ap::ap_error::make_assertion(st>=1&&st<=2, "PSpline2BuildPeriodic: incorrect spline type!");
    ap::ap_error::make_assertion(pt>=0&&pt<=2, "PSpline2BuildPeriodic: incorrect parameterization type!");
    ap::ap_error::make_assertion(n>=3, "PSpline2BuildPeriodic: N<3!");
    
    //
    // Prepare
    //
    p.n = n;
    p.periodic = true;
    tmp.setlength(n+1);
    xyp.setlength(n+1, 2);
    ap::vmove(&xyp(0, 0), xyp.getstride(), &xy(0, 0), xy.getstride(), ap::vlen(0,n-1));
    ap::vmove(&xyp(0, 1), xyp.getstride(), &xy(0, 1), xy.getstride(), ap::vlen(0,n-1));
    ap::vmove(&xyp(n, 0), 1, &xy(0, 0), 1, ap::vlen(0,1));
    
    //
    // Build parameterization, check that all parameters are distinct
    //
    pspline2par(xyp, n+1, pt, p.p);
    ap::ap_error::make_assertion(apservaredistinct(p.p, n+1), "PSpline2BuildPeriodic: consequent (or first and last) points are too close!");
    
    //
    // Build splines
    //
    if( st==1 )
    {
        ap::vmove(&tmp(0), 1, &xyp(0, 0), xyp.getstride(), ap::vlen(0,n));
        spline1dbuildcatmullrom(p.p, tmp, n+1, -1, 0.0, p.x);
        ap::vmove(&tmp(0), 1, &xyp(0, 1), xyp.getstride(), ap::vlen(0,n));
        spline1dbuildcatmullrom(p.p, tmp, n+1, -1, 0.0, p.y);
    }
    if( st==2 )
    {
        ap::vmove(&tmp(0), 1, &xyp(0, 0), xyp.getstride(), ap::vlen(0,n));
        spline1dbuildcubic(p.p, tmp, n+1, -1, 0.0, -1, 0.0, p.x);
        ap::vmove(&tmp(0), 1, &xyp(0, 1), xyp.getstride(), ap::vlen(0,n));
        spline1dbuildcubic(p.p, tmp, n+1, -1, 0.0, -1, 0.0, p.y);
    }
}


/*************************************************************************
This  function  builds  periodic  3-dimensional  parametric  spline  which
starts at (X[0],Y[0],Z[0]), goes through all points to (X[N-1],Y[N-1],Z[N-1])
and then back to (X[0],Y[0],Z[0]).

Same as PSpline2Build() function, but for 3D, so we  won't  duplicate  its
description here.

  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline3buildperiodic(ap::real_2d_array xy,
     int n,
     int st,
     int pt,
     pspline3interpolant& p)
{
    ap::real_2d_array xyp;
    ap::real_1d_array tmp;
    double v;
    int i;

    ap::ap_error::make_assertion(st>=1&&st<=2, "PSpline3BuildPeriodic: incorrect spline type!");
    ap::ap_error::make_assertion(pt>=0&&pt<=2, "PSpline3BuildPeriodic: incorrect parameterization type!");
    ap::ap_error::make_assertion(n>=3, "PSpline3BuildPeriodic: N<3!");
    
    //
    // Prepare
    //
    p.n = n;
    p.periodic = true;
    tmp.setlength(n+1);
    xyp.setlength(n+1, 3);
    ap::vmove(&xyp(0, 0), xyp.getstride(), &xy(0, 0), xy.getstride(), ap::vlen(0,n-1));
    ap::vmove(&xyp(0, 1), xyp.getstride(), &xy(0, 1), xy.getstride(), ap::vlen(0,n-1));
    ap::vmove(&xyp(0, 2), xyp.getstride(), &xy(0, 2), xy.getstride(), ap::vlen(0,n-1));
    ap::vmove(&xyp(n, 0), 1, &xy(0, 0), 1, ap::vlen(0,2));
    
    //
    // Build parameterization, check that all parameters are distinct
    //
    pspline3par(xyp, n+1, pt, p.p);
    ap::ap_error::make_assertion(apservaredistinct(p.p, n+1), "PSplineBuild2Periodic: consequent (or first and last) points are too close!");
    
    //
    // Build splines
    //
    if( st==1 )
    {
        ap::vmove(&tmp(0), 1, &xyp(0, 0), xyp.getstride(), ap::vlen(0,n));
        spline1dbuildcatmullrom(p.p, tmp, n+1, -1, 0.0, p.x);
        ap::vmove(&tmp(0), 1, &xyp(0, 1), xyp.getstride(), ap::vlen(0,n));
        spline1dbuildcatmullrom(p.p, tmp, n+1, -1, 0.0, p.y);
        ap::vmove(&tmp(0), 1, &xyp(0, 2), xyp.getstride(), ap::vlen(0,n));
        spline1dbuildcatmullrom(p.p, tmp, n+1, -1, 0.0, p.z);
    }
    if( st==2 )
    {
        ap::vmove(&tmp(0), 1, &xyp(0, 0), xyp.getstride(), ap::vlen(0,n));
        spline1dbuildcubic(p.p, tmp, n+1, -1, 0.0, -1, 0.0, p.x);
        ap::vmove(&tmp(0), 1, &xyp(0, 1), xyp.getstride(), ap::vlen(0,n));
        spline1dbuildcubic(p.p, tmp, n+1, -1, 0.0, -1, 0.0, p.y);
        ap::vmove(&tmp(0), 1, &xyp(0, 2), xyp.getstride(), ap::vlen(0,n));
        spline1dbuildcubic(p.p, tmp, n+1, -1, 0.0, -1, 0.0, p.z);
    }
}


/*************************************************************************
This function returns vector of parameter values correspoding to points.

I.e. for P created from (X[0],Y[0])...(X[N-1],Y[N-1]) and U=TValues(P)  we
have
    (X[0],Y[0]) = PSpline2Calc(P,U[0]),
    (X[1],Y[1]) = PSpline2Calc(P,U[1]),
    (X[2],Y[2]) = PSpline2Calc(P,U[2]),
    ...

INPUT PARAMETERS:
    P   -   parametric spline interpolant

OUTPUT PARAMETERS:
    N   -   array size
    T   -   array[0..N-1]


NOTES:
* for non-periodic splines U[0]=0, U[0]<U[1]<...<U[N-1], U[N-1]=1
* for periodic splines     U[0]=0, U[0]<U[1]<...<U[N-1], U[N-1]<1

  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline2parametervalues(const pspline2interpolant& p,
     int& n,
     ap::real_1d_array& t)
{

    ap::ap_error::make_assertion(p.n>=2, "PSpline2ParameterValues: internal error!");
    n = p.n;
    t.setlength(n);
    ap::vmove(&t(0), 1, &p.p(0), 1, ap::vlen(0,n-1));
    t(0) = 0;
    if( !p.periodic )
    {
        t(n-1) = 1;
    }
}


/*************************************************************************
This function returns vector of parameter values correspoding to points.

Same as PSpline2ParameterValues(), but for 3D.

  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline3parametervalues(const pspline3interpolant& p,
     int& n,
     ap::real_1d_array& t)
{

    ap::ap_error::make_assertion(p.n>=2, "PSpline3ParameterValues: internal error!");
    n = p.n;
    t.setlength(n);
    ap::vmove(&t(0), 1, &p.p(0), 1, ap::vlen(0,n-1));
    t(0) = 0;
    if( !p.periodic )
    {
        t(n-1) = 1;
    }
}


/*************************************************************************
This function  calculates  the value of the parametric spline for a  given
value of parameter T

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X   -   X-position
    Y   -   Y-position


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline2calc(const pspline2interpolant& p,
     double t,
     double& x,
     double& y)
{

    if( p.periodic )
    {
        t = t-ap::ifloor(t);
    }
    x = spline1dcalc(p.x, t);
    y = spline1dcalc(p.y, t);
}


/*************************************************************************
This function  calculates  the value of the parametric spline for a  given
value of parameter T.

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X   -   X-position
    Y   -   Y-position
    Z   -   Z-position


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline3calc(const pspline3interpolant& p,
     double t,
     double& x,
     double& y,
     double& z)
{

    if( p.periodic )
    {
        t = t-ap::ifloor(t);
    }
    x = spline1dcalc(p.x, t);
    y = spline1dcalc(p.y, t);
    z = spline1dcalc(p.z, t);
}


/*************************************************************************
This function  calculates  tangent vector for a given value of parameter T

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X    -   X-component of tangent vector (normalized)
    Y    -   Y-component of tangent vector (normalized)
    
NOTE:
    X^2+Y^2 is either 1 (for non-zero tangent vector) or 0.


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline2tangent(const pspline2interpolant& p,
     double t,
     double& x,
     double& y)
{
    double v;
    double v0;
    double v1;

    if( p.periodic )
    {
        t = t-ap::ifloor(t);
    }
    pspline2diff(p, t, v0, x, v1, y);
    if( ap::fp_neq(x,0)||ap::fp_neq(y,0) )
    {
        
        //
        // this code is a bit more complex than X^2+Y^2 to avoid
        // overflow for large values of X and Y.
        //
        v = safepythag2(x, y);
        x = x/v;
        y = y/v;
    }
}


/*************************************************************************
This function  calculates  tangent vector for a given value of parameter T

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X    -   X-component of tangent vector (normalized)
    Y    -   Y-component of tangent vector (normalized)
    Z    -   Z-component of tangent vector (normalized)

NOTE:
    X^2+Y^2+Z^2 is either 1 (for non-zero tangent vector) or 0.


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline3tangent(const pspline3interpolant& p,
     double t,
     double& x,
     double& y,
     double& z)
{
    double v;
    double v0;
    double v1;
    double v2;

    if( p.periodic )
    {
        t = t-ap::ifloor(t);
    }
    pspline3diff(p, t, v0, x, v1, y, v2, z);
    if( ap::fp_neq(x,0)||ap::fp_neq(y,0)||ap::fp_neq(z,0) )
    {
        v = safepythag3(x, y, z);
        x = x/v;
        y = y/v;
        z = z/v;
    }
}


/*************************************************************************
This function calculates derivative, i.e. it returns (dX/dT,dY/dT).

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X   -   X-value
    DX  -   X-derivative
    Y   -   Y-value
    DY  -   Y-derivative


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline2diff(const pspline2interpolant& p,
     double t,
     double& x,
     double& dx,
     double& y,
     double& dy)
{
    double d2s;

    if( p.periodic )
    {
        t = t-ap::ifloor(t);
    }
    spline1ddiff(p.x, t, x, dx, d2s);
    spline1ddiff(p.y, t, y, dy, d2s);
}


/*************************************************************************
This function calculates derivative, i.e. it returns (dX/dT,dY/dT,dZ/dT).

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X   -   X-value
    DX  -   X-derivative
    Y   -   Y-value
    DY  -   Y-derivative
    Z   -   Z-value
    DZ  -   Z-derivative


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline3diff(const pspline3interpolant& p,
     double t,
     double& x,
     double& dx,
     double& y,
     double& dy,
     double& z,
     double& dz)
{
    double d2s;

    if( p.periodic )
    {
        t = t-ap::ifloor(t);
    }
    spline1ddiff(p.x, t, x, dx, d2s);
    spline1ddiff(p.y, t, y, dy, d2s);
    spline1ddiff(p.z, t, z, dz, d2s);
}


/*************************************************************************
This function calculates first and second derivative with respect to T.

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X   -   X-value
    DX  -   derivative
    D2X -   second derivative
    Y   -   Y-value
    DY  -   derivative
    D2Y -   second derivative


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline2diff2(const pspline2interpolant& p,
     double t,
     double& x,
     double& dx,
     double& d2x,
     double& y,
     double& dy,
     double& d2y)
{

    if( p.periodic )
    {
        t = t-ap::ifloor(t);
    }
    spline1ddiff(p.x, t, x, dx, d2x);
    spline1ddiff(p.y, t, y, dy, d2y);
}


/*************************************************************************
This function calculates first and second derivative with respect to T.

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    T   -   point:
            * T in [0,1] corresponds to interval spanned by points
            * for non-periodic splines T<0 (or T>1) correspond to parts of
              the curve before the first (after the last) point
            * for periodic splines T<0 (or T>1) are projected  into  [0,1]
              by making T=T-floor(T).

OUTPUT PARAMETERS:
    X   -   X-value
    DX  -   derivative
    D2X -   second derivative
    Y   -   Y-value
    DY  -   derivative
    D2Y -   second derivative
    Z   -   Z-value
    DZ  -   derivative
    D2Z -   second derivative


  -- ALGLIB PROJECT --
     Copyright 28.05.2010 by Bochkanov Sergey
*************************************************************************/
void pspline3diff2(const pspline3interpolant& p,
     double t,
     double& x,
     double& dx,
     double& d2x,
     double& y,
     double& dy,
     double& d2y,
     double& z,
     double& dz,
     double& d2z)
{

    if( p.periodic )
    {
        t = t-ap::ifloor(t);
    }
    spline1ddiff(p.x, t, x, dx, d2x);
    spline1ddiff(p.y, t, y, dy, d2y);
    spline1ddiff(p.z, t, z, dz, d2z);
}


/*************************************************************************
This function  calculates  arc length, i.e. length of  curve  between  t=a
and t=b.

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    A,B -   parameter values corresponding to arc ends:
            * B>A will result in positive length returned
            * B<A will result in negative length returned

RESULT:
    length of arc starting at T=A and ending at T=B.


  -- ALGLIB PROJECT --
     Copyright 30.05.2010 by Bochkanov Sergey
*************************************************************************/
double pspline2arclength(const pspline2interpolant& p, double a, double b)
{
    double result;
    autogkstate state;
    autogkreport rep;
    double sx;
    double dsx;
    double d2sx;
    double sy;
    double dsy;
    double d2sy;

    autogksmooth(a, b, state);
    while(autogkiteration(state))
    {
        spline1ddiff(p.x, state.x, sx, dsx, d2sx);
        spline1ddiff(p.y, state.x, sy, dsy, d2sy);
        state.f = safepythag2(dsx, dsy);
    }
    autogkresults(state, result, rep);
    ap::ap_error::make_assertion(rep.terminationtype>0, "PSpline2ArcLength: internal error!");
    return result;
}


/*************************************************************************
This function  calculates  arc length, i.e. length of  curve  between  t=a
and t=b.

INPUT PARAMETERS:
    P   -   parametric spline interpolant
    A,B -   parameter values corresponding to arc ends:
            * B>A will result in positive length returned
            * B<A will result in negative length returned

RESULT:
    length of arc starting at T=A and ending at T=B.


  -- ALGLIB PROJECT --
     Copyright 30.05.2010 by Bochkanov Sergey
*************************************************************************/
double pspline3arclength(const pspline3interpolant& p, double a, double b)
{
    double result;
    autogkstate state;
    autogkreport rep;
    double sx;
    double dsx;
    double d2sx;
    double sy;
    double dsy;
    double d2sy;
    double sz;
    double dsz;
    double d2sz;

    autogksmooth(a, b, state);
    while(autogkiteration(state))
    {
        spline1ddiff(p.x, state.x, sx, dsx, d2sx);
        spline1ddiff(p.y, state.x, sy, dsy, d2sy);
        spline1ddiff(p.z, state.x, sz, dsz, d2sz);
        state.f = safepythag3(dsx, dsy, dsz);
    }
    autogkresults(state, result, rep);
    ap::ap_error::make_assertion(rep.terminationtype>0, "PSpline3ArcLength: internal error!");
    return result;
}


/*************************************************************************
Builds non-periodic parameterization for 2-dimensional spline
*************************************************************************/
static void pspline2par(const ap::real_2d_array& xy,
     int n,
     int pt,
     ap::real_1d_array& p)
{
    double v;
    int i;

    ap::ap_error::make_assertion(pt>=0&&pt<=2, "PSpline2Par: internal error!");
    
    //
    // Build parameterization:
    // * fill by non-normalized values
    // * normalize them so we have P[0]=0, P[N-1]=1.
    //
    p.setlength(n);
    if( pt==0 )
    {
        for(i = 0; i <= n-1; i++)
        {
            p(i) = i;
        }
    }
    if( pt==1 )
    {
        p(0) = 0;
        for(i = 1; i <= n-1; i++)
        {
            p(i) = p(i-1)+safepythag2(xy(i,0)-xy(i-1,0), xy(i,1)-xy(i-1,1));
        }
    }
    if( pt==2 )
    {
        p(0) = 0;
        for(i = 1; i <= n-1; i++)
        {
            p(i) = p(i-1)+sqrt(safepythag2(xy(i,0)-xy(i-1,0), xy(i,1)-xy(i-1,1)));
        }
    }
    v = 1/p(n-1);
    ap::vmul(&p(0), 1, ap::vlen(0,n-1), v);
}


/*************************************************************************
Builds non-periodic parameterization for 3-dimensional spline
*************************************************************************/
static void pspline3par(const ap::real_2d_array& xy,
     int n,
     int pt,
     ap::real_1d_array& p)
{
    double v;
    int i;

    ap::ap_error::make_assertion(pt>=0&&pt<=2, "PSpline3Par: internal error!");
    
    //
    // Build parameterization:
    // * fill by non-normalized values
    // * normalize them so we have P[0]=0, P[N-1]=1.
    //
    p.setlength(n);
    if( pt==0 )
    {
        for(i = 0; i <= n-1; i++)
        {
            p(i) = i;
        }
    }
    if( pt==1 )
    {
        p(0) = 0;
        for(i = 1; i <= n-1; i++)
        {
            p(i) = p(i-1)+safepythag3(xy(i,0)-xy(i-1,0), xy(i,1)-xy(i-1,1), xy(i,2)-xy(i-1,2));
        }
    }
    if( pt==2 )
    {
        p(0) = 0;
        for(i = 1; i <= n-1; i++)
        {
            p(i) = p(i-1)+sqrt(safepythag3(xy(i,0)-xy(i-1,0), xy(i,1)-xy(i-1,1), xy(i,2)-xy(i-1,2)));
        }
    }
    v = 1/p(n-1);
    ap::vmul(&p(0), 1, ap::vlen(0,n-1), v);
}




