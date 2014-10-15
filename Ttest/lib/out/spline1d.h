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

#ifndef _spline1d_h
#define _spline1d_h

#include "ap.h"
#include "ialglib.h"

#include "spline3.h"
#include "blas.h"
#include "reflections.h"
#include "creflections.h"
#include "hqrnd.h"
#include "matgen.h"
#include "ablasf.h"
#include "ablas.h"
#include "trfac.h"
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "matinv.h"
#include "hblas.h"
#include "sblas.h"
#include "ortfac.h"
#include "rotations.h"
#include "bdsvd.h"
#include "svd.h"
#include "xblas.h"
#include "densesolver.h"
#include "linmin.h"
#include "minlbfgs.h"
#include "minlm.h"
#include "lsfit.h"
#include "apserv.h"


/*************************************************************************
1-dimensional spline inteprolant
*************************************************************************/
struct spline1dinterpolant
{
    bool periodic;
    int n;
    int k;
    ap::real_1d_array x;
    ap::real_1d_array c;
};


/*************************************************************************
Spline fitting report:
    TaskRCond       reciprocal of task's condition number
    RMSError        RMS error
    AvgError        average error
    AvgRelError     average relative error (for non-zero Y[I])
    MaxError        maximum error
*************************************************************************/
struct spline1dfitreport
{
    double taskrcond;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double maxerror;
};




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
     spline1dinterpolant& c);


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
     spline1dinterpolant& c);


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
     spline1dinterpolant& c);


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
     spline1dinterpolant& c);


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
     spline1dinterpolant& c);


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
     spline1dfitreport& rep);


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
     spline1dfitreport& rep);


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
     spline1dfitreport& rep);


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
     spline1dfitreport& rep);


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
double spline1dcalc(const spline1dinterpolant& c, double x);


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
     double& d2s);


/*************************************************************************
This subroutine makes the copy of the spline.

INPUT PARAMETERS:
    C   -   spline interpolant.

Result:
    CC  -   spline copy

  -- ALGLIB PROJECT --
     Copyright 29.06.2007 by Bochkanov Sergey
*************************************************************************/
void spline1dcopy(const spline1dinterpolant& c, spline1dinterpolant& cc);


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
     ap::real_2d_array& tbl);


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
void spline1dlintransx(spline1dinterpolant& c, double a, double b);


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
void spline1dlintransy(spline1dinterpolant& c, double a, double b);


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
double spline1dintegrate(const spline1dinterpolant& c, double x);


#endif

