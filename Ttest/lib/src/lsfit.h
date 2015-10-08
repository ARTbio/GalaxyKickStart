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

#ifndef _lsfit_h
#define _lsfit_h

#include "ap.h"
#include "ialglib.h"

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


/*************************************************************************
Least squares fitting report:
    TaskRCond       reciprocal of task's condition number
    RMSError        RMS error
    AvgError        average error
    AvgRelError     average relative error (for non-zero Y[I])
    MaxError        maximum error
*************************************************************************/
struct lsfitreport
{
    double taskrcond;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double maxerror;
};


struct lsfitstate
{
    int n;
    int m;
    int k;
    double epsf;
    double epsx;
    int maxits;
    double stpmax;
    ap::real_2d_array taskx;
    ap::real_1d_array tasky;
    ap::real_1d_array w;
    bool cheapfg;
    bool havehess;
    bool needf;
    bool needfg;
    bool needfgh;
    int pointindex;
    ap::real_1d_array x;
    ap::real_1d_array c;
    double f;
    ap::real_1d_array g;
    ap::real_2d_array h;
    int repterminationtype;
    double reprmserror;
    double repavgerror;
    double repavgrelerror;
    double repmaxerror;
    minlmstate optstate;
    minlmreport optrep;
    ap::rcommstate rstate;
};




/*************************************************************************
Weighted linear least squares fitting.

QR decomposition is used to reduce task to MxM, then triangular solver  or
SVD-based solver is used depending on condition number of the  system.  It
allows to maximize speed and retain decent accuracy.

INPUT PARAMETERS:
    Y       -   array[0..N-1] Function values in  N  points.
    W       -   array[0..N-1]  Weights  corresponding to function  values.
                Each summand in square  sum  of  approximation  deviations
                from  given  values  is  multiplied  by  the   square   of
                corresponding weight.
    FMatrix -   a table of basis functions values, array[0..N-1, 0..M-1].
                FMatrix[I, J] - value of J-th basis function in I-th point.
    N       -   number of points used. N>=1.
    M       -   number of basis functions, M>=1.

OUTPUT PARAMETERS:
    Info    -   error code:
                * -4    internal SVD decomposition subroutine failed (very
                        rare and for degenerate systems only)
                * -1    incorrect N/M were specified
                *  1    task is solved
    C       -   decomposition coefficients, array[0..M-1]
    Rep     -   fitting report. Following fields are set:
                * Rep.TaskRCond     reciprocal of condition number
                * RMSError          rms error on the (X,Y).
                * AvgError          average error on the (X,Y).
                * AvgRelError       average relative error on the non-zero Y
                * MaxError          maximum error
                                    NON-WEIGHTED ERRORS ARE CALCULATED

SEE ALSO
    LSFitLinear
    LSFitLinearC
    LSFitLinearWC

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitlinearw(const ap::real_1d_array& y,
     const ap::real_1d_array& w,
     const ap::real_2d_array& fmatrix,
     int n,
     int m,
     int& info,
     ap::real_1d_array& c,
     lsfitreport& rep);


/*************************************************************************
Weighted constained linear least squares fitting.

This  is  variation  of LSFitLinearW(), which searchs for min|A*x=b| given
that  K  additional  constaints  C*x=bc are satisfied. It reduces original
task to modified one: min|B*y-d| WITHOUT constraints,  then LSFitLinearW()
is called.

INPUT PARAMETERS:
    Y       -   array[0..N-1] Function values in  N  points.
    W       -   array[0..N-1]  Weights  corresponding to function  values.
                Each summand in square  sum  of  approximation  deviations
                from  given  values  is  multiplied  by  the   square   of
                corresponding weight.
    FMatrix -   a table of basis functions values, array[0..N-1, 0..M-1].
                FMatrix[I,J] - value of J-th basis function in I-th point.
    CMatrix -   a table of constaints, array[0..K-1,0..M].
                I-th row of CMatrix corresponds to I-th linear constraint:
                CMatrix[I,0]*C[0] + ... + CMatrix[I,M-1]*C[M-1] = CMatrix[I,M]
    N       -   number of points used. N>=1.
    M       -   number of basis functions, M>=1.
    K       -   number of constraints, 0 <= K < M
                K=0 corresponds to absence of constraints.

OUTPUT PARAMETERS:
    Info    -   error code:
                * -4    internal SVD decomposition subroutine failed (very
                        rare and for degenerate systems only)
                * -3    either   too   many  constraints  (M   or   more),
                        degenerate  constraints   (some   constraints  are
                        repetead twice) or inconsistent  constraints  were
                        specified.
                * -1    incorrect N/M/K were specified
                *  1    task is solved
    C       -   decomposition coefficients, array[0..M-1]
    Rep     -   fitting report. Following fields are set:
                * RMSError          rms error on the (X,Y).
                * AvgError          average error on the (X,Y).
                * AvgRelError       average relative error on the non-zero Y
                * MaxError          maximum error
                                    NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroitine doesn't calculate task's condition number for K<>0.

SEE ALSO
    LSFitLinear
    LSFitLinearC
    LSFitLinearWC

  -- ALGLIB --
     Copyright 07.09.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitlinearwc(ap::real_1d_array y,
     const ap::real_1d_array& w,
     const ap::real_2d_array& fmatrix,
     ap::real_2d_array cmatrix,
     int n,
     int m,
     int k,
     int& info,
     ap::real_1d_array& c,
     lsfitreport& rep);


/*************************************************************************
Linear least squares fitting, without weights.

See LSFitLinearW for more information.

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitlinear(const ap::real_1d_array& y,
     const ap::real_2d_array& fmatrix,
     int n,
     int m,
     int& info,
     ap::real_1d_array& c,
     lsfitreport& rep);


/*************************************************************************
Constained linear least squares fitting, without weights.

See LSFitLinearWC() for more information.

  -- ALGLIB --
     Copyright 07.09.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitlinearc(ap::real_1d_array y,
     const ap::real_2d_array& fmatrix,
     const ap::real_2d_array& cmatrix,
     int n,
     int m,
     int k,
     int& info,
     ap::real_1d_array& c,
     lsfitreport& rep);


/*************************************************************************
Weighted nonlinear least squares fitting using gradient and Hessian.

Nonlinear task min(F(c)) is solved, where

    F(c) = (w[0]*(f(x[0],c)-y[0]))^2 + ... + (w[n-1]*(f(x[n-1],c)-y[n-1]))^2,
    
    * N is a number of points,
    * M is a dimension of a space points belong to,
    * K is a dimension of a space of parameters being fitted,
    * w is an N-dimensional vector of weight coefficients,
    * x is a set of N points, each of them is an M-dimensional vector,
    * c is a K-dimensional vector of parameters being fitted
    
This subroutine uses only f(x[i],c) and its gradient.
    
INPUT PARAMETERS:
    X       -   array[0..N-1,0..M-1], points (one row = one point)
    Y       -   array[0..N-1], function values.
    W       -   weights, array[0..N-1]
    C       -   array[0..K-1], initial approximation to the solution,
    N       -   number of points, N>1
    M       -   dimension of space
    K       -   number of parameters being fitted
    CheapFG -   boolean flag, which is:
                * True  if both function and gradient calculation complexity
                        are less than O(M^2).  An improved  algorithm  can
                        be  used  which corresponds  to  FGJ  scheme  from
                        MINLM unit.
                * False otherwise.
                        Standard Jacibian-bases  Levenberg-Marquardt  algo
                        will be used (FJ scheme).

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state between subsequent
                calls  of   LSFitNonlinearIteration.   Used  for  reverse
                communication.  This  structure   should   be  passed  to
                LSFitNonlinearIteration subroutine.

See also:
    LSFitNonlinearIteration
    LSFitNonlinearResults
    LSFitNonlinearFG (fitting without weights)
    LSFitNonlinearWFGH (fitting using Hessian)
    LSFitNonlinearFGH (fitting using Hessian, without weights)


  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitnonlinearwfg(const ap::real_2d_array& x,
     const ap::real_1d_array& y,
     const ap::real_1d_array& w,
     const ap::real_1d_array& c,
     int n,
     int m,
     int k,
     bool cheapfg,
     lsfitstate& state);


/*************************************************************************
Nonlinear least squares fitting, no individual weights.
See LSFitNonlinearWFG for more information.

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitnonlinearfg(const ap::real_2d_array& x,
     const ap::real_1d_array& y,
     const ap::real_1d_array& c,
     int n,
     int m,
     int k,
     bool cheapfg,
     lsfitstate& state);


/*************************************************************************
Weighted nonlinear least squares fitting using gradient/Hessian.

Nonlinear task min(F(c)) is solved, where

    F(c) = (w[0]*(f(x[0],c)-y[0]))^2 + ... + (w[n-1]*(f(x[n-1],c)-y[n-1]))^2,

    * N is a number of points,
    * M is a dimension of a space points belong to,
    * K is a dimension of a space of parameters being fitted,
    * w is an N-dimensional vector of weight coefficients,
    * x is a set of N points, each of them is an M-dimensional vector,
    * c is a K-dimensional vector of parameters being fitted

This subroutine uses f(x[i],c), its gradient and its Hessian.

See LSFitNonlinearWFG() subroutine for information about function
parameters.

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitnonlinearwfgh(const ap::real_2d_array& x,
     const ap::real_1d_array& y,
     const ap::real_1d_array& w,
     const ap::real_1d_array& c,
     int n,
     int m,
     int k,
     lsfitstate& state);


/*************************************************************************
Nonlinear least squares fitting using gradient/Hessian without  individual
weights. See LSFitNonlinearWFGH() for more information.


  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitnonlinearfgh(const ap::real_2d_array& x,
     const ap::real_1d_array& y,
     const ap::real_1d_array& c,
     int n,
     int m,
     int k,
     lsfitstate& state);


/*************************************************************************
Stopping conditions for nonlinear least squares fitting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be initialized
                with LSFitNonLinearCreate???()
    EpsF    -   stopping criterion. Algorithm stops if
                |F(k+1)-F(k)| <= EpsF*max{|F(k)|, |F(k+1)|, 1}
    EpsX    -   stopping criterion. Algorithm stops if
                |X(k+1)-X(k)| <= EpsX*(1+|X(k)|)
    MaxIts  -   stopping criterion. Algorithm stops after MaxIts iterations.
                MaxIts=0 means no stopping criterion.

NOTE

Passing EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to automatic
stopping criterion selection (according to the scheme used by MINLM unit).


  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitnonlinearsetcond(lsfitstate& state,
     double epsf,
     double epsx,
     int maxits);


/*************************************************************************
This function sets maximum step length

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be
                initialized with LSFitNonLinearCreate???()
    StpMax  -   maximum step length, >=0. Set StpMax to 0.0,  if you don't
                want to limit step length.

Use this subroutine when you optimize target function which contains exp()
or  other  fast  growing  functions,  and optimization algorithm makes too
large  steps  which  leads  to overflow. This function allows us to reject
steps  that  are  too  large  (and  therefore  expose  us  to the possible
overflow) without actually calculating function value at the x+stp*d.

NOTE: non-zero StpMax leads to moderate  performance  degradation  because
intermediate  step  of  preconditioned L-BFGS optimization is incompatible
with limits on step size.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void lsfitnonlinearsetstpmax(lsfitstate& state, double stpmax);


/*************************************************************************
Nonlinear least squares fitting. Algorithm iteration.

Called after inialization of the State structure with  LSFitNonlinearXXX()
subroutine. See HTML docs for examples.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between  subsequent
                calls and which is used for reverse communication. Must be
                initialized with LSFitNonlinearXXX() call first.

RESULT
1. If subroutine returned False, iterative algorithm has converged.
2. If subroutine returned True, then if:
* if State.NeedF=True,      function value F(X,C) is required
* if State.NeedFG=True,     function value F(X,C) and gradient  dF/dC(X,C)
                            are required
* if State.NeedFGH=True     function value F(X,C), gradient dF/dC(X,C) and
                            Hessian are required

One and only one of this fields can be set at time.

Function, its gradient and Hessian are calculated at  (X,C),  where  X  is
stored in State.X[0..M-1] and C is stored in State.C[0..K-1].

Results are stored:
* function value            -   in State.F
* gradient                  -   in State.G[0..K-1]
* Hessian                   -   in State.H[0..K-1,0..K-1]

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
bool lsfitnonlineariteration(lsfitstate& state);


/*************************************************************************
Nonlinear least squares fitting results.

Called after LSFitNonlinearIteration() returned False.

INPUT PARAMETERS:
    State   -   algorithm state (used by LSFitNonlinearIteration).

OUTPUT PARAMETERS:
    Info    -   completetion code:
                    * -1    incorrect parameters were specified
                    *  1    relative function improvement is no more than
                            EpsF.
                    *  2    relative step is no more than EpsX.
                    *  4    gradient norm is no more than EpsG
                    *  5    MaxIts steps was taken
    C       -   array[0..K-1], solution
    Rep     -   optimization report. Following fields are set:
                * Rep.TerminationType completetion code:
                * RMSError          rms error on the (X,Y).
                * AvgError          average error on the (X,Y).
                * AvgRelError       average relative error on the non-zero Y
                * MaxError          maximum error
                                    NON-WEIGHTED ERRORS ARE CALCULATED


  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************/
void lsfitnonlinearresults(const lsfitstate& state,
     int& info,
     ap::real_1d_array& c,
     lsfitreport& rep);


void lsfitscalexy(ap::real_1d_array& x,
     ap::real_1d_array& y,
     int n,
     ap::real_1d_array& xc,
     ap::real_1d_array& yc,
     const ap::integer_1d_array& dc,
     int k,
     double& xa,
     double& xb,
     double& sa,
     double& sb,
     ap::real_1d_array& xoriginal,
     ap::real_1d_array& yoriginal);


#endif

