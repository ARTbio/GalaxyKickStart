/*************************************************************************
Copyright (c) 2009, Sergey Bochkanov (ALGLIB project).

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

#ifndef _minlm_h
#define _minlm_h

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


struct minlmstate
{
    bool wrongparams;
    int n;
    int m;
    double epsg;
    double epsf;
    double epsx;
    int maxits;
    bool xrep;
    double stpmax;
    int flags;
    int usermode;
    ap::real_1d_array x;
    double f;
    ap::real_1d_array fi;
    ap::real_2d_array j;
    ap::real_2d_array h;
    ap::real_1d_array g;
    bool needf;
    bool needfg;
    bool needfgh;
    bool needfij;
    bool xupdated;
    minlbfgsstate internalstate;
    minlbfgsreport internalrep;
    ap::real_1d_array xprec;
    ap::real_1d_array xbase;
    ap::real_1d_array xdir;
    ap::real_1d_array gbase;
    ap::real_1d_array xprev;
    double fprev;
    ap::real_2d_array rawmodel;
    ap::real_2d_array model;
    ap::real_1d_array work;
    ap::rcommstate rstate;
    int repiterationscount;
    int repterminationtype;
    int repnfunc;
    int repnjac;
    int repngrad;
    int repnhess;
    int repncholesky;
    int solverinfo;
    densesolverreport solverrep;
    int invinfo;
    matinvreport invrep;
};


struct minlmreport
{
    int iterationscount;
    int terminationtype;
    int nfunc;
    int njac;
    int ngrad;
    int nhess;
    int ncholesky;
};




/*************************************************************************
    LEVENBERG-MARQUARDT-LIKE METHOD FOR NON-LINEAR OPTIMIZATION

Optimization using function gradient and Hessian.  Algorithm -  Levenberg-
Marquardt   modification   with   L-BFGS   pre-optimization  and  internal
pre-conditioned L-BFGS optimization after each Levenberg-Marquardt step.

Function F has general form (not "sum-of-squares"):

    F = F(x[0], ..., x[n-1])

EXAMPLE

See HTML-documentation.

INPUT PARAMETERS:
    N       -   dimension, N>1
    X       -   initial solution, array[0..N-1]

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state between subsequent
                calls of MinLMIteration. Used for reverse communication.
                This structure should be passed to MinLMIteration subroutine.

See also MinLMIteration, MinLMResults.

NOTES:

1. you may tune stopping conditions with MinLMSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLMSetStpMax() function to bound algorithm's steps.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatefgh(const int& n,
     const ap::real_1d_array& x,
     minlmstate& state);


/*************************************************************************
    LEVENBERG-MARQUARDT-LIKE METHOD FOR NON-LINEAR OPTIMIZATION

Optimization using function gradient and Jacobian.  Algorithm -  Levenberg-
Marquardt   modification   with   L-BFGS   pre-optimization  and  internal
pre-conditioned L-BFGS optimization after each Levenberg-Marquardt step.

Function F is represented as sum of squares:

    F = f[0]^2(x[0],...,x[n-1]) + ... + f[m-1]^2(x[0],...,x[n-1])

EXAMPLE

See HTML-documentation.

INPUT PARAMETERS:
    N       -   dimension, N>1
    M       -   number of functions f[i]
    X       -   initial solution, array[0..N-1]

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state between subsequent
                calls of MinLMIteration. Used for reverse communication.
                This structure should be passed to MinLMIteration subroutine.

See also MinLMIteration, MinLMResults.

NOTES:

1. you may tune stopping conditions with MinLMSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLMSetStpMax() function to bound algorithm's steps.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatefgj(const int& n,
     const int& m,
     const ap::real_1d_array& x,
     minlmstate& state);


/*************************************************************************
    CLASSIC LEVENBERG-MARQUARDT METHOD FOR NON-LINEAR OPTIMIZATION

Optimization using Jacobi matrix. Algorithm  -  classic Levenberg-Marquardt
method.

Function F is represented as sum of squares:

    F = f[0]^2(x[0],...,x[n-1]) + ... + f[m-1]^2(x[0],...,x[n-1])

EXAMPLE

See HTML-documentation.

INPUT PARAMETERS:
    N       -   dimension, N>1
    M       -   number of functions f[i]
    X       -   initial solution, array[0..N-1]

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state between subsequent
                calls of MinLMIteration. Used for reverse communication.
                This structure should be passed to MinLMIteration subroutine.

See also MinLMIteration, MinLMResults.

NOTES:

1. you may tune stopping conditions with MinLMSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLMSetStpMax() function to bound algorithm's steps.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatefj(const int& n,
     const int& m,
     const ap::real_1d_array& x,
     minlmstate& state);


/*************************************************************************
This function sets stopping conditions for Levenberg-Marquardt optimization
algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be initialized
                with MinLMCreate???()
    EpsG    -   >=0
                The  subroutine  finishes  its  work   if   the  condition
                ||G||<EpsG is satisfied, where ||.|| means Euclidian norm,
                G - gradient.
    EpsF    -   >=0
                The  subroutine  finishes  its work if on k+1-th iteration
                the  condition  |F(k+1)-F(k)|<=EpsF*max{|F(k)|,|F(k+1)|,1}
                is satisfied.
    EpsX    -   >=0
                The subroutine finishes its work if  on  k+1-th  iteration
                the condition |X(k+1)-X(k)| <= EpsX is fulfilled.
    MaxIts  -   maximum number of iterations. If MaxIts=0, the  number  of
                iterations   is    unlimited.   Only   Levenberg-Marquardt
                iterations  are  counted  (L-BFGS/CG  iterations  are  NOT
                counted  because their cost is very low copared to that of
                LM).

Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to
automatic stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlmsetcond(minlmstate& state,
     double epsg,
     double epsf,
     double epsx,
     int maxits);


/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be
                initialized with MinLMCreate???()
    NeedXRep-   whether iteration reports are needed or not

Usually  algorithm  returns  from  MinLMIteration()  only  when  it  needs
function/gradient/Hessian. However, with this function we can let it  stop
after  each  iteration  (one iteration may include  more than one function
evaluation), which is indicated by XUpdated field.

Both Levenberg-Marquardt and L-BFGS iterations are reported.


  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlmsetxrep(minlmstate& state, bool needxrep);


/*************************************************************************
This function sets maximum step length

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be
                initialized with MinCGCreate???()
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
void minlmsetstpmax(minlmstate& state, double stpmax);


/*************************************************************************
One Levenberg-Marquardt iteration.

Called after inialization of State structure with MinLMXXX subroutine.
See HTML docs for examples.

Input parameters:
    State   -   structure which stores algorithm state between subsequent
                calls and which is used for reverse communication. Must be
                initialized with MinLMXXX call first.

If subroutine returned False, iterative algorithm has converged.

If subroutine returned True, then:
* if State.NeedF=True,      -   function value F at State.X[0..N-1]
                                is required
* if State.NeedFG=True      -   function value F and gradient G
                                are required
* if State.NeedFiJ=True     -   function vector f[i] and Jacobi matrix J
                                are required
* if State.NeedFGH=True     -   function value F, gradient G and Hesian H
                                are required
* if State.XUpdated=True    -   algorithm reports about new iteration,
                                State.X contains current point,
                                State.F contains function value.

One and only one of this fields can be set at time.

Results are stored:
* function value            -   in MinLMState.F
* gradient                  -   in MinLMState.G[0..N-1]
* Jacobi matrix             -   in MinLMState.J[0..M-1,0..N-1]
* Hessian                   -   in MinLMState.H[0..N-1,0..N-1]

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
bool minlmiteration(minlmstate& state);


/*************************************************************************
Levenberg-Marquardt algorithm results

Called after MinLMIteration returned False.

Input parameters:
    State   -   algorithm state (used by MinLMIteration).

Output parameters:
    X       -   array[0..N-1], solution
    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    * -1    incorrect parameters were specified
                    *  1    relative function improvement is no more than
                            EpsF.
                    *  2    relative step is no more than EpsX.
                    *  4    gradient is no more than EpsG.
                    *  5    MaxIts steps was taken
                    *  7    stopping conditions are too stringent,
                            further improvement is impossible
                * Rep.IterationsCount contains iterations count
                * Rep.NFunc     - number of function calculations
                * Rep.NJac      - number of Jacobi matrix calculations
                * Rep.NGrad     - number of gradient calculations
                * Rep.NHess     - number of Hessian calculations
                * Rep.NCholesky - number of Cholesky decomposition calculations

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmresults(const minlmstate& state,
     ap::real_1d_array& x,
     minlmreport& rep);


#endif

