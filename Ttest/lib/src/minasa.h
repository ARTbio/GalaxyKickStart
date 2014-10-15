/*************************************************************************
Copyright (c) 2010, Sergey Bochkanov (ALGLIB project).

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

#ifndef _minasa_h
#define _minasa_h

#include "ap.h"
#include "ialglib.h"

#include "linmin.h"


struct minasastate
{
    int n;
    double epsg;
    double epsf;
    double epsx;
    int maxits;
    bool xrep;
    double stpmax;
    int cgtype;
    int k;
    int nfev;
    int mcstage;
    ap::real_1d_array bndl;
    ap::real_1d_array bndu;
    int curalgo;
    int acount;
    double mu;
    double finit;
    double dginit;
    ap::real_1d_array ak;
    ap::real_1d_array xk;
    ap::real_1d_array dk;
    ap::real_1d_array an;
    ap::real_1d_array xn;
    ap::real_1d_array dn;
    ap::real_1d_array d;
    double fold;
    double stp;
    ap::real_1d_array work;
    ap::real_1d_array yk;
    ap::real_1d_array gc;
    ap::real_1d_array x;
    double f;
    ap::real_1d_array g;
    bool needfg;
    bool xupdated;
    ap::rcommstate rstate;
    int repiterationscount;
    int repnfev;
    int repterminationtype;
    int debugrestartscount;
    linminstate lstate;
    double betahs;
    double betady;
};


struct minasareport
{
    int iterationscount;
    int nfev;
    int terminationtype;
    int activeconstraints;
};




/*************************************************************************
              NONLINEAR BOUND CONSTRAINED OPTIMIZATION USING
                               MODIFIED
                   WILLIAM W. HAGER AND HONGCHAO ZHANG
                         ACTIVE SET ALGORITHM

The  subroutine  minimizes  function  F(x)  of  N  arguments  with   bound
constraints: BndL[i] <= x[i] <= BndU[i]

This method is  globally  convergent  as  long  as  grad(f)  is  Lipschitz
continuous on a level set: L = { x : f(x)<=f(x0) }.

INPUT PARAMETERS:
    N       -   problem dimension. N>0
    X       -   initial solution approximation, array[0..N-1].
    BndL    -   lower bounds, array[0..N-1].
                all elements MUST be specified,  i.e.  all  variables  are
                bounded. However, if some (all) variables  are  unbounded,
                you may specify very small number as bound: -1000,  -1.0E6
                or -1.0E300, or something like that.
    BndU    -   upper bounds, array[0..N-1].
                all elements MUST be specified,  i.e.  all  variables  are
                bounded. However, if some (all) variables  are  unbounded,
                you may specify very large number as bound: +1000,  +1.0E6
                or +1.0E300, or something like that.
    EpsG    -   positive number which  defines  a  precision  of  search.  The
                subroutine finishes its work if the condition ||G|| < EpsG  is
                satisfied, where ||.|| means Euclidian norm, G - gradient, X -
                current approximation.
    EpsF    -   positive number which  defines  a  precision  of  search.  The
                subroutine finishes its work if on iteration  number  k+1  the
                condition |F(k+1)-F(k)| <= EpsF*max{|F(k)|, |F(k+1)|, 1}    is
                satisfied.
    EpsX    -   positive number which  defines  a  precision  of  search.  The
                subroutine finishes its work if on iteration number k+1    the
                condition |X(k+1)-X(k)| <= EpsX is fulfilled.
    MaxIts  -   maximum number of iterations. If MaxIts=0, the number of
                iterations is unlimited.

OUTPUT PARAMETERS:
    State - structure used for reverse communication.

This function  initializes  State   structure  with  default  optimization
parameters (stopping conditions, step size, etc.).  Use  MinASASet??????()
functions to tune optimization parameters.

After   all   optimization   parameters   are   tuned,   you   should  use
MinASAIteration() function to advance algorithm iterations.

NOTES:

1. you may tune stopping conditions with MinASASetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinASASetStpMax() function to bound algorithm's steps.

  -- ALGLIB --
     Copyright 25.03.2010 by Bochkanov Sergey
*************************************************************************/
void minasacreate(int n,
     const ap::real_1d_array& x,
     const ap::real_1d_array& bndl,
     const ap::real_1d_array& bndu,
     minasastate& state);


/*************************************************************************
This function sets stopping conditions for the ASA optimization algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be initialized
                with MinASACreate()
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
                iterations is unlimited.

Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to
automatic stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minasasetcond(minasastate& state,
     double epsg,
     double epsf,
     double epsx,
     int maxits);


/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be
                initialized with MinASACreate()
    NeedXRep-   whether iteration reports are needed or not

Usually  algorithm  returns from  MinASAIteration()  only  when  it  needs
function/gradient. However, with this function we can let  it  stop  after
each  iteration  (one  iteration  may  include   more  than  one  function
evaluation), which is indicated by XUpdated field.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minasasetxrep(minasastate& state, bool needxrep);


/*************************************************************************
This function sets optimization algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be
                initialized with MinASACreate()
    UAType  -   algorithm type:
                * -1    automatic selection of the best algorithm
                * 0     DY (Dai and Yuan) algorithm
                * 1     Hybrid DY-HS algorithm

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minasasetalgorithm(minasastate& state, int algotype);


/*************************************************************************
This function sets maximum step length

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be
                initialized with MinCGCreate()
    StpMax  -   maximum step length, >=0. Set StpMax to 0.0,  if you don't
                want to limit step length.

Use this subroutine when you optimize target function which contains exp()
or  other  fast  growing  functions,  and optimization algorithm makes too
large  steps  which  leads  to overflow. This function allows us to reject
steps  that  are  too  large  (and  therefore  expose  us  to the possible
overflow) without actually calculating function value at the x+stp*d.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minasasetstpmax(minasastate& state, double stpmax);


/*************************************************************************
One ASA iteration

Called after initialization with MinASACreate.
See HTML documentation for examples.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be initialized
                with MinASACreate.
RESULT:
* if function returned False, iterative proces has converged.
  Use MinLBFGSResults() to obtain optimization results.
* if subroutine returned True, then, depending on structure fields, we
  have one of the following situations


=== FUNC/GRAD REQUEST ===
State.NeedFG is True => function value/gradient are needed.
Caller should calculate function value State.F and gradient
State.G[0..N-1] at State.X[0..N-1] and call MinLBFGSIteration() again.

=== NEW INTERATION IS REPORTED ===
State.XUpdated is True => one more iteration was made.
State.X contains current position, State.F contains function value at X.
You can read info from these fields, but never modify  them  because  they
contain the only copy of optimization algorithm state.

One and only one of these fields (NeedFG, XUpdated) is true on return. New
iterations are reported only when reports  are  explicitly  turned  on  by
MinLBFGSSetXRep() function, so if you never called it, you can expect that
NeedFG is always True.


  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
bool minasaiteration(minasastate& state);


/*************************************************************************
Conjugate gradient results

Called after MinASA returned False.

INPUT PARAMETERS:
    State   -   algorithm state (used by MinASAIteration).

OUTPUT PARAMETERS:
    X       -   array[0..N-1], solution
    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    * -2    rounding errors prevent further improvement.
                            X contains best point found.
                    * -1    incorrect parameters were specified
                    *  1    relative function improvement is no more than
                            EpsF.
                    *  2    relative step is no more than EpsX.
                    *  4    gradient norm is no more than EpsG
                    *  5    MaxIts steps was taken
                    *  7    stopping conditions are too stringent,
                            further improvement is impossible
                * Rep.IterationsCount contains iterations count
                * NFEV countains number of function calculations
                * ActiveConstraints contains number of active constraints

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void minasaresults(const minasastate& state,
     ap::real_1d_array& x,
     minasareport& rep);


#endif

