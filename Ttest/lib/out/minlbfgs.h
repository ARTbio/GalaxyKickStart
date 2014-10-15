/*************************************************************************
Copyright (c) 2007-2008, Sergey Bochkanov (ALGLIB project).

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

#ifndef _minlbfgs_h
#define _minlbfgs_h

#include "ap.h"
#include "ialglib.h"

#include "linmin.h"


struct minlbfgsstate
{
    int n;
    int m;
    double epsg;
    double epsf;
    double epsx;
    int maxits;
    int flags;
    bool xrep;
    double stpmax;
    int nfev;
    int mcstage;
    int k;
    int q;
    int p;
    ap::real_1d_array rho;
    ap::real_2d_array y;
    ap::real_2d_array s;
    ap::real_1d_array theta;
    ap::real_1d_array d;
    double stp;
    ap::real_1d_array work;
    double fold;
    double gammak;
    ap::real_1d_array x;
    double f;
    ap::real_1d_array g;
    bool needfg;
    bool xupdated;
    ap::rcommstate rstate;
    int repiterationscount;
    int repnfev;
    int repterminationtype;
    linminstate lstate;
};


struct minlbfgsreport
{
    int iterationscount;
    int nfev;
    int terminationtype;
};




/*************************************************************************
        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION

The subroutine minimizes function F(x) of N arguments by  using  a  quasi-
Newton method (LBFGS scheme) which is optimized to use  a  minimum  amount
of memory.

The subroutine generates the approximation of an inverse Hessian matrix by
using information about the last M steps of the algorithm  (instead of N).
It lessens a required amount of memory from a value  of  order  N^2  to  a
value of order 2*N*M.

INPUT PARAMETERS:
    N       -   problem dimension. N>0
    M       -   number of corrections in the BFGS scheme of Hessian
                approximation update. Recommended value:  3<=M<=7. The smaller
                value causes worse convergence, the bigger will  not  cause  a
                considerably better convergence, but will cause a fall in  the
                performance. M<=N.
    X       -   initial solution approximation, array[0..N-1].

OUTPUT PARAMETERS:
    State   -   structure used for reverse communication.
    
This function  initializes  State   structure  with  default  optimization
parameters (stopping conditions, step size, etc.). Use MinLBFGSSet??????()
functions to tune optimization parameters.

After   all   optimization   parameters   are   tuned,   you   should  use
MinLBFGSIteration() function to advance algorithm iterations.

NOTES:

1. you may tune stopping conditions with MinLBFGSSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLBFGSSetStpMax() function to bound algorithm's  steps.  However,
   L-BFGS rarely needs such a tuning.


  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgscreate(int n,
     int m,
     const ap::real_1d_array& x,
     minlbfgsstate& state);


/*************************************************************************
This function sets stopping conditions for L-BFGS optimization algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be initialized
                with MinLBFGSCreate()
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
void minlbfgssetcond(minlbfgsstate& state,
     double epsg,
     double epsf,
     double epsx,
     int maxits);


/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be
                initialized with MinLBFGSCreate()
    NeedXRep-   whether iteration reports are needed or not

Usually algorithm returns  from  MinLBFGSIteration()  only when  it  needs
function/gradient/ (which is indicated by NeedFG field. However, with this
function we can let it  stop  after  each  iteration  (one  iteration  may
include more than one function evaluation), which is indicated by XUpdated
field.


  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetxrep(minlbfgsstate& state, bool needxrep);


/*************************************************************************
This function sets maximum step length

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be
                initialized with MinLBFGSCreate()
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
void minlbfgssetstpmax(minlbfgsstate& state, double stpmax);


/*************************************************************************
Extended subroutine for internal use only.

Accepts additional parameters:

    Flags - additional settings:
            * Flags = 0     means no additional settings
            * Flags = 1     "do not allocate memory". used when solving
                            a many subsequent tasks with  same N/M  values.
                            First  call MUST  be without this flag bit set,
                            subsequent  calls   of   MinLBFGS   with   same
                            MinLBFGSState structure can set Flags to 1.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgscreatex(int n,
     int m,
     const ap::real_1d_array& x,
     int flags,
     minlbfgsstate& state);


/*************************************************************************
L-BFGS iterations

Called after initialization with MinLBFGSCreate() function.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be initialized
                with MinLBFGSCreate()

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
bool minlbfgsiteration(minlbfgsstate& state);


/*************************************************************************
L-BFGS algorithm results

Called after MinLBFGSIteration() returned False.

INPUT PARAMETERS:
    State   -   algorithm state (used by MinLBFGSIteration).

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

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgsresults(const minlbfgsstate& state,
     ap::real_1d_array& x,
     minlbfgsreport& rep);


#endif

