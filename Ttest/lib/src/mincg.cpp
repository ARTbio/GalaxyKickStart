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

#include <stdafx.h>
#include "mincg.h"

static void clearrequestfields(mincgstate& state);

/*************************************************************************
        NONLINEAR CONJUGATE GRADIENT METHOD

The subroutine minimizes function F(x) of N arguments by using one of  the
nonlinear conjugate gradient methods.

These CG methods are globally convergent (even on non-convex functions) as
long as grad(f) is Lipschitz continuous in  a  some  neighborhood  of  the
L = { x : f(x)<=f(x0) }.

INPUT PARAMETERS:
    N       -   problem dimension. N>0
    X       -   initial solution approximation, array[0..N-1].
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

See also MinCGIteration, MinCGResults

NOTE:

Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to
automatic stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 25.03.2010 by Bochkanov Sergey
*************************************************************************/
void mincgcreate(int n, const ap::real_1d_array& x, mincgstate& state)
{

    ap::ap_error::make_assertion(n>=1, "MinCGCreate: N too small!");
    
    //
    // Initialize
    //
    state.n = n;
    mincgsetcond(state, double(0), double(0), double(0), 0);
    mincgsetxrep(state, false);
    mincgsetstpmax(state, double(0));
    mincgsetcgtype(state, -1);
    state.xk.setlength(n);
    state.dk.setlength(n);
    state.xn.setlength(n);
    state.dn.setlength(n);
    state.x.setlength(n);
    state.d.setlength(n);
    state.g.setlength(n);
    state.work.setlength(n);
    state.yk.setlength(n);
    
    //
    // Prepare first run
    //
    ap::vmove(&state.x(0), 1, &x(0), 1, ap::vlen(0,n-1));
    state.rstate.ia.setbounds(0, 2);
    state.rstate.ra.setbounds(0, 2);
    state.rstate.stage = -1;
}


/*************************************************************************
This function sets stopping conditions for CG optimization algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be initialized
                with MinCGCreate()
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
void mincgsetcond(mincgstate& state,
     double epsg,
     double epsf,
     double epsx,
     int maxits)
{

    ap::ap_error::make_assertion(ap::fp_greater_eq(epsg,0), "MinCGSetCond: negative EpsG!");
    ap::ap_error::make_assertion(ap::fp_greater_eq(epsf,0), "MinCGSetCond: negative EpsF!");
    ap::ap_error::make_assertion(ap::fp_greater_eq(epsx,0), "MinCGSetCond: negative EpsX!");
    ap::ap_error::make_assertion(maxits>=0, "MinCGSetCond: negative MaxIts!");
    if( ap::fp_eq(epsg,0)&&ap::fp_eq(epsf,0)&&ap::fp_eq(epsx,0)&&maxits==0 )
    {
        epsx = 1.0E-6;
    }
    state.epsg = epsg;
    state.epsf = epsf;
    state.epsx = epsx;
    state.maxits = maxits;
}


/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be
                initialized with MinCGCreate()
    NeedXRep-   whether iteration reports are needed or not

Usually  algorithm  returns  from  MinCGIteration()  only  when  it  needs
function/gradient. However, with this function we can let  it  stop  after
each  iteration  (one  iteration  may  include   more  than  one  function
evaluation), which is indicated by XUpdated field.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsetxrep(mincgstate& state, bool needxrep)
{

    state.xrep = needxrep;
}


/*************************************************************************
This function sets CG algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be
                initialized with MinCGCreate()
    CGType  -   algorithm type:
                * -1    automatic selection of the best algorithm
                * 0     DY (Dai and Yuan) algorithm
                * 1     Hybrid DY-HS algorithm

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsetcgtype(mincgstate& state, int cgtype)
{

    ap::ap_error::make_assertion(cgtype>=-1&&cgtype<=1, "MinCGSetCGType: incorrect CGType!");
    if( cgtype==-1 )
    {
        cgtype = 1;
    }
    state.cgtype = cgtype;
}


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
void mincgsetstpmax(mincgstate& state, double stpmax)
{

    ap::ap_error::make_assertion(ap::fp_greater_eq(stpmax,0), "MinCGSetStpMax: StpMax<0!");
    state.stpmax = stpmax;
}


/*************************************************************************
One conjugate gradient iteration

Called after initialization with MinCG.
See HTML documentation for examples.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be initialized
                with MinCG.

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
     Copyright 20.04.2009 by Bochkanov Sergey
*************************************************************************/
bool mincgiteration(mincgstate& state)
{
    bool result;
    int n;
    int i;
    double betak;
    double v;
    double vv;
    int mcinfo;

    
    //
    // Reverse communication preparations
    // I know it looks ugly, but it works the same way
    // anywhere from C++ to Python.
    //
    // This code initializes locals by:
    // * random values determined during code
    //   generation - on first subroutine call
    // * values from previous call - on subsequent calls
    //
    if( state.rstate.stage>=0 )
    {
        n = state.rstate.ia(0);
        i = state.rstate.ia(1);
        mcinfo = state.rstate.ia(2);
        betak = state.rstate.ra(0);
        v = state.rstate.ra(1);
        vv = state.rstate.ra(2);
    }
    else
    {
        n = -983;
        i = -989;
        mcinfo = -834;
        betak = 900;
        v = -287;
        vv = 364;
    }
    if( state.rstate.stage==0 )
    {
        goto lbl_0;
    }
    if( state.rstate.stage==1 )
    {
        goto lbl_1;
    }
    if( state.rstate.stage==2 )
    {
        goto lbl_2;
    }
    if( state.rstate.stage==3 )
    {
        goto lbl_3;
    }
    
    //
    // Routine body
    //
    
    //
    // Prepare
    //
    n = state.n;
    state.repterminationtype = 0;
    state.repiterationscount = 0;
    state.repnfev = 0;
    state.debugrestartscount = 0;
    
    //
    // Calculate F/G, initialize algorithm
    //
    clearrequestfields(state);
    state.needfg = true;
    state.rstate.stage = 0;
    goto lbl_rcomm;
lbl_0:
    if( !state.xrep )
    {
        goto lbl_4;
    }
    clearrequestfields(state);
    state.xupdated = true;
    state.rstate.stage = 1;
    goto lbl_rcomm;
lbl_1:
lbl_4:
    v = ap::vdotproduct(&state.g(0), 1, &state.g(0), 1, ap::vlen(0,n-1));
    v = sqrt(v);
    if( ap::fp_eq(v,0) )
    {
        state.repterminationtype = 4;
        result = false;
        return result;
    }
    state.repnfev = 1;
    state.k = 0;
    state.fold = state.f;
    ap::vmove(&state.xk(0), 1, &state.x(0), 1, ap::vlen(0,n-1));
    ap::vmoveneg(&state.dk(0), 1, &state.g(0), 1, ap::vlen(0,n-1));
    
    //
    // Main cycle
    //
lbl_6:
    if( false )
    {
        goto lbl_7;
    }
    
    //
    // Store G[k] for later calculation of Y[k]
    //
    ap::vmoveneg(&state.yk(0), 1, &state.g(0), 1, ap::vlen(0,n-1));
    
    //
    // Calculate X(k+1): minimize F(x+alpha*d)
    //
    ap::vmove(&state.d(0), 1, &state.dk(0), 1, ap::vlen(0,n-1));
    ap::vmove(&state.x(0), 1, &state.xk(0), 1, ap::vlen(0,n-1));
    state.mcstage = 0;
    state.stp = 1.0;
    linminnormalized(state.d, state.stp, n);
    mcsrch(n, state.x, state.f, state.g, state.d, state.stp, state.stpmax, mcinfo, state.nfev, state.work, state.lstate, state.mcstage);
lbl_8:
    if( state.mcstage==0 )
    {
        goto lbl_9;
    }
    clearrequestfields(state);
    state.needfg = true;
    state.rstate.stage = 2;
    goto lbl_rcomm;
lbl_2:
    mcsrch(n, state.x, state.f, state.g, state.d, state.stp, state.stpmax, mcinfo, state.nfev, state.work, state.lstate, state.mcstage);
    goto lbl_8;
lbl_9:
    if( !state.xrep )
    {
        goto lbl_10;
    }
    clearrequestfields(state);
    state.xupdated = true;
    state.rstate.stage = 3;
    goto lbl_rcomm;
lbl_3:
lbl_10:
    ap::vmove(&state.xn(0), 1, &state.x(0), 1, ap::vlen(0,n-1));
    if( mcinfo==1 )
    {
        
        //
        // Standard Wolfe conditions hold
        // Calculate Y[K] and BetaK
        //
        ap::vadd(&state.yk(0), 1, &state.g(0), 1, ap::vlen(0,n-1));
        vv = ap::vdotproduct(&state.yk(0), 1, &state.dk(0), 1, ap::vlen(0,n-1));
        v = ap::vdotproduct(&state.g(0), 1, &state.g(0), 1, ap::vlen(0,n-1));
        state.betady = v/vv;
        v = ap::vdotproduct(&state.g(0), 1, &state.yk(0), 1, ap::vlen(0,n-1));
        state.betahs = v/vv;
        if( state.cgtype==0 )
        {
            betak = state.betady;
        }
        if( state.cgtype==1 )
        {
            betak = ap::maxreal(double(0), ap::minreal(state.betady, state.betahs));
        }
    }
    else
    {
        
        //
        // Something is wrong (may be function is too wild or too flat).
        //
        // We'll set BetaK=0, which will restart CG algorithm.
        // We can stop later (during normal checks) if stopping conditions are met.
        //
        betak = 0;
        state.debugrestartscount = state.debugrestartscount+1;
    }
    
    //
    // Calculate D(k+1)
    //
    ap::vmoveneg(&state.dn(0), 1, &state.g(0), 1, ap::vlen(0,n-1));
    ap::vadd(&state.dn(0), 1, &state.dk(0), 1, ap::vlen(0,n-1), betak);
    
    //
    // Update information and Hessian.
    // Check stopping conditions.
    //
    state.repnfev = state.repnfev+state.nfev;
    state.repiterationscount = state.repiterationscount+1;
    if( state.repiterationscount>=state.maxits&&state.maxits>0 )
    {
        
        //
        // Too many iterations
        //
        state.repterminationtype = 5;
        result = false;
        return result;
    }
    v = ap::vdotproduct(&state.g(0), 1, &state.g(0), 1, ap::vlen(0,n-1));
    if( ap::fp_less_eq(sqrt(v),state.epsg) )
    {
        
        //
        // Gradient is small enough
        //
        state.repterminationtype = 4;
        result = false;
        return result;
    }
    if( ap::fp_less_eq(state.fold-state.f,state.epsf*ap::maxreal(fabs(state.fold), ap::maxreal(fabs(state.f), 1.0))) )
    {
        
        //
        // F(k+1)-F(k) is small enough
        //
        state.repterminationtype = 1;
        result = false;
        return result;
    }
    v = ap::vdotproduct(&state.d(0), 1, &state.d(0), 1, ap::vlen(0,n-1));
    if( ap::fp_less_eq(sqrt(v)*state.stp,state.epsx) )
    {
        
        //
        // X(k+1)-X(k) is small enough
        //
        state.repterminationtype = 2;
        result = false;
        return result;
    }
    
    //
    // Shift Xk/Dk, update other information
    //
    ap::vmove(&state.xk(0), 1, &state.xn(0), 1, ap::vlen(0,n-1));
    ap::vmove(&state.dk(0), 1, &state.dn(0), 1, ap::vlen(0,n-1));
    state.fold = state.f;
    state.k = state.k+1;
    goto lbl_6;
lbl_7:
    result = false;
    return result;
    
    //
    // Saving state
    //
lbl_rcomm:
    result = true;
    state.rstate.ia(0) = n;
    state.rstate.ia(1) = i;
    state.rstate.ia(2) = mcinfo;
    state.rstate.ra(0) = betak;
    state.rstate.ra(1) = v;
    state.rstate.ra(2) = vv;
    return result;
}


/*************************************************************************
Conjugate gradient results

Called after MinCG returned False.

INPUT PARAMETERS:
    State   -   algorithm state (used by MinCGIteration).

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
     Copyright 20.04.2009 by Bochkanov Sergey
*************************************************************************/
void mincgresults(const mincgstate& state,
     ap::real_1d_array& x,
     mincgreport& rep)
{

    x.setbounds(0, state.n-1);
    ap::vmove(&x(0), 1, &state.xn(0), 1, ap::vlen(0,state.n-1));
    rep.iterationscount = state.repiterationscount;
    rep.nfev = state.repnfev;
    rep.terminationtype = state.repterminationtype;
}


/*************************************************************************
Clears request fileds (to be sure that we don't forgot to clear something)
*************************************************************************/
static void clearrequestfields(mincgstate& state)
{

    state.needfg = false;
    state.xupdated = false;
}




