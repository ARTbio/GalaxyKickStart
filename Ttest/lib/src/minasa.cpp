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
#include "minasa.h"

static const int n1 = 2;
static const int n2 = 2;
static const double stpmin = 1.0E-300;
static const double gpaftol = 0.0001;
static const double gpadecay = 0.5;
static const double asarho = 0.5;

static double asaboundval(double x, double b1, double b2);
static double asaboundedantigradnorm(const minasastate& state);
static double asaginorm(const minasastate& state);
static double asad1norm(const minasastate& state);
static bool asauisempty(const minasastate& state);
static bool asawanttounstick(const minasastate& state);
static void clearrequestfields(minasastate& state);

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
     minasastate& state)
{
    int i;

    ap::ap_error::make_assertion(n>=1, "MinASA: N too small!");
    for(i = 0; i <= n-1; i++)
    {
        ap::ap_error::make_assertion(ap::fp_less_eq(bndl(i),bndu(i)), "MinASA: inconsistent bounds!");
        ap::ap_error::make_assertion(ap::fp_less_eq(bndl(i),x(i)), "MinASA: infeasible X!");
        ap::ap_error::make_assertion(ap::fp_less_eq(x(i),bndu(i)), "MinASA: infeasible X!");
    }
    
    //
    // Initialize
    //
    state.n = n;
    minasasetcond(state, double(0), double(0), double(0), 0);
    minasasetxrep(state, false);
    minasasetstpmax(state, double(0));
    minasasetalgorithm(state, -1);
    state.bndl.setlength(n);
    state.bndu.setlength(n);
    state.ak.setlength(n);
    state.xk.setlength(n);
    state.dk.setlength(n);
    state.an.setlength(n);
    state.xn.setlength(n);
    state.dn.setlength(n);
    state.x.setlength(n);
    state.d.setlength(n);
    state.g.setlength(n);
    state.gc.setlength(n);
    state.work.setlength(n);
    state.yk.setlength(n);
    ap::vmove(&state.bndl(0), 1, &bndl(0), 1, ap::vlen(0,n-1));
    ap::vmove(&state.bndu(0), 1, &bndu(0), 1, ap::vlen(0,n-1));
    
    //
    // Prepare first run
    //
    ap::vmove(&state.x(0), 1, &x(0), 1, ap::vlen(0,n-1));
    state.rstate.ia.setbounds(0, 3);
    state.rstate.ba.setbounds(0, 1);
    state.rstate.ra.setbounds(0, 2);
    state.rstate.stage = -1;
}


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
     int maxits)
{

    ap::ap_error::make_assertion(ap::fp_greater_eq(epsg,0), "MinASASetCond: negative EpsG!");
    ap::ap_error::make_assertion(ap::fp_greater_eq(epsf,0), "MinASASetCond: negative EpsF!");
    ap::ap_error::make_assertion(ap::fp_greater_eq(epsx,0), "MinASASetCond: negative EpsX!");
    ap::ap_error::make_assertion(maxits>=0, "MinASASetCond: negative MaxIts!");
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
                initialized with MinASACreate()
    NeedXRep-   whether iteration reports are needed or not

Usually  algorithm  returns from  MinASAIteration()  only  when  it  needs
function/gradient. However, with this function we can let  it  stop  after
each  iteration  (one  iteration  may  include   more  than  one  function
evaluation), which is indicated by XUpdated field.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minasasetxrep(minasastate& state, bool needxrep)
{

    state.xrep = needxrep;
}


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
void minasasetalgorithm(minasastate& state, int algotype)
{

    ap::ap_error::make_assertion(algotype>=-1&&algotype<=1, "MinASASetAlgorithm: incorrect AlgoType!");
    if( algotype==-1 )
    {
        algotype = 1;
    }
    state.cgtype = algotype;
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
void minasasetstpmax(minasastate& state, double stpmax)
{

    ap::ap_error::make_assertion(ap::fp_greater_eq(stpmax,0), "MinASASetStpMax: StpMax<0!");
    state.stpmax = stpmax;
}


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
bool minasaiteration(minasastate& state)
{
    bool result;
    int n;
    int i;
    double betak;
    double v;
    double vv;
    int mcinfo;
    bool b;
    bool stepfound;
    int diffcnt;

    
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
        diffcnt = state.rstate.ia(3);
        b = state.rstate.ba(0);
        stepfound = state.rstate.ba(1);
        betak = state.rstate.ra(0);
        v = state.rstate.ra(1);
        vv = state.rstate.ra(2);
    }
    else
    {
        n = -983;
        i = -989;
        mcinfo = -834;
        diffcnt = 900;
        b = true;
        stepfound = false;
        betak = 214;
        v = -338;
        vv = -686;
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
    if( state.rstate.stage==4 )
    {
        goto lbl_4;
    }
    if( state.rstate.stage==5 )
    {
        goto lbl_5;
    }
    if( state.rstate.stage==6 )
    {
        goto lbl_6;
    }
    if( state.rstate.stage==7 )
    {
        goto lbl_7;
    }
    if( state.rstate.stage==8 )
    {
        goto lbl_8;
    }
    if( state.rstate.stage==9 )
    {
        goto lbl_9;
    }
    if( state.rstate.stage==10 )
    {
        goto lbl_10;
    }
    if( state.rstate.stage==11 )
    {
        goto lbl_11;
    }
    if( state.rstate.stage==12 )
    {
        goto lbl_12;
    }
    if( state.rstate.stage==13 )
    {
        goto lbl_13;
    }
    if( state.rstate.stage==14 )
    {
        goto lbl_14;
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
    state.cgtype = 1;
    ap::vmove(&state.xk(0), 1, &state.x(0), 1, ap::vlen(0,n-1));
    for(i = 0; i <= n-1; i++)
    {
        if( ap::fp_eq(state.xk(i),state.bndl(i))||ap::fp_eq(state.xk(i),state.bndu(i)) )
        {
            state.ak(i) = 0;
        }
        else
        {
            state.ak(i) = 1;
        }
    }
    state.mu = 0.1;
    state.curalgo = 0;
    
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
        goto lbl_15;
    }
    
    //
    // progress report
    //
    clearrequestfields(state);
    state.xupdated = true;
    state.rstate.stage = 1;
    goto lbl_rcomm;
lbl_1:
lbl_15:
    if( ap::fp_less_eq(asaboundedantigradnorm(state),state.epsg) )
    {
        state.repterminationtype = 4;
        result = false;
        return result;
    }
    state.repnfev = state.repnfev+1;
    
    //
    // Main cycle
    //
    // At the beginning of new iteration:
    // * CurAlgo stores current algorithm selector
    // * State.XK, State.F and State.G store current X/F/G
    // * State.AK stores current set of active constraints
    //
lbl_17:
    if( false )
    {
        goto lbl_18;
    }
    
    //
    // GPA algorithm
    //
    if( state.curalgo!=0 )
    {
        goto lbl_19;
    }
    state.k = 0;
    state.acount = 0;
lbl_21:
    if( false )
    {
        goto lbl_22;
    }
    
    //
    // Determine Dk = proj(xk - gk)-xk
    //
    for(i = 0; i <= n-1; i++)
    {
        state.d(i) = asaboundval(state.xk(i)-state.g(i), state.bndl(i), state.bndu(i))-state.xk(i);
    }
    
    //
    // Armijo line search.
    // * exact search with alpha=1 is tried first,
    //   'exact' means that we evaluate f() EXACTLY at
    //   bound(x-g,bndl,bndu), without intermediate floating
    //   point operations.
    // * alpha<1 are tried if explicit search wasn't successful
    // Result is placed into XN.
    //
    // Two types of search are needed because we can't
    // just use second type with alpha=1 because in finite
    // precision arithmetics (x1-x0)+x0 may differ from x1.
    // So while x1 is correctly bounded (it lie EXACTLY on
    // boundary, if it is active), (x1-x0)+x0 may be
    // not bounded.
    //
    v = ap::vdotproduct(&state.d(0), 1, &state.g(0), 1, ap::vlen(0,n-1));
    state.dginit = v;
    state.finit = state.f;
    if( !(ap::fp_less_eq(asad1norm(state),state.stpmax)||ap::fp_eq(state.stpmax,0)) )
    {
        goto lbl_23;
    }
    
    //
    // Try alpha=1 step first
    //
    for(i = 0; i <= n-1; i++)
    {
        state.x(i) = asaboundval(state.xk(i)-state.g(i), state.bndl(i), state.bndu(i));
    }
    clearrequestfields(state);
    state.needfg = true;
    state.rstate.stage = 2;
    goto lbl_rcomm;
lbl_2:
    state.repnfev = state.repnfev+1;
    stepfound = ap::fp_less_eq(state.f,state.finit+gpaftol*state.dginit);
    goto lbl_24;
lbl_23:
    stepfound = false;
lbl_24:
    if( !stepfound )
    {
        goto lbl_25;
    }
    
    //
    // we are at the boundary(ies)
    //
    ap::vmove(&state.xn(0), 1, &state.x(0), 1, ap::vlen(0,n-1));
    state.stp = 1;
    goto lbl_26;
lbl_25:
    
    //
    // alpha=1 is too large, try smaller values
    //
    state.stp = 1;
    linminnormalized(state.d, state.stp, n);
    state.dginit = state.dginit/state.stp;
    state.stp = gpadecay*state.stp;
    if( ap::fp_greater(state.stpmax,0) )
    {
        state.stp = ap::minreal(state.stp, state.stpmax);
    }
lbl_27:
    if( false )
    {
        goto lbl_28;
    }
    v = state.stp;
    ap::vmove(&state.x(0), 1, &state.xk(0), 1, ap::vlen(0,n-1));
    ap::vadd(&state.x(0), 1, &state.d(0), 1, ap::vlen(0,n-1), v);
    clearrequestfields(state);
    state.needfg = true;
    state.rstate.stage = 3;
    goto lbl_rcomm;
lbl_3:
    state.repnfev = state.repnfev+1;
    if( ap::fp_less_eq(state.stp,stpmin) )
    {
        goto lbl_28;
    }
    if( ap::fp_less_eq(state.f,state.finit+state.stp*gpaftol*state.dginit) )
    {
        goto lbl_28;
    }
    state.stp = state.stp*gpadecay;
    goto lbl_27;
lbl_28:
    ap::vmove(&state.xn(0), 1, &state.x(0), 1, ap::vlen(0,n-1));
lbl_26:
    state.repiterationscount = state.repiterationscount+1;
    if( !state.xrep )
    {
        goto lbl_29;
    }
    
    //
    // progress report
    //
    clearrequestfields(state);
    state.xupdated = true;
    state.rstate.stage = 4;
    goto lbl_rcomm;
lbl_4:
lbl_29:
    
    //
    // Calculate new set of active constraints.
    // Reset counter if active set was changed.
    // Prepare for the new iteration
    //
    for(i = 0; i <= n-1; i++)
    {
        if( ap::fp_eq(state.xn(i),state.bndl(i))||ap::fp_eq(state.xn(i),state.bndu(i)) )
        {
            state.an(i) = 0;
        }
        else
        {
            state.an(i) = 1;
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        if( ap::fp_neq(state.ak(i),state.an(i)) )
        {
            state.acount = -1;
            break;
        }
    }
    state.acount = state.acount+1;
    ap::vmove(&state.xk(0), 1, &state.xn(0), 1, ap::vlen(0,n-1));
    ap::vmove(&state.ak(0), 1, &state.an(0), 1, ap::vlen(0,n-1));
    
    //
    // Stopping conditions
    //
    if( !(state.repiterationscount>=state.maxits&&state.maxits>0) )
    {
        goto lbl_31;
    }
    
    //
    // Too many iterations
    //
    state.repterminationtype = 5;
    if( !state.xrep )
    {
        goto lbl_33;
    }
    clearrequestfields(state);
    state.xupdated = true;
    state.rstate.stage = 5;
    goto lbl_rcomm;
lbl_5:
lbl_33:
    result = false;
    return result;
lbl_31:
    if( ap::fp_greater(asaboundedantigradnorm(state),state.epsg) )
    {
        goto lbl_35;
    }
    
    //
    // Gradient is small enough
    //
    state.repterminationtype = 4;
    if( !state.xrep )
    {
        goto lbl_37;
    }
    clearrequestfields(state);
    state.xupdated = true;
    state.rstate.stage = 6;
    goto lbl_rcomm;
lbl_6:
lbl_37:
    result = false;
    return result;
lbl_35:
    v = ap::vdotproduct(&state.d(0), 1, &state.d(0), 1, ap::vlen(0,n-1));
    if( ap::fp_greater(sqrt(v)*state.stp,state.epsx) )
    {
        goto lbl_39;
    }
    
    //
    // Step size is too small, no further improvement is
    // possible
    //
    state.repterminationtype = 2;
    if( !state.xrep )
    {
        goto lbl_41;
    }
    clearrequestfields(state);
    state.xupdated = true;
    state.rstate.stage = 7;
    goto lbl_rcomm;
lbl_7:
lbl_41:
    result = false;
    return result;
lbl_39:
    if( ap::fp_greater(state.finit-state.f,state.epsf*ap::maxreal(fabs(state.finit), ap::maxreal(fabs(state.f), 1.0))) )
    {
        goto lbl_43;
    }
    
    //
    // F(k+1)-F(k) is small enough
    //
    state.repterminationtype = 1;
    if( !state.xrep )
    {
        goto lbl_45;
    }
    clearrequestfields(state);
    state.xupdated = true;
    state.rstate.stage = 8;
    goto lbl_rcomm;
lbl_8:
lbl_45:
    result = false;
    return result;
lbl_43:
    
    //
    // Decide - should we switch algorithm or not
    //
    if( asauisempty(state) )
    {
        if( ap::fp_greater_eq(asaginorm(state),state.mu*asad1norm(state)) )
        {
            state.curalgo = 1;
            goto lbl_22;
        }
        else
        {
            state.mu = state.mu*asarho;
        }
    }
    else
    {
        if( state.acount==n1 )
        {
            if( ap::fp_greater_eq(asaginorm(state),state.mu*asad1norm(state)) )
            {
                state.curalgo = 1;
                goto lbl_22;
            }
        }
    }
    
    //
    // Next iteration
    //
    state.k = state.k+1;
    goto lbl_21;
lbl_22:
lbl_19:
    
    //
    // CG algorithm
    //
    if( state.curalgo!=1 )
    {
        goto lbl_47;
    }
    
    //
    // first, check that there are non-active constraints.
    // move to GPA algorithm, if all constraints are active
    //
    b = true;
    for(i = 0; i <= n-1; i++)
    {
        if( ap::fp_neq(state.ak(i),0) )
        {
            b = false;
            break;
        }
    }
    if( b )
    {
        state.curalgo = 0;
        goto lbl_17;
    }
    
    //
    // CG iterations
    //
    state.fold = state.f;
    ap::vmove(&state.xk(0), 1, &state.x(0), 1, ap::vlen(0,n-1));
    for(i = 0; i <= n-1; i++)
    {
        state.dk(i) = -state.g(i)*state.ak(i);
        state.gc(i) = state.g(i)*state.ak(i);
    }
lbl_49:
    if( false )
    {
        goto lbl_50;
    }
    
    //
    // Store G[k] for later calculation of Y[k]
    //
    for(i = 0; i <= n-1; i++)
    {
        state.yk(i) = -state.gc(i);
    }
    
    //
    // Make a CG step in direction given by DK[]:
    // * calculate step. Step projection into feasible set
    //   is used. It has several benefits: a) step may be
    //   found with usual line search, b) multiple constraints
    //   may be activated with one step, c) activated constraints
    //   are detected in a natural way - just compare x[i] with
    //   bounds
    // * update active set, set B to True, if there
    //   were changes in the set.
    //
    ap::vmove(&state.d(0), 1, &state.dk(0), 1, ap::vlen(0,n-1));
    ap::vmove(&state.xn(0), 1, &state.xk(0), 1, ap::vlen(0,n-1));
    state.mcstage = 0;
    state.stp = 1;
    linminnormalized(state.d, state.stp, n);
    mcsrch(n, state.xn, state.f, state.gc, state.d, state.stp, state.stpmax, mcinfo, state.nfev, state.work, state.lstate, state.mcstage);
lbl_51:
    if( state.mcstage==0 )
    {
        goto lbl_52;
    }
    
    //
    // preprocess data: bound State.XN so it belongs to the
    // feasible set and store it in the State.X
    //
    for(i = 0; i <= n-1; i++)
    {
        state.x(i) = asaboundval(state.xn(i), state.bndl(i), state.bndu(i));
    }
    
    //
    // RComm
    //
    clearrequestfields(state);
    state.needfg = true;
    state.rstate.stage = 9;
    goto lbl_rcomm;
lbl_9:
    
    //
    // postprocess data: zero components of G corresponding to
    // the active constraints
    //
    for(i = 0; i <= n-1; i++)
    {
        if( ap::fp_eq(state.x(i),state.bndl(i))||ap::fp_eq(state.x(i),state.bndu(i)) )
        {
            state.gc(i) = 0;
        }
        else
        {
            state.gc(i) = state.g(i);
        }
    }
    mcsrch(n, state.xn, state.f, state.gc, state.d, state.stp, state.stpmax, mcinfo, state.nfev, state.work, state.lstate, state.mcstage);
    goto lbl_51;
lbl_52:
    diffcnt = 0;
    for(i = 0; i <= n-1; i++)
    {
        
        //
        // XN contains unprojected result, project it,
        // save copy to X (will be used for progress reporting)
        //
        state.xn(i) = asaboundval(state.xn(i), state.bndl(i), state.bndu(i));
        
        //
        // update active set
        //
        if( ap::fp_eq(state.xn(i),state.bndl(i))||ap::fp_eq(state.xn(i),state.bndu(i)) )
        {
            state.an(i) = 0;
        }
        else
        {
            state.an(i) = 1;
        }
        if( ap::fp_neq(state.an(i),state.ak(i)) )
        {
            diffcnt = diffcnt+1;
        }
        state.ak(i) = state.an(i);
    }
    ap::vmove(&state.xk(0), 1, &state.xn(0), 1, ap::vlen(0,n-1));
    state.repnfev = state.repnfev+state.nfev;
    state.repiterationscount = state.repiterationscount+1;
    if( !state.xrep )
    {
        goto lbl_53;
    }
    
    //
    // progress report
    //
    clearrequestfields(state);
    state.xupdated = true;
    state.rstate.stage = 10;
    goto lbl_rcomm;
lbl_10:
lbl_53:
    
    //
    // Check stopping conditions.
    //
    if( ap::fp_greater(asaboundedantigradnorm(state),state.epsg) )
    {
        goto lbl_55;
    }
    
    //
    // Gradient is small enough
    //
    state.repterminationtype = 4;
    if( !state.xrep )
    {
        goto lbl_57;
    }
    clearrequestfields(state);
    state.xupdated = true;
    state.rstate.stage = 11;
    goto lbl_rcomm;
lbl_11:
lbl_57:
    result = false;
    return result;
lbl_55:
    if( !(state.repiterationscount>=state.maxits&&state.maxits>0) )
    {
        goto lbl_59;
    }
    
    //
    // Too many iterations
    //
    state.repterminationtype = 5;
    if( !state.xrep )
    {
        goto lbl_61;
    }
    clearrequestfields(state);
    state.xupdated = true;
    state.rstate.stage = 12;
    goto lbl_rcomm;
lbl_12:
lbl_61:
    result = false;
    return result;
lbl_59:
    if( !(ap::fp_greater_eq(asaginorm(state),state.mu*asad1norm(state))&&diffcnt==0) )
    {
        goto lbl_63;
    }
    
    //
    // These conditions are explicitly or implicitly
    // related to the current step size and influenced
    // by changes in the active constraints.
    //
    // For these reasons they are checked only when we don't
    // want to 'unstick' at the end of the iteration and there
    // were no changes in the active set.
    //
    // NOTE: consition |G|>=Mu*|D1| must be exactly opposite
    // to the condition used to switch back to GPA. At least
    // one inequality must be strict, otherwise infinite cycle
    // may occur when |G|=Mu*|D1| (we DON'T test stopping
    // conditions and we DON'T switch to GPA, so we cycle
    // indefinitely).
    //
    if( ap::fp_greater(state.fold-state.f,state.epsf*ap::maxreal(fabs(state.fold), ap::maxreal(fabs(state.f), 1.0))) )
    {
        goto lbl_65;
    }
    
    //
    // F(k+1)-F(k) is small enough
    //
    state.repterminationtype = 1;
    if( !state.xrep )
    {
        goto lbl_67;
    }
    clearrequestfields(state);
    state.xupdated = true;
    state.rstate.stage = 13;
    goto lbl_rcomm;
lbl_13:
lbl_67:
    result = false;
    return result;
lbl_65:
    v = ap::vdotproduct(&state.d(0), 1, &state.d(0), 1, ap::vlen(0,n-1));
    if( ap::fp_greater(sqrt(v)*state.stp,state.epsx) )
    {
        goto lbl_69;
    }
    
    //
    // X(k+1)-X(k) is small enough
    //
    state.repterminationtype = 2;
    if( !state.xrep )
    {
        goto lbl_71;
    }
    clearrequestfields(state);
    state.xupdated = true;
    state.rstate.stage = 14;
    goto lbl_rcomm;
lbl_14:
lbl_71:
    result = false;
    return result;
lbl_69:
lbl_63:
    
    //
    // Check conditions for switching
    //
    if( ap::fp_less(asaginorm(state),state.mu*asad1norm(state)) )
    {
        state.curalgo = 0;
        goto lbl_50;
    }
    if( diffcnt>0 )
    {
        if( asauisempty(state)||diffcnt>=n2 )
        {
            state.curalgo = 1;
        }
        else
        {
            state.curalgo = 0;
        }
        goto lbl_50;
    }
    
    //
    // Calculate D(k+1)
    //
    // Line search may result in:
    // * maximum feasible step being taken (already processed)
    // * point satisfying Wolfe conditions
    // * some kind of error (CG is restarted by assigning 0.0 to Beta)
    //
    if( mcinfo==1 )
    {
        
        //
        // Standard Wolfe conditions are satisfied:
        // * calculate Y[K] and BetaK
        //
        ap::vadd(&state.yk(0), 1, &state.gc(0), 1, ap::vlen(0,n-1));
        vv = ap::vdotproduct(&state.yk(0), 1, &state.dk(0), 1, ap::vlen(0,n-1));
        v = ap::vdotproduct(&state.gc(0), 1, &state.gc(0), 1, ap::vlen(0,n-1));
        state.betady = v/vv;
        v = ap::vdotproduct(&state.gc(0), 1, &state.yk(0), 1, ap::vlen(0,n-1));
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
    ap::vmoveneg(&state.dn(0), 1, &state.gc(0), 1, ap::vlen(0,n-1));
    ap::vadd(&state.dn(0), 1, &state.dk(0), 1, ap::vlen(0,n-1), betak);
    ap::vmove(&state.dk(0), 1, &state.dn(0), 1, ap::vlen(0,n-1));
    
    //
    // update other information
    //
    state.fold = state.f;
    state.k = state.k+1;
    goto lbl_49;
lbl_50:
lbl_47:
    goto lbl_17;
lbl_18:
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
    state.rstate.ia(3) = diffcnt;
    state.rstate.ba(0) = b;
    state.rstate.ba(1) = stepfound;
    state.rstate.ra(0) = betak;
    state.rstate.ra(1) = v;
    state.rstate.ra(2) = vv;
    return result;
}


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
     minasareport& rep)
{
    int i;

    x.setbounds(0, state.n-1);
    ap::vmove(&x(0), 1, &state.x(0), 1, ap::vlen(0,state.n-1));
    rep.iterationscount = state.repiterationscount;
    rep.nfev = state.repnfev;
    rep.terminationtype = state.repterminationtype;
    rep.activeconstraints = 0;
    for(i = 0; i <= state.n-1; i++)
    {
        if( ap::fp_eq(state.ak(i),0) )
        {
            rep.activeconstraints = rep.activeconstraints+1;
        }
    }
}


/*************************************************************************
'bound' value: map X to [B1,B2]

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
static double asaboundval(double x, double b1, double b2)
{
    double result;

    if( ap::fp_less_eq(x,b1) )
    {
        result = b1;
        return result;
    }
    if( ap::fp_greater_eq(x,b2) )
    {
        result = b2;
        return result;
    }
    result = x;
    return result;
}


/*************************************************************************
Returns norm of bounded anti-gradient.

Bounded antigradient is a vector obtained from  anti-gradient  by  zeroing
components which point outwards:
    result = norm(v)
    v[i]=0     if ((-g[i]<0)and(x[i]=bndl[i])) or
                  ((-g[i]>0)and(x[i]=bndu[i]))
    v[i]=-g[i] otherwise

This function may be used to check a stopping criterion.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
static double asaboundedantigradnorm(const minasastate& state)
{
    double result;
    int i;
    double v;

    result = 0;
    for(i = 0; i <= state.n-1; i++)
    {
        v = -state.g(i);
        if( ap::fp_eq(state.x(i),state.bndl(i))&&ap::fp_less(-state.g(i),0) )
        {
            v = 0;
        }
        if( ap::fp_eq(state.x(i),state.bndu(i))&&ap::fp_greater(-state.g(i),0) )
        {
            v = 0;
        }
        result = result+ap::sqr(v);
    }
    result = sqrt(result);
    return result;
}


/*************************************************************************
Returns norm of GI(x).

GI(x) is  a  gradient  vector  whose  components  associated  with  active
constraints are zeroed. It  differs  from  bounded  anti-gradient  because
components  of   GI(x)   are   zeroed  independently  of  sign(g[i]),  and
anti-gradient's components are zeroed with respect to both constraint  and
sign.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
static double asaginorm(const minasastate& state)
{
    double result;
    int i;
    double v;

    result = 0;
    for(i = 0; i <= state.n-1; i++)
    {
        if( ap::fp_neq(state.x(i),state.bndl(i))&&ap::fp_neq(state.x(i),state.bndu(i)) )
        {
            result = result+ap::sqr(state.g(i));
        }
    }
    result = sqrt(result);
    return result;
}


/*************************************************************************
Returns norm(D1(State.X))

For a meaning of D1 see 'NEW ACTIVE SET ALGORITHM FOR BOX CONSTRAINED
OPTIMIZATION' by WILLIAM W. HAGER AND HONGCHAO ZHANG.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
static double asad1norm(const minasastate& state)
{
    double result;
    int i;

    result = 0;
    for(i = 0; i <= state.n-1; i++)
    {
        result = result+ap::sqr(asaboundval(state.x(i)-state.g(i), state.bndl(i), state.bndu(i))-state.x(i));
    }
    result = sqrt(result);
    return result;
}


/*************************************************************************
Returns True, if U set is empty.

* State.X is used as point,
* State.G - as gradient,
* D is calculated within function (because State.D may have different
  meaning depending on current optimization algorithm)

For a meaning of U see 'NEW ACTIVE SET ALGORITHM FOR BOX CONSTRAINED
OPTIMIZATION' by WILLIAM W. HAGER AND HONGCHAO ZHANG.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
static bool asauisempty(const minasastate& state)
{
    bool result;
    int i;
    double d;
    double d2;
    double d32;

    d = asad1norm(state);
    d2 = sqrt(d);
    d32 = d*d2;
    result = true;
    for(i = 0; i <= state.n-1; i++)
    {
        if( ap::fp_greater_eq(fabs(state.g(i)),d2)&&ap::fp_greater_eq(ap::minreal(state.x(i)-state.bndl(i), state.bndu(i)-state.x(i)),d32) )
        {
            result = false;
            return result;
        }
    }
    return result;
}


/*************************************************************************
Returns True, if optimizer "want  to  unstick"  from  one  of  the  active
constraints, i.e. there is such active constraint with index I that either
lower bound is active and g[i]<0, or upper bound is active and g[i]>0.

State.X is used as current point, State.X - as gradient.
  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
static bool asawanttounstick(const minasastate& state)
{
    bool result;
    int i;

    result = false;
    for(i = 0; i <= state.n-1; i++)
    {
        if( ap::fp_eq(state.x(i),state.bndl(i))&&ap::fp_less(state.g(i),0) )
        {
            result = true;
        }
        if( ap::fp_eq(state.x(i),state.bndu(i))&&ap::fp_greater(state.g(i),0) )
        {
            result = true;
        }
        if( result )
        {
            return result;
        }
    }
    return result;
}


/*************************************************************************
Clears request fileds (to be sure that we don't forgot to clear something)
*************************************************************************/
static void clearrequestfields(minasastate& state)
{

    state.needfg = false;
    state.xupdated = false;
}




