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

#include <stdafx.h>
#include "minlbfgs.h"

static void clearrequestfields(minlbfgsstate& state);

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
     minlbfgsstate& state)
{

    minlbfgscreatex(n, m, x, 0, state);
}


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
     int maxits)
{

    ap::ap_error::make_assertion(ap::fp_greater_eq(epsg,0), "MinLBFGSSetCond: negative EpsG!");
    ap::ap_error::make_assertion(ap::fp_greater_eq(epsf,0), "MinLBFGSSetCond: negative EpsF!");
    ap::ap_error::make_assertion(ap::fp_greater_eq(epsx,0), "MinLBFGSSetCond: negative EpsX!");
    ap::ap_error::make_assertion(maxits>=0, "MinLBFGSSetCond: negative MaxIts!");
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
void minlbfgssetxrep(minlbfgsstate& state, bool needxrep)
{

    state.xrep = needxrep;
}


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
void minlbfgssetstpmax(minlbfgsstate& state, double stpmax)
{

    ap::ap_error::make_assertion(ap::fp_greater_eq(stpmax,0), "MinLBFGSSetStpMax: StpMax<0!");
    state.stpmax = stpmax;
}


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
     minlbfgsstate& state)
{
    bool allocatemem;

    ap::ap_error::make_assertion(n>=1, "MinLBFGS: N too small!");
    ap::ap_error::make_assertion(m>=1, "MinLBFGS: M too small!");
    ap::ap_error::make_assertion(m<=n, "MinLBFGS: M too large!");
    
    //
    // Initialize
    //
    state.n = n;
    state.m = m;
    state.flags = flags;
    allocatemem = flags%2==0;
    flags = flags/2;
    if( allocatemem )
    {
        state.rho.setbounds(0, m-1);
        state.theta.setbounds(0, m-1);
        state.y.setbounds(0, m-1, 0, n-1);
        state.s.setbounds(0, m-1, 0, n-1);
        state.d.setbounds(0, n-1);
        state.x.setbounds(0, n-1);
        state.g.setbounds(0, n-1);
        state.work.setbounds(0, n-1);
    }
    minlbfgssetcond(state, double(0), double(0), double(0), 0);
    minlbfgssetxrep(state, false);
    minlbfgssetstpmax(state, double(0));
    
    //
    // Prepare first run
    //
    state.k = 0;
    ap::vmove(&state.x(0), 1, &x(0), 1, ap::vlen(0,n-1));
    state.rstate.ia.setbounds(0, 6);
    state.rstate.ra.setbounds(0, 4);
    state.rstate.stage = -1;
}


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
bool minlbfgsiteration(minlbfgsstate& state)
{
    bool result;
    int n;
    int m;
    int maxits;
    double epsf;
    double epsg;
    double epsx;
    int i;
    int j;
    int ic;
    int mcinfo;
    double v;
    double vv;

    
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
        m = state.rstate.ia(1);
        maxits = state.rstate.ia(2);
        i = state.rstate.ia(3);
        j = state.rstate.ia(4);
        ic = state.rstate.ia(5);
        mcinfo = state.rstate.ia(6);
        epsf = state.rstate.ra(0);
        epsg = state.rstate.ra(1);
        epsx = state.rstate.ra(2);
        v = state.rstate.ra(3);
        vv = state.rstate.ra(4);
    }
    else
    {
        n = -983;
        m = -989;
        maxits = -834;
        i = 900;
        j = -287;
        ic = 364;
        mcinfo = 214;
        epsf = -338;
        epsg = -686;
        epsx = 912;
        v = 585;
        vv = 497;
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
    // Unload frequently used variables from State structure
    // (just for typing convinience)
    //
    n = state.n;
    m = state.m;
    epsg = state.epsg;
    epsf = state.epsf;
    epsx = state.epsx;
    maxits = state.maxits;
    state.repterminationtype = 0;
    state.repiterationscount = 0;
    state.repnfev = 0;
    
    //
    // Calculate F/G at the initial point
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
    state.repnfev = 1;
    state.fold = state.f;
    v = ap::vdotproduct(&state.g(0), 1, &state.g(0), 1, ap::vlen(0,n-1));
    v = sqrt(v);
    if( ap::fp_less_eq(v,epsg) )
    {
        state.repterminationtype = 4;
        result = false;
        return result;
    }
    
    //
    // Choose initial step
    //
    if( ap::fp_eq(state.stpmax,0) )
    {
        state.stp = ap::minreal(1.0/v, double(1));
    }
    else
    {
        state.stp = ap::minreal(1.0/v, state.stpmax);
    }
    ap::vmoveneg(&state.d(0), 1, &state.g(0), 1, ap::vlen(0,n-1));
    
    //
    // Main cycle
    //
lbl_6:
    if( false )
    {
        goto lbl_7;
    }
    
    //
    // Main cycle: prepare to 1-D line search
    //
    state.p = state.k%m;
    state.q = ap::minint(state.k, m-1);
    
    //
    // Store X[k], G[k]
    //
    ap::vmoveneg(&state.s(state.p, 0), 1, &state.x(0), 1, ap::vlen(0,n-1));
    ap::vmoveneg(&state.y(state.p, 0), 1, &state.g(0), 1, ap::vlen(0,n-1));
    
    //
    // Minimize F(x+alpha*d)
    // Calculate S[k], Y[k]
    //
    state.mcstage = 0;
    if( state.k!=0 )
    {
        state.stp = 1.0;
    }
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
    
    //
    // report
    //
    clearrequestfields(state);
    state.xupdated = true;
    state.rstate.stage = 3;
    goto lbl_rcomm;
lbl_3:
lbl_10:
    state.repnfev = state.repnfev+state.nfev;
    state.repiterationscount = state.repiterationscount+1;
    ap::vadd(&state.s(state.p, 0), 1, &state.x(0), 1, ap::vlen(0,n-1));
    ap::vadd(&state.y(state.p, 0), 1, &state.g(0), 1, ap::vlen(0,n-1));
    
    //
    // Stopping conditions
    //
    if( state.repiterationscount>=maxits&&maxits>0 )
    {
        
        //
        // Too many iterations
        //
        state.repterminationtype = 5;
        result = false;
        return result;
    }
    v = ap::vdotproduct(&state.g(0), 1, &state.g(0), 1, ap::vlen(0,n-1));
    if( ap::fp_less_eq(sqrt(v),epsg) )
    {
        
        //
        // Gradient is small enough
        //
        state.repterminationtype = 4;
        result = false;
        return result;
    }
    if( ap::fp_less_eq(state.fold-state.f,epsf*ap::maxreal(fabs(state.fold), ap::maxreal(fabs(state.f), 1.0))) )
    {
        
        //
        // F(k+1)-F(k) is small enough
        //
        state.repterminationtype = 1;
        result = false;
        return result;
    }
    v = ap::vdotproduct(&state.s(state.p, 0), 1, &state.s(state.p, 0), 1, ap::vlen(0,n-1));
    if( ap::fp_less_eq(sqrt(v),epsx) )
    {
        
        //
        // X(k+1)-X(k) is small enough
        //
        state.repterminationtype = 2;
        result = false;
        return result;
    }
    
    //
    // If Wolfe conditions are satisfied, we can update
    // limited memory model.
    //
    // However, if conditions are not satisfied (NFEV limit is met,
    // function is too wild, ...), we'll skip L-BFGS update
    //
    if( mcinfo!=1 )
    {
        
        //
        // Skip update.
        //
        // In such cases we'll initialize search direction by
        // antigradient vector, because it  leads to more
        // transparent code with less number of special cases
        //
        state.fold = state.f;
        ap::vmoveneg(&state.d(0), 1, &state.g(0), 1, ap::vlen(0,n-1));
    }
    else
    {
        
        //
        // Calculate Rho[k], GammaK
        //
        v = ap::vdotproduct(&state.y(state.p, 0), 1, &state.s(state.p, 0), 1, ap::vlen(0,n-1));
        vv = ap::vdotproduct(&state.y(state.p, 0), 1, &state.y(state.p, 0), 1, ap::vlen(0,n-1));
        if( ap::fp_eq(v,0)||ap::fp_eq(vv,0) )
        {
            
            //
            // Rounding errors make further iterations impossible.
            //
            state.repterminationtype = -2;
            result = false;
            return result;
        }
        state.rho(state.p) = 1/v;
        state.gammak = v/vv;
        
        //
        //  Calculate d(k+1) = -H(k+1)*g(k+1)
        //
        //  for I:=K downto K-Q do
        //      V = s(i)^T * work(iteration:I)
        //      theta(i) = V
        //      work(iteration:I+1) = work(iteration:I) - V*Rho(i)*y(i)
        //  work(last iteration) = H0*work(last iteration)
        //  for I:=K-Q to K do
        //      V = y(i)^T*work(iteration:I)
        //      work(iteration:I+1) = work(iteration:I) +(-V+theta(i))*Rho(i)*s(i)
        //
        //  NOW WORK CONTAINS d(k+1)
        //
        ap::vmove(&state.work(0), 1, &state.g(0), 1, ap::vlen(0,n-1));
        for(i = state.k; i >= state.k-state.q; i--)
        {
            ic = i%m;
            v = ap::vdotproduct(&state.s(ic, 0), 1, &state.work(0), 1, ap::vlen(0,n-1));
            state.theta(ic) = v;
            vv = v*state.rho(ic);
            ap::vsub(&state.work(0), 1, &state.y(ic, 0), 1, ap::vlen(0,n-1), vv);
        }
        v = state.gammak;
        ap::vmul(&state.work(0), 1, ap::vlen(0,n-1), v);
        for(i = state.k-state.q; i <= state.k; i++)
        {
            ic = i%m;
            v = ap::vdotproduct(&state.y(ic, 0), 1, &state.work(0), 1, ap::vlen(0,n-1));
            vv = state.rho(ic)*(-v+state.theta(ic));
            ap::vadd(&state.work(0), 1, &state.s(ic, 0), 1, ap::vlen(0,n-1), vv);
        }
        ap::vmoveneg(&state.d(0), 1, &state.work(0), 1, ap::vlen(0,n-1));
        
        //
        // Next step
        //
        state.fold = state.f;
        state.k = state.k+1;
    }
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
    state.rstate.ia(1) = m;
    state.rstate.ia(2) = maxits;
    state.rstate.ia(3) = i;
    state.rstate.ia(4) = j;
    state.rstate.ia(5) = ic;
    state.rstate.ia(6) = mcinfo;
    state.rstate.ra(0) = epsf;
    state.rstate.ra(1) = epsg;
    state.rstate.ra(2) = epsx;
    state.rstate.ra(3) = v;
    state.rstate.ra(4) = vv;
    return result;
}


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
     minlbfgsreport& rep)
{

    x.setbounds(0, state.n-1);
    ap::vmove(&x(0), 1, &state.x(0), 1, ap::vlen(0,state.n-1));
    rep.iterationscount = state.repiterationscount;
    rep.nfev = state.repnfev;
    rep.terminationtype = state.repterminationtype;
}


/*************************************************************************
Clears request fileds (to be sure that we don't forgot to clear something)
*************************************************************************/
static void clearrequestfields(minlbfgsstate& state)
{

    state.needfg = false;
    state.xupdated = false;
}




