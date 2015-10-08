/*************************************************************************
Copyright 2009 by Sergey Bochkanov (ALGLIB project).

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
#include "odesolver.h"

static const double odesolvermaxgrow = 3.0;
static const double odesolvermaxshrink = 10.0;

static void odesolverinit(int solvertype,
     const ap::real_1d_array& y,
     int n,
     const ap::real_1d_array& x,
     int m,
     double eps,
     double h,
     odesolverstate& state);

/*************************************************************************
Cash-Karp adaptive ODE solver.

This subroutine solves ODE  Y'=f(Y,x)  with  initial  conditions  Y(xs)=Ys
(here Y may be single variable or vector of N variables).

INPUT PARAMETERS:
    Y       -   initial conditions, array[0..N-1].
                contains values of Y[] at X[0]
    N       -   system size
    X       -   points at which Y should be tabulated, array[0..M-1]
                integrations starts at X[0], ends at X[M-1],  intermediate
                values at X[i] are returned too.
                SHOULD BE ORDERED BY ASCENDING OR BY DESCENDING!!!!
    M       -   number of intermediate points + first point + last point:
                * M>2 means that you need both Y(X[M-1]) and M-2 values at
                  intermediate points
                * M=2 means that you want just to integrate from  X[0]  to
                  X[1] and don't interested in intermediate values.
                * M=1 means that you don't want to integrate :)
                  it is degenerate case, but it will be handled correctly.
                * M<1 means error
    Eps     -   tolerance (absolute/relative error on each  step  will  be
                less than Eps). When passing:
                * Eps>0, it means desired ABSOLUTE error
                * Eps<0, it means desired RELATIVE error.  Relative errors
                  are calculated with respect to maximum values of  Y seen
                  so far. Be careful to use this criterion  when  starting
                  from Y[] that are close to zero.
    H       -   initial  step  lenth,  it  will  be adjusted automatically
                after the first  step.  If  H=0,  step  will  be  selected
                automatically  (usualy  it  will  be  equal  to  0.001  of
                min(x[i]-x[j])).

OUTPUT PARAMETERS
    State   -   structure which stores algorithm state between  subsequent
                calls of OdeSolverIteration. Used for reverse communication.
                This structure should be passed  to the OdeSolverIteration
                subroutine.

SEE ALSO
    AutoGKSmoothW, AutoGKSingular, AutoGKIteration, AutoGKResults.


  -- ALGLIB --
     Copyright 01.09.2009 by Bochkanov Sergey
*************************************************************************/
void odesolverrkck(const ap::real_1d_array& y,
     int n,
     const ap::real_1d_array& x,
     int m,
     double eps,
     double h,
     odesolverstate& state)
{

    odesolverinit(0, y, n, x, m, eps, h, state);
}


/*************************************************************************
One iteration of ODE solver.

Called after inialization of State structure with OdeSolverXXX subroutine.
See HTML docs for examples.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between subsequent
                calls and which is used for reverse communication. Must be
                initialized with OdeSolverXXX() call first.

If subroutine returned False, algorithm have finished its work.
If subroutine returned True, then user should:
* calculate F(State.X, State.Y)
* store it in State.DY
Here State.X is real, State.Y and State.DY are arrays[0..N-1] of reals.

  -- ALGLIB --
     Copyright 01.09.2009 by Bochkanov Sergey
*************************************************************************/
bool odesolveriteration(odesolverstate& state)
{
    bool result;
    int n;
    int m;
    int i;
    int j;
    int k;
    double xc;
    double v;
    double h;
    double h2;
    bool gridpoint;
    double err;
    double maxgrowpow;
    int klimit;

    
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
        i = state.rstate.ia(2);
        j = state.rstate.ia(3);
        k = state.rstate.ia(4);
        klimit = state.rstate.ia(5);
        gridpoint = state.rstate.ba(0);
        xc = state.rstate.ra(0);
        v = state.rstate.ra(1);
        h = state.rstate.ra(2);
        h2 = state.rstate.ra(3);
        err = state.rstate.ra(4);
        maxgrowpow = state.rstate.ra(5);
    }
    else
    {
        n = -983;
        m = -989;
        i = -834;
        j = 900;
        k = -287;
        klimit = 364;
        gridpoint = false;
        xc = -338;
        v = -686;
        h = 912;
        h2 = 585;
        err = 497;
        maxgrowpow = -271;
    }
    if( state.rstate.stage==0 )
    {
        goto lbl_0;
    }
    
    //
    // Routine body
    //
    
    //
    // prepare
    //
    if( state.repterminationtype!=0 )
    {
        result = false;
        return result;
    }
    n = state.n;
    m = state.m;
    h = state.h;
    state.y.setlength(n);
    state.dy.setlength(n);
    maxgrowpow = pow(odesolvermaxgrow, double(5));
    state.repnfev = 0;
    
    //
    // some preliminary checks for internal errors
    // after this we assume that H>0 and M>1
    //
    ap::ap_error::make_assertion(ap::fp_greater(state.h,0), "ODESolver: internal error");
    ap::ap_error::make_assertion(m>1, "ODESolverIteration: internal error");
    
    //
    // choose solver
    //
    if( state.solvertype!=0 )
    {
        goto lbl_1;
    }
    
    //
    // Cask-Karp solver
    // Prepare coefficients table.
    // Check it for errors
    //
    state.rka.setlength(6);
    state.rka(0) = 0;
    state.rka(1) = double(1)/double(5);
    state.rka(2) = double(3)/double(10);
    state.rka(3) = double(3)/double(5);
    state.rka(4) = 1;
    state.rka(5) = double(7)/double(8);
    state.rkb.setlength(6, 5);
    state.rkb(1,0) = double(1)/double(5);
    state.rkb(2,0) = double(3)/double(40);
    state.rkb(2,1) = double(9)/double(40);
    state.rkb(3,0) = double(3)/double(10);
    state.rkb(3,1) = -double(9)/double(10);
    state.rkb(3,2) = double(6)/double(5);
    state.rkb(4,0) = -double(11)/double(54);
    state.rkb(4,1) = double(5)/double(2);
    state.rkb(4,2) = -double(70)/double(27);
    state.rkb(4,3) = double(35)/double(27);
    state.rkb(5,0) = double(1631)/double(55296);
    state.rkb(5,1) = double(175)/double(512);
    state.rkb(5,2) = double(575)/double(13824);
    state.rkb(5,3) = double(44275)/double(110592);
    state.rkb(5,4) = double(253)/double(4096);
    state.rkc.setlength(6);
    state.rkc(0) = double(37)/double(378);
    state.rkc(1) = 0;
    state.rkc(2) = double(250)/double(621);
    state.rkc(3) = double(125)/double(594);
    state.rkc(4) = 0;
    state.rkc(5) = double(512)/double(1771);
    state.rkcs.setlength(6);
    state.rkcs(0) = double(2825)/double(27648);
    state.rkcs(1) = 0;
    state.rkcs(2) = double(18575)/double(48384);
    state.rkcs(3) = double(13525)/double(55296);
    state.rkcs(4) = double(277)/double(14336);
    state.rkcs(5) = double(1)/double(4);
    state.rkk.setlength(6, n);
    
    //
    // Main cycle consists of two iterations:
    // * outer where we travel from X[i-1] to X[i]
    // * inner where we travel inside [X[i-1],X[i]]
    //
    state.ytbl.setlength(m, n);
    state.escale.setlength(n);
    state.yn.setlength(n);
    state.yns.setlength(n);
    xc = state.xg(0);
    ap::vmove(&state.ytbl(0, 0), 1, &state.yc(0), 1, ap::vlen(0,n-1));
    for(j = 0; j <= n-1; j++)
    {
        state.escale(j) = 0;
    }
    i = 1;
lbl_3:
    if( i>m-1 )
    {
        goto lbl_5;
    }
    
    //
    // begin inner iteration
    //
lbl_6:
    if( false )
    {
        goto lbl_7;
    }
    
    //
    // truncate step if needed (beyond right boundary).
    // determine should we store X or not
    //
    if( ap::fp_greater_eq(xc+h,state.xg(i)) )
    {
        h = state.xg(i)-xc;
        gridpoint = true;
    }
    else
    {
        gridpoint = false;
    }
    
    //
    // Update error scale maximums
    //
    // These maximums are initialized by zeros,
    // then updated every iterations.
    //
    for(j = 0; j <= n-1; j++)
    {
        state.escale(j) = ap::maxreal(state.escale(j), fabs(state.yc(j)));
    }
    
    //
    // make one step:
    // 1. calculate all info needed to do step
    // 2. update errors scale maximums using values/derivatives
    //    obtained during (1)
    //
    // Take into account that we use scaling of X to reduce task
    // to the form where x[0] < x[1] < ... < x[n-1]. So X is
    // replaced by x=xscale*t, and dy/dx=f(y,x) is replaced
    // by dy/dt=xscale*f(y,xscale*t).
    //
    ap::vmove(&state.yn(0), 1, &state.yc(0), 1, ap::vlen(0,n-1));
    ap::vmove(&state.yns(0), 1, &state.yc(0), 1, ap::vlen(0,n-1));
    k = 0;
lbl_8:
    if( k>5 )
    {
        goto lbl_10;
    }
    
    //
    // prepare data for the next update of YN/YNS
    //
    state.x = state.xscale*(xc+state.rka(k)*h);
    ap::vmove(&state.y(0), 1, &state.yc(0), 1, ap::vlen(0,n-1));
    for(j = 0; j <= k-1; j++)
    {
        v = state.rkb(k,j);
        ap::vadd(&state.y(0), 1, &state.rkk(j, 0), 1, ap::vlen(0,n-1), v);
    }
    state.rstate.stage = 0;
    goto lbl_rcomm;
lbl_0:
    state.repnfev = state.repnfev+1;
    v = h*state.xscale;
    ap::vmove(&state.rkk(k, 0), 1, &state.dy(0), 1, ap::vlen(0,n-1), v);
    
    //
    // update YN/YNS
    //
    v = state.rkc(k);
    ap::vadd(&state.yn(0), 1, &state.rkk(k, 0), 1, ap::vlen(0,n-1), v);
    v = state.rkcs(k);
    ap::vadd(&state.yns(0), 1, &state.rkk(k, 0), 1, ap::vlen(0,n-1), v);
    k = k+1;
    goto lbl_8;
lbl_10:
    
    //
    // estimate error
    //
    err = 0;
    for(j = 0; j <= n-1; j++)
    {
        if( !state.fraceps )
        {
            
            //
            // absolute error is estimated
            //
            err = ap::maxreal(err, fabs(state.yn(j)-state.yns(j)));
        }
        else
        {
            
            //
            // Relative error is estimated
            //
            v = state.escale(j);
            if( ap::fp_eq(v,0) )
            {
                v = 1;
            }
            err = ap::maxreal(err, fabs(state.yn(j)-state.yns(j))/v);
        }
    }
    
    //
    // calculate new step, restart if necessary
    //
    if( ap::fp_less_eq(maxgrowpow*err,state.eps) )
    {
        h2 = odesolvermaxgrow*h;
    }
    else
    {
        h2 = h*pow(state.eps/err, 0.2);
    }
    if( ap::fp_less(h2,h/odesolvermaxshrink) )
    {
        h2 = h/odesolvermaxshrink;
    }
    if( ap::fp_greater(err,state.eps) )
    {
        h = h2;
        goto lbl_6;
    }
    
    //
    // advance position
    //
    xc = xc+h;
    ap::vmove(&state.yc(0), 1, &state.yn(0), 1, ap::vlen(0,n-1));
    
    //
    // update H
    //
    h = h2;
    
    //
    // break on grid point
    //
    if( gridpoint )
    {
        goto lbl_7;
    }
    goto lbl_6;
lbl_7:
    
    //
    // save result
    //
    ap::vmove(&state.ytbl(i, 0), 1, &state.yc(0), 1, ap::vlen(0,n-1));
    i = i+1;
    goto lbl_3;
lbl_5:
    state.repterminationtype = 1;
    result = false;
    return result;
lbl_1:
    result = false;
    return result;
    
    //
    // Saving state
    //
lbl_rcomm:
    result = true;
    state.rstate.ia(0) = n;
    state.rstate.ia(1) = m;
    state.rstate.ia(2) = i;
    state.rstate.ia(3) = j;
    state.rstate.ia(4) = k;
    state.rstate.ia(5) = klimit;
    state.rstate.ba(0) = gridpoint;
    state.rstate.ra(0) = xc;
    state.rstate.ra(1) = v;
    state.rstate.ra(2) = h;
    state.rstate.ra(3) = h2;
    state.rstate.ra(4) = err;
    state.rstate.ra(5) = maxgrowpow;
    return result;
}


/*************************************************************************
ODE solver results

Called after OdeSolverIteration returned False.

INPUT PARAMETERS:
    State   -   algorithm state (used by OdeSolverIteration).

OUTPUT PARAMETERS:
    M       -   number of tabulated values, M>=1
    XTbl    -   array[0..M-1], values of X
    YTbl    -   array[0..M-1,0..N-1], values of Y in X[i]
    Rep     -   solver report:
                * Rep.TerminationType completetion code:
                    * -2    X is not ordered  by  ascending/descending  or
                            there are non-distinct X[],  i.e.  X[i]=X[i+1]
                    * -1    incorrect parameters were specified
                    *  1    task has been solved
                * Rep.NFEV contains number of function calculations

  -- ALGLIB --
     Copyright 01.09.2009 by Bochkanov Sergey
*************************************************************************/
void odesolverresults(const odesolverstate& state,
     int& m,
     ap::real_1d_array& xtbl,
     ap::real_2d_array& ytbl,
     odesolverreport& rep)
{
    double v;
    int i;

    rep.terminationtype = state.repterminationtype;
    if( rep.terminationtype>0 )
    {
        m = state.m;
        rep.nfev = state.repnfev;
        xtbl.setlength(state.m);
        v = state.xscale;
        ap::vmove(&xtbl(0), 1, &state.xg(0), 1, ap::vlen(0,state.m-1), v);
        ytbl.setlength(state.m, state.n);
        for(i = 0; i <= state.m-1; i++)
        {
            ap::vmove(&ytbl(i, 0), 1, &state.ytbl(i, 0), 1, ap::vlen(0,state.n-1));
        }
    }
    else
    {
        rep.nfev = 0;
    }
}


/*************************************************************************
Internal initialization subroutine
*************************************************************************/
static void odesolverinit(int solvertype,
     const ap::real_1d_array& y,
     int n,
     const ap::real_1d_array& x,
     int m,
     double eps,
     double h,
     odesolverstate& state)
{
    int i;
    double v;

    
    //
    // Prepare RComm
    //
    state.rstate.ia.setbounds(0, 5);
    state.rstate.ba.setbounds(0, 0);
    state.rstate.ra.setbounds(0, 5);
    state.rstate.stage = -1;
    
    //
    // check parameters.
    //
    if( n<=0||m<1||ap::fp_eq(eps,0) )
    {
        state.repterminationtype = -1;
        return;
    }
    if( ap::fp_less(h,0) )
    {
        h = -h;
    }
    
    //
    // quick exit if necessary.
    // after this block we assume that M>1
    //
    if( m==1 )
    {
        state.repnfev = 0;
        state.repterminationtype = 1;
        state.ytbl.setlength(1, n);
        ap::vmove(&state.ytbl(0, 0), 1, &y(0), 1, ap::vlen(0,n-1));
        state.xg.setlength(m);
        ap::vmove(&state.xg(0), 1, &x(0), 1, ap::vlen(0,m-1));
        return;
    }
    
    //
    // check again: correct order of X[]
    //
    if( ap::fp_eq(x(1),x(0)) )
    {
        state.repterminationtype = -2;
        return;
    }
    for(i = 1; i <= m-1; i++)
    {
        if( ap::fp_greater(x(1),x(0))&&ap::fp_less_eq(x(i),x(i-1))||ap::fp_less(x(1),x(0))&&ap::fp_greater_eq(x(i),x(i-1)) )
        {
            state.repterminationtype = -2;
            return;
        }
    }
    
    //
    // auto-select H if necessary
    //
    if( ap::fp_eq(h,0) )
    {
        v = fabs(x(1)-x(0));
        for(i = 2; i <= m-1; i++)
        {
            v = ap::minreal(v, fabs(x(i)-x(i-1)));
        }
        h = 0.001*v;
    }
    
    //
    // store parameters
    //
    state.n = n;
    state.m = m;
    state.h = h;
    state.eps = fabs(eps);
    state.fraceps = ap::fp_less(eps,0);
    state.xg.setlength(m);
    ap::vmove(&state.xg(0), 1, &x(0), 1, ap::vlen(0,m-1));
    if( ap::fp_greater(x(1),x(0)) )
    {
        state.xscale = 1;
    }
    else
    {
        state.xscale = -1;
        ap::vmul(&state.xg(0), 1, ap::vlen(0,m-1), -1);
    }
    state.yc.setlength(n);
    ap::vmove(&state.yc(0), 1, &y(0), 1, ap::vlen(0,n-1));
    state.solvertype = solvertype;
    state.repterminationtype = 0;
}




