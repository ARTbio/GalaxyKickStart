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

#ifndef _odesolver_h
#define _odesolver_h

#include "ap.h"
#include "ialglib.h"

struct odesolverstate
{
    int n;
    int m;
    double xscale;
    double h;
    double eps;
    bool fraceps;
    ap::real_1d_array yc;
    ap::real_1d_array escale;
    ap::real_1d_array xg;
    int solvertype;
    double x;
    ap::real_1d_array y;
    ap::real_1d_array dy;
    ap::real_2d_array ytbl;
    int repterminationtype;
    int repnfev;
    ap::real_1d_array yn;
    ap::real_1d_array yns;
    ap::real_1d_array rka;
    ap::real_1d_array rkc;
    ap::real_1d_array rkcs;
    ap::real_2d_array rkb;
    ap::real_2d_array rkk;
    ap::rcommstate rstate;
};


struct odesolverreport
{
    int nfev;
    int terminationtype;
};




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
     odesolverstate& state);


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
bool odesolveriteration(odesolverstate& state);


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
     odesolverreport& rep);


#endif

