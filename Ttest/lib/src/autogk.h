/*************************************************************************
Copyright (c) 2005-2009, Sergey Bochkanov (ALGLIB project).

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

#ifndef _autogk_h
#define _autogk_h

#include "ap.h"
#include "ialglib.h"

#include "tsort.h"
#include "hblas.h"
#include "reflections.h"
#include "creflections.h"
#include "sblas.h"
#include "ablasf.h"
#include "ablas.h"
#include "ortfac.h"
#include "blas.h"
#include "rotations.h"
#include "hsschur.h"
#include "evd.h"
#include "gammafunc.h"
#include "gq.h"
#include "gkq.h"


/*************************************************************************
Integration report:
* TerminationType = completetion code:
    * -5    non-convergence of Gauss-Kronrod nodes
            calculation subroutine.
    * -1    incorrect parameters were specified
    *  1    OK
* Rep.NFEV countains number of function calculations
* Rep.NIntervals contains number of intervals [a,b]
  was partitioned into.
*************************************************************************/
struct autogkreport
{
    int terminationtype;
    int nfev;
    int nintervals;
};


struct autogkinternalstate
{
    double a;
    double b;
    double eps;
    double xwidth;
    double x;
    double f;
    int info;
    double r;
    ap::real_2d_array heap;
    int heapsize;
    int heapwidth;
    int heapused;
    double sumerr;
    double sumabs;
    ap::real_1d_array qn;
    ap::real_1d_array wg;
    ap::real_1d_array wk;
    ap::real_1d_array wr;
    int n;
    ap::rcommstate rstate;
};


/*************************************************************************
This structure stores internal state of the integration algorithm  between
subsequent calls of the AutoGKIteration() subroutine.
*************************************************************************/
struct autogkstate
{
    double a;
    double b;
    double alpha;
    double beta;
    double xwidth;
    double x;
    double xminusa;
    double bminusx;
    double f;
    int wrappermode;
    autogkinternalstate internalstate;
    ap::rcommstate rstate;
    double v;
    int terminationtype;
    int nfev;
    int nintervals;
};




/*************************************************************************
Integration of a smooth function F(x) on a finite interval [a,b].

Fast-convergent algorithm based on a Gauss-Kronrod formula is used. Result
is calculated with accuracy close to the machine precision.

Algorithm works well only with smooth integrands.  It  may  be  used  with
continuous non-smooth integrands, but with  less  performance.

It should never be used with integrands which have integrable singularities
at lower or upper limits - algorithm may crash. Use AutoGKSingular in such
cases.

INPUT PARAMETERS:
    A, B    -   interval boundaries (A<B, A=B or A>B)
    
OUTPUT PARAMETERS
    State   -   structure which stores algorithm state between  subsequent
                calls of AutoGKIteration.  Used for reverse communication.
                This structure should be  passed  to  the  AutoGKIteration
                subroutine.

SEE ALSO
    AutoGKSmoothW, AutoGKSingular, AutoGKIteration, AutoGKResults.
    

  -- ALGLIB --
     Copyright 06.05.2009 by Bochkanov Sergey
*************************************************************************/
void autogksmooth(double a, double b, autogkstate& state);


/*************************************************************************
Integration of a smooth function F(x) on a finite interval [a,b].

This subroutine is same as AutoGKSmooth(), but it guarantees that interval
[a,b] is partitioned into subintervals which have width at most XWidth.

Subroutine  can  be  used  when  integrating nearly-constant function with
narrow "bumps" (about XWidth wide). If "bumps" are too narrow, AutoGKSmooth
subroutine can overlook them.

INPUT PARAMETERS:
    A, B    -   interval boundaries (A<B, A=B or A>B)

OUTPUT PARAMETERS
    State   -   structure which stores algorithm state between  subsequent
                calls of AutoGKIteration.  Used for reverse communication.
                This structure should be  passed  to  the  AutoGKIteration
                subroutine.

SEE ALSO
    AutoGKSmooth, AutoGKSingular, AutoGKIteration, AutoGKResults.


  -- ALGLIB --
     Copyright 06.05.2009 by Bochkanov Sergey
*************************************************************************/
void autogksmoothw(double a, double b, double xwidth, autogkstate& state);


/*************************************************************************
Integration on a finite interval [A,B].
Integrand have integrable singularities at A/B.

F(X) must diverge as "(x-A)^alpha" at A, as "(B-x)^beta" at B,  with known
alpha/beta (alpha>-1, beta>-1).  If alpha/beta  are  not known,  estimates
from below can be used (but these estimates should be greater than -1 too).

One  of  alpha/beta variables (or even both alpha/beta) may be equal to 0,
which means than function F(x) is non-singular at A/B. Anyway (singular at
bounds or not), function F(x) is supposed to be continuous on (A,B).

Fast-convergent algorithm based on a Gauss-Kronrod formula is used. Result
is calculated with accuracy close to the machine precision.

INPUT PARAMETERS:
    A, B    -   interval boundaries (A<B, A=B or A>B)
    Alpha   -   power-law coefficient of the F(x) at A,
                Alpha>-1
    Beta    -   power-law coefficient of the F(x) at B,
                Beta>-1

OUTPUT PARAMETERS
    State   -   structure which stores algorithm state between  subsequent
                calls of AutoGKIteration.  Used for reverse communication.
                This structure should be  passed  to  the  AutoGKIteration
                subroutine.

SEE ALSO
    AutoGKSmooth, AutoGKSmoothW, AutoGKIteration, AutoGKResults.


  -- ALGLIB --
     Copyright 06.05.2009 by Bochkanov Sergey
*************************************************************************/
void autogksingular(double a,
     double b,
     double alpha,
     double beta,
     autogkstate& state);


/*************************************************************************
One step of adaptive integration process.

Called after initialization with one of AutoGKXXX subroutines.
See HTML documentation for examples.

Input parameters:
    State   -   structure which stores algorithm state between  calls  and
                which  is  used  for   reverse   communication.   Must  be
                initialized with one of AutoGKXXX subroutines.

If suborutine returned False, iterative proces has converged. If subroutine
returned True, caller should calculate function value State.F  at  State.X
and call AutoGKIteration again.

NOTE:

When integrating "difficult" functions with integrable singularities like

    F(x) = (x-A)^alpha * (B-x)^beta

subroutine may require the value of F at points which are too close to A/B.
Sometimes to calculate integral with high enough precision we  may need to
calculate F(A+delta) when delta is less than machine  epsilon.  In  finite
precision arithmetics A+delta will be effectively equal to A,  so  we  may
find us in situation when  we  are  trying  to  calculate  something  like
1/sqrt(1-1).

To avoid  such  situations,  AutoGKIteration  subroutine  fills  not  only
State.X  field,  but  also   State.XMinusA   (which  equals  to  X-A)  and
State.BMinusX  (which  equals to B-X) fields.  If X is too close to A or B
(X-A<0.001*A, or B-X<0.001*B, for example) use  these  fields  instead  of
State.X


  -- ALGLIB --
     Copyright 07.05.2009 by Bochkanov Sergey
*************************************************************************/
bool autogkiteration(autogkstate& state);


/*************************************************************************
Adaptive integration results

Called after AutoGKIteration returned False.

Input parameters:
    State   -   algorithm state (used by AutoGKIteration).

Output parameters:
    V       -   integral(f(x)dx,a,b)
    Rep     -   optimization report (see AutoGKReport description)

  -- ALGLIB --
     Copyright 14.11.2007 by Bochkanov Sergey
*************************************************************************/
void autogkresults(const autogkstate& state, double& v, autogkreport& rep);


#endif

