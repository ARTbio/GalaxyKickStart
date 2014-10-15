/*************************************************************************
ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
JORGE J. MORE', DAVID J. THUENTE

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
#include "linmin.h"

static const double ftol = 0.001;
static const double xtol = 100*ap::machineepsilon;
static const double gtol = 0.3;
static const int maxfev = 20;
static const double stpmin = 1.0E-50;
static const double defstpmax = 1.0E+50;

static void mcstep(double& stx,
     double& fx,
     double& dx,
     double& sty,
     double& fy,
     double& dy,
     double& stp,
     const double& fp,
     const double& dp,
     bool& brackt,
     const double& stmin,
     const double& stmax,
     int& info);

/*************************************************************************
Normalizes direction/step pair: makes |D|=1, scales Stp.
If |D|=0, it returns, leavind D/Stp unchanged.

  -- ALGLIB --
     Copyright 01.04.2010 by Bochkanov Sergey
*************************************************************************/
void linminnormalized(ap::real_1d_array& d, double& stp, int n)
{
    double mx;
    double s;
    int i;

    
    //
    // first, scale D to avoid underflow/overflow durng squaring
    //
    mx = 0;
    for(i = 0; i <= n-1; i++)
    {
        mx = ap::maxreal(mx, fabs(d(i)));
    }
    if( ap::fp_eq(mx,0) )
    {
        return;
    }
    s = 1/mx;
    ap::vmul(&d(0), 1, ap::vlen(0,n-1), s);
    stp = stp/s;
    
    //
    // normalize D
    //
    s = ap::vdotproduct(&d(0), 1, &d(0), 1, ap::vlen(0,n-1));
    s = 1/sqrt(s);
    ap::vmul(&d(0), 1, ap::vlen(0,n-1), s);
    stp = stp/s;
}


/*************************************************************************
THE  PURPOSE  OF  MCSRCH  IS  TO  FIND A STEP WHICH SATISFIES A SUFFICIENT
DECREASE CONDITION AND A CURVATURE CONDITION.

AT EACH STAGE THE SUBROUTINE  UPDATES  AN  INTERVAL  OF  UNCERTAINTY  WITH
ENDPOINTS  STX  AND  STY.  THE INTERVAL OF UNCERTAINTY IS INITIALLY CHOSEN
SO THAT IT CONTAINS A MINIMIZER OF THE MODIFIED FUNCTION

    F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).

IF  A STEP  IS OBTAINED FOR  WHICH THE MODIFIED FUNCTION HAS A NONPOSITIVE
FUNCTION  VALUE  AND  NONNEGATIVE  DERIVATIVE,   THEN   THE   INTERVAL  OF
UNCERTAINTY IS CHOSEN SO THAT IT CONTAINS A MINIMIZER OF F(X+STP*S).

THE  ALGORITHM  IS  DESIGNED TO FIND A STEP WHICH SATISFIES THE SUFFICIENT
DECREASE CONDITION

    F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S),

AND THE CURVATURE CONDITION

    ABS(GRADF(X+STP*S)'S)) .LE. GTOL*ABS(GRADF(X)'S).

IF  FTOL  IS  LESS  THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION IS BOUNDED
BELOW,  THEN  THERE  IS  ALWAYS  A  STEP  WHICH SATISFIES BOTH CONDITIONS.
IF  NO  STEP  CAN BE FOUND  WHICH  SATISFIES  BOTH  CONDITIONS,  THEN  THE
ALGORITHM  USUALLY STOPS  WHEN  ROUNDING ERRORS  PREVENT FURTHER PROGRESS.
IN THIS CASE STP ONLY SATISFIES THE SUFFICIENT DECREASE CONDITION.

PARAMETERS DESCRIPRION

N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF VARIABLES.

X IS  AN  ARRAY  OF  LENGTH N. ON INPUT IT MUST CONTAIN THE BASE POINT FOR
THE LINE SEARCH. ON OUTPUT IT CONTAINS X+STP*S.

F IS  A  VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F AT X. ON OUTPUT
IT CONTAINS THE VALUE OF F AT X + STP*S.

G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE GRADIENT OF F AT X.
ON OUTPUT IT CONTAINS THE GRADIENT OF F AT X + STP*S.

S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE SEARCH DIRECTION.

STP  IS  A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN INITIAL ESTIMATE
OF A SATISFACTORY STEP. ON OUTPUT STP CONTAINS THE FINAL ESTIMATE.

FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. TERMINATION OCCURS WHEN THE
SUFFICIENT DECREASE CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE
SATISFIED.

XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS WHEN THE RELATIVE
WIDTH OF THE INTERVAL OF UNCERTAINTY IS AT MOST XTOL.

STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH SPECIFY LOWER  AND
UPPER BOUNDS FOR THE STEP.

MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION OCCURS WHEN THE
NUMBER OF CALLS TO FCN IS AT LEAST MAXFEV BY THE END OF AN ITERATION.

INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
    INFO = 0  IMPROPER INPUT PARAMETERS.

    INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
              DIRECTIONAL DERIVATIVE CONDITION HOLD.

    INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
              IS AT MOST XTOL.

    INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.

    INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.

    INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.

    INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
              THERE MAY NOT BE A STEP WHICH SATISFIES THE
              SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
              TOLERANCES MAY BE TOO SMALL.

NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF CALLS TO FCN.

WA IS A WORK ARRAY OF LENGTH N.

ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
JORGE J. MORE', DAVID J. THUENTE
*************************************************************************/
void mcsrch(const int& n,
     ap::real_1d_array& x,
     double& f,
     ap::real_1d_array& g,
     const ap::real_1d_array& s,
     double& stp,
     double stpmax,
     int& info,
     int& nfev,
     ap::real_1d_array& wa,
     linminstate& state,
     int& stage)
{
    double v;
    double p5;
    double p66;
    double zero;

    
    //
    // init
    //
    p5 = 0.5;
    p66 = 0.66;
    state.xtrapf = 4.0;
    zero = 0;
    if( ap::fp_eq(stpmax,0) )
    {
        stpmax = defstpmax;
    }
    if( ap::fp_less(stp,stpmin) )
    {
        stp = stpmin;
    }
    if( ap::fp_greater(stp,stpmax) )
    {
        stp = stpmax;
    }
    
    //
    // Main cycle
    //
    while(true)
    {
        if( stage==0 )
        {
            
            //
            // NEXT
            //
            stage = 2;
            continue;
        }
        if( stage==2 )
        {
            state.infoc = 1;
            info = 0;
            
            //
            //     CHECK THE INPUT PARAMETERS FOR ERRORS.
            //
            if( n<=0||ap::fp_less_eq(stp,0)||ap::fp_less(ftol,0)||ap::fp_less(gtol,zero)||ap::fp_less(xtol,zero)||ap::fp_less(stpmin,zero)||ap::fp_less(stpmax,stpmin)||maxfev<=0 )
            {
                stage = 0;
                return;
            }
            
            //
            //     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
            //     AND CHECK THAT S IS A DESCENT DIRECTION.
            //
            v = ap::vdotproduct(&g(0), 1, &s(0), 1, ap::vlen(0,n-1));
            state.dginit = v;
            if( ap::fp_greater_eq(state.dginit,0) )
            {
                stage = 0;
                return;
            }
            
            //
            //     INITIALIZE LOCAL VARIABLES.
            //
            state.brackt = false;
            state.stage1 = true;
            nfev = 0;
            state.finit = f;
            state.dgtest = ftol*state.dginit;
            state.width = stpmax-stpmin;
            state.width1 = state.width/p5;
            ap::vmove(&wa(0), 1, &x(0), 1, ap::vlen(0,n-1));
            
            //
            //     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
            //     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
            //     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
            //     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
            //     THE INTERVAL OF UNCERTAINTY.
            //     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
            //     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
            //
            state.stx = 0;
            state.fx = state.finit;
            state.dgx = state.dginit;
            state.sty = 0;
            state.fy = state.finit;
            state.dgy = state.dginit;
            
            //
            // NEXT
            //
            stage = 3;
            continue;
        }
        if( stage==3 )
        {
            
            //
            //     START OF ITERATION.
            //
            //     SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
            //     TO THE PRESENT INTERVAL OF UNCERTAINTY.
            //
            if( state.brackt )
            {
                if( ap::fp_less(state.stx,state.sty) )
                {
                    state.stmin = state.stx;
                    state.stmax = state.sty;
                }
                else
                {
                    state.stmin = state.sty;
                    state.stmax = state.stx;
                }
            }
            else
            {
                state.stmin = state.stx;
                state.stmax = stp+state.xtrapf*(stp-state.stx);
            }
            
            //
            //        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
            //
            if( ap::fp_greater(stp,stpmax) )
            {
                stp = stpmax;
            }
            if( ap::fp_less(stp,stpmin) )
            {
                stp = stpmin;
            }
            
            //
            //        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
            //        STP BE THE LOWEST POINT OBTAINED SO FAR.
            //
            if( state.brackt&&(ap::fp_less_eq(stp,state.stmin)||ap::fp_greater_eq(stp,state.stmax))||nfev>=maxfev-1||state.infoc==0||state.brackt&&ap::fp_less_eq(state.stmax-state.stmin,xtol*state.stmax) )
            {
                stp = state.stx;
            }
            
            //
            //        EVALUATE THE FUNCTION AND GRADIENT AT STP
            //        AND COMPUTE THE DIRECTIONAL DERIVATIVE.
            //
            ap::vmove(&x(0), 1, &wa(0), 1, ap::vlen(0,n-1));
            ap::vadd(&x(0), 1, &s(0), 1, ap::vlen(0,n-1), stp);
            
            //
            // NEXT
            //
            stage = 4;
            return;
        }
        if( stage==4 )
        {
            info = 0;
            nfev = nfev+1;
            v = ap::vdotproduct(&g(0), 1, &s(0), 1, ap::vlen(0,n-1));
            state.dg = v;
            state.ftest1 = state.finit+stp*state.dgtest;
            
            //
            //        TEST FOR CONVERGENCE.
            //
            if( state.brackt&&(ap::fp_less_eq(stp,state.stmin)||ap::fp_greater_eq(stp,state.stmax))||state.infoc==0 )
            {
                info = 6;
            }
            if( ap::fp_eq(stp,stpmax)&&ap::fp_less_eq(f,state.ftest1)&&ap::fp_less_eq(state.dg,state.dgtest) )
            {
                info = 5;
            }
            if( ap::fp_eq(stp,stpmin)&&(ap::fp_greater(f,state.ftest1)||ap::fp_greater_eq(state.dg,state.dgtest)) )
            {
                info = 4;
            }
            if( nfev>=maxfev )
            {
                info = 3;
            }
            if( state.brackt&&ap::fp_less_eq(state.stmax-state.stmin,xtol*state.stmax) )
            {
                info = 2;
            }
            if( ap::fp_less_eq(f,state.ftest1)&&ap::fp_less_eq(fabs(state.dg),-gtol*state.dginit) )
            {
                info = 1;
            }
            
            //
            //        CHECK FOR TERMINATION.
            //
            if( info!=0 )
            {
                stage = 0;
                return;
            }
            
            //
            //        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
            //        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
            //
            if( state.stage1&&ap::fp_less_eq(f,state.ftest1)&&ap::fp_greater_eq(state.dg,ap::minreal(ftol, gtol)*state.dginit) )
            {
                state.stage1 = false;
            }
            
            //
            //        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
            //        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
            //        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
            //        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
            //        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
            //
            if( state.stage1&&ap::fp_less_eq(f,state.fx)&&ap::fp_greater(f,state.ftest1) )
            {
                
                //
                //           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
                //
                state.fm = f-stp*state.dgtest;
                state.fxm = state.fx-state.stx*state.dgtest;
                state.fym = state.fy-state.sty*state.dgtest;
                state.dgm = state.dg-state.dgtest;
                state.dgxm = state.dgx-state.dgtest;
                state.dgym = state.dgy-state.dgtest;
                
                //
                //           CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
                //           AND TO COMPUTE THE NEW STEP.
                //
                mcstep(state.stx, state.fxm, state.dgxm, state.sty, state.fym, state.dgym, stp, state.fm, state.dgm, state.brackt, state.stmin, state.stmax, state.infoc);
                
                //
                //           RESET THE FUNCTION AND GRADIENT VALUES FOR F.
                //
                state.fx = state.fxm+state.stx*state.dgtest;
                state.fy = state.fym+state.sty*state.dgtest;
                state.dgx = state.dgxm+state.dgtest;
                state.dgy = state.dgym+state.dgtest;
            }
            else
            {
                
                //
                //           CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
                //           AND TO COMPUTE THE NEW STEP.
                //
                mcstep(state.stx, state.fx, state.dgx, state.sty, state.fy, state.dgy, stp, f, state.dg, state.brackt, state.stmin, state.stmax, state.infoc);
            }
            
            //
            //        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
            //        INTERVAL OF UNCERTAINTY.
            //
            if( state.brackt )
            {
                if( ap::fp_greater_eq(fabs(state.sty-state.stx),p66*state.width1) )
                {
                    stp = state.stx+p5*(state.sty-state.stx);
                }
                state.width1 = state.width;
                state.width = fabs(state.sty-state.stx);
            }
            
            //
            //  NEXT.
            //
            stage = 3;
            continue;
        }
    }
}


static void mcstep(double& stx,
     double& fx,
     double& dx,
     double& sty,
     double& fy,
     double& dy,
     double& stp,
     const double& fp,
     const double& dp,
     bool& brackt,
     const double& stmin,
     const double& stmax,
     int& info)
{
    bool bound;
    double gamma;
    double p;
    double q;
    double r;
    double s;
    double sgnd;
    double stpc;
    double stpf;
    double stpq;
    double theta;

    info = 0;
    
    //
    //     CHECK THE INPUT PARAMETERS FOR ERRORS.
    //
    if( brackt&&(ap::fp_less_eq(stp,ap::minreal(stx, sty))||ap::fp_greater_eq(stp,ap::maxreal(stx, sty)))||ap::fp_greater_eq(dx*(stp-stx),0)||ap::fp_less(stmax,stmin) )
    {
        return;
    }
    
    //
    //     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.
    //
    sgnd = dp*(dx/fabs(dx));
    
    //
    //     FIRST CASE. A HIGHER FUNCTION VALUE.
    //     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
    //     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
    //     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.
    //
    if( ap::fp_greater(fp,fx) )
    {
        info = 1;
        bound = true;
        theta = 3*(fx-fp)/(stp-stx)+dx+dp;
        s = ap::maxreal(fabs(theta), ap::maxreal(fabs(dx), fabs(dp)));
        gamma = s*sqrt(ap::sqr(theta/s)-dx/s*(dp/s));
        if( ap::fp_less(stp,stx) )
        {
            gamma = -gamma;
        }
        p = gamma-dx+theta;
        q = gamma-dx+gamma+dp;
        r = p/q;
        stpc = stx+r*(stp-stx);
        stpq = stx+dx/((fx-fp)/(stp-stx)+dx)/2*(stp-stx);
        if( ap::fp_less(fabs(stpc-stx),fabs(stpq-stx)) )
        {
            stpf = stpc;
        }
        else
        {
            stpf = stpc+(stpq-stpc)/2;
        }
        brackt = true;
    }
    else
    {
        if( ap::fp_less(sgnd,0) )
        {
            
            //
            //     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
            //     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
            //     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
            //     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
            //
            info = 2;
            bound = false;
            theta = 3*(fx-fp)/(stp-stx)+dx+dp;
            s = ap::maxreal(fabs(theta), ap::maxreal(fabs(dx), fabs(dp)));
            gamma = s*sqrt(ap::sqr(theta/s)-dx/s*(dp/s));
            if( ap::fp_greater(stp,stx) )
            {
                gamma = -gamma;
            }
            p = gamma-dp+theta;
            q = gamma-dp+gamma+dx;
            r = p/q;
            stpc = stp+r*(stx-stp);
            stpq = stp+dp/(dp-dx)*(stx-stp);
            if( ap::fp_greater(fabs(stpc-stp),fabs(stpq-stp)) )
            {
                stpf = stpc;
            }
            else
            {
                stpf = stpq;
            }
            brackt = true;
        }
        else
        {
            if( ap::fp_less(fabs(dp),fabs(dx)) )
            {
                
                //
                //     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
                //     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
                //     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
                //     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
                //     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
                //     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
                //     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
                //     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.
                //
                info = 3;
                bound = true;
                theta = 3*(fx-fp)/(stp-stx)+dx+dp;
                s = ap::maxreal(fabs(theta), ap::maxreal(fabs(dx), fabs(dp)));
                
                //
                //        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
                //        TO INFINITY IN THE DIRECTION OF THE STEP.
                //
                gamma = s*sqrt(ap::maxreal(double(0), ap::sqr(theta/s)-dx/s*(dp/s)));
                if( ap::fp_greater(stp,stx) )
                {
                    gamma = -gamma;
                }
                p = gamma-dp+theta;
                q = gamma+(dx-dp)+gamma;
                r = p/q;
                if( ap::fp_less(r,0)&&ap::fp_neq(gamma,0) )
                {
                    stpc = stp+r*(stx-stp);
                }
                else
                {
                    if( ap::fp_greater(stp,stx) )
                    {
                        stpc = stmax;
                    }
                    else
                    {
                        stpc = stmin;
                    }
                }
                stpq = stp+dp/(dp-dx)*(stx-stp);
                if( brackt )
                {
                    if( ap::fp_less(fabs(stp-stpc),fabs(stp-stpq)) )
                    {
                        stpf = stpc;
                    }
                    else
                    {
                        stpf = stpq;
                    }
                }
                else
                {
                    if( ap::fp_greater(fabs(stp-stpc),fabs(stp-stpq)) )
                    {
                        stpf = stpc;
                    }
                    else
                    {
                        stpf = stpq;
                    }
                }
            }
            else
            {
                
                //
                //     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
                //     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
                //     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
                //     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
                //
                info = 4;
                bound = false;
                if( brackt )
                {
                    theta = 3*(fp-fy)/(sty-stp)+dy+dp;
                    s = ap::maxreal(fabs(theta), ap::maxreal(fabs(dy), fabs(dp)));
                    gamma = s*sqrt(ap::sqr(theta/s)-dy/s*(dp/s));
                    if( ap::fp_greater(stp,sty) )
                    {
                        gamma = -gamma;
                    }
                    p = gamma-dp+theta;
                    q = gamma-dp+gamma+dy;
                    r = p/q;
                    stpc = stp+r*(sty-stp);
                    stpf = stpc;
                }
                else
                {
                    if( ap::fp_greater(stp,stx) )
                    {
                        stpf = stmax;
                    }
                    else
                    {
                        stpf = stmin;
                    }
                }
            }
        }
    }
    
    //
    //     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
    //     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.
    //
    if( ap::fp_greater(fp,fx) )
    {
        sty = stp;
        fy = fp;
        dy = dp;
    }
    else
    {
        if( ap::fp_less(sgnd,0.0) )
        {
            sty = stx;
            fy = fx;
            dy = dx;
        }
        stx = stp;
        fx = fp;
        dx = dp;
    }
    
    //
    //     COMPUTE THE NEW STEP AND SAFEGUARD IT.
    //
    stpf = ap::minreal(stmax, stpf);
    stpf = ap::maxreal(stmin, stpf);
    stp = stpf;
    if( brackt&&bound )
    {
        if( ap::fp_greater(sty,stx) )
        {
            stp = ap::minreal(stx+0.66*(sty-stx), stp);
        }
        else
        {
            stp = ap::maxreal(stx+0.66*(sty-stx), stp);
        }
    }
}




