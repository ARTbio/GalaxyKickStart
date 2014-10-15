/*************************************************************************
Copyright (c) 2008, Sergey Bochkanov (ALGLIB project).

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
#include "logit.h"

static const double xtol = 100*ap::machineepsilon;
static const double ftol = 0.0001;
static const double gtol = 0.3;
static const int maxfev = 20;
static const double stpmin = 1.0E-2;
static const double stpmax = 1.0E5;
static const int logitvnum = 6;

static void mnliexp(ap::real_1d_array& w, const ap::real_1d_array& x);
static void mnlallerrors(logitmodel& lm,
     const ap::real_2d_array& xy,
     int npoints,
     double& relcls,
     double& avgce,
     double& rms,
     double& avg,
     double& avgrel);
static void mnlmcsrch(const int& n,
     ap::real_1d_array& x,
     double& f,
     ap::real_1d_array& g,
     const ap::real_1d_array& s,
     double& stp,
     int& info,
     int& nfev,
     ap::real_1d_array& wa,
     logitmcstate& state,
     int& stage);
static void mnlmcstep(double& stx,
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
This subroutine trains logit model.

INPUT PARAMETERS:
    XY          -   training set, array[0..NPoints-1,0..NVars]
                    First NVars columns store values of independent
                    variables, next column stores number of class (from 0
                    to NClasses-1) which dataset element belongs to. Fractional
                    values are rounded to nearest integer.
    NPoints     -   training set size, NPoints>=1
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   number of classes, NClasses>=2

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<NVars+2, NVars<1, NClasses<2).
                    *  1, if task has been solved
    LM          -   model built
    Rep         -   training report

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
void mnltrainh(const ap::real_2d_array& xy,
     int npoints,
     int nvars,
     int nclasses,
     int& info,
     logitmodel& lm,
     mnlreport& rep)
{
    int i;
    int j;
    int k;
    int ssize;
    bool allsame;
    int offs;
    double threshold;
    double wminstep;
    double decay;
    int wdim;
    int expoffs;
    double v;
    double s;
    multilayerperceptron network;
    int nin;
    int nout;
    int wcount;
    double e;
    ap::real_1d_array g;
    ap::real_2d_array h;
    bool spd;
    ap::real_1d_array x;
    ap::real_1d_array y;
    ap::real_1d_array wbase;
    double wstep;
    ap::real_1d_array wdir;
    ap::real_1d_array work;
    int mcstage;
    logitmcstate mcstate;
    int mcinfo;
    int mcnfev;
    int solverinfo;
    densesolverreport solverrep;

    threshold = 1000*ap::machineepsilon;
    wminstep = 0.001;
    decay = 0.001;
    
    //
    // Test for inputs
    //
    if( npoints<nvars+2||nvars<1||nclasses<2 )
    {
        info = -1;
        return;
    }
    for(i = 0; i <= npoints-1; i++)
    {
        if( ap::round(xy(i,nvars))<0||ap::round(xy(i,nvars))>=nclasses )
        {
            info = -2;
            return;
        }
    }
    info = 1;
    
    //
    // Initialize data
    //
    rep.ngrad = 0;
    rep.nhess = 0;
    
    //
    // Allocate array
    //
    wdim = (nvars+1)*(nclasses-1);
    offs = 5;
    expoffs = offs+wdim;
    ssize = 5+(nvars+1)*(nclasses-1)+nclasses;
    lm.w.setbounds(0, ssize-1);
    lm.w(0) = ssize;
    lm.w(1) = logitvnum;
    lm.w(2) = nvars;
    lm.w(3) = nclasses;
    lm.w(4) = offs;
    
    //
    // Degenerate case: all outputs are equal
    //
    allsame = true;
    for(i = 1; i <= npoints-1; i++)
    {
        if( ap::round(xy(i,nvars))!=ap::round(xy(i-1,nvars)) )
        {
            allsame = false;
        }
    }
    if( allsame )
    {
        for(i = 0; i <= (nvars+1)*(nclasses-1)-1; i++)
        {
            lm.w(offs+i) = 0;
        }
        v = -2*log(ap::minrealnumber);
        k = ap::round(xy(0,nvars));
        if( k==nclasses-1 )
        {
            for(i = 0; i <= nclasses-2; i++)
            {
                lm.w(offs+i*(nvars+1)+nvars) = -v;
            }
        }
        else
        {
            for(i = 0; i <= nclasses-2; i++)
            {
                if( i==k )
                {
                    lm.w(offs+i*(nvars+1)+nvars) = +v;
                }
                else
                {
                    lm.w(offs+i*(nvars+1)+nvars) = 0;
                }
            }
        }
        return;
    }
    
    //
    // General case.
    // Prepare task and network. Allocate space.
    //
    mlpcreatec0(nvars, nclasses, network);
    mlpinitpreprocessor(network, xy, npoints);
    mlpproperties(network, nin, nout, wcount);
    for(i = 0; i <= wcount-1; i++)
    {
        network.weights(i) = (2*ap::randomreal()-1)/nvars;
    }
    g.setbounds(0, wcount-1);
    h.setbounds(0, wcount-1, 0, wcount-1);
    wbase.setbounds(0, wcount-1);
    wdir.setbounds(0, wcount-1);
    work.setbounds(0, wcount-1);
    
    //
    // First stage: optimize in gradient direction.
    //
    for(k = 0; k <= wcount/3+10; k++)
    {
        
        //
        // Calculate gradient in starting point
        //
        mlpgradnbatch(network, xy, npoints, e, g);
        v = ap::vdotproduct(&network.weights(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1));
        e = e+0.5*decay*v;
        ap::vadd(&g(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1), decay);
        rep.ngrad = rep.ngrad+1;
        
        //
        // Setup optimization scheme
        //
        ap::vmoveneg(&wdir(0), 1, &g(0), 1, ap::vlen(0,wcount-1));
        v = ap::vdotproduct(&wdir(0), 1, &wdir(0), 1, ap::vlen(0,wcount-1));
        wstep = sqrt(v);
        v = 1/sqrt(v);
        ap::vmul(&wdir(0), 1, ap::vlen(0,wcount-1), v);
        mcstage = 0;
        mnlmcsrch(wcount, network.weights, e, g, wdir, wstep, mcinfo, mcnfev, work, mcstate, mcstage);
        while(mcstage!=0)
        {
            mlpgradnbatch(network, xy, npoints, e, g);
            v = ap::vdotproduct(&network.weights(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1));
            e = e+0.5*decay*v;
            ap::vadd(&g(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1), decay);
            rep.ngrad = rep.ngrad+1;
            mnlmcsrch(wcount, network.weights, e, g, wdir, wstep, mcinfo, mcnfev, work, mcstate, mcstage);
        }
    }
    
    //
    // Second stage: use Hessian when we are close to the minimum
    //
    while(true)
    {
        
        //
        // Calculate and update E/G/H
        //
        mlphessiannbatch(network, xy, npoints, e, g, h);
        v = ap::vdotproduct(&network.weights(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1));
        e = e+0.5*decay*v;
        ap::vadd(&g(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1), decay);
        for(k = 0; k <= wcount-1; k++)
        {
            h(k,k) = h(k,k)+decay;
        }
        rep.nhess = rep.nhess+1;
        
        //
        // Select step direction
        // NOTE: it is important to use lower-triangle Cholesky
        // factorization since it is much faster than higher-triangle version.
        //
        spd = spdmatrixcholesky(h, wcount, false);
        spdmatrixcholeskysolve(h, wcount, false, g, solverinfo, solverrep, wdir);
        spd = solverinfo>0;
        if( spd )
        {
            
            //
            // H is positive definite.
            // Step in Newton direction.
            //
            ap::vmul(&wdir(0), 1, ap::vlen(0,wcount-1), -1);
            spd = true;
        }
        else
        {
            
            //
            // H is indefinite.
            // Step in gradient direction.
            //
            ap::vmoveneg(&wdir(0), 1, &g(0), 1, ap::vlen(0,wcount-1));
            spd = false;
        }
        
        //
        // Optimize in WDir direction
        //
        v = ap::vdotproduct(&wdir(0), 1, &wdir(0), 1, ap::vlen(0,wcount-1));
        wstep = sqrt(v);
        v = 1/sqrt(v);
        ap::vmul(&wdir(0), 1, ap::vlen(0,wcount-1), v);
        mcstage = 0;
        mnlmcsrch(wcount, network.weights, e, g, wdir, wstep, mcinfo, mcnfev, work, mcstate, mcstage);
        while(mcstage!=0)
        {
            mlpgradnbatch(network, xy, npoints, e, g);
            v = ap::vdotproduct(&network.weights(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1));
            e = e+0.5*decay*v;
            ap::vadd(&g(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1), decay);
            rep.ngrad = rep.ngrad+1;
            mnlmcsrch(wcount, network.weights, e, g, wdir, wstep, mcinfo, mcnfev, work, mcstate, mcstage);
        }
        if( spd&&(mcinfo==2||mcinfo==4||mcinfo==6) )
        {
            break;
        }
    }
    
    //
    // Convert from NN format to MNL format
    //
    ap::vmove(&lm.w(offs), 1, &network.weights(0), 1, ap::vlen(offs,offs+wcount-1));
    for(k = 0; k <= nvars-1; k++)
    {
        for(i = 0; i <= nclasses-2; i++)
        {
            s = network.columnsigmas(k);
            if( ap::fp_eq(s,0) )
            {
                s = 1;
            }
            j = offs+(nvars+1)*i;
            v = lm.w(j+k);
            lm.w(j+k) = v/s;
            lm.w(j+nvars) = lm.w(j+nvars)+v*network.columnmeans(k)/s;
        }
    }
    for(k = 0; k <= nclasses-2; k++)
    {
        lm.w(offs+(nvars+1)*k+nvars) = -lm.w(offs+(nvars+1)*k+nvars);
    }
}


/*************************************************************************
Procesing

INPUT PARAMETERS:
    LM      -   logit model, passed by non-constant reference
                (some fields of structure are used as temporaries
                when calculating model output).
    X       -   input vector,  array[0..NVars-1].

OUTPUT PARAMETERS:
    Y       -   result, array[0..NClasses-1]
                Vector of posterior probabilities for classification task.
                Subroutine does not allocate memory for this vector, it is
                responsibility of a caller to allocate it. Array  must  be
                at least [0..NClasses-1].

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
void mnlprocess(logitmodel& lm,
     const ap::real_1d_array& x,
     ap::real_1d_array& y)
{
    int nvars;
    int nclasses;
    int offs;
    int i;
    int i1;
    double s;

    ap::ap_error::make_assertion(ap::fp_eq(lm.w(1),logitvnum), "MNLProcess: unexpected model version");
    nvars = ap::round(lm.w(2));
    nclasses = ap::round(lm.w(3));
    offs = ap::round(lm.w(4));
    mnliexp(lm.w, x);
    s = 0;
    i1 = offs+(nvars+1)*(nclasses-1);
    for(i = i1; i <= i1+nclasses-1; i++)
    {
        s = s+lm.w(i);
    }
    for(i = 0; i <= nclasses-1; i++)
    {
        y(i) = lm.w(i1+i)/s;
    }
}


/*************************************************************************
Unpacks coefficients of logit model. Logit model have form:

    P(class=i) = S(i) / (S(0) + S(1) + ... +S(M-1))
          S(i) = Exp(A[i,0]*X[0] + ... + A[i,N-1]*X[N-1] + A[i,N]), when i<M-1
        S(M-1) = 1

INPUT PARAMETERS:
    LM          -   logit model in ALGLIB format

OUTPUT PARAMETERS:
    V           -   coefficients, array[0..NClasses-2,0..NVars]
    NVars       -   number of independent variables
    NClasses    -   number of classes

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
void mnlunpack(const logitmodel& lm,
     ap::real_2d_array& a,
     int& nvars,
     int& nclasses)
{
    int offs;
    int i;

    ap::ap_error::make_assertion(ap::fp_eq(lm.w(1),logitvnum), "MNLUnpack: unexpected model version");
    nvars = ap::round(lm.w(2));
    nclasses = ap::round(lm.w(3));
    offs = ap::round(lm.w(4));
    a.setbounds(0, nclasses-2, 0, nvars);
    for(i = 0; i <= nclasses-2; i++)
    {
        ap::vmove(&a(i, 0), 1, &lm.w(offs+i*(nvars+1)), 1, ap::vlen(0,nvars));
    }
}


/*************************************************************************
"Packs" coefficients and creates logit model in ALGLIB format (MNLUnpack
reversed).

INPUT PARAMETERS:
    A           -   model (see MNLUnpack)
    NVars       -   number of independent variables
    NClasses    -   number of classes

OUTPUT PARAMETERS:
    LM          -   logit model.

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
void mnlpack(const ap::real_2d_array& a,
     int nvars,
     int nclasses,
     logitmodel& lm)
{
    int offs;
    int i;
    int wdim;
    int ssize;

    wdim = (nvars+1)*(nclasses-1);
    offs = 5;
    ssize = 5+(nvars+1)*(nclasses-1)+nclasses;
    lm.w.setbounds(0, ssize-1);
    lm.w(0) = ssize;
    lm.w(1) = logitvnum;
    lm.w(2) = nvars;
    lm.w(3) = nclasses;
    lm.w(4) = offs;
    for(i = 0; i <= nclasses-2; i++)
    {
        ap::vmove(&lm.w(offs+i*(nvars+1)), 1, &a(i, 0), 1, ap::vlen(offs+i*(nvars+1),offs+i*(nvars+1)+nvars));
    }
}


/*************************************************************************
Copying of LogitModel strucure

INPUT PARAMETERS:
    LM1 -   original

OUTPUT PARAMETERS:
    LM2 -   copy

  -- ALGLIB --
     Copyright 15.03.2009 by Bochkanov Sergey
*************************************************************************/
void mnlcopy(const logitmodel& lm1, logitmodel& lm2)
{
    int k;

    k = ap::round(lm1.w(0));
    lm2.w.setbounds(0, k-1);
    ap::vmove(&lm2.w(0), 1, &lm1.w(0), 1, ap::vlen(0,k-1));
}


/*************************************************************************
Serialization of LogitModel strucure

INPUT PARAMETERS:
    LM      -   original

OUTPUT PARAMETERS:
    RA      -   array of real numbers which stores model,
                array[0..RLen-1]
    RLen    -   RA lenght

  -- ALGLIB --
     Copyright 15.03.2009 by Bochkanov Sergey
*************************************************************************/
void mnlserialize(const logitmodel& lm, ap::real_1d_array& ra, int& rlen)
{

    rlen = ap::round(lm.w(0))+1;
    ra.setbounds(0, rlen-1);
    ra(0) = logitvnum;
    ap::vmove(&ra(1), 1, &lm.w(0), 1, ap::vlen(1,rlen-1));
}


/*************************************************************************
Unserialization of LogitModel strucure

INPUT PARAMETERS:
    RA      -   real array which stores model

OUTPUT PARAMETERS:
    LM      -   restored model

  -- ALGLIB --
     Copyright 15.03.2009 by Bochkanov Sergey
*************************************************************************/
void mnlunserialize(const ap::real_1d_array& ra, logitmodel& lm)
{

    ap::ap_error::make_assertion(ap::round(ra(0))==logitvnum, "MNLUnserialize: incorrect array!");
    lm.w.setbounds(0, ap::round(ra(1))-1);
    ap::vmove(&lm.w(0), 1, &ra(1), 1, ap::vlen(0,ap::round(ra(1))-1));
}


/*************************************************************************
Average cross-entropy (in bits per element) on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    CrossEntropy/(NPoints*ln(2)).

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
double mnlavgce(logitmodel& lm, const ap::real_2d_array& xy, int npoints)
{
    double result;
    int nvars;
    int nclasses;
    int i;
    ap::real_1d_array workx;
    ap::real_1d_array worky;

    ap::ap_error::make_assertion(ap::fp_eq(lm.w(1),logitvnum), "MNLClsError: unexpected model version");
    nvars = ap::round(lm.w(2));
    nclasses = ap::round(lm.w(3));
    workx.setbounds(0, nvars-1);
    worky.setbounds(0, nclasses-1);
    result = 0;
    for(i = 0; i <= npoints-1; i++)
    {
        ap::ap_error::make_assertion(ap::round(xy(i,nvars))>=0&&ap::round(xy(i,nvars))<nclasses, "MNLAvgCE: incorrect class number!");
        
        //
        // Process
        //
        ap::vmove(&workx(0), 1, &xy(i, 0), 1, ap::vlen(0,nvars-1));
        mnlprocess(lm, workx, worky);
        if( ap::fp_greater(worky(ap::round(xy(i,nvars))),0) )
        {
            result = result-log(worky(ap::round(xy(i,nvars))));
        }
        else
        {
            result = result-log(ap::minrealnumber);
        }
    }
    result = result/(npoints*log(double(2)));
    return result;
}


/*************************************************************************
Relative classification error on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    percent of incorrectly classified cases.

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
double mnlrelclserror(logitmodel& lm,
     const ap::real_2d_array& xy,
     int npoints)
{
    double result;

    result = double(mnlclserror(lm, xy, npoints))/double(npoints);
    return result;
}


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error (error when estimating posterior probabilities).

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double mnlrmserror(logitmodel& lm, const ap::real_2d_array& xy, int npoints)
{
    double result;
    double relcls;
    double avgce;
    double rms;
    double avg;
    double avgrel;

    ap::ap_error::make_assertion(ap::round(lm.w(1))==logitvnum, "MNLRMSError: Incorrect MNL version!");
    mnlallerrors(lm, xy, npoints, relcls, avgce, rms, avg, avgrel);
    result = rms;
    return result;
}


/*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average error (error when estimating posterior probabilities).

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double mnlavgerror(logitmodel& lm, const ap::real_2d_array& xy, int npoints)
{
    double result;
    double relcls;
    double avgce;
    double rms;
    double avg;
    double avgrel;

    ap::ap_error::make_assertion(ap::round(lm.w(1))==logitvnum, "MNLRMSError: Incorrect MNL version!");
    mnlallerrors(lm, xy, npoints, relcls, avgce, rms, avg, avgrel);
    result = avg;
    return result;
}


/*************************************************************************
Average relative error on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average relative error (error when estimating posterior probabilities).

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double mnlavgrelerror(logitmodel& lm, const ap::real_2d_array& xy, int ssize)
{
    double result;
    double relcls;
    double avgce;
    double rms;
    double avg;
    double avgrel;

    ap::ap_error::make_assertion(ap::round(lm.w(1))==logitvnum, "MNLRMSError: Incorrect MNL version!");
    mnlallerrors(lm, xy, ssize, relcls, avgce, rms, avg, avgrel);
    result = avgrel;
    return result;
}


/*************************************************************************
Classification error on test set = MNLRelClsError*NPoints

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
int mnlclserror(logitmodel& lm, const ap::real_2d_array& xy, int npoints)
{
    int result;
    int nvars;
    int nclasses;
    int i;
    int j;
    ap::real_1d_array workx;
    ap::real_1d_array worky;
    int nmax;

    ap::ap_error::make_assertion(ap::fp_eq(lm.w(1),logitvnum), "MNLClsError: unexpected model version");
    nvars = ap::round(lm.w(2));
    nclasses = ap::round(lm.w(3));
    workx.setbounds(0, nvars-1);
    worky.setbounds(0, nclasses-1);
    result = 0;
    for(i = 0; i <= npoints-1; i++)
    {
        
        //
        // Process
        //
        ap::vmove(&workx(0), 1, &xy(i, 0), 1, ap::vlen(0,nvars-1));
        mnlprocess(lm, workx, worky);
        
        //
        // Logit version of the answer
        //
        nmax = 0;
        for(j = 0; j <= nclasses-1; j++)
        {
            if( ap::fp_greater(worky(j),worky(nmax)) )
            {
                nmax = j;
            }
        }
        
        //
        // compare
        //
        if( nmax!=ap::round(xy(i,nvars)) )
        {
            result = result+1;
        }
    }
    return result;
}


/*************************************************************************
Internal subroutine. Places exponents of the anti-overflow shifted
internal linear outputs into the service part of the W array.
*************************************************************************/
static void mnliexp(ap::real_1d_array& w, const ap::real_1d_array& x)
{
    int nvars;
    int nclasses;
    int offs;
    int i;
    int i1;
    double v;
    double mx;

    ap::ap_error::make_assertion(ap::fp_eq(w(1),logitvnum), "LOGIT: unexpected model version");
    nvars = ap::round(w(2));
    nclasses = ap::round(w(3));
    offs = ap::round(w(4));
    i1 = offs+(nvars+1)*(nclasses-1);
    for(i = 0; i <= nclasses-2; i++)
    {
        v = ap::vdotproduct(&w(offs+i*(nvars+1)), 1, &x(0), 1, ap::vlen(offs+i*(nvars+1),offs+i*(nvars+1)+nvars-1));
        w(i1+i) = v+w(offs+i*(nvars+1)+nvars);
    }
    w(i1+nclasses-1) = 0;
    mx = 0;
    for(i = i1; i <= i1+nclasses-1; i++)
    {
        mx = ap::maxreal(mx, w(i));
    }
    for(i = i1; i <= i1+nclasses-1; i++)
    {
        w(i) = exp(w(i)-mx);
    }
}


/*************************************************************************
Calculation of all types of errors

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
static void mnlallerrors(logitmodel& lm,
     const ap::real_2d_array& xy,
     int npoints,
     double& relcls,
     double& avgce,
     double& rms,
     double& avg,
     double& avgrel)
{
    int nvars;
    int nclasses;
    int i;
    ap::real_1d_array buf;
    ap::real_1d_array workx;
    ap::real_1d_array y;
    ap::real_1d_array dy;

    ap::ap_error::make_assertion(ap::round(lm.w(1))==logitvnum, "MNL unit: Incorrect MNL version!");
    nvars = ap::round(lm.w(2));
    nclasses = ap::round(lm.w(3));
    workx.setbounds(0, nvars-1);
    y.setbounds(0, nclasses-1);
    dy.setbounds(0, 0);
    dserrallocate(nclasses, buf);
    for(i = 0; i <= npoints-1; i++)
    {
        ap::vmove(&workx(0), 1, &xy(i, 0), 1, ap::vlen(0,nvars-1));
        mnlprocess(lm, workx, y);
        dy(0) = xy(i,nvars);
        dserraccumulate(buf, y, dy);
    }
    dserrfinish(buf);
    relcls = buf(0);
    avgce = buf(1);
    rms = buf(2);
    avg = buf(3);
    avgrel = buf(4);
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
static void mnlmcsrch(const int& n,
     ap::real_1d_array& x,
     double& f,
     ap::real_1d_array& g,
     const ap::real_1d_array& s,
     double& stp,
     int& info,
     int& nfev,
     ap::real_1d_array& wa,
     logitmcstate& state,
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
                mnlmcstep(state.stx, state.fxm, state.dgxm, state.sty, state.fym, state.dgym, stp, state.fm, state.dgm, state.brackt, state.stmin, state.stmax, state.infoc);
                
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
                mnlmcstep(state.stx, state.fx, state.dgx, state.sty, state.fy, state.dgy, stp, f, state.dg, state.brackt, state.stmin, state.stmax, state.infoc);
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


static void mnlmcstep(double& stx,
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




