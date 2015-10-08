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
#include "mlptrain.h"

static const double mindecay = 0.001;

static void mlpkfoldcvgeneral(const multilayerperceptron& n,
     const ap::real_2d_array& xy,
     int npoints,
     double decay,
     int restarts,
     int foldscount,
     bool lmalgorithm,
     double wstep,
     int maxits,
     int& info,
     mlpreport& rep,
     mlpcvreport& cvrep);
static void mlpkfoldsplit(const ap::real_2d_array& xy,
     int npoints,
     int nclasses,
     int foldscount,
     bool stratifiedsplits,
     ap::integer_1d_array& folds);

/*************************************************************************
Neural network training  using  modified  Levenberg-Marquardt  with  exact
Hessian calculation and regularization. Subroutine trains  neural  network
with restarts from random positions. Algorithm is well  suited  for  small
and medium scale problems (hundreds of weights).

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay constant, >=0.001
                    Decay term 'Decay*||Weights||^2' is added to error
                    function.
                    If you don't know what Decay to choose, use 0.001.
    Restarts    -   number of restarts from random position, >0.
                    If you don't know what Restarts to choose, use 2.

OUTPUT PARAMETERS:
    Network     -   trained neural network.
    Info        -   return code:
                    * -9, if internal matrix inverse subroutine failed
                    * -2, if there is a point with class number
                          outside of [0..NOut-1].
                    * -1, if wrong parameters specified
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void mlptrainlm(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int npoints,
     double decay,
     int restarts,
     int& info,
     mlpreport& rep)
{
    int nin;
    int nout;
    int wcount;
    double lmftol;
    double lmsteptol;
    int i;
    int k;
    double v;
    double e;
    double enew;
    double xnorm2;
    double stepnorm;
    ap::real_1d_array g;
    ap::real_1d_array d;
    ap::real_2d_array h;
    ap::real_2d_array hmod;
    ap::real_2d_array z;
    bool spd;
    double nu;
    double lambda;
    double lambdaup;
    double lambdadown;
    minlbfgsreport internalrep;
    minlbfgsstate state;
    ap::real_1d_array x;
    ap::real_1d_array y;
    ap::real_1d_array wbase;
    ap::real_1d_array wdir;
    ap::real_1d_array wt;
    ap::real_1d_array wx;
    int pass;
    ap::real_1d_array wbest;
    double ebest;
    int invinfo;
    matinvreport invrep;
    int solverinfo;
    densesolverreport solverrep;

    mlpproperties(network, nin, nout, wcount);
    lambdaup = 10;
    lambdadown = 0.3;
    lmftol = 0.001;
    lmsteptol = 0.001;
    
    //
    // Test for inputs
    //
    if( npoints<=0||restarts<1 )
    {
        info = -1;
        return;
    }
    if( mlpissoftmax(network) )
    {
        for(i = 0; i <= npoints-1; i++)
        {
            if( ap::round(xy(i,nin))<0||ap::round(xy(i,nin))>=nout )
            {
                info = -2;
                return;
            }
        }
    }
    decay = ap::maxreal(decay, mindecay);
    info = 2;
    
    //
    // Initialize data
    //
    rep.ngrad = 0;
    rep.nhess = 0;
    rep.ncholesky = 0;
    
    //
    // General case.
    // Prepare task and network. Allocate space.
    //
    mlpinitpreprocessor(network, xy, npoints);
    g.setbounds(0, wcount-1);
    h.setbounds(0, wcount-1, 0, wcount-1);
    hmod.setbounds(0, wcount-1, 0, wcount-1);
    wbase.setbounds(0, wcount-1);
    wdir.setbounds(0, wcount-1);
    wbest.setbounds(0, wcount-1);
    wt.setbounds(0, wcount-1);
    wx.setbounds(0, wcount-1);
    ebest = ap::maxrealnumber;
    
    //
    // Multiple passes
    //
    for(pass = 1; pass <= restarts; pass++)
    {
        
        //
        // Initialize weights
        //
        mlprandomize(network);
        
        //
        // First stage of the hybrid algorithm: LBFGS
        //
        ap::vmove(&wbase(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1));
        minlbfgscreate(wcount, ap::minint(wcount, 5), wbase, state);
        minlbfgssetcond(state, double(0), double(0), double(0), ap::maxint(25, wcount));
        while(minlbfgsiteration(state))
        {
            
            //
            // gradient
            //
            ap::vmove(&network.weights(0), 1, &state.x(0), 1, ap::vlen(0,wcount-1));
            mlpgradbatch(network, xy, npoints, state.f, state.g);
            
            //
            // weight decay
            //
            v = ap::vdotproduct(&network.weights(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1));
            state.f = state.f+0.5*decay*v;
            ap::vadd(&state.g(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1), decay);
            
            //
            // next iteration
            //
            rep.ngrad = rep.ngrad+1;
        }
        minlbfgsresults(state, wbase, internalrep);
        ap::vmove(&network.weights(0), 1, &wbase(0), 1, ap::vlen(0,wcount-1));
        
        //
        // Second stage of the hybrid algorithm: LM
        //
        // Initialize H with identity matrix,
        // G with gradient,
        // E with regularized error.
        //
        mlphessianbatch(network, xy, npoints, e, g, h);
        v = ap::vdotproduct(&network.weights(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1));
        e = e+0.5*decay*v;
        ap::vadd(&g(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1), decay);
        for(k = 0; k <= wcount-1; k++)
        {
            h(k,k) = h(k,k)+decay;
        }
        rep.nhess = rep.nhess+1;
        lambda = 0.001;
        nu = 2;
        while(true)
        {
            
            //
            // 1. HMod = H+lambda*I
            // 2. Try to solve (H+Lambda*I)*dx = -g.
            //    Increase lambda if left part is not positive definite.
            //
            for(i = 0; i <= wcount-1; i++)
            {
                ap::vmove(&hmod(i, 0), 1, &h(i, 0), 1, ap::vlen(0,wcount-1));
                hmod(i,i) = hmod(i,i)+lambda;
            }
            spd = spdmatrixcholesky(hmod, wcount, true);
            rep.ncholesky = rep.ncholesky+1;
            if( !spd )
            {
                lambda = lambda*lambdaup*nu;
                nu = nu*2;
                continue;
            }
            spdmatrixcholeskysolve(hmod, wcount, true, g, solverinfo, solverrep, wdir);
            if( solverinfo<0 )
            {
                lambda = lambda*lambdaup*nu;
                nu = nu*2;
                continue;
            }
            ap::vmul(&wdir(0), 1, ap::vlen(0,wcount-1), -1);
            
            //
            // Lambda found.
            // 1. Save old w in WBase
            // 1. Test some stopping criterions
            // 2. If error(w+wdir)>error(w), increase lambda
            //
            ap::vadd(&network.weights(0), 1, &wdir(0), 1, ap::vlen(0,wcount-1));
            xnorm2 = ap::vdotproduct(&network.weights(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1));
            stepnorm = ap::vdotproduct(&wdir(0), 1, &wdir(0), 1, ap::vlen(0,wcount-1));
            stepnorm = sqrt(stepnorm);
            enew = mlperror(network, xy, npoints)+0.5*decay*xnorm2;
            if( ap::fp_less(stepnorm,lmsteptol*(1+sqrt(xnorm2))) )
            {
                break;
            }
            if( ap::fp_greater(enew,e) )
            {
                lambda = lambda*lambdaup*nu;
                nu = nu*2;
                continue;
            }
            
            //
            // Optimize using inv(cholesky(H)) as preconditioner
            //
            rmatrixtrinverse(hmod, wcount, true, false, invinfo, invrep);
            if( invinfo<=0 )
            {
                
                //
                // if matrix can't be inverted then exit with errors
                // TODO: make WCount steps in direction suggested by HMod
                //
                info = -9;
                return;
            }
            ap::vmove(&wbase(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1));
            for(i = 0; i <= wcount-1; i++)
            {
                wt(i) = 0;
            }
            minlbfgscreatex(wcount, wcount, wt, 1, state);
            minlbfgssetcond(state, double(0), double(0), double(0), 5);
            while(minlbfgsiteration(state))
            {
                
                //
                // gradient
                //
                for(i = 0; i <= wcount-1; i++)
                {
                    v = ap::vdotproduct(&state.x(i), 1, &hmod(i, i), 1, ap::vlen(i,wcount-1));
                    network.weights(i) = wbase(i)+v;
                }
                mlpgradbatch(network, xy, npoints, state.f, g);
                for(i = 0; i <= wcount-1; i++)
                {
                    state.g(i) = 0;
                }
                for(i = 0; i <= wcount-1; i++)
                {
                    v = g(i);
                    ap::vadd(&state.g(i), 1, &hmod(i, i), 1, ap::vlen(i,wcount-1), v);
                }
                
                //
                // weight decay
                // grad(x'*x) = A'*(x0+A*t)
                //
                v = ap::vdotproduct(&network.weights(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1));
                state.f = state.f+0.5*decay*v;
                for(i = 0; i <= wcount-1; i++)
                {
                    v = decay*network.weights(i);
                    ap::vadd(&state.g(i), 1, &hmod(i, i), 1, ap::vlen(i,wcount-1), v);
                }
                
                //
                // next iteration
                //
                rep.ngrad = rep.ngrad+1;
            }
            minlbfgsresults(state, wt, internalrep);
            
            //
            // Accept new position.
            // Calculate Hessian
            //
            for(i = 0; i <= wcount-1; i++)
            {
                v = ap::vdotproduct(&wt(i), 1, &hmod(i, i), 1, ap::vlen(i,wcount-1));
                network.weights(i) = wbase(i)+v;
            }
            mlphessianbatch(network, xy, npoints, e, g, h);
            v = ap::vdotproduct(&network.weights(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1));
            e = e+0.5*decay*v;
            ap::vadd(&g(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1), decay);
            for(k = 0; k <= wcount-1; k++)
            {
                h(k,k) = h(k,k)+decay;
            }
            rep.nhess = rep.nhess+1;
            
            //
            // Update lambda
            //
            lambda = lambda*lambdadown;
            nu = 2;
        }
        
        //
        // update WBest
        //
        v = ap::vdotproduct(&network.weights(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1));
        e = 0.5*decay*v+mlperror(network, xy, npoints);
        if( ap::fp_less(e,ebest) )
        {
            ebest = e;
            ap::vmove(&wbest(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1));
        }
    }
    
    //
    // copy WBest to output
    //
    ap::vmove(&network.weights(0), 1, &wbest(0), 1, ap::vlen(0,wcount-1));
}


/*************************************************************************
Neural  network  training  using  L-BFGS  algorithm  with  regularization.
Subroutine  trains  neural  network  with  restarts from random positions.
Algorithm  is  well  suited  for  problems  of  any dimensionality (memory
requirements and step complexity are linear by weights number).

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay constant, >=0.001
                    Decay term 'Decay*||Weights||^2' is added to error
                    function.
                    If you don't know what Decay to choose, use 0.001.
    Restarts    -   number of restarts from random position, >0.
                    If you don't know what Restarts to choose, use 2.
    WStep       -   stopping criterion. Algorithm stops if  step  size  is
                    less than WStep. Recommended value - 0.01.  Zero  step
                    size means stopping after MaxIts iterations.
    MaxIts      -   stopping   criterion.  Algorithm  stops  after  MaxIts
                    iterations (NOT gradient  calculations).  Zero  MaxIts
                    means stopping when step is sufficiently small.

OUTPUT PARAMETERS:
    Network     -   trained neural network.
    Info        -   return code:
                    * -8, if both WStep=0 and MaxIts=0
                    * -2, if there is a point with class number
                          outside of [0..NOut-1].
                    * -1, if wrong parameters specified
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report

  -- ALGLIB --
     Copyright 09.12.2007 by Bochkanov Sergey
*************************************************************************/
void mlptrainlbfgs(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int npoints,
     double decay,
     int restarts,
     double wstep,
     int maxits,
     int& info,
     mlpreport& rep)
{
    int i;
    int pass;
    int nin;
    int nout;
    int wcount;
    ap::real_1d_array w;
    ap::real_1d_array wbest;
    double e;
    double v;
    double ebest;
    minlbfgsreport internalrep;
    minlbfgsstate state;

    
    //
    // Test inputs, parse flags, read network geometry
    //
    if( ap::fp_eq(wstep,0)&&maxits==0 )
    {
        info = -8;
        return;
    }
    if( npoints<=0||restarts<1||ap::fp_less(wstep,0)||maxits<0 )
    {
        info = -1;
        return;
    }
    mlpproperties(network, nin, nout, wcount);
    if( mlpissoftmax(network) )
    {
        for(i = 0; i <= npoints-1; i++)
        {
            if( ap::round(xy(i,nin))<0||ap::round(xy(i,nin))>=nout )
            {
                info = -2;
                return;
            }
        }
    }
    decay = ap::maxreal(decay, mindecay);
    info = 2;
    
    //
    // Prepare
    //
    mlpinitpreprocessor(network, xy, npoints);
    w.setbounds(0, wcount-1);
    wbest.setbounds(0, wcount-1);
    ebest = ap::maxrealnumber;
    
    //
    // Multiple starts
    //
    rep.ncholesky = 0;
    rep.nhess = 0;
    rep.ngrad = 0;
    for(pass = 1; pass <= restarts; pass++)
    {
        
        //
        // Process
        //
        mlprandomize(network);
        ap::vmove(&w(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1));
        minlbfgscreate(wcount, ap::minint(wcount, 10), w, state);
        minlbfgssetcond(state, 0.0, 0.0, wstep, maxits);
        while(minlbfgsiteration(state))
        {
            ap::vmove(&network.weights(0), 1, &state.x(0), 1, ap::vlen(0,wcount-1));
            mlpgradnbatch(network, xy, npoints, state.f, state.g);
            v = ap::vdotproduct(&network.weights(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1));
            state.f = state.f+0.5*decay*v;
            ap::vadd(&state.g(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1), decay);
            rep.ngrad = rep.ngrad+1;
        }
        minlbfgsresults(state, w, internalrep);
        ap::vmove(&network.weights(0), 1, &w(0), 1, ap::vlen(0,wcount-1));
        
        //
        // Compare with best
        //
        v = ap::vdotproduct(&network.weights(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1));
        e = mlperrorn(network, xy, npoints)+0.5*decay*v;
        if( ap::fp_less(e,ebest) )
        {
            ap::vmove(&wbest(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1));
            ebest = e;
        }
    }
    
    //
    // The best network
    //
    ap::vmove(&network.weights(0), 1, &wbest(0), 1, ap::vlen(0,wcount-1));
}


/*************************************************************************
Neural network training using early stopping (base algorithm - L-BFGS with
regularization).

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry
    TrnXY       -   training set
    TrnSize     -   training set size
    ValXY       -   validation set
    ValSize     -   validation set size
    Decay       -   weight decay constant, >=0.001
                    Decay term 'Decay*||Weights||^2' is added to error
                    function.
                    If you don't know what Decay to choose, use 0.001.
    Restarts    -   number of restarts from random position, >0.
                    If you don't know what Restarts to choose, use 2.

OUTPUT PARAMETERS:
    Network     -   trained neural network.
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NOut-1].
                    * -1, if wrong parameters specified
                          (NPoints<0, Restarts<1, ...).
                    *  2, task has been solved, stopping  criterion  met -
                          sufficiently small step size.  Not expected  (we
                          use  EARLY  stopping)  but  possible  and not an
                          error.
                    *  6, task has been solved, stopping  criterion  met -
                          increasing of validation set error.
    Rep         -   training report

NOTE:

Algorithm stops if validation set error increases for  a  long  enough  or
step size is small enought  (there  are  task  where  validation  set  may
decrease for eternity). In any case solution returned corresponds  to  the
minimum of validation set error.

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void mlptraines(multilayerperceptron& network,
     const ap::real_2d_array& trnxy,
     int trnsize,
     const ap::real_2d_array& valxy,
     int valsize,
     double decay,
     int restarts,
     int& info,
     mlpreport& rep)
{
    int i;
    int pass;
    int nin;
    int nout;
    int wcount;
    ap::real_1d_array w;
    ap::real_1d_array wbest;
    double e;
    double v;
    double ebest;
    ap::real_1d_array wfinal;
    double efinal;
    int itbest;
    minlbfgsreport internalrep;
    minlbfgsstate state;
    double wstep;

    wstep = 0.001;
    
    //
    // Test inputs, parse flags, read network geometry
    //
    if( trnsize<=0||valsize<=0||restarts<1||ap::fp_less(decay,0) )
    {
        info = -1;
        return;
    }
    mlpproperties(network, nin, nout, wcount);
    if( mlpissoftmax(network) )
    {
        for(i = 0; i <= trnsize-1; i++)
        {
            if( ap::round(trnxy(i,nin))<0||ap::round(trnxy(i,nin))>=nout )
            {
                info = -2;
                return;
            }
        }
        for(i = 0; i <= valsize-1; i++)
        {
            if( ap::round(valxy(i,nin))<0||ap::round(valxy(i,nin))>=nout )
            {
                info = -2;
                return;
            }
        }
    }
    info = 2;
    
    //
    // Prepare
    //
    mlpinitpreprocessor(network, trnxy, trnsize);
    w.setbounds(0, wcount-1);
    wbest.setbounds(0, wcount-1);
    wfinal.setbounds(0, wcount-1);
    efinal = ap::maxrealnumber;
    for(i = 0; i <= wcount-1; i++)
    {
        wfinal(i) = 0;
    }
    
    //
    // Multiple starts
    //
    rep.ncholesky = 0;
    rep.nhess = 0;
    rep.ngrad = 0;
    for(pass = 1; pass <= restarts; pass++)
    {
        
        //
        // Process
        //
        mlprandomize(network);
        ebest = mlperror(network, valxy, valsize);
        ap::vmove(&wbest(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1));
        itbest = 0;
        ap::vmove(&w(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1));
        minlbfgscreate(wcount, ap::minint(wcount, 10), w, state);
        minlbfgssetcond(state, 0.0, 0.0, wstep, 0);
        minlbfgssetxrep(state, true);
        while(minlbfgsiteration(state))
        {
            
            //
            // Calculate gradient
            //
            ap::vmove(&network.weights(0), 1, &state.x(0), 1, ap::vlen(0,wcount-1));
            mlpgradnbatch(network, trnxy, trnsize, state.f, state.g);
            v = ap::vdotproduct(&network.weights(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1));
            state.f = state.f+0.5*decay*v;
            ap::vadd(&state.g(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1), decay);
            rep.ngrad = rep.ngrad+1;
            
            //
            // Validation set
            //
            if( state.xupdated )
            {
                ap::vmove(&network.weights(0), 1, &w(0), 1, ap::vlen(0,wcount-1));
                e = mlperror(network, valxy, valsize);
                if( ap::fp_less(e,ebest) )
                {
                    ebest = e;
                    ap::vmove(&wbest(0), 1, &network.weights(0), 1, ap::vlen(0,wcount-1));
                    itbest = internalrep.iterationscount;
                }
                if( internalrep.iterationscount>30&&ap::fp_greater(internalrep.iterationscount,1.5*itbest) )
                {
                    info = 6;
                    break;
                }
            }
        }
        minlbfgsresults(state, w, internalrep);
        
        //
        // Compare with final answer
        //
        if( ap::fp_less(ebest,efinal) )
        {
            ap::vmove(&wfinal(0), 1, &wbest(0), 1, ap::vlen(0,wcount-1));
            efinal = ebest;
        }
    }
    
    //
    // The best network
    //
    ap::vmove(&network.weights(0), 1, &wfinal(0), 1, ap::vlen(0,wcount-1));
}


/*************************************************************************
Cross-validation estimate of generalization error.

Base algorithm - L-BFGS.

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry.   Network is
                    not changed during cross-validation -  it is used only
                    as a representative of its architecture.
    XY          -   training set.
    SSize       -   training set size
    Decay       -   weight  decay, same as in MLPTrainLBFGS
    Restarts    -   number of restarts, >0.
                    restarts are counted for each partition separately, so
                    total number of restarts will be Restarts*FoldsCount.
    WStep       -   stopping criterion, same as in MLPTrainLBFGS
    MaxIts      -   stopping criterion, same as in MLPTrainLBFGS
    FoldsCount  -   number of folds in k-fold cross-validation,
                    2<=FoldsCount<=SSize.
                    recommended value: 10.

OUTPUT PARAMETERS:
    Info        -   return code, same as in MLPTrainLBFGS
    Rep         -   report, same as in MLPTrainLM/MLPTrainLBFGS
    CVRep       -   generalization error estimates

  -- ALGLIB --
     Copyright 09.12.2007 by Bochkanov Sergey
*************************************************************************/
void mlpkfoldcvlbfgs(const multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int npoints,
     double decay,
     int restarts,
     double wstep,
     int maxits,
     int foldscount,
     int& info,
     mlpreport& rep,
     mlpcvreport& cvrep)
{

    mlpkfoldcvgeneral(network, xy, npoints, decay, restarts, foldscount, false, wstep, maxits, info, rep, cvrep);
}


/*************************************************************************
Cross-validation estimate of generalization error.

Base algorithm - Levenberg-Marquardt.

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry.   Network is
                    not changed during cross-validation -  it is used only
                    as a representative of its architecture.
    XY          -   training set.
    SSize       -   training set size
    Decay       -   weight  decay, same as in MLPTrainLBFGS
    Restarts    -   number of restarts, >0.
                    restarts are counted for each partition separately, so
                    total number of restarts will be Restarts*FoldsCount.
    FoldsCount  -   number of folds in k-fold cross-validation,
                    2<=FoldsCount<=SSize.
                    recommended value: 10.

OUTPUT PARAMETERS:
    Info        -   return code, same as in MLPTrainLBFGS
    Rep         -   report, same as in MLPTrainLM/MLPTrainLBFGS
    CVRep       -   generalization error estimates

  -- ALGLIB --
     Copyright 09.12.2007 by Bochkanov Sergey
*************************************************************************/
void mlpkfoldcvlm(const multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int npoints,
     double decay,
     int restarts,
     int foldscount,
     int& info,
     mlpreport& rep,
     mlpcvreport& cvrep)
{

    mlpkfoldcvgeneral(network, xy, npoints, decay, restarts, foldscount, true, 0.0, 0, info, rep, cvrep);
}


/*************************************************************************
Internal cross-validation subroutine
*************************************************************************/
static void mlpkfoldcvgeneral(const multilayerperceptron& n,
     const ap::real_2d_array& xy,
     int npoints,
     double decay,
     int restarts,
     int foldscount,
     bool lmalgorithm,
     double wstep,
     int maxits,
     int& info,
     mlpreport& rep,
     mlpcvreport& cvrep)
{
    int i;
    int fold;
    int j;
    int k;
    multilayerperceptron network;
    int nin;
    int nout;
    int rowlen;
    int wcount;
    int nclasses;
    int tssize;
    int cvssize;
    ap::real_2d_array cvset;
    ap::real_2d_array testset;
    ap::integer_1d_array folds;
    int relcnt;
    mlpreport internalrep;
    ap::real_1d_array x;
    ap::real_1d_array y;

    
    //
    // Read network geometry, test parameters
    //
    mlpproperties(n, nin, nout, wcount);
    if( mlpissoftmax(n) )
    {
        nclasses = nout;
        rowlen = nin+1;
    }
    else
    {
        nclasses = -nout;
        rowlen = nin+nout;
    }
    if( npoints<=0||foldscount<2||foldscount>npoints )
    {
        info = -1;
        return;
    }
    mlpcopy(n, network);
    
    //
    // K-fold out cross-validation.
    // First, estimate generalization error
    //
    testset.setbounds(0, npoints-1, 0, rowlen-1);
    cvset.setbounds(0, npoints-1, 0, rowlen-1);
    x.setbounds(0, nin-1);
    y.setbounds(0, nout-1);
    mlpkfoldsplit(xy, npoints, nclasses, foldscount, false, folds);
    cvrep.relclserror = 0;
    cvrep.avgce = 0;
    cvrep.rmserror = 0;
    cvrep.avgerror = 0;
    cvrep.avgrelerror = 0;
    rep.ngrad = 0;
    rep.nhess = 0;
    rep.ncholesky = 0;
    relcnt = 0;
    for(fold = 0; fold <= foldscount-1; fold++)
    {
        
        //
        // Separate set
        //
        tssize = 0;
        cvssize = 0;
        for(i = 0; i <= npoints-1; i++)
        {
            if( folds(i)==fold )
            {
                ap::vmove(&testset(tssize, 0), 1, &xy(i, 0), 1, ap::vlen(0,rowlen-1));
                tssize = tssize+1;
            }
            else
            {
                ap::vmove(&cvset(cvssize, 0), 1, &xy(i, 0), 1, ap::vlen(0,rowlen-1));
                cvssize = cvssize+1;
            }
        }
        
        //
        // Train on CV training set
        //
        if( lmalgorithm )
        {
            mlptrainlm(network, cvset, cvssize, decay, restarts, info, internalrep);
        }
        else
        {
            mlptrainlbfgs(network, cvset, cvssize, decay, restarts, wstep, maxits, info, internalrep);
        }
        if( info<0 )
        {
            cvrep.relclserror = 0;
            cvrep.avgce = 0;
            cvrep.rmserror = 0;
            cvrep.avgerror = 0;
            cvrep.avgrelerror = 0;
            return;
        }
        rep.ngrad = rep.ngrad+internalrep.ngrad;
        rep.nhess = rep.nhess+internalrep.nhess;
        rep.ncholesky = rep.ncholesky+internalrep.ncholesky;
        
        //
        // Estimate error using CV test set
        //
        if( mlpissoftmax(network) )
        {
            
            //
            // classification-only code
            //
            cvrep.relclserror = cvrep.relclserror+mlpclserror(network, testset, tssize);
            cvrep.avgce = cvrep.avgce+mlperrorn(network, testset, tssize);
        }
        for(i = 0; i <= tssize-1; i++)
        {
            ap::vmove(&x(0), 1, &testset(i, 0), 1, ap::vlen(0,nin-1));
            mlpprocess(network, x, y);
            if( mlpissoftmax(network) )
            {
                
                //
                // Classification-specific code
                //
                k = ap::round(testset(i,nin));
                for(j = 0; j <= nout-1; j++)
                {
                    if( j==k )
                    {
                        cvrep.rmserror = cvrep.rmserror+ap::sqr(y(j)-1);
                        cvrep.avgerror = cvrep.avgerror+fabs(y(j)-1);
                        cvrep.avgrelerror = cvrep.avgrelerror+fabs(y(j)-1);
                        relcnt = relcnt+1;
                    }
                    else
                    {
                        cvrep.rmserror = cvrep.rmserror+ap::sqr(y(j));
                        cvrep.avgerror = cvrep.avgerror+fabs(y(j));
                    }
                }
            }
            else
            {
                
                //
                // Regression-specific code
                //
                for(j = 0; j <= nout-1; j++)
                {
                    cvrep.rmserror = cvrep.rmserror+ap::sqr(y(j)-testset(i,nin+j));
                    cvrep.avgerror = cvrep.avgerror+fabs(y(j)-testset(i,nin+j));
                    if( ap::fp_neq(testset(i,nin+j),0) )
                    {
                        cvrep.avgrelerror = cvrep.avgrelerror+fabs((y(j)-testset(i,nin+j))/testset(i,nin+j));
                        relcnt = relcnt+1;
                    }
                }
            }
        }
    }
    if( mlpissoftmax(network) )
    {
        cvrep.relclserror = cvrep.relclserror/npoints;
        cvrep.avgce = cvrep.avgce/(log(double(2))*npoints);
    }
    cvrep.rmserror = sqrt(cvrep.rmserror/(npoints*nout));
    cvrep.avgerror = cvrep.avgerror/(npoints*nout);
    cvrep.avgrelerror = cvrep.avgrelerror/relcnt;
    info = 1;
}


/*************************************************************************
Subroutine prepares K-fold split of the training set.

NOTES:
    "NClasses>0" means that we have classification task.
    "NClasses<0" means regression task with -NClasses real outputs.
*************************************************************************/
static void mlpkfoldsplit(const ap::real_2d_array& xy,
     int npoints,
     int nclasses,
     int foldscount,
     bool stratifiedsplits,
     ap::integer_1d_array& folds)
{
    int i;
    int j;
    int k;

    
    //
    // test parameters
    //
    ap::ap_error::make_assertion(npoints>0, "MLPKFoldSplit: wrong NPoints!");
    ap::ap_error::make_assertion(nclasses>1||nclasses<0, "MLPKFoldSplit: wrong NClasses!");
    ap::ap_error::make_assertion(foldscount>=2&&foldscount<=npoints, "MLPKFoldSplit: wrong FoldsCount!");
    ap::ap_error::make_assertion(!stratifiedsplits, "MLPKFoldSplit: stratified splits are not supported!");
    
    //
    // Folds
    //
    folds.setbounds(0, npoints-1);
    for(i = 0; i <= npoints-1; i++)
    {
        folds(i) = i*foldscount/npoints;
    }
    for(i = 0; i <= npoints-2; i++)
    {
        j = i+ap::randominteger(npoints-i);
        if( j!=i )
        {
            k = folds(i);
            folds(i) = folds(j);
            folds(j) = k;
        }
    }
}




