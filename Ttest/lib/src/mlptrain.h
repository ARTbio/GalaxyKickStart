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

#ifndef _mlptrain_h
#define _mlptrain_h

#include "ap.h"
#include "ialglib.h"

#include "mlpbase.h"
#include "reflections.h"
#include "creflections.h"
#include "hqrnd.h"
#include "matgen.h"
#include "ablasf.h"
#include "ablas.h"
#include "trfac.h"
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "matinv.h"
#include "linmin.h"
#include "minlbfgs.h"
#include "hblas.h"
#include "sblas.h"
#include "ortfac.h"
#include "blas.h"
#include "rotations.h"
#include "bdsvd.h"
#include "svd.h"
#include "xblas.h"
#include "densesolver.h"


/*************************************************************************
Training report:
    * NGrad     - number of gradient calculations
    * NHess     - number of Hessian calculations
    * NCholesky - number of Cholesky decompositions
*************************************************************************/
struct mlpreport
{
    int ngrad;
    int nhess;
    int ncholesky;
};


/*************************************************************************
Cross-validation estimates of generalization error
*************************************************************************/
struct mlpcvreport
{
    double relclserror;
    double avgce;
    double rmserror;
    double avgerror;
    double avgrelerror;
};




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
     mlpreport& rep);


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
     mlpreport& rep);


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
     mlpreport& rep);


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
     mlpcvreport& cvrep);


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
     mlpcvreport& cvrep);


#endif

