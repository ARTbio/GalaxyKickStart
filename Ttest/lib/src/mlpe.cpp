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
#include "mlpe.h"

static const int mlpntotaloffset = 3;
static const int mlpevnum = 9;

static void mlpeallerrors(mlpensemble& ensemble,
     const ap::real_2d_array& xy,
     int npoints,
     double& relcls,
     double& avgce,
     double& rms,
     double& avg,
     double& avgrel);
static void mlpebagginginternal(mlpensemble& ensemble,
     const ap::real_2d_array& xy,
     int npoints,
     double decay,
     int restarts,
     double wstep,
     int maxits,
     bool lmalgorithm,
     int& info,
     mlpreport& rep,
     mlpcvreport& ooberrors);

/*************************************************************************
Like MLPCreate0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreate0(int nin, int nout, int ensemblesize, mlpensemble& ensemble)
{
    multilayerperceptron net;

    mlpcreate0(nin, nout, net);
    mlpecreatefromnetwork(net, ensemblesize, ensemble);
}


/*************************************************************************
Like MLPCreate1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreate1(int nin,
     int nhid,
     int nout,
     int ensemblesize,
     mlpensemble& ensemble)
{
    multilayerperceptron net;

    mlpcreate1(nin, nhid, nout, net);
    mlpecreatefromnetwork(net, ensemblesize, ensemble);
}


/*************************************************************************
Like MLPCreate2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreate2(int nin,
     int nhid1,
     int nhid2,
     int nout,
     int ensemblesize,
     mlpensemble& ensemble)
{
    multilayerperceptron net;

    mlpcreate2(nin, nhid1, nhid2, nout, net);
    mlpecreatefromnetwork(net, ensemblesize, ensemble);
}


/*************************************************************************
Like MLPCreateB0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreateb0(int nin,
     int nout,
     double b,
     double d,
     int ensemblesize,
     mlpensemble& ensemble)
{
    multilayerperceptron net;

    mlpcreateb0(nin, nout, b, d, net);
    mlpecreatefromnetwork(net, ensemblesize, ensemble);
}


/*************************************************************************
Like MLPCreateB1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreateb1(int nin,
     int nhid,
     int nout,
     double b,
     double d,
     int ensemblesize,
     mlpensemble& ensemble)
{
    multilayerperceptron net;

    mlpcreateb1(nin, nhid, nout, b, d, net);
    mlpecreatefromnetwork(net, ensemblesize, ensemble);
}


/*************************************************************************
Like MLPCreateB2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreateb2(int nin,
     int nhid1,
     int nhid2,
     int nout,
     double b,
     double d,
     int ensemblesize,
     mlpensemble& ensemble)
{
    multilayerperceptron net;

    mlpcreateb2(nin, nhid1, nhid2, nout, b, d, net);
    mlpecreatefromnetwork(net, ensemblesize, ensemble);
}


/*************************************************************************
Like MLPCreateR0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreater0(int nin,
     int nout,
     double a,
     double b,
     int ensemblesize,
     mlpensemble& ensemble)
{
    multilayerperceptron net;

    mlpcreater0(nin, nout, a, b, net);
    mlpecreatefromnetwork(net, ensemblesize, ensemble);
}


/*************************************************************************
Like MLPCreateR1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreater1(int nin,
     int nhid,
     int nout,
     double a,
     double b,
     int ensemblesize,
     mlpensemble& ensemble)
{
    multilayerperceptron net;

    mlpcreater1(nin, nhid, nout, a, b, net);
    mlpecreatefromnetwork(net, ensemblesize, ensemble);
}


/*************************************************************************
Like MLPCreateR2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreater2(int nin,
     int nhid1,
     int nhid2,
     int nout,
     double a,
     double b,
     int ensemblesize,
     mlpensemble& ensemble)
{
    multilayerperceptron net;

    mlpcreater2(nin, nhid1, nhid2, nout, a, b, net);
    mlpecreatefromnetwork(net, ensemblesize, ensemble);
}


/*************************************************************************
Like MLPCreateC0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatec0(int nin, int nout, int ensemblesize, mlpensemble& ensemble)
{
    multilayerperceptron net;

    mlpcreatec0(nin, nout, net);
    mlpecreatefromnetwork(net, ensemblesize, ensemble);
}


/*************************************************************************
Like MLPCreateC1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatec1(int nin,
     int nhid,
     int nout,
     int ensemblesize,
     mlpensemble& ensemble)
{
    multilayerperceptron net;

    mlpcreatec1(nin, nhid, nout, net);
    mlpecreatefromnetwork(net, ensemblesize, ensemble);
}


/*************************************************************************
Like MLPCreateC2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatec2(int nin,
     int nhid1,
     int nhid2,
     int nout,
     int ensemblesize,
     mlpensemble& ensemble)
{
    multilayerperceptron net;

    mlpcreatec2(nin, nhid1, nhid2, nout, net);
    mlpecreatefromnetwork(net, ensemblesize, ensemble);
}


/*************************************************************************
Creates ensemble from network. Only network geometry is copied.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatefromnetwork(const multilayerperceptron& network,
     int ensemblesize,
     mlpensemble& ensemble)
{
    int i;
    int ccount;

    ap::ap_error::make_assertion(ensemblesize>0, "MLPECreate: incorrect ensemble size!");
    
    //
    // network properties
    //
    mlpproperties(network, ensemble.nin, ensemble.nout, ensemble.wcount);
    if( mlpissoftmax(network) )
    {
        ccount = ensemble.nin;
    }
    else
    {
        ccount = ensemble.nin+ensemble.nout;
    }
    ensemble.postprocessing = false;
    ensemble.issoftmax = mlpissoftmax(network);
    ensemble.ensemblesize = ensemblesize;
    
    //
    // structure information
    //
    ensemble.structinfo.setbounds(0, network.structinfo(0)-1);
    for(i = 0; i <= network.structinfo(0)-1; i++)
    {
        ensemble.structinfo(i) = network.structinfo(i);
    }
    
    //
    // weights, means, sigmas
    //
    ensemble.weights.setbounds(0, ensemblesize*ensemble.wcount-1);
    ensemble.columnmeans.setbounds(0, ensemblesize*ccount-1);
    ensemble.columnsigmas.setbounds(0, ensemblesize*ccount-1);
    for(i = 0; i <= ensemblesize*ensemble.wcount-1; i++)
    {
        ensemble.weights(i) = ap::randomreal()-0.5;
    }
    for(i = 0; i <= ensemblesize-1; i++)
    {
        ap::vmove(&ensemble.columnmeans(i*ccount), 1, &network.columnmeans(0), 1, ap::vlen(i*ccount,(i+1)*ccount-1));
        ap::vmove(&ensemble.columnsigmas(i*ccount), 1, &network.columnsigmas(0), 1, ap::vlen(i*ccount,(i+1)*ccount-1));
    }
    
    //
    // serialized part
    //
    mlpserialize(network, ensemble.serializedmlp, ensemble.serializedlen);
    
    //
    // temporaries, internal buffers
    //
    ensemble.tmpweights.setbounds(0, ensemble.wcount-1);
    ensemble.tmpmeans.setbounds(0, ccount-1);
    ensemble.tmpsigmas.setbounds(0, ccount-1);
    ensemble.neurons.setbounds(0, ensemble.structinfo(mlpntotaloffset)-1);
    ensemble.dfdnet.setbounds(0, ensemble.structinfo(mlpntotaloffset)-1);
    ensemble.y.setbounds(0, ensemble.nout-1);
}


/*************************************************************************
Copying of MLPEnsemble strucure

INPUT PARAMETERS:
    Ensemble1 -   original

OUTPUT PARAMETERS:
    Ensemble2 -   copy

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecopy(const mlpensemble& ensemble1, mlpensemble& ensemble2)
{
    int i;
    int ssize;
    int ccount;
    int ntotal;

    
    //
    // Unload info
    //
    ssize = ensemble1.structinfo(0);
    if( ensemble1.issoftmax )
    {
        ccount = ensemble1.nin;
    }
    else
    {
        ccount = ensemble1.nin+ensemble1.nout;
    }
    ntotal = ensemble1.structinfo(mlpntotaloffset);
    
    //
    // Allocate space
    //
    ensemble2.structinfo.setbounds(0, ssize-1);
    ensemble2.weights.setbounds(0, ensemble1.ensemblesize*ensemble1.wcount-1);
    ensemble2.columnmeans.setbounds(0, ensemble1.ensemblesize*ccount-1);
    ensemble2.columnsigmas.setbounds(0, ensemble1.ensemblesize*ccount-1);
    ensemble2.tmpweights.setbounds(0, ensemble1.wcount-1);
    ensemble2.tmpmeans.setbounds(0, ccount-1);
    ensemble2.tmpsigmas.setbounds(0, ccount-1);
    ensemble2.serializedmlp.setbounds(0, ensemble1.serializedlen-1);
    ensemble2.neurons.setbounds(0, ntotal-1);
    ensemble2.dfdnet.setbounds(0, ntotal-1);
    ensemble2.y.setbounds(0, ensemble1.nout-1);
    
    //
    // Copy
    //
    ensemble2.nin = ensemble1.nin;
    ensemble2.nout = ensemble1.nout;
    ensemble2.wcount = ensemble1.wcount;
    ensemble2.ensemblesize = ensemble1.ensemblesize;
    ensemble2.issoftmax = ensemble1.issoftmax;
    ensemble2.postprocessing = ensemble1.postprocessing;
    ensemble2.serializedlen = ensemble1.serializedlen;
    for(i = 0; i <= ssize-1; i++)
    {
        ensemble2.structinfo(i) = ensemble1.structinfo(i);
    }
    ap::vmove(&ensemble2.weights(0), 1, &ensemble1.weights(0), 1, ap::vlen(0,ensemble1.ensemblesize*ensemble1.wcount-1));
    ap::vmove(&ensemble2.columnmeans(0), 1, &ensemble1.columnmeans(0), 1, ap::vlen(0,ensemble1.ensemblesize*ccount-1));
    ap::vmove(&ensemble2.columnsigmas(0), 1, &ensemble1.columnsigmas(0), 1, ap::vlen(0,ensemble1.ensemblesize*ccount-1));
    ap::vmove(&ensemble2.serializedmlp(0), 1, &ensemble1.serializedmlp(0), 1, ap::vlen(0,ensemble1.serializedlen-1));
}


/*************************************************************************
Serialization of MLPEnsemble strucure

INPUT PARAMETERS:
    Ensemble-   original

OUTPUT PARAMETERS:
    RA      -   array of real numbers which stores ensemble,
                array[0..RLen-1]
    RLen    -   RA lenght

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpeserialize(mlpensemble& ensemble, ap::real_1d_array& ra, int& rlen)
{
    int i;
    int ssize;
    int ntotal;
    int ccount;
    int hsize;
    int offs;

    hsize = 13;
    ssize = ensemble.structinfo(0);
    if( ensemble.issoftmax )
    {
        ccount = ensemble.nin;
    }
    else
    {
        ccount = ensemble.nin+ensemble.nout;
    }
    ntotal = ensemble.structinfo(mlpntotaloffset);
    rlen = hsize+ssize+ensemble.ensemblesize*ensemble.wcount+2*ccount*ensemble.ensemblesize+ensemble.serializedlen;
    
    //
    //  RA format:
    //  [0]     RLen
    //  [1]     Version (MLPEVNum)
    //  [2]     EnsembleSize
    //  [3]     NIn
    //  [4]     NOut
    //  [5]     WCount
    //  [6]     IsSoftmax 0/1
    //  [7]     PostProcessing 0/1
    //  [8]     sizeof(StructInfo)
    //  [9]     NTotal (sizeof(Neurons), sizeof(DFDNET))
    //  [10]    CCount (sizeof(ColumnMeans), sizeof(ColumnSigmas))
    //  [11]    data offset
    //  [12]    SerializedLen
    //
    //  [..]    StructInfo
    //  [..]    Weights
    //  [..]    ColumnMeans
    //  [..]    ColumnSigmas
    //
    ra.setbounds(0, rlen-1);
    ra(0) = rlen;
    ra(1) = mlpevnum;
    ra(2) = ensemble.ensemblesize;
    ra(3) = ensemble.nin;
    ra(4) = ensemble.nout;
    ra(5) = ensemble.wcount;
    if( ensemble.issoftmax )
    {
        ra(6) = 1;
    }
    else
    {
        ra(6) = 0;
    }
    if( ensemble.postprocessing )
    {
        ra(7) = 1;
    }
    else
    {
        ra(7) = 9;
    }
    ra(8) = ssize;
    ra(9) = ntotal;
    ra(10) = ccount;
    ra(11) = hsize;
    ra(12) = ensemble.serializedlen;
    offs = hsize;
    for(i = offs; i <= offs+ssize-1; i++)
    {
        ra(i) = ensemble.structinfo(i-offs);
    }
    offs = offs+ssize;
    ap::vmove(&ra(offs), 1, &ensemble.weights(0), 1, ap::vlen(offs,offs+ensemble.ensemblesize*ensemble.wcount-1));
    offs = offs+ensemble.ensemblesize*ensemble.wcount;
    ap::vmove(&ra(offs), 1, &ensemble.columnmeans(0), 1, ap::vlen(offs,offs+ensemble.ensemblesize*ccount-1));
    offs = offs+ensemble.ensemblesize*ccount;
    ap::vmove(&ra(offs), 1, &ensemble.columnsigmas(0), 1, ap::vlen(offs,offs+ensemble.ensemblesize*ccount-1));
    offs = offs+ensemble.ensemblesize*ccount;
    ap::vmove(&ra(offs), 1, &ensemble.serializedmlp(0), 1, ap::vlen(offs,offs+ensemble.serializedlen-1));
    offs = offs+ensemble.serializedlen;
}


/*************************************************************************
Unserialization of MLPEnsemble strucure

INPUT PARAMETERS:
    RA      -   real array which stores ensemble

OUTPUT PARAMETERS:
    Ensemble-   restored structure

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpeunserialize(const ap::real_1d_array& ra, mlpensemble& ensemble)
{
    int i;
    int ssize;
    int ntotal;
    int ccount;
    int hsize;
    int offs;

    ap::ap_error::make_assertion(ap::round(ra(1))==mlpevnum, "MLPEUnserialize: incorrect array!");
    
    //
    // load info
    //
    hsize = 13;
    ensemble.ensemblesize = ap::round(ra(2));
    ensemble.nin = ap::round(ra(3));
    ensemble.nout = ap::round(ra(4));
    ensemble.wcount = ap::round(ra(5));
    ensemble.issoftmax = ap::round(ra(6))==1;
    ensemble.postprocessing = ap::round(ra(7))==1;
    ssize = ap::round(ra(8));
    ntotal = ap::round(ra(9));
    ccount = ap::round(ra(10));
    offs = ap::round(ra(11));
    ensemble.serializedlen = ap::round(ra(12));
    
    //
    //  Allocate arrays
    //
    ensemble.structinfo.setbounds(0, ssize-1);
    ensemble.weights.setbounds(0, ensemble.ensemblesize*ensemble.wcount-1);
    ensemble.columnmeans.setbounds(0, ensemble.ensemblesize*ccount-1);
    ensemble.columnsigmas.setbounds(0, ensemble.ensemblesize*ccount-1);
    ensemble.tmpweights.setbounds(0, ensemble.wcount-1);
    ensemble.tmpmeans.setbounds(0, ccount-1);
    ensemble.tmpsigmas.setbounds(0, ccount-1);
    ensemble.neurons.setbounds(0, ntotal-1);
    ensemble.dfdnet.setbounds(0, ntotal-1);
    ensemble.serializedmlp.setbounds(0, ensemble.serializedlen-1);
    ensemble.y.setbounds(0, ensemble.nout-1);
    
    //
    // load data
    //
    for(i = offs; i <= offs+ssize-1; i++)
    {
        ensemble.structinfo(i-offs) = ap::round(ra(i));
    }
    offs = offs+ssize;
    ap::vmove(&ensemble.weights(0), 1, &ra(offs), 1, ap::vlen(0,ensemble.ensemblesize*ensemble.wcount-1));
    offs = offs+ensemble.ensemblesize*ensemble.wcount;
    ap::vmove(&ensemble.columnmeans(0), 1, &ra(offs), 1, ap::vlen(0,ensemble.ensemblesize*ccount-1));
    offs = offs+ensemble.ensemblesize*ccount;
    ap::vmove(&ensemble.columnsigmas(0), 1, &ra(offs), 1, ap::vlen(0,ensemble.ensemblesize*ccount-1));
    offs = offs+ensemble.ensemblesize*ccount;
    ap::vmove(&ensemble.serializedmlp(0), 1, &ra(offs), 1, ap::vlen(0,ensemble.serializedlen-1));
    offs = offs+ensemble.serializedlen;
}


/*************************************************************************
Randomization of MLP ensemble

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlperandomize(mlpensemble& ensemble)
{
    int i;

    for(i = 0; i <= ensemble.ensemblesize*ensemble.wcount-1; i++)
    {
        ensemble.weights(i) = ap::randomreal()-0.5;
    }
}


/*************************************************************************
Return ensemble properties (number of inputs and outputs).

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpeproperties(const mlpensemble& ensemble, int& nin, int& nout)
{

    nin = ensemble.nin;
    nout = ensemble.nout;
}


/*************************************************************************
Return normalization type (whether ensemble is SOFTMAX-normalized or not).

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
bool mlpeissoftmax(const mlpensemble& ensemble)
{
    bool result;

    result = ensemble.issoftmax;
    return result;
}


/*************************************************************************
Procesing

INPUT PARAMETERS:
    Ensemble-   neural networks ensemble
    X       -   input vector,  array[0..NIn-1].

OUTPUT PARAMETERS:
    Y       -   result. Regression estimate when solving regression  task,
                vector of posterior probabilities for classification task.
                Subroutine does not allocate memory for this vector, it is
                responsibility of a caller to allocate it. Array  must  be
                at least [0..NOut-1].

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpeprocess(mlpensemble& ensemble,
     const ap::real_1d_array& x,
     ap::real_1d_array& y)
{
    int i;
    int es;
    int wc;
    int cc;
    double v;

    es = ensemble.ensemblesize;
    wc = ensemble.wcount;
    if( ensemble.issoftmax )
    {
        cc = ensemble.nin;
    }
    else
    {
        cc = ensemble.nin+ensemble.nout;
    }
    v = double(1)/double(es);
    for(i = 0; i <= ensemble.nout-1; i++)
    {
        y(i) = 0;
    }
    for(i = 0; i <= es-1; i++)
    {
        ap::vmove(&ensemble.tmpweights(0), 1, &ensemble.weights(i*wc), 1, ap::vlen(0,wc-1));
        ap::vmove(&ensemble.tmpmeans(0), 1, &ensemble.columnmeans(i*cc), 1, ap::vlen(0,cc-1));
        ap::vmove(&ensemble.tmpsigmas(0), 1, &ensemble.columnsigmas(i*cc), 1, ap::vlen(0,cc-1));
        mlpinternalprocessvector(ensemble.structinfo, ensemble.tmpweights, ensemble.tmpmeans, ensemble.tmpsigmas, ensemble.neurons, ensemble.dfdnet, x, ensemble.y);
        ap::vadd(&y(0), 1, &ensemble.y(0), 1, ap::vlen(0,ensemble.nout-1), v);
    }
}


/*************************************************************************
Relative classification error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    percent of incorrectly classified cases.
    Works both for classifier betwork and for regression networks which
are used as classifiers.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlperelclserror(mlpensemble& ensemble,
     const ap::real_2d_array& xy,
     int npoints)
{
    double result;
    double relcls;
    double avgce;
    double rms;
    double avg;
    double avgrel;

    mlpeallerrors(ensemble, xy, npoints, relcls, avgce, rms, avg, avgrel);
    result = relcls;
    return result;
}


/*************************************************************************
Average cross-entropy (in bits per element) on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    CrossEntropy/(NPoints*LN(2)).
    Zero if ensemble solves regression task.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlpeavgce(mlpensemble& ensemble,
     const ap::real_2d_array& xy,
     int npoints)
{
    double result;
    double relcls;
    double avgce;
    double rms;
    double avg;
    double avgrel;

    mlpeallerrors(ensemble, xy, npoints, relcls, avgce, rms, avg, avgrel);
    result = avgce;
    return result;
}


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.
    Its meaning for regression task is obvious. As for classification task
RMS error means error when estimating posterior probabilities.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlpermserror(mlpensemble& ensemble,
     const ap::real_2d_array& xy,
     int npoints)
{
    double result;
    double relcls;
    double avgce;
    double rms;
    double avg;
    double avgrel;

    mlpeallerrors(ensemble, xy, npoints, relcls, avgce, rms, avg, avgrel);
    result = rms;
    return result;
}


/*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for classification task
it means average error when estimating posterior probabilities.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlpeavgerror(mlpensemble& ensemble,
     const ap::real_2d_array& xy,
     int npoints)
{
    double result;
    double relcls;
    double avgce;
    double rms;
    double avg;
    double avgrel;

    mlpeallerrors(ensemble, xy, npoints, relcls, avgce, rms, avg, avgrel);
    result = avg;
    return result;
}


/*************************************************************************
Average relative error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for classification task
it means average relative error when estimating posterior probabilities.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlpeavgrelerror(mlpensemble& ensemble,
     const ap::real_2d_array& xy,
     int npoints)
{
    double result;
    double relcls;
    double avgce;
    double rms;
    double avg;
    double avgrel;

    mlpeallerrors(ensemble, xy, npoints, relcls, avgce, rms, avg, avgrel);
    result = avgrel;
    return result;
}


/*************************************************************************
Training neural networks ensemble using  bootstrap  aggregating (bagging).
Modified Levenberg-Marquardt algorithm is used as base training method.

INPUT PARAMETERS:
    Ensemble    -   model with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay coefficient, >=0.001
    Restarts    -   restarts, >0.

OUTPUT PARAMETERS:
    Ensemble    -   trained model
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report.
    OOBErrors   -   out-of-bag generalization error estimate

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpebagginglm(mlpensemble& ensemble,
     const ap::real_2d_array& xy,
     int npoints,
     double decay,
     int restarts,
     int& info,
     mlpreport& rep,
     mlpcvreport& ooberrors)
{

    mlpebagginginternal(ensemble, xy, npoints, decay, restarts, 0.0, 0, true, info, rep, ooberrors);
}


/*************************************************************************
Training neural networks ensemble using  bootstrap  aggregating (bagging).
L-BFGS algorithm is used as base training method.

INPUT PARAMETERS:
    Ensemble    -   model with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay coefficient, >=0.001
    Restarts    -   restarts, >0.
    WStep       -   stopping criterion, same as in MLPTrainLBFGS
    MaxIts      -   stopping criterion, same as in MLPTrainLBFGS

OUTPUT PARAMETERS:
    Ensemble    -   trained model
    Info        -   return code:
                    * -8, if both WStep=0 and MaxIts=0
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report.
    OOBErrors   -   out-of-bag generalization error estimate

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpebagginglbfgs(mlpensemble& ensemble,
     const ap::real_2d_array& xy,
     int npoints,
     double decay,
     int restarts,
     double wstep,
     int maxits,
     int& info,
     mlpreport& rep,
     mlpcvreport& ooberrors)
{

    mlpebagginginternal(ensemble, xy, npoints, decay, restarts, wstep, maxits, false, info, rep, ooberrors);
}


/*************************************************************************
Training neural networks ensemble using early stopping.

INPUT PARAMETERS:
    Ensemble    -   model with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay coefficient, >=0.001
    Restarts    -   restarts, >0.

OUTPUT PARAMETERS:
    Ensemble    -   trained model
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<0, Restarts<1).
                    *  6, if task has been solved.
    Rep         -   training report.
    OOBErrors   -   out-of-bag generalization error estimate

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void mlpetraines(mlpensemble& ensemble,
     const ap::real_2d_array& xy,
     int npoints,
     double decay,
     int restarts,
     int& info,
     mlpreport& rep)
{
    int i;
    int k;
    int ccount;
    int pcount;
    ap::real_2d_array trnxy;
    ap::real_2d_array valxy;
    int trnsize;
    int valsize;
    multilayerperceptron network;
    int tmpinfo;
    mlpreport tmprep;

    if( npoints<2||restarts<1||ap::fp_less(decay,0) )
    {
        info = -1;
        return;
    }
    if( ensemble.issoftmax )
    {
        for(i = 0; i <= npoints-1; i++)
        {
            if( ap::round(xy(i,ensemble.nin))<0||ap::round(xy(i,ensemble.nin))>=ensemble.nout )
            {
                info = -2;
                return;
            }
        }
    }
    info = 6;
    
    //
    // allocate
    //
    if( ensemble.issoftmax )
    {
        ccount = ensemble.nin+1;
        pcount = ensemble.nin;
    }
    else
    {
        ccount = ensemble.nin+ensemble.nout;
        pcount = ensemble.nin+ensemble.nout;
    }
    trnxy.setbounds(0, npoints-1, 0, ccount-1);
    valxy.setbounds(0, npoints-1, 0, ccount-1);
    mlpunserialize(ensemble.serializedmlp, network);
    rep.ngrad = 0;
    rep.nhess = 0;
    rep.ncholesky = 0;
    
    //
    // train networks
    //
    for(k = 0; k <= ensemble.ensemblesize-1; k++)
    {
        
        //
        // Split set
        //
        do
        {
            trnsize = 0;
            valsize = 0;
            for(i = 0; i <= npoints-1; i++)
            {
                if( ap::fp_less(ap::randomreal(),0.66) )
                {
                    
                    //
                    // Assign sample to training set
                    //
                    ap::vmove(&trnxy(trnsize, 0), 1, &xy(i, 0), 1, ap::vlen(0,ccount-1));
                    trnsize = trnsize+1;
                }
                else
                {
                    
                    //
                    // Assign sample to validation set
                    //
                    ap::vmove(&valxy(valsize, 0), 1, &xy(i, 0), 1, ap::vlen(0,ccount-1));
                    valsize = valsize+1;
                }
            }
        }
        while(!(trnsize!=0&&valsize!=0));
        
        //
        // Train
        //
        mlptraines(network, trnxy, trnsize, valxy, valsize, decay, restarts, tmpinfo, tmprep);
        if( tmpinfo<0 )
        {
            info = tmpinfo;
            return;
        }
        
        //
        // save results
        //
        ap::vmove(&ensemble.weights(k*ensemble.wcount), 1, &network.weights(0), 1, ap::vlen(k*ensemble.wcount,(k+1)*ensemble.wcount-1));
        ap::vmove(&ensemble.columnmeans(k*pcount), 1, &network.columnmeans(0), 1, ap::vlen(k*pcount,(k+1)*pcount-1));
        ap::vmove(&ensemble.columnsigmas(k*pcount), 1, &network.columnsigmas(0), 1, ap::vlen(k*pcount,(k+1)*pcount-1));
        rep.ngrad = rep.ngrad+tmprep.ngrad;
        rep.nhess = rep.nhess+tmprep.nhess;
        rep.ncholesky = rep.ncholesky+tmprep.ncholesky;
    }
}


/*************************************************************************
Calculation of all types of errors

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
static void mlpeallerrors(mlpensemble& ensemble,
     const ap::real_2d_array& xy,
     int npoints,
     double& relcls,
     double& avgce,
     double& rms,
     double& avg,
     double& avgrel)
{
    int i;
    ap::real_1d_array buf;
    ap::real_1d_array workx;
    ap::real_1d_array y;
    ap::real_1d_array dy;

    workx.setbounds(0, ensemble.nin-1);
    y.setbounds(0, ensemble.nout-1);
    if( ensemble.issoftmax )
    {
        dy.setbounds(0, 0);
        dserrallocate(ensemble.nout, buf);
    }
    else
    {
        dy.setbounds(0, ensemble.nout-1);
        dserrallocate(-ensemble.nout, buf);
    }
    for(i = 0; i <= npoints-1; i++)
    {
        ap::vmove(&workx(0), 1, &xy(i, 0), 1, ap::vlen(0,ensemble.nin-1));
        mlpeprocess(ensemble, workx, y);
        if( ensemble.issoftmax )
        {
            dy(0) = xy(i,ensemble.nin);
        }
        else
        {
            ap::vmove(&dy(0), 1, &xy(i, ensemble.nin), 1, ap::vlen(0,ensemble.nout-1));
        }
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
Internal bagging subroutine.

  -- ALGLIB --
     Copyright 19.02.2009 by Bochkanov Sergey
*************************************************************************/
static void mlpebagginginternal(mlpensemble& ensemble,
     const ap::real_2d_array& xy,
     int npoints,
     double decay,
     int restarts,
     double wstep,
     int maxits,
     bool lmalgorithm,
     int& info,
     mlpreport& rep,
     mlpcvreport& ooberrors)
{
    ap::real_2d_array xys;
    ap::boolean_1d_array s;
    ap::real_2d_array oobbuf;
    ap::integer_1d_array oobcntbuf;
    ap::real_1d_array x;
    ap::real_1d_array y;
    ap::real_1d_array dy;
    ap::real_1d_array dsbuf;
    int nin;
    int nout;
    int ccnt;
    int pcnt;
    int i;
    int j;
    int k;
    double v;
    mlpreport tmprep;
    multilayerperceptron network;

    
    //
    // Test for inputs
    //
    if( !lmalgorithm&&ap::fp_eq(wstep,0)&&maxits==0 )
    {
        info = -8;
        return;
    }
    if( npoints<=0||restarts<1||ap::fp_less(wstep,0)||maxits<0 )
    {
        info = -1;
        return;
    }
    if( ensemble.issoftmax )
    {
        for(i = 0; i <= npoints-1; i++)
        {
            if( ap::round(xy(i,ensemble.nin))<0||ap::round(xy(i,ensemble.nin))>=ensemble.nout )
            {
                info = -2;
                return;
            }
        }
    }
    
    //
    // allocate temporaries
    //
    info = 2;
    rep.ngrad = 0;
    rep.nhess = 0;
    rep.ncholesky = 0;
    ooberrors.relclserror = 0;
    ooberrors.avgce = 0;
    ooberrors.rmserror = 0;
    ooberrors.avgerror = 0;
    ooberrors.avgrelerror = 0;
    nin = ensemble.nin;
    nout = ensemble.nout;
    if( ensemble.issoftmax )
    {
        ccnt = nin+1;
        pcnt = nin;
    }
    else
    {
        ccnt = nin+nout;
        pcnt = nin+nout;
    }
    xys.setbounds(0, npoints-1, 0, ccnt-1);
    s.setbounds(0, npoints-1);
    oobbuf.setbounds(0, npoints-1, 0, nout-1);
    oobcntbuf.setbounds(0, npoints-1);
    x.setbounds(0, nin-1);
    y.setbounds(0, nout-1);
    if( ensemble.issoftmax )
    {
        dy.setbounds(0, 0);
    }
    else
    {
        dy.setbounds(0, nout-1);
    }
    for(i = 0; i <= npoints-1; i++)
    {
        for(j = 0; j <= nout-1; j++)
        {
            oobbuf(i,j) = 0;
        }
    }
    for(i = 0; i <= npoints-1; i++)
    {
        oobcntbuf(i) = 0;
    }
    mlpunserialize(ensemble.serializedmlp, network);
    
    //
    // main bagging cycle
    //
    for(k = 0; k <= ensemble.ensemblesize-1; k++)
    {
        
        //
        // prepare dataset
        //
        for(i = 0; i <= npoints-1; i++)
        {
            s(i) = false;
        }
        for(i = 0; i <= npoints-1; i++)
        {
            j = ap::randominteger(npoints);
            s(j) = true;
            ap::vmove(&xys(i, 0), 1, &xy(j, 0), 1, ap::vlen(0,ccnt-1));
        }
        
        //
        // train
        //
        if( lmalgorithm )
        {
            mlptrainlm(network, xys, npoints, decay, restarts, info, tmprep);
        }
        else
        {
            mlptrainlbfgs(network, xys, npoints, decay, restarts, wstep, maxits, info, tmprep);
        }
        if( info<0 )
        {
            return;
        }
        
        //
        // save results
        //
        rep.ngrad = rep.ngrad+tmprep.ngrad;
        rep.nhess = rep.nhess+tmprep.nhess;
        rep.ncholesky = rep.ncholesky+tmprep.ncholesky;
        ap::vmove(&ensemble.weights(k*ensemble.wcount), 1, &network.weights(0), 1, ap::vlen(k*ensemble.wcount,(k+1)*ensemble.wcount-1));
        ap::vmove(&ensemble.columnmeans(k*pcnt), 1, &network.columnmeans(0), 1, ap::vlen(k*pcnt,(k+1)*pcnt-1));
        ap::vmove(&ensemble.columnsigmas(k*pcnt), 1, &network.columnsigmas(0), 1, ap::vlen(k*pcnt,(k+1)*pcnt-1));
        
        //
        // OOB estimates
        //
        for(i = 0; i <= npoints-1; i++)
        {
            if( !s(i) )
            {
                ap::vmove(&x(0), 1, &xy(i, 0), 1, ap::vlen(0,nin-1));
                mlpprocess(network, x, y);
                ap::vadd(&oobbuf(i, 0), 1, &y(0), 1, ap::vlen(0,nout-1));
                oobcntbuf(i) = oobcntbuf(i)+1;
            }
        }
    }
    
    //
    // OOB estimates
    //
    if( ensemble.issoftmax )
    {
        dserrallocate(nout, dsbuf);
    }
    else
    {
        dserrallocate(-nout, dsbuf);
    }
    for(i = 0; i <= npoints-1; i++)
    {
        if( oobcntbuf(i)!=0 )
        {
            v = double(1)/double(oobcntbuf(i));
            ap::vmove(&y(0), 1, &oobbuf(i, 0), 1, ap::vlen(0,nout-1), v);
            if( ensemble.issoftmax )
            {
                dy(0) = xy(i,nin);
            }
            else
            {
                ap::vmove(&dy(0), 1, &xy(i, nin), 1, ap::vlen(0,nout-1), v);
            }
            dserraccumulate(dsbuf, y, dy);
        }
    }
    dserrfinish(dsbuf);
    ooberrors.relclserror = dsbuf(0);
    ooberrors.avgce = dsbuf(1);
    ooberrors.rmserror = dsbuf(2);
    ooberrors.avgerror = dsbuf(3);
    ooberrors.avgrelerror = dsbuf(4);
}




