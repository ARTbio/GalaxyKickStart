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
#include "mlpbase.h"

static const int mlpvnum = 7;
static const int nfieldwidth = 4;
static const int chunksize = 32;

static void addinputlayer(int ncount,
     ap::integer_1d_array& lsizes,
     ap::integer_1d_array& ltypes,
     ap::integer_1d_array& lconnfirst,
     ap::integer_1d_array& lconnlast,
     int& lastproc);
static void addbiasedsummatorlayer(int ncount,
     ap::integer_1d_array& lsizes,
     ap::integer_1d_array& ltypes,
     ap::integer_1d_array& lconnfirst,
     ap::integer_1d_array& lconnlast,
     int& lastproc);
static void addactivationlayer(int functype,
     ap::integer_1d_array& lsizes,
     ap::integer_1d_array& ltypes,
     ap::integer_1d_array& lconnfirst,
     ap::integer_1d_array& lconnlast,
     int& lastproc);
static void addzerolayer(ap::integer_1d_array& lsizes,
     ap::integer_1d_array& ltypes,
     ap::integer_1d_array& lconnfirst,
     ap::integer_1d_array& lconnlast,
     int& lastproc);
static void mlpcreate(int nin,
     int nout,
     const ap::integer_1d_array& lsizes,
     const ap::integer_1d_array& ltypes,
     const ap::integer_1d_array& lconnfirst,
     const ap::integer_1d_array& lconnlast,
     int layerscount,
     bool isclsnet,
     multilayerperceptron& network);
static void mlpactivationfunction(double net,
     int k,
     double& f,
     double& df,
     double& d2f);
static void mlphessianbatchinternal(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int ssize,
     bool naturalerr,
     double& e,
     ap::real_1d_array& grad,
     ap::real_2d_array& h);
static void mlpinternalcalculategradient(multilayerperceptron& network,
     const ap::real_1d_array& neurons,
     const ap::real_1d_array& weights,
     ap::real_1d_array& derror,
     ap::real_1d_array& grad,
     bool naturalerrorfunc);
static void mlpchunkedgradient(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int cstart,
     int csize,
     double& e,
     ap::real_1d_array& grad,
     bool naturalerrorfunc);
static double safecrossentropy(double t, double z);

/*************************************************************************
Creates  neural  network  with  NIn  inputs,  NOut outputs, without hidden
layers, with linear output layer. Network weights are  filled  with  small
random values.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreate0(int nin, int nout, multilayerperceptron& network)
{
    ap::integer_1d_array lsizes;
    ap::integer_1d_array ltypes;
    ap::integer_1d_array lconnfirst;
    ap::integer_1d_array lconnlast;
    int layerscount;
    int lastproc;

    layerscount = 1+2;
    
    //
    // Allocate arrays
    //
    lsizes.setbounds(0, layerscount-1);
    ltypes.setbounds(0, layerscount-1);
    lconnfirst.setbounds(0, layerscount-1);
    lconnlast.setbounds(0, layerscount-1);
    
    //
    // Layers
    //
    addinputlayer(nin, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nout, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    
    //
    // Create
    //
    mlpcreate(nin, nout, lsizes, ltypes, lconnfirst, lconnlast, layerscount, false, network);
}


/*************************************************************************
Same  as  MLPCreate0,  but  with  one  hidden  layer  (NHid  neurons) with
non-linear activation function. Output layer is linear.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreate1(int nin, int nhid, int nout, multilayerperceptron& network)
{
    ap::integer_1d_array lsizes;
    ap::integer_1d_array ltypes;
    ap::integer_1d_array lconnfirst;
    ap::integer_1d_array lconnlast;
    int layerscount;
    int lastproc;

    layerscount = 1+3+2;
    
    //
    // Allocate arrays
    //
    lsizes.setbounds(0, layerscount-1);
    ltypes.setbounds(0, layerscount-1);
    lconnfirst.setbounds(0, layerscount-1);
    lconnlast.setbounds(0, layerscount-1);
    
    //
    // Layers
    //
    addinputlayer(nin, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nhid, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addactivationlayer(1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nout, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    
    //
    // Create
    //
    mlpcreate(nin, nout, lsizes, ltypes, lconnfirst, lconnlast, layerscount, false, network);
}


/*************************************************************************
Same as MLPCreate0, but with two hidden layers (NHid1 and  NHid2  neurons)
with non-linear activation function. Output layer is linear.
 $ALL

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreate2(int nin,
     int nhid1,
     int nhid2,
     int nout,
     multilayerperceptron& network)
{
    ap::integer_1d_array lsizes;
    ap::integer_1d_array ltypes;
    ap::integer_1d_array lconnfirst;
    ap::integer_1d_array lconnlast;
    int layerscount;
    int lastproc;

    layerscount = 1+3+3+2;
    
    //
    // Allocate arrays
    //
    lsizes.setbounds(0, layerscount-1);
    ltypes.setbounds(0, layerscount-1);
    lconnfirst.setbounds(0, layerscount-1);
    lconnlast.setbounds(0, layerscount-1);
    
    //
    // Layers
    //
    addinputlayer(nin, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nhid1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addactivationlayer(1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nhid2, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addactivationlayer(1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nout, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    
    //
    // Create
    //
    mlpcreate(nin, nout, lsizes, ltypes, lconnfirst, lconnlast, layerscount, false, network);
}


/*************************************************************************
Creates  neural  network  with  NIn  inputs,  NOut outputs, without hidden
layers with non-linear output layer. Network weights are filled with small
random values.

Activation function of the output layer takes values:

    (B, +INF), if D>=0

or

    (-INF, B), if D<0.


  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreateb0(int nin,
     int nout,
     double b,
     double d,
     multilayerperceptron& network)
{
    ap::integer_1d_array lsizes;
    ap::integer_1d_array ltypes;
    ap::integer_1d_array lconnfirst;
    ap::integer_1d_array lconnlast;
    int layerscount;
    int lastproc;
    int i;

    layerscount = 1+3;
    if( ap::fp_greater_eq(d,0) )
    {
        d = 1;
    }
    else
    {
        d = -1;
    }
    
    //
    // Allocate arrays
    //
    lsizes.setbounds(0, layerscount-1);
    ltypes.setbounds(0, layerscount-1);
    lconnfirst.setbounds(0, layerscount-1);
    lconnlast.setbounds(0, layerscount-1);
    
    //
    // Layers
    //
    addinputlayer(nin, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nout, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addactivationlayer(3, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    
    //
    // Create
    //
    mlpcreate(nin, nout, lsizes, ltypes, lconnfirst, lconnlast, layerscount, false, network);
    
    //
    // Turn on ouputs shift/scaling.
    //
    for(i = nin; i <= nin+nout-1; i++)
    {
        network.columnmeans(i) = b;
        network.columnsigmas(i) = d;
    }
}


/*************************************************************************
Same as MLPCreateB0 but with non-linear hidden layer.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreateb1(int nin,
     int nhid,
     int nout,
     double b,
     double d,
     multilayerperceptron& network)
{
    ap::integer_1d_array lsizes;
    ap::integer_1d_array ltypes;
    ap::integer_1d_array lconnfirst;
    ap::integer_1d_array lconnlast;
    int layerscount;
    int lastproc;
    int i;

    layerscount = 1+3+3;
    if( ap::fp_greater_eq(d,0) )
    {
        d = 1;
    }
    else
    {
        d = -1;
    }
    
    //
    // Allocate arrays
    //
    lsizes.setbounds(0, layerscount-1);
    ltypes.setbounds(0, layerscount-1);
    lconnfirst.setbounds(0, layerscount-1);
    lconnlast.setbounds(0, layerscount-1);
    
    //
    // Layers
    //
    addinputlayer(nin, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nhid, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addactivationlayer(1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nout, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addactivationlayer(3, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    
    //
    // Create
    //
    mlpcreate(nin, nout, lsizes, ltypes, lconnfirst, lconnlast, layerscount, false, network);
    
    //
    // Turn on ouputs shift/scaling.
    //
    for(i = nin; i <= nin+nout-1; i++)
    {
        network.columnmeans(i) = b;
        network.columnsigmas(i) = d;
    }
}


/*************************************************************************
Same as MLPCreateB0 but with two non-linear hidden layers.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreateb2(int nin,
     int nhid1,
     int nhid2,
     int nout,
     double b,
     double d,
     multilayerperceptron& network)
{
    ap::integer_1d_array lsizes;
    ap::integer_1d_array ltypes;
    ap::integer_1d_array lconnfirst;
    ap::integer_1d_array lconnlast;
    int layerscount;
    int lastproc;
    int i;

    layerscount = 1+3+3+3;
    if( ap::fp_greater_eq(d,0) )
    {
        d = 1;
    }
    else
    {
        d = -1;
    }
    
    //
    // Allocate arrays
    //
    lsizes.setbounds(0, layerscount-1);
    ltypes.setbounds(0, layerscount-1);
    lconnfirst.setbounds(0, layerscount-1);
    lconnlast.setbounds(0, layerscount-1);
    
    //
    // Layers
    //
    addinputlayer(nin, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nhid1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addactivationlayer(1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nhid2, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addactivationlayer(1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nout, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addactivationlayer(3, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    
    //
    // Create
    //
    mlpcreate(nin, nout, lsizes, ltypes, lconnfirst, lconnlast, layerscount, false, network);
    
    //
    // Turn on ouputs shift/scaling.
    //
    for(i = nin; i <= nin+nout-1; i++)
    {
        network.columnmeans(i) = b;
        network.columnsigmas(i) = d;
    }
}


/*************************************************************************
Creates  neural  network  with  NIn  inputs,  NOut outputs, without hidden
layers with non-linear output layer. Network weights are filled with small
random values. Activation function of the output layer takes values [A,B].

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreater0(int nin,
     int nout,
     double a,
     double b,
     multilayerperceptron& network)
{
    ap::integer_1d_array lsizes;
    ap::integer_1d_array ltypes;
    ap::integer_1d_array lconnfirst;
    ap::integer_1d_array lconnlast;
    int layerscount;
    int lastproc;
    int i;

    layerscount = 1+3;
    
    //
    // Allocate arrays
    //
    lsizes.setbounds(0, layerscount-1);
    ltypes.setbounds(0, layerscount-1);
    lconnfirst.setbounds(0, layerscount-1);
    lconnlast.setbounds(0, layerscount-1);
    
    //
    // Layers
    //
    addinputlayer(nin, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nout, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addactivationlayer(1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    
    //
    // Create
    //
    mlpcreate(nin, nout, lsizes, ltypes, lconnfirst, lconnlast, layerscount, false, network);
    
    //
    // Turn on outputs shift/scaling.
    //
    for(i = nin; i <= nin+nout-1; i++)
    {
        network.columnmeans(i) = 0.5*(a+b);
        network.columnsigmas(i) = 0.5*(a-b);
    }
}


/*************************************************************************
Same as MLPCreateR0, but with non-linear hidden layer.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreater1(int nin,
     int nhid,
     int nout,
     double a,
     double b,
     multilayerperceptron& network)
{
    ap::integer_1d_array lsizes;
    ap::integer_1d_array ltypes;
    ap::integer_1d_array lconnfirst;
    ap::integer_1d_array lconnlast;
    int layerscount;
    int lastproc;
    int i;

    layerscount = 1+3+3;
    
    //
    // Allocate arrays
    //
    lsizes.setbounds(0, layerscount-1);
    ltypes.setbounds(0, layerscount-1);
    lconnfirst.setbounds(0, layerscount-1);
    lconnlast.setbounds(0, layerscount-1);
    
    //
    // Layers
    //
    addinputlayer(nin, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nhid, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addactivationlayer(1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nout, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addactivationlayer(1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    
    //
    // Create
    //
    mlpcreate(nin, nout, lsizes, ltypes, lconnfirst, lconnlast, layerscount, false, network);
    
    //
    // Turn on outputs shift/scaling.
    //
    for(i = nin; i <= nin+nout-1; i++)
    {
        network.columnmeans(i) = 0.5*(a+b);
        network.columnsigmas(i) = 0.5*(a-b);
    }
}


/*************************************************************************
Same as MLPCreateR0, but with two non-linear hidden layers.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreater2(int nin,
     int nhid1,
     int nhid2,
     int nout,
     double a,
     double b,
     multilayerperceptron& network)
{
    ap::integer_1d_array lsizes;
    ap::integer_1d_array ltypes;
    ap::integer_1d_array lconnfirst;
    ap::integer_1d_array lconnlast;
    int layerscount;
    int lastproc;
    int i;

    layerscount = 1+3+3+3;
    
    //
    // Allocate arrays
    //
    lsizes.setbounds(0, layerscount-1);
    ltypes.setbounds(0, layerscount-1);
    lconnfirst.setbounds(0, layerscount-1);
    lconnlast.setbounds(0, layerscount-1);
    
    //
    // Layers
    //
    addinputlayer(nin, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nhid1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addactivationlayer(1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nhid2, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addactivationlayer(1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nout, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addactivationlayer(1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    
    //
    // Create
    //
    mlpcreate(nin, nout, lsizes, ltypes, lconnfirst, lconnlast, layerscount, false, network);
    
    //
    // Turn on outputs shift/scaling.
    //
    for(i = nin; i <= nin+nout-1; i++)
    {
        network.columnmeans(i) = 0.5*(a+b);
        network.columnsigmas(i) = 0.5*(a-b);
    }
}


/*************************************************************************
Creates classifier network with NIn  inputs  and  NOut  possible  classes.
Network contains no hidden layers and linear output  layer  with  SOFTMAX-
normalization  (so  outputs  sums  up  to  1.0  and  converge to posterior
probabilities).

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreatec0(int nin, int nout, multilayerperceptron& network)
{
    ap::integer_1d_array lsizes;
    ap::integer_1d_array ltypes;
    ap::integer_1d_array lconnfirst;
    ap::integer_1d_array lconnlast;
    int layerscount;
    int lastproc;

    ap::ap_error::make_assertion(nout>=2, "MLPCreateC0: NOut<2!");
    layerscount = 1+2+1;
    
    //
    // Allocate arrays
    //
    lsizes.setbounds(0, layerscount-1);
    ltypes.setbounds(0, layerscount-1);
    lconnfirst.setbounds(0, layerscount-1);
    lconnlast.setbounds(0, layerscount-1);
    
    //
    // Layers
    //
    addinputlayer(nin, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nout-1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addzerolayer(lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    
    //
    // Create
    //
    mlpcreate(nin, nout, lsizes, ltypes, lconnfirst, lconnlast, layerscount, true, network);
}


/*************************************************************************
Same as MLPCreateC0, but with one non-linear hidden layer.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreatec1(int nin, int nhid, int nout, multilayerperceptron& network)
{
    ap::integer_1d_array lsizes;
    ap::integer_1d_array ltypes;
    ap::integer_1d_array lconnfirst;
    ap::integer_1d_array lconnlast;
    int layerscount;
    int lastproc;

    ap::ap_error::make_assertion(nout>=2, "MLPCreateC1: NOut<2!");
    layerscount = 1+3+2+1;
    
    //
    // Allocate arrays
    //
    lsizes.setbounds(0, layerscount-1);
    ltypes.setbounds(0, layerscount-1);
    lconnfirst.setbounds(0, layerscount-1);
    lconnlast.setbounds(0, layerscount-1);
    
    //
    // Layers
    //
    addinputlayer(nin, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nhid, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addactivationlayer(1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nout-1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addzerolayer(lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    
    //
    // Create
    //
    mlpcreate(nin, nout, lsizes, ltypes, lconnfirst, lconnlast, layerscount, true, network);
}


/*************************************************************************
Same as MLPCreateC0, but with two non-linear hidden layers.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreatec2(int nin,
     int nhid1,
     int nhid2,
     int nout,
     multilayerperceptron& network)
{
    ap::integer_1d_array lsizes;
    ap::integer_1d_array ltypes;
    ap::integer_1d_array lconnfirst;
    ap::integer_1d_array lconnlast;
    int layerscount;
    int lastproc;

    ap::ap_error::make_assertion(nout>=2, "MLPCreateC2: NOut<2!");
    layerscount = 1+3+3+2+1;
    
    //
    // Allocate arrays
    //
    lsizes.setbounds(0, layerscount-1);
    ltypes.setbounds(0, layerscount-1);
    lconnfirst.setbounds(0, layerscount-1);
    lconnlast.setbounds(0, layerscount-1);
    
    //
    // Layers
    //
    addinputlayer(nin, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nhid1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addactivationlayer(1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nhid2, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addactivationlayer(1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addbiasedsummatorlayer(nout-1, lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    addzerolayer(lsizes, ltypes, lconnfirst, lconnlast, lastproc);
    
    //
    // Create
    //
    mlpcreate(nin, nout, lsizes, ltypes, lconnfirst, lconnlast, layerscount, true, network);
}


/*************************************************************************
Copying of neural network

INPUT PARAMETERS:
    Network1 -   original

OUTPUT PARAMETERS:
    Network2 -   copy

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcopy(const multilayerperceptron& network1,
     multilayerperceptron& network2)
{
    int i;
    int ssize;
    int ntotal;
    int nin;
    int nout;
    int wcount;

    
    //
    // Unload info
    //
    ssize = network1.structinfo(0);
    nin = network1.structinfo(1);
    nout = network1.structinfo(2);
    ntotal = network1.structinfo(3);
    wcount = network1.structinfo(4);
    
    //
    // Allocate space
    //
    network2.structinfo.setbounds(0, ssize-1);
    network2.weights.setbounds(0, wcount-1);
    if( mlpissoftmax(network1) )
    {
        network2.columnmeans.setbounds(0, nin-1);
        network2.columnsigmas.setbounds(0, nin-1);
    }
    else
    {
        network2.columnmeans.setbounds(0, nin+nout-1);
        network2.columnsigmas.setbounds(0, nin+nout-1);
    }
    network2.neurons.setbounds(0, ntotal-1);
    network2.chunks.setbounds(0, 3*ntotal, 0, chunksize-1);
    network2.nwbuf.setbounds(0, ap::maxint(wcount, 2*nout)-1);
    network2.dfdnet.setbounds(0, ntotal-1);
    network2.x.setbounds(0, nin-1);
    network2.y.setbounds(0, nout-1);
    network2.derror.setbounds(0, ntotal-1);
    
    //
    // Copy
    //
    for(i = 0; i <= ssize-1; i++)
    {
        network2.structinfo(i) = network1.structinfo(i);
    }
    ap::vmove(&network2.weights(0), 1, &network1.weights(0), 1, ap::vlen(0,wcount-1));
    if( mlpissoftmax(network1) )
    {
        ap::vmove(&network2.columnmeans(0), 1, &network1.columnmeans(0), 1, ap::vlen(0,nin-1));
        ap::vmove(&network2.columnsigmas(0), 1, &network1.columnsigmas(0), 1, ap::vlen(0,nin-1));
    }
    else
    {
        ap::vmove(&network2.columnmeans(0), 1, &network1.columnmeans(0), 1, ap::vlen(0,nin+nout-1));
        ap::vmove(&network2.columnsigmas(0), 1, &network1.columnsigmas(0), 1, ap::vlen(0,nin+nout-1));
    }
    ap::vmove(&network2.neurons(0), 1, &network1.neurons(0), 1, ap::vlen(0,ntotal-1));
    ap::vmove(&network2.dfdnet(0), 1, &network1.dfdnet(0), 1, ap::vlen(0,ntotal-1));
    ap::vmove(&network2.x(0), 1, &network1.x(0), 1, ap::vlen(0,nin-1));
    ap::vmove(&network2.y(0), 1, &network1.y(0), 1, ap::vlen(0,nout-1));
    ap::vmove(&network2.derror(0), 1, &network1.derror(0), 1, ap::vlen(0,ntotal-1));
}


/*************************************************************************
Serialization of MultiLayerPerceptron strucure

INPUT PARAMETERS:
    Network -   original

OUTPUT PARAMETERS:
    RA      -   array of real numbers which stores network,
                array[0..RLen-1]
    RLen    -   RA lenght

  -- ALGLIB --
     Copyright 29.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpserialize(const multilayerperceptron& network,
     ap::real_1d_array& ra,
     int& rlen)
{
    int i;
    int ssize;
    int ntotal;
    int nin;
    int nout;
    int wcount;
    int sigmalen;
    int offs;

    
    //
    // Unload info
    //
    ssize = network.structinfo(0);
    nin = network.structinfo(1);
    nout = network.structinfo(2);
    ntotal = network.structinfo(3);
    wcount = network.structinfo(4);
    if( mlpissoftmax(network) )
    {
        sigmalen = nin;
    }
    else
    {
        sigmalen = nin+nout;
    }
    
    //
    //  RA format:
    //      LEN         DESRC.
    //      1           RLen
    //      1           version (MLPVNum)
    //      1           StructInfo size
    //      SSize       StructInfo
    //      WCount      Weights
    //      SigmaLen    ColumnMeans
    //      SigmaLen    ColumnSigmas
    //
    rlen = 3+ssize+wcount+2*sigmalen;
    ra.setbounds(0, rlen-1);
    ra(0) = rlen;
    ra(1) = mlpvnum;
    ra(2) = ssize;
    offs = 3;
    for(i = 0; i <= ssize-1; i++)
    {
        ra(offs+i) = network.structinfo(i);
    }
    offs = offs+ssize;
    ap::vmove(&ra(offs), 1, &network.weights(0), 1, ap::vlen(offs,offs+wcount-1));
    offs = offs+wcount;
    ap::vmove(&ra(offs), 1, &network.columnmeans(0), 1, ap::vlen(offs,offs+sigmalen-1));
    offs = offs+sigmalen;
    ap::vmove(&ra(offs), 1, &network.columnsigmas(0), 1, ap::vlen(offs,offs+sigmalen-1));
    offs = offs+sigmalen;
}


/*************************************************************************
Unserialization of MultiLayerPerceptron strucure

INPUT PARAMETERS:
    RA      -   real array which stores network

OUTPUT PARAMETERS:
    Network -   restored network

  -- ALGLIB --
     Copyright 29.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpunserialize(const ap::real_1d_array& ra,
     multilayerperceptron& network)
{
    int i;
    int ssize;
    int ntotal;
    int nin;
    int nout;
    int wcount;
    int sigmalen;
    int offs;

    ap::ap_error::make_assertion(ap::round(ra(1))==mlpvnum, "MLPUnserialize: incorrect array!");
    
    //
    // Unload StructInfo from IA
    //
    offs = 3;
    ssize = ap::round(ra(2));
    network.structinfo.setbounds(0, ssize-1);
    for(i = 0; i <= ssize-1; i++)
    {
        network.structinfo(i) = ap::round(ra(offs+i));
    }
    offs = offs+ssize;
    
    //
    // Unload info from StructInfo
    //
    ssize = network.structinfo(0);
    nin = network.structinfo(1);
    nout = network.structinfo(2);
    ntotal = network.structinfo(3);
    wcount = network.structinfo(4);
    if( network.structinfo(6)==0 )
    {
        sigmalen = nin+nout;
    }
    else
    {
        sigmalen = nin;
    }
    
    //
    // Allocate space for other fields
    //
    network.weights.setbounds(0, wcount-1);
    network.columnmeans.setbounds(0, sigmalen-1);
    network.columnsigmas.setbounds(0, sigmalen-1);
    network.neurons.setbounds(0, ntotal-1);
    network.chunks.setbounds(0, 3*ntotal, 0, chunksize-1);
    network.nwbuf.setbounds(0, ap::maxint(wcount, 2*nout)-1);
    network.dfdnet.setbounds(0, ntotal-1);
    network.x.setbounds(0, nin-1);
    network.y.setbounds(0, nout-1);
    network.derror.setbounds(0, ntotal-1);
    
    //
    // Copy parameters from RA
    //
    ap::vmove(&network.weights(0), 1, &ra(offs), 1, ap::vlen(0,wcount-1));
    offs = offs+wcount;
    ap::vmove(&network.columnmeans(0), 1, &ra(offs), 1, ap::vlen(0,sigmalen-1));
    offs = offs+sigmalen;
    ap::vmove(&network.columnsigmas(0), 1, &ra(offs), 1, ap::vlen(0,sigmalen-1));
    offs = offs+sigmalen;
}


/*************************************************************************
Randomization of neural network weights

  -- ALGLIB --
     Copyright 06.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlprandomize(multilayerperceptron& network)
{
    int i;
    int nin;
    int nout;
    int wcount;

    mlpproperties(network, nin, nout, wcount);
    for(i = 0; i <= wcount-1; i++)
    {
        network.weights(i) = ap::randomreal()-0.5;
    }
}


/*************************************************************************
Randomization of neural network weights and standartisator

  -- ALGLIB --
     Copyright 10.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlprandomizefull(multilayerperceptron& network)
{
    int i;
    int nin;
    int nout;
    int wcount;
    int ntotal;
    int istart;
    int offs;
    int ntype;

    mlpproperties(network, nin, nout, wcount);
    ntotal = network.structinfo(3);
    istart = network.structinfo(5);
    
    //
    // Process network
    //
    for(i = 0; i <= wcount-1; i++)
    {
        network.weights(i) = ap::randomreal()-0.5;
    }
    for(i = 0; i <= nin-1; i++)
    {
        network.columnmeans(i) = 2*ap::randomreal()-1;
        network.columnsigmas(i) = 1.5*ap::randomreal()+0.5;
    }
    if( !mlpissoftmax(network) )
    {
        for(i = 0; i <= nout-1; i++)
        {
            offs = istart+(ntotal-nout+i)*nfieldwidth;
            ntype = network.structinfo(offs+0);
            if( ntype==0 )
            {
                
                //
                // Shifts are changed only for linear outputs neurons
                //
                network.columnmeans(nin+i) = 2*ap::randomreal()-1;
            }
            if( ntype==0||ntype==3 )
            {
                
                //
                // Scales are changed only for linear or bounded outputs neurons.
                // Note that scale randomization preserves sign.
                //
                network.columnsigmas(nin+i) = ap::sign(network.columnsigmas(nin+i))*(1.5*ap::randomreal()+0.5);
            }
        }
    }
}


/*************************************************************************
Internal subroutine.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpinitpreprocessor(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int ssize)
{
    int i;
    int j;
    int jmax;
    int nin;
    int nout;
    int wcount;
    int ntotal;
    int istart;
    int offs;
    int ntype;
    ap::real_1d_array means;
    ap::real_1d_array sigmas;
    double s;

    mlpproperties(network, nin, nout, wcount);
    ntotal = network.structinfo(3);
    istart = network.structinfo(5);
    
    //
    // Means/Sigmas
    //
    if( mlpissoftmax(network) )
    {
        jmax = nin-1;
    }
    else
    {
        jmax = nin+nout-1;
    }
    means.setbounds(0, jmax);
    sigmas.setbounds(0, jmax);
    for(j = 0; j <= jmax; j++)
    {
        means(j) = 0;
        for(i = 0; i <= ssize-1; i++)
        {
            means(j) = means(j)+xy(i,j);
        }
        means(j) = means(j)/ssize;
        sigmas(j) = 0;
        for(i = 0; i <= ssize-1; i++)
        {
            sigmas(j) = sigmas(j)+ap::sqr(xy(i,j)-means(j));
        }
        sigmas(j) = sqrt(sigmas(j)/ssize);
    }
    
    //
    // Inputs
    //
    for(i = 0; i <= nin-1; i++)
    {
        network.columnmeans(i) = means(i);
        network.columnsigmas(i) = sigmas(i);
        if( ap::fp_eq(network.columnsigmas(i),0) )
        {
            network.columnsigmas(i) = 1;
        }
    }
    
    //
    // Outputs
    //
    if( !mlpissoftmax(network) )
    {
        for(i = 0; i <= nout-1; i++)
        {
            offs = istart+(ntotal-nout+i)*nfieldwidth;
            ntype = network.structinfo(offs+0);
            
            //
            // Linear outputs
            //
            if( ntype==0 )
            {
                network.columnmeans(nin+i) = means(nin+i);
                network.columnsigmas(nin+i) = sigmas(nin+i);
                if( ap::fp_eq(network.columnsigmas(nin+i),0) )
                {
                    network.columnsigmas(nin+i) = 1;
                }
            }
            
            //
            // Bounded outputs (half-interval)
            //
            if( ntype==3 )
            {
                s = means(nin+i)-network.columnmeans(nin+i);
                if( ap::fp_eq(s,0) )
                {
                    s = ap::sign(network.columnsigmas(nin+i));
                }
                if( ap::fp_eq(s,0) )
                {
                    s = 1.0;
                }
                network.columnsigmas(nin+i) = ap::sign(network.columnsigmas(nin+i))*fabs(s);
                if( ap::fp_eq(network.columnsigmas(nin+i),0) )
                {
                    network.columnsigmas(nin+i) = 1;
                }
            }
        }
    }
}


/*************************************************************************
Returns information about initialized network: number of inputs, outputs,
weights.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpproperties(const multilayerperceptron& network,
     int& nin,
     int& nout,
     int& wcount)
{

    nin = network.structinfo(1);
    nout = network.structinfo(2);
    wcount = network.structinfo(4);
}


/*************************************************************************
Tells whether network is SOFTMAX-normalized (i.e. classifier) or not.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
bool mlpissoftmax(const multilayerperceptron& network)
{
    bool result;

    result = network.structinfo(6)==1;
    return result;
}


/*************************************************************************
Procesing

INPUT PARAMETERS:
    Network -   neural network
    X       -   input vector,  array[0..NIn-1].

OUTPUT PARAMETERS:
    Y       -   result. Regression estimate when solving regression  task,
                vector of posterior probabilities for classification task.
                Subroutine does not allocate memory for this vector, it is
                responsibility of a caller to allocate it. Array  must  be
                at least [0..NOut-1].

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpprocess(multilayerperceptron& network,
     const ap::real_1d_array& x,
     ap::real_1d_array& y)
{

    mlpinternalprocessvector(network.structinfo, network.weights, network.columnmeans, network.columnsigmas, network.neurons, network.dfdnet, x, y);
}


/*************************************************************************
Error function for neural network, internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
double mlperror(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int ssize)
{
    double result;
    int i;
    int k;
    int nin;
    int nout;
    int wcount;
    double e;

    mlpproperties(network, nin, nout, wcount);
    result = 0;
    for(i = 0; i <= ssize-1; i++)
    {
        ap::vmove(&network.x(0), 1, &xy(i, 0), 1, ap::vlen(0,nin-1));
        mlpprocess(network, network.x, network.y);
        if( mlpissoftmax(network) )
        {
            
            //
            // class labels outputs
            //
            k = ap::round(xy(i,nin));
            if( k>=0&&k<nout )
            {
                network.y(k) = network.y(k)-1;
            }
        }
        else
        {
            
            //
            // real outputs
            //
            ap::vsub(&network.y(0), 1, &xy(i, nin), 1, ap::vlen(0,nout-1));
        }
        e = ap::vdotproduct(&network.y(0), 1, &network.y(0), 1, ap::vlen(0,nout-1));
        result = result+e/2;
    }
    return result;
}


/*************************************************************************
Natural error function for neural network, internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
double mlperrorn(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int ssize)
{
    double result;
    int i;
    int k;
    int nin;
    int nout;
    int wcount;
    double e;

    mlpproperties(network, nin, nout, wcount);
    result = 0;
    for(i = 0; i <= ssize-1; i++)
    {
        
        //
        // Process vector
        //
        ap::vmove(&network.x(0), 1, &xy(i, 0), 1, ap::vlen(0,nin-1));
        mlpprocess(network, network.x, network.y);
        
        //
        // Update error function
        //
        if( network.structinfo(6)==0 )
        {
            
            //
            // Least squares error function
            //
            ap::vsub(&network.y(0), 1, &xy(i, nin), 1, ap::vlen(0,nout-1));
            e = ap::vdotproduct(&network.y(0), 1, &network.y(0), 1, ap::vlen(0,nout-1));
            result = result+e/2;
        }
        else
        {
            
            //
            // Cross-entropy error function
            //
            k = ap::round(xy(i,nin));
            if( k>=0&&k<nout )
            {
                result = result+safecrossentropy(double(1), network.y(k));
            }
        }
    }
    return result;
}


/*************************************************************************
Classification error

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
int mlpclserror(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int ssize)
{
    int result;
    int i;
    int j;
    int nin;
    int nout;
    int wcount;
    ap::real_1d_array workx;
    ap::real_1d_array worky;
    int nn;
    int ns;
    int nmax;

    mlpproperties(network, nin, nout, wcount);
    workx.setbounds(0, nin-1);
    worky.setbounds(0, nout-1);
    result = 0;
    for(i = 0; i <= ssize-1; i++)
    {
        
        //
        // Process
        //
        ap::vmove(&workx(0), 1, &xy(i, 0), 1, ap::vlen(0,nin-1));
        mlpprocess(network, workx, worky);
        
        //
        // Network version of the answer
        //
        nmax = 0;
        for(j = 0; j <= nout-1; j++)
        {
            if( ap::fp_greater(worky(j),worky(nmax)) )
            {
                nmax = j;
            }
        }
        nn = nmax;
        
        //
        // Right answer
        //
        if( mlpissoftmax(network) )
        {
            ns = ap::round(xy(i,nin));
        }
        else
        {
            nmax = 0;
            for(j = 0; j <= nout-1; j++)
            {
                if( ap::fp_greater(xy(i,nin+j),xy(i,nin+nmax)) )
                {
                    nmax = j;
                }
            }
            ns = nmax;
        }
        
        //
        // compare
        //
        if( nn!=ns )
        {
            result = result+1;
        }
    }
    return result;
}


/*************************************************************************
Relative classification error on the test set

INPUT PARAMETERS:
    Network -   network
    XY      -   test set
    NPoints -   test set size

RESULT:
    percent of incorrectly classified cases. Works both for
    classifier networks and general purpose networks used as
    classifiers.

  -- ALGLIB --
     Copyright 25.12.2008 by Bochkanov Sergey
*************************************************************************/
double mlprelclserror(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int npoints)
{
    double result;

    result = double(mlpclserror(network, xy, npoints))/double(npoints);
    return result;
}


/*************************************************************************
Average cross-entropy (in bits per element) on the test set

INPUT PARAMETERS:
    Network -   neural network
    XY      -   test set
    NPoints -   test set size

RESULT:
    CrossEntropy/(NPoints*LN(2)).
    Zero if network solves regression task.

  -- ALGLIB --
     Copyright 08.01.2009 by Bochkanov Sergey
*************************************************************************/
double mlpavgce(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int npoints)
{
    double result;
    int nin;
    int nout;
    int wcount;

    if( mlpissoftmax(network) )
    {
        mlpproperties(network, nin, nout, wcount);
        result = mlperrorn(network, xy, npoints)/(npoints*log(double(2)));
    }
    else
    {
        result = 0;
    }
    return result;
}


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    Network -   neural network
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.
    Its meaning for regression task is obvious. As for
    classification task, RMS error means error when estimating posterior
    probabilities.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
double mlprmserror(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int npoints)
{
    double result;
    int nin;
    int nout;
    int wcount;

    mlpproperties(network, nin, nout, wcount);
    result = sqrt(2*mlperror(network, xy, npoints)/(npoints*nout));
    return result;
}


/*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    Network -   neural network
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for
    classification task, it means average error when estimating posterior
    probabilities.

  -- ALGLIB --
     Copyright 11.03.2008 by Bochkanov Sergey
*************************************************************************/
double mlpavgerror(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int npoints)
{
    double result;
    int i;
    int j;
    int k;
    int nin;
    int nout;
    int wcount;

    mlpproperties(network, nin, nout, wcount);
    result = 0;
    for(i = 0; i <= npoints-1; i++)
    {
        ap::vmove(&network.x(0), 1, &xy(i, 0), 1, ap::vlen(0,nin-1));
        mlpprocess(network, network.x, network.y);
        if( mlpissoftmax(network) )
        {
            
            //
            // class labels
            //
            k = ap::round(xy(i,nin));
            for(j = 0; j <= nout-1; j++)
            {
                if( j==k )
                {
                    result = result+fabs(1-network.y(j));
                }
                else
                {
                    result = result+fabs(network.y(j));
                }
            }
        }
        else
        {
            
            //
            // real outputs
            //
            for(j = 0; j <= nout-1; j++)
            {
                result = result+fabs(xy(i,nin+j)-network.y(j));
            }
        }
    }
    result = result/(npoints*nout);
    return result;
}


/*************************************************************************
Average relative error on the test set

INPUT PARAMETERS:
    Network -   neural network
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for
    classification task, it means average relative error when estimating
    posterior probability of belonging to the correct class.

  -- ALGLIB --
     Copyright 11.03.2008 by Bochkanov Sergey
*************************************************************************/
double mlpavgrelerror(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int npoints)
{
    double result;
    int i;
    int j;
    int k;
    int lk;
    int nin;
    int nout;
    int wcount;

    mlpproperties(network, nin, nout, wcount);
    result = 0;
    k = 0;
    for(i = 0; i <= npoints-1; i++)
    {
        ap::vmove(&network.x(0), 1, &xy(i, 0), 1, ap::vlen(0,nin-1));
        mlpprocess(network, network.x, network.y);
        if( mlpissoftmax(network) )
        {
            
            //
            // class labels
            //
            lk = ap::round(xy(i,nin));
            for(j = 0; j <= nout-1; j++)
            {
                if( j==lk )
                {
                    result = result+fabs(1-network.y(j));
                    k = k+1;
                }
            }
        }
        else
        {
            
            //
            // real outputs
            //
            for(j = 0; j <= nout-1; j++)
            {
                if( ap::fp_neq(xy(i,nin+j),0) )
                {
                    result = result+fabs(xy(i,nin+j)-network.y(j))/fabs(xy(i,nin+j));
                    k = k+1;
                }
            }
        }
    }
    if( k!=0 )
    {
        result = result/k;
    }
    return result;
}


/*************************************************************************
Gradient calculation. Internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgrad(multilayerperceptron& network,
     const ap::real_1d_array& x,
     const ap::real_1d_array& desiredy,
     double& e,
     ap::real_1d_array& grad)
{
    int i;
    int nout;
    int ntotal;

    
    //
    // Prepare dError/dOut, internal structures
    //
    mlpprocess(network, x, network.y);
    nout = network.structinfo(2);
    ntotal = network.structinfo(3);
    e = 0;
    for(i = 0; i <= ntotal-1; i++)
    {
        network.derror(i) = 0;
    }
    for(i = 0; i <= nout-1; i++)
    {
        network.derror(ntotal-nout+i) = network.y(i)-desiredy(i);
        e = e+ap::sqr(network.y(i)-desiredy(i))/2;
    }
    
    //
    // gradient
    //
    mlpinternalcalculategradient(network, network.neurons, network.weights, network.derror, grad, false);
}


/*************************************************************************
Gradient calculation (natural error function). Internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgradn(multilayerperceptron& network,
     const ap::real_1d_array& x,
     const ap::real_1d_array& desiredy,
     double& e,
     ap::real_1d_array& grad)
{
    double s;
    int i;
    int nout;
    int ntotal;

    
    //
    // Prepare dError/dOut, internal structures
    //
    mlpprocess(network, x, network.y);
    nout = network.structinfo(2);
    ntotal = network.structinfo(3);
    for(i = 0; i <= ntotal-1; i++)
    {
        network.derror(i) = 0;
    }
    e = 0;
    if( network.structinfo(6)==0 )
    {
        
        //
        // Regression network, least squares
        //
        for(i = 0; i <= nout-1; i++)
        {
            network.derror(ntotal-nout+i) = network.y(i)-desiredy(i);
            e = e+ap::sqr(network.y(i)-desiredy(i))/2;
        }
    }
    else
    {
        
        //
        // Classification network, cross-entropy
        //
        s = 0;
        for(i = 0; i <= nout-1; i++)
        {
            s = s+desiredy(i);
        }
        for(i = 0; i <= nout-1; i++)
        {
            network.derror(ntotal-nout+i) = s*network.y(i)-desiredy(i);
            e = e+safecrossentropy(desiredy(i), network.y(i));
        }
    }
    
    //
    // gradient
    //
    mlpinternalcalculategradient(network, network.neurons, network.weights, network.derror, grad, true);
}


/*************************************************************************
Batch gradient calculation. Internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgradbatch(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int ssize,
     double& e,
     ap::real_1d_array& grad)
{
    int i;
    int nin;
    int nout;
    int wcount;

    mlpproperties(network, nin, nout, wcount);
    for(i = 0; i <= wcount-1; i++)
    {
        grad(i) = 0;
    }
    e = 0;
    i = 0;
    while(i<=ssize-1)
    {
        mlpchunkedgradient(network, xy, i, ap::minint(ssize, i+chunksize)-i, e, grad, false);
        i = i+chunksize;
    }
}


/*************************************************************************
Batch gradient calculation (natural error function). Internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgradnbatch(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int ssize,
     double& e,
     ap::real_1d_array& grad)
{
    int i;
    int nin;
    int nout;
    int wcount;

    mlpproperties(network, nin, nout, wcount);
    for(i = 0; i <= wcount-1; i++)
    {
        grad(i) = 0;
    }
    e = 0;
    i = 0;
    while(i<=ssize-1)
    {
        mlpchunkedgradient(network, xy, i, ap::minint(ssize, i+chunksize)-i, e, grad, true);
        i = i+chunksize;
    }
}


/*************************************************************************
Batch Hessian calculation (natural error function) using R-algorithm.
Internal subroutine.

  -- ALGLIB --
     Copyright 26.01.2008 by Bochkanov Sergey.
     
     Hessian calculation based on R-algorithm described in
     "Fast Exact Multiplication by the Hessian",
     B. A. Pearlmutter,
     Neural Computation, 1994.
*************************************************************************/
void mlphessiannbatch(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int ssize,
     double& e,
     ap::real_1d_array& grad,
     ap::real_2d_array& h)
{

    mlphessianbatchinternal(network, xy, ssize, true, e, grad, h);
}


/*************************************************************************
Batch Hessian calculation using R-algorithm.
Internal subroutine.

  -- ALGLIB --
     Copyright 26.01.2008 by Bochkanov Sergey.

     Hessian calculation based on R-algorithm described in
     "Fast Exact Multiplication by the Hessian",
     B. A. Pearlmutter,
     Neural Computation, 1994.
*************************************************************************/
void mlphessianbatch(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int ssize,
     double& e,
     ap::real_1d_array& grad,
     ap::real_2d_array& h)
{

    mlphessianbatchinternal(network, xy, ssize, false, e, grad, h);
}


/*************************************************************************
Internal subroutine, shouldn't be called by user.
*************************************************************************/
void mlpinternalprocessvector(const ap::integer_1d_array& structinfo,
     const ap::real_1d_array& weights,
     const ap::real_1d_array& columnmeans,
     const ap::real_1d_array& columnsigmas,
     ap::real_1d_array& neurons,
     ap::real_1d_array& dfdnet,
     const ap::real_1d_array& x,
     ap::real_1d_array& y)
{
    int i;
    int n1;
    int n2;
    int w1;
    int w2;
    int ntotal;
    int nin;
    int nout;
    int istart;
    int offs;
    double net;
    double f;
    double df;
    double d2f;
    double mx;
    bool perr;

    
    //
    // Read network geometry
    //
    nin = structinfo(1);
    nout = structinfo(2);
    ntotal = structinfo(3);
    istart = structinfo(5);
    
    //
    // Inputs standartisation and putting in the network
    //
    for(i = 0; i <= nin-1; i++)
    {
        if( ap::fp_neq(columnsigmas(i),0) )
        {
            neurons(i) = (x(i)-columnmeans(i))/columnsigmas(i);
        }
        else
        {
            neurons(i) = x(i)-columnmeans(i);
        }
    }
    
    //
    // Process network
    //
    for(i = 0; i <= ntotal-1; i++)
    {
        offs = istart+i*nfieldwidth;
        if( structinfo(offs+0)>0 )
        {
            
            //
            // Activation function
            //
            mlpactivationfunction(neurons(structinfo(offs+2)), structinfo(offs+0), f, df, d2f);
            neurons(i) = f;
            dfdnet(i) = df;
        }
        if( structinfo(offs+0)==0 )
        {
            
            //
            // Adaptive summator
            //
            n1 = structinfo(offs+2);
            n2 = n1+structinfo(offs+1)-1;
            w1 = structinfo(offs+3);
            w2 = w1+structinfo(offs+1)-1;
            net = ap::vdotproduct(&weights(w1), 1, &neurons(n1), 1, ap::vlen(w1,w2));
            neurons(i) = net;
            dfdnet(i) = 1.0;
        }
        if( structinfo(offs+0)<0 )
        {
            perr = true;
            if( structinfo(offs+0)==-2 )
            {
                
                //
                // input neuron, left unchanged
                //
                perr = false;
            }
            if( structinfo(offs+0)==-3 )
            {
                
                //
                // "-1" neuron
                //
                neurons(i) = -1;
                perr = false;
            }
            if( structinfo(offs+0)==-4 )
            {
                
                //
                // "0" neuron
                //
                neurons(i) = 0;
                perr = false;
            }
            ap::ap_error::make_assertion(!perr, "MLPInternalProcessVector: internal error - unknown neuron type!");
        }
    }
    
    //
    // Extract result
    //
    ap::vmove(&y(0), 1, &neurons(ntotal-nout), 1, ap::vlen(0,nout-1));
    
    //
    // Softmax post-processing or standardisation if needed
    //
    ap::ap_error::make_assertion(structinfo(6)==0||structinfo(6)==1, "MLPInternalProcessVector: unknown normalization type!");
    if( structinfo(6)==1 )
    {
        
        //
        // Softmax
        //
        mx = y(0);
        for(i = 1; i <= nout-1; i++)
        {
            mx = ap::maxreal(mx, y(i));
        }
        net = 0;
        for(i = 0; i <= nout-1; i++)
        {
            y(i) = exp(y(i)-mx);
            net = net+y(i);
        }
        for(i = 0; i <= nout-1; i++)
        {
            y(i) = y(i)/net;
        }
    }
    else
    {
        
        //
        // Standardisation
        //
        for(i = 0; i <= nout-1; i++)
        {
            y(i) = y(i)*columnsigmas(nin+i)+columnmeans(nin+i);
        }
    }
}


/*************************************************************************
Internal subroutine: adding new input layer to network
*************************************************************************/
static void addinputlayer(int ncount,
     ap::integer_1d_array& lsizes,
     ap::integer_1d_array& ltypes,
     ap::integer_1d_array& lconnfirst,
     ap::integer_1d_array& lconnlast,
     int& lastproc)
{

    lsizes(0) = ncount;
    ltypes(0) = -2;
    lconnfirst(0) = 0;
    lconnlast(0) = 0;
    lastproc = 0;
}


/*************************************************************************
Internal subroutine: adding new summator layer to network
*************************************************************************/
static void addbiasedsummatorlayer(int ncount,
     ap::integer_1d_array& lsizes,
     ap::integer_1d_array& ltypes,
     ap::integer_1d_array& lconnfirst,
     ap::integer_1d_array& lconnlast,
     int& lastproc)
{

    lsizes(lastproc+1) = 1;
    ltypes(lastproc+1) = -3;
    lconnfirst(lastproc+1) = 0;
    lconnlast(lastproc+1) = 0;
    lsizes(lastproc+2) = ncount;
    ltypes(lastproc+2) = 0;
    lconnfirst(lastproc+2) = lastproc;
    lconnlast(lastproc+2) = lastproc+1;
    lastproc = lastproc+2;
}


/*************************************************************************
Internal subroutine: adding new summator layer to network
*************************************************************************/
static void addactivationlayer(int functype,
     ap::integer_1d_array& lsizes,
     ap::integer_1d_array& ltypes,
     ap::integer_1d_array& lconnfirst,
     ap::integer_1d_array& lconnlast,
     int& lastproc)
{

    ap::ap_error::make_assertion(functype>0, "AddActivationLayer: incorrect function type");
    lsizes(lastproc+1) = lsizes(lastproc);
    ltypes(lastproc+1) = functype;
    lconnfirst(lastproc+1) = lastproc;
    lconnlast(lastproc+1) = lastproc;
    lastproc = lastproc+1;
}


/*************************************************************************
Internal subroutine: adding new zero layer to network
*************************************************************************/
static void addzerolayer(ap::integer_1d_array& lsizes,
     ap::integer_1d_array& ltypes,
     ap::integer_1d_array& lconnfirst,
     ap::integer_1d_array& lconnlast,
     int& lastproc)
{

    lsizes(lastproc+1) = 1;
    ltypes(lastproc+1) = -4;
    lconnfirst(lastproc+1) = 0;
    lconnlast(lastproc+1) = 0;
    lastproc = lastproc+1;
}


/*************************************************************************
Internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
static void mlpcreate(int nin,
     int nout,
     const ap::integer_1d_array& lsizes,
     const ap::integer_1d_array& ltypes,
     const ap::integer_1d_array& lconnfirst,
     const ap::integer_1d_array& lconnlast,
     int layerscount,
     bool isclsnet,
     multilayerperceptron& network)
{
    int i;
    int j;
    int ssize;
    int ntotal;
    int wcount;
    int offs;
    int nprocessed;
    int wallocated;
    ap::integer_1d_array localtemp;
    ap::integer_1d_array lnfirst;
    ap::integer_1d_array lnsyn;

    
    //
    // Check
    //
    ap::ap_error::make_assertion(layerscount>0, "MLPCreate: wrong parameters!");
    ap::ap_error::make_assertion(ltypes(0)==-2, "MLPCreate: wrong LTypes[0] (must be -2)!");
    for(i = 0; i <= layerscount-1; i++)
    {
        ap::ap_error::make_assertion(lsizes(i)>0, "MLPCreate: wrong LSizes!");
        ap::ap_error::make_assertion(lconnfirst(i)>=0&&(lconnfirst(i)<i||i==0), "MLPCreate: wrong LConnFirst!");
        ap::ap_error::make_assertion(lconnlast(i)>=lconnfirst(i)&&(lconnlast(i)<i||i==0), "MLPCreate: wrong LConnLast!");
    }
    
    //
    // Build network geometry
    //
    lnfirst.setbounds(0, layerscount-1);
    lnsyn.setbounds(0, layerscount-1);
    ntotal = 0;
    wcount = 0;
    for(i = 0; i <= layerscount-1; i++)
    {
        
        //
        // Analyze connections.
        // This code must throw an assertion in case of unknown LTypes[I]
        //
        lnsyn(i) = -1;
        if( ltypes(i)>=0 )
        {
            lnsyn(i) = 0;
            for(j = lconnfirst(i); j <= lconnlast(i); j++)
            {
                lnsyn(i) = lnsyn(i)+lsizes(j);
            }
        }
        else
        {
            if( ltypes(i)==-2||ltypes(i)==-3||ltypes(i)==-4 )
            {
                lnsyn(i) = 0;
            }
        }
        ap::ap_error::make_assertion(lnsyn(i)>=0, "MLPCreate: internal error #0!");
        
        //
        // Other info
        //
        lnfirst(i) = ntotal;
        ntotal = ntotal+lsizes(i);
        if( ltypes(i)==0 )
        {
            wcount = wcount+lnsyn(i)*lsizes(i);
        }
    }
    ssize = 7+ntotal*nfieldwidth;
    
    //
    // Allocate
    //
    network.structinfo.setbounds(0, ssize-1);
    network.weights.setbounds(0, wcount-1);
    if( isclsnet )
    {
        network.columnmeans.setbounds(0, nin-1);
        network.columnsigmas.setbounds(0, nin-1);
    }
    else
    {
        network.columnmeans.setbounds(0, nin+nout-1);
        network.columnsigmas.setbounds(0, nin+nout-1);
    }
    network.neurons.setbounds(0, ntotal-1);
    network.chunks.setbounds(0, 3*ntotal, 0, chunksize-1);
    network.nwbuf.setbounds(0, ap::maxint(wcount, 2*nout)-1);
    network.dfdnet.setbounds(0, ntotal-1);
    network.x.setbounds(0, nin-1);
    network.y.setbounds(0, nout-1);
    network.derror.setbounds(0, ntotal-1);
    
    //
    // Fill structure: global info
    //
    network.structinfo(0) = ssize;
    network.structinfo(1) = nin;
    network.structinfo(2) = nout;
    network.structinfo(3) = ntotal;
    network.structinfo(4) = wcount;
    network.structinfo(5) = 7;
    if( isclsnet )
    {
        network.structinfo(6) = 1;
    }
    else
    {
        network.structinfo(6) = 0;
    }
    
    //
    // Fill structure: neuron connections
    //
    nprocessed = 0;
    wallocated = 0;
    for(i = 0; i <= layerscount-1; i++)
    {
        for(j = 0; j <= lsizes(i)-1; j++)
        {
            offs = network.structinfo(5)+nprocessed*nfieldwidth;
            network.structinfo(offs+0) = ltypes(i);
            if( ltypes(i)==0 )
            {
                
                //
                // Adaptive summator:
                // * connections with weights to previous neurons
                //
                network.structinfo(offs+1) = lnsyn(i);
                network.structinfo(offs+2) = lnfirst(lconnfirst(i));
                network.structinfo(offs+3) = wallocated;
                wallocated = wallocated+lnsyn(i);
                nprocessed = nprocessed+1;
            }
            if( ltypes(i)>0 )
            {
                
                //
                // Activation layer:
                // * each neuron connected to one (only one) of previous neurons.
                // * no weights
                //
                network.structinfo(offs+1) = 1;
                network.structinfo(offs+2) = lnfirst(lconnfirst(i))+j;
                network.structinfo(offs+3) = -1;
                nprocessed = nprocessed+1;
            }
            if( ltypes(i)==-2||ltypes(i)==-3||ltypes(i)==-4 )
            {
                nprocessed = nprocessed+1;
            }
        }
    }
    ap::ap_error::make_assertion(wallocated==wcount, "MLPCreate: internal error #1!");
    ap::ap_error::make_assertion(nprocessed==ntotal, "MLPCreate: internal error #2!");
    
    //
    // Fill weights by small random values
    // Initialize means and sigmas
    //
    for(i = 0; i <= wcount-1; i++)
    {
        network.weights(i) = ap::randomreal()-0.5;
    }
    for(i = 0; i <= nin-1; i++)
    {
        network.columnmeans(i) = 0;
        network.columnsigmas(i) = 1;
    }
    if( !isclsnet )
    {
        for(i = 0; i <= nout-1; i++)
        {
            network.columnmeans(nin+i) = 0;
            network.columnsigmas(nin+i) = 1;
        }
    }
}


/*************************************************************************
Internal subroutine

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
static void mlpactivationfunction(double net,
     int k,
     double& f,
     double& df,
     double& d2f)
{
    double net2;
    double arg;
    double root;
    double r;

    f = 0;
    df = 0;
    if( k==1 )
    {
        
        //
        // TanH activation function
        //
        if( ap::fp_less(fabs(net),100) )
        {
            f = tanh(net);
        }
        else
        {
            f = ap::sign(net);
        }
        df = 1-ap::sqr(f);
        d2f = -2*f*df;
        return;
    }
    if( k==3 )
    {
        
        //
        // EX activation function
        //
        if( ap::fp_greater_eq(net,0) )
        {
            net2 = net*net;
            arg = net2+1;
            root = sqrt(arg);
            f = net+root;
            r = net/root;
            df = 1+r;
            d2f = (root-net*r)/arg;
        }
        else
        {
            f = exp(net);
            df = f;
            d2f = f;
        }
        return;
    }
    if( k==2 )
    {
        f = exp(-ap::sqr(net));
        df = -2*net*f;
        d2f = -2*(f+df*net);
        return;
    }
}


/*************************************************************************
Internal subroutine for Hessian calculation.

WARNING!!! Unspeakable math far beyong human capabilities :)
*************************************************************************/
static void mlphessianbatchinternal(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int ssize,
     bool naturalerr,
     double& e,
     ap::real_1d_array& grad,
     ap::real_2d_array& h)
{
    int nin;
    int nout;
    int wcount;
    int ntotal;
    int istart;
    int i;
    int j;
    int k;
    int kl;
    int offs;
    int n1;
    int n2;
    int w1;
    int w2;
    double s;
    double t;
    double v;
    double et;
    bool bflag;
    double f;
    double df;
    double d2f;
    double deidyj;
    double mx;
    double q;
    double z;
    double s2;
    double expi;
    double expj;
    ap::real_1d_array x;
    ap::real_1d_array desiredy;
    ap::real_1d_array gt;
    ap::real_1d_array zeros;
    ap::real_2d_array rx;
    ap::real_2d_array ry;
    ap::real_2d_array rdx;
    ap::real_2d_array rdy;

    mlpproperties(network, nin, nout, wcount);
    ntotal = network.structinfo(3);
    istart = network.structinfo(5);
    
    //
    // Prepare
    //
    x.setbounds(0, nin-1);
    desiredy.setbounds(0, nout-1);
    zeros.setbounds(0, wcount-1);
    gt.setbounds(0, wcount-1);
    rx.setbounds(0, ntotal+nout-1, 0, wcount-1);
    ry.setbounds(0, ntotal+nout-1, 0, wcount-1);
    rdx.setbounds(0, ntotal+nout-1, 0, wcount-1);
    rdy.setbounds(0, ntotal+nout-1, 0, wcount-1);
    e = 0;
    for(i = 0; i <= wcount-1; i++)
    {
        zeros(i) = 0;
    }
    ap::vmove(&grad(0), 1, &zeros(0), 1, ap::vlen(0,wcount-1));
    for(i = 0; i <= wcount-1; i++)
    {
        ap::vmove(&h(i, 0), 1, &zeros(0), 1, ap::vlen(0,wcount-1));
    }
    
    //
    // Process
    //
    for(k = 0; k <= ssize-1; k++)
    {
        
        //
        // Process vector with MLPGradN.
        // Now Neurons, DFDNET and DError contains results of the last run.
        //
        ap::vmove(&x(0), 1, &xy(k, 0), 1, ap::vlen(0,nin-1));
        if( mlpissoftmax(network) )
        {
            
            //
            // class labels outputs
            //
            kl = ap::round(xy(k,nin));
            for(i = 0; i <= nout-1; i++)
            {
                if( i==kl )
                {
                    desiredy(i) = 1;
                }
                else
                {
                    desiredy(i) = 0;
                }
            }
        }
        else
        {
            
            //
            // real outputs
            //
            ap::vmove(&desiredy(0), 1, &xy(k, nin), 1, ap::vlen(0,nout-1));
        }
        if( naturalerr )
        {
            mlpgradn(network, x, desiredy, et, gt);
        }
        else
        {
            mlpgrad(network, x, desiredy, et, gt);
        }
        
        //
        // grad, error
        //
        e = e+et;
        ap::vadd(&grad(0), 1, &gt(0), 1, ap::vlen(0,wcount-1));
        
        //
        // Hessian.
        // Forward pass of the R-algorithm
        //
        for(i = 0; i <= ntotal-1; i++)
        {
            offs = istart+i*nfieldwidth;
            ap::vmove(&rx(i, 0), 1, &zeros(0), 1, ap::vlen(0,wcount-1));
            ap::vmove(&ry(i, 0), 1, &zeros(0), 1, ap::vlen(0,wcount-1));
            if( network.structinfo(offs+0)>0 )
            {
                
                //
                // Activation function
                //
                n1 = network.structinfo(offs+2);
                ap::vmove(&rx(i, 0), 1, &ry(n1, 0), 1, ap::vlen(0,wcount-1));
                v = network.dfdnet(i);
                ap::vmove(&ry(i, 0), 1, &rx(i, 0), 1, ap::vlen(0,wcount-1), v);
            }
            if( network.structinfo(offs+0)==0 )
            {
                
                //
                // Adaptive summator
                //
                n1 = network.structinfo(offs+2);
                n2 = n1+network.structinfo(offs+1)-1;
                w1 = network.structinfo(offs+3);
                w2 = w1+network.structinfo(offs+1)-1;
                for(j = n1; j <= n2; j++)
                {
                    v = network.weights(w1+j-n1);
                    ap::vadd(&rx(i, 0), 1, &ry(j, 0), 1, ap::vlen(0,wcount-1), v);
                    rx(i,w1+j-n1) = rx(i,w1+j-n1)+network.neurons(j);
                }
                ap::vmove(&ry(i, 0), 1, &rx(i, 0), 1, ap::vlen(0,wcount-1));
            }
            if( network.structinfo(offs+0)<0 )
            {
                bflag = true;
                if( network.structinfo(offs+0)==-2 )
                {
                    
                    //
                    // input neuron, left unchanged
                    //
                    bflag = false;
                }
                if( network.structinfo(offs+0)==-3 )
                {
                    
                    //
                    // "-1" neuron, left unchanged
                    //
                    bflag = false;
                }
                if( network.structinfo(offs+0)==-4 )
                {
                    
                    //
                    // "0" neuron, left unchanged
                    //
                    bflag = false;
                }
                ap::ap_error::make_assertion(!bflag, "MLPHessianNBatch: internal error - unknown neuron type!");
            }
        }
        
        //
        // Hessian. Backward pass of the R-algorithm.
        //
        // Stage 1. Initialize RDY
        //
        for(i = 0; i <= ntotal+nout-1; i++)
        {
            ap::vmove(&rdy(i, 0), 1, &zeros(0), 1, ap::vlen(0,wcount-1));
        }
        if( network.structinfo(6)==0 )
        {
            
            //
            // Standardisation.
            //
            // In context of the Hessian calculation standardisation
            // is considered as additional layer with weightless
            // activation function:
            //
            // F(NET) := Sigma*NET
            //
            // So we add one more layer to forward pass, and
            // make forward/backward pass through this layer.
            //
            for(i = 0; i <= nout-1; i++)
            {
                n1 = ntotal-nout+i;
                n2 = ntotal+i;
                
                //
                // Forward pass from N1 to N2
                //
                ap::vmove(&rx(n2, 0), 1, &ry(n1, 0), 1, ap::vlen(0,wcount-1));
                v = network.columnsigmas(nin+i);
                ap::vmove(&ry(n2, 0), 1, &rx(n2, 0), 1, ap::vlen(0,wcount-1), v);
                
                //
                // Initialization of RDY
                //
                ap::vmove(&rdy(n2, 0), 1, &ry(n2, 0), 1, ap::vlen(0,wcount-1));
                
                //
                // Backward pass from N2 to N1:
                // 1. Calculate R(dE/dX).
                // 2. No R(dE/dWij) is needed since weight of activation neuron
                //    is fixed to 1. So we can update R(dE/dY) for
                //    the connected neuron (note that Vij=0, Wij=1)
                //
                df = network.columnsigmas(nin+i);
                ap::vmove(&rdx(n2, 0), 1, &rdy(n2, 0), 1, ap::vlen(0,wcount-1), df);
                ap::vadd(&rdy(n1, 0), 1, &rdx(n2, 0), 1, ap::vlen(0,wcount-1));
            }
        }
        else
        {
            
            //
            // Softmax.
            //
            // Initialize RDY using generalized expression for ei'(yi)
            // (see expression (9) from p. 5 of "Fast Exact Multiplication by the Hessian").
            //
            // When we are working with softmax network, generalized
            // expression for ei'(yi) is used because softmax
            // normalization leads to ei, which depends on all y's
            //
            if( naturalerr )
            {
                
                //
                // softmax + cross-entropy.
                // We have:
                //
                // S = sum(exp(yk)),
                // ei = sum(trn)*exp(yi)/S-trn_i
                //
                // j=i:   d(ei)/d(yj) = T*exp(yi)*(S-exp(yi))/S^2
                // j<>i:  d(ei)/d(yj) = -T*exp(yi)*exp(yj)/S^2
                //
                t = 0;
                for(i = 0; i <= nout-1; i++)
                {
                    t = t+desiredy(i);
                }
                mx = network.neurons(ntotal-nout);
                for(i = 0; i <= nout-1; i++)
                {
                    mx = ap::maxreal(mx, network.neurons(ntotal-nout+i));
                }
                s = 0;
                for(i = 0; i <= nout-1; i++)
                {
                    network.nwbuf(i) = exp(network.neurons(ntotal-nout+i)-mx);
                    s = s+network.nwbuf(i);
                }
                for(i = 0; i <= nout-1; i++)
                {
                    for(j = 0; j <= nout-1; j++)
                    {
                        if( j==i )
                        {
                            deidyj = t*network.nwbuf(i)*(s-network.nwbuf(i))/ap::sqr(s);
                            ap::vadd(&rdy(ntotal-nout+i, 0), 1, &ry(ntotal-nout+i, 0), 1, ap::vlen(0,wcount-1), deidyj);
                        }
                        else
                        {
                            deidyj = -t*network.nwbuf(i)*network.nwbuf(j)/ap::sqr(s);
                            ap::vadd(&rdy(ntotal-nout+i, 0), 1, &ry(ntotal-nout+j, 0), 1, ap::vlen(0,wcount-1), deidyj);
                        }
                    }
                }
            }
            else
            {
                
                //
                // For a softmax + squared error we have expression
                // far beyond human imagination so we dont even try
                // to comment on it. Just enjoy the code...
                //
                // P.S. That's why "natural error" is called "natural" -
                // compact beatiful expressions, fast code....
                //
                mx = network.neurons(ntotal-nout);
                for(i = 0; i <= nout-1; i++)
                {
                    mx = ap::maxreal(mx, network.neurons(ntotal-nout+i));
                }
                s = 0;
                s2 = 0;
                for(i = 0; i <= nout-1; i++)
                {
                    network.nwbuf(i) = exp(network.neurons(ntotal-nout+i)-mx);
                    s = s+network.nwbuf(i);
                    s2 = s2+ap::sqr(network.nwbuf(i));
                }
                q = 0;
                for(i = 0; i <= nout-1; i++)
                {
                    q = q+(network.y(i)-desiredy(i))*network.nwbuf(i);
                }
                for(i = 0; i <= nout-1; i++)
                {
                    z = -q+(network.y(i)-desiredy(i))*s;
                    expi = network.nwbuf(i);
                    for(j = 0; j <= nout-1; j++)
                    {
                        expj = network.nwbuf(j);
                        if( j==i )
                        {
                            deidyj = expi/ap::sqr(s)*((z+expi)*(s-2*expi)/s+expi*s2/ap::sqr(s));
                        }
                        else
                        {
                            deidyj = expi*expj/ap::sqr(s)*(s2/ap::sqr(s)-2*z/s-(expi+expj)/s+(network.y(i)-desiredy(i))-(network.y(j)-desiredy(j)));
                        }
                        ap::vadd(&rdy(ntotal-nout+i, 0), 1, &ry(ntotal-nout+j, 0), 1, ap::vlen(0,wcount-1), deidyj);
                    }
                }
            }
        }
        
        //
        // Hessian. Backward pass of the R-algorithm
        //
        // Stage 2. Process.
        //
        for(i = ntotal-1; i >= 0; i--)
        {
            
            //
            // Possible variants:
            // 1. Activation function
            // 2. Adaptive summator
            // 3. Special neuron
            //
            offs = istart+i*nfieldwidth;
            if( network.structinfo(offs+0)>0 )
            {
                n1 = network.structinfo(offs+2);
                
                //
                // First, calculate R(dE/dX).
                //
                mlpactivationfunction(network.neurons(n1), network.structinfo(offs+0), f, df, d2f);
                v = d2f*network.derror(i);
                ap::vmove(&rdx(i, 0), 1, &rdy(i, 0), 1, ap::vlen(0,wcount-1), df);
                ap::vadd(&rdx(i, 0), 1, &rx(i, 0), 1, ap::vlen(0,wcount-1), v);
                
                //
                // No R(dE/dWij) is needed since weight of activation neuron
                // is fixed to 1.
                //
                // So we can update R(dE/dY) for the connected neuron.
                // (note that Vij=0, Wij=1)
                //
                ap::vadd(&rdy(n1, 0), 1, &rdx(i, 0), 1, ap::vlen(0,wcount-1));
            }
            if( network.structinfo(offs+0)==0 )
            {
                
                //
                // Adaptive summator
                //
                n1 = network.structinfo(offs+2);
                n2 = n1+network.structinfo(offs+1)-1;
                w1 = network.structinfo(offs+3);
                w2 = w1+network.structinfo(offs+1)-1;
                
                //
                // First, calculate R(dE/dX).
                //
                ap::vmove(&rdx(i, 0), 1, &rdy(i, 0), 1, ap::vlen(0,wcount-1));
                
                //
                // Then, calculate R(dE/dWij)
                //
                for(j = w1; j <= w2; j++)
                {
                    v = network.neurons(n1+j-w1);
                    ap::vadd(&h(j, 0), 1, &rdx(i, 0), 1, ap::vlen(0,wcount-1), v);
                    v = network.derror(i);
                    ap::vadd(&h(j, 0), 1, &ry(n1+j-w1, 0), 1, ap::vlen(0,wcount-1), v);
                }
                
                //
                // And finally, update R(dE/dY) for connected neurons.
                //
                for(j = w1; j <= w2; j++)
                {
                    v = network.weights(j);
                    ap::vadd(&rdy(n1+j-w1, 0), 1, &rdx(i, 0), 1, ap::vlen(0,wcount-1), v);
                    rdy(n1+j-w1,j) = rdy(n1+j-w1,j)+network.derror(i);
                }
            }
            if( network.structinfo(offs+0)<0 )
            {
                bflag = false;
                if( network.structinfo(offs+0)==-2||network.structinfo(offs+0)==-3||network.structinfo(offs+0)==-4 )
                {
                    
                    //
                    // Special neuron type, no back-propagation required
                    //
                    bflag = true;
                }
                ap::ap_error::make_assertion(bflag, "MLPHessianNBatch: unknown neuron type!");
            }
        }
    }
}


/*************************************************************************
Internal subroutine

Network must be processed by MLPProcess on X
*************************************************************************/
static void mlpinternalcalculategradient(multilayerperceptron& network,
     const ap::real_1d_array& neurons,
     const ap::real_1d_array& weights,
     ap::real_1d_array& derror,
     ap::real_1d_array& grad,
     bool naturalerrorfunc)
{
    int i;
    int n1;
    int n2;
    int w1;
    int w2;
    int ntotal;
    int istart;
    int nin;
    int nout;
    int offs;
    double dedf;
    double dfdnet;
    double v;
    double fown;
    double deown;
    double net;
    double mx;
    bool bflag;

    
    //
    // Read network geometry
    //
    nin = network.structinfo(1);
    nout = network.structinfo(2);
    ntotal = network.structinfo(3);
    istart = network.structinfo(5);
    
    //
    // Pre-processing of dError/dOut:
    // from dError/dOut(normalized) to dError/dOut(non-normalized)
    //
    ap::ap_error::make_assertion(network.structinfo(6)==0||network.structinfo(6)==1, "MLPInternalCalculateGradient: unknown normalization type!");
    if( network.structinfo(6)==1 )
    {
        
        //
        // Softmax
        //
        if( !naturalerrorfunc )
        {
            mx = network.neurons(ntotal-nout);
            for(i = 0; i <= nout-1; i++)
            {
                mx = ap::maxreal(mx, network.neurons(ntotal-nout+i));
            }
            net = 0;
            for(i = 0; i <= nout-1; i++)
            {
                network.nwbuf(i) = exp(network.neurons(ntotal-nout+i)-mx);
                net = net+network.nwbuf(i);
            }
            v = ap::vdotproduct(&network.derror(ntotal-nout), 1, &network.nwbuf(0), 1, ap::vlen(ntotal-nout,ntotal-1));
            for(i = 0; i <= nout-1; i++)
            {
                fown = network.nwbuf(i);
                deown = network.derror(ntotal-nout+i);
                network.nwbuf(nout+i) = (-v+deown*fown+deown*(net-fown))*fown/ap::sqr(net);
            }
            for(i = 0; i <= nout-1; i++)
            {
                network.derror(ntotal-nout+i) = network.nwbuf(nout+i);
            }
        }
    }
    else
    {
        
        //
        // Un-standardisation
        //
        for(i = 0; i <= nout-1; i++)
        {
            network.derror(ntotal-nout+i) = network.derror(ntotal-nout+i)*network.columnsigmas(nin+i);
        }
    }
    
    //
    // Backpropagation
    //
    for(i = ntotal-1; i >= 0; i--)
    {
        
        //
        // Extract info
        //
        offs = istart+i*nfieldwidth;
        if( network.structinfo(offs+0)>0 )
        {
            
            //
            // Activation function
            //
            dedf = network.derror(i);
            dfdnet = network.dfdnet(i);
            derror(network.structinfo(offs+2)) = derror(network.structinfo(offs+2))+dedf*dfdnet;
        }
        if( network.structinfo(offs+0)==0 )
        {
            
            //
            // Adaptive summator
            //
            n1 = network.structinfo(offs+2);
            n2 = n1+network.structinfo(offs+1)-1;
            w1 = network.structinfo(offs+3);
            w2 = w1+network.structinfo(offs+1)-1;
            dedf = network.derror(i);
            dfdnet = 1.0;
            v = dedf*dfdnet;
            ap::vmove(&grad(w1), 1, &neurons(n1), 1, ap::vlen(w1,w2), v);
            ap::vadd(&derror(n1), 1, &weights(w1), 1, ap::vlen(n1,n2), v);
        }
        if( network.structinfo(offs+0)<0 )
        {
            bflag = false;
            if( network.structinfo(offs+0)==-2||network.structinfo(offs+0)==-3||network.structinfo(offs+0)==-4 )
            {
                
                //
                // Special neuron type, no back-propagation required
                //
                bflag = true;
            }
            ap::ap_error::make_assertion(bflag, "MLPInternalCalculateGradient: unknown neuron type!");
        }
    }
}


/*************************************************************************
Internal subroutine, chunked gradient
*************************************************************************/
static void mlpchunkedgradient(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int cstart,
     int csize,
     double& e,
     ap::real_1d_array& grad,
     bool naturalerrorfunc)
{
    int i;
    int j;
    int k;
    int kl;
    int n1;
    int n2;
    int w1;
    int w2;
    int c1;
    int c2;
    int ntotal;
    int nin;
    int nout;
    int offs;
    double f;
    double df;
    double d2f;
    double v;
    double s;
    double fown;
    double deown;
    double net;
    double lnnet;
    double mx;
    bool bflag;
    int istart;
    int ineurons;
    int idfdnet;
    int iderror;
    int izeros;

    
    //
    // Read network geometry, prepare data
    //
    nin = network.structinfo(1);
    nout = network.structinfo(2);
    ntotal = network.structinfo(3);
    istart = network.structinfo(5);
    c1 = cstart;
    c2 = cstart+csize-1;
    ineurons = 0;
    idfdnet = ntotal;
    iderror = 2*ntotal;
    izeros = 3*ntotal;
    for(j = 0; j <= csize-1; j++)
    {
        network.chunks(izeros,j) = 0;
    }
    
    //
    // Forward pass:
    // 1. Load inputs from XY to Chunks[0:NIn-1,0:CSize-1]
    // 2. Forward pass
    //
    for(i = 0; i <= nin-1; i++)
    {
        for(j = 0; j <= csize-1; j++)
        {
            if( ap::fp_neq(network.columnsigmas(i),0) )
            {
                network.chunks(i,j) = (xy(c1+j,i)-network.columnmeans(i))/network.columnsigmas(i);
            }
            else
            {
                network.chunks(i,j) = xy(c1+j,i)-network.columnmeans(i);
            }
        }
    }
    for(i = 0; i <= ntotal-1; i++)
    {
        offs = istart+i*nfieldwidth;
        if( network.structinfo(offs+0)>0 )
        {
            
            //
            // Activation function:
            // * calculate F vector, F(i) = F(NET(i))
            //
            n1 = network.structinfo(offs+2);
            ap::vmove(&network.chunks(i, 0), 1, &network.chunks(n1, 0), 1, ap::vlen(0,csize-1));
            for(j = 0; j <= csize-1; j++)
            {
                mlpactivationfunction(network.chunks(i,j), network.structinfo(offs+0), f, df, d2f);
                network.chunks(i,j) = f;
                network.chunks(idfdnet+i,j) = df;
            }
        }
        if( network.structinfo(offs+0)==0 )
        {
            
            //
            // Adaptive summator:
            // * calculate NET vector, NET(i) = SUM(W(j,i)*Neurons(j),j=N1..N2)
            //
            n1 = network.structinfo(offs+2);
            n2 = n1+network.structinfo(offs+1)-1;
            w1 = network.structinfo(offs+3);
            w2 = w1+network.structinfo(offs+1)-1;
            ap::vmove(&network.chunks(i, 0), 1, &network.chunks(izeros, 0), 1, ap::vlen(0,csize-1));
            for(j = n1; j <= n2; j++)
            {
                v = network.weights(w1+j-n1);
                ap::vadd(&network.chunks(i, 0), 1, &network.chunks(j, 0), 1, ap::vlen(0,csize-1), v);
            }
        }
        if( network.structinfo(offs+0)<0 )
        {
            bflag = false;
            if( network.structinfo(offs+0)==-2 )
            {
                
                //
                // input neuron, left unchanged
                //
                bflag = true;
            }
            if( network.structinfo(offs+0)==-3 )
            {
                
                //
                // "-1" neuron
                //
                for(k = 0; k <= csize-1; k++)
                {
                    network.chunks(i,k) = -1;
                }
                bflag = true;
            }
            if( network.structinfo(offs+0)==-4 )
            {
                
                //
                // "0" neuron
                //
                for(k = 0; k <= csize-1; k++)
                {
                    network.chunks(i,k) = 0;
                }
                bflag = true;
            }
            ap::ap_error::make_assertion(bflag, "MLPChunkedGradient: internal error - unknown neuron type!");
        }
    }
    
    //
    // Post-processing, error, dError/dOut
    //
    for(i = 0; i <= ntotal-1; i++)
    {
        ap::vmove(&network.chunks(iderror+i, 0), 1, &network.chunks(izeros, 0), 1, ap::vlen(0,csize-1));
    }
    ap::ap_error::make_assertion(network.structinfo(6)==0||network.structinfo(6)==1, "MLPChunkedGradient: unknown normalization type!");
    if( network.structinfo(6)==1 )
    {
        
        //
        // Softmax output, classification network.
        //
        // For each K = 0..CSize-1 do:
        // 1. place exp(outputs[k]) to NWBuf[0:NOut-1]
        // 2. place sum(exp(..)) to NET
        // 3. calculate dError/dOut and place it to the second block of Chunks
        //
        for(k = 0; k <= csize-1; k++)
        {
            
            //
            // Normalize
            //
            mx = network.chunks(ntotal-nout,k);
            for(i = 1; i <= nout-1; i++)
            {
                mx = ap::maxreal(mx, network.chunks(ntotal-nout+i,k));
            }
            net = 0;
            for(i = 0; i <= nout-1; i++)
            {
                network.nwbuf(i) = exp(network.chunks(ntotal-nout+i,k)-mx);
                net = net+network.nwbuf(i);
            }
            
            //
            // Calculate error function and dError/dOut
            //
            if( naturalerrorfunc )
            {
                
                //
                // Natural error func.
                //
                //
                s = 1;
                lnnet = log(net);
                kl = ap::round(xy(cstart+k,nin));
                for(i = 0; i <= nout-1; i++)
                {
                    if( i==kl )
                    {
                        v = 1;
                    }
                    else
                    {
                        v = 0;
                    }
                    network.chunks(iderror+ntotal-nout+i,k) = s*network.nwbuf(i)/net-v;
                    e = e+safecrossentropy(v, network.nwbuf(i)/net);
                }
            }
            else
            {
                
                //
                // Least squares error func
                // Error, dError/dOut(normalized)
                //
                kl = ap::round(xy(cstart+k,nin));
                for(i = 0; i <= nout-1; i++)
                {
                    if( i==kl )
                    {
                        v = network.nwbuf(i)/net-1;
                    }
                    else
                    {
                        v = network.nwbuf(i)/net;
                    }
                    network.nwbuf(nout+i) = v;
                    e = e+ap::sqr(v)/2;
                }
                
                //
                // From dError/dOut(normalized) to dError/dOut(non-normalized)
                //
                v = ap::vdotproduct(&network.nwbuf(nout), 1, &network.nwbuf(0), 1, ap::vlen(nout,2*nout-1));
                for(i = 0; i <= nout-1; i++)
                {
                    fown = network.nwbuf(i);
                    deown = network.nwbuf(nout+i);
                    network.chunks(iderror+ntotal-nout+i,k) = (-v+deown*fown+deown*(net-fown))*fown/ap::sqr(net);
                }
            }
        }
    }
    else
    {
        
        //
        // Normal output, regression network
        //
        // For each K = 0..CSize-1 do:
        // 1. calculate dError/dOut and place it to the second block of Chunks
        //
        for(i = 0; i <= nout-1; i++)
        {
            for(j = 0; j <= csize-1; j++)
            {
                v = network.chunks(ntotal-nout+i,j)*network.columnsigmas(nin+i)+network.columnmeans(nin+i)-xy(cstart+j,nin+i);
                network.chunks(iderror+ntotal-nout+i,j) = v*network.columnsigmas(nin+i);
                e = e+ap::sqr(v)/2;
            }
        }
    }
    
    //
    // Backpropagation
    //
    for(i = ntotal-1; i >= 0; i--)
    {
        
        //
        // Extract info
        //
        offs = istart+i*nfieldwidth;
        if( network.structinfo(offs+0)>0 )
        {
            
            //
            // Activation function
            //
            n1 = network.structinfo(offs+2);
            for(k = 0; k <= csize-1; k++)
            {
                network.chunks(iderror+i,k) = network.chunks(iderror+i,k)*network.chunks(idfdnet+i,k);
            }
            ap::vadd(&network.chunks(iderror+n1, 0), 1, &network.chunks(iderror+i, 0), 1, ap::vlen(0,csize-1));
        }
        if( network.structinfo(offs+0)==0 )
        {
            
            //
            // "Normal" activation function
            //
            n1 = network.structinfo(offs+2);
            n2 = n1+network.structinfo(offs+1)-1;
            w1 = network.structinfo(offs+3);
            w2 = w1+network.structinfo(offs+1)-1;
            for(j = w1; j <= w2; j++)
            {
                v = ap::vdotproduct(&network.chunks(n1+j-w1, 0), 1, &network.chunks(iderror+i, 0), 1, ap::vlen(0,csize-1));
                grad(j) = grad(j)+v;
            }
            for(j = n1; j <= n2; j++)
            {
                v = network.weights(w1+j-n1);
                ap::vadd(&network.chunks(iderror+j, 0), 1, &network.chunks(iderror+i, 0), 1, ap::vlen(0,csize-1), v);
            }
        }
        if( network.structinfo(offs+0)<0 )
        {
            bflag = false;
            if( network.structinfo(offs+0)==-2||network.structinfo(offs+0)==-3||network.structinfo(offs+0)==-4 )
            {
                
                //
                // Special neuron type, no back-propagation required
                //
                bflag = true;
            }
            ap::ap_error::make_assertion(bflag, "MLPInternalCalculateGradient: unknown neuron type!");
        }
    }
}


/*************************************************************************
Returns T*Ln(T/Z), guarded against overflow/underflow.
Internal subroutine.
*************************************************************************/
static double safecrossentropy(double t, double z)
{
    double result;
    double r;

    if( ap::fp_eq(t,0) )
    {
        result = 0;
    }
    else
    {
        if( ap::fp_greater(fabs(z),1) )
        {
            
            //
            // Shouldn't be the case with softmax,
            // but we just want to be sure.
            //
            if( ap::fp_eq(t/z,0) )
            {
                r = ap::minrealnumber;
            }
            else
            {
                r = t/z;
            }
        }
        else
        {
            
            //
            // Normal case
            //
            if( ap::fp_eq(z,0)||ap::fp_greater_eq(fabs(t),ap::maxrealnumber*fabs(z)) )
            {
                r = ap::maxrealnumber;
            }
            else
            {
                r = t/z;
            }
        }
        result = t*log(r);
    }
    return result;
}




