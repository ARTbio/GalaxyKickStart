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

#ifndef _mlpbase_h
#define _mlpbase_h

#include "ap.h"
#include "ialglib.h"

struct multilayerperceptron
{
    ap::integer_1d_array structinfo;
    ap::real_1d_array weights;
    ap::real_1d_array columnmeans;
    ap::real_1d_array columnsigmas;
    ap::real_1d_array neurons;
    ap::real_1d_array dfdnet;
    ap::real_1d_array derror;
    ap::real_1d_array x;
    ap::real_1d_array y;
    ap::real_2d_array chunks;
    ap::real_1d_array nwbuf;
};




/*************************************************************************
Creates  neural  network  with  NIn  inputs,  NOut outputs, without hidden
layers, with linear output layer. Network weights are  filled  with  small
random values.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreate0(int nin, int nout, multilayerperceptron& network);


/*************************************************************************
Same  as  MLPCreate0,  but  with  one  hidden  layer  (NHid  neurons) with
non-linear activation function. Output layer is linear.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreate1(int nin, int nhid, int nout, multilayerperceptron& network);


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
     multilayerperceptron& network);


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
     multilayerperceptron& network);


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
     multilayerperceptron& network);


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
     multilayerperceptron& network);


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
     multilayerperceptron& network);


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
     multilayerperceptron& network);


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
     multilayerperceptron& network);


/*************************************************************************
Creates classifier network with NIn  inputs  and  NOut  possible  classes.
Network contains no hidden layers and linear output  layer  with  SOFTMAX-
normalization  (so  outputs  sums  up  to  1.0  and  converge to posterior
probabilities).

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreatec0(int nin, int nout, multilayerperceptron& network);


/*************************************************************************
Same as MLPCreateC0, but with one non-linear hidden layer.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreatec1(int nin, int nhid, int nout, multilayerperceptron& network);


/*************************************************************************
Same as MLPCreateC0, but with two non-linear hidden layers.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreatec2(int nin,
     int nhid1,
     int nhid2,
     int nout,
     multilayerperceptron& network);


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
     multilayerperceptron& network2);


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
     int& rlen);


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
     multilayerperceptron& network);


/*************************************************************************
Randomization of neural network weights

  -- ALGLIB --
     Copyright 06.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlprandomize(multilayerperceptron& network);


/*************************************************************************
Randomization of neural network weights and standartisator

  -- ALGLIB --
     Copyright 10.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlprandomizefull(multilayerperceptron& network);


/*************************************************************************
Internal subroutine.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpinitpreprocessor(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int ssize);


/*************************************************************************
Returns information about initialized network: number of inputs, outputs,
weights.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpproperties(const multilayerperceptron& network,
     int& nin,
     int& nout,
     int& wcount);


/*************************************************************************
Tells whether network is SOFTMAX-normalized (i.e. classifier) or not.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
bool mlpissoftmax(const multilayerperceptron& network);


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
     ap::real_1d_array& y);


/*************************************************************************
Error function for neural network, internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
double mlperror(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int ssize);


/*************************************************************************
Natural error function for neural network, internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
double mlperrorn(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int ssize);


/*************************************************************************
Classification error

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
int mlpclserror(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int ssize);


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
     int npoints);


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
     int npoints);


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
     int npoints);


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
     int npoints);


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
     int npoints);


/*************************************************************************
Gradient calculation. Internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgrad(multilayerperceptron& network,
     const ap::real_1d_array& x,
     const ap::real_1d_array& desiredy,
     double& e,
     ap::real_1d_array& grad);


/*************************************************************************
Gradient calculation (natural error function). Internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgradn(multilayerperceptron& network,
     const ap::real_1d_array& x,
     const ap::real_1d_array& desiredy,
     double& e,
     ap::real_1d_array& grad);


/*************************************************************************
Batch gradient calculation. Internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgradbatch(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int ssize,
     double& e,
     ap::real_1d_array& grad);


/*************************************************************************
Batch gradient calculation (natural error function). Internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgradnbatch(multilayerperceptron& network,
     const ap::real_2d_array& xy,
     int ssize,
     double& e,
     ap::real_1d_array& grad);


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
     ap::real_2d_array& h);


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
     ap::real_2d_array& h);


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
     ap::real_1d_array& y);


#endif

