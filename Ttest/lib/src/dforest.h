/*************************************************************************
Copyright (c) 2009, Sergey Bochkanov (ALGLIB project).

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

#ifndef _dforest_h
#define _dforest_h

#include "ap.h"
#include "ialglib.h"

#include "tsort.h"
#include "descriptivestatistics.h"
#include "bdss.h"


struct decisionforest
{
    int nvars;
    int nclasses;
    int ntrees;
    int bufsize;
    ap::real_1d_array trees;
};


struct dfreport
{
    double relclserror;
    double avgce;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double oobrelclserror;
    double oobavgce;
    double oobrmserror;
    double oobavgerror;
    double oobavgrelerror;
};


struct dfinternalbuffers
{
    ap::real_1d_array treebuf;
    ap::integer_1d_array idxbuf;
    ap::real_1d_array tmpbufr;
    ap::real_1d_array tmpbufr2;
    ap::integer_1d_array tmpbufi;
    ap::integer_1d_array classibuf;
    ap::integer_1d_array varpool;
    ap::boolean_1d_array evsbin;
    ap::real_1d_array evssplits;
};




/*************************************************************************
This subroutine builds random decision forest.

INPUT PARAMETERS:
    XY          -   training set
    NPoints     -   training set size, NPoints>=1
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   task type:
                    * NClasses=1 - regression task with one
                                   dependent variable
                    * NClasses>1 - classification task with
                                   NClasses classes.
    NTrees      -   number of trees in a forest, NTrees>=1.
                    recommended values: 50-100.
    R           -   percent of a training set used to build
                    individual trees. 0<R<=1.
                    recommended values: 0.1 <= R <= 0.66.

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<1, NVars<1, NClasses<1, NTrees<1, R<=0
                          or R>1).
                    *  1, if task has been solved
    DF          -   model built
    Rep         -   training report, contains error on a training set
                    and out-of-bag estimates of generalization error.

  -- ALGLIB --
     Copyright 19.02.2009 by Bochkanov Sergey
*************************************************************************/
void dfbuildrandomdecisionforest(const ap::real_2d_array& xy,
     int npoints,
     int nvars,
     int nclasses,
     int ntrees,
     double r,
     int& info,
     decisionforest& df,
     dfreport& rep);


void dfbuildinternal(const ap::real_2d_array& xy,
     int npoints,
     int nvars,
     int nclasses,
     int ntrees,
     int samplesize,
     int nfeatures,
     int flags,
     int& info,
     decisionforest& df,
     dfreport& rep);


/*************************************************************************
Procesing

INPUT PARAMETERS:
    DF      -   decision forest model
    X       -   input vector,  array[0..NVars-1].

OUTPUT PARAMETERS:
    Y       -   result. Regression estimate when solving regression  task,
                vector of posterior probabilities for classification task.
                Subroutine does not allocate memory for this vector, it is
                responsibility of a caller to allocate it. Array  must  be
                at least [0..NClasses-1].

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
void dfprocess(const decisionforest& df,
     const ap::real_1d_array& x,
     ap::real_1d_array& y);


/*************************************************************************
Relative classification error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    percent of incorrectly classified cases.
    Zero if model solves regression task.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfrelclserror(const decisionforest& df,
     const ap::real_2d_array& xy,
     int npoints);


/*************************************************************************
Average cross-entropy (in bits per element) on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    CrossEntropy/(NPoints*LN(2)).
    Zero if model solves regression task.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfavgce(const decisionforest& df,
     const ap::real_2d_array& xy,
     int npoints);


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.
    Its meaning for regression task is obvious. As for
    classification task, RMS error means error when estimating posterior
    probabilities.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfrmserror(const decisionforest& df,
     const ap::real_2d_array& xy,
     int npoints);


/*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for
    classification task, it means average error when estimating posterior
    probabilities.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfavgerror(const decisionforest& df,
     const ap::real_2d_array& xy,
     int npoints);


/*************************************************************************
Average relative error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for
    classification task, it means average relative error when estimating
    posterior probability of belonging to the correct class.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfavgrelerror(const decisionforest& df,
     const ap::real_2d_array& xy,
     int npoints);


/*************************************************************************
Copying of DecisionForest strucure

INPUT PARAMETERS:
    DF1 -   original

OUTPUT PARAMETERS:
    DF2 -   copy

  -- ALGLIB --
     Copyright 13.02.2009 by Bochkanov Sergey
*************************************************************************/
void dfcopy(const decisionforest& df1, decisionforest& df2);


/*************************************************************************
Serialization of DecisionForest strucure

INPUT PARAMETERS:
    DF      -   original

OUTPUT PARAMETERS:
    RA      -   array of real numbers which stores decision forest,
                array[0..RLen-1]
    RLen    -   RA lenght

  -- ALGLIB --
     Copyright 13.02.2009 by Bochkanov Sergey
*************************************************************************/
void dfserialize(const decisionforest& df, ap::real_1d_array& ra, int& rlen);


/*************************************************************************
Unserialization of DecisionForest strucure

INPUT PARAMETERS:
    RA      -   real array which stores decision forest

OUTPUT PARAMETERS:
    DF      -   restored structure

  -- ALGLIB --
     Copyright 13.02.2009 by Bochkanov Sergey
*************************************************************************/
void dfunserialize(const ap::real_1d_array& ra, decisionforest& df);


#endif

