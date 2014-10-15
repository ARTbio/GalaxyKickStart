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

#ifndef _logit_h
#define _logit_h

#include "ap.h"
#include "ialglib.h"

#include "descriptivestatistics.h"
#include "mlpbase.h"
#include "hblas.h"
#include "reflections.h"
#include "creflections.h"
#include "sblas.h"
#include "ablasf.h"
#include "ablas.h"
#include "ortfac.h"
#include "blas.h"
#include "rotations.h"
#include "bdsvd.h"
#include "svd.h"
#include "hqrnd.h"
#include "matgen.h"
#include "trfac.h"
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "xblas.h"
#include "densesolver.h"
#include "tsort.h"
#include "bdss.h"


struct logitmodel
{
    ap::real_1d_array w;
};


struct logitmcstate
{
    bool brackt;
    bool stage1;
    int infoc;
    double dg;
    double dgm;
    double dginit;
    double dgtest;
    double dgx;
    double dgxm;
    double dgy;
    double dgym;
    double finit;
    double ftest1;
    double fm;
    double fx;
    double fxm;
    double fy;
    double fym;
    double stx;
    double sty;
    double stmin;
    double stmax;
    double width;
    double width1;
    double xtrapf;
};


/*************************************************************************
MNLReport structure contains information about training process:
* NGrad     -   number of gradient calculations
* NHess     -   number of Hessian calculations
*************************************************************************/
struct mnlreport
{
    int ngrad;
    int nhess;
};




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
     mnlreport& rep);


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
     ap::real_1d_array& y);


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
     int& nclasses);


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
     logitmodel& lm);


/*************************************************************************
Copying of LogitModel strucure

INPUT PARAMETERS:
    LM1 -   original

OUTPUT PARAMETERS:
    LM2 -   copy

  -- ALGLIB --
     Copyright 15.03.2009 by Bochkanov Sergey
*************************************************************************/
void mnlcopy(const logitmodel& lm1, logitmodel& lm2);


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
void mnlserialize(const logitmodel& lm, ap::real_1d_array& ra, int& rlen);


/*************************************************************************
Unserialization of LogitModel strucure

INPUT PARAMETERS:
    RA      -   real array which stores model

OUTPUT PARAMETERS:
    LM      -   restored model

  -- ALGLIB --
     Copyright 15.03.2009 by Bochkanov Sergey
*************************************************************************/
void mnlunserialize(const ap::real_1d_array& ra, logitmodel& lm);


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
double mnlavgce(logitmodel& lm, const ap::real_2d_array& xy, int npoints);


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
     int npoints);


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
double mnlrmserror(logitmodel& lm, const ap::real_2d_array& xy, int npoints);


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
double mnlavgerror(logitmodel& lm, const ap::real_2d_array& xy, int npoints);


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
double mnlavgrelerror(logitmodel& lm, const ap::real_2d_array& xy, int ssize);


/*************************************************************************
Classification error on test set = MNLRelClsError*NPoints

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
int mnlclserror(logitmodel& lm, const ap::real_2d_array& xy, int npoints);


#endif

