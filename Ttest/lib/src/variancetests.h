/*************************************************************************
Copyright (c) 2007, Sergey Bochkanov (ALGLIB project).

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

#ifndef _variancetests_h
#define _variancetests_h

#include "ap.h"
#include "ialglib.h"

#include "gammafunc.h"
#include "normaldistr.h"
#include "ibetaf.h"
#include "fdistr.h"
#include "igammaf.h"
#include "chisquaredistr.h"


/*************************************************************************
Two-sample F-test

This test checks three hypotheses about dispersions of the given  samples.
The following tests are performed:
    * two-tailed test (null hypothesis - the dispersions are equal)
    * left-tailed test (null hypothesis  -  the  dispersion  of  the first
      sample is greater than or equal to  the  dispersion  of  the  second
      sample).
    * right-tailed test (null hypothesis - the  dispersion  of  the  first
      sample is less than or equal to the dispersion of the second sample)

The test is based on the following assumptions:
    * the given samples have normal distributions
    * the samples are independent.

Input parameters:
    X   -   sample 1. Array whose index goes from 0 to N-1.
    N   -   sample size.
    Y   -   sample 2. Array whose index goes from 0 to M-1.
    M   -   sample size.

Output parameters:
    BothTails   -   p-value for two-tailed test.
                    If BothTails is less than the given significance level
                    the null hypothesis is rejected.
    LeftTail    -   p-value for left-tailed test.
                    If LeftTail is less than the given significance level,
                    the null hypothesis is rejected.
    RightTail   -   p-value for right-tailed test.
                    If RightTail is less than the given significance level
                    the null hypothesis is rejected.

  -- ALGLIB --
     Copyright 19.09.2006 by Bochkanov Sergey
*************************************************************************/
void ftest(const ap::real_1d_array& x,
     int n,
     const ap::real_1d_array& y,
     int m,
     double& bothtails,
     double& lefttail,
     double& righttail);


/*************************************************************************
One-sample chi-square test

This test checks three hypotheses about the dispersion of the given sample
The following tests are performed:
    * two-tailed test (null hypothesis - the dispersion equals  the  given
      number)
    * left-tailed test (null hypothesis - the dispersion is  greater  than
      or equal to the given number)
    * right-tailed test (null hypothesis  -  dispersion is  less  than  or
      equal to the given number).

Test is based on the following assumptions:
    * the given sample has a normal distribution.

Input parameters:
    X           -   sample 1. Array whose index goes from 0 to N-1.
    N           -   size of the sample.
    Variance    -   dispersion value to compare with.

Output parameters:
    BothTails   -   p-value for two-tailed test.
                    If BothTails is less than the given significance level
                    the null hypothesis is rejected.
    LeftTail    -   p-value for left-tailed test.
                    If LeftTail is less than the given significance level,
                    the null hypothesis is rejected.
    RightTail   -   p-value for right-tailed test.
                    If RightTail is less than the given significance level
                    the null hypothesis is rejected.

  -- ALGLIB --
     Copyright 19.09.2006 by Bochkanov Sergey
*************************************************************************/
void onesamplevariancetest(const ap::real_1d_array& x,
     int n,
     double variance,
     double& bothtails,
     double& lefttail,
     double& righttail);


#endif

