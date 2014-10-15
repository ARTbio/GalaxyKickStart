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

#ifndef _mannwhitneyu_h
#define _mannwhitneyu_h

#include "ap.h"
#include "ialglib.h"

/*************************************************************************
Mann-Whitney U-test

This test checks hypotheses about whether X  and  Y  are  samples  of  two
continuous distributions of the same shape  and  same  median  or  whether
their medians are different.

The following tests are performed:
    * two-tailed test (null hypothesis - the medians are equal)
    * left-tailed test (null hypothesis - the median of the  first  sample
      is greater than or equal to the median of the second sample)
    * right-tailed test (null hypothesis - the median of the first  sample
      is less than or equal to the median of the second sample).

Requirements:
    * the samples are independent
    * X and Y are continuous distributions (or discrete distributions well-
      approximating continuous distributions)
    * distributions of X and Y have the  same  shape.  The  only  possible
      difference is their position (i.e. the value of the median)
    * the number of elements in each sample is not less than 5
    * the scale of measurement should be ordinal, interval or ratio  (i.e.
      the test could not be applied to nominal variables).

The test is non-parametric and doesn't require distributions to be normal.

Input parameters:
    X   -   sample 1. Array whose index goes from 0 to N-1.
    N   -   size of the sample. N>=5
    Y   -   sample 2. Array whose index goes from 0 to M-1.
    M   -   size of the sample. M>=5

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

To calculate p-values, special approximation is used. This method lets  us
calculate p-values with satisfactory  accuracy  in  interval  [0.0001, 1].
There is no approximation outside the [0.0001, 1] interval. Therefore,  if
the significance level outlies this interval, the test returns 0.0001.

Relative precision of approximation of p-value:

N          M          Max.err.   Rms.err.
5..10      N..10      1.4e-02    6.0e-04
5..10      N..100     2.2e-02    5.3e-06
10..15     N..15      1.0e-02    3.2e-04
10..15     N..100     1.0e-02    2.2e-05
15..100    N..100     6.1e-03    2.7e-06

For N,M>100 accuracy checks weren't put into  practice,  but  taking  into
account characteristics of asymptotic approximation used, precision should
not be sharply different from the values for interval [5, 100].

  -- ALGLIB --
     Copyright 09.04.2007 by Bochkanov Sergey
*************************************************************************/
void mannwhitneyutest(const ap::real_1d_array& x,
     int n,
     const ap::real_1d_array& y,
     int m,
     double& bothtails,
     double& lefttail,
     double& righttail);


#endif

