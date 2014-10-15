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

#ifndef _jarquebera_h
#define _jarquebera_h

#include "ap.h"
#include "ialglib.h"

/*************************************************************************
Jarque-Bera test

This test checks hypotheses about the fact that a  given  sample  X  is  a
sample of normal random variable.

Requirements:
    * the number of elements in the sample is not less than 5.

Input parameters:
    X   -   sample. Array whose index goes from 0 to N-1.
    N   -   size of the sample. N>=5

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

Accuracy of the approximation used (5<=N<=1951):

p-value  	    relative error (5<=N<=1951)
[1, 0.1]            < 1%
[0.1, 0.01]         < 2%
[0.01, 0.001]       < 6%
[0.001, 0]          wasn't measured

For N>1951 accuracy wasn't measured but it shouldn't be sharply  different
from table values.

  -- ALGLIB --
     Copyright 09.04.2007 by Bochkanov Sergey
*************************************************************************/
void jarqueberatest(const ap::real_1d_array& x, int n, double& p);


#endif

