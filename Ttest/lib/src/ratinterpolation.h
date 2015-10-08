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

#ifndef _ratinterpolation_h
#define _ratinterpolation_h

#include "ap.h"
#include "ialglib.h"

#include "tsort.h"


/*************************************************************************
Rational barycentric interpolation without poles

The subroutine constructs the rational interpolating function without real
poles. It should be noted that the barycentric weights of the  interpolant
constructed are independent of the values of the given function.

Input parameters:
    X   -   interpolation nodes, array[0..N-1].
    N   -   number of nodes, N>0.
    D   -   order of the interpolation scheme, 0 <= D <= N-1.

Output parameters:
    W   -   array of the barycentric weights which  can  be  used  in  the
            BarycentricInterpolate subroutine. Array[0..N-1]

Note:
    this algorithm always succeeds and calculates the weights  with  close
    to machine precision.

  -- ALGLIB PROJECT --
     Copyright 17.06.2007 by Bochkanov Sergey
*************************************************************************/
void buildfloaterhormannrationalinterpolant(ap::real_1d_array x,
     int n,
     int d,
     ap::real_1d_array& w);


#endif

