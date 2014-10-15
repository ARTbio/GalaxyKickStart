/*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from FORTRAN to
      pseudocode.

See subroutines comments for additional copyrights.

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

#ifndef _estnorm_h
#define _estnorm_h

#include "ap.h"
#include "ialglib.h"

/*************************************************************************
Matrix norm estimation

The algorithm estimates the 1-norm of square matrix A  on  the  assumption
that the multiplication of matrix  A  by  the  vector  is  available  (the
iterative method is used). It is recommended to use this algorithm  if  it
is hard  to  calculate  matrix  elements  explicitly  (for  example,  when
estimating the inverse matrix norm).

The algorithm uses back communication for multiplying the  vector  by  the
matrix.  If  KASE=0  after  returning from a subroutine, its execution was
completed successfully, otherwise it is required to multiply the  returned
vector by matrix A and call the subroutine again.

The DemoIterativeEstimateNorm subroutine shows a simple example.

Parameters:
    N       -   size of matrix A.
    V       -   vector.   It is initialized by the subroutine on the first
                call. It is then passed into it on repeated calls.
    X       -   if KASE<>0, it contains the vector to be replaced by:
                    A * X,      if KASE=1
                    A^T * X,    if KASE=2
                Array whose index ranges within [1..N].
    ISGN    -   vector. It is initialized by the subroutine on  the  first
                call. It is then passed into it on repeated calls.
    EST     -   if KASE=0, it contains the lower boundary of the matrix
                norm estimate.
    KASE    -   on the first call, it should be equal to 0. After the last
                return, it is equal to 0 (EST contains the  matrix  norm),
                on intermediate returns it can be equal to 1 or 2 depending
                on the operation to be performed on vector X.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************/
void iterativeestimate1norm(int n,
     ap::real_1d_array& v,
     ap::real_1d_array& x,
     ap::integer_1d_array& isgn,
     double& est,
     int& kase);


/*************************************************************************
Example of usage of an IterativeEstimateNorm subroutine

Input parameters:
    A   -   matrix.
            Array whose indexes range within [1..N, 1..N].

Return:
    Matrix norm estimated by the subroutine.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************/
double demoiterativeestimate1norm(const ap::real_2d_array& a, int n);


#endif

