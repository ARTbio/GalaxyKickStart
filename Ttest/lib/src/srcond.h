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

#ifndef _srcond_h
#define _srcond_h

#include "ap.h"
#include "ialglib.h"

#include "ldlt.h"
#include "ssolve.h"
#include "estnorm.h"


/*************************************************************************
Condition number estimate of a symmetric matrix

The algorithm calculates a lower bound of the condition number. In this
case, the algorithm does not return a lower bound of the condition number,
but an inverse number (to avoid an overflow in case of a singular matrix).

It should be noted that 1-norm and inf-norm condition numbers of symmetric
matrices are equal, so the algorithm doesn't take into account the
differences between these types of norms.

Input parameters:
    A       -   symmetric definite matrix which is given by its upper or
                lower triangle depending on IsUpper.
                Array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format.

Result:
    1/LowerBound(cond(A))
*************************************************************************/
double smatrixrcond(const ap::real_2d_array& a, int n, bool isupper);


/*************************************************************************
Condition number estimate of a matrix given by LDLT-decomposition

The algorithm calculates a lower bound of the condition number. In this
case, the algorithm does not return a lower bound of the condition number,
but an inverse number (to avoid an overflow in case of a singular matrix).

It should be noted that 1-norm and inf-norm condition numbers of symmetric
matrices are equal, so the algorithm doesn't take into account the
differences between these types of norms.

Input parameters:
    L       -   LDLT-decomposition of matrix A given by the upper or lower
                triangle depending on IsUpper.
                Output of SMatrixLDLT subroutine.
    Pivots  -   table of permutations which were made during LDLT-decomposition,
                Output of SMatrixLDLT subroutine.
    N       -   size of matrix A.
    IsUpper -   storage format.

Result:
    1/LowerBound(cond(A))
*************************************************************************/
double smatrixldltrcond(const ap::real_2d_array& l,
     const ap::integer_1d_array& pivots,
     int n,
     bool isupper);


double rcondsymmetric(ap::real_2d_array a, int n, bool isupper);


double rcondldlt(const ap::real_2d_array& l,
     const ap::integer_1d_array& pivots,
     int n,
     bool isupper);


void internalldltrcond(const ap::real_2d_array& l,
     const ap::integer_1d_array& pivots,
     int n,
     bool isupper,
     bool isnormprovided,
     double anorm,
     double& rcond);


#endif

