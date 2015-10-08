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

#ifndef _ssolve_h
#define _ssolve_h

#include "ap.h"
#include "ialglib.h"

#include "ldlt.h"


/*************************************************************************
Solving  a system  of linear equations  with a system matrix  given by its
LDLT decomposition

The algorithm solves systems with a square matrix only.

Input parameters:
    A       -   LDLT decomposition of the matrix (the result of the
                SMatrixLDLT subroutine).
    Pivots  -   row permutation table (the result of the SMatrixLDLT subroutine).
    B       -   right side of a system.
                Array whose index ranges within [0..N-1].
    N       -   size of matrix A.
    IsUpper -   points to the triangle of matrix A in which the LDLT
                decomposition is stored.
                If IsUpper=True, the decomposition has the form of U*D*U',
                matrix U is stored in the upper triangle of  matrix A  (in
                that case, the lower triangle isn't used and isn't changed
                by the subroutine).
                Similarly, if IsUpper=False, the decomposition has the form
                of L*D*L' and the lower triangle stores matrix L.

Output parameters:
    X       -   solution of a system.
                Array whose index ranges within [0..N-1].

Result:
    True, if the matrix is not singular. X contains the solution.
    False, if the matrix is singular (the determinant of matrix D is equal
to 0). In this case, X doesn't contain a solution.
*************************************************************************/
bool smatrixldltsolve(const ap::real_2d_array& a,
     const ap::integer_1d_array& pivots,
     ap::real_1d_array b,
     int n,
     bool isupper,
     ap::real_1d_array& x);


/*************************************************************************
Solving a system of linear equations with a symmetric system matrix

Input parameters:
    A       -   system matrix (upper or lower triangle).
                Array whose indexes range within [0..N-1, 0..N-1].
    B       -   right side of a system.
                Array whose index ranges within [0..N-1].
    N       -   size of matrix A.
    IsUpper -   If IsUpper = True, A contains the upper triangle,
                otherwise A contains the lower triangle.

Output parameters:
    X       -   solution of a system.
                Array whose index ranges within [0..N-1].

Result:
    True, if the matrix is not singular. X contains the solution.
    False, if the matrix is singular (the determinant of the matrix is equal
to 0). In this case, X doesn't contain a solution.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************/
bool smatrixsolve(ap::real_2d_array a,
     const ap::real_1d_array& b,
     int n,
     bool isupper,
     ap::real_1d_array& x);


bool solvesystemldlt(const ap::real_2d_array& a,
     const ap::integer_1d_array& pivots,
     ap::real_1d_array b,
     int n,
     bool isupper,
     ap::real_1d_array& x);


bool solvesymmetricsystem(ap::real_2d_array a,
     ap::real_1d_array b,
     int n,
     bool isupper,
     ap::real_1d_array& x);


#endif

