/*************************************************************************
Copyright (c) 2005-2007, Sergey Bochkanov (ALGLIB project).

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

#ifndef _sdet_h
#define _sdet_h

#include "ap.h"
#include "ialglib.h"

#include "ldlt.h"


/*************************************************************************
Determinant calculation of the matrix given by LDLT decomposition.

Input parameters:
    A       -   LDLT-decomposition of the matrix,
                output of subroutine SMatrixLDLT.
    Pivots  -   table of permutations which were made during
                LDLT decomposition, output of subroutine SMatrixLDLT.
    N       -   size of matrix A.
    IsUpper -   matrix storage format. The value is equal to the input
                parameter of subroutine SMatrixLDLT.

Result:
    matrix determinant.

  -- ALGLIB --
     Copyright 2005-2008 by Bochkanov Sergey
*************************************************************************/
double smatrixldltdet(const ap::real_2d_array& a,
     const ap::integer_1d_array& pivots,
     int n,
     bool isupper);


/*************************************************************************
Determinant calculation of the symmetric matrix

Input parameters:
    A       -   matrix. Array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   if IsUpper = True, then symmetric matrix A is given by its
                upper triangle, and the lower triangle isn’t used by
                subroutine. Similarly, if IsUpper = False, then A is given
                by its lower triangle.

Result:
    determinant of matrix A.

  -- ALGLIB --
     Copyright 2005-2008 by Bochkanov Sergey
*************************************************************************/
double smatrixdet(ap::real_2d_array a, int n, bool isupper);


double determinantldlt(const ap::real_2d_array& a,
     const ap::integer_1d_array& pivots,
     int n,
     bool isupper);


double determinantsymmetric(ap::real_2d_array a, int n, bool isupper);


#endif

