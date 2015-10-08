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

#ifndef _sinverse_h
#define _sinverse_h

#include "ap.h"
#include "ialglib.h"

#include "sblas.h"
#include "ldlt.h"


/*************************************************************************
Inversion of a symmetric indefinite matrix

The algorithm gets an LDLT-decomposition as an input, generates matrix A^-1
and saves the lower or upper triangle of an inverse matrix depending on the
input (U*D*U' or L*D*L').

Input parameters:
    A       -   LDLT-decomposition of the matrix,
                Output of subroutine SMatrixLDLT.
    N       -   size of matrix A.
    IsUpper -   storage format. If IsUpper = True, then the symmetric matrix
                is given as decomposition A = U*D*U' and this decomposition
                is stored in the upper triangle of matrix A and on the main
                diagonal, and the lower triangle of matrix A is not used.
    Pivots  -   a table of permutations, output of subroutine SMatrixLDLT.

Output parameters:
    A       -   inverse of the matrix, whose LDLT-decomposition was stored
                in matrix A as a subroutine input.
                Array with elements [0..N-1, 0..N-1].
                If IsUpper = True, then A contains the upper triangle of
                matrix A^-1, and the elements below the main diagonal are
                not used nor changed. The same applies if IsUpper = False.

Result:
    True, if the matrix is not singular.
    False, if the matrix is singular and could not be inverted.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     March 31, 1993
*************************************************************************/
bool smatrixldltinverse(ap::real_2d_array& a,
     const ap::integer_1d_array& pivots,
     int n,
     bool isupper);


/*************************************************************************
Inversion of a symmetric indefinite matrix

Given a lower or upper triangle of matrix A, the algorithm generates
matrix A^-1 and saves the lower or upper triangle depending on the input.

Input parameters:
    A       -   matrix to be inverted (upper or lower triangle).
                Array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format. If IsUpper = True, then the upper
                triangle of matrix A is given, otherwise the lower
                triangle is given.

Output parameters:
    A       -   inverse of matrix A.
                Array with elements [0..N-1, 0..N-1].
                If IsUpper = True, then A contains the upper triangle of
                matrix A^-1, and the elements below the main diagonal are
                not used nor changed.
                The same applies if IsUpper = False.

Result:
    True, if the matrix is not singular.
    False, if the matrix is singular and could not be inverted.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     March 31, 1993
*************************************************************************/
bool smatrixinverse(ap::real_2d_array& a, int n, bool isupper);


bool inverseldlt(ap::real_2d_array& a,
     const ap::integer_1d_array& pivots,
     int n,
     bool isupper);


bool inversesymmetricindefinite(ap::real_2d_array& a, int n, bool isupper);


#endif

