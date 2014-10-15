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

#ifndef _ldlt_h
#define _ldlt_h

#include "ap.h"
#include "ialglib.h"

/*************************************************************************
LDLTDecomposition of a symmetric matrix

The algorithm represents a symmetric matrix (which is not necessarily
positive definite) as A=L*D*L' or A = U*D*U', where D is a block-diagonal
matrix with blocks 1x1 or 2x2, matrix L (matrix U) is a product of lower
(upper) triangular matrices with unit diagonal and permutation matrices.

Input parameters:
    A       -   factorized matrix, array with elements [0..N-1, 0..N-1].
                If IsUpper – True, then the upper triangle contains
                elements of symmetric matrix A, and the lower triangle is
                not used.
                The same applies if IsUpper = False.
    N       -   size of factorized matrix.
    IsUpper -   parameter which shows a method of matrix definition (lower
                or upper triangle).

Output parameters:
    A       -   matrices D and U, if IsUpper = True, or L, if IsUpper = False,
                in compact form, replacing the upper (lower) triangle of
                matrix A. In that case, the elements under (over) the main
                diagonal are not used nor modified.
    Pivots  -   tables of performed permutations (see below).

If IsUpper = True, then A = U*D*U', U = P(n)*U(n)*...*P(k)*U(k), where
P(k) is the permutation matrix, U(k) - upper triangular matrix with its
unit main diagonal and k decreases from n with step s which is equal to
1 or 2 (according to the size of the blocks of matrix D).

        (   I    v    0   )   k-s+1
U(k) =  (   0    I    0   )   s
        (   0    0    I   )   n-k-1
           k-s+1 s   n-k-1

If Pivots[k]>=0, then s=1, P(k) - permutation of rows k and Pivots[k], the
vectorv forming matrix U(k) is stored in elements A(0:k-1,k), D(k) replaces
A(k,k). If Pivots[k]=Pivots[k-1]<0 then s=2, P(k) - permutation of rows k-1
and N+Pivots[k-1], the vector v forming matrix U(k) is stored in elements
A(0:k-1,k:k+1), the upper triangle of block D(k) is stored in A(k,k),
A(k,k+1) and A(k+1,k+1).

If IsUpper = False, then A = L*D*L', L=P(0)*L(0)*...*P(k)*L(k), where P(k)
is the permutation matrix, L(k) – lower triangular matrix with unit main
diagonal and k decreases from 1 with step s which is equal to 1 or 2
(according to the size of the blocks of matrix D).

        (   I    0     0   )  k-1
L(k) =  (   0    I     0   )  s
        (   0    v     I   )  n-k-s+1
           k-1   s  n-k-s+1

If Pivots[k]>=0 then s=1, P(k) – permutation of rows k and Pivots[k], the
vector v forming matrix L(k) is stored in elements A(k+1:n-1,k), D(k)
replaces A(k,k). If Pivots[k]=Pivots[k+1]<0 then s=2, P(k) - permutation
of rows k+1 and N+Pivots[k+1], the vector v forming matrix L(k) is stored
in elements A(k+2:n-1,k:k+1), the lower triangle of block D(k) is stored in
A(k,k), A(k+1,k) and A(k+1,k+1).

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     June 30, 1999
*************************************************************************/
void smatrixldlt(ap::real_2d_array& a,
     int n,
     bool isupper,
     ap::integer_1d_array& pivots);


void ldltdecomposition(ap::real_2d_array& a,
     int n,
     bool isupper,
     ap::integer_1d_array& pivots);


#endif

