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

#ifndef _hsschur_h
#define _hsschur_h

#include "ap.h"
#include "ialglib.h"

#include "blas.h"
#include "reflections.h"
#include "rotations.h"


/*************************************************************************
Subroutine performing  the  Schur  decomposition  of  a  matrix  in  upper
Hessenberg form using the QR algorithm with multiple shifts.

The  source matrix  H  is  represented as  S'*H*S = T, where H - matrix in
upper Hessenberg form,  S - orthogonal matrix (Schur vectors),   T - upper
quasi-triangular matrix (with blocks of sizes  1x1  and  2x2  on  the main
diagonal).

Input parameters:
    H   -   matrix to be decomposed.
            Array whose indexes range within [1..N, 1..N].
    N   -   size of H, N>=0.


Output parameters:
    H   –   contains the matrix T.
            Array whose indexes range within [1..N, 1..N].
            All elements below the blocks on the main diagonal are equal
            to 0.
    S   -   contains Schur vectors.
            Array whose indexes range within [1..N, 1..N].

Note 1:
    The block structure of matrix T could be easily recognized: since  all
    the elements  below  the blocks are zeros, the elements a[i+1,i] which
    are equal to 0 show the block border.

Note 2:
    the algorithm  performance  depends  on  the  value  of  the  internal
    parameter NS of InternalSchurDecomposition  subroutine  which  defines
    the number of shifts in the QR algorithm (analog of  the  block  width
    in block matrix algorithms in linear algebra). If you require  maximum
    performance  on  your  machine,  it  is  recommended  to  adjust  this
    parameter manually.

Result:
    True, if the algorithm has converged and the parameters H and S contain
        the result.
    False, if the algorithm has not converged.

Algorithm implemented on the basis of subroutine DHSEQR (LAPACK 3.0 library).
*************************************************************************/
bool upperhessenbergschurdecomposition(ap::real_2d_array& h,
     int n,
     ap::real_2d_array& s);


void internalschurdecomposition(ap::real_2d_array& h,
     int n,
     int tneeded,
     int zneeded,
     ap::real_1d_array& wr,
     ap::real_1d_array& wi,
     ap::real_2d_array& z,
     int& info);


#endif

