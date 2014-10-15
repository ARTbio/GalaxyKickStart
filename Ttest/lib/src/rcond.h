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

#ifndef _rcond_h
#define _rcond_h

#include "ap.h"
#include "ialglib.h"

#include "reflections.h"
#include "creflections.h"
#include "hqrnd.h"
#include "matgen.h"
#include "ablasf.h"
#include "ablas.h"
#include "trfac.h"
#include "trlinsolve.h"
#include "safesolve.h"


/*************************************************************************
Estimate of a matrix condition number (1-norm)

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N   -   size of matrix A.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************/
double rmatrixrcond1(ap::real_2d_array a, int n);


/*************************************************************************
Estimate of a matrix condition number (infinity-norm).

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N   -   size of matrix A.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************/
double rmatrixrcondinf(ap::real_2d_array a, int n);


/*************************************************************************
Condition number estimate of a symmetric positive definite matrix.

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

It should be noted that 1-norm and inf-norm of condition numbers of symmetric
matrices are equal, so the algorithm doesn't take into account the
differences between these types of norms.

Input parameters:
    A       -   symmetric positive definite matrix which is given by its
                upper or lower triangle depending on the value of
                IsUpper. Array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format.

Result:
    1/LowerBound(cond(A)), if matrix A is positive definite,
   -1, if matrix A is not positive definite, and its condition number
    could not be found by this algorithm.

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************/
double spdmatrixrcond(ap::real_2d_array a, int n, bool isupper);


/*************************************************************************
Triangular matrix: estimate of a condition number (1-norm)

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    A       -   matrix. Array[0..N-1, 0..N-1].
    N       -   size of A.
    IsUpper -   True, if the matrix is upper triangular.
    IsUnit  -   True, if the matrix has a unit diagonal.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************/
double rmatrixtrrcond1(const ap::real_2d_array& a,
     int n,
     bool isupper,
     bool isunit);


/*************************************************************************
Triangular matrix: estimate of a matrix condition number (infinity-norm).

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N   -   size of matrix A.
    IsUpper -   True, if the matrix is upper triangular.
    IsUnit  -   True, if the matrix has a unit diagonal.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************/
double rmatrixtrrcondinf(const ap::real_2d_array& a,
     int n,
     bool isupper,
     bool isunit);


/*************************************************************************
Condition number estimate of a Hermitian positive definite matrix.

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

It should be noted that 1-norm and inf-norm of condition numbers of symmetric
matrices are equal, so the algorithm doesn't take into account the
differences between these types of norms.

Input parameters:
    A       -   Hermitian positive definite matrix which is given by its
                upper or lower triangle depending on the value of
                IsUpper. Array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format.

Result:
    1/LowerBound(cond(A)), if matrix A is positive definite,
   -1, if matrix A is not positive definite, and its condition number
    could not be found by this algorithm.

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************/
double hpdmatrixrcond(ap::complex_2d_array a, int n, bool isupper);


/*************************************************************************
Estimate of a matrix condition number (1-norm)

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N   -   size of matrix A.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************/
double cmatrixrcond1(ap::complex_2d_array a, int n);


/*************************************************************************
Estimate of a matrix condition number (infinity-norm).

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N   -   size of matrix A.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************/
double cmatrixrcondinf(ap::complex_2d_array a, int n);


/*************************************************************************
Estimate of the condition number of a matrix given by its LU decomposition (1-norm)

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    LUA         -   LU decomposition of a matrix in compact form. Output of
                    the RMatrixLU subroutine.
    N           -   size of matrix A.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************/
double rmatrixlurcond1(const ap::real_2d_array& lua, int n);


/*************************************************************************
Estimate of the condition number of a matrix given by its LU decomposition
(infinity norm).

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    LUA     -   LU decomposition of a matrix in compact form. Output of
                the RMatrixLU subroutine.
    N       -   size of matrix A.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************/
double rmatrixlurcondinf(const ap::real_2d_array& lua, int n);


/*************************************************************************
Condition number estimate of a symmetric positive definite matrix given by
Cholesky decomposition.

The algorithm calculates a lower bound of the condition number. In this
case, the algorithm does not return a lower bound of the condition number,
but an inverse number (to avoid an overflow in case of a singular matrix).

It should be noted that 1-norm and inf-norm condition numbers of symmetric
matrices are equal, so the algorithm doesn't take into account the
differences between these types of norms.

Input parameters:
    CD  - Cholesky decomposition of matrix A,
          output of SMatrixCholesky subroutine.
    N   - size of matrix A.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************/
double spdmatrixcholeskyrcond(const ap::real_2d_array& a, int n, bool isupper);


/*************************************************************************
Condition number estimate of a Hermitian positive definite matrix given by
Cholesky decomposition.

The algorithm calculates a lower bound of the condition number. In this
case, the algorithm does not return a lower bound of the condition number,
but an inverse number (to avoid an overflow in case of a singular matrix).

It should be noted that 1-norm and inf-norm condition numbers of symmetric
matrices are equal, so the algorithm doesn't take into account the
differences between these types of norms.

Input parameters:
    CD  - Cholesky decomposition of matrix A,
          output of SMatrixCholesky subroutine.
    N   - size of matrix A.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************/
double hpdmatrixcholeskyrcond(const ap::complex_2d_array& a,
     int n,
     bool isupper);


/*************************************************************************
Estimate of the condition number of a matrix given by its LU decomposition (1-norm)

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    LUA         -   LU decomposition of a matrix in compact form. Output of
                    the CMatrixLU subroutine.
    N           -   size of matrix A.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************/
double cmatrixlurcond1(const ap::complex_2d_array& lua, int n);


/*************************************************************************
Estimate of the condition number of a matrix given by its LU decomposition
(infinity norm).

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    LUA     -   LU decomposition of a matrix in compact form. Output of
                the CMatrixLU subroutine.
    N       -   size of matrix A.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************/
double cmatrixlurcondinf(const ap::complex_2d_array& lua, int n);


/*************************************************************************
Triangular matrix: estimate of a condition number (1-norm)

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    A       -   matrix. Array[0..N-1, 0..N-1].
    N       -   size of A.
    IsUpper -   True, if the matrix is upper triangular.
    IsUnit  -   True, if the matrix has a unit diagonal.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************/
double cmatrixtrrcond1(const ap::complex_2d_array& a,
     int n,
     bool isupper,
     bool isunit);


/*************************************************************************
Triangular matrix: estimate of a matrix condition number (infinity-norm).

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N   -   size of matrix A.
    IsUpper -   True, if the matrix is upper triangular.
    IsUnit  -   True, if the matrix has a unit diagonal.

Result: 1/LowerBound(cond(A))

NOTE:
    if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
    0.0 is returned in such cases.
*************************************************************************/
double cmatrixtrrcondinf(const ap::complex_2d_array& a,
     int n,
     bool isupper,
     bool isunit);


/*************************************************************************
Threshold for rcond: matrices with condition number beyond this  threshold
are considered singular.

Threshold must be far enough from underflow, at least Sqr(Threshold)  must
be greater than underflow.
*************************************************************************/
double rcondthreshold();


#endif

