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

#include <stdafx.h>
#include "rcond.h"

static void rmatrixrcondtrinternal(const ap::real_2d_array& a,
     int n,
     bool isupper,
     bool isunit,
     bool onenorm,
     double anorm,
     double& rc);
static void cmatrixrcondtrinternal(const ap::complex_2d_array& a,
     const int& n,
     bool isupper,
     bool isunit,
     bool onenorm,
     double anorm,
     double& rc);
static void spdmatrixrcondcholeskyinternal(const ap::real_2d_array& cha,
     int n,
     bool isupper,
     bool isnormprovided,
     double anorm,
     double& rc);
static void hpdmatrixrcondcholeskyinternal(const ap::complex_2d_array& cha,
     int n,
     bool isupper,
     bool isnormprovided,
     double anorm,
     double& rc);
static void rmatrixrcondluinternal(const ap::real_2d_array& lua,
     int n,
     bool onenorm,
     bool isanormprovided,
     double anorm,
     double& rc);
static void cmatrixrcondluinternal(const ap::complex_2d_array& lua,
     const int& n,
     bool onenorm,
     bool isanormprovided,
     double anorm,
     double& rc);
static void rmatrixestimatenorm(int n,
     ap::real_1d_array& v,
     ap::real_1d_array& x,
     ap::integer_1d_array& isgn,
     double& est,
     int& kase);
static void cmatrixestimatenorm(const int& n,
     ap::complex_1d_array& v,
     ap::complex_1d_array& x,
     double& est,
     int& kase,
     ap::integer_1d_array& isave,
     ap::real_1d_array& rsave);
static double internalcomplexrcondscsum1(const ap::complex_1d_array& x, int n);
static int internalcomplexrcondicmax1(const ap::complex_1d_array& x, int n);
static void internalcomplexrcondsaveall(ap::integer_1d_array& isave,
     ap::real_1d_array& rsave,
     int& i,
     int& iter,
     int& j,
     int& jlast,
     int& jump,
     double& absxi,
     double& altsgn,
     double& estold,
     double& temp);
static void internalcomplexrcondloadall(ap::integer_1d_array& isave,
     ap::real_1d_array& rsave,
     int& i,
     int& iter,
     int& j,
     int& jlast,
     int& jump,
     double& absxi,
     double& altsgn,
     double& estold,
     double& temp);

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
double rmatrixrcond1(ap::real_2d_array a, int n)
{
    double result;
    int i;
    int j;
    double v;
    double nrm;
    ap::integer_1d_array pivots;
    ap::real_1d_array t;

    ap::ap_error::make_assertion(n>=1, "RMatrixRCond1: N<1!");
    t.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        t(i) = 0;
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            t(j) = t(j)+fabs(a(i,j));
        }
    }
    nrm = 0;
    for(i = 0; i <= n-1; i++)
    {
        nrm = ap::maxreal(nrm, t(i));
    }
    rmatrixlu(a, n, n, pivots);
    rmatrixrcondluinternal(a, n, true, true, nrm, v);
    result = v;
    return result;
}


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
double rmatrixrcondinf(ap::real_2d_array a, int n)
{
    double result;
    int i;
    int j;
    double v;
    double nrm;
    ap::integer_1d_array pivots;

    ap::ap_error::make_assertion(n>=1, "RMatrixRCondInf: N<1!");
    nrm = 0;
    for(i = 0; i <= n-1; i++)
    {
        v = 0;
        for(j = 0; j <= n-1; j++)
        {
            v = v+fabs(a(i,j));
        }
        nrm = ap::maxreal(nrm, v);
    }
    rmatrixlu(a, n, n, pivots);
    rmatrixrcondluinternal(a, n, false, true, nrm, v);
    result = v;
    return result;
}


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
double spdmatrixrcond(ap::real_2d_array a, int n, bool isupper)
{
    double result;
    int i;
    int j;
    int j1;
    int j2;
    double v;
    double nrm;
    ap::real_1d_array t;

    t.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        t(i) = 0;
    }
    for(i = 0; i <= n-1; i++)
    {
        if( isupper )
        {
            j1 = i;
            j2 = n-1;
        }
        else
        {
            j1 = 0;
            j2 = i;
        }
        for(j = j1; j <= j2; j++)
        {
            if( i==j )
            {
                t(i) = t(i)+fabs(a(i,i));
            }
            else
            {
                t(i) = t(i)+fabs(a(i,j));
                t(j) = t(j)+fabs(a(i,j));
            }
        }
    }
    nrm = 0;
    for(i = 0; i <= n-1; i++)
    {
        nrm = ap::maxreal(nrm, t(i));
    }
    if( spdmatrixcholesky(a, n, isupper) )
    {
        spdmatrixrcondcholeskyinternal(a, n, isupper, true, nrm, v);
        result = v;
    }
    else
    {
        result = -1;
    }
    return result;
}


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
     bool isunit)
{
    double result;
    int i;
    int j;
    double v;
    double nrm;
    ap::integer_1d_array pivots;
    ap::real_1d_array t;
    int j1;
    int j2;

    ap::ap_error::make_assertion(n>=1, "RMatrixTRRCond1: N<1!");
    t.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        t(i) = 0;
    }
    for(i = 0; i <= n-1; i++)
    {
        if( isupper )
        {
            j1 = i+1;
            j2 = n-1;
        }
        else
        {
            j1 = 0;
            j2 = i-1;
        }
        for(j = j1; j <= j2; j++)
        {
            t(j) = t(j)+fabs(a(i,j));
        }
        if( isunit )
        {
            t(i) = t(i)+1;
        }
        else
        {
            t(i) = t(i)+fabs(a(i,i));
        }
    }
    nrm = 0;
    for(i = 0; i <= n-1; i++)
    {
        nrm = ap::maxreal(nrm, t(i));
    }
    rmatrixrcondtrinternal(a, n, isupper, isunit, true, nrm, v);
    result = v;
    return result;
}


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
     bool isunit)
{
    double result;
    int i;
    int j;
    double v;
    double nrm;
    ap::integer_1d_array pivots;
    int j1;
    int j2;

    ap::ap_error::make_assertion(n>=1, "RMatrixTRRCondInf: N<1!");
    nrm = 0;
    for(i = 0; i <= n-1; i++)
    {
        if( isupper )
        {
            j1 = i+1;
            j2 = n-1;
        }
        else
        {
            j1 = 0;
            j2 = i-1;
        }
        v = 0;
        for(j = j1; j <= j2; j++)
        {
            v = v+fabs(a(i,j));
        }
        if( isunit )
        {
            v = v+1;
        }
        else
        {
            v = v+fabs(a(i,i));
        }
        nrm = ap::maxreal(nrm, v);
    }
    rmatrixrcondtrinternal(a, n, isupper, isunit, false, nrm, v);
    result = v;
    return result;
}


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
double hpdmatrixrcond(ap::complex_2d_array a, int n, bool isupper)
{
    double result;
    int i;
    int j;
    int j1;
    int j2;
    double v;
    double nrm;
    ap::real_1d_array t;

    t.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        t(i) = 0;
    }
    for(i = 0; i <= n-1; i++)
    {
        if( isupper )
        {
            j1 = i;
            j2 = n-1;
        }
        else
        {
            j1 = 0;
            j2 = i;
        }
        for(j = j1; j <= j2; j++)
        {
            if( i==j )
            {
                t(i) = t(i)+ap::abscomplex(a(i,i));
            }
            else
            {
                t(i) = t(i)+ap::abscomplex(a(i,j));
                t(j) = t(j)+ap::abscomplex(a(i,j));
            }
        }
    }
    nrm = 0;
    for(i = 0; i <= n-1; i++)
    {
        nrm = ap::maxreal(nrm, t(i));
    }
    if( hpdmatrixcholesky(a, n, isupper) )
    {
        hpdmatrixrcondcholeskyinternal(a, n, isupper, true, nrm, v);
        result = v;
    }
    else
    {
        result = -1;
    }
    return result;
}


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
double cmatrixrcond1(ap::complex_2d_array a, int n)
{
    double result;
    int i;
    int j;
    double v;
    double nrm;
    ap::integer_1d_array pivots;
    ap::real_1d_array t;

    ap::ap_error::make_assertion(n>=1, "CMatrixRCond1: N<1!");
    t.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        t(i) = 0;
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            t(j) = t(j)+ap::abscomplex(a(i,j));
        }
    }
    nrm = 0;
    for(i = 0; i <= n-1; i++)
    {
        nrm = ap::maxreal(nrm, t(i));
    }
    cmatrixlu(a, n, n, pivots);
    cmatrixrcondluinternal(a, n, true, true, nrm, v);
    result = v;
    return result;
}


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
double cmatrixrcondinf(ap::complex_2d_array a, int n)
{
    double result;
    int i;
    int j;
    double v;
    double nrm;
    ap::integer_1d_array pivots;

    ap::ap_error::make_assertion(n>=1, "CMatrixRCondInf: N<1!");
    nrm = 0;
    for(i = 0; i <= n-1; i++)
    {
        v = 0;
        for(j = 0; j <= n-1; j++)
        {
            v = v+ap::abscomplex(a(i,j));
        }
        nrm = ap::maxreal(nrm, v);
    }
    cmatrixlu(a, n, n, pivots);
    cmatrixrcondluinternal(a, n, false, true, nrm, v);
    result = v;
    return result;
}


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
double rmatrixlurcond1(const ap::real_2d_array& lua, int n)
{
    double result;
    double v;

    rmatrixrcondluinternal(lua, n, true, false, double(0), v);
    result = v;
    return result;
}


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
double rmatrixlurcondinf(const ap::real_2d_array& lua, int n)
{
    double result;
    double v;

    rmatrixrcondluinternal(lua, n, false, false, double(0), v);
    result = v;
    return result;
}


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
double spdmatrixcholeskyrcond(const ap::real_2d_array& a,
     int n,
     bool isupper)
{
    double result;
    double v;

    spdmatrixrcondcholeskyinternal(a, n, isupper, false, double(0), v);
    result = v;
    return result;
}


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
     bool isupper)
{
    double result;
    double v;

    hpdmatrixrcondcholeskyinternal(a, n, isupper, false, double(0), v);
    result = v;
    return result;
}


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
double cmatrixlurcond1(const ap::complex_2d_array& lua, int n)
{
    double result;
    double v;

    ap::ap_error::make_assertion(n>=1, "CMatrixLURCond1: N<1!");
    cmatrixrcondluinternal(lua, n, true, false, 0.0, v);
    result = v;
    return result;
}


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
double cmatrixlurcondinf(const ap::complex_2d_array& lua, int n)
{
    double result;
    double v;

    ap::ap_error::make_assertion(n>=1, "CMatrixLURCondInf: N<1!");
    cmatrixrcondluinternal(lua, n, false, false, 0.0, v);
    result = v;
    return result;
}


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
     bool isunit)
{
    double result;
    int i;
    int j;
    double v;
    double nrm;
    ap::integer_1d_array pivots;
    ap::real_1d_array t;
    int j1;
    int j2;

    ap::ap_error::make_assertion(n>=1, "RMatrixTRRCond1: N<1!");
    t.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        t(i) = 0;
    }
    for(i = 0; i <= n-1; i++)
    {
        if( isupper )
        {
            j1 = i+1;
            j2 = n-1;
        }
        else
        {
            j1 = 0;
            j2 = i-1;
        }
        for(j = j1; j <= j2; j++)
        {
            t(j) = t(j)+ap::abscomplex(a(i,j));
        }
        if( isunit )
        {
            t(i) = t(i)+1;
        }
        else
        {
            t(i) = t(i)+ap::abscomplex(a(i,i));
        }
    }
    nrm = 0;
    for(i = 0; i <= n-1; i++)
    {
        nrm = ap::maxreal(nrm, t(i));
    }
    cmatrixrcondtrinternal(a, n, isupper, isunit, true, nrm, v);
    result = v;
    return result;
}


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
     bool isunit)
{
    double result;
    int i;
    int j;
    double v;
    double nrm;
    ap::integer_1d_array pivots;
    int j1;
    int j2;

    ap::ap_error::make_assertion(n>=1, "RMatrixTRRCondInf: N<1!");
    nrm = 0;
    for(i = 0; i <= n-1; i++)
    {
        if( isupper )
        {
            j1 = i+1;
            j2 = n-1;
        }
        else
        {
            j1 = 0;
            j2 = i-1;
        }
        v = 0;
        for(j = j1; j <= j2; j++)
        {
            v = v+ap::abscomplex(a(i,j));
        }
        if( isunit )
        {
            v = v+1;
        }
        else
        {
            v = v+ap::abscomplex(a(i,i));
        }
        nrm = ap::maxreal(nrm, v);
    }
    cmatrixrcondtrinternal(a, n, isupper, isunit, false, nrm, v);
    result = v;
    return result;
}


/*************************************************************************
Threshold for rcond: matrices with condition number beyond this  threshold
are considered singular.

Threshold must be far enough from underflow, at least Sqr(Threshold)  must
be greater than underflow.
*************************************************************************/
double rcondthreshold()
{
    double result;

    result = sqrt(sqrt(ap::minrealnumber));
    return result;
}


/*************************************************************************
Internal subroutine for condition number estimation

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************/
static void rmatrixrcondtrinternal(const ap::real_2d_array& a,
     int n,
     bool isupper,
     bool isunit,
     bool onenorm,
     double anorm,
     double& rc)
{
    ap::real_1d_array ex;
    ap::real_1d_array ev;
    ap::integer_1d_array iwork;
    ap::real_1d_array tmp;
    double v;
    int i;
    int j;
    int kase;
    int kase1;
    int j1;
    int j2;
    double ainvnm;
    double maxgrowth;
    double s;
    bool mupper;
    bool mtrans;
    bool munit;

    
    //
    // RC=0 if something happens
    //
    rc = 0;
    
    //
    // init
    //
    if( onenorm )
    {
        kase1 = 1;
    }
    else
    {
        kase1 = 2;
    }
    mupper = true;
    mtrans = true;
    munit = true;
    iwork.setlength(n+1);
    tmp.setlength(n);
    
    //
    // prepare parameters for triangular solver
    //
    maxgrowth = 1/rcondthreshold();
    s = 0;
    for(i = 0; i <= n-1; i++)
    {
        if( isupper )
        {
            j1 = i+1;
            j2 = n-1;
        }
        else
        {
            j1 = 0;
            j2 = i-1;
        }
        for(j = j1; j <= j2; j++)
        {
            s = ap::maxreal(s, fabs(a(i,j)));
        }
        if( isunit )
        {
            s = ap::maxreal(s, double(1));
        }
        else
        {
            s = ap::maxreal(s, fabs(a(i,i)));
        }
    }
    if( ap::fp_eq(s,0) )
    {
        s = 1;
    }
    s = 1/s;
    
    //
    // Scale according to S
    //
    anorm = anorm*s;
    
    //
    // Quick return if possible
    // We assume that ANORM<>0 after this block
    //
    if( ap::fp_eq(anorm,0) )
    {
        return;
    }
    if( n==1 )
    {
        rc = 1;
        return;
    }
    
    //
    // Estimate the norm of inv(A).
    //
    ainvnm = 0;
    kase = 0;
    while(true)
    {
        rmatrixestimatenorm(n, ev, ex, iwork, ainvnm, kase);
        if( kase==0 )
        {
            break;
        }
        
        //
        // from 1-based array to 0-based
        //
        for(i = 0; i <= n-1; i++)
        {
            ex(i) = ex(i+1);
        }
        
        //
        // multiply by inv(A) or inv(A')
        //
        if( kase==kase1 )
        {
            
            //
            // multiply by inv(A)
            //
            if( !rmatrixscaledtrsafesolve(a, s, n, ex, isupper, 0, isunit, maxgrowth) )
            {
                return;
            }
        }
        else
        {
            
            //
            // multiply by inv(A')
            //
            if( !rmatrixscaledtrsafesolve(a, s, n, ex, isupper, 1, isunit, maxgrowth) )
            {
                return;
            }
        }
        
        //
        // from 0-based array to 1-based
        //
        for(i = n-1; i >= 0; i--)
        {
            ex(i+1) = ex(i);
        }
    }
    
    //
    // Compute the estimate of the reciprocal condition number.
    //
    if( ap::fp_neq(ainvnm,0) )
    {
        rc = 1/ainvnm;
        rc = rc/anorm;
        if( ap::fp_less(rc,rcondthreshold()) )
        {
            rc = 0;
        }
    }
}


/*************************************************************************
Condition number estimation

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     March 31, 1993
*************************************************************************/
static void cmatrixrcondtrinternal(const ap::complex_2d_array& a,
     const int& n,
     bool isupper,
     bool isunit,
     bool onenorm,
     double anorm,
     double& rc)
{
    ap::complex_1d_array ex;
    ap::complex_1d_array cwork2;
    ap::complex_1d_array cwork3;
    ap::complex_1d_array cwork4;
    ap::integer_1d_array isave;
    ap::real_1d_array rsave;
    int kase;
    int kase1;
    double ainvnm;
    ap::complex v;
    int i;
    int j;
    int j1;
    int j2;
    double s;
    double maxgrowth;

    
    //
    // RC=0 if something happens
    //
    rc = 0;
    
    //
    // init
    //
    if( n<=0 )
    {
        return;
    }
    if( n==0 )
    {
        rc = 1;
        return;
    }
    cwork2.setlength(n+1);
    
    //
    // prepare parameters for triangular solver
    //
    maxgrowth = 1/rcondthreshold();
    s = 0;
    for(i = 0; i <= n-1; i++)
    {
        if( isupper )
        {
            j1 = i+1;
            j2 = n-1;
        }
        else
        {
            j1 = 0;
            j2 = i-1;
        }
        for(j = j1; j <= j2; j++)
        {
            s = ap::maxreal(s, ap::abscomplex(a(i,j)));
        }
        if( isunit )
        {
            s = ap::maxreal(s, double(1));
        }
        else
        {
            s = ap::maxreal(s, ap::abscomplex(a(i,i)));
        }
    }
    if( ap::fp_eq(s,0) )
    {
        s = 1;
    }
    s = 1/s;
    
    //
    // Scale according to S
    //
    anorm = anorm*s;
    
    //
    // Quick return if possible
    //
    if( ap::fp_eq(anorm,0) )
    {
        return;
    }
    
    //
    // Estimate the norm of inv(A).
    //
    ainvnm = 0;
    if( onenorm )
    {
        kase1 = 1;
    }
    else
    {
        kase1 = 2;
    }
    kase = 0;
    while(true)
    {
        cmatrixestimatenorm(n, cwork4, ex, ainvnm, kase, isave, rsave);
        if( kase==0 )
        {
            break;
        }
        
        //
        // From 1-based to 0-based
        //
        for(i = 0; i <= n-1; i++)
        {
            ex(i) = ex(i+1);
        }
        
        //
        // multiply by inv(A) or inv(A')
        //
        if( kase==kase1 )
        {
            
            //
            // multiply by inv(A)
            //
            if( !cmatrixscaledtrsafesolve(a, s, n, ex, isupper, 0, isunit, maxgrowth) )
            {
                return;
            }
        }
        else
        {
            
            //
            // multiply by inv(A')
            //
            if( !cmatrixscaledtrsafesolve(a, s, n, ex, isupper, 2, isunit, maxgrowth) )
            {
                return;
            }
        }
        
        //
        // from 0-based to 1-based
        //
        for(i = n-1; i >= 0; i--)
        {
            ex(i+1) = ex(i);
        }
    }
    
    //
    // Compute the estimate of the reciprocal condition number.
    //
    if( ap::fp_neq(ainvnm,0) )
    {
        rc = 1/ainvnm;
        rc = rc/anorm;
        if( ap::fp_less(rc,rcondthreshold()) )
        {
            rc = 0;
        }
    }
}


/*************************************************************************
Internal subroutine for condition number estimation

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************/
static void spdmatrixrcondcholeskyinternal(const ap::real_2d_array& cha,
     int n,
     bool isupper,
     bool isnormprovided,
     double anorm,
     double& rc)
{
    int i;
    int j;
    int kase;
    double ainvnm;
    ap::real_1d_array ex;
    ap::real_1d_array ev;
    ap::real_1d_array tmp;
    ap::integer_1d_array iwork;
    double sa;
    double v;
    double maxgrowth;

    ap::ap_error::make_assertion(n>=1, "");
    tmp.setlength(n);
    
    //
    // RC=0 if something happens
    //
    rc = 0;
    
    //
    // prepare parameters for triangular solver
    //
    maxgrowth = 1/rcondthreshold();
    sa = 0;
    if( isupper )
    {
        for(i = 0; i <= n-1; i++)
        {
            for(j = i; j <= n-1; j++)
            {
                sa = ap::maxreal(sa, ap::abscomplex(cha(i,j)));
            }
        }
    }
    else
    {
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= i; j++)
            {
                sa = ap::maxreal(sa, ap::abscomplex(cha(i,j)));
            }
        }
    }
    if( ap::fp_eq(sa,0) )
    {
        sa = 1;
    }
    sa = 1/sa;
    
    //
    // Estimate the norm of A.
    //
    if( !isnormprovided )
    {
        kase = 0;
        anorm = 0;
        while(true)
        {
            rmatrixestimatenorm(n, ev, ex, iwork, anorm, kase);
            if( kase==0 )
            {
                break;
            }
            if( isupper )
            {
                
                //
                // Multiply by U
                //
                for(i = 1; i <= n; i++)
                {
                    v = ap::vdotproduct(&cha(i-1, i-1), 1, &ex(i), 1, ap::vlen(i-1,n-1));
                    ex(i) = v;
                }
                ap::vmul(&ex(1), 1, ap::vlen(1,n), sa);
                
                //
                // Multiply by U'
                //
                for(i = 0; i <= n-1; i++)
                {
                    tmp(i) = 0;
                }
                for(i = 0; i <= n-1; i++)
                {
                    v = ex(i+1);
                    ap::vadd(&tmp(i), 1, &cha(i, i), 1, ap::vlen(i,n-1), v);
                }
                ap::vmove(&ex(1), 1, &tmp(0), 1, ap::vlen(1,n));
                ap::vmul(&ex(1), 1, ap::vlen(1,n), sa);
            }
            else
            {
                
                //
                // Multiply by L'
                //
                for(i = 0; i <= n-1; i++)
                {
                    tmp(i) = 0;
                }
                for(i = 0; i <= n-1; i++)
                {
                    v = ex(i+1);
                    ap::vadd(&tmp(0), 1, &cha(i, 0), 1, ap::vlen(0,i), v);
                }
                ap::vmove(&ex(1), 1, &tmp(0), 1, ap::vlen(1,n));
                ap::vmul(&ex(1), 1, ap::vlen(1,n), sa);
                
                //
                // Multiply by L
                //
                for(i = n; i >= 1; i--)
                {
                    v = ap::vdotproduct(&cha(i-1, 0), 1, &ex(1), 1, ap::vlen(0,i-1));
                    ex(i) = v;
                }
                ap::vmul(&ex(1), 1, ap::vlen(1,n), sa);
            }
        }
    }
    
    //
    // Quick return if possible
    //
    if( ap::fp_eq(anorm,0) )
    {
        return;
    }
    if( n==1 )
    {
        rc = 1;
        return;
    }
    
    //
    // Estimate the 1-norm of inv(A).
    //
    kase = 0;
    while(true)
    {
        rmatrixestimatenorm(n, ev, ex, iwork, ainvnm, kase);
        if( kase==0 )
        {
            break;
        }
        for(i = 0; i <= n-1; i++)
        {
            ex(i) = ex(i+1);
        }
        if( isupper )
        {
            
            //
            // Multiply by inv(U').
            //
            if( !rmatrixscaledtrsafesolve(cha, sa, n, ex, isupper, 1, false, maxgrowth) )
            {
                return;
            }
            
            //
            // Multiply by inv(U).
            //
            if( !rmatrixscaledtrsafesolve(cha, sa, n, ex, isupper, 0, false, maxgrowth) )
            {
                return;
            }
        }
        else
        {
            
            //
            // Multiply by inv(L).
            //
            if( !rmatrixscaledtrsafesolve(cha, sa, n, ex, isupper, 0, false, maxgrowth) )
            {
                return;
            }
            
            //
            // Multiply by inv(L').
            //
            if( !rmatrixscaledtrsafesolve(cha, sa, n, ex, isupper, 1, false, maxgrowth) )
            {
                return;
            }
        }
        for(i = n-1; i >= 0; i--)
        {
            ex(i+1) = ex(i);
        }
    }
    
    //
    // Compute the estimate of the reciprocal condition number.
    //
    if( ap::fp_neq(ainvnm,0) )
    {
        v = 1/ainvnm;
        rc = v/anorm;
        if( ap::fp_less(rc,rcondthreshold()) )
        {
            rc = 0;
        }
    }
}


/*************************************************************************
Internal subroutine for condition number estimation

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************/
static void hpdmatrixrcondcholeskyinternal(const ap::complex_2d_array& cha,
     int n,
     bool isupper,
     bool isnormprovided,
     double anorm,
     double& rc)
{
    ap::integer_1d_array isave;
    ap::real_1d_array rsave;
    ap::complex_1d_array ex;
    ap::complex_1d_array ev;
    ap::complex_1d_array tmp;
    int kase;
    double ainvnm;
    ap::complex v;
    int i;
    int j;
    double sa;
    double maxgrowth;

    ap::ap_error::make_assertion(n>=1, "");
    tmp.setlength(n);
    
    //
    // RC=0 if something happens
    //
    rc = 0;
    
    //
    // prepare parameters for triangular solver
    //
    maxgrowth = 1/rcondthreshold();
    sa = 0;
    if( isupper )
    {
        for(i = 0; i <= n-1; i++)
        {
            for(j = i; j <= n-1; j++)
            {
                sa = ap::maxreal(sa, ap::abscomplex(cha(i,j)));
            }
        }
    }
    else
    {
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= i; j++)
            {
                sa = ap::maxreal(sa, ap::abscomplex(cha(i,j)));
            }
        }
    }
    if( ap::fp_eq(sa,0) )
    {
        sa = 1;
    }
    sa = 1/sa;
    
    //
    // Estimate the norm of A
    //
    if( !isnormprovided )
    {
        anorm = 0;
        kase = 0;
        while(true)
        {
            cmatrixestimatenorm(n, ev, ex, anorm, kase, isave, rsave);
            if( kase==0 )
            {
                break;
            }
            if( isupper )
            {
                
                //
                // Multiply by U
                //
                for(i = 1; i <= n; i++)
                {
                    v = ap::vdotproduct(&cha(i-1, i-1), 1, "N", &ex(i), 1, "N", ap::vlen(i-1,n-1));
                    ex(i) = v;
                }
                ap::vmul(&ex(1), 1, ap::vlen(1,n), sa);
                
                //
                // Multiply by U'
                //
                for(i = 0; i <= n-1; i++)
                {
                    tmp(i) = 0;
                }
                for(i = 0; i <= n-1; i++)
                {
                    v = ex(i+1);
                    ap::vadd(&tmp(i), 1, &cha(i, i), 1, "Conj", ap::vlen(i,n-1), v);
                }
                ap::vmove(&ex(1), 1, &tmp(0), 1, "N", ap::vlen(1,n));
                ap::vmul(&ex(1), 1, ap::vlen(1,n), sa);
            }
            else
            {
                
                //
                // Multiply by L'
                //
                for(i = 0; i <= n-1; i++)
                {
                    tmp(i) = 0;
                }
                for(i = 0; i <= n-1; i++)
                {
                    v = ex(i+1);
                    ap::vadd(&tmp(0), 1, &cha(i, 0), 1, "Conj", ap::vlen(0,i), v);
                }
                ap::vmove(&ex(1), 1, &tmp(0), 1, "N", ap::vlen(1,n));
                ap::vmul(&ex(1), 1, ap::vlen(1,n), sa);
                
                //
                // Multiply by L
                //
                for(i = n; i >= 1; i--)
                {
                    v = ap::vdotproduct(&cha(i-1, 0), 1, "N", &ex(1), 1, "N", ap::vlen(0,i-1));
                    ex(i) = v;
                }
                ap::vmul(&ex(1), 1, ap::vlen(1,n), sa);
            }
        }
    }
    
    //
    // Quick return if possible
    // After this block we assume that ANORM<>0
    //
    if( ap::fp_eq(anorm,0) )
    {
        return;
    }
    if( n==1 )
    {
        rc = 1;
        return;
    }
    
    //
    // Estimate the norm of inv(A).
    //
    ainvnm = 0;
    kase = 0;
    while(true)
    {
        cmatrixestimatenorm(n, ev, ex, ainvnm, kase, isave, rsave);
        if( kase==0 )
        {
            break;
        }
        for(i = 0; i <= n-1; i++)
        {
            ex(i) = ex(i+1);
        }
        if( isupper )
        {
            
            //
            // Multiply by inv(U').
            //
            if( !cmatrixscaledtrsafesolve(cha, sa, n, ex, isupper, 2, false, maxgrowth) )
            {
                return;
            }
            
            //
            // Multiply by inv(U).
            //
            if( !cmatrixscaledtrsafesolve(cha, sa, n, ex, isupper, 0, false, maxgrowth) )
            {
                return;
            }
        }
        else
        {
            
            //
            // Multiply by inv(L).
            //
            if( !cmatrixscaledtrsafesolve(cha, sa, n, ex, isupper, 0, false, maxgrowth) )
            {
                return;
            }
            
            //
            // Multiply by inv(L').
            //
            if( !cmatrixscaledtrsafesolve(cha, sa, n, ex, isupper, 2, false, maxgrowth) )
            {
                return;
            }
        }
        for(i = n-1; i >= 0; i--)
        {
            ex(i+1) = ex(i);
        }
    }
    
    //
    // Compute the estimate of the reciprocal condition number.
    //
    if( ap::fp_neq(ainvnm,0) )
    {
        rc = 1/ainvnm;
        rc = rc/anorm;
        if( ap::fp_less(rc,rcondthreshold()) )
        {
            rc = 0;
        }
    }
}


/*************************************************************************
Internal subroutine for condition number estimation

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************/
static void rmatrixrcondluinternal(const ap::real_2d_array& lua,
     int n,
     bool onenorm,
     bool isanormprovided,
     double anorm,
     double& rc)
{
    ap::real_1d_array ex;
    ap::real_1d_array ev;
    ap::integer_1d_array iwork;
    ap::real_1d_array tmp;
    double v;
    int i;
    int j;
    int kase;
    int kase1;
    double ainvnm;
    double maxgrowth;
    double su;
    double sl;
    bool mupper;
    bool mtrans;
    bool munit;

    
    //
    // RC=0 if something happens
    //
    rc = 0;
    
    //
    // init
    //
    if( onenorm )
    {
        kase1 = 1;
    }
    else
    {
        kase1 = 2;
    }
    mupper = true;
    mtrans = true;
    munit = true;
    iwork.setlength(n+1);
    tmp.setlength(n);
    
    //
    // prepare parameters for triangular solver
    //
    maxgrowth = 1/rcondthreshold();
    su = 0;
    sl = 1;
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= i-1; j++)
        {
            sl = ap::maxreal(sl, fabs(lua(i,j)));
        }
        for(j = i; j <= n-1; j++)
        {
            su = ap::maxreal(su, fabs(lua(i,j)));
        }
    }
    if( ap::fp_eq(su,0) )
    {
        su = 1;
    }
    su = 1/su;
    sl = 1/sl;
    
    //
    // Estimate the norm of A.
    //
    if( !isanormprovided )
    {
        kase = 0;
        anorm = 0;
        while(true)
        {
            rmatrixestimatenorm(n, ev, ex, iwork, anorm, kase);
            if( kase==0 )
            {
                break;
            }
            if( kase==kase1 )
            {
                
                //
                // Multiply by U
                //
                for(i = 1; i <= n; i++)
                {
                    v = ap::vdotproduct(&lua(i-1, i-1), 1, &ex(i), 1, ap::vlen(i-1,n-1));
                    ex(i) = v;
                }
                
                //
                // Multiply by L
                //
                for(i = n; i >= 1; i--)
                {
                    if( i>1 )
                    {
                        v = ap::vdotproduct(&lua(i-1, 0), 1, &ex(1), 1, ap::vlen(0,i-2));
                    }
                    else
                    {
                        v = 0;
                    }
                    ex(i) = ex(i)+v;
                }
            }
            else
            {
                
                //
                // Multiply by L'
                //
                for(i = 0; i <= n-1; i++)
                {
                    tmp(i) = 0;
                }
                for(i = 0; i <= n-1; i++)
                {
                    v = ex(i+1);
                    if( i>=1 )
                    {
                        ap::vadd(&tmp(0), 1, &lua(i, 0), 1, ap::vlen(0,i-1), v);
                    }
                    tmp(i) = tmp(i)+v;
                }
                ap::vmove(&ex(1), 1, &tmp(0), 1, ap::vlen(1,n));
                
                //
                // Multiply by U'
                //
                for(i = 0; i <= n-1; i++)
                {
                    tmp(i) = 0;
                }
                for(i = 0; i <= n-1; i++)
                {
                    v = ex(i+1);
                    ap::vadd(&tmp(i), 1, &lua(i, i), 1, ap::vlen(i,n-1), v);
                }
                ap::vmove(&ex(1), 1, &tmp(0), 1, ap::vlen(1,n));
            }
        }
    }
    
    //
    // Scale according to SU/SL
    //
    anorm = anorm*su*sl;
    
    //
    // Quick return if possible
    // We assume that ANORM<>0 after this block
    //
    if( ap::fp_eq(anorm,0) )
    {
        return;
    }
    if( n==1 )
    {
        rc = 1;
        return;
    }
    
    //
    // Estimate the norm of inv(A).
    //
    ainvnm = 0;
    kase = 0;
    while(true)
    {
        rmatrixestimatenorm(n, ev, ex, iwork, ainvnm, kase);
        if( kase==0 )
        {
            break;
        }
        
        //
        // from 1-based array to 0-based
        //
        for(i = 0; i <= n-1; i++)
        {
            ex(i) = ex(i+1);
        }
        
        //
        // multiply by inv(A) or inv(A')
        //
        if( kase==kase1 )
        {
            
            //
            // Multiply by inv(L).
            //
            if( !rmatrixscaledtrsafesolve(lua, sl, n, ex, !mupper, 0, munit, maxgrowth) )
            {
                return;
            }
            
            //
            // Multiply by inv(U).
            //
            if( !rmatrixscaledtrsafesolve(lua, su, n, ex, mupper, 0, !munit, maxgrowth) )
            {
                return;
            }
        }
        else
        {
            
            //
            // Multiply by inv(U').
            //
            if( !rmatrixscaledtrsafesolve(lua, su, n, ex, mupper, 1, !munit, maxgrowth) )
            {
                return;
            }
            
            //
            // Multiply by inv(L').
            //
            if( !rmatrixscaledtrsafesolve(lua, sl, n, ex, !mupper, 1, munit, maxgrowth) )
            {
                return;
            }
        }
        
        //
        // from 0-based array to 1-based
        //
        for(i = n-1; i >= 0; i--)
        {
            ex(i+1) = ex(i);
        }
    }
    
    //
    // Compute the estimate of the reciprocal condition number.
    //
    if( ap::fp_neq(ainvnm,0) )
    {
        rc = 1/ainvnm;
        rc = rc/anorm;
        if( ap::fp_less(rc,rcondthreshold()) )
        {
            rc = 0;
        }
    }
}


/*************************************************************************
Condition number estimation

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     March 31, 1993
*************************************************************************/
static void cmatrixrcondluinternal(const ap::complex_2d_array& lua,
     const int& n,
     bool onenorm,
     bool isanormprovided,
     double anorm,
     double& rc)
{
    ap::complex_1d_array ex;
    ap::complex_1d_array cwork2;
    ap::complex_1d_array cwork3;
    ap::complex_1d_array cwork4;
    ap::integer_1d_array isave;
    ap::real_1d_array rsave;
    int kase;
    int kase1;
    double ainvnm;
    ap::complex v;
    int i;
    int j;
    double su;
    double sl;
    double maxgrowth;

    if( n<=0 )
    {
        return;
    }
    cwork2.setlength(n+1);
    rc = 0;
    if( n==0 )
    {
        rc = 1;
        return;
    }
    
    //
    // prepare parameters for triangular solver
    //
    maxgrowth = 1/rcondthreshold();
    su = 0;
    sl = 1;
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= i-1; j++)
        {
            sl = ap::maxreal(sl, ap::abscomplex(lua(i,j)));
        }
        for(j = i; j <= n-1; j++)
        {
            su = ap::maxreal(su, ap::abscomplex(lua(i,j)));
        }
    }
    if( ap::fp_eq(su,0) )
    {
        su = 1;
    }
    su = 1/su;
    sl = 1/sl;
    
    //
    // Estimate the norm of SU*SL*A.
    //
    if( !isanormprovided )
    {
        anorm = 0;
        if( onenorm )
        {
            kase1 = 1;
        }
        else
        {
            kase1 = 2;
        }
        kase = 0;
        do
        {
            cmatrixestimatenorm(n, cwork4, ex, anorm, kase, isave, rsave);
            if( kase!=0 )
            {
                if( kase==kase1 )
                {
                    
                    //
                    // Multiply by U
                    //
                    for(i = 1; i <= n; i++)
                    {
                        v = ap::vdotproduct(&lua(i-1, i-1), 1, "N", &ex(i), 1, "N", ap::vlen(i-1,n-1));
                        ex(i) = v;
                    }
                    
                    //
                    // Multiply by L
                    //
                    for(i = n; i >= 1; i--)
                    {
                        v = 0;
                        if( i>1 )
                        {
                            v = ap::vdotproduct(&lua(i-1, 0), 1, "N", &ex(1), 1, "N", ap::vlen(0,i-2));
                        }
                        ex(i) = v+ex(i);
                    }
                }
                else
                {
                    
                    //
                    // Multiply by L'
                    //
                    for(i = 1; i <= n; i++)
                    {
                        cwork2(i) = 0;
                    }
                    for(i = 1; i <= n; i++)
                    {
                        v = ex(i);
                        if( i>1 )
                        {
                            ap::vadd(&cwork2(1), 1, &lua(i-1, 0), 1, "Conj", ap::vlen(1,i-1), v);
                        }
                        cwork2(i) = cwork2(i)+v;
                    }
                    
                    //
                    // Multiply by U'
                    //
                    for(i = 1; i <= n; i++)
                    {
                        ex(i) = 0;
                    }
                    for(i = 1; i <= n; i++)
                    {
                        v = cwork2(i);
                        ap::vadd(&ex(i), 1, &lua(i-1, i-1), 1, "Conj", ap::vlen(i,n), v);
                    }
                }
            }
        }
        while(kase!=0);
    }
    
    //
    // Scale according to SU/SL
    //
    anorm = anorm*su*sl;
    
    //
    // Quick return if possible
    //
    if( ap::fp_eq(anorm,0) )
    {
        return;
    }
    
    //
    // Estimate the norm of inv(A).
    //
    ainvnm = 0;
    if( onenorm )
    {
        kase1 = 1;
    }
    else
    {
        kase1 = 2;
    }
    kase = 0;
    while(true)
    {
        cmatrixestimatenorm(n, cwork4, ex, ainvnm, kase, isave, rsave);
        if( kase==0 )
        {
            break;
        }
        
        //
        // From 1-based to 0-based
        //
        for(i = 0; i <= n-1; i++)
        {
            ex(i) = ex(i+1);
        }
        
        //
        // multiply by inv(A) or inv(A')
        //
        if( kase==kase1 )
        {
            
            //
            // Multiply by inv(L).
            //
            if( !cmatrixscaledtrsafesolve(lua, sl, n, ex, false, 0, true, maxgrowth) )
            {
                rc = 0;
                return;
            }
            
            //
            // Multiply by inv(U).
            //
            if( !cmatrixscaledtrsafesolve(lua, su, n, ex, true, 0, false, maxgrowth) )
            {
                rc = 0;
                return;
            }
        }
        else
        {
            
            //
            // Multiply by inv(U').
            //
            if( !cmatrixscaledtrsafesolve(lua, su, n, ex, true, 2, false, maxgrowth) )
            {
                rc = 0;
                return;
            }
            
            //
            // Multiply by inv(L').
            //
            if( !cmatrixscaledtrsafesolve(lua, sl, n, ex, false, 2, true, maxgrowth) )
            {
                rc = 0;
                return;
            }
        }
        
        //
        // from 0-based to 1-based
        //
        for(i = n-1; i >= 0; i--)
        {
            ex(i+1) = ex(i);
        }
    }
    
    //
    // Compute the estimate of the reciprocal condition number.
    //
    if( ap::fp_neq(ainvnm,0) )
    {
        rc = 1/ainvnm;
        rc = rc/anorm;
        if( ap::fp_less(rc,rcondthreshold()) )
        {
            rc = 0;
        }
    }
}


/*************************************************************************
Internal subroutine for matrix norm estimation

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************/
static void rmatrixestimatenorm(int n,
     ap::real_1d_array& v,
     ap::real_1d_array& x,
     ap::integer_1d_array& isgn,
     double& est,
     int& kase)
{
    int itmax;
    int i;
    double t;
    bool flg;
    int positer;
    int posj;
    int posjlast;
    int posjump;
    int posaltsgn;
    int posestold;
    int postemp;

    itmax = 5;
    posaltsgn = n+1;
    posestold = n+2;
    postemp = n+3;
    positer = n+1;
    posj = n+2;
    posjlast = n+3;
    posjump = n+4;
    if( kase==0 )
    {
        v.setlength(n+4);
        x.setlength(n+1);
        isgn.setlength(n+5);
        t = double(1)/double(n);
        for(i = 1; i <= n; i++)
        {
            x(i) = t;
        }
        kase = 1;
        isgn(posjump) = 1;
        return;
    }
    
    //
    //     ................ ENTRY   (JUMP = 1)
    //     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
    //
    if( isgn(posjump)==1 )
    {
        if( n==1 )
        {
            v(1) = x(1);
            est = fabs(v(1));
            kase = 0;
            return;
        }
        est = 0;
        for(i = 1; i <= n; i++)
        {
            est = est+fabs(x(i));
        }
        for(i = 1; i <= n; i++)
        {
            if( ap::fp_greater_eq(x(i),0) )
            {
                x(i) = 1;
            }
            else
            {
                x(i) = -1;
            }
            isgn(i) = ap::sign(x(i));
        }
        kase = 2;
        isgn(posjump) = 2;
        return;
    }
    
    //
    //     ................ ENTRY   (JUMP = 2)
    //     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X.
    //
    if( isgn(posjump)==2 )
    {
        isgn(posj) = 1;
        for(i = 2; i <= n; i++)
        {
            if( ap::fp_greater(fabs(x(i)),fabs(x(isgn(posj)))) )
            {
                isgn(posj) = i;
            }
        }
        isgn(positer) = 2;
        
        //
        // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
        //
        for(i = 1; i <= n; i++)
        {
            x(i) = 0;
        }
        x(isgn(posj)) = 1;
        kase = 1;
        isgn(posjump) = 3;
        return;
    }
    
    //
    //     ................ ENTRY   (JUMP = 3)
    //     X HAS BEEN OVERWRITTEN BY A*X.
    //
    if( isgn(posjump)==3 )
    {
        ap::vmove(&v(1), 1, &x(1), 1, ap::vlen(1,n));
        v(posestold) = est;
        est = 0;
        for(i = 1; i <= n; i++)
        {
            est = est+fabs(v(i));
        }
        flg = false;
        for(i = 1; i <= n; i++)
        {
            if( ap::fp_greater_eq(x(i),0)&&isgn(i)<0||ap::fp_less(x(i),0)&&isgn(i)>=0 )
            {
                flg = true;
            }
        }
        
        //
        // REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
        // OR MAY BE CYCLING.
        //
        if( !flg||ap::fp_less_eq(est,v(posestold)) )
        {
            v(posaltsgn) = 1;
            for(i = 1; i <= n; i++)
            {
                x(i) = v(posaltsgn)*(1+double(i-1)/double(n-1));
                v(posaltsgn) = -v(posaltsgn);
            }
            kase = 1;
            isgn(posjump) = 5;
            return;
        }
        for(i = 1; i <= n; i++)
        {
            if( ap::fp_greater_eq(x(i),0) )
            {
                x(i) = 1;
                isgn(i) = 1;
            }
            else
            {
                x(i) = -1;
                isgn(i) = -1;
            }
        }
        kase = 2;
        isgn(posjump) = 4;
        return;
    }
    
    //
    //     ................ ENTRY   (JUMP = 4)
    //     X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X.
    //
    if( isgn(posjump)==4 )
    {
        isgn(posjlast) = isgn(posj);
        isgn(posj) = 1;
        for(i = 2; i <= n; i++)
        {
            if( ap::fp_greater(fabs(x(i)),fabs(x(isgn(posj)))) )
            {
                isgn(posj) = i;
            }
        }
        if( ap::fp_neq(x(isgn(posjlast)),fabs(x(isgn(posj))))&&isgn(positer)<itmax )
        {
            isgn(positer) = isgn(positer)+1;
            for(i = 1; i <= n; i++)
            {
                x(i) = 0;
            }
            x(isgn(posj)) = 1;
            kase = 1;
            isgn(posjump) = 3;
            return;
        }
        
        //
        // ITERATION COMPLETE.  FINAL STAGE.
        //
        v(posaltsgn) = 1;
        for(i = 1; i <= n; i++)
        {
            x(i) = v(posaltsgn)*(1+double(i-1)/double(n-1));
            v(posaltsgn) = -v(posaltsgn);
        }
        kase = 1;
        isgn(posjump) = 5;
        return;
    }
    
    //
    //     ................ ENTRY   (JUMP = 5)
    //     X HAS BEEN OVERWRITTEN BY A*X.
    //
    if( isgn(posjump)==5 )
    {
        v(postemp) = 0;
        for(i = 1; i <= n; i++)
        {
            v(postemp) = v(postemp)+fabs(x(i));
        }
        v(postemp) = 2*v(postemp)/(3*n);
        if( ap::fp_greater(v(postemp),est) )
        {
            ap::vmove(&v(1), 1, &x(1), 1, ap::vlen(1,n));
            est = v(postemp);
        }
        kase = 0;
        return;
    }
}


static void cmatrixestimatenorm(const int& n,
     ap::complex_1d_array& v,
     ap::complex_1d_array& x,
     double& est,
     int& kase,
     ap::integer_1d_array& isave,
     ap::real_1d_array& rsave)
{
    int itmax;
    int i;
    int iter;
    int j;
    int jlast;
    int jump;
    double absxi;
    double altsgn;
    double estold;
    double safmin;
    double temp;

    
    //
    //Executable Statements ..
    //
    itmax = 5;
    safmin = ap::minrealnumber;
    if( kase==0 )
    {
        v.setlength(n+1);
        x.setlength(n+1);
        isave.setlength(5);
        rsave.setlength(4);
        for(i = 1; i <= n; i++)
        {
            x(i) = double(1)/double(n);
        }
        kase = 1;
        jump = 1;
        internalcomplexrcondsaveall(isave, rsave, i, iter, j, jlast, jump, absxi, altsgn, estold, temp);
        return;
    }
    internalcomplexrcondloadall(isave, rsave, i, iter, j, jlast, jump, absxi, altsgn, estold, temp);
    
    //
    // ENTRY   (JUMP = 1)
    // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
    //
    if( jump==1 )
    {
        if( n==1 )
        {
            v(1) = x(1);
            est = ap::abscomplex(v(1));
            kase = 0;
            internalcomplexrcondsaveall(isave, rsave, i, iter, j, jlast, jump, absxi, altsgn, estold, temp);
            return;
        }
        est = internalcomplexrcondscsum1(x, n);
        for(i = 1; i <= n; i++)
        {
            absxi = ap::abscomplex(x(i));
            if( ap::fp_greater(absxi,safmin) )
            {
                x(i) = x(i)/absxi;
            }
            else
            {
                x(i) = 1;
            }
        }
        kase = 2;
        jump = 2;
        internalcomplexrcondsaveall(isave, rsave, i, iter, j, jlast, jump, absxi, altsgn, estold, temp);
        return;
    }
    
    //
    // ENTRY   (JUMP = 2)
    // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
    //
    if( jump==2 )
    {
        j = internalcomplexrcondicmax1(x, n);
        iter = 2;
        
        //
        // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
        //
        for(i = 1; i <= n; i++)
        {
            x(i) = 0;
        }
        x(j) = 1;
        kase = 1;
        jump = 3;
        internalcomplexrcondsaveall(isave, rsave, i, iter, j, jlast, jump, absxi, altsgn, estold, temp);
        return;
    }
    
    //
    // ENTRY   (JUMP = 3)
    // X HAS BEEN OVERWRITTEN BY A*X.
    //
    if( jump==3 )
    {
        ap::vmove(&v(1), 1, &x(1), 1, "N", ap::vlen(1,n));
        estold = est;
        est = internalcomplexrcondscsum1(v, n);
        
        //
        // TEST FOR CYCLING.
        //
        if( ap::fp_less_eq(est,estold) )
        {
            
            //
            // ITERATION COMPLETE.  FINAL STAGE.
            //
            altsgn = 1;
            for(i = 1; i <= n; i++)
            {
                x(i) = altsgn*(1+double(i-1)/double(n-1));
                altsgn = -altsgn;
            }
            kase = 1;
            jump = 5;
            internalcomplexrcondsaveall(isave, rsave, i, iter, j, jlast, jump, absxi, altsgn, estold, temp);
            return;
        }
        for(i = 1; i <= n; i++)
        {
            absxi = ap::abscomplex(x(i));
            if( ap::fp_greater(absxi,safmin) )
            {
                x(i) = x(i)/absxi;
            }
            else
            {
                x(i) = 1;
            }
        }
        kase = 2;
        jump = 4;
        internalcomplexrcondsaveall(isave, rsave, i, iter, j, jlast, jump, absxi, altsgn, estold, temp);
        return;
    }
    
    //
    // ENTRY   (JUMP = 4)
    // X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
    //
    if( jump==4 )
    {
        jlast = j;
        j = internalcomplexrcondicmax1(x, n);
        if( ap::fp_neq(ap::abscomplex(x(jlast)),ap::abscomplex(x(j)))&&iter<itmax )
        {
            iter = iter+1;
            
            //
            // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
            //
            for(i = 1; i <= n; i++)
            {
                x(i) = 0;
            }
            x(j) = 1;
            kase = 1;
            jump = 3;
            internalcomplexrcondsaveall(isave, rsave, i, iter, j, jlast, jump, absxi, altsgn, estold, temp);
            return;
        }
        
        //
        // ITERATION COMPLETE.  FINAL STAGE.
        //
        altsgn = 1;
        for(i = 1; i <= n; i++)
        {
            x(i) = altsgn*(1+double(i-1)/double(n-1));
            altsgn = -altsgn;
        }
        kase = 1;
        jump = 5;
        internalcomplexrcondsaveall(isave, rsave, i, iter, j, jlast, jump, absxi, altsgn, estold, temp);
        return;
    }
    
    //
    // ENTRY   (JUMP = 5)
    // X HAS BEEN OVERWRITTEN BY A*X.
    //
    if( jump==5 )
    {
        temp = 2*(internalcomplexrcondscsum1(x, n)/(3*n));
        if( ap::fp_greater(temp,est) )
        {
            ap::vmove(&v(1), 1, &x(1), 1, "N", ap::vlen(1,n));
            est = temp;
        }
        kase = 0;
        internalcomplexrcondsaveall(isave, rsave, i, iter, j, jlast, jump, absxi, altsgn, estold, temp);
        return;
    }
}


static double internalcomplexrcondscsum1(const ap::complex_1d_array& x,
     int n)
{
    double result;
    int i;

    result = 0;
    for(i = 1; i <= n; i++)
    {
        result = result+ap::abscomplex(x(i));
    }
    return result;
}


static int internalcomplexrcondicmax1(const ap::complex_1d_array& x, int n)
{
    int result;
    int i;
    double m;

    result = 1;
    m = ap::abscomplex(x(1));
    for(i = 2; i <= n; i++)
    {
        if( ap::fp_greater(ap::abscomplex(x(i)),m) )
        {
            result = i;
            m = ap::abscomplex(x(i));
        }
    }
    return result;
}


static void internalcomplexrcondsaveall(ap::integer_1d_array& isave,
     ap::real_1d_array& rsave,
     int& i,
     int& iter,
     int& j,
     int& jlast,
     int& jump,
     double& absxi,
     double& altsgn,
     double& estold,
     double& temp)
{

    isave(0) = i;
    isave(1) = iter;
    isave(2) = j;
    isave(3) = jlast;
    isave(4) = jump;
    rsave(0) = absxi;
    rsave(1) = altsgn;
    rsave(2) = estold;
    rsave(3) = temp;
}


static void internalcomplexrcondloadall(ap::integer_1d_array& isave,
     ap::real_1d_array& rsave,
     int& i,
     int& iter,
     int& j,
     int& jlast,
     int& jump,
     double& absxi,
     double& altsgn,
     double& estold,
     double& temp)
{

    i = isave(0);
    iter = isave(1);
    j = isave(2);
    jlast = isave(3);
    jump = isave(4);
    absxi = rsave(0);
    altsgn = rsave(1);
    estold = rsave(2);
    temp = rsave(3);
}




