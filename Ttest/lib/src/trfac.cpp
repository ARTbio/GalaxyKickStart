/*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee. All rights reserved.

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
#include "trfac.h"

static void cmatrixluprec(ap::complex_2d_array& a,
     int offs,
     int m,
     int n,
     ap::integer_1d_array& pivots,
     ap::complex_1d_array& tmp);
static void rmatrixluprec(ap::real_2d_array& a,
     int offs,
     int m,
     int n,
     ap::integer_1d_array& pivots,
     ap::real_1d_array& tmp);
static void cmatrixplurec(ap::complex_2d_array& a,
     int offs,
     int m,
     int n,
     ap::integer_1d_array& pivots,
     ap::complex_1d_array& tmp);
static void rmatrixplurec(ap::real_2d_array& a,
     int offs,
     int m,
     int n,
     ap::integer_1d_array& pivots,
     ap::real_1d_array& tmp);
static void cmatrixlup2(ap::complex_2d_array& a,
     int offs,
     int m,
     int n,
     ap::integer_1d_array& pivots,
     ap::complex_1d_array& tmp);
static void rmatrixlup2(ap::real_2d_array& a,
     int offs,
     int m,
     int n,
     ap::integer_1d_array& pivots,
     ap::real_1d_array& tmp);
static void cmatrixplu2(ap::complex_2d_array& a,
     int offs,
     int m,
     int n,
     ap::integer_1d_array& pivots,
     ap::complex_1d_array& tmp);
static void rmatrixplu2(ap::real_2d_array& a,
     int offs,
     int m,
     int n,
     ap::integer_1d_array& pivots,
     ap::real_1d_array& tmp);
static bool hpdmatrixcholeskyrec(ap::complex_2d_array& a,
     int offs,
     int n,
     bool isupper,
     ap::complex_1d_array& tmp);
static bool spdmatrixcholeskyrec(ap::real_2d_array& a,
     int offs,
     int n,
     bool isupper,
     ap::real_1d_array& tmp);
static bool hpdmatrixcholesky2(ap::complex_2d_array& aaa,
     int offs,
     int n,
     bool isupper,
     ap::complex_1d_array& tmp);
static bool spdmatrixcholesky2(ap::real_2d_array& aaa,
     int offs,
     int n,
     bool isupper,
     ap::real_1d_array& tmp);

/*************************************************************************
LU decomposition of a general real matrix with row pivoting

A is represented as A = P*L*U, where:
* L is lower unitriangular matrix
* U is upper triangular matrix
* P = P0*P1*...*PK, K=min(M,N)-1,
  Pi - permutation matrix for I and Pivots[I]

This is cache-oblivous implementation of LU decomposition.
It is optimized for square matrices. As for rectangular matrices:
* best case - M>>N
* worst case - N>>M, small M, large N, matrix does not fit in CPU cache

INPUT PARAMETERS:
    A       -   array[0..M-1, 0..N-1].
    M       -   number of rows in matrix A.
    N       -   number of columns in matrix A.


OUTPUT PARAMETERS:
    A       -   matrices L and U in compact form:
                * L is stored under main diagonal
                * U is stored on and above main diagonal
    Pivots  -   permutation matrix in compact form.
                array[0..Min(M-1,N-1)].

  -- ALGLIB routine --
     10.01.2010
     Bochkanov Sergey
*************************************************************************/
void rmatrixlu(ap::real_2d_array& a,
     int m,
     int n,
     ap::integer_1d_array& pivots)
{

    ap::ap_error::make_assertion(m>0, "RMatrixLU: incorrect M!");
    ap::ap_error::make_assertion(n>0, "RMatrixLU: incorrect N!");
    rmatrixplu(a, m, n, pivots);
}


/*************************************************************************
LU decomposition of a general complex matrix with row pivoting

A is represented as A = P*L*U, where:
* L is lower unitriangular matrix
* U is upper triangular matrix
* P = P0*P1*...*PK, K=min(M,N)-1,
  Pi - permutation matrix for I and Pivots[I]

This is cache-oblivous implementation of LU decomposition. It is optimized
for square matrices. As for rectangular matrices:
* best case - M>>N
* worst case - N>>M, small M, large N, matrix does not fit in CPU cache

INPUT PARAMETERS:
    A       -   array[0..M-1, 0..N-1].
    M       -   number of rows in matrix A.
    N       -   number of columns in matrix A.


OUTPUT PARAMETERS:
    A       -   matrices L and U in compact form:
                * L is stored under main diagonal
                * U is stored on and above main diagonal
    Pivots  -   permutation matrix in compact form.
                array[0..Min(M-1,N-1)].

  -- ALGLIB routine --
     10.01.2010
     Bochkanov Sergey
*************************************************************************/
void cmatrixlu(ap::complex_2d_array& a,
     int m,
     int n,
     ap::integer_1d_array& pivots)
{

    ap::ap_error::make_assertion(m>0, "CMatrixLU: incorrect M!");
    ap::ap_error::make_assertion(n>0, "CMatrixLU: incorrect N!");
    cmatrixplu(a, m, n, pivots);
}


/*************************************************************************
Cache-oblivious Cholesky decomposition

The algorithm computes Cholesky decomposition  of  a  Hermitian  positive-
definite matrix. The result of an algorithm is a representation  of  A  as
A=U'*U  or A=L*L' (here X' detones conj(X^T)).

INPUT PARAMETERS:
    A       -   upper or lower triangle of a factorized matrix.
                array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   if IsUpper=True, then A contains an upper triangle of
                a symmetric matrix, otherwise A contains a lower one.

OUTPUT PARAMETERS:
    A       -   the result of factorization. If IsUpper=True, then
                the upper triangle contains matrix U, so that A = U'*U,
                and the elements below the main diagonal are not modified.
                Similarly, if IsUpper = False.

RESULT:
    If  the  matrix  is  positive-definite,  the  function  returns  True.
    Otherwise, the function returns False. Contents of A is not determined
    in such case.

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
bool hpdmatrixcholesky(ap::complex_2d_array& a, int n, bool isupper)
{
    bool result;
    ap::complex_1d_array tmp;

    if( n<1 )
    {
        result = false;
        return result;
    }
    tmp.setlength(2*n);
    result = hpdmatrixcholeskyrec(a, 0, n, isupper, tmp);
    return result;
}


/*************************************************************************
Cache-oblivious Cholesky decomposition

The algorithm computes Cholesky decomposition  of  a  symmetric  positive-
definite matrix. The result of an algorithm is a representation  of  A  as
A=U^T*U  or A=L*L^T

INPUT PARAMETERS:
    A       -   upper or lower triangle of a factorized matrix.
                array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   if IsUpper=True, then A contains an upper triangle of
                a symmetric matrix, otherwise A contains a lower one.

OUTPUT PARAMETERS:
    A       -   the result of factorization. If IsUpper=True, then
                the upper triangle contains matrix U, so that A = U^T*U,
                and the elements below the main diagonal are not modified.
                Similarly, if IsUpper = False.

RESULT:
    If  the  matrix  is  positive-definite,  the  function  returns  True.
    Otherwise, the function returns False. Contents of A is not determined
    in such case.

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
bool spdmatrixcholesky(ap::real_2d_array& a, int n, bool isupper)
{
    bool result;
    ap::real_1d_array tmp;

    if( n<1 )
    {
        result = false;
        return result;
    }
    tmp.setlength(2*n);
    result = spdmatrixcholeskyrec(a, 0, n, isupper, tmp);
    return result;
}


void rmatrixlup(ap::real_2d_array& a,
     int m,
     int n,
     ap::integer_1d_array& pivots)
{
    ap::real_1d_array tmp;
    int i;
    int j;
    double mx;
    double v;

    
    //
    // Internal LU decomposition subroutine.
    // Never call it directly.
    //
    ap::ap_error::make_assertion(m>0, "RMatrixLUP: incorrect M!");
    ap::ap_error::make_assertion(n>0, "RMatrixLUP: incorrect N!");
    
    //
    // Scale matrix to avoid overflows,
    // decompose it, then scale back.
    //
    mx = 0;
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            mx = ap::maxreal(mx, fabs(a(i,j)));
        }
    }
    if( ap::fp_neq(mx,0) )
    {
        v = 1/mx;
        for(i = 0; i <= m-1; i++)
        {
            ap::vmul(&a(i, 0), 1, ap::vlen(0,n-1), v);
        }
    }
    pivots.setlength(ap::minint(m, n));
    tmp.setlength(2*ap::maxint(m, n));
    rmatrixluprec(a, 0, m, n, pivots, tmp);
    if( ap::fp_neq(mx,0) )
    {
        v = mx;
        for(i = 0; i <= m-1; i++)
        {
            ap::vmul(&a(i, 0), 1, ap::vlen(0,ap::minint(i, n-1)), v);
        }
    }
}


void cmatrixlup(ap::complex_2d_array& a,
     int m,
     int n,
     ap::integer_1d_array& pivots)
{
    ap::complex_1d_array tmp;
    int i;
    int j;
    double mx;
    double v;

    
    //
    // Internal LU decomposition subroutine.
    // Never call it directly.
    //
    ap::ap_error::make_assertion(m>0, "CMatrixLUP: incorrect M!");
    ap::ap_error::make_assertion(n>0, "CMatrixLUP: incorrect N!");
    
    //
    // Scale matrix to avoid overflows,
    // decompose it, then scale back.
    //
    mx = 0;
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            mx = ap::maxreal(mx, ap::abscomplex(a(i,j)));
        }
    }
    if( ap::fp_neq(mx,0) )
    {
        v = 1/mx;
        for(i = 0; i <= m-1; i++)
        {
            ap::vmul(&a(i, 0), 1, ap::vlen(0,n-1), v);
        }
    }
    pivots.setlength(ap::minint(m, n));
    tmp.setlength(2*ap::maxint(m, n));
    cmatrixluprec(a, 0, m, n, pivots, tmp);
    if( ap::fp_neq(mx,0) )
    {
        v = mx;
        for(i = 0; i <= m-1; i++)
        {
            ap::vmul(&a(i, 0), 1, ap::vlen(0,ap::minint(i, n-1)), v);
        }
    }
}


void rmatrixplu(ap::real_2d_array& a,
     int m,
     int n,
     ap::integer_1d_array& pivots)
{
    ap::real_1d_array tmp;
    int i;
    int j;
    double mx;
    double v;

    
    //
    // Internal LU decomposition subroutine.
    // Never call it directly.
    //
    ap::ap_error::make_assertion(m>0, "RMatrixPLU: incorrect M!");
    ap::ap_error::make_assertion(n>0, "RMatrixPLU: incorrect N!");
    tmp.setlength(2*ap::maxint(m, n));
    pivots.setlength(ap::minint(m, n));
    
    //
    // Scale matrix to avoid overflows,
    // decompose it, then scale back.
    //
    mx = 0;
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            mx = ap::maxreal(mx, fabs(a(i,j)));
        }
    }
    if( ap::fp_neq(mx,0) )
    {
        v = 1/mx;
        for(i = 0; i <= m-1; i++)
        {
            ap::vmul(&a(i, 0), 1, ap::vlen(0,n-1), v);
        }
    }
    rmatrixplurec(a, 0, m, n, pivots, tmp);
    if( ap::fp_neq(mx,0) )
    {
        v = mx;
        for(i = 0; i <= ap::minint(m, n)-1; i++)
        {
            ap::vmul(&a(i, i), 1, ap::vlen(i,n-1), v);
        }
    }
}


void cmatrixplu(ap::complex_2d_array& a,
     int m,
     int n,
     ap::integer_1d_array& pivots)
{
    ap::complex_1d_array tmp;
    int i;
    int j;
    double mx;
    ap::complex v;

    
    //
    // Internal LU decomposition subroutine.
    // Never call it directly.
    //
    ap::ap_error::make_assertion(m>0, "CMatrixPLU: incorrect M!");
    ap::ap_error::make_assertion(n>0, "CMatrixPLU: incorrect N!");
    tmp.setlength(2*ap::maxint(m, n));
    pivots.setlength(ap::minint(m, n));
    
    //
    // Scale matrix to avoid overflows,
    // decompose it, then scale back.
    //
    mx = 0;
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            mx = ap::maxreal(mx, ap::abscomplex(a(i,j)));
        }
    }
    if( ap::fp_neq(mx,0) )
    {
        v = 1/mx;
        for(i = 0; i <= m-1; i++)
        {
            ap::vmul(&a(i, 0), 1, ap::vlen(0,n-1), v);
        }
    }
    cmatrixplurec(a, 0, m, n, pivots, tmp);
    if( ap::fp_neq(mx,0) )
    {
        v = mx;
        for(i = 0; i <= ap::minint(m, n)-1; i++)
        {
            ap::vmul(&a(i, i), 1, ap::vlen(i,n-1), v);
        }
    }
}


/*************************************************************************
Recurrent complex LU subroutine.
Never call it directly.

  -- ALGLIB routine --
     04.01.2010
     Bochkanov Sergey
*************************************************************************/
static void cmatrixluprec(ap::complex_2d_array& a,
     int offs,
     int m,
     int n,
     ap::integer_1d_array& pivots,
     ap::complex_1d_array& tmp)
{
    int i;
    int m1;
    int m2;

    
    //
    // Kernel case
    //
    if( ap::minint(m, n)<=ablascomplexblocksize(a) )
    {
        cmatrixlup2(a, offs, m, n, pivots, tmp);
        return;
    }
    
    //
    // Preliminary step, make N>=M
    //
    //     ( A1 )
    // A = (    ), where A1 is square
    //     ( A2 )
    //
    // Factorize A1, update A2
    //
    if( m>n )
    {
        cmatrixluprec(a, offs, n, n, pivots, tmp);
        for(i = 0; i <= n-1; i++)
        {
            ap::vmove(&tmp(0), 1, &a(offs+n, offs+i), a.getstride(), "N", ap::vlen(0,m-n-1));
            ap::vmove(&a(offs+n, offs+i), a.getstride(), &a(offs+n, pivots(offs+i)), a.getstride(), "N", ap::vlen(offs+n,offs+m-1));
            ap::vmove(&a(offs+n, pivots(offs+i)), a.getstride(), &tmp(0), 1, "N", ap::vlen(offs+n,offs+m-1));
        }
        cmatrixrighttrsm(m-n, n, a, offs, offs, true, true, 0, a, offs+n, offs);
        return;
    }
    
    //
    // Non-kernel case
    //
    ablascomplexsplitlength(a, m, m1, m2);
    cmatrixluprec(a, offs, m1, n, pivots, tmp);
    if( m2>0 )
    {
        for(i = 0; i <= m1-1; i++)
        {
            if( offs+i!=pivots(offs+i) )
            {
                ap::vmove(&tmp(0), 1, &a(offs+m1, offs+i), a.getstride(), "N", ap::vlen(0,m2-1));
                ap::vmove(&a(offs+m1, offs+i), a.getstride(), &a(offs+m1, pivots(offs+i)), a.getstride(), "N", ap::vlen(offs+m1,offs+m-1));
                ap::vmove(&a(offs+m1, pivots(offs+i)), a.getstride(), &tmp(0), 1, "N", ap::vlen(offs+m1,offs+m-1));
            }
        }
        cmatrixrighttrsm(m2, m1, a, offs, offs, true, true, 0, a, offs+m1, offs);
        cmatrixgemm(m-m1, n-m1, m1, -1.0, a, offs+m1, offs, 0, a, offs, offs+m1, 0, +1.0, a, offs+m1, offs+m1);
        cmatrixluprec(a, offs+m1, m-m1, n-m1, pivots, tmp);
        for(i = 0; i <= m2-1; i++)
        {
            if( offs+m1+i!=pivots(offs+m1+i) )
            {
                ap::vmove(&tmp(0), 1, &a(offs, offs+m1+i), a.getstride(), "N", ap::vlen(0,m1-1));
                ap::vmove(&a(offs, offs+m1+i), a.getstride(), &a(offs, pivots(offs+m1+i)), a.getstride(), "N", ap::vlen(offs,offs+m1-1));
                ap::vmove(&a(offs, pivots(offs+m1+i)), a.getstride(), &tmp(0), 1, "N", ap::vlen(offs,offs+m1-1));
            }
        }
    }
}


/*************************************************************************
Recurrent real LU subroutine.
Never call it directly.

  -- ALGLIB routine --
     04.01.2010
     Bochkanov Sergey
*************************************************************************/
static void rmatrixluprec(ap::real_2d_array& a,
     int offs,
     int m,
     int n,
     ap::integer_1d_array& pivots,
     ap::real_1d_array& tmp)
{
    int i;
    int m1;
    int m2;

    
    //
    // Kernel case
    //
    if( ap::minint(m, n)<=ablasblocksize(a) )
    {
        rmatrixlup2(a, offs, m, n, pivots, tmp);
        return;
    }
    
    //
    // Preliminary step, make N>=M
    //
    //     ( A1 )
    // A = (    ), where A1 is square
    //     ( A2 )
    //
    // Factorize A1, update A2
    //
    if( m>n )
    {
        rmatrixluprec(a, offs, n, n, pivots, tmp);
        for(i = 0; i <= n-1; i++)
        {
            if( offs+i!=pivots(offs+i) )
            {
                ap::vmove(&tmp(0), 1, &a(offs+n, offs+i), a.getstride(), ap::vlen(0,m-n-1));
                ap::vmove(&a(offs+n, offs+i), a.getstride(), &a(offs+n, pivots(offs+i)), a.getstride(), ap::vlen(offs+n,offs+m-1));
                ap::vmove(&a(offs+n, pivots(offs+i)), a.getstride(), &tmp(0), 1, ap::vlen(offs+n,offs+m-1));
            }
        }
        rmatrixrighttrsm(m-n, n, a, offs, offs, true, true, 0, a, offs+n, offs);
        return;
    }
    
    //
    // Non-kernel case
    //
    ablassplitlength(a, m, m1, m2);
    rmatrixluprec(a, offs, m1, n, pivots, tmp);
    if( m2>0 )
    {
        for(i = 0; i <= m1-1; i++)
        {
            if( offs+i!=pivots(offs+i) )
            {
                ap::vmove(&tmp(0), 1, &a(offs+m1, offs+i), a.getstride(), ap::vlen(0,m2-1));
                ap::vmove(&a(offs+m1, offs+i), a.getstride(), &a(offs+m1, pivots(offs+i)), a.getstride(), ap::vlen(offs+m1,offs+m-1));
                ap::vmove(&a(offs+m1, pivots(offs+i)), a.getstride(), &tmp(0), 1, ap::vlen(offs+m1,offs+m-1));
            }
        }
        rmatrixrighttrsm(m2, m1, a, offs, offs, true, true, 0, a, offs+m1, offs);
        rmatrixgemm(m-m1, n-m1, m1, -1.0, a, offs+m1, offs, 0, a, offs, offs+m1, 0, +1.0, a, offs+m1, offs+m1);
        rmatrixluprec(a, offs+m1, m-m1, n-m1, pivots, tmp);
        for(i = 0; i <= m2-1; i++)
        {
            if( offs+m1+i!=pivots(offs+m1+i) )
            {
                ap::vmove(&tmp(0), 1, &a(offs, offs+m1+i), a.getstride(), ap::vlen(0,m1-1));
                ap::vmove(&a(offs, offs+m1+i), a.getstride(), &a(offs, pivots(offs+m1+i)), a.getstride(), ap::vlen(offs,offs+m1-1));
                ap::vmove(&a(offs, pivots(offs+m1+i)), a.getstride(), &tmp(0), 1, ap::vlen(offs,offs+m1-1));
            }
        }
    }
}


/*************************************************************************
Recurrent complex LU subroutine.
Never call it directly.

  -- ALGLIB routine --
     04.01.2010
     Bochkanov Sergey
*************************************************************************/
static void cmatrixplurec(ap::complex_2d_array& a,
     int offs,
     int m,
     int n,
     ap::integer_1d_array& pivots,
     ap::complex_1d_array& tmp)
{
    int i;
    int n1;
    int n2;

    
    //
    // Kernel case
    //
    if( ap::minint(m, n)<=ablascomplexblocksize(a) )
    {
        cmatrixplu2(a, offs, m, n, pivots, tmp);
        return;
    }
    
    //
    // Preliminary step, make M>=N.
    //
    // A = (A1 A2), where A1 is square
    // Factorize A1, update A2
    //
    if( n>m )
    {
        cmatrixplurec(a, offs, m, m, pivots, tmp);
        for(i = 0; i <= m-1; i++)
        {
            ap::vmove(&tmp(0), 1, &a(offs+i, offs+m), 1, "N", ap::vlen(0,n-m-1));
            ap::vmove(&a(offs+i, offs+m), 1, &a(pivots(offs+i), offs+m), 1, "N", ap::vlen(offs+m,offs+n-1));
            ap::vmove(&a(pivots(offs+i), offs+m), 1, &tmp(0), 1, "N", ap::vlen(offs+m,offs+n-1));
        }
        cmatrixlefttrsm(m, n-m, a, offs, offs, false, true, 0, a, offs, offs+m);
        return;
    }
    
    //
    // Non-kernel case
    //
    ablascomplexsplitlength(a, n, n1, n2);
    cmatrixplurec(a, offs, m, n1, pivots, tmp);
    if( n2>0 )
    {
        for(i = 0; i <= n1-1; i++)
        {
            if( offs+i!=pivots(offs+i) )
            {
                ap::vmove(&tmp(0), 1, &a(offs+i, offs+n1), 1, "N", ap::vlen(0,n2-1));
                ap::vmove(&a(offs+i, offs+n1), 1, &a(pivots(offs+i), offs+n1), 1, "N", ap::vlen(offs+n1,offs+n-1));
                ap::vmove(&a(pivots(offs+i), offs+n1), 1, &tmp(0), 1, "N", ap::vlen(offs+n1,offs+n-1));
            }
        }
        cmatrixlefttrsm(n1, n2, a, offs, offs, false, true, 0, a, offs, offs+n1);
        cmatrixgemm(m-n1, n-n1, n1, -1.0, a, offs+n1, offs, 0, a, offs, offs+n1, 0, +1.0, a, offs+n1, offs+n1);
        cmatrixplurec(a, offs+n1, m-n1, n-n1, pivots, tmp);
        for(i = 0; i <= n2-1; i++)
        {
            if( offs+n1+i!=pivots(offs+n1+i) )
            {
                ap::vmove(&tmp(0), 1, &a(offs+n1+i, offs), 1, "N", ap::vlen(0,n1-1));
                ap::vmove(&a(offs+n1+i, offs), 1, &a(pivots(offs+n1+i), offs), 1, "N", ap::vlen(offs,offs+n1-1));
                ap::vmove(&a(pivots(offs+n1+i), offs), 1, &tmp(0), 1, "N", ap::vlen(offs,offs+n1-1));
            }
        }
    }
}


/*************************************************************************
Recurrent real LU subroutine.
Never call it directly.

  -- ALGLIB routine --
     04.01.2010
     Bochkanov Sergey
*************************************************************************/
static void rmatrixplurec(ap::real_2d_array& a,
     int offs,
     int m,
     int n,
     ap::integer_1d_array& pivots,
     ap::real_1d_array& tmp)
{
    int i;
    int n1;
    int n2;

    
    //
    // Kernel case
    //
    if( ap::minint(m, n)<=ablasblocksize(a) )
    {
        rmatrixplu2(a, offs, m, n, pivots, tmp);
        return;
    }
    
    //
    // Preliminary step, make M>=N.
    //
    // A = (A1 A2), where A1 is square
    // Factorize A1, update A2
    //
    if( n>m )
    {
        rmatrixplurec(a, offs, m, m, pivots, tmp);
        for(i = 0; i <= m-1; i++)
        {
            ap::vmove(&tmp(0), 1, &a(offs+i, offs+m), 1, ap::vlen(0,n-m-1));
            ap::vmove(&a(offs+i, offs+m), 1, &a(pivots(offs+i), offs+m), 1, ap::vlen(offs+m,offs+n-1));
            ap::vmove(&a(pivots(offs+i), offs+m), 1, &tmp(0), 1, ap::vlen(offs+m,offs+n-1));
        }
        rmatrixlefttrsm(m, n-m, a, offs, offs, false, true, 0, a, offs, offs+m);
        return;
    }
    
    //
    // Non-kernel case
    //
    ablassplitlength(a, n, n1, n2);
    rmatrixplurec(a, offs, m, n1, pivots, tmp);
    if( n2>0 )
    {
        for(i = 0; i <= n1-1; i++)
        {
            if( offs+i!=pivots(offs+i) )
            {
                ap::vmove(&tmp(0), 1, &a(offs+i, offs+n1), 1, ap::vlen(0,n2-1));
                ap::vmove(&a(offs+i, offs+n1), 1, &a(pivots(offs+i), offs+n1), 1, ap::vlen(offs+n1,offs+n-1));
                ap::vmove(&a(pivots(offs+i), offs+n1), 1, &tmp(0), 1, ap::vlen(offs+n1,offs+n-1));
            }
        }
        rmatrixlefttrsm(n1, n2, a, offs, offs, false, true, 0, a, offs, offs+n1);
        rmatrixgemm(m-n1, n-n1, n1, -1.0, a, offs+n1, offs, 0, a, offs, offs+n1, 0, +1.0, a, offs+n1, offs+n1);
        rmatrixplurec(a, offs+n1, m-n1, n-n1, pivots, tmp);
        for(i = 0; i <= n2-1; i++)
        {
            if( offs+n1+i!=pivots(offs+n1+i) )
            {
                ap::vmove(&tmp(0), 1, &a(offs+n1+i, offs), 1, ap::vlen(0,n1-1));
                ap::vmove(&a(offs+n1+i, offs), 1, &a(pivots(offs+n1+i), offs), 1, ap::vlen(offs,offs+n1-1));
                ap::vmove(&a(pivots(offs+n1+i), offs), 1, &tmp(0), 1, ap::vlen(offs,offs+n1-1));
            }
        }
    }
}


/*************************************************************************
Complex LUP kernel

  -- ALGLIB routine --
     10.01.2010
     Bochkanov Sergey
*************************************************************************/
static void cmatrixlup2(ap::complex_2d_array& a,
     int offs,
     int m,
     int n,
     ap::integer_1d_array& pivots,
     ap::complex_1d_array& tmp)
{
    int i;
    int j;
    int jp;
    ap::complex s;

    
    //
    // Quick return if possible
    //
    if( m==0||n==0 )
    {
        return;
    }
    
    //
    // main cycle
    //
    for(j = 0; j <= ap::minint(m-1, n-1); j++)
    {
        
        //
        // Find pivot, swap columns
        //
        jp = j;
        for(i = j+1; i <= n-1; i++)
        {
            if( ap::fp_greater(ap::abscomplex(a(offs+j,offs+i)),ap::abscomplex(a(offs+j,offs+jp))) )
            {
                jp = i;
            }
        }
        pivots(offs+j) = offs+jp;
        if( jp!=j )
        {
            ap::vmove(&tmp(0), 1, &a(offs, offs+j), a.getstride(), "N", ap::vlen(0,m-1));
            ap::vmove(&a(offs, offs+j), a.getstride(), &a(offs, offs+jp), a.getstride(), "N", ap::vlen(offs,offs+m-1));
            ap::vmove(&a(offs, offs+jp), a.getstride(), &tmp(0), 1, "N", ap::vlen(offs,offs+m-1));
        }
        
        //
        // LU decomposition of 1x(N-J) matrix
        //
        if( a(offs+j,offs+j)!=0&&j+1<=n-1 )
        {
            s = 1/a(offs+j,offs+j);
            ap::vmul(&a(offs+j, offs+j+1), 1, ap::vlen(offs+j+1,offs+n-1), s);
        }
        
        //
        // Update trailing (M-J-1)x(N-J-1) matrix
        //
        if( j<ap::minint(m-1, n-1) )
        {
            ap::vmove(&tmp(0), 1, &a(offs+j+1, offs+j), a.getstride(), "N", ap::vlen(0,m-j-2));
            ap::vmoveneg(&tmp(m), 1, &a(offs+j, offs+j+1), 1, "N", ap::vlen(m,m+n-j-2));
            cmatrixrank1(m-j-1, n-j-1, a, offs+j+1, offs+j+1, tmp, 0, tmp, m);
        }
    }
}


/*************************************************************************
Real LUP kernel

  -- ALGLIB routine --
     10.01.2010
     Bochkanov Sergey
*************************************************************************/
static void rmatrixlup2(ap::real_2d_array& a,
     int offs,
     int m,
     int n,
     ap::integer_1d_array& pivots,
     ap::real_1d_array& tmp)
{
    int i;
    int j;
    int jp;
    double s;

    
    //
    // Quick return if possible
    //
    if( m==0||n==0 )
    {
        return;
    }
    
    //
    // main cycle
    //
    for(j = 0; j <= ap::minint(m-1, n-1); j++)
    {
        
        //
        // Find pivot, swap columns
        //
        jp = j;
        for(i = j+1; i <= n-1; i++)
        {
            if( ap::fp_greater(fabs(a(offs+j,offs+i)),fabs(a(offs+j,offs+jp))) )
            {
                jp = i;
            }
        }
        pivots(offs+j) = offs+jp;
        if( jp!=j )
        {
            ap::vmove(&tmp(0), 1, &a(offs, offs+j), a.getstride(), ap::vlen(0,m-1));
            ap::vmove(&a(offs, offs+j), a.getstride(), &a(offs, offs+jp), a.getstride(), ap::vlen(offs,offs+m-1));
            ap::vmove(&a(offs, offs+jp), a.getstride(), &tmp(0), 1, ap::vlen(offs,offs+m-1));
        }
        
        //
        // LU decomposition of 1x(N-J) matrix
        //
        if( ap::fp_neq(a(offs+j,offs+j),0)&&j+1<=n-1 )
        {
            s = 1/a(offs+j,offs+j);
            ap::vmul(&a(offs+j, offs+j+1), 1, ap::vlen(offs+j+1,offs+n-1), s);
        }
        
        //
        // Update trailing (M-J-1)x(N-J-1) matrix
        //
        if( j<ap::minint(m-1, n-1) )
        {
            ap::vmove(&tmp(0), 1, &a(offs+j+1, offs+j), a.getstride(), ap::vlen(0,m-j-2));
            ap::vmoveneg(&tmp(m), 1, &a(offs+j, offs+j+1), 1, ap::vlen(m,m+n-j-2));
            rmatrixrank1(m-j-1, n-j-1, a, offs+j+1, offs+j+1, tmp, 0, tmp, m);
        }
    }
}


/*************************************************************************
Complex PLU kernel

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     June 30, 1992
*************************************************************************/
static void cmatrixplu2(ap::complex_2d_array& a,
     int offs,
     int m,
     int n,
     ap::integer_1d_array& pivots,
     ap::complex_1d_array& tmp)
{
    int i;
    int j;
    int jp;
    ap::complex s;

    
    //
    // Quick return if possible
    //
    if( m==0||n==0 )
    {
        return;
    }
    for(j = 0; j <= ap::minint(m-1, n-1); j++)
    {
        
        //
        // Find pivot and test for singularity.
        //
        jp = j;
        for(i = j+1; i <= m-1; i++)
        {
            if( ap::fp_greater(ap::abscomplex(a(offs+i,offs+j)),ap::abscomplex(a(offs+jp,offs+j))) )
            {
                jp = i;
            }
        }
        pivots(offs+j) = offs+jp;
        if( a(offs+jp,offs+j)!=0 )
        {
            
            //
            //Apply the interchange to rows
            //
            if( jp!=j )
            {
                for(i = 0; i <= n-1; i++)
                {
                    s = a(offs+j,offs+i);
                    a(offs+j,offs+i) = a(offs+jp,offs+i);
                    a(offs+jp,offs+i) = s;
                }
            }
            
            //
            //Compute elements J+1:M of J-th column.
            //
            if( j+1<=m-1 )
            {
                s = 1/a(offs+j,offs+j);
                ap::vmul(&a(offs+j+1, offs+j), a.getstride(), ap::vlen(offs+j+1,offs+m-1), s);
            }
        }
        if( j<ap::minint(m, n)-1 )
        {
            
            //
            //Update trailing submatrix.
            //
            ap::vmove(&tmp(0), 1, &a(offs+j+1, offs+j), a.getstride(), "N", ap::vlen(0,m-j-2));
            ap::vmoveneg(&tmp(m), 1, &a(offs+j, offs+j+1), 1, "N", ap::vlen(m,m+n-j-2));
            cmatrixrank1(m-j-1, n-j-1, a, offs+j+1, offs+j+1, tmp, 0, tmp, m);
        }
    }
}


/*************************************************************************
Real PLU kernel

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     June 30, 1992
*************************************************************************/
static void rmatrixplu2(ap::real_2d_array& a,
     int offs,
     int m,
     int n,
     ap::integer_1d_array& pivots,
     ap::real_1d_array& tmp)
{
    int i;
    int j;
    int jp;
    double s;

    
    //
    // Quick return if possible
    //
    if( m==0||n==0 )
    {
        return;
    }
    for(j = 0; j <= ap::minint(m-1, n-1); j++)
    {
        
        //
        // Find pivot and test for singularity.
        //
        jp = j;
        for(i = j+1; i <= m-1; i++)
        {
            if( ap::fp_greater(fabs(a(offs+i,offs+j)),fabs(a(offs+jp,offs+j))) )
            {
                jp = i;
            }
        }
        pivots(offs+j) = offs+jp;
        if( ap::fp_neq(a(offs+jp,offs+j),0) )
        {
            
            //
            //Apply the interchange to rows
            //
            if( jp!=j )
            {
                for(i = 0; i <= n-1; i++)
                {
                    s = a(offs+j,offs+i);
                    a(offs+j,offs+i) = a(offs+jp,offs+i);
                    a(offs+jp,offs+i) = s;
                }
            }
            
            //
            //Compute elements J+1:M of J-th column.
            //
            if( j+1<=m-1 )
            {
                s = 1/a(offs+j,offs+j);
                ap::vmul(&a(offs+j+1, offs+j), a.getstride(), ap::vlen(offs+j+1,offs+m-1), s);
            }
        }
        if( j<ap::minint(m, n)-1 )
        {
            
            //
            //Update trailing submatrix.
            //
            ap::vmove(&tmp(0), 1, &a(offs+j+1, offs+j), a.getstride(), ap::vlen(0,m-j-2));
            ap::vmoveneg(&tmp(m), 1, &a(offs+j, offs+j+1), 1, ap::vlen(m,m+n-j-2));
            rmatrixrank1(m-j-1, n-j-1, a, offs+j+1, offs+j+1, tmp, 0, tmp, m);
        }
    }
}


/*************************************************************************
Recursive computational subroutine for HPDMatrixCholesky

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
static bool hpdmatrixcholeskyrec(ap::complex_2d_array& a,
     int offs,
     int n,
     bool isupper,
     ap::complex_1d_array& tmp)
{
    bool result;
    int n1;
    int n2;

    
    //
    // check N
    //
    if( n<1 )
    {
        result = false;
        return result;
    }
    
    //
    // special cases
    //
    if( n==1 )
    {
        if( ap::fp_greater(a(offs,offs).x,0) )
        {
            a(offs,offs) = sqrt(a(offs,offs).x);
            result = true;
        }
        else
        {
            result = false;
        }
        return result;
    }
    if( n<=ablascomplexblocksize(a) )
    {
        result = hpdmatrixcholesky2(a, offs, n, isupper, tmp);
        return result;
    }
    
    //
    // general case: split task in cache-oblivious manner
    //
    result = true;
    ablascomplexsplitlength(a, n, n1, n2);
    result = hpdmatrixcholeskyrec(a, offs, n1, isupper, tmp);
    if( !result )
    {
        return result;
    }
    if( n2>0 )
    {
        if( isupper )
        {
            cmatrixlefttrsm(n1, n2, a, offs, offs, isupper, false, 2, a, offs, offs+n1);
            cmatrixsyrk(n2, n1, -1.0, a, offs, offs+n1, 2, +1.0, a, offs+n1, offs+n1, isupper);
        }
        else
        {
            cmatrixrighttrsm(n2, n1, a, offs, offs, isupper, false, 2, a, offs+n1, offs);
            cmatrixsyrk(n2, n1, -1.0, a, offs+n1, offs, 0, +1.0, a, offs+n1, offs+n1, isupper);
        }
        result = hpdmatrixcholeskyrec(a, offs+n1, n2, isupper, tmp);
        if( !result )
        {
            return result;
        }
    }
    return result;
}


/*************************************************************************
Recursive computational subroutine for SPDMatrixCholesky

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
static bool spdmatrixcholeskyrec(ap::real_2d_array& a,
     int offs,
     int n,
     bool isupper,
     ap::real_1d_array& tmp)
{
    bool result;
    int n1;
    int n2;

    
    //
    // check N
    //
    if( n<1 )
    {
        result = false;
        return result;
    }
    
    //
    // special cases
    //
    if( n==1 )
    {
        if( ap::fp_greater(a(offs,offs),0) )
        {
            a(offs,offs) = sqrt(a(offs,offs));
            result = true;
        }
        else
        {
            result = false;
        }
        return result;
    }
    if( n<=ablasblocksize(a) )
    {
        result = spdmatrixcholesky2(a, offs, n, isupper, tmp);
        return result;
    }
    
    //
    // general case: split task in cache-oblivious manner
    //
    result = true;
    ablassplitlength(a, n, n1, n2);
    result = spdmatrixcholeskyrec(a, offs, n1, isupper, tmp);
    if( !result )
    {
        return result;
    }
    if( n2>0 )
    {
        if( isupper )
        {
            rmatrixlefttrsm(n1, n2, a, offs, offs, isupper, false, 1, a, offs, offs+n1);
            rmatrixsyrk(n2, n1, -1.0, a, offs, offs+n1, 1, +1.0, a, offs+n1, offs+n1, isupper);
        }
        else
        {
            rmatrixrighttrsm(n2, n1, a, offs, offs, isupper, false, 1, a, offs+n1, offs);
            rmatrixsyrk(n2, n1, -1.0, a, offs+n1, offs, 0, +1.0, a, offs+n1, offs+n1, isupper);
        }
        result = spdmatrixcholeskyrec(a, offs+n1, n2, isupper, tmp);
        if( !result )
        {
            return result;
        }
    }
    return result;
}


/*************************************************************************
Level-2 Hermitian Cholesky subroutine.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************/
static bool hpdmatrixcholesky2(ap::complex_2d_array& aaa,
     int offs,
     int n,
     bool isupper,
     ap::complex_1d_array& tmp)
{
    bool result;
    int i;
    int j;
    int k;
    int j1;
    int j2;
    double ajj;
    ap::complex v;
    double r;

    result = true;
    if( n<0 )
    {
        result = false;
        return result;
    }
    
    //
    // Quick return if possible
    //
    if( n==0 )
    {
        return result;
    }
    if( isupper )
    {
        
        //
        // Compute the Cholesky factorization A = U'*U.
        //
        for(j = 0; j <= n-1; j++)
        {
            
            //
            // Compute U(J,J) and test for non-positive-definiteness.
            //
            v = ap::vdotproduct(&aaa(offs, offs+j), aaa.getstride(), "Conj", &aaa(offs, offs+j), aaa.getstride(), "N", ap::vlen(offs,offs+j-1));
            ajj = (aaa(offs+j,offs+j)-v).x;
            if( ap::fp_less_eq(ajj,0) )
            {
                aaa(offs+j,offs+j) = ajj;
                result = false;
                return result;
            }
            ajj = sqrt(ajj);
            aaa(offs+j,offs+j) = ajj;
            
            //
            // Compute elements J+1:N-1 of row J.
            //
            if( j<n-1 )
            {
                if( j>0 )
                {
                    ap::vmoveneg(&tmp(0), 1, &aaa(offs, offs+j), aaa.getstride(), "Conj", ap::vlen(0,j-1));
                    cmatrixmv(n-j-1, j, aaa, offs, offs+j+1, 1, tmp, 0, tmp, n);
                    ap::vadd(&aaa(offs+j, offs+j+1), 1, &tmp(n), 1, "N", ap::vlen(offs+j+1,offs+n-1));
                }
                r = 1/ajj;
                ap::vmul(&aaa(offs+j, offs+j+1), 1, ap::vlen(offs+j+1,offs+n-1), r);
            }
        }
    }
    else
    {
        
        //
        // Compute the Cholesky factorization A = L*L'.
        //
        for(j = 0; j <= n-1; j++)
        {
            
            //
            // Compute L(J+1,J+1) and test for non-positive-definiteness.
            //
            v = ap::vdotproduct(&aaa(offs+j, offs), 1, "Conj", &aaa(offs+j, offs), 1, "N", ap::vlen(offs,offs+j-1));
            ajj = (aaa(offs+j,offs+j)-v).x;
            if( ap::fp_less_eq(ajj,0) )
            {
                aaa(offs+j,offs+j) = ajj;
                result = false;
                return result;
            }
            ajj = sqrt(ajj);
            aaa(offs+j,offs+j) = ajj;
            
            //
            // Compute elements J+1:N of column J.
            //
            if( j<n-1 )
            {
                if( j>0 )
                {
                    ap::vmove(&tmp(0), 1, &aaa(offs+j, offs), 1, "Conj", ap::vlen(0,j-1));
                    cmatrixmv(n-j-1, j, aaa, offs+j+1, offs, 0, tmp, 0, tmp, n);
                    for(i = 0; i <= n-j-2; i++)
                    {
                        aaa(offs+j+1+i,offs+j) = (aaa(offs+j+1+i,offs+j)-tmp(n+i))/ajj;
                    }
                }
                else
                {
                    for(i = 0; i <= n-j-2; i++)
                    {
                        aaa(offs+j+1+i,offs+j) = aaa(offs+j+1+i,offs+j)/ajj;
                    }
                }
            }
        }
    }
    return result;
}


/*************************************************************************
Level-2 Cholesky subroutine

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************/
static bool spdmatrixcholesky2(ap::real_2d_array& aaa,
     int offs,
     int n,
     bool isupper,
     ap::real_1d_array& tmp)
{
    bool result;
    int i;
    int j;
    int k;
    int j1;
    int j2;
    double ajj;
    double v;
    double r;

    result = true;
    if( n<0 )
    {
        result = false;
        return result;
    }
    
    //
    // Quick return if possible
    //
    if( n==0 )
    {
        return result;
    }
    if( isupper )
    {
        
        //
        // Compute the Cholesky factorization A = U'*U.
        //
        for(j = 0; j <= n-1; j++)
        {
            
            //
            // Compute U(J,J) and test for non-positive-definiteness.
            //
            v = ap::vdotproduct(&aaa(offs, offs+j), aaa.getstride(), &aaa(offs, offs+j), aaa.getstride(), ap::vlen(offs,offs+j-1));
            ajj = aaa(offs+j,offs+j)-v;
            if( ap::fp_less_eq(ajj,0) )
            {
                aaa(offs+j,offs+j) = ajj;
                result = false;
                return result;
            }
            ajj = sqrt(ajj);
            aaa(offs+j,offs+j) = ajj;
            
            //
            // Compute elements J+1:N-1 of row J.
            //
            if( j<n-1 )
            {
                if( j>0 )
                {
                    ap::vmoveneg(&tmp(0), 1, &aaa(offs, offs+j), aaa.getstride(), ap::vlen(0,j-1));
                    rmatrixmv(n-j-1, j, aaa, offs, offs+j+1, 1, tmp, 0, tmp, n);
                    ap::vadd(&aaa(offs+j, offs+j+1), 1, &tmp(n), 1, ap::vlen(offs+j+1,offs+n-1));
                }
                r = 1/ajj;
                ap::vmul(&aaa(offs+j, offs+j+1), 1, ap::vlen(offs+j+1,offs+n-1), r);
            }
        }
    }
    else
    {
        
        //
        // Compute the Cholesky factorization A = L*L'.
        //
        for(j = 0; j <= n-1; j++)
        {
            
            //
            // Compute L(J+1,J+1) and test for non-positive-definiteness.
            //
            v = ap::vdotproduct(&aaa(offs+j, offs), 1, &aaa(offs+j, offs), 1, ap::vlen(offs,offs+j-1));
            ajj = aaa(offs+j,offs+j)-v;
            if( ap::fp_less_eq(ajj,0) )
            {
                aaa(offs+j,offs+j) = ajj;
                result = false;
                return result;
            }
            ajj = sqrt(ajj);
            aaa(offs+j,offs+j) = ajj;
            
            //
            // Compute elements J+1:N of column J.
            //
            if( j<n-1 )
            {
                if( j>0 )
                {
                    ap::vmove(&tmp(0), 1, &aaa(offs+j, offs), 1, ap::vlen(0,j-1));
                    rmatrixmv(n-j-1, j, aaa, offs+j+1, offs, 0, tmp, 0, tmp, n);
                    for(i = 0; i <= n-j-2; i++)
                    {
                        aaa(offs+j+1+i,offs+j) = (aaa(offs+j+1+i,offs+j)-tmp(n+i))/ajj;
                    }
                }
                else
                {
                    for(i = 0; i <= n-j-2; i++)
                    {
                        aaa(offs+j+1+i,offs+j) = aaa(offs+j+1+i,offs+j)/ajj;
                    }
                }
            }
        }
    }
    return result;
}




