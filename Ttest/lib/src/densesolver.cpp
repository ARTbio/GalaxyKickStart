/*************************************************************************
Copyright (c) 2007-2008, Sergey Bochkanov (ALGLIB project).

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
#include "densesolver.h"

static void rmatrixlusolveinternal(const ap::real_2d_array& lua,
     const ap::integer_1d_array& p,
     const double& scalea,
     int n,
     const ap::real_2d_array& a,
     bool havea,
     const ap::real_2d_array& b,
     int m,
     int& info,
     densesolverreport& rep,
     ap::real_2d_array& x);
static void spdmatrixcholeskysolveinternal(const ap::real_2d_array& cha,
     const double& sqrtscalea,
     int n,
     bool isupper,
     const ap::real_2d_array& a,
     bool havea,
     const ap::real_2d_array& b,
     int m,
     int& info,
     densesolverreport& rep,
     ap::real_2d_array& x);
static void cmatrixlusolveinternal(const ap::complex_2d_array& lua,
     const ap::integer_1d_array& p,
     const double& scalea,
     int n,
     const ap::complex_2d_array& a,
     bool havea,
     const ap::complex_2d_array& b,
     int m,
     int& info,
     densesolverreport& rep,
     ap::complex_2d_array& x);
static void hpdmatrixcholeskysolveinternal(const ap::complex_2d_array& cha,
     const double& sqrtscalea,
     int n,
     bool isupper,
     const ap::complex_2d_array& a,
     bool havea,
     const ap::complex_2d_array& b,
     int m,
     int& info,
     densesolverreport& rep,
     ap::complex_2d_array& x);
static int densesolverrfsmax(int n, double r1, double rinf);
static int densesolverrfsmaxv2(int n, double r2);
static void rbasiclusolve(const ap::real_2d_array& lua,
     const ap::integer_1d_array& p,
     double scalea,
     int n,
     ap::real_1d_array& xb,
     ap::real_1d_array& tmp);
static void spdbasiccholeskysolve(const ap::real_2d_array& cha,
     double sqrtscalea,
     int n,
     bool isupper,
     ap::real_1d_array& xb,
     ap::real_1d_array& tmp);
static void cbasiclusolve(const ap::complex_2d_array& lua,
     const ap::integer_1d_array& p,
     double scalea,
     int n,
     ap::complex_1d_array& xb,
     ap::complex_1d_array& tmp);
static void hpdbasiccholeskysolve(const ap::complex_2d_array& cha,
     double sqrtscalea,
     int n,
     bool isupper,
     ap::complex_1d_array& xb,
     ap::complex_1d_array& tmp);

/*************************************************************************
Dense solver.

This  subroutine  solves  a  system  A*x=b,  where A is NxN non-denegerate
real matrix, x and b are vectors.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^3) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is singular, or VERY close to singular.
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   solver report, see below for more info
    X       -   array[0..N-1], it contains:
                * solution of A*x=b if A is non-singular (well-conditioned
                  or ill-conditioned, but not very close to singular)
                * zeros,  if  A  is  singular  or  VERY  close to singular
                  (in this case Info=-3).

SOLVER REPORT

Subroutine sets following fields of the Rep structure:
* R1        reciprocal of condition number: 1/cond(A), 1-norm.
* RInf      reciprocal of condition number: 1/cond(A), inf-norm.

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixsolve(const ap::real_2d_array& a,
     int n,
     const ap::real_1d_array& b,
     int& info,
     densesolverreport& rep,
     ap::real_1d_array& x)
{
    ap::real_2d_array bm;
    ap::real_2d_array xm;

    if( n<=0 )
    {
        info = -1;
        return;
    }
    bm.setlength(n, 1);
    ap::vmove(&bm(0, 0), bm.getstride(), &b(0), 1, ap::vlen(0,n-1));
    rmatrixsolvem(a, n, bm, 1, true, info, rep, xm);
    x.setlength(n);
    ap::vmove(&x(0), 1, &xm(0, 0), xm.getstride(), ap::vlen(0,n-1));
}


/*************************************************************************
Dense solver.

Similar to RMatrixSolve() but solves task with multiple right parts (where
b and x are NxM matrices).

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* optional iterative refinement
* O(N^3+M*N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size
    RFS     -   iterative refinement switch:
                * True - refinement is used.
                  Less performance, more precision.
                * False - refinement is not used.
                  More performance, less precision.

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixsolvem(const ap::real_2d_array& a,
     int n,
     const ap::real_2d_array& b,
     int m,
     bool rfs,
     int& info,
     densesolverreport& rep,
     ap::real_2d_array& x)
{
    ap::real_2d_array da;
    ap::real_2d_array emptya;
    ap::integer_1d_array p;
    double scalea;
    int i;
    int j;

    
    //
    // prepare: check inputs, allocate space...
    //
    if( n<=0||m<=0 )
    {
        info = -1;
        return;
    }
    da.setlength(n, n);
    
    //
    // 1. scale matrix, max(|A[i,j]|)
    // 2. factorize scaled matrix
    // 3. solve
    //
    scalea = 0;
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            scalea = ap::maxreal(scalea, fabs(a(i,j)));
        }
    }
    if( ap::fp_eq(scalea,0) )
    {
        scalea = 1;
    }
    scalea = 1/scalea;
    for(i = 0; i <= n-1; i++)
    {
        ap::vmove(&da(i, 0), 1, &a(i, 0), 1, ap::vlen(0,n-1));
    }
    rmatrixlu(da, n, n, p);
    if( rfs )
    {
        rmatrixlusolveinternal(da, p, scalea, n, a, true, b, m, info, rep, x);
    }
    else
    {
        rmatrixlusolveinternal(da, p, scalea, n, emptya, false, b, m, info, rep, x);
    }
}


/*************************************************************************
Dense solver.

This  subroutine  solves  a  system  A*X=B,  where A is NxN non-denegerate
real matrix given by its LU decomposition, X and B are NxM real matrices.

Algorithm features:
* automatic detection of degenerate cases
* O(N^2) complexity
* condition number estimation

No iterative refinement  is provided because exact form of original matrix
is not known to subroutine. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve
    
  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixlusolve(const ap::real_2d_array& lua,
     const ap::integer_1d_array& p,
     int n,
     const ap::real_1d_array& b,
     int& info,
     densesolverreport& rep,
     ap::real_1d_array& x)
{
    ap::real_2d_array bm;
    ap::real_2d_array xm;

    if( n<=0 )
    {
        info = -1;
        return;
    }
    bm.setlength(n, 1);
    ap::vmove(&bm(0, 0), bm.getstride(), &b(0), 1, ap::vlen(0,n-1));
    rmatrixlusolvem(lua, p, n, bm, 1, info, rep, xm);
    x.setlength(n);
    ap::vmove(&x(0), 1, &xm(0, 0), xm.getstride(), ap::vlen(0,n-1));
}


/*************************************************************************
Dense solver.

Similar to RMatrixLUSolve() but solves task with multiple right parts
(where b and x are NxM matrices).

Algorithm features:
* automatic detection of degenerate cases
* O(M*N^2) complexity
* condition number estimation

No iterative refinement  is provided because exact form of original matrix
is not known to subroutine. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixlusolvem(const ap::real_2d_array& lua,
     const ap::integer_1d_array& p,
     int n,
     const ap::real_2d_array& b,
     int m,
     int& info,
     densesolverreport& rep,
     ap::real_2d_array& x)
{
    ap::real_2d_array emptya;
    int i;
    int j;
    double scalea;

    
    //
    // prepare: check inputs, allocate space...
    //
    if( n<=0||m<=0 )
    {
        info = -1;
        return;
    }
    
    //
    // 1. scale matrix, max(|U[i,j]|)
    //    we assume that LU is in its normal form, i.e. |L[i,j]|<=1
    // 2. solve
    //
    scalea = 0;
    for(i = 0; i <= n-1; i++)
    {
        for(j = i; j <= n-1; j++)
        {
            scalea = ap::maxreal(scalea, fabs(lua(i,j)));
        }
    }
    if( ap::fp_eq(scalea,0) )
    {
        scalea = 1;
    }
    scalea = 1/scalea;
    rmatrixlusolveinternal(lua, p, scalea, n, emptya, false, b, m, info, rep, x);
}


/*************************************************************************
Dense solver.

This  subroutine  solves  a  system  A*x=b,  where BOTH ORIGINAL A AND ITS
LU DECOMPOSITION ARE KNOWN. You can use it if for some  reasons  you  have
both A and its LU decomposition.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolveM
    Rep     -   same as in RMatrixSolveM
    X       -   same as in RMatrixSolveM

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixmixedsolve(const ap::real_2d_array& a,
     const ap::real_2d_array& lua,
     const ap::integer_1d_array& p,
     int n,
     const ap::real_1d_array& b,
     int& info,
     densesolverreport& rep,
     ap::real_1d_array& x)
{
    ap::real_2d_array bm;
    ap::real_2d_array xm;

    if( n<=0 )
    {
        info = -1;
        return;
    }
    bm.setlength(n, 1);
    ap::vmove(&bm(0, 0), bm.getstride(), &b(0), 1, ap::vlen(0,n-1));
    rmatrixmixedsolvem(a, lua, p, n, bm, 1, info, rep, xm);
    x.setlength(n);
    ap::vmove(&x(0), 1, &xm(0, 0), xm.getstride(), ap::vlen(0,n-1));
}


/*************************************************************************
Dense solver.

Similar to RMatrixMixedSolve() but  solves task with multiple right  parts
(where b and x are NxM matrices).

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(M*N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolveM
    Rep     -   same as in RMatrixSolveM
    X       -   same as in RMatrixSolveM

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixmixedsolvem(const ap::real_2d_array& a,
     const ap::real_2d_array& lua,
     const ap::integer_1d_array& p,
     int n,
     const ap::real_2d_array& b,
     int m,
     int& info,
     densesolverreport& rep,
     ap::real_2d_array& x)
{
    double scalea;
    int i;
    int j;

    
    //
    // prepare: check inputs, allocate space...
    //
    if( n<=0||m<=0 )
    {
        info = -1;
        return;
    }
    
    //
    // 1. scale matrix, max(|A[i,j]|)
    // 2. factorize scaled matrix
    // 3. solve
    //
    scalea = 0;
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            scalea = ap::maxreal(scalea, fabs(a(i,j)));
        }
    }
    if( ap::fp_eq(scalea,0) )
    {
        scalea = 1;
    }
    scalea = 1/scalea;
    rmatrixlusolveinternal(lua, p, scalea, n, a, true, b, m, info, rep, x);
}


/*************************************************************************
Dense solver. Same as RMatrixSolveM(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^3+M*N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size
    RFS     -   iterative refinement switch:
                * True - refinement is used.
                  Less performance, more precision.
                * False - refinement is not used.
                  More performance, less precision.

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixsolvem(const ap::complex_2d_array& a,
     int n,
     const ap::complex_2d_array& b,
     int m,
     bool rfs,
     int& info,
     densesolverreport& rep,
     ap::complex_2d_array& x)
{
    ap::complex_2d_array da;
    ap::complex_2d_array emptya;
    ap::integer_1d_array p;
    double scalea;
    int i;
    int j;

    
    //
    // prepare: check inputs, allocate space...
    //
    if( n<=0||m<=0 )
    {
        info = -1;
        return;
    }
    da.setlength(n, n);
    
    //
    // 1. scale matrix, max(|A[i,j]|)
    // 2. factorize scaled matrix
    // 3. solve
    //
    scalea = 0;
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            scalea = ap::maxreal(scalea, ap::abscomplex(a(i,j)));
        }
    }
    if( ap::fp_eq(scalea,0) )
    {
        scalea = 1;
    }
    scalea = 1/scalea;
    for(i = 0; i <= n-1; i++)
    {
        ap::vmove(&da(i, 0), 1, &a(i, 0), 1, "N", ap::vlen(0,n-1));
    }
    cmatrixlu(da, n, n, p);
    if( rfs )
    {
        cmatrixlusolveinternal(da, p, scalea, n, a, true, b, m, info, rep, x);
    }
    else
    {
        cmatrixlusolveinternal(da, p, scalea, n, emptya, false, b, m, info, rep, x);
    }
}


/*************************************************************************
Dense solver. Same as RMatrixSolve(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^3) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixsolve(const ap::complex_2d_array& a,
     int n,
     const ap::complex_1d_array& b,
     int& info,
     densesolverreport& rep,
     ap::complex_1d_array& x)
{
    ap::complex_2d_array bm;
    ap::complex_2d_array xm;

    if( n<=0 )
    {
        info = -1;
        return;
    }
    bm.setlength(n, 1);
    ap::vmove(&bm(0, 0), bm.getstride(), &b(0), 1, "N", ap::vlen(0,n-1));
    cmatrixsolvem(a, n, bm, 1, true, info, rep, xm);
    x.setlength(n);
    ap::vmove(&x(0), 1, &xm(0, 0), xm.getstride(), "N", ap::vlen(0,n-1));
}


/*************************************************************************
Dense solver. Same as RMatrixLUSolveM(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* O(M*N^2) complexity
* condition number estimation

No iterative refinement  is provided because exact form of original matrix
is not known to subroutine. Use CMatrixSolve or CMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixlusolvem(const ap::complex_2d_array& lua,
     const ap::integer_1d_array& p,
     int n,
     const ap::complex_2d_array& b,
     int m,
     int& info,
     densesolverreport& rep,
     ap::complex_2d_array& x)
{
    ap::complex_2d_array emptya;
    int i;
    int j;
    double scalea;

    
    //
    // prepare: check inputs, allocate space...
    //
    if( n<=0||m<=0 )
    {
        info = -1;
        return;
    }
    
    //
    // 1. scale matrix, max(|U[i,j]|)
    //    we assume that LU is in its normal form, i.e. |L[i,j]|<=1
    // 2. solve
    //
    scalea = 0;
    for(i = 0; i <= n-1; i++)
    {
        for(j = i; j <= n-1; j++)
        {
            scalea = ap::maxreal(scalea, ap::abscomplex(lua(i,j)));
        }
    }
    if( ap::fp_eq(scalea,0) )
    {
        scalea = 1;
    }
    scalea = 1/scalea;
    cmatrixlusolveinternal(lua, p, scalea, n, emptya, false, b, m, info, rep, x);
}


/*************************************************************************
Dense solver. Same as RMatrixLUSolve(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* O(N^2) complexity
* condition number estimation

No iterative refinement is provided because exact form of original matrix
is not known to subroutine. Use CMatrixSolve or CMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
    P       -   array[0..N-1], pivots array, CMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixlusolve(const ap::complex_2d_array& lua,
     const ap::integer_1d_array& p,
     int n,
     const ap::complex_1d_array& b,
     int& info,
     densesolverreport& rep,
     ap::complex_1d_array& x)
{
    ap::complex_2d_array bm;
    ap::complex_2d_array xm;

    if( n<=0 )
    {
        info = -1;
        return;
    }
    bm.setlength(n, 1);
    ap::vmove(&bm(0, 0), bm.getstride(), &b(0), 1, "N", ap::vlen(0,n-1));
    cmatrixlusolvem(lua, p, n, bm, 1, info, rep, xm);
    x.setlength(n);
    ap::vmove(&x(0), 1, &xm(0, 0), xm.getstride(), "N", ap::vlen(0,n-1));
}


/*************************************************************************
Dense solver. Same as RMatrixMixedSolveM(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(M*N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
    P       -   array[0..N-1], pivots array, CMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolveM
    Rep     -   same as in RMatrixSolveM
    X       -   same as in RMatrixSolveM

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixmixedsolvem(const ap::complex_2d_array& a,
     const ap::complex_2d_array& lua,
     const ap::integer_1d_array& p,
     int n,
     const ap::complex_2d_array& b,
     int m,
     int& info,
     densesolverreport& rep,
     ap::complex_2d_array& x)
{
    double scalea;
    int i;
    int j;

    
    //
    // prepare: check inputs, allocate space...
    //
    if( n<=0||m<=0 )
    {
        info = -1;
        return;
    }
    
    //
    // 1. scale matrix, max(|A[i,j]|)
    // 2. factorize scaled matrix
    // 3. solve
    //
    scalea = 0;
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            scalea = ap::maxreal(scalea, ap::abscomplex(a(i,j)));
        }
    }
    if( ap::fp_eq(scalea,0) )
    {
        scalea = 1;
    }
    scalea = 1/scalea;
    cmatrixlusolveinternal(lua, p, scalea, n, a, true, b, m, info, rep, x);
}


/*************************************************************************
Dense solver. Same as RMatrixMixedSolve(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
    P       -   array[0..N-1], pivots array, CMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolveM
    Rep     -   same as in RMatrixSolveM
    X       -   same as in RMatrixSolveM

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixmixedsolve(const ap::complex_2d_array& a,
     const ap::complex_2d_array& lua,
     const ap::integer_1d_array& p,
     int n,
     const ap::complex_1d_array& b,
     int& info,
     densesolverreport& rep,
     ap::complex_1d_array& x)
{
    ap::complex_2d_array bm;
    ap::complex_2d_array xm;

    if( n<=0 )
    {
        info = -1;
        return;
    }
    bm.setlength(n, 1);
    ap::vmove(&bm(0, 0), bm.getstride(), &b(0), 1, "N", ap::vlen(0,n-1));
    cmatrixmixedsolvem(a, lua, p, n, bm, 1, info, rep, xm);
    x.setlength(n);
    ap::vmove(&x(0), 1, &xm(0, 0), xm.getstride(), "N", ap::vlen(0,n-1));
}


/*************************************************************************
Dense solver. Same as RMatrixSolveM(), but for symmetric positive definite
matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* O(N^3+M*N^2) complexity
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve.
                Returns -3 for non-SPD matrices.
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void spdmatrixsolvem(const ap::real_2d_array& a,
     int n,
     bool isupper,
     const ap::real_2d_array& b,
     int m,
     int& info,
     densesolverreport& rep,
     ap::real_2d_array& x)
{
    ap::real_2d_array da;
    double sqrtscalea;
    int i;
    int j;
    int j1;
    int j2;

    
    //
    // prepare: check inputs, allocate space...
    //
    if( n<=0||m<=0 )
    {
        info = -1;
        return;
    }
    da.setlength(n, n);
    
    //
    // 1. scale matrix, max(|A[i,j]|)
    // 2. factorize scaled matrix
    // 3. solve
    //
    sqrtscalea = 0;
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
            sqrtscalea = ap::maxreal(sqrtscalea, fabs(a(i,j)));
        }
    }
    if( ap::fp_eq(sqrtscalea,0) )
    {
        sqrtscalea = 1;
    }
    sqrtscalea = 1/sqrtscalea;
    sqrtscalea = sqrt(sqrtscalea);
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
        ap::vmove(&da(i, j1), 1, &a(i, j1), 1, ap::vlen(j1,j2));
    }
    if( !spdmatrixcholesky(da, n, isupper) )
    {
        x.setlength(n, m);
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                x(i,j) = 0;
            }
        }
        rep.r1 = 0;
        rep.rinf = 0;
        info = -3;
        return;
    }
    info = 1;
    spdmatrixcholeskysolveinternal(da, sqrtscalea, n, isupper, a, true, b, m, info, rep, x);
}


/*************************************************************************
Dense solver. Same as RMatrixSolve(), but for SPD matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* O(N^3) complexity
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
                Returns -3 for non-SPD matrices.
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void spdmatrixsolve(const ap::real_2d_array& a,
     int n,
     bool isupper,
     const ap::real_1d_array& b,
     int& info,
     densesolverreport& rep,
     ap::real_1d_array& x)
{
    ap::real_2d_array bm;
    ap::real_2d_array xm;

    if( n<=0 )
    {
        info = -1;
        return;
    }
    bm.setlength(n, 1);
    ap::vmove(&bm(0, 0), bm.getstride(), &b(0), 1, ap::vlen(0,n-1));
    spdmatrixsolvem(a, n, isupper, bm, 1, info, rep, xm);
    x.setlength(n);
    ap::vmove(&x(0), 1, &xm(0, 0), xm.getstride(), ap::vlen(0,n-1));
}


/*************************************************************************
Dense solver. Same as RMatrixLUSolveM(), but for SPD matrices  represented
by their Cholesky decomposition.

Algorithm features:
* automatic detection of degenerate cases
* O(M*N^2) complexity
* condition number estimation
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                SPDMatrixCholesky result
    N       -   size of CHA
    IsUpper -   what half of CHA is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void spdmatrixcholeskysolvem(const ap::real_2d_array& cha,
     int n,
     bool isupper,
     const ap::real_2d_array& b,
     int m,
     int& info,
     densesolverreport& rep,
     ap::real_2d_array& x)
{
    ap::real_2d_array emptya;
    double sqrtscalea;
    int i;
    int j;
    int j1;
    int j2;

    
    //
    // prepare: check inputs, allocate space...
    //
    if( n<=0||m<=0 )
    {
        info = -1;
        return;
    }
    
    //
    // 1. scale matrix, max(|U[i,j]|)
    // 2. factorize scaled matrix
    // 3. solve
    //
    sqrtscalea = 0;
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
            sqrtscalea = ap::maxreal(sqrtscalea, fabs(cha(i,j)));
        }
    }
    if( ap::fp_eq(sqrtscalea,0) )
    {
        sqrtscalea = 1;
    }
    sqrtscalea = 1/sqrtscalea;
    spdmatrixcholeskysolveinternal(cha, sqrtscalea, n, isupper, emptya, false, b, m, info, rep, x);
}


/*************************************************************************
Dense solver. Same as RMatrixLUSolve(), but for  SPD matrices  represented
by their Cholesky decomposition.

Algorithm features:
* automatic detection of degenerate cases
* O(N^2) complexity
* condition number estimation
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                SPDMatrixCholesky result
    N       -   size of A
    IsUpper -   what half of CHA is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void spdmatrixcholeskysolve(const ap::real_2d_array& cha,
     int n,
     bool isupper,
     const ap::real_1d_array& b,
     int& info,
     densesolverreport& rep,
     ap::real_1d_array& x)
{
    ap::real_2d_array bm;
    ap::real_2d_array xm;

    if( n<=0 )
    {
        info = -1;
        return;
    }
    bm.setlength(n, 1);
    ap::vmove(&bm(0, 0), bm.getstride(), &b(0), 1, ap::vlen(0,n-1));
    spdmatrixcholeskysolvem(cha, n, isupper, bm, 1, info, rep, xm);
    x.setlength(n);
    ap::vmove(&x(0), 1, &xm(0, 0), xm.getstride(), ap::vlen(0,n-1));
}


/*************************************************************************
Dense solver. Same as RMatrixSolveM(), but for Hermitian positive definite
matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* O(N^3+M*N^2) complexity
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve.
                Returns -3 for non-HPD matrices.
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixsolvem(const ap::complex_2d_array& a,
     int n,
     bool isupper,
     const ap::complex_2d_array& b,
     int m,
     int& info,
     densesolverreport& rep,
     ap::complex_2d_array& x)
{
    ap::complex_2d_array da;
    double sqrtscalea;
    int i;
    int j;
    int j1;
    int j2;

    
    //
    // prepare: check inputs, allocate space...
    //
    if( n<=0||m<=0 )
    {
        info = -1;
        return;
    }
    da.setlength(n, n);
    
    //
    // 1. scale matrix, max(|A[i,j]|)
    // 2. factorize scaled matrix
    // 3. solve
    //
    sqrtscalea = 0;
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
            sqrtscalea = ap::maxreal(sqrtscalea, ap::abscomplex(a(i,j)));
        }
    }
    if( ap::fp_eq(sqrtscalea,0) )
    {
        sqrtscalea = 1;
    }
    sqrtscalea = 1/sqrtscalea;
    sqrtscalea = sqrt(sqrtscalea);
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
        ap::vmove(&da(i, j1), 1, &a(i, j1), 1, "N", ap::vlen(j1,j2));
    }
    if( !hpdmatrixcholesky(da, n, isupper) )
    {
        x.setlength(n, m);
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                x(i,j) = 0;
            }
        }
        rep.r1 = 0;
        rep.rinf = 0;
        info = -3;
        return;
    }
    info = 1;
    hpdmatrixcholeskysolveinternal(da, sqrtscalea, n, isupper, a, true, b, m, info, rep, x);
}


/*************************************************************************
Dense solver. Same as RMatrixSolve(),  but for Hermitian positive definite
matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* O(N^3) complexity
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
                Returns -3 for non-HPD matrices.
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixsolve(const ap::complex_2d_array& a,
     int n,
     bool isupper,
     const ap::complex_1d_array& b,
     int& info,
     densesolverreport& rep,
     ap::complex_1d_array& x)
{
    ap::complex_2d_array bm;
    ap::complex_2d_array xm;

    if( n<=0 )
    {
        info = -1;
        return;
    }
    bm.setlength(n, 1);
    ap::vmove(&bm(0, 0), bm.getstride(), &b(0), 1, "N", ap::vlen(0,n-1));
    hpdmatrixsolvem(a, n, isupper, bm, 1, info, rep, xm);
    x.setlength(n);
    ap::vmove(&x(0), 1, &xm(0, 0), xm.getstride(), "N", ap::vlen(0,n-1));
}


/*************************************************************************
Dense solver. Same as RMatrixLUSolveM(), but for HPD matrices  represented
by their Cholesky decomposition.

Algorithm features:
* automatic detection of degenerate cases
* O(M*N^2) complexity
* condition number estimation
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                HPDMatrixCholesky result
    N       -   size of CHA
    IsUpper -   what half of CHA is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixcholeskysolvem(const ap::complex_2d_array& cha,
     int n,
     bool isupper,
     const ap::complex_2d_array& b,
     int m,
     int& info,
     densesolverreport& rep,
     ap::complex_2d_array& x)
{
    ap::complex_2d_array emptya;
    double sqrtscalea;
    int i;
    int j;
    int j1;
    int j2;

    
    //
    // prepare: check inputs, allocate space...
    //
    if( n<=0||m<=0 )
    {
        info = -1;
        return;
    }
    
    //
    // 1. scale matrix, max(|U[i,j]|)
    // 2. factorize scaled matrix
    // 3. solve
    //
    sqrtscalea = 0;
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
            sqrtscalea = ap::maxreal(sqrtscalea, ap::abscomplex(cha(i,j)));
        }
    }
    if( ap::fp_eq(sqrtscalea,0) )
    {
        sqrtscalea = 1;
    }
    sqrtscalea = 1/sqrtscalea;
    hpdmatrixcholeskysolveinternal(cha, sqrtscalea, n, isupper, emptya, false, b, m, info, rep, x);
}


/*************************************************************************
Dense solver. Same as RMatrixLUSolve(), but for  HPD matrices  represented
by their Cholesky decomposition.

Algorithm features:
* automatic detection of degenerate cases
* O(N^2) complexity
* condition number estimation
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                SPDMatrixCholesky result
    N       -   size of A
    IsUpper -   what half of CHA is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void hpdmatrixcholeskysolve(const ap::complex_2d_array& cha,
     int n,
     bool isupper,
     const ap::complex_1d_array& b,
     int& info,
     densesolverreport& rep,
     ap::complex_1d_array& x)
{
    ap::complex_2d_array bm;
    ap::complex_2d_array xm;

    if( n<=0 )
    {
        info = -1;
        return;
    }
    bm.setlength(n, 1);
    ap::vmove(&bm(0, 0), bm.getstride(), &b(0), 1, "N", ap::vlen(0,n-1));
    hpdmatrixcholeskysolvem(cha, n, isupper, bm, 1, info, rep, xm);
    x.setlength(n);
    ap::vmove(&x(0), 1, &xm(0, 0), xm.getstride(), "N", ap::vlen(0,n-1));
}


/*************************************************************************
Dense solver.

This subroutine finds solution of the linear system A*X=B with non-square,
possibly degenerate A.  System  is  solved in the least squares sense, and
general least squares solution  X = X0 + CX*y  which  minimizes |A*X-B| is
returned. If A is non-degenerate, solution in the  usual sense is returned

Algorithm features:
* automatic detection of degenerate cases
* iterative refinement
* O(N^3) complexity

INPUT PARAMETERS
    A       -   array[0..NRows-1,0..NCols-1], system matrix
    NRows   -   vertical size of A
    NCols   -   horizontal size of A
    B       -   array[0..NCols-1], right part
    Threshold-  a number in [0,1]. Singular values  beyond  Threshold  are
                considered  zero.  Set  it to 0.0, if you don't understand
                what it means, so the solver will choose good value on its
                own.
                
OUTPUT PARAMETERS
    Info    -   return code:
                * -4    SVD subroutine failed
                * -1    if NRows<=0 or NCols<=0 or Threshold<0 was passed
                *  1    if task is solved
    Rep     -   solver report, see below for more info
    X       -   array[0..N-1,0..M-1], it contains:
                * solution of A*X=B if A is non-singular (well-conditioned
                  or ill-conditioned, but not very close to singular)
                * zeros,  if  A  is  singular  or  VERY  close to singular
                  (in this case Info=-3).

SOLVER REPORT

Subroutine sets following fields of the Rep structure:
* R2        reciprocal of condition number: 1/cond(A), 2-norm.
* N         = NCols
* K         dim(Null(A))
* CX        array[0..N-1,0..K-1], kernel of A.
            Columns of CX store such vectors that A*CX[i]=0.

  -- ALGLIB --
     Copyright 24.08.2009 by Bochkanov Sergey
*************************************************************************/
void rmatrixsolvels(const ap::real_2d_array& a,
     int nrows,
     int ncols,
     const ap::real_1d_array& b,
     double threshold,
     int& info,
     densesolverlsreport& rep,
     ap::real_1d_array& x)
{
    ap::real_1d_array sv;
    ap::real_2d_array u;
    ap::real_2d_array vt;
    ap::real_1d_array rp;
    ap::real_1d_array utb;
    ap::real_1d_array sutb;
    ap::real_1d_array tmp;
    ap::real_1d_array ta;
    ap::real_1d_array tx;
    ap::real_1d_array buf;
    ap::real_1d_array w;
    int i;
    int j;
    int nsv;
    int kernelidx;
    double v;
    double verr;
    bool svdfailed;
    bool zeroa;
    int rfs;
    int nrfs;
    bool terminatenexttime;
    bool smallerr;

    if( nrows<=0||ncols<=0||ap::fp_less(threshold,0) )
    {
        info = -1;
        return;
    }
    if( ap::fp_eq(threshold,0) )
    {
        threshold = 1000*ap::machineepsilon;
    }
    
    //
    // Factorize A first
    //
    svdfailed = !rmatrixsvd(a, nrows, ncols, 1, 2, 2, sv, u, vt);
    zeroa = ap::fp_eq(sv(0),0);
    if( svdfailed||zeroa )
    {
        if( svdfailed )
        {
            info = -4;
        }
        else
        {
            info = 1;
        }
        x.setlength(ncols);
        for(i = 0; i <= ncols-1; i++)
        {
            x(i) = 0;
        }
        rep.n = ncols;
        rep.k = ncols;
        rep.cx.setlength(ncols, ncols);
        for(i = 0; i <= ncols-1; i++)
        {
            for(j = 0; j <= ncols-1; j++)
            {
                if( i==j )
                {
                    rep.cx(i,j) = 1;
                }
                else
                {
                    rep.cx(i,j) = 0;
                }
            }
        }
        rep.r2 = 0;
        return;
    }
    nsv = ap::minint(ncols, nrows);
    if( nsv==ncols )
    {
        rep.r2 = sv(nsv-1)/sv(0);
    }
    else
    {
        rep.r2 = 0;
    }
    rep.n = ncols;
    info = 1;
    
    //
    // Iterative refinement of xc combined with solution:
    // 1. xc = 0
    // 2. calculate r = bc-A*xc using extra-precise dot product
    // 3. solve A*y = r
    // 4. update x:=x+r
    // 5. goto 2
    //
    // This cycle is executed until one of two things happens:
    // 1. maximum number of iterations reached
    // 2. last iteration decreased error to the lower limit
    //
    utb.setlength(nsv);
    sutb.setlength(nsv);
    x.setlength(ncols);
    tmp.setlength(ncols);
    ta.setlength(ncols+1);
    tx.setlength(ncols+1);
    buf.setlength(ncols+1);
    for(i = 0; i <= ncols-1; i++)
    {
        x(i) = 0;
    }
    kernelidx = nsv;
    for(i = 0; i <= nsv-1; i++)
    {
        if( ap::fp_less_eq(sv(i),threshold*sv(0)) )
        {
            kernelidx = i;
            break;
        }
    }
    rep.k = ncols-kernelidx;
    nrfs = densesolverrfsmaxv2(ncols, rep.r2);
    terminatenexttime = false;
    rp.setlength(nrows);
    for(rfs = 0; rfs <= nrfs; rfs++)
    {
        if( terminatenexttime )
        {
            break;
        }
        
        //
        // calculate right part
        //
        if( rfs==0 )
        {
            ap::vmove(&rp(0), 1, &b(0), 1, ap::vlen(0,nrows-1));
        }
        else
        {
            smallerr = true;
            for(i = 0; i <= nrows-1; i++)
            {
                ap::vmove(&ta(0), 1, &a(i, 0), 1, ap::vlen(0,ncols-1));
                ta(ncols) = -1;
                ap::vmove(&tx(0), 1, &x(0), 1, ap::vlen(0,ncols-1));
                tx(ncols) = b(i);
                xdot(ta, tx, ncols+1, buf, v, verr);
                rp(i) = -v;
                smallerr = smallerr&&ap::fp_less(fabs(v),4*verr);
            }
            if( smallerr )
            {
                terminatenexttime = true;
            }
        }
        
        //
        // solve A*dx = rp
        //
        for(i = 0; i <= ncols-1; i++)
        {
            tmp(i) = 0;
        }
        for(i = 0; i <= nsv-1; i++)
        {
            utb(i) = 0;
        }
        for(i = 0; i <= nrows-1; i++)
        {
            v = rp(i);
            ap::vadd(&utb(0), 1, &u(i, 0), 1, ap::vlen(0,nsv-1), v);
        }
        for(i = 0; i <= nsv-1; i++)
        {
            if( i<kernelidx )
            {
                sutb(i) = utb(i)/sv(i);
            }
            else
            {
                sutb(i) = 0;
            }
        }
        for(i = 0; i <= nsv-1; i++)
        {
            v = sutb(i);
            ap::vadd(&tmp(0), 1, &vt(i, 0), 1, ap::vlen(0,ncols-1), v);
        }
        
        //
        // update x:  x:=x+dx
        //
        ap::vadd(&x(0), 1, &tmp(0), 1, ap::vlen(0,ncols-1));
    }
    
    //
    // fill CX
    //
    if( rep.k>0 )
    {
        rep.cx.setlength(ncols, rep.k);
        for(i = 0; i <= rep.k-1; i++)
        {
            ap::vmove(&rep.cx(0, i), rep.cx.getstride(), &vt(kernelidx+i, 0), 1, ap::vlen(0,ncols-1));
        }
    }
}


/*************************************************************************
Internal LU solver

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
static void rmatrixlusolveinternal(const ap::real_2d_array& lua,
     const ap::integer_1d_array& p,
     const double& scalea,
     int n,
     const ap::real_2d_array& a,
     bool havea,
     const ap::real_2d_array& b,
     int m,
     int& info,
     densesolverreport& rep,
     ap::real_2d_array& x)
{
    int i;
    int j;
    int k;
    int rfs;
    int nrfs;
    ap::real_1d_array xc;
    ap::real_1d_array y;
    ap::real_1d_array bc;
    ap::real_1d_array xa;
    ap::real_1d_array xb;
    ap::real_1d_array tx;
    double v;
    double verr;
    double mxb;
    double scaleright;
    bool smallerr;
    bool terminatenexttime;

    ap::ap_error::make_assertion(ap::fp_greater(scalea,0), "");
    
    //
    // prepare: check inputs, allocate space...
    //
    if( n<=0||m<=0 )
    {
        info = -1;
        return;
    }
    for(i = 0; i <= n-1; i++)
    {
        if( p(i)>n-1||p(i)<i )
        {
            info = -1;
            return;
        }
    }
    x.setlength(n, m);
    y.setlength(n);
    xc.setlength(n);
    bc.setlength(n);
    tx.setlength(n+1);
    xa.setlength(n+1);
    xb.setlength(n+1);
    
    //
    // estimate condition number, test for near singularity
    //
    rep.r1 = rmatrixlurcond1(lua, n);
    rep.rinf = rmatrixlurcondinf(lua, n);
    if( ap::fp_less(rep.r1,rcondthreshold())||ap::fp_less(rep.rinf,rcondthreshold()) )
    {
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                x(i,j) = 0;
            }
        }
        rep.r1 = 0;
        rep.rinf = 0;
        info = -3;
        return;
    }
    info = 1;
    
    //
    // solve
    //
    for(k = 0; k <= m-1; k++)
    {
        
        //
        // copy B to contiguous storage
        //
        ap::vmove(&bc(0), 1, &b(0, k), b.getstride(), ap::vlen(0,n-1));
        
        //
        // Scale right part:
        // * MX stores max(|Bi|)
        // * ScaleRight stores actual scaling applied to B when solving systems
        //   it is chosen to make |scaleRight*b| close to 1.
        //
        mxb = 0;
        for(i = 0; i <= n-1; i++)
        {
            mxb = ap::maxreal(mxb, fabs(bc(i)));
        }
        if( ap::fp_eq(mxb,0) )
        {
            mxb = 1;
        }
        scaleright = 1/mxb;
        
        //
        // First, non-iterative part of solution process.
        // We use separate code for this task because
        // XDot is quite slow and we want to save time.
        //
        ap::vmove(&xc(0), 1, &bc(0), 1, ap::vlen(0,n-1), scaleright);
        rbasiclusolve(lua, p, scalea, n, xc, tx);
        
        //
        // Iterative refinement of xc:
        // * calculate r = bc-A*xc using extra-precise dot product
        // * solve A*y = r
        // * update x:=x+r
        //
        // This cycle is executed until one of two things happens:
        // 1. maximum number of iterations reached
        // 2. last iteration decreased error to the lower limit
        //
        if( havea )
        {
            nrfs = densesolverrfsmax(n, rep.r1, rep.rinf);
            terminatenexttime = false;
            for(rfs = 0; rfs <= nrfs-1; rfs++)
            {
                if( terminatenexttime )
                {
                    break;
                }
                
                //
                // generate right part
                //
                smallerr = true;
                ap::vmove(&xb(0), 1, &xc(0), 1, ap::vlen(0,n-1));
                for(i = 0; i <= n-1; i++)
                {
                    ap::vmove(&xa(0), 1, &a(i, 0), 1, ap::vlen(0,n-1), scalea);
                    xa(n) = -1;
                    xb(n) = scaleright*bc(i);
                    xdot(xa, xb, n+1, tx, v, verr);
                    y(i) = -v;
                    smallerr = smallerr&&ap::fp_less(fabs(v),4*verr);
                }
                if( smallerr )
                {
                    terminatenexttime = true;
                }
                
                //
                // solve and update
                //
                rbasiclusolve(lua, p, scalea, n, y, tx);
                ap::vadd(&xc(0), 1, &y(0), 1, ap::vlen(0,n-1));
            }
        }
        
        //
        // Store xc.
        // Post-scale result.
        //
        v = scalea*mxb;
        ap::vmove(&x(0, k), x.getstride(), &xc(0), 1, ap::vlen(0,n-1), v);
    }
}


/*************************************************************************
Internal Cholesky solver

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
static void spdmatrixcholeskysolveinternal(const ap::real_2d_array& cha,
     const double& sqrtscalea,
     int n,
     bool isupper,
     const ap::real_2d_array& a,
     bool havea,
     const ap::real_2d_array& b,
     int m,
     int& info,
     densesolverreport& rep,
     ap::real_2d_array& x)
{
    int i;
    int j;
    int k;
    int rfs;
    int nrfs;
    ap::real_1d_array xc;
    ap::real_1d_array y;
    ap::real_1d_array bc;
    ap::real_1d_array xa;
    ap::real_1d_array xb;
    ap::real_1d_array tx;
    double v;
    double verr;
    double mxb;
    double scaleright;
    bool smallerr;
    bool terminatenexttime;

    ap::ap_error::make_assertion(ap::fp_greater(sqrtscalea,0), "");
    
    //
    // prepare: check inputs, allocate space...
    //
    if( n<=0||m<=0 )
    {
        info = -1;
        return;
    }
    x.setlength(n, m);
    y.setlength(n);
    xc.setlength(n);
    bc.setlength(n);
    tx.setlength(n+1);
    xa.setlength(n+1);
    xb.setlength(n+1);
    
    //
    // estimate condition number, test for near singularity
    //
    rep.r1 = spdmatrixcholeskyrcond(cha, n, isupper);
    rep.rinf = rep.r1;
    if( ap::fp_less(rep.r1,rcondthreshold()) )
    {
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                x(i,j) = 0;
            }
        }
        rep.r1 = 0;
        rep.rinf = 0;
        info = -3;
        return;
    }
    info = 1;
    
    //
    // solve
    //
    for(k = 0; k <= m-1; k++)
    {
        
        //
        // copy B to contiguous storage
        //
        ap::vmove(&bc(0), 1, &b(0, k), b.getstride(), ap::vlen(0,n-1));
        
        //
        // Scale right part:
        // * MX stores max(|Bi|)
        // * ScaleRight stores actual scaling applied to B when solving systems
        //   it is chosen to make |scaleRight*b| close to 1.
        //
        mxb = 0;
        for(i = 0; i <= n-1; i++)
        {
            mxb = ap::maxreal(mxb, fabs(bc(i)));
        }
        if( ap::fp_eq(mxb,0) )
        {
            mxb = 1;
        }
        scaleright = 1/mxb;
        
        //
        // First, non-iterative part of solution process.
        // We use separate code for this task because
        // XDot is quite slow and we want to save time.
        //
        ap::vmove(&xc(0), 1, &bc(0), 1, ap::vlen(0,n-1), scaleright);
        spdbasiccholeskysolve(cha, sqrtscalea, n, isupper, xc, tx);
        
        //
        // Store xc.
        // Post-scale result.
        //
        v = ap::sqr(sqrtscalea)*mxb;
        ap::vmove(&x(0, k), x.getstride(), &xc(0), 1, ap::vlen(0,n-1), v);
    }
}


/*************************************************************************
Internal LU solver

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
static void cmatrixlusolveinternal(const ap::complex_2d_array& lua,
     const ap::integer_1d_array& p,
     const double& scalea,
     int n,
     const ap::complex_2d_array& a,
     bool havea,
     const ap::complex_2d_array& b,
     int m,
     int& info,
     densesolverreport& rep,
     ap::complex_2d_array& x)
{
    int i;
    int j;
    int k;
    int rfs;
    int nrfs;
    ap::complex_1d_array xc;
    ap::complex_1d_array y;
    ap::complex_1d_array bc;
    ap::complex_1d_array xa;
    ap::complex_1d_array xb;
    ap::complex_1d_array tx;
    ap::real_1d_array tmpbuf;
    ap::complex v;
    double verr;
    double mxb;
    double scaleright;
    bool smallerr;
    bool terminatenexttime;

    ap::ap_error::make_assertion(ap::fp_greater(scalea,0), "");
    
    //
    // prepare: check inputs, allocate space...
    //
    if( n<=0||m<=0 )
    {
        info = -1;
        return;
    }
    for(i = 0; i <= n-1; i++)
    {
        if( p(i)>n-1||p(i)<i )
        {
            info = -1;
            return;
        }
    }
    x.setlength(n, m);
    y.setlength(n);
    xc.setlength(n);
    bc.setlength(n);
    tx.setlength(n);
    xa.setlength(n+1);
    xb.setlength(n+1);
    tmpbuf.setlength(2*n+2);
    
    //
    // estimate condition number, test for near singularity
    //
    rep.r1 = cmatrixlurcond1(lua, n);
    rep.rinf = cmatrixlurcondinf(lua, n);
    if( ap::fp_less(rep.r1,rcondthreshold())||ap::fp_less(rep.rinf,rcondthreshold()) )
    {
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                x(i,j) = 0;
            }
        }
        rep.r1 = 0;
        rep.rinf = 0;
        info = -3;
        return;
    }
    info = 1;
    
    //
    // solve
    //
    for(k = 0; k <= m-1; k++)
    {
        
        //
        // copy B to contiguous storage
        //
        ap::vmove(&bc(0), 1, &b(0, k), b.getstride(), "N", ap::vlen(0,n-1));
        
        //
        // Scale right part:
        // * MX stores max(|Bi|)
        // * ScaleRight stores actual scaling applied to B when solving systems
        //   it is chosen to make |scaleRight*b| close to 1.
        //
        mxb = 0;
        for(i = 0; i <= n-1; i++)
        {
            mxb = ap::maxreal(mxb, ap::abscomplex(bc(i)));
        }
        if( ap::fp_eq(mxb,0) )
        {
            mxb = 1;
        }
        scaleright = 1/mxb;
        
        //
        // First, non-iterative part of solution process.
        // We use separate code for this task because
        // XDot is quite slow and we want to save time.
        //
        ap::vmove(&xc(0), 1, &bc(0), 1, "N", ap::vlen(0,n-1), scaleright);
        cbasiclusolve(lua, p, scalea, n, xc, tx);
        
        //
        // Iterative refinement of xc:
        // * calculate r = bc-A*xc using extra-precise dot product
        // * solve A*y = r
        // * update x:=x+r
        //
        // This cycle is executed until one of two things happens:
        // 1. maximum number of iterations reached
        // 2. last iteration decreased error to the lower limit
        //
        if( havea )
        {
            nrfs = densesolverrfsmax(n, rep.r1, rep.rinf);
            terminatenexttime = false;
            for(rfs = 0; rfs <= nrfs-1; rfs++)
            {
                if( terminatenexttime )
                {
                    break;
                }
                
                //
                // generate right part
                //
                smallerr = true;
                ap::vmove(&xb(0), 1, &xc(0), 1, "N", ap::vlen(0,n-1));
                for(i = 0; i <= n-1; i++)
                {
                    ap::vmove(&xa(0), 1, &a(i, 0), 1, "N", ap::vlen(0,n-1), scalea);
                    xa(n) = -1;
                    xb(n) = scaleright*bc(i);
                    xcdot(xa, xb, n+1, tmpbuf, v, verr);
                    y(i) = -v;
                    smallerr = smallerr&&ap::fp_less(ap::abscomplex(v),4*verr);
                }
                if( smallerr )
                {
                    terminatenexttime = true;
                }
                
                //
                // solve and update
                //
                cbasiclusolve(lua, p, scalea, n, y, tx);
                ap::vadd(&xc(0), 1, &y(0), 1, "N", ap::vlen(0,n-1));
            }
        }
        
        //
        // Store xc.
        // Post-scale result.
        //
        v = scalea*mxb;
        ap::vmove(&x(0, k), x.getstride(), &xc(0), 1, "N", ap::vlen(0,n-1), v);
    }
}


/*************************************************************************
Internal Cholesky solver

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
static void hpdmatrixcholeskysolveinternal(const ap::complex_2d_array& cha,
     const double& sqrtscalea,
     int n,
     bool isupper,
     const ap::complex_2d_array& a,
     bool havea,
     const ap::complex_2d_array& b,
     int m,
     int& info,
     densesolverreport& rep,
     ap::complex_2d_array& x)
{
    int i;
    int j;
    int k;
    int rfs;
    int nrfs;
    ap::complex_1d_array xc;
    ap::complex_1d_array y;
    ap::complex_1d_array bc;
    ap::complex_1d_array xa;
    ap::complex_1d_array xb;
    ap::complex_1d_array tx;
    double v;
    double verr;
    double mxb;
    double scaleright;
    bool smallerr;
    bool terminatenexttime;

    ap::ap_error::make_assertion(ap::fp_greater(sqrtscalea,0), "");
    
    //
    // prepare: check inputs, allocate space...
    //
    if( n<=0||m<=0 )
    {
        info = -1;
        return;
    }
    x.setlength(n, m);
    y.setlength(n);
    xc.setlength(n);
    bc.setlength(n);
    tx.setlength(n+1);
    xa.setlength(n+1);
    xb.setlength(n+1);
    
    //
    // estimate condition number, test for near singularity
    //
    rep.r1 = hpdmatrixcholeskyrcond(cha, n, isupper);
    rep.rinf = rep.r1;
    if( ap::fp_less(rep.r1,rcondthreshold()) )
    {
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                x(i,j) = 0;
            }
        }
        rep.r1 = 0;
        rep.rinf = 0;
        info = -3;
        return;
    }
    info = 1;
    
    //
    // solve
    //
    for(k = 0; k <= m-1; k++)
    {
        
        //
        // copy B to contiguous storage
        //
        ap::vmove(&bc(0), 1, &b(0, k), b.getstride(), "N", ap::vlen(0,n-1));
        
        //
        // Scale right part:
        // * MX stores max(|Bi|)
        // * ScaleRight stores actual scaling applied to B when solving systems
        //   it is chosen to make |scaleRight*b| close to 1.
        //
        mxb = 0;
        for(i = 0; i <= n-1; i++)
        {
            mxb = ap::maxreal(mxb, ap::abscomplex(bc(i)));
        }
        if( ap::fp_eq(mxb,0) )
        {
            mxb = 1;
        }
        scaleright = 1/mxb;
        
        //
        // First, non-iterative part of solution process.
        // We use separate code for this task because
        // XDot is quite slow and we want to save time.
        //
        ap::vmove(&xc(0), 1, &bc(0), 1, "N", ap::vlen(0,n-1), scaleright);
        hpdbasiccholeskysolve(cha, sqrtscalea, n, isupper, xc, tx);
        
        //
        // Store xc.
        // Post-scale result.
        //
        v = ap::sqr(sqrtscalea)*mxb;
        ap::vmove(&x(0, k), x.getstride(), &xc(0), 1, "N", ap::vlen(0,n-1), v);
    }
}


/*************************************************************************
Internal subroutine.
Returns maximum count of RFS iterations as function of:
1. machine epsilon
2. task size.
3. condition number

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
static int densesolverrfsmax(int n, double r1, double rinf)
{
    int result;

    result = 5;
    return result;
}


/*************************************************************************
Internal subroutine.
Returns maximum count of RFS iterations as function of:
1. machine epsilon
2. task size.
3. norm-2 condition number

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
static int densesolverrfsmaxv2(int n, double r2)
{
    int result;

    result = densesolverrfsmax(n, double(0), double(0));
    return result;
}


/*************************************************************************
Basic LU solver for ScaleA*PLU*x = y.

This subroutine assumes that:
* L is well-scaled, and it is U which needs scaling by ScaleA.
* A=PLU is well-conditioned, so no zero divisions or overflow may occur

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
static void rbasiclusolve(const ap::real_2d_array& lua,
     const ap::integer_1d_array& p,
     double scalea,
     int n,
     ap::real_1d_array& xb,
     ap::real_1d_array& tmp)
{
    int i;
    double v;

    for(i = 0; i <= n-1; i++)
    {
        if( p(i)!=i )
        {
            v = xb(i);
            xb(i) = xb(p(i));
            xb(p(i)) = v;
        }
    }
    for(i = 1; i <= n-1; i++)
    {
        v = ap::vdotproduct(&lua(i, 0), 1, &xb(0), 1, ap::vlen(0,i-1));
        xb(i) = xb(i)-v;
    }
    xb(n-1) = xb(n-1)/(scalea*lua(n-1,n-1));
    for(i = n-2; i >= 0; i--)
    {
        ap::vmove(&tmp(i+1), 1, &lua(i, i+1), 1, ap::vlen(i+1,n-1), scalea);
        v = ap::vdotproduct(&tmp(i+1), 1, &xb(i+1), 1, ap::vlen(i+1,n-1));
        xb(i) = (xb(i)-v)/(scalea*lua(i,i));
    }
}


/*************************************************************************
Basic Cholesky solver for ScaleA*Cholesky(A)'*x = y.

This subroutine assumes that:
* A*ScaleA is well scaled
* A is well-conditioned, so no zero divisions or overflow may occur

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
static void spdbasiccholeskysolve(const ap::real_2d_array& cha,
     double sqrtscalea,
     int n,
     bool isupper,
     ap::real_1d_array& xb,
     ap::real_1d_array& tmp)
{
    int i;
    double v;

    
    //
    // A = L*L' or A=U'*U
    //
    if( isupper )
    {
        
        //
        // Solve U'*y=b first.
        //
        for(i = 0; i <= n-1; i++)
        {
            xb(i) = xb(i)/(sqrtscalea*cha(i,i));
            if( i<n-1 )
            {
                v = xb(i);
                ap::vmove(&tmp(i+1), 1, &cha(i, i+1), 1, ap::vlen(i+1,n-1), sqrtscalea);
                ap::vsub(&xb(i+1), 1, &tmp(i+1), 1, ap::vlen(i+1,n-1), v);
            }
        }
        
        //
        // Solve U*x=y then.
        //
        for(i = n-1; i >= 0; i--)
        {
            if( i<n-1 )
            {
                ap::vmove(&tmp(i+1), 1, &cha(i, i+1), 1, ap::vlen(i+1,n-1), sqrtscalea);
                v = ap::vdotproduct(&tmp(i+1), 1, &xb(i+1), 1, ap::vlen(i+1,n-1));
                xb(i) = xb(i)-v;
            }
            xb(i) = xb(i)/(sqrtscalea*cha(i,i));
        }
    }
    else
    {
        
        //
        // Solve L*y=b first
        //
        for(i = 0; i <= n-1; i++)
        {
            if( i>0 )
            {
                ap::vmove(&tmp(0), 1, &cha(i, 0), 1, ap::vlen(0,i-1), sqrtscalea);
                v = ap::vdotproduct(&tmp(0), 1, &xb(0), 1, ap::vlen(0,i-1));
                xb(i) = xb(i)-v;
            }
            xb(i) = xb(i)/(sqrtscalea*cha(i,i));
        }
        
        //
        // Solve L'*x=y then.
        //
        for(i = n-1; i >= 0; i--)
        {
            xb(i) = xb(i)/(sqrtscalea*cha(i,i));
            if( i>0 )
            {
                v = xb(i);
                ap::vmove(&tmp(0), 1, &cha(i, 0), 1, ap::vlen(0,i-1), sqrtscalea);
                ap::vsub(&xb(0), 1, &tmp(0), 1, ap::vlen(0,i-1), v);
            }
        }
    }
}


/*************************************************************************
Basic LU solver for ScaleA*PLU*x = y.

This subroutine assumes that:
* L is well-scaled, and it is U which needs scaling by ScaleA.
* A=PLU is well-conditioned, so no zero divisions or overflow may occur

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
static void cbasiclusolve(const ap::complex_2d_array& lua,
     const ap::integer_1d_array& p,
     double scalea,
     int n,
     ap::complex_1d_array& xb,
     ap::complex_1d_array& tmp)
{
    int i;
    ap::complex v;

    for(i = 0; i <= n-1; i++)
    {
        if( p(i)!=i )
        {
            v = xb(i);
            xb(i) = xb(p(i));
            xb(p(i)) = v;
        }
    }
    for(i = 1; i <= n-1; i++)
    {
        v = ap::vdotproduct(&lua(i, 0), 1, "N", &xb(0), 1, "N", ap::vlen(0,i-1));
        xb(i) = xb(i)-v;
    }
    xb(n-1) = xb(n-1)/(scalea*lua(n-1,n-1));
    for(i = n-2; i >= 0; i--)
    {
        ap::vmove(&tmp(i+1), 1, &lua(i, i+1), 1, "N", ap::vlen(i+1,n-1), scalea);
        v = ap::vdotproduct(&tmp(i+1), 1, "N", &xb(i+1), 1, "N", ap::vlen(i+1,n-1));
        xb(i) = (xb(i)-v)/(scalea*lua(i,i));
    }
}


/*************************************************************************
Basic Cholesky solver for ScaleA*Cholesky(A)'*x = y.

This subroutine assumes that:
* A*ScaleA is well scaled
* A is well-conditioned, so no zero divisions or overflow may occur

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
static void hpdbasiccholeskysolve(const ap::complex_2d_array& cha,
     double sqrtscalea,
     int n,
     bool isupper,
     ap::complex_1d_array& xb,
     ap::complex_1d_array& tmp)
{
    int i;
    ap::complex v;

    
    //
    // A = L*L' or A=U'*U
    //
    if( isupper )
    {
        
        //
        // Solve U'*y=b first.
        //
        for(i = 0; i <= n-1; i++)
        {
            xb(i) = xb(i)/(sqrtscalea*ap::conj(cha(i,i)));
            if( i<n-1 )
            {
                v = xb(i);
                ap::vmove(&tmp(i+1), 1, &cha(i, i+1), 1, "Conj", ap::vlen(i+1,n-1), sqrtscalea);
                ap::vsub(&xb(i+1), 1, &tmp(i+1), 1, "N", ap::vlen(i+1,n-1), v);
            }
        }
        
        //
        // Solve U*x=y then.
        //
        for(i = n-1; i >= 0; i--)
        {
            if( i<n-1 )
            {
                ap::vmove(&tmp(i+1), 1, &cha(i, i+1), 1, "N", ap::vlen(i+1,n-1), sqrtscalea);
                v = ap::vdotproduct(&tmp(i+1), 1, "N", &xb(i+1), 1, "N", ap::vlen(i+1,n-1));
                xb(i) = xb(i)-v;
            }
            xb(i) = xb(i)/(sqrtscalea*cha(i,i));
        }
    }
    else
    {
        
        //
        // Solve L*y=b first
        //
        for(i = 0; i <= n-1; i++)
        {
            if( i>0 )
            {
                ap::vmove(&tmp(0), 1, &cha(i, 0), 1, "N", ap::vlen(0,i-1), sqrtscalea);
                v = ap::vdotproduct(&tmp(0), 1, "N", &xb(0), 1, "N", ap::vlen(0,i-1));
                xb(i) = xb(i)-v;
            }
            xb(i) = xb(i)/(sqrtscalea*cha(i,i));
        }
        
        //
        // Solve L'*x=y then.
        //
        for(i = n-1; i >= 0; i--)
        {
            xb(i) = xb(i)/(sqrtscalea*ap::conj(cha(i,i)));
            if( i>0 )
            {
                v = xb(i);
                ap::vmove(&tmp(0), 1, &cha(i, 0), 1, "Conj", ap::vlen(0,i-1), sqrtscalea);
                ap::vsub(&xb(0), 1, &tmp(0), 1, "N", ap::vlen(0,i-1), v);
            }
        }
    }
}




