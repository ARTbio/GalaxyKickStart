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
#include "matinv.h"

static void rmatrixtrinverserec(ap::real_2d_array& a,
     int offs,
     int n,
     bool isupper,
     bool isunit,
     ap::real_1d_array& tmp,
     int& info,
     matinvreport& rep);
static void cmatrixtrinverserec(ap::complex_2d_array& a,
     int offs,
     int n,
     bool isupper,
     bool isunit,
     ap::complex_1d_array& tmp,
     int& info,
     matinvreport& rep);
static void rmatrixluinverserec(ap::real_2d_array& a,
     int offs,
     int n,
     ap::real_1d_array& work,
     int& info,
     matinvreport& rep);
static void cmatrixluinverserec(ap::complex_2d_array& a,
     int offs,
     int n,
     ap::complex_1d_array& work,
     int& info,
     matinvreport& rep);
static void spdmatrixcholeskyinverserec(ap::real_2d_array& a,
     int offs,
     int n,
     bool isupper,
     ap::real_1d_array& tmp);
static void hpdmatrixcholeskyinverserec(ap::complex_2d_array& a,
     int offs,
     int n,
     bool isupper,
     ap::complex_1d_array& tmp);

/*************************************************************************
Inversion of a matrix given by its LU decomposition.

INPUT PARAMETERS:
    A       -   LU decomposition of the matrix (output of RMatrixLU subroutine).
    Pivots  -   table of permutations which were made during the LU decomposition
                (the output of RMatrixLU subroutine).
    N       -   size of matrix A.

OUTPUT PARAMETERS:
    Info    -   return code:
                * -3    A is singular, or VERY close to singular.
                        it is filled by zeros in such cases.
                * -1    N<=0 was passed, or incorrect Pivots was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   solver report, see below for more info
    A       -   inverse of matrix A.
                Array whose indexes range within [0..N-1, 0..N-1].

SOLVER REPORT

Subroutine sets following fields of the Rep structure:
* R1        reciprocal of condition number: 1/cond(A), 1-norm.
* RInf      reciprocal of condition number: 1/cond(A), inf-norm.

  -- ALGLIB routine --
     05.02.2010
     Bochkanov Sergey
*************************************************************************/
void rmatrixluinverse(ap::real_2d_array& a,
     const ap::integer_1d_array& pivots,
     int n,
     int& info,
     matinvreport& rep)
{
    ap::real_1d_array work;
    int i;
    int j;
    int k;
    double v;

    info = 1;
    
    //
    // Quick return if possible
    //
    if( n==0 )
    {
        info = -1;
        return;
    }
    for(i = 0; i <= n-1; i++)
    {
        if( pivots(i)>n-1||pivots(i)<i )
        {
            info = -1;
            return;
        }
    }
    
    //
    // calculate condition numbers
    //
    rep.r1 = rmatrixlurcond1(a, n);
    rep.rinf = rmatrixlurcondinf(a, n);
    if( ap::fp_less(rep.r1,rcondthreshold())||ap::fp_less(rep.rinf,rcondthreshold()) )
    {
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                a(i,j) = 0;
            }
        }
        rep.r1 = 0;
        rep.rinf = 0;
        info = -3;
        return;
    }
    
    //
    // Call cache-oblivious code
    //
    work.setlength(n);
    rmatrixluinverserec(a, 0, n, work, info, rep);
    
    //
    // apply permutations
    //
    for(i = 0; i <= n-1; i++)
    {
        for(j = n-2; j >= 0; j--)
        {
            k = pivots(j);
            v = a(i,j);
            a(i,j) = a(i,k);
            a(i,k) = v;
        }
    }
}


/*************************************************************************
Inversion of a general matrix.

Input parameters:
    A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N   -   size of matrix A.

Output parameters:
    Info    -   return code, same as in RMatrixLUInverse
    Rep     -   solver report, same as in RMatrixLUInverse
    A       -   inverse of matrix A, same as in RMatrixLUInverse

Result:
    True, if the matrix is not singular.
    False, if the matrix is singular.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************/
void rmatrixinverse(ap::real_2d_array& a,
     int n,
     int& info,
     matinvreport& rep)
{
    ap::integer_1d_array pivots;

    rmatrixlu(a, n, n, pivots);
    rmatrixluinverse(a, pivots, n, info, rep);
}


/*************************************************************************
Inversion of a matrix given by its LU decomposition.

INPUT PARAMETERS:
    A       -   LU decomposition of the matrix (output of CMatrixLU subroutine).
    Pivots  -   table of permutations which were made during the LU decomposition
                (the output of CMatrixLU subroutine).
    N       -   size of matrix A.

OUTPUT PARAMETERS:
    Info    -   return code, same as in RMatrixLUInverse
    Rep     -   solver report, same as in RMatrixLUInverse
    A       -   inverse of matrix A, same as in RMatrixLUInverse

  -- ALGLIB routine --
     05.02.2010
     Bochkanov Sergey
*************************************************************************/
void cmatrixluinverse(ap::complex_2d_array& a,
     const ap::integer_1d_array& pivots,
     int n,
     int& info,
     matinvreport& rep)
{
    ap::complex_1d_array work;
    int i;
    int j;
    int k;
    ap::complex v;

    info = 1;
    
    //
    // Quick return if possible
    //
    if( n==0 )
    {
        info = -1;
        return;
    }
    for(i = 0; i <= n-1; i++)
    {
        if( pivots(i)>n-1||pivots(i)<i )
        {
            info = -1;
            return;
        }
    }
    
    //
    // calculate condition numbers
    //
    rep.r1 = cmatrixlurcond1(a, n);
    rep.rinf = cmatrixlurcondinf(a, n);
    if( ap::fp_less(rep.r1,rcondthreshold())||ap::fp_less(rep.rinf,rcondthreshold()) )
    {
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                a(i,j) = 0;
            }
        }
        rep.r1 = 0;
        rep.rinf = 0;
        info = -3;
        return;
    }
    
    //
    // Call cache-oblivious code
    //
    work.setlength(n);
    cmatrixluinverserec(a, 0, n, work, info, rep);
    
    //
    // apply permutations
    //
    for(i = 0; i <= n-1; i++)
    {
        for(j = n-2; j >= 0; j--)
        {
            k = pivots(j);
            v = a(i,j);
            a(i,j) = a(i,k);
            a(i,k) = v;
        }
    }
}


/*************************************************************************
Inversion of a general matrix.

Input parameters:
    A   -   matrix, array[0..N-1,0..N-1].
    N   -   size of A.

Output parameters:
    Info    -   return code, same as in RMatrixLUInverse
    Rep     -   solver report, same as in RMatrixLUInverse
    A       -   inverse of matrix A, same as in RMatrixLUInverse

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************/
void cmatrixinverse(ap::complex_2d_array& a,
     int n,
     int& info,
     matinvreport& rep)
{
    ap::integer_1d_array pivots;

    cmatrixlu(a, n, n, pivots);
    cmatrixluinverse(a, pivots, n, info, rep);
}


/*************************************************************************
Inversion of a symmetric positive definite matrix which is given
by Cholesky decomposition.

Input parameters:
    A       -   Cholesky decomposition of the matrix to be inverted:
                A=U’*U or A = L*L'.
                Output of  SPDMatrixCholesky subroutine.
    N       -   size of matrix A.
    IsUpper –   storage format.
                If IsUpper = True, then matrix A is given as A = U'*U
                (matrix contains upper triangle).
                Similarly, if IsUpper = False, then A = L*L'.

Output parameters:
    Info    -   return code, same as in RMatrixLUInverse
    Rep     -   solver report, same as in RMatrixLUInverse
    A       -   inverse of matrix A, same as in RMatrixLUInverse

  -- ALGLIB routine --
     10.02.2010
     Bochkanov Sergey
*************************************************************************/
void spdmatrixcholeskyinverse(ap::real_2d_array& a,
     int n,
     bool isupper,
     int& info,
     matinvreport& rep)
{
    int i;
    int j;
    int k;
    double v;
    double ajj;
    double aii;
    ap::real_1d_array tmp;
    int info2;
    matinvreport rep2;

    if( n<1 )
    {
        info = -1;
        return;
    }
    info = 1;
    
    //
    // calculate condition numbers
    //
    rep.r1 = spdmatrixcholeskyrcond(a, n, isupper);
    rep.rinf = rep.r1;
    if( ap::fp_less(rep.r1,rcondthreshold())||ap::fp_less(rep.rinf,rcondthreshold()) )
    {
        if( isupper )
        {
            for(i = 0; i <= n-1; i++)
            {
                for(j = i; j <= n-1; j++)
                {
                    a(i,j) = 0;
                }
            }
        }
        else
        {
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= i; j++)
                {
                    a(i,j) = 0;
                }
            }
        }
        rep.r1 = 0;
        rep.rinf = 0;
        info = -3;
        return;
    }
    
    //
    // Inverse
    //
    tmp.setlength(n);
    spdmatrixcholeskyinverserec(a, 0, n, isupper, tmp);
}


/*************************************************************************
Inversion of a symmetric positive definite matrix.

Given an upper or lower triangle of a symmetric positive definite matrix,
the algorithm generates matrix A^-1 and saves the upper or lower triangle
depending on the input.

Input parameters:
    A       -   matrix to be inverted (upper or lower triangle).
                Array with elements [0..N-1,0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format.
                If IsUpper = True, then the upper triangle of matrix A is
                given, otherwise the lower triangle is given.

Output parameters:
    Info    -   return code, same as in RMatrixLUInverse
    Rep     -   solver report, same as in RMatrixLUInverse
    A       -   inverse of matrix A, same as in RMatrixLUInverse

  -- ALGLIB routine --
     10.02.2010
     Bochkanov Sergey
*************************************************************************/
void spdmatrixinverse(ap::real_2d_array& a,
     int n,
     bool isupper,
     int& info,
     matinvreport& rep)
{

    if( n<1 )
    {
        info = -1;
        return;
    }
    info = 1;
    if( spdmatrixcholesky(a, n, isupper) )
    {
        spdmatrixcholeskyinverse(a, n, isupper, info, rep);
    }
    else
    {
        info = -3;
    }
}


/*************************************************************************
Inversion of a Hermitian positive definite matrix which is given
by Cholesky decomposition.

Input parameters:
    A       -   Cholesky decomposition of the matrix to be inverted:
                A=U’*U or A = L*L'.
                Output of  HPDMatrixCholesky subroutine.
    N       -   size of matrix A.
    IsUpper –   storage format.
                If IsUpper = True, then matrix A is given as A = U'*U
                (matrix contains upper triangle).
                Similarly, if IsUpper = False, then A = L*L'.

Output parameters:
    Info    -   return code, same as in RMatrixLUInverse
    Rep     -   solver report, same as in RMatrixLUInverse
    A       -   inverse of matrix A, same as in RMatrixLUInverse

  -- ALGLIB routine --
     10.02.2010
     Bochkanov Sergey
*************************************************************************/
void hpdmatrixcholeskyinverse(ap::complex_2d_array& a,
     int n,
     bool isupper,
     int& info,
     matinvreport& rep)
{
    int i;
    int j;
    int info2;
    matinvreport rep2;
    ap::complex_1d_array tmp;
    ap::complex v;

    if( n<1 )
    {
        info = -1;
        return;
    }
    info = 1;
    
    //
    // calculate condition numbers
    //
    rep.r1 = hpdmatrixcholeskyrcond(a, n, isupper);
    rep.rinf = rep.r1;
    if( ap::fp_less(rep.r1,rcondthreshold())||ap::fp_less(rep.rinf,rcondthreshold()) )
    {
        if( isupper )
        {
            for(i = 0; i <= n-1; i++)
            {
                for(j = i; j <= n-1; j++)
                {
                    a(i,j) = 0;
                }
            }
        }
        else
        {
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= i; j++)
                {
                    a(i,j) = 0;
                }
            }
        }
        rep.r1 = 0;
        rep.rinf = 0;
        info = -3;
        return;
    }
    
    //
    // Inverse
    //
    tmp.setlength(n);
    hpdmatrixcholeskyinverserec(a, 0, n, isupper, tmp);
}


/*************************************************************************
Inversion of a Hermitian positive definite matrix.

Given an upper or lower triangle of a Hermitian positive definite matrix,
the algorithm generates matrix A^-1 and saves the upper or lower triangle
depending on the input.

Input parameters:
    A       -   matrix to be inverted (upper or lower triangle).
                Array with elements [0..N-1,0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format.
                If IsUpper = True, then the upper triangle of matrix A is
                given, otherwise the lower triangle is given.

Output parameters:
    Info    -   return code, same as in RMatrixLUInverse
    Rep     -   solver report, same as in RMatrixLUInverse
    A       -   inverse of matrix A, same as in RMatrixLUInverse

  -- ALGLIB routine --
     10.02.2010
     Bochkanov Sergey
*************************************************************************/
void hpdmatrixinverse(ap::complex_2d_array& a,
     int n,
     bool isupper,
     int& info,
     matinvreport& rep)
{

    if( n<1 )
    {
        info = -1;
        return;
    }
    info = 1;
    if( hpdmatrixcholesky(a, n, isupper) )
    {
        hpdmatrixcholeskyinverse(a, n, isupper, info, rep);
    }
    else
    {
        info = -3;
    }
}


/*************************************************************************
Triangular matrix inverse (real)

The subroutine inverts the following types of matrices:
    * upper triangular
    * upper triangular with unit diagonal
    * lower triangular
    * lower triangular with unit diagonal

In case of an upper (lower) triangular matrix,  the  inverse  matrix  will
also be upper (lower) triangular, and after the end of the algorithm,  the
inverse matrix replaces the source matrix. The elements  below (above) the
main diagonal are not changed by the algorithm.

If  the matrix  has a unit diagonal, the inverse matrix also  has  a  unit
diagonal, and the diagonal elements are not passed to the algorithm.

Input parameters:
    A       -   matrix, array[0..N-1, 0..N-1].
    N       -   size of A.
    IsUpper -   True, if the matrix is upper triangular.
    IsUnit  -   True, if the matrix has a unit diagonal.

Output parameters:
    Info    -   same as for RMatrixLUInverse
    Rep     -   same as for RMatrixLUInverse
    A       -   same as for RMatrixLUInverse.

  -- ALGLIB --
     Copyright 05.02.2010 by Bochkanov Sergey
*************************************************************************/
void rmatrixtrinverse(ap::real_2d_array& a,
     int n,
     bool isupper,
     bool isunit,
     int& info,
     matinvreport& rep)
{
    int i;
    int j;
    ap::real_1d_array tmp;

    if( n<1 )
    {
        info = -1;
        return;
    }
    info = 1;
    
    //
    // calculate condition numbers
    //
    rep.r1 = rmatrixtrrcond1(a, n, isupper, isunit);
    rep.rinf = rmatrixtrrcondinf(a, n, isupper, isunit);
    if( ap::fp_less(rep.r1,rcondthreshold())||ap::fp_less(rep.rinf,rcondthreshold()) )
    {
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                a(i,j) = 0;
            }
        }
        rep.r1 = 0;
        rep.rinf = 0;
        info = -3;
        return;
    }
    
    //
    // Invert
    //
    tmp.setlength(n);
    rmatrixtrinverserec(a, 0, n, isupper, isunit, tmp, info, rep);
}


/*************************************************************************
Triangular matrix inverse (complex)

The subroutine inverts the following types of matrices:
    * upper triangular
    * upper triangular with unit diagonal
    * lower triangular
    * lower triangular with unit diagonal

In case of an upper (lower) triangular matrix,  the  inverse  matrix  will
also be upper (lower) triangular, and after the end of the algorithm,  the
inverse matrix replaces the source matrix. The elements  below (above) the
main diagonal are not changed by the algorithm.

If  the matrix  has a unit diagonal, the inverse matrix also  has  a  unit
diagonal, and the diagonal elements are not passed to the algorithm.

Input parameters:
    A       -   matrix, array[0..N-1, 0..N-1].
    N       -   size of A.
    IsUpper -   True, if the matrix is upper triangular.
    IsUnit  -   True, if the matrix has a unit diagonal.

Output parameters:
    Info    -   same as for RMatrixLUInverse
    Rep     -   same as for RMatrixLUInverse
    A       -   same as for RMatrixLUInverse.

  -- ALGLIB --
     Copyright 05.02.2010 by Bochkanov Sergey
*************************************************************************/
void cmatrixtrinverse(ap::complex_2d_array& a,
     int n,
     bool isupper,
     bool isunit,
     int& info,
     matinvreport& rep)
{
    int i;
    int j;
    ap::complex_1d_array tmp;

    if( n<1 )
    {
        info = -1;
        return;
    }
    info = 1;
    
    //
    // calculate condition numbers
    //
    rep.r1 = cmatrixtrrcond1(a, n, isupper, isunit);
    rep.rinf = cmatrixtrrcondinf(a, n, isupper, isunit);
    if( ap::fp_less(rep.r1,rcondthreshold())||ap::fp_less(rep.rinf,rcondthreshold()) )
    {
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                a(i,j) = 0;
            }
        }
        rep.r1 = 0;
        rep.rinf = 0;
        info = -3;
        return;
    }
    
    //
    // Invert
    //
    tmp.setlength(n);
    cmatrixtrinverserec(a, 0, n, isupper, isunit, tmp, info, rep);
}


/*************************************************************************
Triangular matrix inversion, recursive subroutine

  -- ALGLIB --
     05.02.2010, Bochkanov Sergey.
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992.
*************************************************************************/
static void rmatrixtrinverserec(ap::real_2d_array& a,
     int offs,
     int n,
     bool isupper,
     bool isunit,
     ap::real_1d_array& tmp,
     int& info,
     matinvreport& rep)
{
    int n1;
    int n2;
    int i;
    int j;
    double v;
    double ajj;

    if( n<1 )
    {
        info = -1;
        return;
    }
    
    //
    // Base case
    //
    if( n<=ablasblocksize(a) )
    {
        if( isupper )
        {
            
            //
            // Compute inverse of upper triangular matrix.
            //
            for(j = 0; j <= n-1; j++)
            {
                if( !isunit )
                {
                    if( ap::fp_eq(a(offs+j,offs+j),0) )
                    {
                        info = -3;
                        return;
                    }
                    a(offs+j,offs+j) = 1/a(offs+j,offs+j);
                    ajj = -a(offs+j,offs+j);
                }
                else
                {
                    ajj = -1;
                }
                
                //
                // Compute elements 1:j-1 of j-th column.
                //
                if( j>0 )
                {
                    ap::vmove(&tmp(0), 1, &a(offs+0, offs+j), a.getstride(), ap::vlen(0,j-1));
                    for(i = 0; i <= j-1; i++)
                    {
                        if( i<j-1 )
                        {
                            v = ap::vdotproduct(&a(offs+i, offs+i+1), 1, &tmp(i+1), 1, ap::vlen(offs+i+1,offs+j-1));
                        }
                        else
                        {
                            v = 0;
                        }
                        if( !isunit )
                        {
                            a(offs+i,offs+j) = v+a(offs+i,offs+i)*tmp(i);
                        }
                        else
                        {
                            a(offs+i,offs+j) = v+tmp(i);
                        }
                    }
                    ap::vmul(&a(offs+0, offs+j), a.getstride(), ap::vlen(offs+0,offs+j-1), ajj);
                }
            }
        }
        else
        {
            
            //
            // Compute inverse of lower triangular matrix.
            //
            for(j = n-1; j >= 0; j--)
            {
                if( !isunit )
                {
                    if( ap::fp_eq(a(offs+j,offs+j),0) )
                    {
                        info = -3;
                        return;
                    }
                    a(offs+j,offs+j) = 1/a(offs+j,offs+j);
                    ajj = -a(offs+j,offs+j);
                }
                else
                {
                    ajj = -1;
                }
                if( j<n-1 )
                {
                    
                    //
                    // Compute elements j+1:n of j-th column.
                    //
                    ap::vmove(&tmp(j+1), 1, &a(offs+j+1, offs+j), a.getstride(), ap::vlen(j+1,n-1));
                    for(i = j+1; i <= n-1; i++)
                    {
                        if( i>j+1 )
                        {
                            v = ap::vdotproduct(&a(offs+i, offs+j+1), 1, &tmp(j+1), 1, ap::vlen(offs+j+1,offs+i-1));
                        }
                        else
                        {
                            v = 0;
                        }
                        if( !isunit )
                        {
                            a(offs+i,offs+j) = v+a(offs+i,offs+i)*tmp(i);
                        }
                        else
                        {
                            a(offs+i,offs+j) = v+tmp(i);
                        }
                    }
                    ap::vmul(&a(offs+j+1, offs+j), a.getstride(), ap::vlen(offs+j+1,offs+n-1), ajj);
                }
            }
        }
        return;
    }
    
    //
    // Recursive case
    //
    ablassplitlength(a, n, n1, n2);
    if( n2>0 )
    {
        if( isupper )
        {
            for(i = 0; i <= n1-1; i++)
            {
                ap::vmul(&a(offs+i, offs+n1), 1, ap::vlen(offs+n1,offs+n-1), -1);
            }
            rmatrixlefttrsm(n1, n2, a, offs, offs, isupper, isunit, 0, a, offs, offs+n1);
            rmatrixrighttrsm(n1, n2, a, offs+n1, offs+n1, isupper, isunit, 0, a, offs, offs+n1);
        }
        else
        {
            for(i = 0; i <= n2-1; i++)
            {
                ap::vmul(&a(offs+n1+i, offs), 1, ap::vlen(offs,offs+n1-1), -1);
            }
            rmatrixrighttrsm(n2, n1, a, offs, offs, isupper, isunit, 0, a, offs+n1, offs);
            rmatrixlefttrsm(n2, n1, a, offs+n1, offs+n1, isupper, isunit, 0, a, offs+n1, offs);
        }
        rmatrixtrinverserec(a, offs+n1, n2, isupper, isunit, tmp, info, rep);
    }
    rmatrixtrinverserec(a, offs, n1, isupper, isunit, tmp, info, rep);
}


/*************************************************************************
Triangular matrix inversion, recursive subroutine

  -- ALGLIB --
     05.02.2010, Bochkanov Sergey.
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992.
*************************************************************************/
static void cmatrixtrinverserec(ap::complex_2d_array& a,
     int offs,
     int n,
     bool isupper,
     bool isunit,
     ap::complex_1d_array& tmp,
     int& info,
     matinvreport& rep)
{
    int n1;
    int n2;
    int i;
    int j;
    ap::complex v;
    ap::complex ajj;

    if( n<1 )
    {
        info = -1;
        return;
    }
    
    //
    // Base case
    //
    if( n<=ablascomplexblocksize(a) )
    {
        if( isupper )
        {
            
            //
            // Compute inverse of upper triangular matrix.
            //
            for(j = 0; j <= n-1; j++)
            {
                if( !isunit )
                {
                    if( a(offs+j,offs+j)==0 )
                    {
                        info = -3;
                        return;
                    }
                    a(offs+j,offs+j) = 1/a(offs+j,offs+j);
                    ajj = -a(offs+j,offs+j);
                }
                else
                {
                    ajj = -1;
                }
                
                //
                // Compute elements 1:j-1 of j-th column.
                //
                if( j>0 )
                {
                    ap::vmove(&tmp(0), 1, &a(offs+0, offs+j), a.getstride(), "N", ap::vlen(0,j-1));
                    for(i = 0; i <= j-1; i++)
                    {
                        if( i<j-1 )
                        {
                            v = ap::vdotproduct(&a(offs+i, offs+i+1), 1, "N", &tmp(i+1), 1, "N", ap::vlen(offs+i+1,offs+j-1));
                        }
                        else
                        {
                            v = 0;
                        }
                        if( !isunit )
                        {
                            a(offs+i,offs+j) = v+a(offs+i,offs+i)*tmp(i);
                        }
                        else
                        {
                            a(offs+i,offs+j) = v+tmp(i);
                        }
                    }
                    ap::vmul(&a(offs+0, offs+j), a.getstride(), ap::vlen(offs+0,offs+j-1), ajj);
                }
            }
        }
        else
        {
            
            //
            // Compute inverse of lower triangular matrix.
            //
            for(j = n-1; j >= 0; j--)
            {
                if( !isunit )
                {
                    if( a(offs+j,offs+j)==0 )
                    {
                        info = -3;
                        return;
                    }
                    a(offs+j,offs+j) = 1/a(offs+j,offs+j);
                    ajj = -a(offs+j,offs+j);
                }
                else
                {
                    ajj = -1;
                }
                if( j<n-1 )
                {
                    
                    //
                    // Compute elements j+1:n of j-th column.
                    //
                    ap::vmove(&tmp(j+1), 1, &a(offs+j+1, offs+j), a.getstride(), "N", ap::vlen(j+1,n-1));
                    for(i = j+1; i <= n-1; i++)
                    {
                        if( i>j+1 )
                        {
                            v = ap::vdotproduct(&a(offs+i, offs+j+1), 1, "N", &tmp(j+1), 1, "N", ap::vlen(offs+j+1,offs+i-1));
                        }
                        else
                        {
                            v = 0;
                        }
                        if( !isunit )
                        {
                            a(offs+i,offs+j) = v+a(offs+i,offs+i)*tmp(i);
                        }
                        else
                        {
                            a(offs+i,offs+j) = v+tmp(i);
                        }
                    }
                    ap::vmul(&a(offs+j+1, offs+j), a.getstride(), ap::vlen(offs+j+1,offs+n-1), ajj);
                }
            }
        }
        return;
    }
    
    //
    // Recursive case
    //
    ablascomplexsplitlength(a, n, n1, n2);
    if( n2>0 )
    {
        if( isupper )
        {
            for(i = 0; i <= n1-1; i++)
            {
                ap::vmul(&a(offs+i, offs+n1), 1, ap::vlen(offs+n1,offs+n-1), -1);
            }
            cmatrixlefttrsm(n1, n2, a, offs, offs, isupper, isunit, 0, a, offs, offs+n1);
            cmatrixrighttrsm(n1, n2, a, offs+n1, offs+n1, isupper, isunit, 0, a, offs, offs+n1);
        }
        else
        {
            for(i = 0; i <= n2-1; i++)
            {
                ap::vmul(&a(offs+n1+i, offs), 1, ap::vlen(offs,offs+n1-1), -1);
            }
            cmatrixrighttrsm(n2, n1, a, offs, offs, isupper, isunit, 0, a, offs+n1, offs);
            cmatrixlefttrsm(n2, n1, a, offs+n1, offs+n1, isupper, isunit, 0, a, offs+n1, offs);
        }
        cmatrixtrinverserec(a, offs+n1, n2, isupper, isunit, tmp, info, rep);
    }
    cmatrixtrinverserec(a, offs, n1, isupper, isunit, tmp, info, rep);
}


static void rmatrixluinverserec(ap::real_2d_array& a,
     int offs,
     int n,
     ap::real_1d_array& work,
     int& info,
     matinvreport& rep)
{
    int i;
    int iws;
    int j;
    int jb;
    int jj;
    int jp;
    int k;
    double v;
    int n1;
    int n2;

    if( n<1 )
    {
        info = -1;
        return;
    }
    
    //
    // Base case
    //
    if( n<=ablasblocksize(a) )
    {
        
        //
        // Form inv(U)
        //
        rmatrixtrinverserec(a, offs, n, true, false, work, info, rep);
        if( info<=0 )
        {
            return;
        }
        
        //
        // Solve the equation inv(A)*L = inv(U) for inv(A).
        //
        for(j = n-1; j >= 0; j--)
        {
            
            //
            // Copy current column of L to WORK and replace with zeros.
            //
            for(i = j+1; i <= n-1; i++)
            {
                work(i) = a(offs+i,offs+j);
                a(offs+i,offs+j) = 0;
            }
            
            //
            // Compute current column of inv(A).
            //
            if( j<n-1 )
            {
                for(i = 0; i <= n-1; i++)
                {
                    v = ap::vdotproduct(&a(offs+i, offs+j+1), 1, &work(j+1), 1, ap::vlen(offs+j+1,offs+n-1));
                    a(offs+i,offs+j) = a(offs+i,offs+j)-v;
                }
            }
        }
        return;
    }
    
    //
    // Recursive code:
    //
    //         ( L1      )   ( U1  U12 )
    // A    =  (         ) * (         )
    //         ( L12  L2 )   (     U2  )
    //
    //         ( W   X )
    // A^-1 =  (       )
    //         ( Y   Z )
    //
    ablassplitlength(a, n, n1, n2);
    ap::ap_error::make_assertion(n2>0, "LUInverseRec: internal error!");
    
    //
    // X := inv(U1)*U12*inv(U2)
    //
    rmatrixlefttrsm(n1, n2, a, offs, offs, true, false, 0, a, offs, offs+n1);
    rmatrixrighttrsm(n1, n2, a, offs+n1, offs+n1, true, false, 0, a, offs, offs+n1);
    
    //
    // Y := inv(L2)*L12*inv(L1)
    //
    rmatrixlefttrsm(n2, n1, a, offs+n1, offs+n1, false, true, 0, a, offs+n1, offs);
    rmatrixrighttrsm(n2, n1, a, offs, offs, false, true, 0, a, offs+n1, offs);
    
    //
    // W := inv(L1*U1)+X*Y
    //
    rmatrixluinverserec(a, offs, n1, work, info, rep);
    if( info<=0 )
    {
        return;
    }
    rmatrixgemm(n1, n1, n2, 1.0, a, offs, offs+n1, 0, a, offs+n1, offs, 0, 1.0, a, offs, offs);
    
    //
    // X := -X*inv(L2)
    // Y := -inv(U2)*Y
    //
    rmatrixrighttrsm(n1, n2, a, offs+n1, offs+n1, false, true, 0, a, offs, offs+n1);
    for(i = 0; i <= n1-1; i++)
    {
        ap::vmul(&a(offs+i, offs+n1), 1, ap::vlen(offs+n1,offs+n-1), -1);
    }
    rmatrixlefttrsm(n2, n1, a, offs+n1, offs+n1, true, false, 0, a, offs+n1, offs);
    for(i = 0; i <= n2-1; i++)
    {
        ap::vmul(&a(offs+n1+i, offs), 1, ap::vlen(offs,offs+n1-1), -1);
    }
    
    //
    // Z := inv(L2*U2)
    //
    rmatrixluinverserec(a, offs+n1, n2, work, info, rep);
}


static void cmatrixluinverserec(ap::complex_2d_array& a,
     int offs,
     int n,
     ap::complex_1d_array& work,
     int& info,
     matinvreport& rep)
{
    int i;
    int iws;
    int j;
    int jb;
    int jj;
    int jp;
    int k;
    ap::complex v;
    int n1;
    int n2;

    if( n<1 )
    {
        info = -1;
        return;
    }
    
    //
    // Base case
    //
    if( n<=ablascomplexblocksize(a) )
    {
        
        //
        // Form inv(U)
        //
        cmatrixtrinverserec(a, offs, n, true, false, work, info, rep);
        if( info<=0 )
        {
            return;
        }
        
        //
        // Solve the equation inv(A)*L = inv(U) for inv(A).
        //
        for(j = n-1; j >= 0; j--)
        {
            
            //
            // Copy current column of L to WORK and replace with zeros.
            //
            for(i = j+1; i <= n-1; i++)
            {
                work(i) = a(offs+i,offs+j);
                a(offs+i,offs+j) = 0;
            }
            
            //
            // Compute current column of inv(A).
            //
            if( j<n-1 )
            {
                for(i = 0; i <= n-1; i++)
                {
                    v = ap::vdotproduct(&a(offs+i, offs+j+1), 1, "N", &work(j+1), 1, "N", ap::vlen(offs+j+1,offs+n-1));
                    a(offs+i,offs+j) = a(offs+i,offs+j)-v;
                }
            }
        }
        return;
    }
    
    //
    // Recursive code:
    //
    //         ( L1      )   ( U1  U12 )
    // A    =  (         ) * (         )
    //         ( L12  L2 )   (     U2  )
    //
    //         ( W   X )
    // A^-1 =  (       )
    //         ( Y   Z )
    //
    ablascomplexsplitlength(a, n, n1, n2);
    ap::ap_error::make_assertion(n2>0, "LUInverseRec: internal error!");
    
    //
    // X := inv(U1)*U12*inv(U2)
    //
    cmatrixlefttrsm(n1, n2, a, offs, offs, true, false, 0, a, offs, offs+n1);
    cmatrixrighttrsm(n1, n2, a, offs+n1, offs+n1, true, false, 0, a, offs, offs+n1);
    
    //
    // Y := inv(L2)*L12*inv(L1)
    //
    cmatrixlefttrsm(n2, n1, a, offs+n1, offs+n1, false, true, 0, a, offs+n1, offs);
    cmatrixrighttrsm(n2, n1, a, offs, offs, false, true, 0, a, offs+n1, offs);
    
    //
    // W := inv(L1*U1)+X*Y
    //
    cmatrixluinverserec(a, offs, n1, work, info, rep);
    if( info<=0 )
    {
        return;
    }
    cmatrixgemm(n1, n1, n2, 1.0, a, offs, offs+n1, 0, a, offs+n1, offs, 0, 1.0, a, offs, offs);
    
    //
    // X := -X*inv(L2)
    // Y := -inv(U2)*Y
    //
    cmatrixrighttrsm(n1, n2, a, offs+n1, offs+n1, false, true, 0, a, offs, offs+n1);
    for(i = 0; i <= n1-1; i++)
    {
        ap::vmul(&a(offs+i, offs+n1), 1, ap::vlen(offs+n1,offs+n-1), -1);
    }
    cmatrixlefttrsm(n2, n1, a, offs+n1, offs+n1, true, false, 0, a, offs+n1, offs);
    for(i = 0; i <= n2-1; i++)
    {
        ap::vmul(&a(offs+n1+i, offs), 1, ap::vlen(offs,offs+n1-1), -1);
    }
    
    //
    // Z := inv(L2*U2)
    //
    cmatrixluinverserec(a, offs+n1, n2, work, info, rep);
}


/*************************************************************************
Recursive subroutine for SPD inversion.

  -- ALGLIB routine --
     10.02.2010
     Bochkanov Sergey
*************************************************************************/
static void spdmatrixcholeskyinverserec(ap::real_2d_array& a,
     int offs,
     int n,
     bool isupper,
     ap::real_1d_array& tmp)
{
    int i;
    int j;
    double v;
    int n1;
    int n2;
    int info2;
    matinvreport rep2;

    if( n<1 )
    {
        return;
    }
    
    //
    // Base case
    //
    if( n<=ablasblocksize(a) )
    {
        rmatrixtrinverserec(a, offs, n, isupper, false, tmp, info2, rep2);
        if( isupper )
        {
            
            //
            // Compute the product U * U'.
            // NOTE: we never assume that diagonal of U is real
            //
            for(i = 0; i <= n-1; i++)
            {
                if( i==0 )
                {
                    
                    //
                    // 1x1 matrix
                    //
                    a(offs+i,offs+i) = ap::sqr(a(offs+i,offs+i));
                }
                else
                {
                    
                    //
                    // (I+1)x(I+1) matrix,
                    //
                    // ( A11  A12 )   ( A11^H        )   ( A11*A11^H+A12*A12^H  A12*A22^H )
                    // (          ) * (              ) = (                                )
                    // (      A22 )   ( A12^H  A22^H )   ( A22*A12^H            A22*A22^H )
                    //
                    // A11 is IxI, A22 is 1x1.
                    //
                    ap::vmove(&tmp(0), 1, &a(offs, offs+i), a.getstride(), ap::vlen(0,i-1));
                    for(j = 0; j <= i-1; j++)
                    {
                        v = a(offs+j,offs+i);
                        ap::vadd(&a(offs+j, offs+j), 1, &tmp(j), 1, ap::vlen(offs+j,offs+i-1), v);
                    }
                    v = a(offs+i,offs+i);
                    ap::vmul(&a(offs, offs+i), a.getstride(), ap::vlen(offs,offs+i-1), v);
                    a(offs+i,offs+i) = ap::sqr(a(offs+i,offs+i));
                }
            }
        }
        else
        {
            
            //
            // Compute the product L' * L
            // NOTE: we never assume that diagonal of L is real
            //
            for(i = 0; i <= n-1; i++)
            {
                if( i==0 )
                {
                    
                    //
                    // 1x1 matrix
                    //
                    a(offs+i,offs+i) = ap::sqr(a(offs+i,offs+i));
                }
                else
                {
                    
                    //
                    // (I+1)x(I+1) matrix,
                    //
                    // ( A11^H  A21^H )   ( A11      )   ( A11^H*A11+A21^H*A21  A21^H*A22 )
                    // (              ) * (          ) = (                                )
                    // (        A22^H )   ( A21  A22 )   ( A22^H*A21            A22^H*A22 )
                    //
                    // A11 is IxI, A22 is 1x1.
                    //
                    ap::vmove(&tmp(0), 1, &a(offs+i, offs), 1, ap::vlen(0,i-1));
                    for(j = 0; j <= i-1; j++)
                    {
                        v = a(offs+i,offs+j);
                        ap::vadd(&a(offs+j, offs), 1, &tmp(0), 1, ap::vlen(offs,offs+j), v);
                    }
                    v = a(offs+i,offs+i);
                    ap::vmul(&a(offs+i, offs), 1, ap::vlen(offs,offs+i-1), v);
                    a(offs+i,offs+i) = ap::sqr(a(offs+i,offs+i));
                }
            }
        }
        return;
    }
    
    //
    // Recursive code: triangular factor inversion merged with
    // UU' or L'L multiplication
    //
    ablassplitlength(a, n, n1, n2);
    
    //
    // form off-diagonal block of trangular inverse
    //
    if( isupper )
    {
        for(i = 0; i <= n1-1; i++)
        {
            ap::vmul(&a(offs+i, offs+n1), 1, ap::vlen(offs+n1,offs+n-1), -1);
        }
        rmatrixlefttrsm(n1, n2, a, offs, offs, isupper, false, 0, a, offs, offs+n1);
        rmatrixrighttrsm(n1, n2, a, offs+n1, offs+n1, isupper, false, 0, a, offs, offs+n1);
    }
    else
    {
        for(i = 0; i <= n2-1; i++)
        {
            ap::vmul(&a(offs+n1+i, offs), 1, ap::vlen(offs,offs+n1-1), -1);
        }
        rmatrixrighttrsm(n2, n1, a, offs, offs, isupper, false, 0, a, offs+n1, offs);
        rmatrixlefttrsm(n2, n1, a, offs+n1, offs+n1, isupper, false, 0, a, offs+n1, offs);
    }
    
    //
    // invert first diagonal block
    //
    spdmatrixcholeskyinverserec(a, offs, n1, isupper, tmp);
    
    //
    // update first diagonal block with off-diagonal block,
    // update off-diagonal block
    //
    if( isupper )
    {
        rmatrixsyrk(n1, n2, 1.0, a, offs, offs+n1, 0, 1.0, a, offs, offs, isupper);
        rmatrixrighttrsm(n1, n2, a, offs+n1, offs+n1, isupper, false, 1, a, offs, offs+n1);
    }
    else
    {
        rmatrixsyrk(n1, n2, 1.0, a, offs+n1, offs, 1, 1.0, a, offs, offs, isupper);
        rmatrixlefttrsm(n2, n1, a, offs+n1, offs+n1, isupper, false, 1, a, offs+n1, offs);
    }
    
    //
    // invert second diagonal block
    //
    spdmatrixcholeskyinverserec(a, offs+n1, n2, isupper, tmp);
}


/*************************************************************************
Recursive subroutine for HPD inversion.

  -- ALGLIB routine --
     10.02.2010
     Bochkanov Sergey
*************************************************************************/
static void hpdmatrixcholeskyinverserec(ap::complex_2d_array& a,
     int offs,
     int n,
     bool isupper,
     ap::complex_1d_array& tmp)
{
    int i;
    int j;
    ap::complex v;
    int n1;
    int n2;
    int info2;
    matinvreport rep2;

    if( n<1 )
    {
        return;
    }
    
    //
    // Base case
    //
    if( n<=ablascomplexblocksize(a) )
    {
        cmatrixtrinverserec(a, offs, n, isupper, false, tmp, info2, rep2);
        if( isupper )
        {
            
            //
            // Compute the product U * U'.
            // NOTE: we never assume that diagonal of U is real
            //
            for(i = 0; i <= n-1; i++)
            {
                if( i==0 )
                {
                    
                    //
                    // 1x1 matrix
                    //
                    a(offs+i,offs+i) = ap::sqr(a(offs+i,offs+i).x)+ap::sqr(a(offs+i,offs+i).y);
                }
                else
                {
                    
                    //
                    // (I+1)x(I+1) matrix,
                    //
                    // ( A11  A12 )   ( A11^H        )   ( A11*A11^H+A12*A12^H  A12*A22^H )
                    // (          ) * (              ) = (                                )
                    // (      A22 )   ( A12^H  A22^H )   ( A22*A12^H            A22*A22^H )
                    //
                    // A11 is IxI, A22 is 1x1.
                    //
                    ap::vmove(&tmp(0), 1, &a(offs, offs+i), a.getstride(), "Conj", ap::vlen(0,i-1));
                    for(j = 0; j <= i-1; j++)
                    {
                        v = a(offs+j,offs+i);
                        ap::vadd(&a(offs+j, offs+j), 1, &tmp(j), 1, "N", ap::vlen(offs+j,offs+i-1), v);
                    }
                    v = ap::conj(a(offs+i,offs+i));
                    ap::vmul(&a(offs, offs+i), a.getstride(), ap::vlen(offs,offs+i-1), v);
                    a(offs+i,offs+i) = ap::sqr(a(offs+i,offs+i).x)+ap::sqr(a(offs+i,offs+i).y);
                }
            }
        }
        else
        {
            
            //
            // Compute the product L' * L
            // NOTE: we never assume that diagonal of L is real
            //
            for(i = 0; i <= n-1; i++)
            {
                if( i==0 )
                {
                    
                    //
                    // 1x1 matrix
                    //
                    a(offs+i,offs+i) = ap::sqr(a(offs+i,offs+i).x)+ap::sqr(a(offs+i,offs+i).y);
                }
                else
                {
                    
                    //
                    // (I+1)x(I+1) matrix,
                    //
                    // ( A11^H  A21^H )   ( A11      )   ( A11^H*A11+A21^H*A21  A21^H*A22 )
                    // (              ) * (          ) = (                                )
                    // (        A22^H )   ( A21  A22 )   ( A22^H*A21            A22^H*A22 )
                    //
                    // A11 is IxI, A22 is 1x1.
                    //
                    ap::vmove(&tmp(0), 1, &a(offs+i, offs), 1, "N", ap::vlen(0,i-1));
                    for(j = 0; j <= i-1; j++)
                    {
                        v = ap::conj(a(offs+i,offs+j));
                        ap::vadd(&a(offs+j, offs), 1, &tmp(0), 1, "N", ap::vlen(offs,offs+j), v);
                    }
                    v = ap::conj(a(offs+i,offs+i));
                    ap::vmul(&a(offs+i, offs), 1, ap::vlen(offs,offs+i-1), v);
                    a(offs+i,offs+i) = ap::sqr(a(offs+i,offs+i).x)+ap::sqr(a(offs+i,offs+i).y);
                }
            }
        }
        return;
    }
    
    //
    // Recursive code: triangular factor inversion merged with
    // UU' or L'L multiplication
    //
    ablascomplexsplitlength(a, n, n1, n2);
    
    //
    // form off-diagonal block of trangular inverse
    //
    if( isupper )
    {
        for(i = 0; i <= n1-1; i++)
        {
            ap::vmul(&a(offs+i, offs+n1), 1, ap::vlen(offs+n1,offs+n-1), -1);
        }
        cmatrixlefttrsm(n1, n2, a, offs, offs, isupper, false, 0, a, offs, offs+n1);
        cmatrixrighttrsm(n1, n2, a, offs+n1, offs+n1, isupper, false, 0, a, offs, offs+n1);
    }
    else
    {
        for(i = 0; i <= n2-1; i++)
        {
            ap::vmul(&a(offs+n1+i, offs), 1, ap::vlen(offs,offs+n1-1), -1);
        }
        cmatrixrighttrsm(n2, n1, a, offs, offs, isupper, false, 0, a, offs+n1, offs);
        cmatrixlefttrsm(n2, n1, a, offs+n1, offs+n1, isupper, false, 0, a, offs+n1, offs);
    }
    
    //
    // invert first diagonal block
    //
    hpdmatrixcholeskyinverserec(a, offs, n1, isupper, tmp);
    
    //
    // update first diagonal block with off-diagonal block,
    // update off-diagonal block
    //
    if( isupper )
    {
        cmatrixsyrk(n1, n2, 1.0, a, offs, offs+n1, 0, 1.0, a, offs, offs, isupper);
        cmatrixrighttrsm(n1, n2, a, offs+n1, offs+n1, isupper, false, 2, a, offs, offs+n1);
    }
    else
    {
        cmatrixsyrk(n1, n2, 1.0, a, offs+n1, offs, 2, 1.0, a, offs, offs, isupper);
        cmatrixlefttrsm(n2, n1, a, offs+n1, offs+n1, isupper, false, 2, a, offs+n1, offs);
    }
    
    //
    // invert second diagonal block
    //
    hpdmatrixcholeskyinverserec(a, offs+n1, n2, isupper, tmp);
}




