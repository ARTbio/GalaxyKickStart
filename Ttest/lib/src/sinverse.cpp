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
#include "sinverse.h"

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
     bool isupper)
{
    bool result;
    ap::real_1d_array work;
    ap::real_1d_array work2;
    int i;
    int k;
    int kp;
    int kstep;
    double ak;
    double akkp1;
    double akp1;
    double d;
    double t;
    double temp;
    int km1;
    int kp1;
    int l;
    int i1;
    int i2;
    double v;

    work.setbounds(1, n);
    work2.setbounds(1, n);
    result = true;
    
    //
    // Quick return if possible
    //
    if( n==0 )
    {
        return result;
    }
    
    //
    // Check that the diagonal matrix D is nonsingular.
    //
    for(i = 0; i <= n-1; i++)
    {
        if( pivots(i)>=0&&ap::fp_eq(a(i,i),0) )
        {
            result = false;
            return result;
        }
    }
    if( isupper )
    {
        
        //
        // Compute inv(A) from the factorization A = U*D*U'.
        //
        // K+1 is the main loop index, increasing from 1 to N in steps of
        // 1 or 2, depending on the size of the diagonal blocks.
        //
        k = 0;
        while(k<=n-1)
        {
            if( pivots(k)>=0 )
            {
                
                //
                // 1 x 1 diagonal block
                //
                // Invert the diagonal block.
                //
                a(k,k) = 1/a(k,k);
                
                //
                // Compute column K+1 of the inverse.
                //
                if( k>0 )
                {
                    ap::vmove(&work(1), 1, &a(0, k), a.getstride(), ap::vlen(1,k));
                    symmetricmatrixvectormultiply(a, isupper, 1-1, k+1-1-1, work, double(-1), work2);
                    ap::vmove(&a(0, k), a.getstride(), &work2(1), 1, ap::vlen(0,k-1));
                    v = ap::vdotproduct(&work2(1), 1, &work(1), 1, ap::vlen(1,k));
                    a(k,k) = a(k,k)-v;
                }
                kstep = 1;
            }
            else
            {
                
                //
                // 2 x 2 diagonal block
                //
                // Invert the diagonal block.
                //
                t = fabs(a(k,k+1));
                ak = a(k,k)/t;
                akp1 = a(k+1,k+1)/t;
                akkp1 = a(k,k+1)/t;
                d = t*(ak*akp1-1);
                a(k,k) = akp1/d;
                a(k+1,k+1) = ak/d;
                a(k,k+1) = -akkp1/d;
                
                //
                // Compute columns K+1 and K+1+1 of the inverse.
                //
                if( k>0 )
                {
                    ap::vmove(&work(1), 1, &a(0, k), a.getstride(), ap::vlen(1,k));
                    symmetricmatrixvectormultiply(a, isupper, 0, k-1, work, double(-1), work2);
                    ap::vmove(&a(0, k), a.getstride(), &work2(1), 1, ap::vlen(0,k-1));
                    v = ap::vdotproduct(&work(1), 1, &work2(1), 1, ap::vlen(1,k));
                    a(k,k) = a(k,k)-v;
                    v = ap::vdotproduct(&a(0, k), a.getstride(), &a(0, k+1), a.getstride(), ap::vlen(0,k-1));
                    a(k,k+1) = a(k,k+1)-v;
                    ap::vmove(&work(1), 1, &a(0, k+1), a.getstride(), ap::vlen(1,k));
                    symmetricmatrixvectormultiply(a, isupper, 0, k-1, work, double(-1), work2);
                    ap::vmove(&a(0, k+1), a.getstride(), &work2(1), 1, ap::vlen(0,k-1));
                    v = ap::vdotproduct(&work(1), 1, &work2(1), 1, ap::vlen(1,k));
                    a(k+1,k+1) = a(k+1,k+1)-v;
                }
                kstep = 2;
            }
            if( pivots(k)>=0 )
            {
                kp = pivots(k);
            }
            else
            {
                kp = n+pivots(k);
            }
            if( kp!=k )
            {
                
                //
                // Interchange rows and columns K and KP in the leading
                // submatrix
                //
                ap::vmove(&work(1), 1, &a(0, k), a.getstride(), ap::vlen(1,kp));
                ap::vmove(&a(0, k), a.getstride(), &a(0, kp), a.getstride(), ap::vlen(0,kp-1));
                ap::vmove(&a(0, kp), a.getstride(), &work(1), 1, ap::vlen(0,kp-1));
                ap::vmove(&work(1), 1, &a(kp+1, k), a.getstride(), ap::vlen(1,k-1-kp));
                ap::vmove(&a(kp+1, k), a.getstride(), &a(kp, kp+1), 1, ap::vlen(kp+1,k-1));
                ap::vmove(&a(kp, kp+1), 1, &work(1), 1, ap::vlen(kp+1,k-1));
                temp = a(k,k);
                a(k,k) = a(kp,kp);
                a(kp,kp) = temp;
                if( kstep==2 )
                {
                    temp = a(k,k+1);
                    a(k,k+1) = a(kp,k+1);
                    a(kp,k+1) = temp;
                }
            }
            k = k+kstep;
        }
    }
    else
    {
        
        //
        // Compute inv(A) from the factorization A = L*D*L'.
        //
        // K is the main loop index, increasing from 0 to N-1 in steps of
        // 1 or 2, depending on the size of the diagonal blocks.
        //
        k = n-1;
        while(k>=0)
        {
            if( pivots(k)>=0 )
            {
                
                //
                // 1 x 1 diagonal block
                //
                // Invert the diagonal block.
                //
                a(k,k) = 1/a(k,k);
                
                //
                // Compute column K+1 of the inverse.
                //
                if( k<n-1 )
                {
                    ap::vmove(&work(1), 1, &a(k+1, k), a.getstride(), ap::vlen(1,n-k-1));
                    symmetricmatrixvectormultiply(a, isupper, k+1, n-1, work, double(-1), work2);
                    ap::vmove(&a(k+1, k), a.getstride(), &work2(1), 1, ap::vlen(k+1,n-1));
                    v = ap::vdotproduct(&work(1), 1, &work2(1), 1, ap::vlen(1,n-k-1));
                    a(k,k) = a(k,k)-v;
                }
                kstep = 1;
            }
            else
            {
                
                //
                // 2 x 2 diagonal block
                //
                // Invert the diagonal block.
                //
                t = fabs(a(k,k-1));
                ak = a(k-1,k-1)/t;
                akp1 = a(k,k)/t;
                akkp1 = a(k,k-1)/t;
                d = t*(ak*akp1-1);
                a(k-1,k-1) = akp1/d;
                a(k,k) = ak/d;
                a(k,k-1) = -akkp1/d;
                
                //
                // Compute columns K+1-1 and K+1 of the inverse.
                //
                if( k<n-1 )
                {
                    ap::vmove(&work(1), 1, &a(k+1, k), a.getstride(), ap::vlen(1,n-k-1));
                    symmetricmatrixvectormultiply(a, isupper, k+1, n-1, work, double(-1), work2);
                    ap::vmove(&a(k+1, k), a.getstride(), &work2(1), 1, ap::vlen(k+1,n-1));
                    v = ap::vdotproduct(&work(1), 1, &work2(1), 1, ap::vlen(1,n-k-1));
                    a(k,k) = a(k,k)-v;
                    v = ap::vdotproduct(&a(k+1, k), a.getstride(), &a(k+1, k-1), a.getstride(), ap::vlen(k+1,n-1));
                    a(k,k-1) = a(k,k-1)-v;
                    ap::vmove(&work(1), 1, &a(k+1, k-1), a.getstride(), ap::vlen(1,n-k-1));
                    symmetricmatrixvectormultiply(a, isupper, k+1, n-1, work, double(-1), work2);
                    ap::vmove(&a(k+1, k-1), a.getstride(), &work2(1), 1, ap::vlen(k+1,n-1));
                    v = ap::vdotproduct(&work(1), 1, &work2(1), 1, ap::vlen(1,n-k-1));
                    a(k-1,k-1) = a(k-1,k-1)-v;
                }
                kstep = 2;
            }
            if( pivots(k)>=0 )
            {
                kp = pivots(k);
            }
            else
            {
                kp = pivots(k)+n;
            }
            if( kp!=k )
            {
                
                //
                // Interchange rows and columns K and KP
                //
                if( kp<n-1 )
                {
                    ap::vmove(&work(1), 1, &a(kp+1, k), a.getstride(), ap::vlen(1,n-kp-1));
                    ap::vmove(&a(kp+1, k), a.getstride(), &a(kp+1, kp), a.getstride(), ap::vlen(kp+1,n-1));
                    ap::vmove(&a(kp+1, kp), a.getstride(), &work(1), 1, ap::vlen(kp+1,n-1));
                }
                ap::vmove(&work(1), 1, &a(k+1, k), a.getstride(), ap::vlen(1,kp-k-1));
                ap::vmove(&a(k+1, k), a.getstride(), &a(kp, k+1), 1, ap::vlen(k+1,kp-1));
                ap::vmove(&a(kp, k+1), 1, &work(1), 1, ap::vlen(k+1,kp-1));
                temp = a(k,k);
                a(k,k) = a(kp,kp);
                a(kp,kp) = temp;
                if( kstep==2 )
                {
                    temp = a(k,k-1);
                    a(k,k-1) = a(kp,k-1);
                    a(kp,k-1) = temp;
                }
            }
            k = k-kstep;
        }
    }
    return result;
}


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
bool smatrixinverse(ap::real_2d_array& a, int n, bool isupper)
{
    bool result;
    ap::integer_1d_array pivots;

    smatrixldlt(a, n, isupper, pivots);
    result = smatrixldltinverse(a, pivots, n, isupper);
    return result;
}


bool inverseldlt(ap::real_2d_array& a,
     const ap::integer_1d_array& pivots,
     int n,
     bool isupper)
{
    bool result;
    ap::real_1d_array work;
    ap::real_1d_array work2;
    int i;
    int k;
    int kp;
    int kstep;
    double ak;
    double akkp1;
    double akp1;
    double d;
    double t;
    double temp;
    int km1;
    int kp1;
    int l;
    int i1;
    int i2;
    double v;

    work.setbounds(1, n);
    work2.setbounds(1, n);
    result = true;
    
    //
    // Quick return if possible
    //
    if( n==0 )
    {
        return result;
    }
    
    //
    // Check that the diagonal matrix D is nonsingular.
    //
    for(i = 1; i <= n; i++)
    {
        if( pivots(i)>0&&ap::fp_eq(a(i,i),0) )
        {
            result = false;
            return result;
        }
    }
    if( isupper )
    {
        
        //
        // Compute inv(A) from the factorization A = U*D*U'.
        //
        // K is the main loop index, increasing from 1 to N in steps of
        // 1 or 2, depending on the size of the diagonal blocks.
        //
        k = 1;
        while(k<=n)
        {
            if( pivots(k)>0 )
            {
                
                //
                // 1 x 1 diagonal block
                //
                // Invert the diagonal block.
                //
                a(k,k) = 1/a(k,k);
                
                //
                // Compute column K of the inverse.
                //
                if( k>1 )
                {
                    km1 = k-1;
                    ap::vmove(&work(1), 1, &a(1, k), a.getstride(), ap::vlen(1,km1));
                    symmetricmatrixvectormultiply(a, isupper, 1, k-1, work, double(-1), work2);
                    ap::vmove(&a(1, k), a.getstride(), &work2(1), 1, ap::vlen(1,km1));
                    v = ap::vdotproduct(&work2(1), 1, &work(1), 1, ap::vlen(1,km1));
                    a(k,k) = a(k,k)-v;
                }
                kstep = 1;
            }
            else
            {
                
                //
                // 2 x 2 diagonal block
                //
                // Invert the diagonal block.
                //
                t = fabs(a(k,k+1));
                ak = a(k,k)/t;
                akp1 = a(k+1,k+1)/t;
                akkp1 = a(k,k+1)/t;
                d = t*(ak*akp1-1);
                a(k,k) = akp1/d;
                a(k+1,k+1) = ak/d;
                a(k,k+1) = -akkp1/d;
                
                //
                // Compute columns K and K+1 of the inverse.
                //
                if( k>1 )
                {
                    km1 = k-1;
                    kp1 = k+1;
                    ap::vmove(&work(1), 1, &a(1, k), a.getstride(), ap::vlen(1,km1));
                    symmetricmatrixvectormultiply(a, isupper, 1, k-1, work, double(-1), work2);
                    ap::vmove(&a(1, k), a.getstride(), &work2(1), 1, ap::vlen(1,km1));
                    v = ap::vdotproduct(&work(1), 1, &work2(1), 1, ap::vlen(1,km1));
                    a(k,k) = a(k,k)-v;
                    v = ap::vdotproduct(&a(1, k), a.getstride(), &a(1, kp1), a.getstride(), ap::vlen(1,km1));
                    a(k,k+1) = a(k,k+1)-v;
                    ap::vmove(&work(1), 1, &a(1, kp1), a.getstride(), ap::vlen(1,km1));
                    symmetricmatrixvectormultiply(a, isupper, 1, k-1, work, double(-1), work2);
                    ap::vmove(&a(1, kp1), a.getstride(), &work2(1), 1, ap::vlen(1,km1));
                    v = ap::vdotproduct(&work(1), 1, &work2(1), 1, ap::vlen(1,km1));
                    a(k+1,k+1) = a(k+1,k+1)-v;
                }
                kstep = 2;
            }
            kp = abs(pivots(k));
            if( kp!=k )
            {
                
                //
                // Interchange rows and columns K and KP in the leading
                // submatrix A(1:k+1,1:k+1)
                //
                l = kp-1;
                ap::vmove(&work(1), 1, &a(1, k), a.getstride(), ap::vlen(1,l));
                ap::vmove(&a(1, k), a.getstride(), &a(1, kp), a.getstride(), ap::vlen(1,l));
                ap::vmove(&a(1, kp), a.getstride(), &work(1), 1, ap::vlen(1,l));
                l = k-kp-1;
                i1 = kp+1;
                i2 = k-1;
                ap::vmove(&work(1), 1, &a(i1, k), a.getstride(), ap::vlen(1,l));
                ap::vmove(&a(i1, k), a.getstride(), &a(kp, i1), 1, ap::vlen(i1,i2));
                ap::vmove(&a(kp, i1), 1, &work(1), 1, ap::vlen(i1,i2));
                temp = a(k,k);
                a(k,k) = a(kp,kp);
                a(kp,kp) = temp;
                if( kstep==2 )
                {
                    temp = a(k,k+1);
                    a(k,k+1) = a(kp,k+1);
                    a(kp,k+1) = temp;
                }
            }
            k = k+kstep;
        }
    }
    else
    {
        
        //
        // Compute inv(A) from the factorization A = L*D*L'.
        //
        // K is the main loop index, increasing from 1 to N in steps of
        // 1 or 2, depending on the size of the diagonal blocks.
        //
        k = n;
        while(k>=1)
        {
            if( pivots(k)>0 )
            {
                
                //
                // 1 x 1 diagonal block
                //
                // Invert the diagonal block.
                //
                a(k,k) = 1/a(k,k);
                
                //
                // Compute column K of the inverse.
                //
                if( k<n )
                {
                    kp1 = k+1;
                    km1 = k-1;
                    l = n-k;
                    ap::vmove(&work(1), 1, &a(kp1, k), a.getstride(), ap::vlen(1,l));
                    symmetricmatrixvectormultiply(a, isupper, k+1, n, work, double(-1), work2);
                    ap::vmove(&a(kp1, k), a.getstride(), &work2(1), 1, ap::vlen(kp1,n));
                    v = ap::vdotproduct(&work(1), 1, &work2(1), 1, ap::vlen(1,l));
                    a(k,k) = a(k,k)-v;
                }
                kstep = 1;
            }
            else
            {
                
                //
                // 2 x 2 diagonal block
                //
                // Invert the diagonal block.
                //
                t = fabs(a(k,k-1));
                ak = a(k-1,k-1)/t;
                akp1 = a(k,k)/t;
                akkp1 = a(k,k-1)/t;
                d = t*(ak*akp1-1);
                a(k-1,k-1) = akp1/d;
                a(k,k) = ak/d;
                a(k,k-1) = -akkp1/d;
                
                //
                // Compute columns K-1 and K of the inverse.
                //
                if( k<n )
                {
                    kp1 = k+1;
                    km1 = k-1;
                    l = n-k;
                    ap::vmove(&work(1), 1, &a(kp1, k), a.getstride(), ap::vlen(1,l));
                    symmetricmatrixvectormultiply(a, isupper, k+1, n, work, double(-1), work2);
                    ap::vmove(&a(kp1, k), a.getstride(), &work2(1), 1, ap::vlen(kp1,n));
                    v = ap::vdotproduct(&work(1), 1, &work2(1), 1, ap::vlen(1,l));
                    a(k,k) = a(k,k)-v;
                    v = ap::vdotproduct(&a(kp1, k), a.getstride(), &a(kp1, km1), a.getstride(), ap::vlen(kp1,n));
                    a(k,k-1) = a(k,k-1)-v;
                    ap::vmove(&work(1), 1, &a(kp1, km1), a.getstride(), ap::vlen(1,l));
                    symmetricmatrixvectormultiply(a, isupper, k+1, n, work, double(-1), work2);
                    ap::vmove(&a(kp1, km1), a.getstride(), &work2(1), 1, ap::vlen(kp1,n));
                    v = ap::vdotproduct(&work(1), 1, &work2(1), 1, ap::vlen(1,l));
                    a(k-1,k-1) = a(k-1,k-1)-v;
                }
                kstep = 2;
            }
            kp = abs(pivots(k));
            if( kp!=k )
            {
                
                //
                // Interchange rows and columns K and KP in the trailing
                // submatrix A(k-1:n,k-1:n)
                //
                if( kp<n )
                {
                    l = n-kp;
                    kp1 = kp+1;
                    ap::vmove(&work(1), 1, &a(kp1, k), a.getstride(), ap::vlen(1,l));
                    ap::vmove(&a(kp1, k), a.getstride(), &a(kp1, kp), a.getstride(), ap::vlen(kp1,n));
                    ap::vmove(&a(kp1, kp), a.getstride(), &work(1), 1, ap::vlen(kp1,n));
                }
                l = kp-k-1;
                i1 = k+1;
                i2 = kp-1;
                ap::vmove(&work(1), 1, &a(i1, k), a.getstride(), ap::vlen(1,l));
                ap::vmove(&a(i1, k), a.getstride(), &a(kp, i1), 1, ap::vlen(i1,i2));
                ap::vmove(&a(kp, i1), 1, &work(1), 1, ap::vlen(i1,i2));
                temp = a(k,k);
                a(k,k) = a(kp,kp);
                a(kp,kp) = temp;
                if( kstep==2 )
                {
                    temp = a(k,k-1);
                    a(k,k-1) = a(kp,k-1);
                    a(kp,k-1) = temp;
                }
            }
            k = k-kstep;
        }
    }
    return result;
}


bool inversesymmetricindefinite(ap::real_2d_array& a, int n, bool isupper)
{
    bool result;
    ap::integer_1d_array pivots;

    ldltdecomposition(a, n, isupper, pivots);
    result = inverseldlt(a, pivots, n, isupper);
    return result;
}




