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
#include "ssolve.h"

/*************************************************************************
Solving  a system  of linear equations  with a system matrix  given by its
LDLT decomposition

The algorithm solves systems with a square matrix only.

Input parameters:
    A       -   LDLT decomposition of the matrix (the result of the
                SMatrixLDLT subroutine).
    Pivots  -   row permutation table (the result of the SMatrixLDLT subroutine).
    B       -   right side of a system.
                Array whose index ranges within [0..N-1].
    N       -   size of matrix A.
    IsUpper -   points to the triangle of matrix A in which the LDLT
                decomposition is stored.
                If IsUpper=True, the decomposition has the form of U*D*U',
                matrix U is stored in the upper triangle of  matrix A  (in
                that case, the lower triangle isn't used and isn't changed
                by the subroutine).
                Similarly, if IsUpper=False, the decomposition has the form
                of L*D*L' and the lower triangle stores matrix L.

Output parameters:
    X       -   solution of a system.
                Array whose index ranges within [0..N-1].

Result:
    True, if the matrix is not singular. X contains the solution.
    False, if the matrix is singular (the determinant of matrix D is equal
to 0). In this case, X doesn't contain a solution.
*************************************************************************/
bool smatrixldltsolve(const ap::real_2d_array& a,
     const ap::integer_1d_array& pivots,
     ap::real_1d_array b,
     int n,
     bool isupper,
     ap::real_1d_array& x)
{
    bool result;
    int i;
    int k;
    int kp;
    double ak;
    double akm1;
    double akm1k;
    double bk;
    double bkm1;
    double denom;
    double v;

    
    //
    // Quick return if possible
    //
    result = true;
    if( n==0 )
    {
        return result;
    }
    
    //
    // Check that the diagonal matrix D is nonsingular
    //
    if( isupper )
    {
        
        //
        // Upper triangular storage: examine D from bottom to top
        //
        for(i = n-1; i >= 0; i--)
        {
            if( pivots(i)>=0&&ap::fp_eq(a(i,i),0) )
            {
                result = false;
                return result;
            }
        }
    }
    else
    {
        
        //
        // Lower triangular storage: examine D from top to bottom.
        //
        for(i = 0; i <= n-1; i++)
        {
            if( pivots(i)>=0&&ap::fp_eq(a(i,i),0) )
            {
                result = false;
                return result;
            }
        }
    }
    
    //
    // Solve Ax = b
    //
    if( isupper )
    {
        
        //
        // Solve A*X = B, where A = U*D*U'.
        //
        // First solve U*D*X = B, overwriting B with X.
        //
        // K+1 is the main loop index, decreasing from N to 1 in steps of
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
                // Interchange rows K+1 and IPIV(K+1).
                //
                kp = pivots(k);
                if( kp!=k )
                {
                    v = b(k);
                    b(k) = b(kp);
                    b(kp) = v;
                }
                
                //
                // Multiply by inv(U(K+1)), where U(K+1) is the transformation
                // stored in column K+1 of A.
                //
                v = b(k);
                ap::vsub(&b(0), 1, &a(0, k), a.getstride(), ap::vlen(0,k-1), v);
                
                //
                // Multiply by the inverse of the diagonal block.
                //
                b(k) = b(k)/a(k,k);
                k = k-1;
            }
            else
            {
                
                //
                // 2 x 2 diagonal block
                //
                // Interchange rows K+1-1 and -IPIV(K+1).
                //
                kp = pivots(k)+n;
                if( kp!=k-1 )
                {
                    v = b(k-1);
                    b(k-1) = b(kp);
                    b(kp) = v;
                }
                
                //
                // Multiply by inv(U(K+1)), where U(K+1) is the transformation
                // stored in columns K+1-1 and K+1 of A.
                //
                v = b(k);
                ap::vsub(&b(0), 1, &a(0, k), a.getstride(), ap::vlen(0,k-2), v);
                v = b(k-1);
                ap::vsub(&b(0), 1, &a(0, k-1), a.getstride(), ap::vlen(0,k-2), v);
                
                //
                // Multiply by the inverse of the diagonal block.
                //
                akm1k = a(k-1,k);
                akm1 = a(k-1,k-1)/akm1k;
                ak = a(k,k)/akm1k;
                denom = akm1*ak-1;
                bkm1 = b(k-1)/akm1k;
                bk = b(k)/akm1k;
                b(k-1) = (ak*bkm1-bk)/denom;
                b(k) = (akm1*bk-bkm1)/denom;
                k = k-2;
            }
        }
        
        //
        // Next solve U'*X = B, overwriting B with X.
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
                // Multiply by inv(U'(K+1)), where U(K+1) is the transformation
                // stored in column K+1 of A.
                //
                v = ap::vdotproduct(&b(0), 1, &a(0, k), a.getstride(), ap::vlen(0,k-1));
                b(k) = b(k)-v;
                
                //
                // Interchange rows K+1 and IPIV(K+1).
                //
                kp = pivots(k);
                if( kp!=k )
                {
                    v = b(k);
                    b(k) = b(kp);
                    b(kp) = v;
                }
                k = k+1;
            }
            else
            {
                
                //
                // 2 x 2 diagonal block
                //
                // Multiply by inv(U'(K+1+1)), where U(K+1+1) is the transformation
                // stored in columns K+1 and K+1+1 of A.
                //
                v = ap::vdotproduct(&b(0), 1, &a(0, k), a.getstride(), ap::vlen(0,k-1));
                b(k) = b(k)-v;
                v = ap::vdotproduct(&b(0), 1, &a(0, k+1), a.getstride(), ap::vlen(0,k-1));
                b(k+1) = b(k+1)-v;
                
                //
                // Interchange rows K+1 and -IPIV(K+1).
                //
                kp = pivots(k)+n;
                if( kp!=k )
                {
                    v = b(k);
                    b(k) = b(kp);
                    b(kp) = v;
                }
                k = k+2;
            }
        }
    }
    else
    {
        
        //
        // Solve A*X = B, where A = L*D*L'.
        //
        // First solve L*D*X = B, overwriting B with X.
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
                // Interchange rows K+1 and IPIV(K+1).
                //
                kp = pivots(k);
                if( kp!=k )
                {
                    v = b(k);
                    b(k) = b(kp);
                    b(kp) = v;
                }
                
                //
                // Multiply by inv(L(K+1)), where L(K+1) is the transformation
                // stored in column K+1 of A.
                //
                if( k+1<n )
                {
                    v = b(k);
                    ap::vsub(&b(k+1), 1, &a(k+1, k), a.getstride(), ap::vlen(k+1,n-1), v);
                }
                
                //
                // Multiply by the inverse of the diagonal block.
                //
                b(k) = b(k)/a(k,k);
                k = k+1;
            }
            else
            {
                
                //
                // 2 x 2 diagonal block
                //
                // Interchange rows K+1+1 and -IPIV(K+1).
                //
                kp = pivots(k)+n;
                if( kp!=k+1 )
                {
                    v = b(k+1);
                    b(k+1) = b(kp);
                    b(kp) = v;
                }
                
                //
                // Multiply by inv(L(K+1)), where L(K+1) is the transformation
                // stored in columns K+1 and K+1+1 of A.
                //
                if( k+1<n-1 )
                {
                    v = b(k);
                    ap::vsub(&b(k+2), 1, &a(k+2, k), a.getstride(), ap::vlen(k+2,n-1), v);
                    v = b(k+1);
                    ap::vsub(&b(k+2), 1, &a(k+2, k+1), a.getstride(), ap::vlen(k+2,n-1), v);
                }
                
                //
                // Multiply by the inverse of the diagonal block.
                //
                akm1k = a(k+1,k);
                akm1 = a(k,k)/akm1k;
                ak = a(k+1,k+1)/akm1k;
                denom = akm1*ak-1;
                bkm1 = b(k)/akm1k;
                bk = b(k+1)/akm1k;
                b(k) = (ak*bkm1-bk)/denom;
                b(k+1) = (akm1*bk-bkm1)/denom;
                k = k+2;
            }
        }
        
        //
        // Next solve L'*X = B, overwriting B with X.
        //
        // K+1 is the main loop index, decreasing from N to 1 in steps of
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
                // Multiply by inv(L'(K+1)), where L(K+1) is the transformation
                // stored in column K+1 of A.
                //
                if( k+1<n )
                {
                    v = ap::vdotproduct(&b(k+1), 1, &a(k+1, k), a.getstride(), ap::vlen(k+1,n-1));
                    b(k) = b(k)-v;
                }
                
                //
                // Interchange rows K+1 and IPIV(K+1).
                //
                kp = pivots(k);
                if( kp!=k )
                {
                    v = b(k);
                    b(k) = b(kp);
                    b(kp) = v;
                }
                k = k-1;
            }
            else
            {
                
                //
                // 2 x 2 diagonal block
                //
                // Multiply by inv(L'(K+1-1)), where L(K+1-1) is the transformation
                // stored in columns K+1-1 and K+1 of A.
                //
                if( k+1<n )
                {
                    v = ap::vdotproduct(&b(k+1), 1, &a(k+1, k), a.getstride(), ap::vlen(k+1,n-1));
                    b(k) = b(k)-v;
                    v = ap::vdotproduct(&b(k+1), 1, &a(k+1, k-1), a.getstride(), ap::vlen(k+1,n-1));
                    b(k-1) = b(k-1)-v;
                }
                
                //
                // Interchange rows K+1 and -IPIV(K+1).
                //
                kp = pivots(k)+n;
                if( kp!=k )
                {
                    v = b(k);
                    b(k) = b(kp);
                    b(kp) = v;
                }
                k = k-2;
            }
        }
    }
    x.setbounds(0, n-1);
    ap::vmove(&x(0), 1, &b(0), 1, ap::vlen(0,n-1));
    return result;
}


/*************************************************************************
Solving a system of linear equations with a symmetric system matrix

Input parameters:
    A       -   system matrix (upper or lower triangle).
                Array whose indexes range within [0..N-1, 0..N-1].
    B       -   right side of a system.
                Array whose index ranges within [0..N-1].
    N       -   size of matrix A.
    IsUpper -   If IsUpper = True, A contains the upper triangle,
                otherwise A contains the lower triangle.

Output parameters:
    X       -   solution of a system.
                Array whose index ranges within [0..N-1].

Result:
    True, if the matrix is not singular. X contains the solution.
    False, if the matrix is singular (the determinant of the matrix is equal
to 0). In this case, X doesn't contain a solution.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************/
bool smatrixsolve(ap::real_2d_array a,
     const ap::real_1d_array& b,
     int n,
     bool isupper,
     ap::real_1d_array& x)
{
    bool result;
    ap::integer_1d_array pivots;

    smatrixldlt(a, n, isupper, pivots);
    result = smatrixldltsolve(a, pivots, b, n, isupper, x);
    return result;
}


bool solvesystemldlt(const ap::real_2d_array& a,
     const ap::integer_1d_array& pivots,
     ap::real_1d_array b,
     int n,
     bool isupper,
     ap::real_1d_array& x)
{
    bool result;
    int i;
    int k;
    int kp;
    int km1;
    int km2;
    int kp1;
    int kp2;
    double ak;
    double akm1;
    double akm1k;
    double bk;
    double bkm1;
    double denom;
    double v;

    
    //
    // Quick return if possible
    //
    result = true;
    if( n==0 )
    {
        return result;
    }
    
    //
    // Check that the diagonal matrix D is nonsingular
    //
    if( isupper )
    {
        
        //
        // Upper triangular storage: examine D from bottom to top
        //
        for(i = n; i >= 1; i--)
        {
            if( pivots(i)>0&&ap::fp_eq(a(i,i),0) )
            {
                result = false;
                return result;
            }
        }
    }
    else
    {
        
        //
        // Lower triangular storage: examine D from top to bottom.
        //
        for(i = 1; i <= n; i++)
        {
            if( pivots(i)>0&&ap::fp_eq(a(i,i),0) )
            {
                result = false;
                return result;
            }
        }
    }
    
    //
    // Solve Ax = b
    //
    if( isupper )
    {
        
        //
        // Solve A*X = B, where A = U*D*U'.
        //
        // First solve U*D*X = B, overwriting B with X.
        //
        // K is the main loop index, decreasing from N to 1 in steps of
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
                // Interchange rows K and IPIV(K).
                //
                kp = pivots(k);
                if( kp!=k )
                {
                    v = b(k);
                    b(k) = b(kp);
                    b(kp) = v;
                }
                
                //
                // Multiply by inv(U(K)), where U(K) is the transformation
                // stored in column K of A.
                //
                km1 = k-1;
                v = b(k);
                ap::vsub(&b(1), 1, &a(1, k), a.getstride(), ap::vlen(1,km1), v);
                
                //
                // Multiply by the inverse of the diagonal block.
                //
                b(k) = b(k)/a(k,k);
                k = k-1;
            }
            else
            {
                
                //
                // 2 x 2 diagonal block
                //
                // Interchange rows K-1 and -IPIV(K).
                //
                kp = -pivots(k);
                if( kp!=k-1 )
                {
                    v = b(k-1);
                    b(k-1) = b(kp);
                    b(kp) = v;
                }
                
                //
                // Multiply by inv(U(K)), where U(K) is the transformation
                // stored in columns K-1 and K of A.
                //
                km2 = k-2;
                km1 = k-1;
                v = b(k);
                ap::vsub(&b(1), 1, &a(1, k), a.getstride(), ap::vlen(1,km2), v);
                v = b(k-1);
                ap::vsub(&b(1), 1, &a(1, km1), a.getstride(), ap::vlen(1,km2), v);
                
                //
                // Multiply by the inverse of the diagonal block.
                //
                akm1k = a(k-1,k);
                akm1 = a(k-1,k-1)/akm1k;
                ak = a(k,k)/akm1k;
                denom = akm1*ak-1;
                bkm1 = b(k-1)/akm1k;
                bk = b(k)/akm1k;
                b(k-1) = (ak*bkm1-bk)/denom;
                b(k) = (akm1*bk-bkm1)/denom;
                k = k-2;
            }
        }
        
        //
        // Next solve U'*X = B, overwriting B with X.
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
                // Multiply by inv(U'(K)), where U(K) is the transformation
                // stored in column K of A.
                //
                km1 = k-1;
                v = ap::vdotproduct(&b(1), 1, &a(1, k), a.getstride(), ap::vlen(1,km1));
                b(k) = b(k)-v;
                
                //
                // Interchange rows K and IPIV(K).
                //
                kp = pivots(k);
                if( kp!=k )
                {
                    v = b(k);
                    b(k) = b(kp);
                    b(kp) = v;
                }
                k = k+1;
            }
            else
            {
                
                //
                // 2 x 2 diagonal block
                //
                // Multiply by inv(U'(K+1)), where U(K+1) is the transformation
                // stored in columns K and K+1 of A.
                //
                km1 = k-1;
                kp1 = k+1;
                v = ap::vdotproduct(&b(1), 1, &a(1, k), a.getstride(), ap::vlen(1,km1));
                b(k) = b(k)-v;
                v = ap::vdotproduct(&b(1), 1, &a(1, kp1), a.getstride(), ap::vlen(1,km1));
                b(k+1) = b(k+1)-v;
                
                //
                // Interchange rows K and -IPIV(K).
                //
                kp = -pivots(k);
                if( kp!=k )
                {
                    v = b(k);
                    b(k) = b(kp);
                    b(kp) = v;
                }
                k = k+2;
            }
        }
    }
    else
    {
        
        //
        // Solve A*X = B, where A = L*D*L'.
        //
        // First solve L*D*X = B, overwriting B with X.
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
                // Interchange rows K and IPIV(K).
                //
                kp = pivots(k);
                if( kp!=k )
                {
                    v = b(k);
                    b(k) = b(kp);
                    b(kp) = v;
                }
                
                //
                // Multiply by inv(L(K)), where L(K) is the transformation
                // stored in column K of A.
                //
                if( k<n )
                {
                    kp1 = k+1;
                    v = b(k);
                    ap::vsub(&b(kp1), 1, &a(kp1, k), a.getstride(), ap::vlen(kp1,n), v);
                }
                
                //
                // Multiply by the inverse of the diagonal block.
                //
                b(k) = b(k)/a(k,k);
                k = k+1;
            }
            else
            {
                
                //
                // 2 x 2 diagonal block
                //
                // Interchange rows K+1 and -IPIV(K).
                //
                kp = -pivots(k);
                if( kp!=k+1 )
                {
                    v = b(k+1);
                    b(k+1) = b(kp);
                    b(kp) = v;
                }
                
                //
                // Multiply by inv(L(K)), where L(K) is the transformation
                // stored in columns K and K+1 of A.
                //
                if( k<n-1 )
                {
                    kp1 = k+1;
                    kp2 = k+2;
                    v = b(k);
                    ap::vsub(&b(kp2), 1, &a(kp2, k), a.getstride(), ap::vlen(kp2,n), v);
                    v = b(k+1);
                    ap::vsub(&b(kp2), 1, &a(kp2, kp1), a.getstride(), ap::vlen(kp2,n), v);
                }
                
                //
                // Multiply by the inverse of the diagonal block.
                //
                akm1k = a(k+1,k);
                akm1 = a(k,k)/akm1k;
                ak = a(k+1,k+1)/akm1k;
                denom = akm1*ak-1;
                bkm1 = b(k)/akm1k;
                bk = b(k+1)/akm1k;
                b(k) = (ak*bkm1-bk)/denom;
                b(k+1) = (akm1*bk-bkm1)/denom;
                k = k+2;
            }
        }
        
        //
        // Next solve L'*X = B, overwriting B with X.
        //
        // K is the main loop index, decreasing from N to 1 in steps of
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
                // Multiply by inv(L'(K)), where L(K) is the transformation
                // stored in column K of A.
                //
                if( k<n )
                {
                    kp1 = k+1;
                    v = ap::vdotproduct(&b(kp1), 1, &a(kp1, k), a.getstride(), ap::vlen(kp1,n));
                    b(k) = b(k)-v;
                }
                
                //
                // Interchange rows K and IPIV(K).
                //
                kp = pivots(k);
                if( kp!=k )
                {
                    v = b(k);
                    b(k) = b(kp);
                    b(kp) = v;
                }
                k = k-1;
            }
            else
            {
                
                //
                // 2 x 2 diagonal block
                //
                // Multiply by inv(L'(K-1)), where L(K-1) is the transformation
                // stored in columns K-1 and K of A.
                //
                if( k<n )
                {
                    kp1 = k+1;
                    km1 = k-1;
                    v = ap::vdotproduct(&b(kp1), 1, &a(kp1, k), a.getstride(), ap::vlen(kp1,n));
                    b(k) = b(k)-v;
                    v = ap::vdotproduct(&b(kp1), 1, &a(kp1, km1), a.getstride(), ap::vlen(kp1,n));
                    b(k-1) = b(k-1)-v;
                }
                
                //
                // Interchange rows K and -IPIV(K).
                //
                kp = -pivots(k);
                if( kp!=k )
                {
                    v = b(k);
                    b(k) = b(kp);
                    b(kp) = v;
                }
                k = k-2;
            }
        }
    }
    x.setbounds(1, n);
    ap::vmove(&x(1), 1, &b(1), 1, ap::vlen(1,n));
    return result;
}


bool solvesymmetricsystem(ap::real_2d_array a,
     ap::real_1d_array b,
     int n,
     bool isupper,
     ap::real_1d_array& x)
{
    bool result;
    ap::integer_1d_array pivots;

    ldltdecomposition(a, n, isupper, pivots);
    result = solvesystemldlt(a, pivots, b, n, isupper, x);
    return result;
}




