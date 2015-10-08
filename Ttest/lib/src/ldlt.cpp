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
#include "ldlt.h"

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
     ap::integer_1d_array& pivots)
{
    int i;
    int imax;
    int j;
    int jmax;
    int k;
    int kk;
    int kp;
    int kstep;
    double absakk;
    double alpha;
    double colmax;
    double d11;
    double d12;
    double d21;
    double d22;
    double r1;
    double rowmax;
    double t;
    double wk;
    double wkm1;
    double wkp1;
    int ii;
    int i1;
    int i2;
    double vv;
    ap::real_1d_array temp;

    pivots.setbounds(0, n-1);
    temp.setbounds(0, n-1);
    
    //
    // Initialize ALPHA for use in choosing pivot block size.
    //
    alpha = (1+sqrt(double(17)))/8;
    if( isupper )
    {
        
        //
        // Factorize A as U*D*U' using the upper triangle of A
        //
        //
        // K is the main loop index, decreasing from N to 1 in steps of
        // 1 or 2
        //
        k = n-1;
        while(k>=0)
        {
            kstep = 1;
            
            //
            // Determine rows and columns to be interchanged and whether
            // a 1-by-1 or 2-by-2 pivot block will be used
            //
            absakk = fabs(a(k,k));
            
            //
            // IMAX is the row-index of the largest off-diagonal element in
            // column K+1, and COLMAX is its absolute value
            //
            if( k>0 )
            {
                imax = 1;
                for(ii = 2; ii <= k; ii++)
                {
                    if( ap::fp_greater(fabs(a(ii-1,k)),fabs(a(imax-1,k))) )
                    {
                        imax = ii;
                    }
                }
                colmax = fabs(a(imax-1,k));
            }
            else
            {
                colmax = 0;
            }
            if( ap::fp_eq(ap::maxreal(absakk, colmax),0) )
            {
                
                //
                // Column K is zero
                //
                kp = k;
            }
            else
            {
                if( ap::fp_greater_eq(absakk,alpha*colmax) )
                {
                    
                    //
                    // no interchange, use 1-by-1 pivot block
                    //
                    kp = k;
                }
                else
                {
                    
                    //
                    // JMAX is the column-index of the largest off-diagonal
                    // element in row IMAX, and ROWMAX is its absolute value
                    //
                    jmax = imax+1;
                    for(ii = imax+2; ii <= k+1; ii++)
                    {
                        if( ap::fp_greater(fabs(a(imax-1,ii-1)),fabs(a(imax-1,jmax-1))) )
                        {
                            jmax = ii;
                        }
                    }
                    rowmax = fabs(a(imax-1,jmax-1));
                    if( imax>1 )
                    {
                        jmax = 1;
                        for(ii = 2; ii <= imax-1; ii++)
                        {
                            if( ap::fp_greater(fabs(a(ii-1,imax-1)),fabs(a(jmax-1,imax-1))) )
                            {
                                jmax = ii;
                            }
                        }
                        rowmax = ap::maxreal(rowmax, fabs(a(jmax-1,imax-1)));
                    }
                    vv = colmax/rowmax;
                    if( ap::fp_greater_eq(absakk,alpha*colmax*vv) )
                    {
                        
                        //
                        // no interchange, use 1-by-1 pivot block
                        //
                        kp = k;
                    }
                    else
                    {
                        if( ap::fp_greater_eq(fabs(a(imax-1,imax-1)),alpha*rowmax) )
                        {
                            
                            //
                            // interchange rows and columns K and IMAX, use 1-by-1
                            // pivot block
                            //
                            kp = imax-1;
                        }
                        else
                        {
                            
                            //
                            // interchange rows and columns K-1 and IMAX, use 2-by-2
                            // pivot block
                            //
                            kp = imax-1;
                            kstep = 2;
                        }
                    }
                }
                kk = k+1-kstep;
                if( kp+1!=kk+1 )
                {
                    
                    //
                    // Interchange rows and columns KK and KP+1 in the leading
                    // submatrix A(0:K,0:K)
                    //
                    ap::vmove(&temp(0), 1, &a(0, kk), a.getstride(), ap::vlen(0,kp-1));
                    ap::vmove(&a(0, kk), a.getstride(), &a(0, kp), a.getstride(), ap::vlen(0,kp-1));
                    ap::vmove(&a(0, kp), a.getstride(), &temp(0), 1, ap::vlen(0,kp-1));
                    ap::vmove(&temp(kp+1), 1, &a(kp+1, kk), a.getstride(), ap::vlen(kp+1,kk-1));
                    ap::vmove(&a(kp+1, kk), a.getstride(), &a(kp, kp+1), 1, ap::vlen(kp+1,kk-1));
                    ap::vmove(&a(kp, kp+1), 1, &temp(kp+1), 1, ap::vlen(kp+1,kk-1));
                    t = a(kk,kk);
                    a(kk,kk) = a(kp,kp);
                    a(kp,kp) = t;
                    if( kstep==2 )
                    {
                        t = a(k-1,k);
                        a(k-1,k) = a(kp,k);
                        a(kp,k) = t;
                    }
                }
                
                //
                // Update the leading submatrix
                //
                if( kstep==1 )
                {
                    
                    //
                    // 1-by-1 pivot block D(k): column k now holds
                    //
                    // W(k) = U(k)*D(k)
                    //
                    // where U(k) is the k-th column of U
                    //
                    // Perform a rank-1 update of A(1:k-1,1:k-1) as
                    //
                    // A := A - U(k)*D(k)*U(k)' = A - W(k)*1/D(k)*W(k)'
                    //
                    r1 = 1/a(k,k);
                    for(i = 0; i <= k-1; i++)
                    {
                        vv = -r1*a(i,k);
                        ap::vadd(&a(i, i), 1, &a(i, k), a.getstride(), ap::vlen(i,k-1), vv);
                    }
                    
                    //
                    // Store U(K+1) in column K+1
                    //
                    ap::vmul(&a(0, k), a.getstride(), ap::vlen(0,k-1), r1);
                }
                else
                {
                    
                    //
                    // 2-by-2 pivot block D(k): columns k and k-1 now hold
                    //
                    // ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
                    //
                    // where U(k) and U(k-1) are the k-th and (k-1)-th columns
                    // of U
                    //
                    // Perform a rank-2 update of A(1:k-2,1:k-2) as
                    //
                    // A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )'
                    //    = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )'
                    //
                    if( k>1 )
                    {
                        d12 = a(k-1,k);
                        d22 = a(k-1,k-1)/d12;
                        d11 = a(k,k)/d12;
                        t = 1/(d11*d22-1);
                        d12 = t/d12;
                        for(j = k-2; j >= 0; j--)
                        {
                            wkm1 = d12*(d11*a(j,k-1)-a(j,k));
                            wk = d12*(d22*a(j,k)-a(j,k-1));
                            ap::vsub(&a(0, j), a.getstride(), &a(0, k), a.getstride(), ap::vlen(0,j), wk);
                            ap::vsub(&a(0, j), a.getstride(), &a(0, k-1), a.getstride(), ap::vlen(0,j), wkm1);
                            a(j,k) = wk;
                            a(j,k-1) = wkm1;
                        }
                    }
                }
            }
            
            //
            // Store details of the interchanges in IPIV
            //
            if( kstep==1 )
            {
                pivots(k) = kp;
            }
            else
            {
                pivots(k) = kp-n;
                pivots(k-1) = kp-n;
            }
            
            //
            // Decrease K+1 and return to the start of the main loop
            //
            k = k-kstep;
        }
    }
    else
    {
        
        //
        // Factorize A as L*D*L' using the lower triangle of A
        //
        // K+1 is the main loop index, increasing from 1 to N in steps of
        // 1 or 2
        //
        k = 0;
        while(k<=n-1)
        {
            kstep = 1;
            
            //
            // Determine rows and columns to be interchanged and whether
            // a 1-by-1 or 2-by-2 pivot block will be used
            //
            absakk = fabs(a(k,k));
            
            //
            // IMAX is the row-index of the largest off-diagonal element in
            // column K+1, and COLMAX is its absolute value
            //
            if( k<n-1 )
            {
                imax = k+1+1;
                for(ii = k+1+2; ii <= n; ii++)
                {
                    if( ap::fp_greater(fabs(a(ii-1,k)),fabs(a(imax-1,k))) )
                    {
                        imax = ii;
                    }
                }
                colmax = fabs(a(imax-1,k));
            }
            else
            {
                colmax = 0;
            }
            if( ap::fp_eq(ap::maxreal(absakk, colmax),0) )
            {
                
                //
                // Column K+1 is zero
                //
                kp = k;
            }
            else
            {
                if( ap::fp_greater_eq(absakk,alpha*colmax) )
                {
                    
                    //
                    // no interchange, use 1-by-1 pivot block
                    //
                    kp = k;
                }
                else
                {
                    
                    //
                    // JMAX is the column-index of the largest off-diagonal
                    // element in row IMAX, and ROWMAX is its absolute value
                    //
                    jmax = k+1;
                    for(ii = k+1+1; ii <= imax-1; ii++)
                    {
                        if( ap::fp_greater(fabs(a(imax-1,ii-1)),fabs(a(imax-1,jmax-1))) )
                        {
                            jmax = ii;
                        }
                    }
                    rowmax = fabs(a(imax-1,jmax-1));
                    if( imax<n )
                    {
                        jmax = imax+1;
                        for(ii = imax+2; ii <= n; ii++)
                        {
                            if( ap::fp_greater(fabs(a(ii-1,imax-1)),fabs(a(jmax-1,imax-1))) )
                            {
                                jmax = ii;
                            }
                        }
                        rowmax = ap::maxreal(rowmax, fabs(a(jmax-1,imax-1)));
                    }
                    vv = colmax/rowmax;
                    if( ap::fp_greater_eq(absakk,alpha*colmax*vv) )
                    {
                        
                        //
                        // no interchange, use 1-by-1 pivot block
                        //
                        kp = k;
                    }
                    else
                    {
                        if( ap::fp_greater_eq(fabs(a(imax-1,imax-1)),alpha*rowmax) )
                        {
                            
                            //
                            // interchange rows and columns K+1 and IMAX, use 1-by-1
                            // pivot block
                            //
                            kp = imax-1;
                        }
                        else
                        {
                            
                            //
                            // interchange rows and columns K+1+1 and IMAX, use 2-by-2
                            // pivot block
                            //
                            kp = imax-1;
                            kstep = 2;
                        }
                    }
                }
                kk = k+kstep-1;
                if( kp!=kk )
                {
                    
                    //
                    //              Interchange rows and columns KK+1 and KP+1 in the trailing
                    //              submatrix A(K+1:n,K+1:n)
                    //
                    if( kp+1<n )
                    {
                        ap::vmove(&temp(kp+1), 1, &a(kp+1, kk), a.getstride(), ap::vlen(kp+1,n-1));
                        ap::vmove(&a(kp+1, kk), a.getstride(), &a(kp+1, kp), a.getstride(), ap::vlen(kp+1,n-1));
                        ap::vmove(&a(kp+1, kp), a.getstride(), &temp(kp+1), 1, ap::vlen(kp+1,n-1));
                    }
                    ap::vmove(&temp(kk+1), 1, &a(kk+1, kk), a.getstride(), ap::vlen(kk+1,kp-1));
                    ap::vmove(&a(kk+1, kk), a.getstride(), &a(kp, kk+1), 1, ap::vlen(kk+1,kp-1));
                    ap::vmove(&a(kp, kk+1), 1, &temp(kk+1), 1, ap::vlen(kk+1,kp-1));
                    t = a(kk,kk);
                    a(kk,kk) = a(kp,kp);
                    a(kp,kp) = t;
                    if( kstep==2 )
                    {
                        t = a(k+1,k);
                        a(k+1,k) = a(kp,k);
                        a(kp,k) = t;
                    }
                }
                
                //
                // Update the trailing submatrix
                //
                if( kstep==1 )
                {
                    
                    //
                    // 1-by-1 pivot block D(K+1): column K+1 now holds
                    //
                    // W(K+1) = L(K+1)*D(K+1)
                    //
                    // where L(K+1) is the K+1-th column of L
                    //
                    if( k+1<n )
                    {
                        
                        //
                        // Perform a rank-1 update of A(K+1+1:n,K+1+1:n) as
                        //
                        // A := A - L(K+1)*D(K+1)*L(K+1)' = A - W(K+1)*(1/D(K+1))*W(K+1)'
                        //
                        d11 = 1/a(k+1-1,k+1-1);
                        for(ii = k+1; ii <= n-1; ii++)
                        {
                            vv = -d11*a(ii,k);
                            ap::vadd(&a(ii, k+1), 1, &a(k+1, k), a.getstride(), ap::vlen(k+1,ii), vv);
                        }
                        
                        //
                        // Store L(K+1) in column K+1
                        //
                        ap::vmul(&a(k+1, k), a.getstride(), ap::vlen(k+1,n-1), d11);
                    }
                }
                else
                {
                    
                    //
                    // 2-by-2 pivot block D(K+1)
                    //
                    if( k<n-2 )
                    {
                        
                        //
                        // Perform a rank-2 update of A(K+1+2:n,K+1+2:n) as
                        //
                        // A := A - ( (A(K+1) A(K+1+1))*D(K+1)**(-1) ) * (A(K+1) A(K+1+1))'
                        //
                        // where L(K+1) and L(K+1+1) are the K+1-th and (K+1+1)-th
                        // columns of L
                        //
                        d21 = a(k+1,k);
                        d11 = a(k+1,k+1)/d21;
                        d22 = a(k,k)/d21;
                        t = 1/(d11*d22-1);
                        d21 = t/d21;
                        for(j = k+2; j <= n-1; j++)
                        {
                            wk = d21*(d11*a(j,k)-a(j,k+1));
                            wkp1 = d21*(d22*a(j,k+1)-a(j,k));
                            ap::vsub(&a(j, j), a.getstride(), &a(j, k), a.getstride(), ap::vlen(j,n-1), wk);
                            ap::vsub(&a(j, j), a.getstride(), &a(j, k+1), a.getstride(), ap::vlen(j,n-1), wkp1);
                            a(j,k) = wk;
                            a(j,k+1) = wkp1;
                        }
                    }
                }
            }
            
            //
            // Store details of the interchanges in IPIV
            //
            if( kstep==1 )
            {
                pivots(k+1-1) = kp+1-1;
            }
            else
            {
                pivots(k+1-1) = kp+1-1-n;
                pivots(k+1+1-1) = kp+1-1-n;
            }
            
            //
            // Increase K+1 and return to the start of the main loop
            //
            k = k+kstep;
        }
    }
}


void ldltdecomposition(ap::real_2d_array& a,
     int n,
     bool isupper,
     ap::integer_1d_array& pivots)
{
    int i;
    int imax;
    int j;
    int jmax;
    int k;
    int kk;
    int kp;
    int kstep;
    double absakk;
    double alpha;
    double colmax;
    double d11;
    double d12;
    double d21;
    double d22;
    double r1;
    double rowmax;
    double t;
    double wk;
    double wkm1;
    double wkp1;
    int ii;
    int i1;
    int i2;
    double vv;
    ap::real_1d_array temp;

    pivots.setbounds(1, n);
    temp.setbounds(1, n);
    
    //
    // Initialize ALPHA for use in choosing pivot block size.
    //
    alpha = (1+sqrt(double(17)))/8;
    if( isupper )
    {
        
        //
        // Factorize A as U*D*U' using the upper triangle of A
        //
        //
        // K is the main loop index, decreasing from N to 1 in steps of
        // 1 or 2
        //
        k = n;
        while(k>=1)
        {
            kstep = 1;
            
            //
            // Determine rows and columns to be interchanged and whether
            // a 1-by-1 or 2-by-2 pivot block will be used
            //
            absakk = fabs(a(k,k));
            
            //
            // IMAX is the row-index of the largest off-diagonal element in
            // column K, and COLMAX is its absolute value
            //
            if( k>1 )
            {
                imax = 1;
                for(ii = 2; ii <= k-1; ii++)
                {
                    if( ap::fp_greater(fabs(a(ii,k)),fabs(a(imax,k))) )
                    {
                        imax = ii;
                    }
                }
                colmax = fabs(a(imax,k));
            }
            else
            {
                colmax = 0;
            }
            if( ap::fp_eq(ap::maxreal(absakk, colmax),0) )
            {
                
                //
                // Column K is zero
                //
                kp = k;
            }
            else
            {
                if( ap::fp_greater_eq(absakk,alpha*colmax) )
                {
                    
                    //
                    // no interchange, use 1-by-1 pivot block
                    //
                    kp = k;
                }
                else
                {
                    
                    //
                    // JMAX is the column-index of the largest off-diagonal
                    // element in row IMAX, and ROWMAX is its absolute value
                    //
                    jmax = imax+1;
                    for(ii = imax+2; ii <= k; ii++)
                    {
                        if( ap::fp_greater(fabs(a(imax,ii)),fabs(a(imax,jmax))) )
                        {
                            jmax = ii;
                        }
                    }
                    rowmax = fabs(a(imax,jmax));
                    if( imax>1 )
                    {
                        jmax = 1;
                        for(ii = 2; ii <= imax-1; ii++)
                        {
                            if( ap::fp_greater(fabs(a(ii,imax)),fabs(a(jmax,imax))) )
                            {
                                jmax = ii;
                            }
                        }
                        rowmax = ap::maxreal(rowmax, fabs(a(jmax,imax)));
                    }
                    vv = colmax/rowmax;
                    if( ap::fp_greater_eq(absakk,alpha*colmax*vv) )
                    {
                        
                        //
                        // no interchange, use 1-by-1 pivot block
                        //
                        kp = k;
                    }
                    else
                    {
                        if( ap::fp_greater_eq(fabs(a(imax,imax)),alpha*rowmax) )
                        {
                            
                            //
                            // interchange rows and columns K and IMAX, use 1-by-1
                            // pivot block
                            //
                            kp = imax;
                        }
                        else
                        {
                            
                            //
                            // interchange rows and columns K-1 and IMAX, use 2-by-2
                            // pivot block
                            //
                            kp = imax;
                            kstep = 2;
                        }
                    }
                }
                kk = k-kstep+1;
                if( kp!=kk )
                {
                    
                    //
                    // Interchange rows and columns KK and KP in the leading
                    // submatrix A(1:k,1:k)
                    //
                    i1 = kp-1;
                    ap::vmove(&temp(1), 1, &a(1, kk), a.getstride(), ap::vlen(1,i1));
                    ap::vmove(&a(1, kk), a.getstride(), &a(1, kp), a.getstride(), ap::vlen(1,i1));
                    ap::vmove(&a(1, kp), a.getstride(), &temp(1), 1, ap::vlen(1,i1));
                    i1 = kp+1;
                    i2 = kk-1;
                    ap::vmove(&temp(i1), 1, &a(i1, kk), a.getstride(), ap::vlen(i1,i2));
                    ap::vmove(&a(i1, kk), a.getstride(), &a(kp, i1), 1, ap::vlen(i1,i2));
                    ap::vmove(&a(kp, i1), 1, &temp(i1), 1, ap::vlen(i1,i2));
                    t = a(kk,kk);
                    a(kk,kk) = a(kp,kp);
                    a(kp,kp) = t;
                    if( kstep==2 )
                    {
                        t = a(k-1,k);
                        a(k-1,k) = a(kp,k);
                        a(kp,k) = t;
                    }
                }
                
                //
                // Update the leading submatrix
                //
                if( kstep==1 )
                {
                    
                    //
                    // 1-by-1 pivot block D(k): column k now holds
                    //
                    // W(k) = U(k)*D(k)
                    //
                    // where U(k) is the k-th column of U
                    //
                    // Perform a rank-1 update of A(1:k-1,1:k-1) as
                    //
                    // A := A - U(k)*D(k)*U(k)' = A - W(k)*1/D(k)*W(k)'
                    //
                    r1 = 1/a(k,k);
                    for(i = 1; i <= k-1; i++)
                    {
                        i2 = k-1;
                        vv = -r1*a(i,k);
                        ap::vadd(&a(i, i), 1, &a(i, k), a.getstride(), ap::vlen(i,i2), vv);
                    }
                    
                    //
                    // Store U(k) in column k
                    //
                    i2 = k-1;
                    ap::vmul(&a(1, k), a.getstride(), ap::vlen(1,i2), r1);
                }
                else
                {
                    
                    //
                    // 2-by-2 pivot block D(k): columns k and k-1 now hold
                    //
                    // ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
                    //
                    // where U(k) and U(k-1) are the k-th and (k-1)-th columns
                    // of U
                    //
                    // Perform a rank-2 update of A(1:k-2,1:k-2) as
                    //
                    // A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )'
                    //    = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )'
                    //
                    if( k>2 )
                    {
                        d12 = a(k-1,k);
                        d22 = a(k-1,k-1)/d12;
                        d11 = a(k,k)/d12;
                        t = 1/(d11*d22-1);
                        d12 = t/d12;
                        for(j = k-2; j >= 1; j--)
                        {
                            wkm1 = d12*(d11*a(j,k-1)-a(j,k));
                            wk = d12*(d22*a(j,k)-a(j,k-1));
                            ap::vsub(&a(1, j), a.getstride(), &a(1, k), a.getstride(), ap::vlen(1,j), wk);
                            i1 = k-1;
                            ap::vsub(&a(1, j), a.getstride(), &a(1, i1), a.getstride(), ap::vlen(1,j), wkm1);
                            a(j,k) = wk;
                            a(j,k-1) = wkm1;
                        }
                    }
                }
            }
            
            //
            // Store details of the interchanges in IPIV
            //
            if( kstep==1 )
            {
                pivots(k) = kp;
            }
            else
            {
                pivots(k) = -kp;
                pivots(k-1) = -kp;
            }
            
            //
            // Decrease K and return to the start of the main loop
            //
            k = k-kstep;
        }
    }
    else
    {
        
        //
        // Factorize A as L*D*L' using the lower triangle of A
        //
        // K is the main loop index, increasing from 1 to N in steps of
        // 1 or 2
        //
        k = 1;
        while(k<=n)
        {
            kstep = 1;
            
            //
            // Determine rows and columns to be interchanged and whether
            // a 1-by-1 or 2-by-2 pivot block will be used
            //
            absakk = fabs(a(k,k));
            
            //
            // IMAX is the row-index of the largest off-diagonal element in
            // column K, and COLMAX is its absolute value
            //
            if( k<n )
            {
                imax = k+1;
                for(ii = k+2; ii <= n; ii++)
                {
                    if( ap::fp_greater(fabs(a(ii,k)),fabs(a(imax,k))) )
                    {
                        imax = ii;
                    }
                }
                colmax = fabs(a(imax,k));
            }
            else
            {
                colmax = 0;
            }
            if( ap::fp_eq(ap::maxreal(absakk, colmax),0) )
            {
                
                //
                // Column K is zero
                //
                kp = k;
            }
            else
            {
                if( ap::fp_greater_eq(absakk,alpha*colmax) )
                {
                    
                    //
                    // no interchange, use 1-by-1 pivot block
                    //
                    kp = k;
                }
                else
                {
                    
                    //
                    // JMAX is the column-index of the largest off-diagonal
                    // element in row IMAX, and ROWMAX is its absolute value
                    //
                    jmax = k;
                    for(ii = k+1; ii <= imax-1; ii++)
                    {
                        if( ap::fp_greater(fabs(a(imax,ii)),fabs(a(imax,jmax))) )
                        {
                            jmax = ii;
                        }
                    }
                    rowmax = fabs(a(imax,jmax));
                    if( imax<n )
                    {
                        jmax = imax+1;
                        for(ii = imax+2; ii <= n; ii++)
                        {
                            if( ap::fp_greater(fabs(a(ii,imax)),fabs(a(jmax,imax))) )
                            {
                                jmax = ii;
                            }
                        }
                        rowmax = ap::maxreal(rowmax, fabs(a(jmax,imax)));
                    }
                    vv = colmax/rowmax;
                    if( ap::fp_greater_eq(absakk,alpha*colmax*vv) )
                    {
                        
                        //
                        // no interchange, use 1-by-1 pivot block
                        //
                        kp = k;
                    }
                    else
                    {
                        if( ap::fp_greater_eq(fabs(a(imax,imax)),alpha*rowmax) )
                        {
                            
                            //
                            // interchange rows and columns K and IMAX, use 1-by-1
                            // pivot block
                            //
                            kp = imax;
                        }
                        else
                        {
                            
                            //
                            // interchange rows and columns K+1 and IMAX, use 2-by-2
                            // pivot block
                            //
                            kp = imax;
                            kstep = 2;
                        }
                    }
                }
                kk = k+kstep-1;
                if( kp!=kk )
                {
                    
                    //
                    //              Interchange rows and columns KK and KP in the trailing
                    //              submatrix A(k:n,k:n)
                    //
                    if( kp<n )
                    {
                        i1 = kp+1;
                        ap::vmove(&temp(i1), 1, &a(i1, kk), a.getstride(), ap::vlen(i1,n));
                        ap::vmove(&a(i1, kk), a.getstride(), &a(i1, kp), a.getstride(), ap::vlen(i1,n));
                        ap::vmove(&a(i1, kp), a.getstride(), &temp(i1), 1, ap::vlen(i1,n));
                    }
                    i1 = kk+1;
                    i2 = kp-1;
                    ap::vmove(&temp(i1), 1, &a(i1, kk), a.getstride(), ap::vlen(i1,i2));
                    ap::vmove(&a(i1, kk), a.getstride(), &a(kp, i1), 1, ap::vlen(i1,i2));
                    ap::vmove(&a(kp, i1), 1, &temp(i1), 1, ap::vlen(i1,i2));
                    t = a(kk,kk);
                    a(kk,kk) = a(kp,kp);
                    a(kp,kp) = t;
                    if( kstep==2 )
                    {
                        t = a(k+1,k);
                        a(k+1,k) = a(kp,k);
                        a(kp,k) = t;
                    }
                }
                
                //
                // Update the trailing submatrix
                //
                if( kstep==1 )
                {
                    
                    //
                    // 1-by-1 pivot block D(k): column k now holds
                    //
                    // W(k) = L(k)*D(k)
                    //
                    // where L(k) is the k-th column of L
                    //
                    if( k<n )
                    {
                        
                        //
                        // Perform a rank-1 update of A(k+1:n,k+1:n) as
                        //
                        // A := A - L(k)*D(k)*L(k)' = A - W(k)*(1/D(k))*W(k)'
                        //
                        d11 = 1/a(k,k);
                        for(ii = k+1; ii <= n; ii++)
                        {
                            i1 = k+1;
                            i2 = ii;
                            vv = -d11*a(ii,k);
                            ap::vadd(&a(ii, i1), 1, &a(i1, k), a.getstride(), ap::vlen(i1,i2), vv);
                        }
                        
                        //
                        // Store L(k) in column K
                        //
                        i1 = k+1;
                        ap::vmul(&a(i1, k), a.getstride(), ap::vlen(i1,n), d11);
                    }
                }
                else
                {
                    
                    //
                    // 2-by-2 pivot block D(k)
                    //
                    if( k<n-1 )
                    {
                        
                        //
                        // Perform a rank-2 update of A(k+2:n,k+2:n) as
                        //
                        // A := A - ( (A(k) A(k+1))*D(k)**(-1) ) * (A(k) A(k+1))'
                        //
                        // where L(k) and L(k+1) are the k-th and (k+1)-th
                        // columns of L
                        //
                        d21 = a(k+1,k);
                        d11 = a(k+1,k+1)/d21;
                        d22 = a(k,k)/d21;
                        t = 1/(d11*d22-1);
                        d21 = t/d21;
                        for(j = k+2; j <= n; j++)
                        {
                            wk = d21*(d11*a(j,k)-a(j,k+1));
                            wkp1 = d21*(d22*a(j,k+1)-a(j,k));
                            ii = k+1;
                            ap::vsub(&a(j, j), a.getstride(), &a(j, k), a.getstride(), ap::vlen(j,n), wk);
                            ap::vsub(&a(j, j), a.getstride(), &a(j, ii), a.getstride(), ap::vlen(j,n), wkp1);
                            a(j,k) = wk;
                            a(j,k+1) = wkp1;
                        }
                    }
                }
            }
            
            //
            // Store details of the interchanges in IPIV
            //
            if( kstep==1 )
            {
                pivots(k) = kp;
            }
            else
            {
                pivots(k) = -kp;
                pivots(k+1) = -kp;
            }
            
            //
            // Increase K and return to the start of the main loop
            //
            k = k+kstep;
        }
    }
}




