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
#include "srcond.h"

/*************************************************************************
Condition number estimate of a symmetric matrix

The algorithm calculates a lower bound of the condition number. In this
case, the algorithm does not return a lower bound of the condition number,
but an inverse number (to avoid an overflow in case of a singular matrix).

It should be noted that 1-norm and inf-norm condition numbers of symmetric
matrices are equal, so the algorithm doesn't take into account the
differences between these types of norms.

Input parameters:
    A       -   symmetric definite matrix which is given by its upper or
                lower triangle depending on IsUpper.
                Array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format.

Result:
    1/LowerBound(cond(A))
*************************************************************************/
double smatrixrcond(const ap::real_2d_array& a, int n, bool isupper)
{
    double result;
    int i;
    int j;
    ap::real_2d_array a1;

    a1.setbounds(1, n, 1, n);
    for(i = 1; i <= n; i++)
    {
        if( isupper )
        {
            for(j = i; j <= n; j++)
            {
                a1(i,j) = a(i-1,j-1);
            }
        }
        else
        {
            for(j = 1; j <= i; j++)
            {
                a1(i,j) = a(i-1,j-1);
            }
        }
    }
    result = rcondsymmetric(a1, n, isupper);
    return result;
}


/*************************************************************************
Condition number estimate of a matrix given by LDLT-decomposition

The algorithm calculates a lower bound of the condition number. In this
case, the algorithm does not return a lower bound of the condition number,
but an inverse number (to avoid an overflow in case of a singular matrix).

It should be noted that 1-norm and inf-norm condition numbers of symmetric
matrices are equal, so the algorithm doesn't take into account the
differences between these types of norms.

Input parameters:
    L       -   LDLT-decomposition of matrix A given by the upper or lower
                triangle depending on IsUpper.
                Output of SMatrixLDLT subroutine.
    Pivots  -   table of permutations which were made during LDLT-decomposition,
                Output of SMatrixLDLT subroutine.
    N       -   size of matrix A.
    IsUpper -   storage format.

Result:
    1/LowerBound(cond(A))
*************************************************************************/
double smatrixldltrcond(const ap::real_2d_array& l,
     const ap::integer_1d_array& pivots,
     int n,
     bool isupper)
{
    double result;
    int i;
    int j;
    ap::real_2d_array l1;
    ap::integer_1d_array p1;

    l1.setbounds(1, n, 1, n);
    for(i = 1; i <= n; i++)
    {
        if( isupper )
        {
            for(j = i; j <= n; j++)
            {
                l1(i,j) = l(i-1,j-1);
            }
        }
        else
        {
            for(j = 1; j <= i; j++)
            {
                l1(i,j) = l(i-1,j-1);
            }
        }
    }
    p1.setbounds(1, n);
    for(i = 1; i <= n; i++)
    {
        if( pivots(i-1)>=0 )
        {
            p1(i) = pivots(i-1)+1;
        }
        else
        {
            p1(i) = -(pivots(i-1)+n+1);
        }
    }
    result = rcondldlt(l1, p1, n, isupper);
    return result;
}


double rcondsymmetric(ap::real_2d_array a, int n, bool isupper)
{
    double result;
    int i;
    int j;
    int im;
    int jm;
    double v;
    double nrm;
    ap::integer_1d_array pivots;

    nrm = 0;
    for(j = 1; j <= n; j++)
    {
        v = 0;
        for(i = 1; i <= n; i++)
        {
            im = i;
            jm = j;
            if( isupper&&j<i )
            {
                im = j;
                jm = i;
            }
            if( !isupper&&j>i )
            {
                im = j;
                jm = i;
            }
            v = v+fabs(a(im,jm));
        }
        nrm = ap::maxreal(nrm, v);
    }
    ldltdecomposition(a, n, isupper, pivots);
    internalldltrcond(a, pivots, n, isupper, true, nrm, v);
    result = v;
    return result;
}


double rcondldlt(const ap::real_2d_array& l,
     const ap::integer_1d_array& pivots,
     int n,
     bool isupper)
{
    double result;
    double v;

    internalldltrcond(l, pivots, n, isupper, false, double(0), v);
    result = v;
    return result;
}


void internalldltrcond(const ap::real_2d_array& l,
     const ap::integer_1d_array& pivots,
     int n,
     bool isupper,
     bool isnormprovided,
     double anorm,
     double& rcond)
{
    int i;
    int kase;
    int k;
    int km1;
    int km2;
    int kp1;
    int kp2;
    double ainvnm;
    ap::real_1d_array work0;
    ap::real_1d_array work1;
    ap::real_1d_array work2;
    ap::integer_1d_array iwork;
    double v;

    ap::ap_error::make_assertion(n>=0, "");
    
    //
    // Check that the diagonal matrix D is nonsingular.
    //
    rcond = 0;
    if( isupper )
    {
        for(i = n; i >= 1; i--)
        {
            if( pivots(i)>0&&ap::fp_eq(l(i,i),0) )
            {
                return;
            }
        }
    }
    else
    {
        for(i = 1; i <= n; i++)
        {
            if( pivots(i)>0&&ap::fp_eq(l(i,i),0) )
            {
                return;
            }
        }
    }
    
    //
    // Estimate the norm of A.
    //
    if( !isnormprovided )
    {
        kase = 0;
        anorm = 0;
        while(true)
        {
            iterativeestimate1norm(n, work1, work0, iwork, anorm, kase);
            if( kase==0 )
            {
                break;
            }
            if( isupper )
            {
                
                //
                // Multiply by U'
                //
                k = n;
                while(k>=1)
                {
                    if( pivots(k)>0 )
                    {
                        
                        //
                        // P(k)
                        //
                        v = work0(k);
                        work0(k) = work0(pivots(k));
                        work0(pivots(k)) = v;
                        
                        //
                        // U(k)
                        //
                        km1 = k-1;
                        v = ap::vdotproduct(&work0(1), 1, &l(1, k), l.getstride(), ap::vlen(1,km1));
                        work0(k) = work0(k)+v;
                        
                        //
                        // Next k
                        //
                        k = k-1;
                    }
                    else
                    {
                        
                        //
                        // P(k)
                        //
                        v = work0(k-1);
                        work0(k-1) = work0(-pivots(k-1));
                        work0(-pivots(k-1)) = v;
                        
                        //
                        // U(k)
                        //
                        km1 = k-1;
                        km2 = k-2;
                        v = ap::vdotproduct(&work0(1), 1, &l(1, km1), l.getstride(), ap::vlen(1,km2));
                        work0(km1) = work0(km1)+v;
                        v = ap::vdotproduct(&work0(1), 1, &l(1, k), l.getstride(), ap::vlen(1,km2));
                        work0(k) = work0(k)+v;
                        
                        //
                        // Next k
                        //
                        k = k-2;
                    }
                }
                
                //
                // Multiply by D
                //
                k = n;
                while(k>=1)
                {
                    if( pivots(k)>0 )
                    {
                        work0(k) = work0(k)*l(k,k);
                        k = k-1;
                    }
                    else
                    {
                        v = work0(k-1);
                        work0(k-1) = l(k-1,k-1)*work0(k-1)+l(k-1,k)*work0(k);
                        work0(k) = l(k-1,k)*v+l(k,k)*work0(k);
                        k = k-2;
                    }
                }
                
                //
                // Multiply by U
                //
                k = 1;
                while(k<=n)
                {
                    if( pivots(k)>0 )
                    {
                        
                        //
                        // U(k)
                        //
                        km1 = k-1;
                        v = work0(k);
                        ap::vadd(&work0(1), 1, &l(1, k), l.getstride(), ap::vlen(1,km1), v);
                        
                        //
                        // P(k)
                        //
                        v = work0(k);
                        work0(k) = work0(pivots(k));
                        work0(pivots(k)) = v;
                        
                        //
                        // Next k
                        //
                        k = k+1;
                    }
                    else
                    {
                        
                        //
                        // U(k)
                        //
                        km1 = k-1;
                        kp1 = k+1;
                        v = work0(k);
                        ap::vadd(&work0(1), 1, &l(1, k), l.getstride(), ap::vlen(1,km1), v);
                        v = work0(kp1);
                        ap::vadd(&work0(1), 1, &l(1, kp1), l.getstride(), ap::vlen(1,km1), v);
                        
                        //
                        // P(k)
                        //
                        v = work0(k);
                        work0(k) = work0(-pivots(k));
                        work0(-pivots(k)) = v;
                        
                        //
                        // Next k
                        //
                        k = k+2;
                    }
                }
            }
            else
            {
                
                //
                // Multiply by L'
                //
                k = 1;
                while(k<=n)
                {
                    if( pivots(k)>0 )
                    {
                        
                        //
                        // P(k)
                        //
                        v = work0(k);
                        work0(k) = work0(pivots(k));
                        work0(pivots(k)) = v;
                        
                        //
                        // L(k)
                        //
                        kp1 = k+1;
                        v = ap::vdotproduct(&work0(kp1), 1, &l(kp1, k), l.getstride(), ap::vlen(kp1,n));
                        work0(k) = work0(k)+v;
                        
                        //
                        // Next k
                        //
                        k = k+1;
                    }
                    else
                    {
                        
                        //
                        // P(k)
                        //
                        v = work0(k+1);
                        work0(k+1) = work0(-pivots(k+1));
                        work0(-pivots(k+1)) = v;
                        
                        //
                        // L(k)
                        //
                        kp1 = k+1;
                        kp2 = k+2;
                        v = ap::vdotproduct(&work0(kp2), 1, &l(kp2, k), l.getstride(), ap::vlen(kp2,n));
                        work0(k) = work0(k)+v;
                        v = ap::vdotproduct(&work0(kp2), 1, &l(kp2, kp1), l.getstride(), ap::vlen(kp2,n));
                        work0(kp1) = work0(kp1)+v;
                        
                        //
                        // Next k
                        //
                        k = k+2;
                    }
                }
                
                //
                // Multiply by D
                //
                k = n;
                while(k>=1)
                {
                    if( pivots(k)>0 )
                    {
                        work0(k) = work0(k)*l(k,k);
                        k = k-1;
                    }
                    else
                    {
                        v = work0(k-1);
                        work0(k-1) = l(k-1,k-1)*work0(k-1)+l(k,k-1)*work0(k);
                        work0(k) = l(k,k-1)*v+l(k,k)*work0(k);
                        k = k-2;
                    }
                }
                
                //
                // Multiply by L
                //
                k = n;
                while(k>=1)
                {
                    if( pivots(k)>0 )
                    {
                        
                        //
                        // L(k)
                        //
                        kp1 = k+1;
                        v = work0(k);
                        ap::vadd(&work0(kp1), 1, &l(kp1, k), l.getstride(), ap::vlen(kp1,n), v);
                        
                        //
                        // P(k)
                        //
                        v = work0(k);
                        work0(k) = work0(pivots(k));
                        work0(pivots(k)) = v;
                        
                        //
                        // Next k
                        //
                        k = k-1;
                    }
                    else
                    {
                        
                        //
                        // L(k)
                        //
                        kp1 = k+1;
                        km1 = k-1;
                        v = work0(k);
                        ap::vadd(&work0(kp1), 1, &l(kp1, k), l.getstride(), ap::vlen(kp1,n), v);
                        v = work0(km1);
                        ap::vadd(&work0(kp1), 1, &l(kp1, km1), l.getstride(), ap::vlen(kp1,n), v);
                        
                        //
                        // P(k)
                        //
                        v = work0(k);
                        work0(k) = work0(-pivots(k));
                        work0(-pivots(k)) = v;
                        
                        //
                        // Next k
                        //
                        k = k-2;
                    }
                }
            }
        }
    }
    
    //
    // Quick return if possible
    //
    rcond = 0;
    if( n==0 )
    {
        rcond = 1;
        return;
    }
    if( ap::fp_eq(anorm,0) )
    {
        return;
    }
    
    //
    // Estimate the 1-norm of inv(A).
    //
    kase = 0;
    while(true)
    {
        iterativeestimate1norm(n, work1, work0, iwork, ainvnm, kase);
        if( kase==0 )
        {
            break;
        }
        solvesystemldlt(l, pivots, work0, n, isupper, work2);
        ap::vmove(&work0(1), 1, &work2(1), 1, ap::vlen(1,n));
    }
    
    //
    // Compute the estimate of the reciprocal condition number.
    //
    if( ap::fp_neq(ainvnm,0) )
    {
        v = 1/ainvnm;
        rcond = v/anorm;
    }
}




