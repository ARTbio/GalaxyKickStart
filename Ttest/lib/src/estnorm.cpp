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
#include "estnorm.h"

/*************************************************************************
Matrix norm estimation

The algorithm estimates the 1-norm of square matrix A  on  the  assumption
that the multiplication of matrix  A  by  the  vector  is  available  (the
iterative method is used). It is recommended to use this algorithm  if  it
is hard  to  calculate  matrix  elements  explicitly  (for  example,  when
estimating the inverse matrix norm).

The algorithm uses back communication for multiplying the  vector  by  the
matrix.  If  KASE=0  after  returning from a subroutine, its execution was
completed successfully, otherwise it is required to multiply the  returned
vector by matrix A and call the subroutine again.

The DemoIterativeEstimateNorm subroutine shows a simple example.

Parameters:
    N       -   size of matrix A.
    V       -   vector.   It is initialized by the subroutine on the first
                call. It is then passed into it on repeated calls.
    X       -   if KASE<>0, it contains the vector to be replaced by:
                    A * X,      if KASE=1
                    A^T * X,    if KASE=2
                Array whose index ranges within [1..N].
    ISGN    -   vector. It is initialized by the subroutine on  the  first
                call. It is then passed into it on repeated calls.
    EST     -   if KASE=0, it contains the lower boundary of the matrix
                norm estimate.
    KASE    -   on the first call, it should be equal to 0. After the last
                return, it is equal to 0 (EST contains the  matrix  norm),
                on intermediate returns it can be equal to 1 or 2 depending
                on the operation to be performed on vector X.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************/
void iterativeestimate1norm(int n,
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
        v.setbounds(1, n+3);
        x.setbounds(1, n);
        isgn.setbounds(1, n+4);
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


/*************************************************************************
Example of usage of an IterativeEstimateNorm subroutine

Input parameters:
    A   -   matrix.
            Array whose indexes range within [1..N, 1..N].

Return:
    Matrix norm estimated by the subroutine.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************/
double demoiterativeestimate1norm(const ap::real_2d_array& a, int n)
{
    double result;
    int i;
    double s;
    ap::real_1d_array x;
    ap::real_1d_array t;
    ap::real_1d_array v;
    ap::integer_1d_array iv;
    int kase;

    kase = 0;
    t.setbounds(1, n);
    iterativeestimate1norm(n, v, x, iv, result, kase);
    while(kase!=0)
    {
        if( kase==1 )
        {
            for(i = 1; i <= n; i++)
            {
                s = ap::vdotproduct(&a(i, 1), 1, &x(1), 1, ap::vlen(1,n));
                t(i) = s;
            }
        }
        else
        {
            for(i = 1; i <= n; i++)
            {
                s = ap::vdotproduct(&a(1, i), a.getstride(), &x(1), 1, ap::vlen(1,n));
                t(i) = s;
            }
        }
        ap::vmove(&x(1), 1, &t(1), 1, ap::vlen(1,n));
        iterativeestimate1norm(n, v, x, iv, result, kase);
    }
    return result;
}




