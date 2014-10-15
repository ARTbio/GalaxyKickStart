/*************************************************************************
This file is a part of ALGLIB project.

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
#include "safesolve.h"

static bool cbasicsolveandupdate(ap::complex alpha,
     ap::complex beta,
     double lnmax,
     double bnorm,
     double maxgrowth,
     double& xnorm,
     ap::complex& x);

/*************************************************************************
Real implementation of CMatrixScaledTRSafeSolve

  -- ALGLIB routine --
     21.01.2010
     Bochkanov Sergey
*************************************************************************/
bool rmatrixscaledtrsafesolve(const ap::real_2d_array& a,
     double sa,
     int n,
     ap::real_1d_array& x,
     bool isupper,
     int trans,
     bool isunit,
     double maxgrowth)
{
    bool result;
    double lnmax;
    double nrmb;
    double nrmx;
    int i;
    ap::complex alpha;
    ap::complex beta;
    double vr;
    ap::complex cx;
    ap::real_1d_array tmp;

    ap::ap_error::make_assertion(n>0, "RMatrixTRSafeSolve: incorrect N!");
    ap::ap_error::make_assertion(trans==0||trans==1, "RMatrixTRSafeSolve: incorrect Trans!");
    result = true;
    lnmax = log(ap::maxrealnumber);
    
    //
    // Quick return if possible
    //
    if( n<=0 )
    {
        return result;
    }
    
    //
    // Load norms: right part and X
    //
    nrmb = 0;
    for(i = 0; i <= n-1; i++)
    {
        nrmb = ap::maxreal(nrmb, fabs(x(i)));
    }
    nrmx = 0;
    
    //
    // Solve
    //
    tmp.setlength(n);
    result = true;
    if( isupper&&trans==0 )
    {
        
        //
        // U*x = b
        //
        for(i = n-1; i >= 0; i--)
        {
            
            //
            // Task is reduced to alpha*x[i] = beta
            //
            if( isunit )
            {
                alpha = sa;
            }
            else
            {
                alpha = a(i,i)*sa;
            }
            if( i<n-1 )
            {
                ap::vmove(&tmp(i+1), 1, &a(i, i+1), 1, ap::vlen(i+1,n-1), sa);
                vr = ap::vdotproduct(&tmp(i+1), 1, &x(i+1), 1, ap::vlen(i+1,n-1));
                beta = x(i)-vr;
            }
            else
            {
                beta = x(i);
            }
            
            //
            // solve alpha*x[i] = beta
            //
            result = cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, nrmx, cx);
            if( !result )
            {
                return result;
            }
            x(i) = cx.x;
        }
        return result;
    }
    if( !isupper&&trans==0 )
    {
        
        //
        // L*x = b
        //
        for(i = 0; i <= n-1; i++)
        {
            
            //
            // Task is reduced to alpha*x[i] = beta
            //
            if( isunit )
            {
                alpha = sa;
            }
            else
            {
                alpha = a(i,i)*sa;
            }
            if( i>0 )
            {
                ap::vmove(&tmp(0), 1, &a(i, 0), 1, ap::vlen(0,i-1), sa);
                vr = ap::vdotproduct(&tmp(0), 1, &x(0), 1, ap::vlen(0,i-1));
                beta = x(i)-vr;
            }
            else
            {
                beta = x(i);
            }
            
            //
            // solve alpha*x[i] = beta
            //
            result = cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, nrmx, cx);
            if( !result )
            {
                return result;
            }
            x(i) = cx.x;
        }
        return result;
    }
    if( isupper&&trans==1 )
    {
        
        //
        // U^T*x = b
        //
        for(i = 0; i <= n-1; i++)
        {
            
            //
            // Task is reduced to alpha*x[i] = beta
            //
            if( isunit )
            {
                alpha = sa;
            }
            else
            {
                alpha = a(i,i)*sa;
            }
            beta = x(i);
            
            //
            // solve alpha*x[i] = beta
            //
            result = cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, nrmx, cx);
            if( !result )
            {
                return result;
            }
            x(i) = cx.x;
            
            //
            // update the rest of right part
            //
            if( i<n-1 )
            {
                vr = cx.x;
                ap::vmove(&tmp(i+1), 1, &a(i, i+1), 1, ap::vlen(i+1,n-1), sa);
                ap::vsub(&x(i+1), 1, &tmp(i+1), 1, ap::vlen(i+1,n-1), vr);
            }
        }
        return result;
    }
    if( !isupper&&trans==1 )
    {
        
        //
        // L^T*x = b
        //
        for(i = n-1; i >= 0; i--)
        {
            
            //
            // Task is reduced to alpha*x[i] = beta
            //
            if( isunit )
            {
                alpha = sa;
            }
            else
            {
                alpha = a(i,i)*sa;
            }
            beta = x(i);
            
            //
            // solve alpha*x[i] = beta
            //
            result = cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, nrmx, cx);
            if( !result )
            {
                return result;
            }
            x(i) = cx.x;
            
            //
            // update the rest of right part
            //
            if( i>0 )
            {
                vr = cx.x;
                ap::vmove(&tmp(0), 1, &a(i, 0), 1, ap::vlen(0,i-1), sa);
                ap::vsub(&x(0), 1, &tmp(0), 1, ap::vlen(0,i-1), vr);
            }
        }
        return result;
    }
    result = false;
    return result;
}


/*************************************************************************
Internal subroutine for safe solution of

    SA*op(A)=b
    
where  A  is  NxN  upper/lower  triangular/unitriangular  matrix, op(A) is
either identity transform, transposition or Hermitian transposition, SA is
a scaling factor such that max(|SA*A[i,j]|) is close to 1.0 in magnutude.

This subroutine  limits  relative  growth  of  solution  (in inf-norm)  by
MaxGrowth,  returning  False  if  growth  exceeds MaxGrowth. Degenerate or
near-degenerate matrices are handled correctly (False is returned) as long
as MaxGrowth is significantly less than MaxRealNumber/norm(b).

  -- ALGLIB routine --
     21.01.2010
     Bochkanov Sergey
*************************************************************************/
bool cmatrixscaledtrsafesolve(const ap::complex_2d_array& a,
     double sa,
     int n,
     ap::complex_1d_array& x,
     bool isupper,
     int trans,
     bool isunit,
     double maxgrowth)
{
    bool result;
    double lnmax;
    double nrmb;
    double nrmx;
    int i;
    ap::complex alpha;
    ap::complex beta;
    ap::complex vc;
    ap::complex_1d_array tmp;

    ap::ap_error::make_assertion(n>0, "CMatrixTRSafeSolve: incorrect N!");
    ap::ap_error::make_assertion(trans==0||trans==1||trans==2, "CMatrixTRSafeSolve: incorrect Trans!");
    result = true;
    lnmax = log(ap::maxrealnumber);
    
    //
    // Quick return if possible
    //
    if( n<=0 )
    {
        return result;
    }
    
    //
    // Load norms: right part and X
    //
    nrmb = 0;
    for(i = 0; i <= n-1; i++)
    {
        nrmb = ap::maxreal(nrmb, ap::abscomplex(x(i)));
    }
    nrmx = 0;
    
    //
    // Solve
    //
    tmp.setlength(n);
    result = true;
    if( isupper&&trans==0 )
    {
        
        //
        // U*x = b
        //
        for(i = n-1; i >= 0; i--)
        {
            
            //
            // Task is reduced to alpha*x[i] = beta
            //
            if( isunit )
            {
                alpha = sa;
            }
            else
            {
                alpha = a(i,i)*sa;
            }
            if( i<n-1 )
            {
                ap::vmove(&tmp(i+1), 1, &a(i, i+1), 1, "N", ap::vlen(i+1,n-1), sa);
                vc = ap::vdotproduct(&tmp(i+1), 1, "N", &x(i+1), 1, "N", ap::vlen(i+1,n-1));
                beta = x(i)-vc;
            }
            else
            {
                beta = x(i);
            }
            
            //
            // solve alpha*x[i] = beta
            //
            result = cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, nrmx, vc);
            if( !result )
            {
                return result;
            }
            x(i) = vc;
        }
        return result;
    }
    if( !isupper&&trans==0 )
    {
        
        //
        // L*x = b
        //
        for(i = 0; i <= n-1; i++)
        {
            
            //
            // Task is reduced to alpha*x[i] = beta
            //
            if( isunit )
            {
                alpha = sa;
            }
            else
            {
                alpha = a(i,i)*sa;
            }
            if( i>0 )
            {
                ap::vmove(&tmp(0), 1, &a(i, 0), 1, "N", ap::vlen(0,i-1), sa);
                vc = ap::vdotproduct(&tmp(0), 1, "N", &x(0), 1, "N", ap::vlen(0,i-1));
                beta = x(i)-vc;
            }
            else
            {
                beta = x(i);
            }
            
            //
            // solve alpha*x[i] = beta
            //
            result = cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, nrmx, vc);
            if( !result )
            {
                return result;
            }
            x(i) = vc;
        }
        return result;
    }
    if( isupper&&trans==1 )
    {
        
        //
        // U^T*x = b
        //
        for(i = 0; i <= n-1; i++)
        {
            
            //
            // Task is reduced to alpha*x[i] = beta
            //
            if( isunit )
            {
                alpha = sa;
            }
            else
            {
                alpha = a(i,i)*sa;
            }
            beta = x(i);
            
            //
            // solve alpha*x[i] = beta
            //
            result = cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, nrmx, vc);
            if( !result )
            {
                return result;
            }
            x(i) = vc;
            
            //
            // update the rest of right part
            //
            if( i<n-1 )
            {
                ap::vmove(&tmp(i+1), 1, &a(i, i+1), 1, "N", ap::vlen(i+1,n-1), sa);
                ap::vsub(&x(i+1), 1, &tmp(i+1), 1, "N", ap::vlen(i+1,n-1), vc);
            }
        }
        return result;
    }
    if( !isupper&&trans==1 )
    {
        
        //
        // L^T*x = b
        //
        for(i = n-1; i >= 0; i--)
        {
            
            //
            // Task is reduced to alpha*x[i] = beta
            //
            if( isunit )
            {
                alpha = sa;
            }
            else
            {
                alpha = a(i,i)*sa;
            }
            beta = x(i);
            
            //
            // solve alpha*x[i] = beta
            //
            result = cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, nrmx, vc);
            if( !result )
            {
                return result;
            }
            x(i) = vc;
            
            //
            // update the rest of right part
            //
            if( i>0 )
            {
                ap::vmove(&tmp(0), 1, &a(i, 0), 1, "N", ap::vlen(0,i-1), sa);
                ap::vsub(&x(0), 1, &tmp(0), 1, "N", ap::vlen(0,i-1), vc);
            }
        }
        return result;
    }
    if( isupper&&trans==2 )
    {
        
        //
        // U^H*x = b
        //
        for(i = 0; i <= n-1; i++)
        {
            
            //
            // Task is reduced to alpha*x[i] = beta
            //
            if( isunit )
            {
                alpha = sa;
            }
            else
            {
                alpha = ap::conj(a(i,i))*sa;
            }
            beta = x(i);
            
            //
            // solve alpha*x[i] = beta
            //
            result = cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, nrmx, vc);
            if( !result )
            {
                return result;
            }
            x(i) = vc;
            
            //
            // update the rest of right part
            //
            if( i<n-1 )
            {
                ap::vmove(&tmp(i+1), 1, &a(i, i+1), 1, "Conj", ap::vlen(i+1,n-1), sa);
                ap::vsub(&x(i+1), 1, &tmp(i+1), 1, "N", ap::vlen(i+1,n-1), vc);
            }
        }
        return result;
    }
    if( !isupper&&trans==2 )
    {
        
        //
        // L^T*x = b
        //
        for(i = n-1; i >= 0; i--)
        {
            
            //
            // Task is reduced to alpha*x[i] = beta
            //
            if( isunit )
            {
                alpha = sa;
            }
            else
            {
                alpha = ap::conj(a(i,i))*sa;
            }
            beta = x(i);
            
            //
            // solve alpha*x[i] = beta
            //
            result = cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, nrmx, vc);
            if( !result )
            {
                return result;
            }
            x(i) = vc;
            
            //
            // update the rest of right part
            //
            if( i>0 )
            {
                ap::vmove(&tmp(0), 1, &a(i, 0), 1, "Conj", ap::vlen(0,i-1), sa);
                ap::vsub(&x(0), 1, &tmp(0), 1, "N", ap::vlen(0,i-1), vc);
            }
        }
        return result;
    }
    result = false;
    return result;
}


/*************************************************************************
complex basic solver-updater for reduced linear system

    alpha*x[i] = beta

solves this equation and updates it in overlfow-safe manner (keeping track
of relative growth of solution).

Parameters:
    Alpha   -   alpha
    Beta    -   beta
    LnMax   -   precomputed Ln(MaxRealNumber)
    BNorm   -   inf-norm of b (right part of original system)
    MaxGrowth-  maximum growth of norm(x) relative to norm(b)
    XNorm   -   inf-norm of other components of X (which are already processed)
                it is updated by CBasicSolveAndUpdate.
    X       -   solution

  -- ALGLIB routine --
     26.01.2009
     Bochkanov Sergey
*************************************************************************/
static bool cbasicsolveandupdate(ap::complex alpha,
     ap::complex beta,
     double lnmax,
     double bnorm,
     double maxgrowth,
     double& xnorm,
     ap::complex& x)
{
    bool result;
    double v;

    result = false;
    if( alpha==0 )
    {
        return result;
    }
    if( beta!=0 )
    {
        
        //
        // alpha*x[i]=beta
        //
        v = log(ap::abscomplex(beta))-log(ap::abscomplex(alpha));
        if( ap::fp_greater(v,lnmax) )
        {
            return result;
        }
        x = beta/alpha;
    }
    else
    {
        
        //
        // alpha*x[i]=0
        //
        x = 0;
    }
    
    //
    // update NrmX, test growth limit
    //
    xnorm = ap::maxreal(xnorm, ap::abscomplex(x));
    if( ap::fp_greater(xnorm,maxgrowth*bnorm) )
    {
        return result;
    }
    result = true;
    return result;
}




