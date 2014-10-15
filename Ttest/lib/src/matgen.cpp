/*************************************************************************
Copyright (c) 2007, Sergey Bochkanov (ALGLIB project).

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
#include "matgen.h"

/*************************************************************************
Generation of a random uniformly distributed (Haar) orthogonal matrix

INPUT PARAMETERS:
    N   -   matrix size, N>=1
    
OUTPUT PARAMETERS:
    A   -   orthogonal NxN matrix, array[0..N-1,0..N-1]

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************/
void rmatrixrndorthogonal(int n, ap::real_2d_array& a)
{
    int i;
    int j;

    ap::ap_error::make_assertion(n>=1, "RMatrixRndOrthogonal: N<1!");
    a.setbounds(0, n-1, 0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            if( i==j )
            {
                a(i,j) = 1;
            }
            else
            {
                a(i,j) = 0;
            }
        }
    }
    rmatrixrndorthogonalfromtheright(a, n, n);
}


/*************************************************************************
Generation of random NxN matrix with given condition number and norm2(A)=1

INPUT PARAMETERS:
    N   -   matrix size
    C   -   condition number (in 2-norm)

OUTPUT PARAMETERS:
    A   -   random matrix with norm2(A)=1 and cond(A)=C

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************/
void rmatrixrndcond(int n, double c, ap::real_2d_array& a)
{
    int i;
    int j;
    double l1;
    double l2;

    ap::ap_error::make_assertion(n>=1&&ap::fp_greater_eq(c,1), "RMatrixRndCond: N<1 or C<1!");
    a.setbounds(0, n-1, 0, n-1);
    if( n==1 )
    {
        
        //
        // special case
        //
        a(0,0) = 2*ap::randominteger(2)-1;
        return;
    }
    l1 = 0;
    l2 = log(1/c);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            a(i,j) = 0;
        }
    }
    a(0,0) = exp(l1);
    for(i = 1; i <= n-2; i++)
    {
        a(i,i) = exp(ap::randomreal()*(l2-l1)+l1);
    }
    a(n-1,n-1) = exp(l2);
    rmatrixrndorthogonalfromtheleft(a, n, n);
    rmatrixrndorthogonalfromtheright(a, n, n);
}


/*************************************************************************
Generation of a random Haar distributed orthogonal complex matrix

INPUT PARAMETERS:
    N   -   matrix size, N>=1

OUTPUT PARAMETERS:
    A   -   orthogonal NxN matrix, array[0..N-1,0..N-1]

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************/
void cmatrixrndorthogonal(int n, ap::complex_2d_array& a)
{
    int i;
    int j;

    ap::ap_error::make_assertion(n>=1, "CMatrixRndOrthogonal: N<1!");
    a.setbounds(0, n-1, 0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            if( i==j )
            {
                a(i,j) = 1;
            }
            else
            {
                a(i,j) = 0;
            }
        }
    }
    cmatrixrndorthogonalfromtheright(a, n, n);
}


/*************************************************************************
Generation of random NxN complex matrix with given condition number C and
norm2(A)=1

INPUT PARAMETERS:
    N   -   matrix size
    C   -   condition number (in 2-norm)

OUTPUT PARAMETERS:
    A   -   random matrix with norm2(A)=1 and cond(A)=C

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************/
void cmatrixrndcond(int n, double c, ap::complex_2d_array& a)
{
    int i;
    int j;
    double l1;
    double l2;
    hqrndstate state;
    ap::complex v;

    ap::ap_error::make_assertion(n>=1&&ap::fp_greater_eq(c,1), "CMatrixRndCond: N<1 or C<1!");
    a.setbounds(0, n-1, 0, n-1);
    if( n==1 )
    {
        
        //
        // special case
        //
        hqrndrandomize(state);
        hqrndunit2(state, v.x, v.y);
        a(0,0) = v;
        return;
    }
    l1 = 0;
    l2 = log(1/c);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            a(i,j) = 0;
        }
    }
    a(0,0) = exp(l1);
    for(i = 1; i <= n-2; i++)
    {
        a(i,i) = exp(ap::randomreal()*(l2-l1)+l1);
    }
    a(n-1,n-1) = exp(l2);
    cmatrixrndorthogonalfromtheleft(a, n, n);
    cmatrixrndorthogonalfromtheright(a, n, n);
}


/*************************************************************************
Generation of random NxN symmetric matrix with given condition number  and
norm2(A)=1

INPUT PARAMETERS:
    N   -   matrix size
    C   -   condition number (in 2-norm)

OUTPUT PARAMETERS:
    A   -   random matrix with norm2(A)=1 and cond(A)=C

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************/
void smatrixrndcond(int n, double c, ap::real_2d_array& a)
{
    int i;
    int j;
    double l1;
    double l2;

    ap::ap_error::make_assertion(n>=1&&ap::fp_greater_eq(c,1), "SMatrixRndCond: N<1 or C<1!");
    a.setbounds(0, n-1, 0, n-1);
    if( n==1 )
    {
        
        //
        // special case
        //
        a(0,0) = 2*ap::randominteger(2)-1;
        return;
    }
    
    //
    // Prepare matrix
    //
    l1 = 0;
    l2 = log(1/c);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            a(i,j) = 0;
        }
    }
    a(0,0) = exp(l1);
    for(i = 1; i <= n-2; i++)
    {
        a(i,i) = (2*ap::randominteger(2)-1)*exp(ap::randomreal()*(l2-l1)+l1);
    }
    a(n-1,n-1) = exp(l2);
    
    //
    // Multiply
    //
    smatrixrndmultiply(a, n);
}


/*************************************************************************
Generation of random NxN symmetric positive definite matrix with given
condition number and norm2(A)=1

INPUT PARAMETERS:
    N   -   matrix size
    C   -   condition number (in 2-norm)

OUTPUT PARAMETERS:
    A   -   random SPD matrix with norm2(A)=1 and cond(A)=C

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************/
void spdmatrixrndcond(int n, double c, ap::real_2d_array& a)
{
    int i;
    int j;
    double l1;
    double l2;

    
    //
    // Special cases
    //
    if( n<=0||ap::fp_less(c,1) )
    {
        return;
    }
    a.setbounds(0, n-1, 0, n-1);
    if( n==1 )
    {
        a(0,0) = 1;
        return;
    }
    
    //
    // Prepare matrix
    //
    l1 = 0;
    l2 = log(1/c);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            a(i,j) = 0;
        }
    }
    a(0,0) = exp(l1);
    for(i = 1; i <= n-2; i++)
    {
        a(i,i) = exp(ap::randomreal()*(l2-l1)+l1);
    }
    a(n-1,n-1) = exp(l2);
    
    //
    // Multiply
    //
    smatrixrndmultiply(a, n);
}


/*************************************************************************
Generation of random NxN Hermitian matrix with given condition number  and
norm2(A)=1

INPUT PARAMETERS:
    N   -   matrix size
    C   -   condition number (in 2-norm)

OUTPUT PARAMETERS:
    A   -   random matrix with norm2(A)=1 and cond(A)=C

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************/
void hmatrixrndcond(int n, double c, ap::complex_2d_array& a)
{
    int i;
    int j;
    double l1;
    double l2;

    ap::ap_error::make_assertion(n>=1&&ap::fp_greater_eq(c,1), "HMatrixRndCond: N<1 or C<1!");
    a.setbounds(0, n-1, 0, n-1);
    if( n==1 )
    {
        
        //
        // special case
        //
        a(0,0) = 2*ap::randominteger(2)-1;
        return;
    }
    
    //
    // Prepare matrix
    //
    l1 = 0;
    l2 = log(1/c);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            a(i,j) = 0;
        }
    }
    a(0,0) = exp(l1);
    for(i = 1; i <= n-2; i++)
    {
        a(i,i) = (2*ap::randominteger(2)-1)*exp(ap::randomreal()*(l2-l1)+l1);
    }
    a(n-1,n-1) = exp(l2);
    
    //
    // Multiply
    //
    hmatrixrndmultiply(a, n);
    
    //
    // post-process to ensure that matrix diagonal is real
    //
    for(i = 0; i <= n-1; i++)
    {
        a(i,i).y = 0;
    }
}


/*************************************************************************
Generation of random NxN Hermitian positive definite matrix with given
condition number and norm2(A)=1

INPUT PARAMETERS:
    N   -   matrix size
    C   -   condition number (in 2-norm)

OUTPUT PARAMETERS:
    A   -   random HPD matrix with norm2(A)=1 and cond(A)=C

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************/
void hpdmatrixrndcond(int n, double c, ap::complex_2d_array& a)
{
    int i;
    int j;
    double l1;
    double l2;

    
    //
    // Special cases
    //
    if( n<=0||ap::fp_less(c,1) )
    {
        return;
    }
    a.setbounds(0, n-1, 0, n-1);
    if( n==1 )
    {
        a(0,0) = 1;
        return;
    }
    
    //
    // Prepare matrix
    //
    l1 = 0;
    l2 = log(1/c);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            a(i,j) = 0;
        }
    }
    a(0,0) = exp(l1);
    for(i = 1; i <= n-2; i++)
    {
        a(i,i) = exp(ap::randomreal()*(l2-l1)+l1);
    }
    a(n-1,n-1) = exp(l2);
    
    //
    // Multiply
    //
    hmatrixrndmultiply(a, n);
    
    //
    // post-process to ensure that matrix diagonal is real
    //
    for(i = 0; i <= n-1; i++)
    {
        a(i,i).y = 0;
    }
}


/*************************************************************************
Multiplication of MxN matrix by NxN random Haar distributed orthogonal matrix

INPUT PARAMETERS:
    A   -   matrix, array[0..M-1, 0..N-1]
    M, N-   matrix size

OUTPUT PARAMETERS:
    A   -   A*Q, where Q is random NxN orthogonal matrix

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************/
void rmatrixrndorthogonalfromtheright(ap::real_2d_array& a, int m, int n)
{
    double tau;
    double lambda;
    int s;
    int i;
    double u1;
    double u2;
    ap::real_1d_array w;
    ap::real_1d_array v;
    hqrndstate state;

    ap::ap_error::make_assertion(n>=1&&m>=1, "RMatrixRndOrthogonalFromTheRight: N<1 or M<1!");
    if( n==1 )
    {
        
        //
        // Special case
        //
        tau = 2*ap::randominteger(2)-1;
        for(i = 0; i <= m-1; i++)
        {
            a(i,0) = a(i,0)*tau;
        }
        return;
    }
    
    //
    // General case.
    // First pass.
    //
    w.setbounds(0, m-1);
    v.setbounds(1, n);
    hqrndrandomize(state);
    for(s = 2; s <= n; s++)
    {
        
        //
        // Prepare random normal v
        //
        do
        {
            i = 1;
            while(i<=s)
            {
                hqrndnormal2(state, u1, u2);
                v(i) = u1;
                if( i+1<=s )
                {
                    v(i+1) = u2;
                }
                i = i+2;
            }
            lambda = ap::vdotproduct(&v(1), 1, &v(1), 1, ap::vlen(1,s));
        }
        while(ap::fp_eq(lambda,0));
        
        //
        // Prepare and apply reflection
        //
        generatereflection(v, s, tau);
        v(1) = 1;
        applyreflectionfromtheright(a, tau, v, 0, m-1, n-s, n-1, w);
    }
    
    //
    // Second pass.
    //
    for(i = 0; i <= n-1; i++)
    {
        tau = 2*ap::randominteger(2)-1;
        ap::vmul(&a(0, i), a.getstride(), ap::vlen(0,m-1), tau);
    }
}


/*************************************************************************
Multiplication of MxN matrix by MxM random Haar distributed orthogonal matrix

INPUT PARAMETERS:
    A   -   matrix, array[0..M-1, 0..N-1]
    M, N-   matrix size

OUTPUT PARAMETERS:
    A   -   Q*A, where Q is random MxM orthogonal matrix

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************/
void rmatrixrndorthogonalfromtheleft(ap::real_2d_array& a, int m, int n)
{
    double tau;
    double lambda;
    int s;
    int i;
    int j;
    double u1;
    double u2;
    ap::real_1d_array w;
    ap::real_1d_array v;
    hqrndstate state;

    ap::ap_error::make_assertion(n>=1&&m>=1, "RMatrixRndOrthogonalFromTheRight: N<1 or M<1!");
    if( m==1 )
    {
        
        //
        // special case
        //
        tau = 2*ap::randominteger(2)-1;
        for(j = 0; j <= n-1; j++)
        {
            a(0,j) = a(0,j)*tau;
        }
        return;
    }
    
    //
    // General case.
    // First pass.
    //
    w.setbounds(0, n-1);
    v.setbounds(1, m);
    hqrndrandomize(state);
    for(s = 2; s <= m; s++)
    {
        
        //
        // Prepare random normal v
        //
        do
        {
            i = 1;
            while(i<=s)
            {
                hqrndnormal2(state, u1, u2);
                v(i) = u1;
                if( i+1<=s )
                {
                    v(i+1) = u2;
                }
                i = i+2;
            }
            lambda = ap::vdotproduct(&v(1), 1, &v(1), 1, ap::vlen(1,s));
        }
        while(ap::fp_eq(lambda,0));
        
        //
        // Prepare and apply reflection
        //
        generatereflection(v, s, tau);
        v(1) = 1;
        applyreflectionfromtheleft(a, tau, v, m-s, m-1, 0, n-1, w);
    }
    
    //
    // Second pass.
    //
    for(i = 0; i <= m-1; i++)
    {
        tau = 2*ap::randominteger(2)-1;
        ap::vmul(&a(i, 0), 1, ap::vlen(0,n-1), tau);
    }
}


/*************************************************************************
Multiplication of MxN complex matrix by NxN random Haar distributed
complex orthogonal matrix

INPUT PARAMETERS:
    A   -   matrix, array[0..M-1, 0..N-1]
    M, N-   matrix size

OUTPUT PARAMETERS:
    A   -   A*Q, where Q is random NxN orthogonal matrix

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************/
void cmatrixrndorthogonalfromtheright(ap::complex_2d_array& a, int m, int n)
{
    ap::complex lambda;
    ap::complex tau;
    int s;
    int i;
    ap::complex_1d_array w;
    ap::complex_1d_array v;
    hqrndstate state;

    ap::ap_error::make_assertion(n>=1&&m>=1, "CMatrixRndOrthogonalFromTheRight: N<1 or M<1!");
    if( n==1 )
    {
        
        //
        // Special case
        //
        hqrndrandomize(state);
        hqrndunit2(state, tau.x, tau.y);
        for(i = 0; i <= m-1; i++)
        {
            a(i,0) = a(i,0)*tau;
        }
        return;
    }
    
    //
    // General case.
    // First pass.
    //
    w.setbounds(0, m-1);
    v.setbounds(1, n);
    hqrndrandomize(state);
    for(s = 2; s <= n; s++)
    {
        
        //
        // Prepare random normal v
        //
        do
        {
            for(i = 1; i <= s; i++)
            {
                hqrndnormal2(state, tau.x, tau.y);
                v(i) = tau;
            }
            lambda = ap::vdotproduct(&v(1), 1, "N", &v(1), 1, "Conj", ap::vlen(1,s));
        }
        while(lambda==0);
        
        //
        // Prepare and apply reflection
        //
        complexgeneratereflection(v, s, tau);
        v(1) = 1;
        complexapplyreflectionfromtheright(a, tau, v, 0, m-1, n-s, n-1, w);
    }
    
    //
    // Second pass.
    //
    for(i = 0; i <= n-1; i++)
    {
        hqrndunit2(state, tau.x, tau.y);
        ap::vmul(&a(0, i), a.getstride(), ap::vlen(0,m-1), tau);
    }
}


/*************************************************************************
Multiplication of MxN complex matrix by MxM random Haar distributed
complex orthogonal matrix

INPUT PARAMETERS:
    A   -   matrix, array[0..M-1, 0..N-1]
    M, N-   matrix size

OUTPUT PARAMETERS:
    A   -   Q*A, where Q is random MxM orthogonal matrix

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************/
void cmatrixrndorthogonalfromtheleft(ap::complex_2d_array& a, int m, int n)
{
    ap::complex tau;
    ap::complex lambda;
    int s;
    int i;
    int j;
    ap::complex_1d_array w;
    ap::complex_1d_array v;
    hqrndstate state;

    ap::ap_error::make_assertion(n>=1&&m>=1, "CMatrixRndOrthogonalFromTheRight: N<1 or M<1!");
    if( m==1 )
    {
        
        //
        // special case
        //
        hqrndrandomize(state);
        hqrndunit2(state, tau.x, tau.y);
        for(j = 0; j <= n-1; j++)
        {
            a(0,j) = a(0,j)*tau;
        }
        return;
    }
    
    //
    // General case.
    // First pass.
    //
    w.setbounds(0, n-1);
    v.setbounds(1, m);
    hqrndrandomize(state);
    for(s = 2; s <= m; s++)
    {
        
        //
        // Prepare random normal v
        //
        do
        {
            for(i = 1; i <= s; i++)
            {
                hqrndnormal2(state, tau.x, tau.y);
                v(i) = tau;
            }
            lambda = ap::vdotproduct(&v(1), 1, "N", &v(1), 1, "Conj", ap::vlen(1,s));
        }
        while(lambda==0);
        
        //
        // Prepare and apply reflection
        //
        complexgeneratereflection(v, s, tau);
        v(1) = 1;
        complexapplyreflectionfromtheleft(a, tau, v, m-s, m-1, 0, n-1, w);
    }
    
    //
    // Second pass.
    //
    for(i = 0; i <= m-1; i++)
    {
        hqrndunit2(state, tau.x, tau.y);
        ap::vmul(&a(i, 0), 1, ap::vlen(0,n-1), tau);
    }
}


/*************************************************************************
Symmetric multiplication of NxN matrix by random Haar distributed
orthogonal  matrix

INPUT PARAMETERS:
    A   -   matrix, array[0..N-1, 0..N-1]
    N   -   matrix size

OUTPUT PARAMETERS:
    A   -   Q'*A*Q, where Q is random NxN orthogonal matrix

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************/
void smatrixrndmultiply(ap::real_2d_array& a, int n)
{
    double tau;
    double lambda;
    int s;
    int i;
    double u1;
    double u2;
    ap::real_1d_array w;
    ap::real_1d_array v;
    hqrndstate state;

    
    //
    // General case.
    //
    w.setbounds(0, n-1);
    v.setbounds(1, n);
    hqrndrandomize(state);
    for(s = 2; s <= n; s++)
    {
        
        //
        // Prepare random normal v
        //
        do
        {
            i = 1;
            while(i<=s)
            {
                hqrndnormal2(state, u1, u2);
                v(i) = u1;
                if( i+1<=s )
                {
                    v(i+1) = u2;
                }
                i = i+2;
            }
            lambda = ap::vdotproduct(&v(1), 1, &v(1), 1, ap::vlen(1,s));
        }
        while(ap::fp_eq(lambda,0));
        
        //
        // Prepare and apply reflection
        //
        generatereflection(v, s, tau);
        v(1) = 1;
        applyreflectionfromtheright(a, tau, v, 0, n-1, n-s, n-1, w);
        applyreflectionfromtheleft(a, tau, v, n-s, n-1, 0, n-1, w);
    }
    
    //
    // Second pass.
    //
    for(i = 0; i <= n-1; i++)
    {
        tau = 2*ap::randominteger(2)-1;
        ap::vmul(&a(0, i), a.getstride(), ap::vlen(0,n-1), tau);
        ap::vmul(&a(i, 0), 1, ap::vlen(0,n-1), tau);
    }
}


/*************************************************************************
Hermitian multiplication of NxN matrix by random Haar distributed
complex orthogonal matrix

INPUT PARAMETERS:
    A   -   matrix, array[0..N-1, 0..N-1]
    N   -   matrix size

OUTPUT PARAMETERS:
    A   -   Q^H*A*Q, where Q is random NxN orthogonal matrix

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************/
void hmatrixrndmultiply(ap::complex_2d_array& a, int n)
{
    ap::complex tau;
    ap::complex lambda;
    int s;
    int i;
    ap::complex_1d_array w;
    ap::complex_1d_array v;
    hqrndstate state;

    
    //
    // General case.
    //
    w.setbounds(0, n-1);
    v.setbounds(1, n);
    hqrndrandomize(state);
    for(s = 2; s <= n; s++)
    {
        
        //
        // Prepare random normal v
        //
        do
        {
            for(i = 1; i <= s; i++)
            {
                hqrndnormal2(state, tau.x, tau.y);
                v(i) = tau;
            }
            lambda = ap::vdotproduct(&v(1), 1, "N", &v(1), 1, "Conj", ap::vlen(1,s));
        }
        while(lambda==0);
        
        //
        // Prepare and apply reflection
        //
        complexgeneratereflection(v, s, tau);
        v(1) = 1;
        complexapplyreflectionfromtheright(a, tau, v, 0, n-1, n-s, n-1, w);
        complexapplyreflectionfromtheleft(a, ap::conj(tau), v, n-s, n-1, 0, n-1, w);
    }
    
    //
    // Second pass.
    //
    for(i = 0; i <= n-1; i++)
    {
        hqrndunit2(state, tau.x, tau.y);
        ap::vmul(&a(0, i), a.getstride(), ap::vlen(0,n-1), tau);
        tau = ap::conj(tau);
        ap::vmul(&a(i, 0), 1, ap::vlen(0,n-1), tau);
    }
}




