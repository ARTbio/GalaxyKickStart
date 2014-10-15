/*************************************************************************
Copyright (c) 2009, Sergey Bochkanov (ALGLIB project).

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
#include "conv.h"

/*************************************************************************
1-dimensional complex convolution.

For given A/B returns conv(A,B) (non-circular). Subroutine can automatically
choose between three implementations: straightforward O(M*N)  formula  for
very small N (or M), overlap-add algorithm for  cases  where  max(M,N)  is
significantly larger than min(M,N), but O(M*N) algorithm is too slow,  and
general FFT-based formula for cases where two previois algorithms are  too
slow.

Algorithm has max(M,N)*log(max(M,N)) complexity for any M/N.

INPUT PARAMETERS
    A   -   array[0..M-1] - complex function to be transformed
    M   -   problem size
    B   -   array[0..N-1] - complex function to be transformed
    N   -   problem size

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..N+M-2].

NOTE:
    It is assumed that A is zero at T<0, B is zero too.  If  one  or  both
functions have non-zero values at negative T's, you  can  still  use  this
subroutine - just shift its result correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convc1d(const ap::complex_1d_array& a,
     int m,
     const ap::complex_1d_array& b,
     int n,
     ap::complex_1d_array& r)
{

    ap::ap_error::make_assertion(n>0&&m>0, "ConvC1D: incorrect N or M!");
    
    //
    // normalize task: make M>=N,
    // so A will be longer that B.
    //
    if( m<n )
    {
        convc1d(b, n, a, m, r);
        return;
    }
    convc1dx(a, m, b, n, false, -1, 0, r);
}


/*************************************************************************
1-dimensional complex non-circular deconvolution (inverse of ConvC1D()).

Algorithm has M*log(M)) complexity for any M (composite or prime).

INPUT PARAMETERS
    A   -   array[0..M-1] - convolved signal, A = conv(R, B)
    M   -   convolved signal length
    B   -   array[0..N-1] - response
    N   -   response length, N<=M

OUTPUT PARAMETERS
    R   -   deconvolved signal. array[0..M-N].

NOTE:
    deconvolution is unstable process and may result in division  by  zero
(if your response function is degenerate, i.e. has zero Fourier coefficient).

NOTE:
    It is assumed that A is zero at T<0, B is zero too.  If  one  or  both
functions have non-zero values at negative T's, you  can  still  use  this
subroutine - just shift its result correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convc1dinv(const ap::complex_1d_array& a,
     int m,
     const ap::complex_1d_array& b,
     int n,
     ap::complex_1d_array& r)
{
    int i;
    int p;
    ap::real_1d_array buf;
    ap::real_1d_array buf2;
    ftplan plan;
    ap::complex c1;
    ap::complex c2;
    ap::complex c3;
    double t;

    ap::ap_error::make_assertion(n>0&&m>0&&n<=m, "ConvC1DInv: incorrect N or M!");
    p = ftbasefindsmooth(m);
    ftbasegeneratecomplexfftplan(p, plan);
    buf.setlength(2*p);
    for(i = 0; i <= m-1; i++)
    {
        buf(2*i+0) = a(i).x;
        buf(2*i+1) = a(i).y;
    }
    for(i = m; i <= p-1; i++)
    {
        buf(2*i+0) = 0;
        buf(2*i+1) = 0;
    }
    buf2.setlength(2*p);
    for(i = 0; i <= n-1; i++)
    {
        buf2(2*i+0) = b(i).x;
        buf2(2*i+1) = b(i).y;
    }
    for(i = n; i <= p-1; i++)
    {
        buf2(2*i+0) = 0;
        buf2(2*i+1) = 0;
    }
    ftbaseexecuteplan(buf, 0, p, plan);
    ftbaseexecuteplan(buf2, 0, p, plan);
    for(i = 0; i <= p-1; i++)
    {
        c1.x = buf(2*i+0);
        c1.y = buf(2*i+1);
        c2.x = buf2(2*i+0);
        c2.y = buf2(2*i+1);
        c3 = c1/c2;
        buf(2*i+0) = c3.x;
        buf(2*i+1) = -c3.y;
    }
    ftbaseexecuteplan(buf, 0, p, plan);
    t = double(1)/double(p);
    r.setlength(m-n+1);
    for(i = 0; i <= m-n; i++)
    {
        r(i).x = +t*buf(2*i+0);
        r(i).y = -t*buf(2*i+1);
    }
}


/*************************************************************************
1-dimensional circular complex convolution.

For given S/R returns conv(S,R) (circular). Algorithm has linearithmic
complexity for any M/N.

IMPORTANT:  normal convolution is commutative,  i.e.   it  is symmetric  -
conv(A,B)=conv(B,A).  Cyclic convolution IS NOT.  One function - S - is  a
signal,  periodic function, and another - R - is a response,  non-periodic
function with limited length.

INPUT PARAMETERS
    S   -   array[0..M-1] - complex periodic signal
    M   -   problem size
    B   -   array[0..N-1] - complex non-periodic response
    N   -   problem size

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..M-1].

NOTE:
    It is assumed that B is zero at T<0. If  it  has  non-zero  values  at
negative T's, you can still use this subroutine - just  shift  its  result
correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convc1dcircular(const ap::complex_1d_array& s,
     int m,
     const ap::complex_1d_array& r,
     int n,
     ap::complex_1d_array& c)
{
    ap::complex_1d_array buf;
    int i1;
    int i2;
    int j2;

    ap::ap_error::make_assertion(n>0&&m>0, "ConvC1DCircular: incorrect N or M!");
    
    //
    // normalize task: make M>=N,
    // so A will be longer (at least - not shorter) that B.
    //
    if( m<n )
    {
        buf.setlength(m);
        for(i1 = 0; i1 <= m-1; i1++)
        {
            buf(i1) = 0;
        }
        i1 = 0;
        while(i1<n)
        {
            i2 = ap::minint(i1+m-1, n-1);
            j2 = i2-i1;
            ap::vadd(&buf(0), 1, &r(i1), 1, "N", ap::vlen(0,j2));
            i1 = i1+m;
        }
        convc1dcircular(s, m, buf, m, c);
        return;
    }
    convc1dx(s, m, r, n, true, -1, 0, c);
}


/*************************************************************************
1-dimensional circular complex deconvolution (inverse of ConvC1DCircular()).

Algorithm has M*log(M)) complexity for any M (composite or prime).

INPUT PARAMETERS
    A   -   array[0..M-1] - convolved periodic signal, A = conv(R, B)
    M   -   convolved signal length
    B   -   array[0..N-1] - non-periodic response
    N   -   response length

OUTPUT PARAMETERS
    R   -   deconvolved signal. array[0..M-1].

NOTE:
    deconvolution is unstable process and may result in division  by  zero
(if your response function is degenerate, i.e. has zero Fourier coefficient).

NOTE:
    It is assumed that B is zero at T<0. If  it  has  non-zero  values  at
negative T's, you can still use this subroutine - just  shift  its  result
correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convc1dcircularinv(const ap::complex_1d_array& a,
     int m,
     const ap::complex_1d_array& b,
     int n,
     ap::complex_1d_array& r)
{
    int i;
    int i1;
    int i2;
    int j2;
    ap::real_1d_array buf;
    ap::real_1d_array buf2;
    ap::complex_1d_array cbuf;
    ftplan plan;
    ap::complex c1;
    ap::complex c2;
    ap::complex c3;
    double t;

    ap::ap_error::make_assertion(n>0&&m>0, "ConvC1DCircularInv: incorrect N or M!");
    
    //
    // normalize task: make M>=N,
    // so A will be longer (at least - not shorter) that B.
    //
    if( m<n )
    {
        cbuf.setlength(m);
        for(i = 0; i <= m-1; i++)
        {
            cbuf(i) = 0;
        }
        i1 = 0;
        while(i1<n)
        {
            i2 = ap::minint(i1+m-1, n-1);
            j2 = i2-i1;
            ap::vadd(&cbuf(0), 1, &b(i1), 1, "N", ap::vlen(0,j2));
            i1 = i1+m;
        }
        convc1dcircularinv(a, m, cbuf, m, r);
        return;
    }
    
    //
    // Task is normalized
    //
    ftbasegeneratecomplexfftplan(m, plan);
    buf.setlength(2*m);
    for(i = 0; i <= m-1; i++)
    {
        buf(2*i+0) = a(i).x;
        buf(2*i+1) = a(i).y;
    }
    buf2.setlength(2*m);
    for(i = 0; i <= n-1; i++)
    {
        buf2(2*i+0) = b(i).x;
        buf2(2*i+1) = b(i).y;
    }
    for(i = n; i <= m-1; i++)
    {
        buf2(2*i+0) = 0;
        buf2(2*i+1) = 0;
    }
    ftbaseexecuteplan(buf, 0, m, plan);
    ftbaseexecuteplan(buf2, 0, m, plan);
    for(i = 0; i <= m-1; i++)
    {
        c1.x = buf(2*i+0);
        c1.y = buf(2*i+1);
        c2.x = buf2(2*i+0);
        c2.y = buf2(2*i+1);
        c3 = c1/c2;
        buf(2*i+0) = c3.x;
        buf(2*i+1) = -c3.y;
    }
    ftbaseexecuteplan(buf, 0, m, plan);
    t = double(1)/double(m);
    r.setlength(m);
    for(i = 0; i <= m-1; i++)
    {
        r(i).x = +t*buf(2*i+0);
        r(i).y = -t*buf(2*i+1);
    }
}


/*************************************************************************
1-dimensional real convolution.

Analogous to ConvC1D(), see ConvC1D() comments for more details.

INPUT PARAMETERS
    A   -   array[0..M-1] - real function to be transformed
    M   -   problem size
    B   -   array[0..N-1] - real function to be transformed
    N   -   problem size

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..N+M-2].

NOTE:
    It is assumed that A is zero at T<0, B is zero too.  If  one  or  both
functions have non-zero values at negative T's, you  can  still  use  this
subroutine - just shift its result correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convr1d(const ap::real_1d_array& a,
     int m,
     const ap::real_1d_array& b,
     int n,
     ap::real_1d_array& r)
{

    ap::ap_error::make_assertion(n>0&&m>0, "ConvR1D: incorrect N or M!");
    
    //
    // normalize task: make M>=N,
    // so A will be longer that B.
    //
    if( m<n )
    {
        convr1d(b, n, a, m, r);
        return;
    }
    convr1dx(a, m, b, n, false, -1, 0, r);
}


/*************************************************************************
1-dimensional real deconvolution (inverse of ConvC1D()).

Algorithm has M*log(M)) complexity for any M (composite or prime).

INPUT PARAMETERS
    A   -   array[0..M-1] - convolved signal, A = conv(R, B)
    M   -   convolved signal length
    B   -   array[0..N-1] - response
    N   -   response length, N<=M

OUTPUT PARAMETERS
    R   -   deconvolved signal. array[0..M-N].

NOTE:
    deconvolution is unstable process and may result in division  by  zero
(if your response function is degenerate, i.e. has zero Fourier coefficient).

NOTE:
    It is assumed that A is zero at T<0, B is zero too.  If  one  or  both
functions have non-zero values at negative T's, you  can  still  use  this
subroutine - just shift its result correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convr1dinv(const ap::real_1d_array& a,
     int m,
     const ap::real_1d_array& b,
     int n,
     ap::real_1d_array& r)
{
    int i;
    int p;
    ap::real_1d_array buf;
    ap::real_1d_array buf2;
    ap::real_1d_array buf3;
    ftplan plan;
    ap::complex c1;
    ap::complex c2;
    ap::complex c3;

    ap::ap_error::make_assertion(n>0&&m>0&&n<=m, "ConvR1DInv: incorrect N or M!");
    p = ftbasefindsmootheven(m);
    buf.setlength(p);
    ap::vmove(&buf(0), 1, &a(0), 1, ap::vlen(0,m-1));
    for(i = m; i <= p-1; i++)
    {
        buf(i) = 0;
    }
    buf2.setlength(p);
    ap::vmove(&buf2(0), 1, &b(0), 1, ap::vlen(0,n-1));
    for(i = n; i <= p-1; i++)
    {
        buf2(i) = 0;
    }
    buf3.setlength(p);
    ftbasegeneratecomplexfftplan(p/2, plan);
    fftr1dinternaleven(buf, p, buf3, plan);
    fftr1dinternaleven(buf2, p, buf3, plan);
    buf(0) = buf(0)/buf2(0);
    buf(1) = buf(1)/buf2(1);
    for(i = 1; i <= p/2-1; i++)
    {
        c1.x = buf(2*i+0);
        c1.y = buf(2*i+1);
        c2.x = buf2(2*i+0);
        c2.y = buf2(2*i+1);
        c3 = c1/c2;
        buf(2*i+0) = c3.x;
        buf(2*i+1) = c3.y;
    }
    fftr1dinvinternaleven(buf, p, buf3, plan);
    r.setlength(m-n+1);
    ap::vmove(&r(0), 1, &buf(0), 1, ap::vlen(0,m-n));
}


/*************************************************************************
1-dimensional circular real convolution.

Analogous to ConvC1DCircular(), see ConvC1DCircular() comments for more details.

INPUT PARAMETERS
    S   -   array[0..M-1] - real signal
    M   -   problem size
    B   -   array[0..N-1] - real response
    N   -   problem size

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..M-1].

NOTE:
    It is assumed that B is zero at T<0. If  it  has  non-zero  values  at
negative T's, you can still use this subroutine - just  shift  its  result
correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convr1dcircular(const ap::real_1d_array& s,
     int m,
     const ap::real_1d_array& r,
     int n,
     ap::real_1d_array& c)
{
    ap::real_1d_array buf;
    int i1;
    int i2;
    int j2;

    ap::ap_error::make_assertion(n>0&&m>0, "ConvC1DCircular: incorrect N or M!");
    
    //
    // normalize task: make M>=N,
    // so A will be longer (at least - not shorter) that B.
    //
    if( m<n )
    {
        buf.setlength(m);
        for(i1 = 0; i1 <= m-1; i1++)
        {
            buf(i1) = 0;
        }
        i1 = 0;
        while(i1<n)
        {
            i2 = ap::minint(i1+m-1, n-1);
            j2 = i2-i1;
            ap::vadd(&buf(0), 1, &r(i1), 1, ap::vlen(0,j2));
            i1 = i1+m;
        }
        convr1dcircular(s, m, buf, m, c);
        return;
    }
    
    //
    // reduce to usual convolution
    //
    convr1dx(s, m, r, n, true, -1, 0, c);
}


/*************************************************************************
1-dimensional complex deconvolution (inverse of ConvC1D()).

Algorithm has M*log(M)) complexity for any M (composite or prime).

INPUT PARAMETERS
    A   -   array[0..M-1] - convolved signal, A = conv(R, B)
    M   -   convolved signal length
    B   -   array[0..N-1] - response
    N   -   response length

OUTPUT PARAMETERS
    R   -   deconvolved signal. array[0..M-N].

NOTE:
    deconvolution is unstable process and may result in division  by  zero
(if your response function is degenerate, i.e. has zero Fourier coefficient).

NOTE:
    It is assumed that B is zero at T<0. If  it  has  non-zero  values  at
negative T's, you can still use this subroutine - just  shift  its  result
correspondingly.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convr1dcircularinv(const ap::real_1d_array& a,
     int m,
     const ap::real_1d_array& b,
     int n,
     ap::real_1d_array& r)
{
    int i;
    int i1;
    int i2;
    int j2;
    ap::real_1d_array buf;
    ap::real_1d_array buf2;
    ap::real_1d_array buf3;
    ap::complex_1d_array cbuf;
    ap::complex_1d_array cbuf2;
    ftplan plan;
    ap::complex c1;
    ap::complex c2;
    ap::complex c3;

    ap::ap_error::make_assertion(n>0&&m>0, "ConvR1DCircularInv: incorrect N or M!");
    
    //
    // normalize task: make M>=N,
    // so A will be longer (at least - not shorter) that B.
    //
    if( m<n )
    {
        buf.setlength(m);
        for(i = 0; i <= m-1; i++)
        {
            buf(i) = 0;
        }
        i1 = 0;
        while(i1<n)
        {
            i2 = ap::minint(i1+m-1, n-1);
            j2 = i2-i1;
            ap::vadd(&buf(0), 1, &b(i1), 1, ap::vlen(0,j2));
            i1 = i1+m;
        }
        convr1dcircularinv(a, m, buf, m, r);
        return;
    }
    
    //
    // Task is normalized
    //
    if( m%2==0 )
    {
        
        //
        // size is even, use fast even-size FFT
        //
        buf.setlength(m);
        ap::vmove(&buf(0), 1, &a(0), 1, ap::vlen(0,m-1));
        buf2.setlength(m);
        ap::vmove(&buf2(0), 1, &b(0), 1, ap::vlen(0,n-1));
        for(i = n; i <= m-1; i++)
        {
            buf2(i) = 0;
        }
        buf3.setlength(m);
        ftbasegeneratecomplexfftplan(m/2, plan);
        fftr1dinternaleven(buf, m, buf3, plan);
        fftr1dinternaleven(buf2, m, buf3, plan);
        buf(0) = buf(0)/buf2(0);
        buf(1) = buf(1)/buf2(1);
        for(i = 1; i <= m/2-1; i++)
        {
            c1.x = buf(2*i+0);
            c1.y = buf(2*i+1);
            c2.x = buf2(2*i+0);
            c2.y = buf2(2*i+1);
            c3 = c1/c2;
            buf(2*i+0) = c3.x;
            buf(2*i+1) = c3.y;
        }
        fftr1dinvinternaleven(buf, m, buf3, plan);
        r.setlength(m);
        ap::vmove(&r(0), 1, &buf(0), 1, ap::vlen(0,m-1));
    }
    else
    {
        
        //
        // odd-size, use general real FFT
        //
        fftr1d(a, m, cbuf);
        buf2.setlength(m);
        ap::vmove(&buf2(0), 1, &b(0), 1, ap::vlen(0,n-1));
        for(i = n; i <= m-1; i++)
        {
            buf2(i) = 0;
        }
        fftr1d(buf2, m, cbuf2);
        for(i = 0; i <= ap::ifloor(double(m)/double(2)); i++)
        {
            cbuf(i) = cbuf(i)/cbuf2(i);
        }
        fftr1dinv(cbuf, m, r);
    }
}


/*************************************************************************
1-dimensional complex convolution.

Extended subroutine which allows to choose convolution algorithm.
Intended for internal use, ALGLIB users should call ConvC1D()/ConvC1DCircular().

INPUT PARAMETERS
    A   -   array[0..M-1] - complex function to be transformed
    M   -   problem size
    B   -   array[0..N-1] - complex function to be transformed
    N   -   problem size, N<=M
    Alg -   algorithm type:
            *-2     auto-select Q for overlap-add
            *-1     auto-select algorithm and parameters
            * 0     straightforward formula for small N's
            * 1     general FFT-based code
            * 2     overlap-add with length Q
    Q   -   length for overlap-add

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..N+M-1].

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convc1dx(const ap::complex_1d_array& a,
     int m,
     const ap::complex_1d_array& b,
     int n,
     bool circular,
     int alg,
     int q,
     ap::complex_1d_array& r)
{
    int i;
    int j;
    int p;
    int ptotal;
    int i1;
    int i2;
    int j1;
    int j2;
    ap::complex_1d_array bbuf;
    ap::complex v;
    double ax;
    double ay;
    double bx;
    double by;
    double t;
    double tx;
    double ty;
    double flopcand;
    double flopbest;
    int algbest;
    ftplan plan;
    ap::real_1d_array buf;
    ap::real_1d_array buf2;

    ap::ap_error::make_assertion(n>0&&m>0, "ConvC1DX: incorrect N or M!");
    ap::ap_error::make_assertion(n<=m, "ConvC1DX: N<M assumption is false!");
    
    //
    // Auto-select
    //
    if( alg==-1||alg==-2 )
    {
        
        //
        // Initial candidate: straightforward implementation.
        //
        // If we want to use auto-fitted overlap-add,
        // flop count is initialized by large real number - to force
        // another algorithm selection
        //
        algbest = 0;
        if( alg==-1 )
        {
            flopbest = 2*m*n;
        }
        else
        {
            flopbest = ap::maxrealnumber;
        }
        
        //
        // Another candidate - generic FFT code
        //
        if( alg==-1 )
        {
            if( circular&&ftbaseissmooth(m) )
            {
                
                //
                // special code for circular convolution of a sequence with a smooth length
                //
                flopcand = 3*ftbasegetflopestimate(m)+6*m;
                if( ap::fp_less(flopcand,flopbest) )
                {
                    algbest = 1;
                    flopbest = flopcand;
                }
            }
            else
            {
                
                //
                // general cyclic/non-cyclic convolution
                //
                p = ftbasefindsmooth(m+n-1);
                flopcand = 3*ftbasegetflopestimate(p)+6*p;
                if( ap::fp_less(flopcand,flopbest) )
                {
                    algbest = 1;
                    flopbest = flopcand;
                }
            }
        }
        
        //
        // Another candidate - overlap-add
        //
        q = 1;
        ptotal = 1;
        while(ptotal<n)
        {
            ptotal = ptotal*2;
        }
        while(ptotal<=m+n-1)
        {
            p = ptotal-n+1;
            flopcand = ap::iceil(double(m)/double(p))*(2*ftbasegetflopestimate(ptotal)+8*ptotal);
            if( ap::fp_less(flopcand,flopbest) )
            {
                flopbest = flopcand;
                algbest = 2;
                q = p;
            }
            ptotal = ptotal*2;
        }
        alg = algbest;
        convc1dx(a, m, b, n, circular, alg, q, r);
        return;
    }
    
    //
    // straightforward formula for
    // circular and non-circular convolutions.
    //
    // Very simple code, no further comments needed.
    //
    if( alg==0 )
    {
        
        //
        // Special case: N=1
        //
        if( n==1 )
        {
            r.setlength(m);
            v = b(0);
            ap::vmove(&r(0), 1, &a(0), 1, "N", ap::vlen(0,m-1), v);
            return;
        }
        
        //
        // use straightforward formula
        //
        if( circular )
        {
            
            //
            // circular convolution
            //
            r.setlength(m);
            v = b(0);
            ap::vmove(&r(0), 1, &a(0), 1, "N", ap::vlen(0,m-1), v);
            for(i = 1; i <= n-1; i++)
            {
                v = b(i);
                i1 = 0;
                i2 = i-1;
                j1 = m-i;
                j2 = m-1;
                ap::vadd(&r(i1), 1, &a(j1), 1, "N", ap::vlen(i1,i2), v);
                i1 = i;
                i2 = m-1;
                j1 = 0;
                j2 = m-i-1;
                ap::vadd(&r(i1), 1, &a(j1), 1, "N", ap::vlen(i1,i2), v);
            }
        }
        else
        {
            
            //
            // non-circular convolution
            //
            r.setlength(m+n-1);
            for(i = 0; i <= m+n-2; i++)
            {
                r(i) = 0;
            }
            for(i = 0; i <= n-1; i++)
            {
                v = b(i);
                ap::vadd(&r(i), 1, &a(0), 1, "N", ap::vlen(i,i+m-1), v);
            }
        }
        return;
    }
    
    //
    // general FFT-based code for
    // circular and non-circular convolutions.
    //
    // First, if convolution is circular, we test whether M is smooth or not.
    // If it is smooth, we just use M-length FFT to calculate convolution.
    // If it is not, we calculate non-circular convolution and wrap it arount.
    //
    // IF convolution is non-circular, we use zero-padding + FFT.
    //
    if( alg==1 )
    {
        if( circular&&ftbaseissmooth(m) )
        {
            
            //
            // special code for circular convolution with smooth M
            //
            ftbasegeneratecomplexfftplan(m, plan);
            buf.setlength(2*m);
            for(i = 0; i <= m-1; i++)
            {
                buf(2*i+0) = a(i).x;
                buf(2*i+1) = a(i).y;
            }
            buf2.setlength(2*m);
            for(i = 0; i <= n-1; i++)
            {
                buf2(2*i+0) = b(i).x;
                buf2(2*i+1) = b(i).y;
            }
            for(i = n; i <= m-1; i++)
            {
                buf2(2*i+0) = 0;
                buf2(2*i+1) = 0;
            }
            ftbaseexecuteplan(buf, 0, m, plan);
            ftbaseexecuteplan(buf2, 0, m, plan);
            for(i = 0; i <= m-1; i++)
            {
                ax = buf(2*i+0);
                ay = buf(2*i+1);
                bx = buf2(2*i+0);
                by = buf2(2*i+1);
                tx = ax*bx-ay*by;
                ty = ax*by+ay*bx;
                buf(2*i+0) = tx;
                buf(2*i+1) = -ty;
            }
            ftbaseexecuteplan(buf, 0, m, plan);
            t = double(1)/double(m);
            r.setlength(m);
            for(i = 0; i <= m-1; i++)
            {
                r(i).x = +t*buf(2*i+0);
                r(i).y = -t*buf(2*i+1);
            }
        }
        else
        {
            
            //
            // M is non-smooth, general code (circular/non-circular):
            // * first part is the same for circular and non-circular
            //   convolutions. zero padding, FFTs, inverse FFTs
            // * second part differs:
            //   * for non-circular convolution we just copy array
            //   * for circular convolution we add array tail to its head
            //
            p = ftbasefindsmooth(m+n-1);
            ftbasegeneratecomplexfftplan(p, plan);
            buf.setlength(2*p);
            for(i = 0; i <= m-1; i++)
            {
                buf(2*i+0) = a(i).x;
                buf(2*i+1) = a(i).y;
            }
            for(i = m; i <= p-1; i++)
            {
                buf(2*i+0) = 0;
                buf(2*i+1) = 0;
            }
            buf2.setlength(2*p);
            for(i = 0; i <= n-1; i++)
            {
                buf2(2*i+0) = b(i).x;
                buf2(2*i+1) = b(i).y;
            }
            for(i = n; i <= p-1; i++)
            {
                buf2(2*i+0) = 0;
                buf2(2*i+1) = 0;
            }
            ftbaseexecuteplan(buf, 0, p, plan);
            ftbaseexecuteplan(buf2, 0, p, plan);
            for(i = 0; i <= p-1; i++)
            {
                ax = buf(2*i+0);
                ay = buf(2*i+1);
                bx = buf2(2*i+0);
                by = buf2(2*i+1);
                tx = ax*bx-ay*by;
                ty = ax*by+ay*bx;
                buf(2*i+0) = tx;
                buf(2*i+1) = -ty;
            }
            ftbaseexecuteplan(buf, 0, p, plan);
            t = double(1)/double(p);
            if( circular )
            {
                
                //
                // circular, add tail to head
                //
                r.setlength(m);
                for(i = 0; i <= m-1; i++)
                {
                    r(i).x = +t*buf(2*i+0);
                    r(i).y = -t*buf(2*i+1);
                }
                for(i = m; i <= m+n-2; i++)
                {
                    r(i-m).x = r(i-m).x+t*buf(2*i+0);
                    r(i-m).y = r(i-m).y-t*buf(2*i+1);
                }
            }
            else
            {
                
                //
                // non-circular, just copy
                //
                r.setlength(m+n-1);
                for(i = 0; i <= m+n-2; i++)
                {
                    r(i).x = +t*buf(2*i+0);
                    r(i).y = -t*buf(2*i+1);
                }
            }
        }
        return;
    }
    
    //
    // overlap-add method for
    // circular and non-circular convolutions.
    //
    // First part of code (separate FFTs of input blocks) is the same
    // for all types of convolution. Second part (overlapping outputs)
    // differs for different types of convolution. We just copy output
    // when convolution is non-circular. We wrap it around, if it is
    // circular.
    //
    if( alg==2 )
    {
        buf.setlength(2*(q+n-1));
        
        //
        // prepare R
        //
        if( circular )
        {
            r.setlength(m);
            for(i = 0; i <= m-1; i++)
            {
                r(i) = 0;
            }
        }
        else
        {
            r.setlength(m+n-1);
            for(i = 0; i <= m+n-2; i++)
            {
                r(i) = 0;
            }
        }
        
        //
        // pre-calculated FFT(B)
        //
        bbuf.setlength(q+n-1);
        ap::vmove(&bbuf(0), 1, &b(0), 1, "N", ap::vlen(0,n-1));
        for(j = n; j <= q+n-2; j++)
        {
            bbuf(j) = 0;
        }
        fftc1d(bbuf, q+n-1);
        
        //
        // prepare FFT plan for chunks of A
        //
        ftbasegeneratecomplexfftplan(q+n-1, plan);
        
        //
        // main overlap-add cycle
        //
        i = 0;
        while(i<=m-1)
        {
            p = ap::minint(q, m-i);
            for(j = 0; j <= p-1; j++)
            {
                buf(2*j+0) = a(i+j).x;
                buf(2*j+1) = a(i+j).y;
            }
            for(j = p; j <= q+n-2; j++)
            {
                buf(2*j+0) = 0;
                buf(2*j+1) = 0;
            }
            ftbaseexecuteplan(buf, 0, q+n-1, plan);
            for(j = 0; j <= q+n-2; j++)
            {
                ax = buf(2*j+0);
                ay = buf(2*j+1);
                bx = bbuf(j).x;
                by = bbuf(j).y;
                tx = ax*bx-ay*by;
                ty = ax*by+ay*bx;
                buf(2*j+0) = tx;
                buf(2*j+1) = -ty;
            }
            ftbaseexecuteplan(buf, 0, q+n-1, plan);
            t = double(1)/double(q+n-1);
            if( circular )
            {
                j1 = ap::minint(i+p+n-2, m-1)-i;
                j2 = j1+1;
            }
            else
            {
                j1 = p+n-2;
                j2 = j1+1;
            }
            for(j = 0; j <= j1; j++)
            {
                r(i+j).x = r(i+j).x+buf(2*j+0)*t;
                r(i+j).y = r(i+j).y-buf(2*j+1)*t;
            }
            for(j = j2; j <= p+n-2; j++)
            {
                r(j-j2).x = r(j-j2).x+buf(2*j+0)*t;
                r(j-j2).y = r(j-j2).y-buf(2*j+1)*t;
            }
            i = i+p;
        }
        return;
    }
}


/*************************************************************************
1-dimensional real convolution.

Extended subroutine which allows to choose convolution algorithm.
Intended for internal use, ALGLIB users should call ConvR1D().

INPUT PARAMETERS
    A   -   array[0..M-1] - complex function to be transformed
    M   -   problem size
    B   -   array[0..N-1] - complex function to be transformed
    N   -   problem size, N<=M
    Alg -   algorithm type:
            *-2     auto-select Q for overlap-add
            *-1     auto-select algorithm and parameters
            * 0     straightforward formula for small N's
            * 1     general FFT-based code
            * 2     overlap-add with length Q
    Q   -   length for overlap-add

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..N+M-1].

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void convr1dx(const ap::real_1d_array& a,
     int m,
     const ap::real_1d_array& b,
     int n,
     bool circular,
     int alg,
     int q,
     ap::real_1d_array& r)
{
    double v;
    int i;
    int j;
    int p;
    int ptotal;
    int i1;
    int i2;
    int j1;
    int j2;
    double ax;
    double ay;
    double bx;
    double by;
    double tx;
    double ty;
    double flopcand;
    double flopbest;
    int algbest;
    ftplan plan;
    ap::real_1d_array buf;
    ap::real_1d_array buf2;
    ap::real_1d_array buf3;

    ap::ap_error::make_assertion(n>0&&m>0, "ConvC1DX: incorrect N or M!");
    ap::ap_error::make_assertion(n<=m, "ConvC1DX: N<M assumption is false!");
    
    //
    // handle special cases
    //
    if( ap::minint(m, n)<=2 )
    {
        alg = 0;
    }
    
    //
    // Auto-select
    //
    if( alg<0 )
    {
        
        //
        // Initial candidate: straightforward implementation.
        //
        // If we want to use auto-fitted overlap-add,
        // flop count is initialized by large real number - to force
        // another algorithm selection
        //
        algbest = 0;
        if( alg==-1 )
        {
            flopbest = 0.15*m*n;
        }
        else
        {
            flopbest = ap::maxrealnumber;
        }
        
        //
        // Another candidate - generic FFT code
        //
        if( alg==-1 )
        {
            if( circular&&ftbaseissmooth(m)&&m%2==0 )
            {
                
                //
                // special code for circular convolution of a sequence with a smooth length
                //
                flopcand = 3*ftbasegetflopestimate(m/2)+double(6*m)/double(2);
                if( ap::fp_less(flopcand,flopbest) )
                {
                    algbest = 1;
                    flopbest = flopcand;
                }
            }
            else
            {
                
                //
                // general cyclic/non-cyclic convolution
                //
                p = ftbasefindsmootheven(m+n-1);
                flopcand = 3*ftbasegetflopestimate(p/2)+double(6*p)/double(2);
                if( ap::fp_less(flopcand,flopbest) )
                {
                    algbest = 1;
                    flopbest = flopcand;
                }
            }
        }
        
        //
        // Another candidate - overlap-add
        //
        q = 1;
        ptotal = 1;
        while(ptotal<n)
        {
            ptotal = ptotal*2;
        }
        while(ptotal<=m+n-1)
        {
            p = ptotal-n+1;
            flopcand = ap::iceil(double(m)/double(p))*(2*ftbasegetflopestimate(ptotal/2)+1*(ptotal/2));
            if( ap::fp_less(flopcand,flopbest) )
            {
                flopbest = flopcand;
                algbest = 2;
                q = p;
            }
            ptotal = ptotal*2;
        }
        alg = algbest;
        convr1dx(a, m, b, n, circular, alg, q, r);
        return;
    }
    
    //
    // straightforward formula for
    // circular and non-circular convolutions.
    //
    // Very simple code, no further comments needed.
    //
    if( alg==0 )
    {
        
        //
        // Special case: N=1
        //
        if( n==1 )
        {
            r.setlength(m);
            v = b(0);
            ap::vmove(&r(0), 1, &a(0), 1, ap::vlen(0,m-1), v);
            return;
        }
        
        //
        // use straightforward formula
        //
        if( circular )
        {
            
            //
            // circular convolution
            //
            r.setlength(m);
            v = b(0);
            ap::vmove(&r(0), 1, &a(0), 1, ap::vlen(0,m-1), v);
            for(i = 1; i <= n-1; i++)
            {
                v = b(i);
                i1 = 0;
                i2 = i-1;
                j1 = m-i;
                j2 = m-1;
                ap::vadd(&r(i1), 1, &a(j1), 1, ap::vlen(i1,i2), v);
                i1 = i;
                i2 = m-1;
                j1 = 0;
                j2 = m-i-1;
                ap::vadd(&r(i1), 1, &a(j1), 1, ap::vlen(i1,i2), v);
            }
        }
        else
        {
            
            //
            // non-circular convolution
            //
            r.setlength(m+n-1);
            for(i = 0; i <= m+n-2; i++)
            {
                r(i) = 0;
            }
            for(i = 0; i <= n-1; i++)
            {
                v = b(i);
                ap::vadd(&r(i), 1, &a(0), 1, ap::vlen(i,i+m-1), v);
            }
        }
        return;
    }
    
    //
    // general FFT-based code for
    // circular and non-circular convolutions.
    //
    // First, if convolution is circular, we test whether M is smooth or not.
    // If it is smooth, we just use M-length FFT to calculate convolution.
    // If it is not, we calculate non-circular convolution and wrap it arount.
    //
    // If convolution is non-circular, we use zero-padding + FFT.
    //
    // We assume that M+N-1>2 - we should call small case code otherwise
    //
    if( alg==1 )
    {
        ap::ap_error::make_assertion(m+n-1>2, "ConvR1DX: internal error!");
        if( circular&&ftbaseissmooth(m)&&m%2==0 )
        {
            
            //
            // special code for circular convolution with smooth even M
            //
            buf.setlength(m);
            ap::vmove(&buf(0), 1, &a(0), 1, ap::vlen(0,m-1));
            buf2.setlength(m);
            ap::vmove(&buf2(0), 1, &b(0), 1, ap::vlen(0,n-1));
            for(i = n; i <= m-1; i++)
            {
                buf2(i) = 0;
            }
            buf3.setlength(m);
            ftbasegeneratecomplexfftplan(m/2, plan);
            fftr1dinternaleven(buf, m, buf3, plan);
            fftr1dinternaleven(buf2, m, buf3, plan);
            buf(0) = buf(0)*buf2(0);
            buf(1) = buf(1)*buf2(1);
            for(i = 1; i <= m/2-1; i++)
            {
                ax = buf(2*i+0);
                ay = buf(2*i+1);
                bx = buf2(2*i+0);
                by = buf2(2*i+1);
                tx = ax*bx-ay*by;
                ty = ax*by+ay*bx;
                buf(2*i+0) = tx;
                buf(2*i+1) = ty;
            }
            fftr1dinvinternaleven(buf, m, buf3, plan);
            r.setlength(m);
            ap::vmove(&r(0), 1, &buf(0), 1, ap::vlen(0,m-1));
        }
        else
        {
            
            //
            // M is non-smooth or non-even, general code (circular/non-circular):
            // * first part is the same for circular and non-circular
            //   convolutions. zero padding, FFTs, inverse FFTs
            // * second part differs:
            //   * for non-circular convolution we just copy array
            //   * for circular convolution we add array tail to its head
            //
            p = ftbasefindsmootheven(m+n-1);
            buf.setlength(p);
            ap::vmove(&buf(0), 1, &a(0), 1, ap::vlen(0,m-1));
            for(i = m; i <= p-1; i++)
            {
                buf(i) = 0;
            }
            buf2.setlength(p);
            ap::vmove(&buf2(0), 1, &b(0), 1, ap::vlen(0,n-1));
            for(i = n; i <= p-1; i++)
            {
                buf2(i) = 0;
            }
            buf3.setlength(p);
            ftbasegeneratecomplexfftplan(p/2, plan);
            fftr1dinternaleven(buf, p, buf3, plan);
            fftr1dinternaleven(buf2, p, buf3, plan);
            buf(0) = buf(0)*buf2(0);
            buf(1) = buf(1)*buf2(1);
            for(i = 1; i <= p/2-1; i++)
            {
                ax = buf(2*i+0);
                ay = buf(2*i+1);
                bx = buf2(2*i+0);
                by = buf2(2*i+1);
                tx = ax*bx-ay*by;
                ty = ax*by+ay*bx;
                buf(2*i+0) = tx;
                buf(2*i+1) = ty;
            }
            fftr1dinvinternaleven(buf, p, buf3, plan);
            if( circular )
            {
                
                //
                // circular, add tail to head
                //
                r.setlength(m);
                ap::vmove(&r(0), 1, &buf(0), 1, ap::vlen(0,m-1));
                if( n>=2 )
                {
                    ap::vadd(&r(0), 1, &buf(m), 1, ap::vlen(0,n-2));
                }
            }
            else
            {
                
                //
                // non-circular, just copy
                //
                r.setlength(m+n-1);
                ap::vmove(&r(0), 1, &buf(0), 1, ap::vlen(0,m+n-2));
            }
        }
        return;
    }
    
    //
    // overlap-add method
    //
    if( alg==2 )
    {
        ap::ap_error::make_assertion((q+n-1)%2==0, "ConvR1DX: internal error!");
        buf.setlength(q+n-1);
        buf2.setlength(q+n-1);
        buf3.setlength(q+n-1);
        ftbasegeneratecomplexfftplan((q+n-1)/2, plan);
        
        //
        // prepare R
        //
        if( circular )
        {
            r.setlength(m);
            for(i = 0; i <= m-1; i++)
            {
                r(i) = 0;
            }
        }
        else
        {
            r.setlength(m+n-1);
            for(i = 0; i <= m+n-2; i++)
            {
                r(i) = 0;
            }
        }
        
        //
        // pre-calculated FFT(B)
        //
        ap::vmove(&buf2(0), 1, &b(0), 1, ap::vlen(0,n-1));
        for(j = n; j <= q+n-2; j++)
        {
            buf2(j) = 0;
        }
        fftr1dinternaleven(buf2, q+n-1, buf3, plan);
        
        //
        // main overlap-add cycle
        //
        i = 0;
        while(i<=m-1)
        {
            p = ap::minint(q, m-i);
            ap::vmove(&buf(0), 1, &a(i), 1, ap::vlen(0,p-1));
            for(j = p; j <= q+n-2; j++)
            {
                buf(j) = 0;
            }
            fftr1dinternaleven(buf, q+n-1, buf3, plan);
            buf(0) = buf(0)*buf2(0);
            buf(1) = buf(1)*buf2(1);
            for(j = 1; j <= (q+n-1)/2-1; j++)
            {
                ax = buf(2*j+0);
                ay = buf(2*j+1);
                bx = buf2(2*j+0);
                by = buf2(2*j+1);
                tx = ax*bx-ay*by;
                ty = ax*by+ay*bx;
                buf(2*j+0) = tx;
                buf(2*j+1) = ty;
            }
            fftr1dinvinternaleven(buf, q+n-1, buf3, plan);
            if( circular )
            {
                j1 = ap::minint(i+p+n-2, m-1)-i;
                j2 = j1+1;
            }
            else
            {
                j1 = p+n-2;
                j2 = j1+1;
            }
            ap::vadd(&r(i), 1, &buf(0), 1, ap::vlen(i,i+j1));
            if( p+n-2>=j2 )
            {
                ap::vadd(&r(0), 1, &buf(j2), 1, ap::vlen(0,p+n-2-j2));
            }
            i = i+p;
        }
        return;
    }
}




