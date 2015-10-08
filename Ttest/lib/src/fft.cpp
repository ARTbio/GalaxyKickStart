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
#include "fft.h"

/*************************************************************************
1-dimensional complex FFT.

Array size N may be arbitrary number (composite or prime).  Composite  N's
are handled with cache-oblivious variation of  a  Cooley-Tukey  algorithm.
Small prime-factors are transformed using hard coded  codelets (similar to
FFTW codelets, but without low-level  optimization),  large  prime-factors
are handled with Bluestein's algorithm.

Fastests transforms are for smooth N's (prime factors are 2, 3,  5  only),
most fast for powers of 2. When N have prime factors  larger  than  these,
but orders of magnitude smaller than N, computations will be about 4 times
slower than for nearby highly composite N's. When N itself is prime, speed
will be 6 times lower.

Algorithm has O(N*logN) complexity for any N (composite or prime).

INPUT PARAMETERS
    A   -   array[0..N-1] - complex function to be transformed
    N   -   problem size
    
OUTPUT PARAMETERS
    A   -   DFT of a input array, array[0..N-1]
            A_out[j] = SUM(A_in[k]*exp(-2*pi*sqrt(-1)*j*k/N), k = 0..N-1)


  -- ALGLIB --
     Copyright 29.05.2009 by Bochkanov Sergey
*************************************************************************/
void fftc1d(ap::complex_1d_array& a, int n)
{
    ftplan plan;
    int i;
    ap::real_1d_array buf;

    ap::ap_error::make_assertion(n>0, "FFTC1D: incorrect N!");
    
    //
    // Special case: N=1, FFT is just identity transform.
    // After this block we assume that N is strictly greater than 1.
    //
    if( n==1 )
    {
        return;
    }
    
    //
    // convert input array to the more convinient format
    //
    buf.setlength(2*n);
    for(i = 0; i <= n-1; i++)
    {
        buf(2*i+0) = a(i).x;
        buf(2*i+1) = a(i).y;
    }
    
    //
    // Generate plan and execute it.
    //
    // Plan is a combination of a successive factorizations of N and
    // precomputed data. It is much like a FFTW plan, but is not stored
    // between subroutine calls and is much simpler.
    //
    ftbasegeneratecomplexfftplan(n, plan);
    ftbaseexecuteplan(buf, 0, n, plan);
    
    //
    // result
    //
    for(i = 0; i <= n-1; i++)
    {
        a(i).x = buf(2*i+0);
        a(i).y = buf(2*i+1);
    }
}


/*************************************************************************
1-dimensional complex inverse FFT.

Array size N may be arbitrary number (composite or prime).  Algorithm  has
O(N*logN) complexity for any N (composite or prime).

See FFTC1D() description for more information about algorithm performance.

INPUT PARAMETERS
    A   -   array[0..N-1] - complex array to be transformed
    N   -   problem size

OUTPUT PARAMETERS
    A   -   inverse DFT of a input array, array[0..N-1]
            A_out[j] = SUM(A_in[k]/N*exp(+2*pi*sqrt(-1)*j*k/N), k = 0..N-1)


  -- ALGLIB --
     Copyright 29.05.2009 by Bochkanov Sergey
*************************************************************************/
void fftc1dinv(ap::complex_1d_array& a, int n)
{
    int i;

    ap::ap_error::make_assertion(n>0, "FFTC1DInv: incorrect N!");
    
    //
    // Inverse DFT can be expressed in terms of the DFT as
    //
    //     invfft(x) = fft(x')'/N
    //
    // here x' means conj(x).
    //
    for(i = 0; i <= n-1; i++)
    {
        a(i).y = -a(i).y;
    }
    fftc1d(a, n);
    for(i = 0; i <= n-1; i++)
    {
        a(i).x = a(i).x/n;
        a(i).y = -a(i).y/n;
    }
}


/*************************************************************************
1-dimensional real FFT.

Algorithm has O(N*logN) complexity for any N (composite or prime).

INPUT PARAMETERS
    A   -   array[0..N-1] - real function to be transformed
    N   -   problem size

OUTPUT PARAMETERS
    F   -   DFT of a input array, array[0..N-1]
            F[j] = SUM(A[k]*exp(-2*pi*sqrt(-1)*j*k/N), k = 0..N-1)

NOTE:
    F[] satisfies symmetry property F[k] = conj(F[N-k]),  so just one half
of  array  is  usually needed. But for convinience subroutine returns full
complex array (with frequencies above N/2), so its result may be  used  by
other FFT-related subroutines.


  -- ALGLIB --
     Copyright 01.06.2009 by Bochkanov Sergey
*************************************************************************/
void fftr1d(const ap::real_1d_array& a, int n, ap::complex_1d_array& f)
{
    int i;
    int n2;
    int idx;
    ap::complex hn;
    ap::complex hmnc;
    ap::complex v;
    ap::real_1d_array buf;
    ftplan plan;

    ap::ap_error::make_assertion(n>0, "FFTR1D: incorrect N!");
    
    //
    // Special cases:
    // * N=1, FFT is just identity transform.
    // * N=2, FFT is simple too
    //
    // After this block we assume that N is strictly greater than 2
    //
    if( n==1 )
    {
        f.setlength(1);
        f(0) = a(0);
        return;
    }
    if( n==2 )
    {
        f.setlength(2);
        f(0).x = a(0)+a(1);
        f(0).y = 0;
        f(1).x = a(0)-a(1);
        f(1).y = 0;
        return;
    }
    
    //
    // Choose between odd-size and even-size FFTs
    //
    if( n%2==0 )
    {
        
        //
        // even-size real FFT, use reduction to the complex task
        //
        n2 = n/2;
        buf.setlength(n);
        ap::vmove(&buf(0), 1, &a(0), 1, ap::vlen(0,n-1));
        ftbasegeneratecomplexfftplan(n2, plan);
        ftbaseexecuteplan(buf, 0, n2, plan);
        f.setlength(n);
        for(i = 0; i <= n2; i++)
        {
            idx = 2*(i%n2);
            hn.x = buf(idx+0);
            hn.y = buf(idx+1);
            idx = 2*((n2-i)%n2);
            hmnc.x = buf(idx+0);
            hmnc.y = -buf(idx+1);
            v.x = -sin(-2*ap::pi()*i/n);
            v.y = cos(-2*ap::pi()*i/n);
            f(i) = hn+hmnc-v*(hn-hmnc);
            f(i).x = 0.5*f(i).x;
            f(i).y = 0.5*f(i).y;
        }
        for(i = n2+1; i <= n-1; i++)
        {
            f(i) = ap::conj(f(n-i));
        }
        return;
    }
    else
    {
        
        //
        // use complex FFT
        //
        f.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            f(i) = a(i);
        }
        fftc1d(f, n);
        return;
    }
}


/*************************************************************************
1-dimensional real inverse FFT.

Algorithm has O(N*logN) complexity for any N (composite or prime).

INPUT PARAMETERS
    F   -   array[0..floor(N/2)] - frequencies from forward real FFT
    N   -   problem size

OUTPUT PARAMETERS
    A   -   inverse DFT of a input array, array[0..N-1]

NOTE:
    F[] should satisfy symmetry property F[k] = conj(F[N-k]), so just  one
half of frequencies array is needed - elements from 0 to floor(N/2).  F[0]
is ALWAYS real. If N is even F[floor(N/2)] is real too. If N is odd,  then
F[floor(N/2)] has no special properties.

Relying on properties noted above, FFTR1DInv subroutine uses only elements
from 0th to floor(N/2)-th. It ignores imaginary part of F[0],  and in case
N is even it ignores imaginary part of F[floor(N/2)] too.  So you can pass
either frequencies array with N elements or reduced array with roughly N/2
elements - subroutine will successfully transform both.


  -- ALGLIB --
     Copyright 01.06.2009 by Bochkanov Sergey
*************************************************************************/
void fftr1dinv(const ap::complex_1d_array& f, int n, ap::real_1d_array& a)
{
    int i;
    ap::real_1d_array h;
    ap::complex_1d_array fh;

    ap::ap_error::make_assertion(n>0, "FFTR1DInv: incorrect N!");
    
    //
    // Special case: N=1, FFT is just identity transform.
    // After this block we assume that N is strictly greater than 1.
    //
    if( n==1 )
    {
        a.setlength(1);
        a(0) = f(0).x;
        return;
    }
    
    //
    // inverse real FFT is reduced to the inverse real FHT,
    // which is reduced to the forward real FHT,
    // which is reduced to the forward real FFT.
    //
    // Don't worry, it is really compact and efficient reduction :)
    //
    h.setlength(n);
    a.setlength(n);
    h(0) = f(0).x;
    for(i = 1; i <= ap::ifloor(double(n)/double(2))-1; i++)
    {
        h(i) = f(i).x-f(i).y;
        h(n-i) = f(i).x+f(i).y;
    }
    if( n%2==0 )
    {
        h(ap::ifloor(double(n)/double(2))) = f(ap::ifloor(double(n)/double(2))).x;
    }
    else
    {
        h(ap::ifloor(double(n)/double(2))) = f(ap::ifloor(double(n)/double(2))).x-f(ap::ifloor(double(n)/double(2))).y;
        h(ap::ifloor(double(n)/double(2))+1) = f(ap::ifloor(double(n)/double(2))).x+f(ap::ifloor(double(n)/double(2))).y;
    }
    fftr1d(h, n, fh);
    for(i = 0; i <= n-1; i++)
    {
        a(i) = (fh(i).x-fh(i).y)/n;
    }
}


/*************************************************************************
Internal subroutine. Never call it directly!


  -- ALGLIB --
     Copyright 01.06.2009 by Bochkanov Sergey
*************************************************************************/
void fftr1dinternaleven(ap::real_1d_array& a,
     int n,
     ap::real_1d_array& buf,
     ftplan& plan)
{
    double x;
    double y;
    int i;
    int n2;
    int idx;
    ap::complex hn;
    ap::complex hmnc;
    ap::complex v;

    ap::ap_error::make_assertion(n>0&&n%2==0, "FFTR1DEvenInplace: incorrect N!");
    
    //
    // Special cases:
    // * N=2
    //
    // After this block we assume that N is strictly greater than 2
    //
    if( n==2 )
    {
        x = a(0)+a(1);
        y = a(0)-a(1);
        a(0) = x;
        a(1) = y;
        return;
    }
    
    //
    // even-size real FFT, use reduction to the complex task
    //
    n2 = n/2;
    ap::vmove(&buf(0), 1, &a(0), 1, ap::vlen(0,n-1));
    ftbaseexecuteplan(buf, 0, n2, plan);
    a(0) = buf(0)+buf(1);
    for(i = 1; i <= n2-1; i++)
    {
        idx = 2*(i%n2);
        hn.x = buf(idx+0);
        hn.y = buf(idx+1);
        idx = 2*(n2-i);
        hmnc.x = buf(idx+0);
        hmnc.y = -buf(idx+1);
        v.x = -sin(-2*ap::pi()*i/n);
        v.y = cos(-2*ap::pi()*i/n);
        v = hn+hmnc-v*(hn-hmnc);
        a(2*i+0) = 0.5*v.x;
        a(2*i+1) = 0.5*v.y;
    }
    a(1) = buf(0)-buf(1);
}


/*************************************************************************
Internal subroutine. Never call it directly!


  -- ALGLIB --
     Copyright 01.06.2009 by Bochkanov Sergey
*************************************************************************/
void fftr1dinvinternaleven(ap::real_1d_array& a,
     int n,
     ap::real_1d_array& buf,
     ftplan& plan)
{
    double x;
    double y;
    double t;
    int i;
    int n2;

    ap::ap_error::make_assertion(n>0&&n%2==0, "FFTR1DInvInternalEven: incorrect N!");
    
    //
    // Special cases:
    // * N=2
    //
    // After this block we assume that N is strictly greater than 2
    //
    if( n==2 )
    {
        x = 0.5*(a(0)+a(1));
        y = 0.5*(a(0)-a(1));
        a(0) = x;
        a(1) = y;
        return;
    }
    
    //
    // inverse real FFT is reduced to the inverse real FHT,
    // which is reduced to the forward real FHT,
    // which is reduced to the forward real FFT.
    //
    // Don't worry, it is really compact and efficient reduction :)
    //
    n2 = n/2;
    buf(0) = a(0);
    for(i = 1; i <= n2-1; i++)
    {
        x = a(2*i+0);
        y = a(2*i+1);
        buf(i) = x-y;
        buf(n-i) = x+y;
    }
    buf(n2) = a(1);
    fftr1dinternaleven(buf, n, a, plan);
    a(0) = buf(0)/n;
    t = double(1)/double(n);
    for(i = 1; i <= n2-1; i++)
    {
        x = buf(2*i+0);
        y = buf(2*i+1);
        a(i) = t*(x-y);
        a(n-i) = t*(x+y);
    }
    a(n2) = buf(1)/n;
}




