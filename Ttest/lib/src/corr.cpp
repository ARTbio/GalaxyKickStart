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
#include "corr.h"

/*************************************************************************
1-dimensional complex cross-correlation.

For given Pattern/Signal returns corr(Pattern,Signal) (non-circular).

Correlation is calculated using reduction to  convolution.  Algorithm with
max(N,N)*log(max(N,N)) complexity is used (see  ConvC1D()  for  more  info
about performance).

IMPORTANT:
    for  historical reasons subroutine accepts its parameters in  reversed
    order: CorrC1D(Signal, Pattern) = Pattern x Signal (using  traditional
    definition of cross-correlation, denoting cross-correlation as "x").

INPUT PARAMETERS
    Signal  -   array[0..N-1] - complex function to be transformed,
                signal containing pattern
    N       -   problem size
    Pattern -   array[0..M-1] - complex function to be transformed,
                pattern to search withing signal
    M       -   problem size

OUTPUT PARAMETERS
    R       -   cross-correlation, array[0..N+M-2]:
                * positive lags are stored in R[0..N-1],
                  R[i] = sum(conj(pattern[j])*signal[i+j]
                * negative lags are stored in R[N..N+M-2],
                  R[N+M-1-i] = sum(conj(pattern[j])*signal[-i+j]

NOTE:
    It is assumed that pattern domain is [0..M-1].  If Pattern is non-zero
on [-K..M-1],  you can still use this subroutine, just shift result by K.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void corrc1d(const ap::complex_1d_array& signal,
     int n,
     const ap::complex_1d_array& pattern,
     int m,
     ap::complex_1d_array& r)
{
    ap::complex_1d_array p;
    ap::complex_1d_array b;
    int i;

    ap::ap_error::make_assertion(n>0&&m>0, "CorrC1D: incorrect N or M!");
    p.setlength(m);
    for(i = 0; i <= m-1; i++)
    {
        p(m-1-i) = ap::conj(pattern(i));
    }
    convc1d(p, m, signal, n, b);
    r.setlength(m+n-1);
    ap::vmove(&r(0), 1, &b(m-1), 1, "N", ap::vlen(0,n-1));
    if( m+n-2>=n )
    {
        ap::vmove(&r(n), 1, &b(0), 1, "N", ap::vlen(n,m+n-2));
    }
}


/*************************************************************************
1-dimensional circular complex cross-correlation.

For given Pattern/Signal returns corr(Pattern,Signal) (circular).
Algorithm has linearithmic complexity for any M/N.

IMPORTANT:
    for  historical reasons subroutine accepts its parameters in  reversed
    order:   CorrC1DCircular(Signal, Pattern) = Pattern x Signal    (using
    traditional definition of cross-correlation, denoting cross-correlation
    as "x").

INPUT PARAMETERS
    Signal  -   array[0..N-1] - complex function to be transformed,
                periodic signal containing pattern
    N       -   problem size
    Pattern -   array[0..M-1] - complex function to be transformed,
                non-periodic pattern to search withing signal
    M       -   problem size

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..M-1].


  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void corrc1dcircular(const ap::complex_1d_array& signal,
     int m,
     const ap::complex_1d_array& pattern,
     int n,
     ap::complex_1d_array& c)
{
    ap::complex_1d_array p;
    ap::complex_1d_array b;
    int i1;
    int i2;
    int i;
    int j2;

    ap::ap_error::make_assertion(n>0&&m>0, "ConvC1DCircular: incorrect N or M!");
    
    //
    // normalize task: make M>=N,
    // so A will be longer (at least - not shorter) that B.
    //
    if( m<n )
    {
        b.setlength(m);
        for(i1 = 0; i1 <= m-1; i1++)
        {
            b(i1) = 0;
        }
        i1 = 0;
        while(i1<n)
        {
            i2 = ap::minint(i1+m-1, n-1);
            j2 = i2-i1;
            ap::vadd(&b(0), 1, &pattern(i1), 1, "N", ap::vlen(0,j2));
            i1 = i1+m;
        }
        corrc1dcircular(signal, m, b, m, c);
        return;
    }
    
    //
    // Task is normalized
    //
    p.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        p(n-1-i) = ap::conj(pattern(i));
    }
    convc1dcircular(signal, m, p, n, b);
    c.setlength(m);
    ap::vmove(&c(0), 1, &b(n-1), 1, "N", ap::vlen(0,m-n));
    if( m-n+1<=m-1 )
    {
        ap::vmove(&c(m-n+1), 1, &b(0), 1, "N", ap::vlen(m-n+1,m-1));
    }
}


/*************************************************************************
1-dimensional real cross-correlation.

For given Pattern/Signal returns corr(Pattern,Signal) (non-circular).

Correlation is calculated using reduction to  convolution.  Algorithm with
max(N,N)*log(max(N,N)) complexity is used (see  ConvC1D()  for  more  info
about performance).

IMPORTANT:
    for  historical reasons subroutine accepts its parameters in  reversed
    order: CorrR1D(Signal, Pattern) = Pattern x Signal (using  traditional
    definition of cross-correlation, denoting cross-correlation as "x").

INPUT PARAMETERS
    Signal  -   array[0..N-1] - real function to be transformed,
                signal containing pattern
    N       -   problem size
    Pattern -   array[0..M-1] - real function to be transformed,
                pattern to search withing signal
    M       -   problem size

OUTPUT PARAMETERS
    R       -   cross-correlation, array[0..N+M-2]:
                * positive lags are stored in R[0..N-1],
                  R[i] = sum(pattern[j]*signal[i+j]
                * negative lags are stored in R[N..N+M-2],
                  R[N+M-1-i] = sum(pattern[j]*signal[-i+j]

NOTE:
    It is assumed that pattern domain is [0..M-1].  If Pattern is non-zero
on [-K..M-1],  you can still use this subroutine, just shift result by K.

  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void corrr1d(const ap::real_1d_array& signal,
     int n,
     const ap::real_1d_array& pattern,
     int m,
     ap::real_1d_array& r)
{
    ap::real_1d_array p;
    ap::real_1d_array b;
    int i;

    ap::ap_error::make_assertion(n>0&&m>0, "CorrR1D: incorrect N or M!");
    p.setlength(m);
    for(i = 0; i <= m-1; i++)
    {
        p(m-1-i) = pattern(i);
    }
    convr1d(p, m, signal, n, b);
    r.setlength(m+n-1);
    ap::vmove(&r(0), 1, &b(m-1), 1, ap::vlen(0,n-1));
    if( m+n-2>=n )
    {
        ap::vmove(&r(n), 1, &b(0), 1, ap::vlen(n,m+n-2));
    }
}


/*************************************************************************
1-dimensional circular real cross-correlation.

For given Pattern/Signal returns corr(Pattern,Signal) (circular).
Algorithm has linearithmic complexity for any M/N.

IMPORTANT:
    for  historical reasons subroutine accepts its parameters in  reversed
    order:   CorrR1DCircular(Signal, Pattern) = Pattern x Signal    (using
    traditional definition of cross-correlation, denoting cross-correlation
    as "x").

INPUT PARAMETERS
    Signal  -   array[0..N-1] - real function to be transformed,
                periodic signal containing pattern
    N       -   problem size
    Pattern -   array[0..M-1] - real function to be transformed,
                non-periodic pattern to search withing signal
    M       -   problem size

OUTPUT PARAMETERS
    R   -   convolution: A*B. array[0..M-1].


  -- ALGLIB --
     Copyright 21.07.2009 by Bochkanov Sergey
*************************************************************************/
void corrr1dcircular(const ap::real_1d_array& signal,
     int m,
     const ap::real_1d_array& pattern,
     int n,
     ap::real_1d_array& c)
{
    ap::real_1d_array p;
    ap::real_1d_array b;
    int i1;
    int i2;
    int i;
    int j2;

    ap::ap_error::make_assertion(n>0&&m>0, "ConvC1DCircular: incorrect N or M!");
    
    //
    // normalize task: make M>=N,
    // so A will be longer (at least - not shorter) that B.
    //
    if( m<n )
    {
        b.setlength(m);
        for(i1 = 0; i1 <= m-1; i1++)
        {
            b(i1) = 0;
        }
        i1 = 0;
        while(i1<n)
        {
            i2 = ap::minint(i1+m-1, n-1);
            j2 = i2-i1;
            ap::vadd(&b(0), 1, &pattern(i1), 1, ap::vlen(0,j2));
            i1 = i1+m;
        }
        corrr1dcircular(signal, m, b, m, c);
        return;
    }
    
    //
    // Task is normalized
    //
    p.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        p(n-1-i) = pattern(i);
    }
    convr1dcircular(signal, m, p, n, b);
    c.setlength(m);
    ap::vmove(&c(0), 1, &b(n-1), 1, ap::vlen(0,m-n));
    if( m-n+1<=m-1 )
    {
        ap::vmove(&c(m-n+1), 1, &b(0), 1, ap::vlen(m-n+1,m-1));
    }
}




