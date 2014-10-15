/*************************************************************************
Copyright (c)
    2007, Sergey Bochkanov (ALGLIB project).
    1988, Pierre L'Ecuyer

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
#include "hqrnd.h"

static const int hqrndmax = 2147483563;
static const int hqrndm1 = 2147483563;
static const int hqrndm2 = 2147483399;
static const int hqrndmagic = 1634357784;

static int hqrndintegerbase(hqrndstate& state);

/*************************************************************************
HQRNDState  initialization  with  random  values  which come from standard
RNG.

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void hqrndrandomize(hqrndstate& state)
{

    hqrndseed(ap::randominteger(hqrndm1), ap::randominteger(hqrndm2), state);
}


/*************************************************************************
HQRNDState initialization with seed values

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void hqrndseed(int s1, int s2, hqrndstate& state)
{

    state.s1 = s1%(hqrndm1-1)+1;
    state.s2 = s2%(hqrndm2-1)+1;
    state.v = double(1)/double(hqrndmax);
    state.magicv = hqrndmagic;
}


/*************************************************************************
This function generates random real number in (0,1),
not including interval boundaries

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
double hqrnduniformr(hqrndstate& state)
{
    double result;

    result = state.v*hqrndintegerbase(state);
    return result;
}


/*************************************************************************
This function generates random integer number in [0, N)

1. N must be less than HQRNDMax-1.
2. State structure must be initialized with HQRNDRandomize() or HQRNDSeed()

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
int hqrnduniformi(int n, hqrndstate& state)
{
    int result;
    int mx;

    
    //
    // Correct handling of N's close to RNDBaseMax
    // (avoiding skewed distributions for RNDBaseMax<>K*N)
    //
    ap::ap_error::make_assertion(n>0, "HQRNDUniformI: N<=0!");
    ap::ap_error::make_assertion(n<hqrndmax-1, "HQRNDUniformI: N>=RNDBaseMax-1!");
    mx = hqrndmax-1-(hqrndmax-1)%n;
    do
    {
        result = hqrndintegerbase(state)-1;
    }
    while(result>=mx);
    result = result%n;
    return result;
}


/*************************************************************************
Random number generator: normal numbers

This function generates one random number from normal distribution.
Its performance is equal to that of HQRNDNormal2()

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
double hqrndnormal(hqrndstate& state)
{
    double result;
    double v1;
    double v2;

    hqrndnormal2(state, v1, v2);
    result = v1;
    return result;
}


/*************************************************************************
Random number generator: random X and Y such that X^2+Y^2=1

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void hqrndunit2(hqrndstate& state, double& x, double& y)
{
    double v;
    double mx;
    double mn;

    do
    {
        hqrndnormal2(state, x, y);
    }
    while(!(ap::fp_neq(x,0)||ap::fp_neq(y,0)));
    mx = ap::maxreal(fabs(x), fabs(y));
    mn = ap::minreal(fabs(x), fabs(y));
    v = mx*sqrt(1+ap::sqr(mn/mx));
    x = x/v;
    y = y/v;
}


/*************************************************************************
Random number generator: normal numbers

This function generates two independent random numbers from normal
distribution. Its performance is equal to that of HQRNDNormal()

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void hqrndnormal2(hqrndstate& state, double& x1, double& x2)
{
    double u;
    double v;
    double s;

    while(true)
    {
        u = 2*hqrnduniformr(state)-1;
        v = 2*hqrnduniformr(state)-1;
        s = ap::sqr(u)+ap::sqr(v);
        if( ap::fp_greater(s,0)&&ap::fp_less(s,1) )
        {
            
            //
            // two Sqrt's instead of one to
            // avoid overflow when S is too small
            //
            s = sqrt(-2*log(s))/sqrt(s);
            x1 = u*s;
            x2 = v*s;
            return;
        }
    }
}


/*************************************************************************
Random number generator: exponential distribution

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 11.08.2007 by Bochkanov Sergey
*************************************************************************/
double hqrndexponential(double lambda, hqrndstate& state)
{
    double result;

    ap::ap_error::make_assertion(ap::fp_greater(lambda,0), "HQRNDExponential: Lambda<=0!");
    result = -log(hqrnduniformr(state))/lambda;
    return result;
}


/*************************************************************************

L'Ecuyer, Efficient and portable combined random number generators
*************************************************************************/
static int hqrndintegerbase(hqrndstate& state)
{
    int result;
    int k;

    ap::ap_error::make_assertion(state.magicv==hqrndmagic, "HQRNDIntegerBase: State is not correctly initialized!");
    k = state.s1/53668;
    state.s1 = 40014*(state.s1-k*53668)-k*12211;
    if( state.s1<0 )
    {
        state.s1 = state.s1+2147483563;
    }
    k = state.s2/52774;
    state.s2 = 40692*(state.s2-k*52774)-k*3791;
    if( state.s2<0 )
    {
        state.s2 = state.s2+2147483399;
    }
    
    //
    // Result
    //
    result = state.s1-state.s2;
    if( result<1 )
    {
        result = result+2147483562;
    }
    return result;
}




