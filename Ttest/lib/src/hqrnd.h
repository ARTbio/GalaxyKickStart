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

#ifndef _hqrnd_h
#define _hqrnd_h

#include "ap.h"
#include "ialglib.h"

/*************************************************************************
Portable high quality random number generator state.
Initialized with HQRNDRandomize() or HQRNDSeed().

Fields:
    S1, S2      -   seed values
    V           -   precomputed value
    MagicV      -   'magic' value used to determine whether State structure
                    was correctly initialized.
*************************************************************************/
struct hqrndstate
{
    int s1;
    int s2;
    double v;
    int magicv;
};




/*************************************************************************
HQRNDState  initialization  with  random  values  which come from standard
RNG.

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void hqrndrandomize(hqrndstate& state);


/*************************************************************************
HQRNDState initialization with seed values

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void hqrndseed(int s1, int s2, hqrndstate& state);


/*************************************************************************
This function generates random real number in (0,1),
not including interval boundaries

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
double hqrnduniformr(hqrndstate& state);


/*************************************************************************
This function generates random integer number in [0, N)

1. N must be less than HQRNDMax-1.
2. State structure must be initialized with HQRNDRandomize() or HQRNDSeed()

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
int hqrnduniformi(int n, hqrndstate& state);


/*************************************************************************
Random number generator: normal numbers

This function generates one random number from normal distribution.
Its performance is equal to that of HQRNDNormal2()

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
double hqrndnormal(hqrndstate& state);


/*************************************************************************
Random number generator: random X and Y such that X^2+Y^2=1

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void hqrndunit2(hqrndstate& state, double& x, double& y);


/*************************************************************************
Random number generator: normal numbers

This function generates two independent random numbers from normal
distribution. Its performance is equal to that of HQRNDNormal()

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void hqrndnormal2(hqrndstate& state, double& x1, double& x2);


/*************************************************************************
Random number generator: exponential distribution

State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

  -- ALGLIB --
     Copyright 11.08.2007 by Bochkanov Sergey
*************************************************************************/
double hqrndexponential(double lambda, hqrndstate& state);


#endif

