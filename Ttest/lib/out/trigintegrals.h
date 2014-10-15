/*************************************************************************
Cephes Math Library Release 2.8:  June, 2000
Copyright by Stephen L. Moshier

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from C to
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

#ifndef _trigintegrals_h
#define _trigintegrals_h

#include "ap.h"
#include "ialglib.h"

/*************************************************************************
Sine and cosine integrals

Evaluates the integrals

                         x
                         -
                        |  cos t - 1
  Ci(x) = eul + ln x +  |  --------- dt,
                        |      t
                       -
                        0
            x
            -
           |  sin t
  Si(x) =  |  ----- dt
           |    t
          -
           0

where eul = 0.57721566490153286061 is Euler's constant.
The integrals are approximated by rational functions.
For x > 8 auxiliary functions f(x) and g(x) are employed
such that

Ci(x) = f(x) sin(x) - g(x) cos(x)
Si(x) = pi/2 - f(x) cos(x) - g(x) sin(x)


ACCURACY:
   Test interval = [0,50].
Absolute error, except relative when > 1:
arithmetic   function   # trials      peak         rms
   IEEE        Si        30000       4.4e-16     7.3e-17
   IEEE        Ci        30000       6.9e-16     5.1e-17

Cephes Math Library Release 2.1:  January, 1989
Copyright 1984, 1987, 1989 by Stephen L. Moshier
*************************************************************************/
void sinecosineintegrals(double x, double& si, double& ci);


/*************************************************************************
Hyperbolic sine and cosine integrals

Approximates the integrals

                           x
                           -
                          | |   cosh t - 1
  Chi(x) = eul + ln x +   |    -----------  dt,
                        | |          t
                         -
                         0

              x
              -
             | |  sinh t
  Shi(x) =   |    ------  dt
           | |       t
            -
            0

where eul = 0.57721566490153286061 is Euler's constant.
The integrals are evaluated by power series for x < 8
and by Chebyshev expansions for x between 8 and 88.
For large x, both functions approach exp(x)/2x.
Arguments greater than 88 in magnitude return MAXNUM.


ACCURACY:

Test interval 0 to 88.
                     Relative error:
arithmetic   function  # trials      peak         rms
   IEEE         Shi      30000       6.9e-16     1.6e-16
       Absolute error, except relative when |Chi| > 1:
   IEEE         Chi      30000       8.4e-16     1.4e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
void hyperbolicsinecosineintegrals(double x, double& shi, double& chi);


#endif

