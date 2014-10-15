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
#include "fht.h"

/*************************************************************************
1-dimensional Fast Hartley Transform.

Algorithm has O(N*logN) complexity for any N (composite or prime).

INPUT PARAMETERS
    A   -   array[0..N-1] - real function to be transformed
    N   -   problem size
    
OUTPUT PARAMETERS
    A   -   FHT of a input array, array[0..N-1],
            A_out[k] = sum(A_in[j]*(cos(2*pi*j*k/N)+sin(2*pi*j*k/N)), j=0..N-1)


  -- ALGLIB --
     Copyright 04.06.2009 by Bochkanov Sergey
*************************************************************************/
void fhtr1d(ap::real_1d_array& a, int n)
{
    ftplan plan;
    int i;
    ap::complex_1d_array fa;

    ap::ap_error::make_assertion(n>0, "FHTR1D: incorrect N!");
    
    //
    // Special case: N=1, FHT is just identity transform.
    // After this block we assume that N is strictly greater than 1.
    //
    if( n==1 )
    {
        return;
    }
    
    //
    // Reduce FHt to real FFT
    //
    fftr1d(a, n, fa);
    for(i = 0; i <= n-1; i++)
    {
        a(i) = fa(i).x-fa(i).y;
    }
}


/*************************************************************************
1-dimensional inverse FHT.

Algorithm has O(N*logN) complexity for any N (composite or prime).

INPUT PARAMETERS
    A   -   array[0..N-1] - complex array to be transformed
    N   -   problem size

OUTPUT PARAMETERS
    A   -   inverse FHT of a input array, array[0..N-1]


  -- ALGLIB --
     Copyright 29.05.2009 by Bochkanov Sergey
*************************************************************************/
void fhtr1dinv(ap::real_1d_array& a, int n)
{
    int i;

    ap::ap_error::make_assertion(n>0, "FHTR1DInv: incorrect N!");
    
    //
    // Special case: N=1, iFHT is just identity transform.
    // After this block we assume that N is strictly greater than 1.
    //
    if( n==1 )
    {
        return;
    }
    
    //
    // Inverse FHT can be expressed in terms of the FHT as
    //
    //     invfht(x) = fht(x)/N
    //
    fhtr1d(a, n);
    for(i = 0; i <= n-1; i++)
    {
        a(i) = a(i)/n;
    }
}




