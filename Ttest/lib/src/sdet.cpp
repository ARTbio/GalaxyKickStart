/*************************************************************************
Copyright (c) 2005-2007, Sergey Bochkanov (ALGLIB project).

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
#include "sdet.h"

/*************************************************************************
Determinant calculation of the matrix given by LDLT decomposition.

Input parameters:
    A       -   LDLT-decomposition of the matrix,
                output of subroutine SMatrixLDLT.
    Pivots  -   table of permutations which were made during
                LDLT decomposition, output of subroutine SMatrixLDLT.
    N       -   size of matrix A.
    IsUpper -   matrix storage format. The value is equal to the input
                parameter of subroutine SMatrixLDLT.

Result:
    matrix determinant.

  -- ALGLIB --
     Copyright 2005-2008 by Bochkanov Sergey
*************************************************************************/
double smatrixldltdet(const ap::real_2d_array& a,
     const ap::integer_1d_array& pivots,
     int n,
     bool isupper)
{
    double result;
    int k;

    result = 1;
    if( isupper )
    {
        k = 0;
        while(k<n)
        {
            if( pivots(k)>=0 )
            {
                result = result*a(k,k);
                k = k+1;
            }
            else
            {
                result = result*(a(k,k)*a(k+1,k+1)-a(k,k+1)*a(k,k+1));
                k = k+2;
            }
        }
    }
    else
    {
        k = n-1;
        while(k>=0)
        {
            if( pivots(k)>=0 )
            {
                result = result*a(k,k);
                k = k-1;
            }
            else
            {
                result = result*(a(k-1,k-1)*a(k,k)-a(k,k-1)*a(k,k-1));
                k = k-2;
            }
        }
    }
    return result;
}


/*************************************************************************
Determinant calculation of the symmetric matrix

Input parameters:
    A       -   matrix. Array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   if IsUpper = True, then symmetric matrix A is given by its
                upper triangle, and the lower triangle isn’t used by
                subroutine. Similarly, if IsUpper = False, then A is given
                by its lower triangle.

Result:
    determinant of matrix A.

  -- ALGLIB --
     Copyright 2005-2008 by Bochkanov Sergey
*************************************************************************/
double smatrixdet(ap::real_2d_array a, int n, bool isupper)
{
    double result;
    ap::integer_1d_array pivots;

    smatrixldlt(a, n, isupper, pivots);
    result = smatrixldltdet(a, pivots, n, isupper);
    return result;
}


double determinantldlt(const ap::real_2d_array& a,
     const ap::integer_1d_array& pivots,
     int n,
     bool isupper)
{
    double result;
    int k;

    result = 1;
    if( isupper )
    {
        k = 1;
        while(k<=n)
        {
            if( pivots(k)>0 )
            {
                result = result*a(k,k);
                k = k+1;
            }
            else
            {
                result = result*(a(k,k)*a(k+1,k+1)-a(k,k+1)*a(k,k+1));
                k = k+2;
            }
        }
    }
    else
    {
        k = n;
        while(k>=1)
        {
            if( pivots(k)>0 )
            {
                result = result*a(k,k);
                k = k-1;
            }
            else
            {
                result = result*(a(k-1,k-1)*a(k,k)-a(k,k-1)*a(k,k-1));
                k = k-2;
            }
        }
    }
    return result;
}


double determinantsymmetric(ap::real_2d_array a, int n, bool isupper)
{
    double result;
    ap::integer_1d_array pivots;

    ldltdecomposition(a, n, isupper, pivots);
    result = determinantldlt(a, pivots, n, isupper);
    return result;
}




