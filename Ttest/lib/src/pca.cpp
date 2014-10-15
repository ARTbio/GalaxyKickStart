/*************************************************************************
Copyright (c) 2008, Sergey Bochkanov (ALGLIB project).

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
#include "pca.h"

/*************************************************************************
Principal components analysis

Subroutine  builds  orthogonal  basis  where  first  axis  corresponds  to
direction with maximum variance, second axis maximizes variance in subspace
orthogonal to first axis and so on.

It should be noted that, unlike LDA, PCA does not use class labels.

INPUT PARAMETERS:
    X           -   dataset, array[0..NPoints-1,0..NVars-1].
                    matrix contains ONLY INDEPENDENT VARIABLES.
    NPoints     -   dataset size, NPoints>=0
    NVars       -   number of independent variables, NVars>=1

бшундмше оюпюлерпш:
    Info        -   return code:
                    * -4, if SVD subroutine haven't converged
                    * -1, if wrong parameters has been passed (NPoints<0,
                          NVars<1)
                    *  1, if task is solved
    S2          -   array[0..NVars-1]. variance values corresponding
                    to basis vectors.
    V           -   array[0..NVars-1,0..NVars-1]
                    matrix, whose columns store basis vectors.

  -- ALGLIB --
     Copyright 25.08.2008 by Bochkanov Sergey
*************************************************************************/
void pcabuildbasis(const ap::real_2d_array& x,
     int npoints,
     int nvars,
     int& info,
     ap::real_1d_array& s2,
     ap::real_2d_array& v)
{
    ap::real_2d_array a;
    ap::real_2d_array u;
    ap::real_2d_array vt;
    ap::real_1d_array m;
    ap::real_1d_array t;
    int i;
    int j;
    double mean;
    double variance;
    double skewness;
    double kurtosis;

    
    //
    // Check input data
    //
    if( npoints<0||nvars<1 )
    {
        info = -1;
        return;
    }
    info = 1;
    
    //
    // Special case: NPoints=0
    //
    if( npoints==0 )
    {
        s2.setbounds(0, nvars-1);
        v.setbounds(0, nvars-1, 0, nvars-1);
        for(i = 0; i <= nvars-1; i++)
        {
            s2(i) = 0;
        }
        for(i = 0; i <= nvars-1; i++)
        {
            for(j = 0; j <= nvars-1; j++)
            {
                if( i==j )
                {
                    v(i,j) = 1;
                }
                else
                {
                    v(i,j) = 0;
                }
            }
        }
        return;
    }
    
    //
    // Calculate means
    //
    m.setbounds(0, nvars-1);
    t.setbounds(0, npoints-1);
    for(j = 0; j <= nvars-1; j++)
    {
        ap::vmove(&t(0), 1, &x(0, j), x.getstride(), ap::vlen(0,npoints-1));
        calculatemoments(t, npoints, mean, variance, skewness, kurtosis);
        m(j) = mean;
    }
    
    //
    // Center, apply SVD, prepare output
    //
    a.setbounds(0, ap::maxint(npoints, nvars)-1, 0, nvars-1);
    for(i = 0; i <= npoints-1; i++)
    {
        ap::vmove(&a(i, 0), 1, &x(i, 0), 1, ap::vlen(0,nvars-1));
        ap::vsub(&a(i, 0), 1, &m(0), 1, ap::vlen(0,nvars-1));
    }
    for(i = npoints; i <= nvars-1; i++)
    {
        for(j = 0; j <= nvars-1; j++)
        {
            a(i,j) = 0;
        }
    }
    if( !rmatrixsvd(a, ap::maxint(npoints, nvars), nvars, 0, 1, 2, s2, u, vt) )
    {
        info = -4;
        return;
    }
    if( npoints!=1 )
    {
        for(i = 0; i <= nvars-1; i++)
        {
            s2(i) = ap::sqr(s2(i))/(npoints-1);
        }
    }
    v.setbounds(0, nvars-1, 0, nvars-1);
    copyandtranspose(vt, 0, nvars-1, 0, nvars-1, v, 0, nvars-1, 0, nvars-1);
}




