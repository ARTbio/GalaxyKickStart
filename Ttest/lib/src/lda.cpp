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
#include "lda.h"

/*************************************************************************
Multiclass Fisher LDA

Subroutine finds coefficients of linear combination which optimally separates
training set on classes.

INPUT PARAMETERS:
    XY          -   training set, array[0..NPoints-1,0..NVars].
                    First NVars columns store values of independent
                    variables, next column stores number of class (from 0
                    to NClasses-1) which dataset element belongs to. Fractional
                    values are rounded to nearest integer.
    NPoints     -   training set size, NPoints>=0
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   number of classes, NClasses>=2


OUTPUT PARAMETERS:
    Info        -   return code:
                    * -4, if internal EVD subroutine hasn't converged
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed (NPoints<0,
                          NVars<1, NClasses<2)
                    *  1, if task has been solved
                    *  2, if there was a multicollinearity in training set,
                          but task has been solved.
    W           -   linear combination coefficients, array[0..NVars-1]

  -- ALGLIB --
     Copyright 31.05.2008 by Bochkanov Sergey
*************************************************************************/
void fisherlda(const ap::real_2d_array& xy,
     int npoints,
     int nvars,
     int nclasses,
     int& info,
     ap::real_1d_array& w)
{
    ap::real_2d_array w2;

    fisherldan(xy, npoints, nvars, nclasses, info, w2);
    if( info>0 )
    {
        w.setbounds(0, nvars-1);
        ap::vmove(&w(0), 1, &w2(0, 0), w2.getstride(), ap::vlen(0,nvars-1));
    }
}


/*************************************************************************
N-dimensional multiclass Fisher LDA

Subroutine finds coefficients of linear combinations which optimally separates
training set on classes. It returns N-dimensional basis whose vector are sorted
by quality of training set separation (in descending order).

INPUT PARAMETERS:
    XY          -   training set, array[0..NPoints-1,0..NVars].
                    First NVars columns store values of independent
                    variables, next column stores number of class (from 0
                    to NClasses-1) which dataset element belongs to. Fractional
                    values are rounded to nearest integer.
    NPoints     -   training set size, NPoints>=0
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   number of classes, NClasses>=2


OUTPUT PARAMETERS:
    Info        -   return code:
                    * -4, if internal EVD subroutine hasn't converged
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed (NPoints<0,
                          NVars<1, NClasses<2)
                    *  1, if task has been solved
                    *  2, if there was a multicollinearity in training set,
                          but task has been solved.
    W           -   basis, array[0..NVars-1,0..NVars-1]
                    columns of matrix stores basis vectors, sorted by
                    quality of training set separation (in descending order)

  -- ALGLIB --
     Copyright 31.05.2008 by Bochkanov Sergey
*************************************************************************/
void fisherldan(const ap::real_2d_array& xy,
     int npoints,
     int nvars,
     int nclasses,
     int& info,
     ap::real_2d_array& w)
{
    int i;
    int j;
    int k;
    int m;
    double v;
    ap::integer_1d_array c;
    ap::real_1d_array mu;
    ap::real_2d_array muc;
    ap::integer_1d_array nc;
    ap::real_2d_array sw;
    ap::real_2d_array st;
    ap::real_2d_array z;
    ap::real_2d_array z2;
    ap::real_2d_array tm;
    ap::real_2d_array sbroot;
    ap::real_2d_array a;
    ap::real_2d_array xyproj;
    ap::real_2d_array wproj;
    ap::real_1d_array tf;
    ap::real_1d_array d;
    ap::real_1d_array d2;
    ap::real_1d_array work;

    
    //
    // Test data
    //
    if( npoints<0||nvars<1||nclasses<2 )
    {
        info = -1;
        return;
    }
    for(i = 0; i <= npoints-1; i++)
    {
        if( ap::round(xy(i,nvars))<0||ap::round(xy(i,nvars))>=nclasses )
        {
            info = -2;
            return;
        }
    }
    info = 1;
    
    //
    // Special case: NPoints<=1
    // Degenerate task.
    //
    if( npoints<=1 )
    {
        info = 2;
        w.setbounds(0, nvars-1, 0, nvars-1);
        for(i = 0; i <= nvars-1; i++)
        {
            for(j = 0; j <= nvars-1; j++)
            {
                if( i==j )
                {
                    w(i,j) = 1;
                }
                else
                {
                    w(i,j) = 0;
                }
            }
        }
        return;
    }
    
    //
    // Prepare temporaries
    //
    tf.setbounds(0, nvars-1);
    work.setbounds(1, ap::maxint(nvars, npoints));
    
    //
    // Convert class labels from reals to integers (just for convenience)
    //
    c.setbounds(0, npoints-1);
    for(i = 0; i <= npoints-1; i++)
    {
        c(i) = ap::round(xy(i,nvars));
    }
    
    //
    // Calculate class sizes and means
    //
    mu.setbounds(0, nvars-1);
    muc.setbounds(0, nclasses-1, 0, nvars-1);
    nc.setbounds(0, nclasses-1);
    for(j = 0; j <= nvars-1; j++)
    {
        mu(j) = 0;
    }
    for(i = 0; i <= nclasses-1; i++)
    {
        nc(i) = 0;
        for(j = 0; j <= nvars-1; j++)
        {
            muc(i,j) = 0;
        }
    }
    for(i = 0; i <= npoints-1; i++)
    {
        ap::vadd(&mu(0), 1, &xy(i, 0), 1, ap::vlen(0,nvars-1));
        ap::vadd(&muc(c(i), 0), 1, &xy(i, 0), 1, ap::vlen(0,nvars-1));
        nc(c(i)) = nc(c(i))+1;
    }
    for(i = 0; i <= nclasses-1; i++)
    {
        v = double(1)/double(nc(i));
        ap::vmul(&muc(i, 0), 1, ap::vlen(0,nvars-1), v);
    }
    v = double(1)/double(npoints);
    ap::vmul(&mu(0), 1, ap::vlen(0,nvars-1), v);
    
    //
    // Create ST matrix
    //
    st.setbounds(0, nvars-1, 0, nvars-1);
    for(i = 0; i <= nvars-1; i++)
    {
        for(j = 0; j <= nvars-1; j++)
        {
            st(i,j) = 0;
        }
    }
    for(k = 0; k <= npoints-1; k++)
    {
        ap::vmove(&tf(0), 1, &xy(k, 0), 1, ap::vlen(0,nvars-1));
        ap::vsub(&tf(0), 1, &mu(0), 1, ap::vlen(0,nvars-1));
        for(i = 0; i <= nvars-1; i++)
        {
            v = tf(i);
            ap::vadd(&st(i, 0), 1, &tf(0), 1, ap::vlen(0,nvars-1), v);
        }
    }
    
    //
    // Create SW matrix
    //
    sw.setbounds(0, nvars-1, 0, nvars-1);
    for(i = 0; i <= nvars-1; i++)
    {
        for(j = 0; j <= nvars-1; j++)
        {
            sw(i,j) = 0;
        }
    }
    for(k = 0; k <= npoints-1; k++)
    {
        ap::vmove(&tf(0), 1, &xy(k, 0), 1, ap::vlen(0,nvars-1));
        ap::vsub(&tf(0), 1, &muc(c(k), 0), 1, ap::vlen(0,nvars-1));
        for(i = 0; i <= nvars-1; i++)
        {
            v = tf(i);
            ap::vadd(&sw(i, 0), 1, &tf(0), 1, ap::vlen(0,nvars-1), v);
        }
    }
    
    //
    // Maximize ratio J=(w'*ST*w)/(w'*SW*w).
    //
    // First, make transition from w to v such that w'*ST*w becomes v'*v:
    //    v  = root(ST)*w = R*w
    //    R  = root(D)*Z'
    //    w  = (root(ST)^-1)*v = RI*v
    //    RI = Z*inv(root(D))
    //    J  = (v'*v)/(v'*(RI'*SW*RI)*v)
    //    ST = Z*D*Z'
    //
    //    so we have
    //
    //    J = (v'*v) / (v'*(inv(root(D))*Z'*SW*Z*inv(root(D)))*v)  =
    //      = (v'*v) / (v'*A*v)
    //
    if( !smatrixevd(st, nvars, 1, true, d, z) )
    {
        info = -4;
        return;
    }
    w.setbounds(0, nvars-1, 0, nvars-1);
    if( ap::fp_less_eq(d(nvars-1),0)||ap::fp_less_eq(d(0),1000*ap::machineepsilon*d(nvars-1)) )
    {
        
        //
        // Special case: D[NVars-1]<=0
        // Degenerate task (all variables takes the same value).
        //
        if( ap::fp_less_eq(d(nvars-1),0) )
        {
            info = 2;
            for(i = 0; i <= nvars-1; i++)
            {
                for(j = 0; j <= nvars-1; j++)
                {
                    if( i==j )
                    {
                        w(i,j) = 1;
                    }
                    else
                    {
                        w(i,j) = 0;
                    }
                }
            }
            return;
        }
        
        //
        // Special case: degenerate ST matrix, multicollinearity found.
        // Since we know ST eigenvalues/vectors we can translate task to
        // non-degenerate form.
        //
        // Let WG is orthogonal basis of the non zero variance subspace
        // of the ST and let WZ is orthogonal basis of the zero variance
        // subspace.
        //
        // Projection on WG allows us to use LDA on reduced M-dimensional
        // subspace, N-M vectors of WZ allows us to update reduced LDA
        // factors to full N-dimensional subspace.
        //
        m = 0;
        for(k = 0; k <= nvars-1; k++)
        {
            if( ap::fp_less_eq(d(k),1000*ap::machineepsilon*d(nvars-1)) )
            {
                m = k+1;
            }
        }
        ap::ap_error::make_assertion(m!=0, "FisherLDAN: internal error #1");
        xyproj.setbounds(0, npoints-1, 0, nvars-m);
        matrixmatrixmultiply(xy, 0, npoints-1, 0, nvars-1, false, z, 0, nvars-1, m, nvars-1, false, 1.0, xyproj, 0, npoints-1, 0, nvars-m-1, 0.0, work);
        for(i = 0; i <= npoints-1; i++)
        {
            xyproj(i,nvars-m) = xy(i,nvars);
        }
        fisherldan(xyproj, npoints, nvars-m, nclasses, info, wproj);
        if( info<0 )
        {
            return;
        }
        matrixmatrixmultiply(z, 0, nvars-1, m, nvars-1, false, wproj, 0, nvars-m-1, 0, nvars-m-1, false, 1.0, w, 0, nvars-1, 0, nvars-m-1, 0.0, work);
        for(k = nvars-m; k <= nvars-1; k++)
        {
            ap::vmove(&w(0, k), w.getstride(), &z(0, k-(nvars-m)), z.getstride(), ap::vlen(0,nvars-1));
        }
        info = 2;
    }
    else
    {
        
        //
        // General case: no multicollinearity
        //
        tm.setbounds(0, nvars-1, 0, nvars-1);
        a.setbounds(0, nvars-1, 0, nvars-1);
        matrixmatrixmultiply(sw, 0, nvars-1, 0, nvars-1, false, z, 0, nvars-1, 0, nvars-1, false, 1.0, tm, 0, nvars-1, 0, nvars-1, 0.0, work);
        matrixmatrixmultiply(z, 0, nvars-1, 0, nvars-1, true, tm, 0, nvars-1, 0, nvars-1, false, 1.0, a, 0, nvars-1, 0, nvars-1, 0.0, work);
        for(i = 0; i <= nvars-1; i++)
        {
            for(j = 0; j <= nvars-1; j++)
            {
                a(i,j) = a(i,j)/sqrt(d(i)*d(j));
            }
        }
        if( !smatrixevd(a, nvars, 1, true, d2, z2) )
        {
            info = -4;
            return;
        }
        for(k = 0; k <= nvars-1; k++)
        {
            for(i = 0; i <= nvars-1; i++)
            {
                tf(i) = z2(i,k)/sqrt(d(i));
            }
            for(i = 0; i <= nvars-1; i++)
            {
                v = ap::vdotproduct(&z(i, 0), 1, &tf(0), 1, ap::vlen(0,nvars-1));
                w(i,k) = v;
            }
        }
    }
    
    //
    // Post-processing:
    // * normalization
    // * converting to non-negative form, if possible
    //
    for(k = 0; k <= nvars-1; k++)
    {
        v = ap::vdotproduct(&w(0, k), w.getstride(), &w(0, k), w.getstride(), ap::vlen(0,nvars-1));
        v = 1/sqrt(v);
        ap::vmul(&w(0, k), w.getstride(), ap::vlen(0,nvars-1), v);
        v = 0;
        for(i = 0; i <= nvars-1; i++)
        {
            v = v+w(i,k);
        }
        if( ap::fp_less(v,0) )
        {
            ap::vmul(&w(0, k), w.getstride(), ap::vlen(0,nvars-1), -1);
        }
    }
}




