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
#include "kmeans.h"

static bool selectcenterpp(const ap::real_2d_array& xy,
     int npoints,
     int nvars,
     ap::real_2d_array& centers,
     ap::boolean_1d_array busycenters,
     int ccnt,
     ap::real_1d_array& d2,
     ap::real_1d_array& p,
     ap::real_1d_array& tmp);

/*************************************************************************
k-means++ clusterization

INPUT PARAMETERS:
    XY          -   dataset, array [0..NPoints-1,0..NVars-1].
    NPoints     -   dataset size, NPoints>=K
    NVars       -   number of variables, NVars>=1
    K           -   desired number of clusters, K>=1
    Restarts    -   number of restarts, Restarts>=1

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -3, if task is degenerate (number of distinct points is
                          less than K)
                    * -1, if incorrect NPoints/NFeatures/K/Restarts was passed
                    *  1, if subroutine finished successfully
    C           -   array[0..NVars-1,0..K-1].matrix whose columns store
                    cluster's centers
    XYC         -   array which contains number of clusters dataset points
                    belong to.

  -- ALGLIB --
     Copyright 21.03.2009 by Bochkanov Sergey
*************************************************************************/
void kmeansgenerate(const ap::real_2d_array& xy,
     int npoints,
     int nvars,
     int k,
     int restarts,
     int& info,
     ap::real_2d_array& c,
     ap::integer_1d_array& xyc)
{
    int i;
    int j;
    ap::real_2d_array ct;
    ap::real_2d_array ctbest;
    ap::integer_1d_array xycbest;
    double e;
    double ebest;
    ap::real_1d_array x;
    ap::real_1d_array tmp;
    ap::real_1d_array d2;
    ap::real_1d_array p;
    ap::integer_1d_array csizes;
    ap::boolean_1d_array cbusy;
    double v;
    int cclosest;
    double dclosest;
    ap::real_1d_array work;
    bool waschanges;
    bool zerosizeclusters;
    int pass;

    
    //
    // Test parameters
    //
    if( npoints<k||nvars<1||k<1||restarts<1 )
    {
        info = -1;
        return;
    }
    
    //
    // TODO: special case K=1
    // TODO: special case K=NPoints
    //
    info = 1;
    
    //
    // Multiple passes of k-means++ algorithm
    //
    ct.setlength(k, nvars);
    ctbest.setlength(k, nvars);
    xyc.setlength(npoints);
    xycbest.setlength(npoints);
    d2.setlength(npoints);
    p.setlength(npoints);
    tmp.setlength(nvars);
    csizes.setlength(k);
    cbusy.setlength(k);
    ebest = ap::maxrealnumber;
    for(pass = 1; pass <= restarts; pass++)
    {
        
        //
        // Select initial centers  using k-means++ algorithm
        // 1. Choose first center at random
        // 2. Choose next centers using their distance from centers already chosen
        //
        // Note that for performance reasons centers are stored in ROWS of CT, not
        // in columns. We'll transpose CT in the end and store it in the C.
        //
        i = ap::randominteger(npoints);
        ap::vmove(&ct(0, 0), 1, &xy(i, 0), 1, ap::vlen(0,nvars-1));
        cbusy(0) = true;
        for(i = 1; i <= k-1; i++)
        {
            cbusy(i) = false;
        }
        if( !selectcenterpp(xy, npoints, nvars, ct, cbusy, k, d2, p, tmp) )
        {
            info = -3;
            return;
        }
        
        //
        // Update centers:
        // 2. update center positions
        //
        while(true)
        {
            
            //
            // fill XYC with center numbers
            //
            waschanges = false;
            for(i = 0; i <= npoints-1; i++)
            {
                cclosest = -1;
                dclosest = ap::maxrealnumber;
                for(j = 0; j <= k-1; j++)
                {
                    ap::vmove(&tmp(0), 1, &xy(i, 0), 1, ap::vlen(0,nvars-1));
                    ap::vsub(&tmp(0), 1, &ct(j, 0), 1, ap::vlen(0,nvars-1));
                    v = ap::vdotproduct(&tmp(0), 1, &tmp(0), 1, ap::vlen(0,nvars-1));
                    if( ap::fp_less(v,dclosest) )
                    {
                        cclosest = j;
                        dclosest = v;
                    }
                }
                if( xyc(i)!=cclosest )
                {
                    waschanges = true;
                }
                xyc(i) = cclosest;
            }
            
            //
            // Update centers
            //
            for(j = 0; j <= k-1; j++)
            {
                csizes(j) = 0;
            }
            for(i = 0; i <= k-1; i++)
            {
                for(j = 0; j <= nvars-1; j++)
                {
                    ct(i,j) = 0;
                }
            }
            for(i = 0; i <= npoints-1; i++)
            {
                csizes(xyc(i)) = csizes(xyc(i))+1;
                ap::vadd(&ct(xyc(i), 0), 1, &xy(i, 0), 1, ap::vlen(0,nvars-1));
            }
            zerosizeclusters = false;
            for(i = 0; i <= k-1; i++)
            {
                cbusy(i) = csizes(i)!=0;
                zerosizeclusters = zerosizeclusters||csizes(i)==0;
            }
            if( zerosizeclusters )
            {
                
                //
                // Some clusters have zero size - rare, but possible.
                // We'll choose new centers for such clusters using k-means++ rule
                // and restart algorithm
                //
                if( !selectcenterpp(xy, npoints, nvars, ct, cbusy, k, d2, p, tmp) )
                {
                    info = -3;
                    return;
                }
                continue;
            }
            for(j = 0; j <= k-1; j++)
            {
                v = double(1)/double(csizes(j));
                ap::vmul(&ct(j, 0), 1, ap::vlen(0,nvars-1), v);
            }
            
            //
            // if nothing has changed during iteration
            //
            if( !waschanges )
            {
                break;
            }
        }
        
        //
        // 3. Calculate E, compare with best centers found so far
        //
        e = 0;
        for(i = 0; i <= npoints-1; i++)
        {
            ap::vmove(&tmp(0), 1, &xy(i, 0), 1, ap::vlen(0,nvars-1));
            ap::vsub(&tmp(0), 1, &ct(xyc(i), 0), 1, ap::vlen(0,nvars-1));
            v = ap::vdotproduct(&tmp(0), 1, &tmp(0), 1, ap::vlen(0,nvars-1));
            e = e+v;
        }
        if( ap::fp_less(e,ebest) )
        {
            
            //
            // store partition.
            //
            ebest = e;
            copymatrix(ct, 0, k-1, 0, nvars-1, ctbest, 0, k-1, 0, nvars-1);
            for(i = 0; i <= npoints-1; i++)
            {
                xycbest(i) = xyc(i);
            }
        }
    }
    
    //
    // Copy and transpose
    //
    c.setbounds(0, nvars-1, 0, k-1);
    copyandtranspose(ctbest, 0, k-1, 0, nvars-1, c, 0, nvars-1, 0, k-1);
    for(i = 0; i <= npoints-1; i++)
    {
        xyc(i) = xycbest(i);
    }
}


/*************************************************************************
Select center for a new cluster using k-means++ rule
*************************************************************************/
static bool selectcenterpp(const ap::real_2d_array& xy,
     int npoints,
     int nvars,
     ap::real_2d_array& centers,
     ap::boolean_1d_array busycenters,
     int ccnt,
     ap::real_1d_array& d2,
     ap::real_1d_array& p,
     ap::real_1d_array& tmp)
{
    bool result;
    int i;
    int j;
    int cc;
    double v;
    double s;

    result = true;
    for(cc = 0; cc <= ccnt-1; cc++)
    {
        if( !busycenters(cc) )
        {
            
            //
            // fill D2
            //
            for(i = 0; i <= npoints-1; i++)
            {
                d2(i) = ap::maxrealnumber;
                for(j = 0; j <= ccnt-1; j++)
                {
                    if( busycenters(j) )
                    {
                        ap::vmove(&tmp(0), 1, &xy(i, 0), 1, ap::vlen(0,nvars-1));
                        ap::vsub(&tmp(0), 1, &centers(j, 0), 1, ap::vlen(0,nvars-1));
                        v = ap::vdotproduct(&tmp(0), 1, &tmp(0), 1, ap::vlen(0,nvars-1));
                        if( ap::fp_less(v,d2(i)) )
                        {
                            d2(i) = v;
                        }
                    }
                }
            }
            
            //
            // calculate P (non-cumulative)
            //
            s = 0;
            for(i = 0; i <= npoints-1; i++)
            {
                s = s+d2(i);
            }
            if( ap::fp_eq(s,0) )
            {
                result = false;
                return result;
            }
            s = 1/s;
            ap::vmove(&p(0), 1, &d2(0), 1, ap::vlen(0,npoints-1), s);
            
            //
            // choose one of points with probability P
            // random number within (0,1) is generated and
            // inverse empirical CDF is used to randomly choose a point.
            //
            s = 0;
            v = ap::randomreal();
            for(i = 0; i <= npoints-1; i++)
            {
                s = s+p(i);
                if( ap::fp_less_eq(v,s)||i==npoints-1 )
                {
                    ap::vmove(&centers(cc, 0), 1, &xy(i, 0), 1, ap::vlen(0,nvars-1));
                    busycenters(cc) = true;
                    break;
                }
            }
        }
    }
    return result;
}




