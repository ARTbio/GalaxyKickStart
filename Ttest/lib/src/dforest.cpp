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
#include "dforest.h"

static const int dfvnum = 8;
static const int innernodewidth = 3;
static const int leafnodewidth = 2;
static const int dfusestrongsplits = 1;
static const int dfuseevs = 2;

static int dfclserror(const decisionforest& df,
     const ap::real_2d_array& xy,
     int npoints);
static void dfprocessinternal(const decisionforest& df,
     int offs,
     const ap::real_1d_array& x,
     ap::real_1d_array& y);
static void dfbuildtree(const ap::real_2d_array& xy,
     int npoints,
     int nvars,
     int nclasses,
     int nfeatures,
     int nvarsinpool,
     int flags,
     dfinternalbuffers& bufs);
static void dfbuildtreerec(const ap::real_2d_array& xy,
     int npoints,
     int nvars,
     int nclasses,
     int nfeatures,
     int nvarsinpool,
     int flags,
     int& numprocessed,
     int idx1,
     int idx2,
     dfinternalbuffers& bufs);
static void dfweakspliti(ap::real_1d_array& x,
     ap::integer_1d_array& y,
     int n,
     int nclasses,
     int& info,
     double& threshold,
     double& e);
static void dfsplitc(ap::real_1d_array& x,
     ap::integer_1d_array& c,
     ap::integer_1d_array& cntbuf,
     int n,
     int nc,
     int flags,
     int& info,
     double& threshold,
     double& e);
static void dfsplitr(ap::real_1d_array& x,
     ap::real_1d_array& y,
     int n,
     int flags,
     int& info,
     double& threshold,
     double& e);

/*************************************************************************
This subroutine builds random decision forest.

INPUT PARAMETERS:
    XY          -   training set
    NPoints     -   training set size, NPoints>=1
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   task type:
                    * NClasses=1 - regression task with one
                                   dependent variable
                    * NClasses>1 - classification task with
                                   NClasses classes.
    NTrees      -   number of trees in a forest, NTrees>=1.
                    recommended values: 50-100.
    R           -   percent of a training set used to build
                    individual trees. 0<R<=1.
                    recommended values: 0.1 <= R <= 0.66.

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<1, NVars<1, NClasses<1, NTrees<1, R<=0
                          or R>1).
                    *  1, if task has been solved
    DF          -   model built
    Rep         -   training report, contains error on a training set
                    and out-of-bag estimates of generalization error.

  -- ALGLIB --
     Copyright 19.02.2009 by Bochkanov Sergey
*************************************************************************/
void dfbuildrandomdecisionforest(const ap::real_2d_array& xy,
     int npoints,
     int nvars,
     int nclasses,
     int ntrees,
     double r,
     int& info,
     decisionforest& df,
     dfreport& rep)
{
    int samplesize;

    if( ap::fp_less_eq(r,0)||ap::fp_greater(r,1) )
    {
        info = -1;
        return;
    }
    samplesize = ap::maxint(ap::round(r*npoints), 1);
    dfbuildinternal(xy, npoints, nvars, nclasses, ntrees, samplesize, ap::maxint(nvars/2, 1), dfusestrongsplits+dfuseevs, info, df, rep);
}


void dfbuildinternal(const ap::real_2d_array& xy,
     int npoints,
     int nvars,
     int nclasses,
     int ntrees,
     int samplesize,
     int nfeatures,
     int flags,
     int& info,
     decisionforest& df,
     dfreport& rep)
{
    int i;
    int j;
    int k;
    int tmpi;
    int lasttreeoffs;
    int offs;
    int ooboffs;
    int treesize;
    int nvarsinpool;
    bool useevs;
    dfinternalbuffers bufs;
    ap::integer_1d_array permbuf;
    ap::real_1d_array oobbuf;
    ap::integer_1d_array oobcntbuf;
    ap::real_2d_array xys;
    ap::real_1d_array x;
    ap::real_1d_array y;
    int oobcnt;
    int oobrelcnt;
    double v;
    double vmin;
    double vmax;
    bool bflag;

    
    //
    // Test for inputs
    //
    if( npoints<1||samplesize<1||samplesize>npoints||nvars<1||nclasses<1||ntrees<1||nfeatures<1 )
    {
        info = -1;
        return;
    }
    if( nclasses>1 )
    {
        for(i = 0; i <= npoints-1; i++)
        {
            if( ap::round(xy(i,nvars))<0||ap::round(xy(i,nvars))>=nclasses )
            {
                info = -2;
                return;
            }
        }
    }
    info = 1;
    
    //
    // Flags
    //
    useevs = flags/dfuseevs%2!=0;
    
    //
    // Allocate data, prepare header
    //
    treesize = 1+innernodewidth*(samplesize-1)+leafnodewidth*samplesize;
    permbuf.setbounds(0, npoints-1);
    bufs.treebuf.setbounds(0, treesize-1);
    bufs.idxbuf.setbounds(0, npoints-1);
    bufs.tmpbufr.setbounds(0, npoints-1);
    bufs.tmpbufr2.setbounds(0, npoints-1);
    bufs.tmpbufi.setbounds(0, npoints-1);
    bufs.varpool.setbounds(0, nvars-1);
    bufs.evsbin.setbounds(0, nvars-1);
    bufs.evssplits.setbounds(0, nvars-1);
    bufs.classibuf.setbounds(0, 2*nclasses-1);
    oobbuf.setbounds(0, nclasses*npoints-1);
    oobcntbuf.setbounds(0, npoints-1);
    df.trees.setbounds(0, ntrees*treesize-1);
    xys.setbounds(0, samplesize-1, 0, nvars);
    x.setbounds(0, nvars-1);
    y.setbounds(0, nclasses-1);
    for(i = 0; i <= npoints-1; i++)
    {
        permbuf(i) = i;
    }
    for(i = 0; i <= npoints*nclasses-1; i++)
    {
        oobbuf(i) = 0;
    }
    for(i = 0; i <= npoints-1; i++)
    {
        oobcntbuf(i) = 0;
    }
    
    //
    // Prepare variable pool and EVS (extended variable selection/splitting) buffers
    // (whether EVS is turned on or not):
    // 1. detect binary variables and pre-calculate splits for them
    // 2. detect variables with non-distinct values and exclude them from pool
    //
    for(i = 0; i <= nvars-1; i++)
    {
        bufs.varpool(i) = i;
    }
    nvarsinpool = nvars;
    if( useevs )
    {
        for(j = 0; j <= nvars-1; j++)
        {
            vmin = xy(0,j);
            vmax = vmin;
            for(i = 0; i <= npoints-1; i++)
            {
                v = xy(i,j);
                vmin = ap::minreal(vmin, v);
                vmax = ap::maxreal(vmax, v);
            }
            if( ap::fp_eq(vmin,vmax) )
            {
                
                //
                // exclude variable from pool
                //
                bufs.varpool(j) = bufs.varpool(nvarsinpool-1);
                bufs.varpool(nvarsinpool-1) = -1;
                nvarsinpool = nvarsinpool-1;
                continue;
            }
            bflag = false;
            for(i = 0; i <= npoints-1; i++)
            {
                v = xy(i,j);
                if( ap::fp_neq(v,vmin)&&ap::fp_neq(v,vmax) )
                {
                    bflag = true;
                    break;
                }
            }
            if( bflag )
            {
                
                //
                // non-binary variable
                //
                bufs.evsbin(j) = false;
            }
            else
            {
                
                //
                // Prepare
                //
                bufs.evsbin(j) = true;
                bufs.evssplits(j) = 0.5*(vmin+vmax);
                if( ap::fp_less_eq(bufs.evssplits(j),vmin) )
                {
                    bufs.evssplits(j) = vmax;
                }
            }
        }
    }
    
    //
    // RANDOM FOREST FORMAT
    // W[0]         -   size of array
    // W[1]         -   version number
    // W[2]         -   NVars
    // W[3]         -   NClasses (1 for regression)
    // W[4]         -   NTrees
    // W[5]         -   trees offset
    //
    //
    // TREE FORMAT
    // W[Offs]      -   size of sub-array
    //     node info:
    // W[K+0]       -   variable number        (-1 for leaf mode)
    // W[K+1]       -   threshold              (class/value for leaf node)
    // W[K+2]       -   ">=" branch index      (absent for leaf node)
    //
    //
    df.nvars = nvars;
    df.nclasses = nclasses;
    df.ntrees = ntrees;
    
    //
    // Build forest
    //
    offs = 0;
    for(i = 0; i <= ntrees-1; i++)
    {
        
        //
        // Prepare sample
        //
        for(k = 0; k <= samplesize-1; k++)
        {
            j = k+ap::randominteger(npoints-k);
            tmpi = permbuf(k);
            permbuf(k) = permbuf(j);
            permbuf(j) = tmpi;
            j = permbuf(k);
            ap::vmove(&xys(k, 0), 1, &xy(j, 0), 1, ap::vlen(0,nvars));
        }
        
        //
        // build tree, copy
        //
        dfbuildtree(xys, samplesize, nvars, nclasses, nfeatures, nvarsinpool, flags, bufs);
        j = ap::round(bufs.treebuf(0));
        ap::vmove(&df.trees(offs), 1, &bufs.treebuf(0), 1, ap::vlen(offs,offs+j-1));
        lasttreeoffs = offs;
        offs = offs+j;
        
        //
        // OOB estimates
        //
        for(k = samplesize; k <= npoints-1; k++)
        {
            for(j = 0; j <= nclasses-1; j++)
            {
                y(j) = 0;
            }
            j = permbuf(k);
            ap::vmove(&x(0), 1, &xy(j, 0), 1, ap::vlen(0,nvars-1));
            dfprocessinternal(df, lasttreeoffs, x, y);
            ap::vadd(&oobbuf(j*nclasses), 1, &y(0), 1, ap::vlen(j*nclasses,(j+1)*nclasses-1));
            oobcntbuf(j) = oobcntbuf(j)+1;
        }
    }
    df.bufsize = offs;
    
    //
    // Normalize OOB results
    //
    for(i = 0; i <= npoints-1; i++)
    {
        if( oobcntbuf(i)!=0 )
        {
            v = double(1)/double(oobcntbuf(i));
            ap::vmul(&oobbuf(i*nclasses), 1, ap::vlen(i*nclasses,i*nclasses+nclasses-1), v);
        }
    }
    
    //
    // Calculate training set estimates
    //
    rep.relclserror = dfrelclserror(df, xy, npoints);
    rep.avgce = dfavgce(df, xy, npoints);
    rep.rmserror = dfrmserror(df, xy, npoints);
    rep.avgerror = dfavgerror(df, xy, npoints);
    rep.avgrelerror = dfavgrelerror(df, xy, npoints);
    
    //
    // Calculate OOB estimates.
    //
    rep.oobrelclserror = 0;
    rep.oobavgce = 0;
    rep.oobrmserror = 0;
    rep.oobavgerror = 0;
    rep.oobavgrelerror = 0;
    oobcnt = 0;
    oobrelcnt = 0;
    for(i = 0; i <= npoints-1; i++)
    {
        if( oobcntbuf(i)!=0 )
        {
            ooboffs = i*nclasses;
            if( nclasses>1 )
            {
                
                //
                // classification-specific code
                //
                k = ap::round(xy(i,nvars));
                tmpi = 0;
                for(j = 1; j <= nclasses-1; j++)
                {
                    if( ap::fp_greater(oobbuf(ooboffs+j),oobbuf(ooboffs+tmpi)) )
                    {
                        tmpi = j;
                    }
                }
                if( tmpi!=k )
                {
                    rep.oobrelclserror = rep.oobrelclserror+1;
                }
                if( ap::fp_neq(oobbuf(ooboffs+k),0) )
                {
                    rep.oobavgce = rep.oobavgce-log(oobbuf(ooboffs+k));
                }
                else
                {
                    rep.oobavgce = rep.oobavgce-log(ap::minrealnumber);
                }
                for(j = 0; j <= nclasses-1; j++)
                {
                    if( j==k )
                    {
                        rep.oobrmserror = rep.oobrmserror+ap::sqr(oobbuf(ooboffs+j)-1);
                        rep.oobavgerror = rep.oobavgerror+fabs(oobbuf(ooboffs+j)-1);
                        rep.oobavgrelerror = rep.oobavgrelerror+fabs(oobbuf(ooboffs+j)-1);
                        oobrelcnt = oobrelcnt+1;
                    }
                    else
                    {
                        rep.oobrmserror = rep.oobrmserror+ap::sqr(oobbuf(ooboffs+j));
                        rep.oobavgerror = rep.oobavgerror+fabs(oobbuf(ooboffs+j));
                    }
                }
            }
            else
            {
                
                //
                // regression-specific code
                //
                rep.oobrmserror = rep.oobrmserror+ap::sqr(oobbuf(ooboffs)-xy(i,nvars));
                rep.oobavgerror = rep.oobavgerror+fabs(oobbuf(ooboffs)-xy(i,nvars));
                if( ap::fp_neq(xy(i,nvars),0) )
                {
                    rep.oobavgrelerror = rep.oobavgrelerror+fabs((oobbuf(ooboffs)-xy(i,nvars))/xy(i,nvars));
                    oobrelcnt = oobrelcnt+1;
                }
            }
            
            //
            // update OOB estimates count.
            //
            oobcnt = oobcnt+1;
        }
    }
    if( oobcnt>0 )
    {
        rep.oobrelclserror = rep.oobrelclserror/oobcnt;
        rep.oobavgce = rep.oobavgce/oobcnt;
        rep.oobrmserror = sqrt(rep.oobrmserror/(oobcnt*nclasses));
        rep.oobavgerror = rep.oobavgerror/(oobcnt*nclasses);
        if( oobrelcnt>0 )
        {
            rep.oobavgrelerror = rep.oobavgrelerror/oobrelcnt;
        }
    }
}


/*************************************************************************
Procesing

INPUT PARAMETERS:
    DF      -   decision forest model
    X       -   input vector,  array[0..NVars-1].

OUTPUT PARAMETERS:
    Y       -   result. Regression estimate when solving regression  task,
                vector of posterior probabilities for classification task.
                Subroutine does not allocate memory for this vector, it is
                responsibility of a caller to allocate it. Array  must  be
                at least [0..NClasses-1].

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
void dfprocess(const decisionforest& df,
     const ap::real_1d_array& x,
     ap::real_1d_array& y)
{
    int offs;
    int i;
    double v;

    
    //
    // Proceed
    //
    offs = 0;
    for(i = 0; i <= df.nclasses-1; i++)
    {
        y(i) = 0;
    }
    for(i = 0; i <= df.ntrees-1; i++)
    {
        
        //
        // Process basic tree
        //
        dfprocessinternal(df, offs, x, y);
        
        //
        // Next tree
        //
        offs = offs+ap::round(df.trees(offs));
    }
    v = double(1)/double(df.ntrees);
    ap::vmul(&y(0), 1, ap::vlen(0,df.nclasses-1), v);
}


/*************************************************************************
Relative classification error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    percent of incorrectly classified cases.
    Zero if model solves regression task.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfrelclserror(const decisionforest& df,
     const ap::real_2d_array& xy,
     int npoints)
{
    double result;

    result = double(dfclserror(df, xy, npoints))/double(npoints);
    return result;
}


/*************************************************************************
Average cross-entropy (in bits per element) on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    CrossEntropy/(NPoints*LN(2)).
    Zero if model solves regression task.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfavgce(const decisionforest& df,
     const ap::real_2d_array& xy,
     int npoints)
{
    double result;
    ap::real_1d_array x;
    ap::real_1d_array y;
    int i;
    int j;
    int k;
    int tmpi;

    x.setbounds(0, df.nvars-1);
    y.setbounds(0, df.nclasses-1);
    result = 0;
    for(i = 0; i <= npoints-1; i++)
    {
        ap::vmove(&x(0), 1, &xy(i, 0), 1, ap::vlen(0,df.nvars-1));
        dfprocess(df, x, y);
        if( df.nclasses>1 )
        {
            
            //
            // classification-specific code
            //
            k = ap::round(xy(i,df.nvars));
            tmpi = 0;
            for(j = 1; j <= df.nclasses-1; j++)
            {
                if( ap::fp_greater(y(j),y(tmpi)) )
                {
                    tmpi = j;
                }
            }
            if( ap::fp_neq(y(k),0) )
            {
                result = result-log(y(k));
            }
            else
            {
                result = result-log(ap::minrealnumber);
            }
        }
    }
    result = result/npoints;
    return result;
}


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.
    Its meaning for regression task is obvious. As for
    classification task, RMS error means error when estimating posterior
    probabilities.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfrmserror(const decisionforest& df,
     const ap::real_2d_array& xy,
     int npoints)
{
    double result;
    ap::real_1d_array x;
    ap::real_1d_array y;
    int i;
    int j;
    int k;
    int tmpi;

    x.setbounds(0, df.nvars-1);
    y.setbounds(0, df.nclasses-1);
    result = 0;
    for(i = 0; i <= npoints-1; i++)
    {
        ap::vmove(&x(0), 1, &xy(i, 0), 1, ap::vlen(0,df.nvars-1));
        dfprocess(df, x, y);
        if( df.nclasses>1 )
        {
            
            //
            // classification-specific code
            //
            k = ap::round(xy(i,df.nvars));
            tmpi = 0;
            for(j = 1; j <= df.nclasses-1; j++)
            {
                if( ap::fp_greater(y(j),y(tmpi)) )
                {
                    tmpi = j;
                }
            }
            for(j = 0; j <= df.nclasses-1; j++)
            {
                if( j==k )
                {
                    result = result+ap::sqr(y(j)-1);
                }
                else
                {
                    result = result+ap::sqr(y(j));
                }
            }
        }
        else
        {
            
            //
            // regression-specific code
            //
            result = result+ap::sqr(y(0)-xy(i,df.nvars));
        }
    }
    result = sqrt(result/(npoints*df.nclasses));
    return result;
}


/*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for
    classification task, it means average error when estimating posterior
    probabilities.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfavgerror(const decisionforest& df,
     const ap::real_2d_array& xy,
     int npoints)
{
    double result;
    ap::real_1d_array x;
    ap::real_1d_array y;
    int i;
    int j;
    int k;

    x.setbounds(0, df.nvars-1);
    y.setbounds(0, df.nclasses-1);
    result = 0;
    for(i = 0; i <= npoints-1; i++)
    {
        ap::vmove(&x(0), 1, &xy(i, 0), 1, ap::vlen(0,df.nvars-1));
        dfprocess(df, x, y);
        if( df.nclasses>1 )
        {
            
            //
            // classification-specific code
            //
            k = ap::round(xy(i,df.nvars));
            for(j = 0; j <= df.nclasses-1; j++)
            {
                if( j==k )
                {
                    result = result+fabs(y(j)-1);
                }
                else
                {
                    result = result+fabs(y(j));
                }
            }
        }
        else
        {
            
            //
            // regression-specific code
            //
            result = result+fabs(y(0)-xy(i,df.nvars));
        }
    }
    result = result/(npoints*df.nclasses);
    return result;
}


/*************************************************************************
Average relative error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for
    classification task, it means average relative error when estimating
    posterior probability of belonging to the correct class.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfavgrelerror(const decisionforest& df,
     const ap::real_2d_array& xy,
     int npoints)
{
    double result;
    ap::real_1d_array x;
    ap::real_1d_array y;
    int relcnt;
    int i;
    int j;
    int k;

    x.setbounds(0, df.nvars-1);
    y.setbounds(0, df.nclasses-1);
    result = 0;
    relcnt = 0;
    for(i = 0; i <= npoints-1; i++)
    {
        ap::vmove(&x(0), 1, &xy(i, 0), 1, ap::vlen(0,df.nvars-1));
        dfprocess(df, x, y);
        if( df.nclasses>1 )
        {
            
            //
            // classification-specific code
            //
            k = ap::round(xy(i,df.nvars));
            for(j = 0; j <= df.nclasses-1; j++)
            {
                if( j==k )
                {
                    result = result+fabs(y(j)-1);
                    relcnt = relcnt+1;
                }
            }
        }
        else
        {
            
            //
            // regression-specific code
            //
            if( ap::fp_neq(xy(i,df.nvars),0) )
            {
                result = result+fabs((y(0)-xy(i,df.nvars))/xy(i,df.nvars));
                relcnt = relcnt+1;
            }
        }
    }
    if( relcnt>0 )
    {
        result = result/relcnt;
    }
    return result;
}


/*************************************************************************
Copying of DecisionForest strucure

INPUT PARAMETERS:
    DF1 -   original

OUTPUT PARAMETERS:
    DF2 -   copy

  -- ALGLIB --
     Copyright 13.02.2009 by Bochkanov Sergey
*************************************************************************/
void dfcopy(const decisionforest& df1, decisionforest& df2)
{

    df2.nvars = df1.nvars;
    df2.nclasses = df1.nclasses;
    df2.ntrees = df1.ntrees;
    df2.bufsize = df1.bufsize;
    df2.trees.setbounds(0, df1.bufsize-1);
    ap::vmove(&df2.trees(0), 1, &df1.trees(0), 1, ap::vlen(0,df1.bufsize-1));
}


/*************************************************************************
Serialization of DecisionForest strucure

INPUT PARAMETERS:
    DF      -   original

OUTPUT PARAMETERS:
    RA      -   array of real numbers which stores decision forest,
                array[0..RLen-1]
    RLen    -   RA lenght

  -- ALGLIB --
     Copyright 13.02.2009 by Bochkanov Sergey
*************************************************************************/
void dfserialize(const decisionforest& df, ap::real_1d_array& ra, int& rlen)
{

    ra.setbounds(0, df.bufsize+5-1);
    ra(0) = dfvnum;
    ra(1) = df.nvars;
    ra(2) = df.nclasses;
    ra(3) = df.ntrees;
    ra(4) = df.bufsize;
    ap::vmove(&ra(5), 1, &df.trees(0), 1, ap::vlen(5,5+df.bufsize-1));
    rlen = 5+df.bufsize;
}


/*************************************************************************
Unserialization of DecisionForest strucure

INPUT PARAMETERS:
    RA      -   real array which stores decision forest

OUTPUT PARAMETERS:
    DF      -   restored structure

  -- ALGLIB --
     Copyright 13.02.2009 by Bochkanov Sergey
*************************************************************************/
void dfunserialize(const ap::real_1d_array& ra, decisionforest& df)
{

    ap::ap_error::make_assertion(ap::round(ra(0))==dfvnum, "DFUnserialize: incorrect array!");
    df.nvars = ap::round(ra(1));
    df.nclasses = ap::round(ra(2));
    df.ntrees = ap::round(ra(3));
    df.bufsize = ap::round(ra(4));
    df.trees.setbounds(0, df.bufsize-1);
    ap::vmove(&df.trees(0), 1, &ra(5), 1, ap::vlen(0,df.bufsize-1));
}


/*************************************************************************
Classification error
*************************************************************************/
static int dfclserror(const decisionforest& df,
     const ap::real_2d_array& xy,
     int npoints)
{
    int result;
    ap::real_1d_array x;
    ap::real_1d_array y;
    int i;
    int j;
    int k;
    int tmpi;

    if( df.nclasses<=1 )
    {
        result = 0;
        return result;
    }
    x.setbounds(0, df.nvars-1);
    y.setbounds(0, df.nclasses-1);
    result = 0;
    for(i = 0; i <= npoints-1; i++)
    {
        ap::vmove(&x(0), 1, &xy(i, 0), 1, ap::vlen(0,df.nvars-1));
        dfprocess(df, x, y);
        k = ap::round(xy(i,df.nvars));
        tmpi = 0;
        for(j = 1; j <= df.nclasses-1; j++)
        {
            if( ap::fp_greater(y(j),y(tmpi)) )
            {
                tmpi = j;
            }
        }
        if( tmpi!=k )
        {
            result = result+1;
        }
    }
    return result;
}


/*************************************************************************
Internal subroutine for processing one decision tree starting at Offs
*************************************************************************/
static void dfprocessinternal(const decisionforest& df,
     int offs,
     const ap::real_1d_array& x,
     ap::real_1d_array& y)
{
    int k;
    int idx;

    
    //
    // Set pointer to the root
    //
    k = offs+1;
    
    //
    // Navigate through the tree
    //
    while(true)
    {
        if( ap::fp_eq(df.trees(k),-1) )
        {
            if( df.nclasses==1 )
            {
                y(0) = y(0)+df.trees(k+1);
            }
            else
            {
                idx = ap::round(df.trees(k+1));
                y(idx) = y(idx)+1;
            }
            break;
        }
        if( ap::fp_less(x(ap::round(df.trees(k))),df.trees(k+1)) )
        {
            k = k+innernodewidth;
        }
        else
        {
            k = offs+ap::round(df.trees(k+2));
        }
    }
}


/*************************************************************************
Builds one decision tree. Just a wrapper for the DFBuildTreeRec.
*************************************************************************/
static void dfbuildtree(const ap::real_2d_array& xy,
     int npoints,
     int nvars,
     int nclasses,
     int nfeatures,
     int nvarsinpool,
     int flags,
     dfinternalbuffers& bufs)
{
    int numprocessed;
    int i;

    ap::ap_error::make_assertion(npoints>0, "");
    
    //
    // Prepare IdxBuf. It stores indices of the training set elements.
    // When training set is being split, contents of IdxBuf is
    // correspondingly reordered so we can know which elements belong
    // to which branch of decision tree.
    //
    for(i = 0; i <= npoints-1; i++)
    {
        bufs.idxbuf(i) = i;
    }
    
    //
    // Recursive procedure
    //
    numprocessed = 1;
    dfbuildtreerec(xy, npoints, nvars, nclasses, nfeatures, nvarsinpool, flags, numprocessed, 0, npoints-1, bufs);
    bufs.treebuf(0) = numprocessed;
}


/*************************************************************************
Builds one decision tree (internal recursive subroutine)

Parameters:
    TreeBuf     -   large enough array, at least TreeSize
    IdxBuf      -   at least NPoints elements
    TmpBufR     -   at least NPoints
    TmpBufR2    -   at least NPoints
    TmpBufI     -   at least NPoints
    TmpBufI2    -   at least NPoints+1
*************************************************************************/
static void dfbuildtreerec(const ap::real_2d_array& xy,
     int npoints,
     int nvars,
     int nclasses,
     int nfeatures,
     int nvarsinpool,
     int flags,
     int& numprocessed,
     int idx1,
     int idx2,
     dfinternalbuffers& bufs)
{
    int i;
    int j;
    int k;
    bool bflag;
    int i1;
    int i2;
    int info;
    double sl;
    double sr;
    double w;
    int idxbest;
    double ebest;
    double tbest;
    int varcur;
    double s;
    double v;
    double v1;
    double v2;
    double threshold;
    int oldnp;
    double currms;
    bool useevs;

    ap::ap_error::make_assertion(npoints>0, "");
    ap::ap_error::make_assertion(idx2>=idx1, "");
    useevs = flags/dfuseevs%2!=0;
    
    //
    // Leaf node
    //
    if( idx2==idx1 )
    {
        bufs.treebuf(numprocessed) = -1;
        bufs.treebuf(numprocessed+1) = xy(bufs.idxbuf(idx1),nvars);
        numprocessed = numprocessed+leafnodewidth;
        return;
    }
    
    //
    // Non-leaf node.
    // Select random variable, prepare split:
    // 1. prepare default solution - no splitting, class at random
    // 2. investigate possible splits, compare with default/best
    //
    idxbest = -1;
    if( nclasses>1 )
    {
        
        //
        // default solution for classification
        //
        for(i = 0; i <= nclasses-1; i++)
        {
            bufs.classibuf(i) = 0;
        }
        s = idx2-idx1+1;
        for(i = idx1; i <= idx2; i++)
        {
            j = ap::round(xy(bufs.idxbuf(i),nvars));
            bufs.classibuf(j) = bufs.classibuf(j)+1;
        }
        ebest = 0;
        for(i = 0; i <= nclasses-1; i++)
        {
            ebest = ebest+bufs.classibuf(i)*ap::sqr(1-bufs.classibuf(i)/s)+(s-bufs.classibuf(i))*ap::sqr(bufs.classibuf(i)/s);
        }
        ebest = sqrt(ebest/(nclasses*(idx2-idx1+1)));
    }
    else
    {
        
        //
        // default solution for regression
        //
        v = 0;
        for(i = idx1; i <= idx2; i++)
        {
            v = v+xy(bufs.idxbuf(i),nvars);
        }
        v = v/(idx2-idx1+1);
        ebest = 0;
        for(i = idx1; i <= idx2; i++)
        {
            ebest = ebest+ap::sqr(xy(bufs.idxbuf(i),nvars)-v);
        }
        ebest = sqrt(ebest/(idx2-idx1+1));
    }
    i = 0;
    while(i<=ap::minint(nfeatures, nvarsinpool)-1)
    {
        
        //
        // select variables from pool
        //
        j = i+ap::randominteger(nvarsinpool-i);
        k = bufs.varpool(i);
        bufs.varpool(i) = bufs.varpool(j);
        bufs.varpool(j) = k;
        varcur = bufs.varpool(i);
        
        //
        // load variable values to working array
        //
        // apply EVS preprocessing: if all variable values are same,
        // variable is excluded from pool.
        //
        // This is necessary for binary pre-splits (see later) to work.
        //
        for(j = idx1; j <= idx2; j++)
        {
            bufs.tmpbufr(j-idx1) = xy(bufs.idxbuf(j),varcur);
        }
        if( useevs )
        {
            bflag = false;
            v = bufs.tmpbufr(0);
            for(j = 0; j <= idx2-idx1; j++)
            {
                if( ap::fp_neq(bufs.tmpbufr(j),v) )
                {
                    bflag = true;
                    break;
                }
            }
            if( !bflag )
            {
                
                //
                // exclude variable from pool,
                // go to the next iteration.
                // I is not increased.
                //
                k = bufs.varpool(i);
                bufs.varpool(i) = bufs.varpool(nvarsinpool-1);
                bufs.varpool(nvarsinpool-1) = k;
                nvarsinpool = nvarsinpool-1;
                continue;
            }
        }
        
        //
        // load labels to working array
        //
        if( nclasses>1 )
        {
            for(j = idx1; j <= idx2; j++)
            {
                bufs.tmpbufi(j-idx1) = ap::round(xy(bufs.idxbuf(j),nvars));
            }
        }
        else
        {
            for(j = idx1; j <= idx2; j++)
            {
                bufs.tmpbufr2(j-idx1) = xy(bufs.idxbuf(j),nvars);
            }
        }
        
        //
        // calculate split
        //
        if( useevs&&bufs.evsbin(varcur) )
        {
            
            //
            // Pre-calculated splits for binary variables.
            // Threshold is already known, just calculate RMS error
            //
            threshold = bufs.evssplits(varcur);
            if( nclasses>1 )
            {
                
                //
                // classification-specific code
                //
                for(j = 0; j <= 2*nclasses-1; j++)
                {
                    bufs.classibuf(j) = 0;
                }
                sl = 0;
                sr = 0;
                for(j = 0; j <= idx2-idx1; j++)
                {
                    k = bufs.tmpbufi(j);
                    if( ap::fp_less(bufs.tmpbufr(j),threshold) )
                    {
                        bufs.classibuf(k) = bufs.classibuf(k)+1;
                        sl = sl+1;
                    }
                    else
                    {
                        bufs.classibuf(k+nclasses) = bufs.classibuf(k+nclasses)+1;
                        sr = sr+1;
                    }
                }
                ap::ap_error::make_assertion(ap::fp_neq(sl,0)&&ap::fp_neq(sr,0), "DFBuildTreeRec: something strange!");
                currms = 0;
                for(j = 0; j <= nclasses-1; j++)
                {
                    w = bufs.classibuf(j);
                    currms = currms+w*ap::sqr(w/sl-1);
                    currms = currms+(sl-w)*ap::sqr(w/sl);
                    w = bufs.classibuf(nclasses+j);
                    currms = currms+w*ap::sqr(w/sr-1);
                    currms = currms+(sr-w)*ap::sqr(w/sr);
                }
                currms = sqrt(currms/(nclasses*(idx2-idx1+1)));
            }
            else
            {
                
                //
                // regression-specific code
                //
                sl = 0;
                sr = 0;
                v1 = 0;
                v2 = 0;
                for(j = 0; j <= idx2-idx1; j++)
                {
                    if( ap::fp_less(bufs.tmpbufr(j),threshold) )
                    {
                        v1 = v1+bufs.tmpbufr2(j);
                        sl = sl+1;
                    }
                    else
                    {
                        v2 = v2+bufs.tmpbufr2(j);
                        sr = sr+1;
                    }
                }
                ap::ap_error::make_assertion(ap::fp_neq(sl,0)&&ap::fp_neq(sr,0), "DFBuildTreeRec: something strange!");
                v1 = v1/sl;
                v2 = v2/sr;
                currms = 0;
                for(j = 0; j <= idx2-idx1; j++)
                {
                    if( ap::fp_less(bufs.tmpbufr(j),threshold) )
                    {
                        currms = currms+ap::sqr(v1-bufs.tmpbufr2(j));
                    }
                    else
                    {
                        currms = currms+ap::sqr(v2-bufs.tmpbufr2(j));
                    }
                }
                currms = sqrt(currms/(idx2-idx1+1));
            }
            info = 1;
        }
        else
        {
            
            //
            // Generic splits
            //
            if( nclasses>1 )
            {
                dfsplitc(bufs.tmpbufr, bufs.tmpbufi, bufs.classibuf, idx2-idx1+1, nclasses, dfusestrongsplits, info, threshold, currms);
            }
            else
            {
                dfsplitr(bufs.tmpbufr, bufs.tmpbufr2, idx2-idx1+1, dfusestrongsplits, info, threshold, currms);
            }
        }
        if( info>0 )
        {
            if( ap::fp_less_eq(currms,ebest) )
            {
                ebest = currms;
                idxbest = varcur;
                tbest = threshold;
            }
        }
        
        //
        // Next iteration
        //
        i = i+1;
    }
    
    //
    // to split or not to split
    //
    if( idxbest<0 )
    {
        
        //
        // All values are same, cannot split.
        //
        bufs.treebuf(numprocessed) = -1;
        if( nclasses>1 )
        {
            
            //
            // Select random class label (randomness allows us to
            // approximate distribution of the classes)
            //
            bufs.treebuf(numprocessed+1) = ap::round(xy(bufs.idxbuf(idx1+ap::randominteger(idx2-idx1+1)),nvars));
        }
        else
        {
            
            //
            // Select average (for regression task).
            //
            v = 0;
            for(i = idx1; i <= idx2; i++)
            {
                v = v+xy(bufs.idxbuf(i),nvars)/(idx2-idx1+1);
            }
            bufs.treebuf(numprocessed+1) = v;
        }
        numprocessed = numprocessed+leafnodewidth;
    }
    else
    {
        
        //
        // we can split
        //
        bufs.treebuf(numprocessed) = idxbest;
        bufs.treebuf(numprocessed+1) = tbest;
        i1 = idx1;
        i2 = idx2;
        while(i1<=i2)
        {
            
            //
            // Reorder indices so that left partition is in [Idx1..I1-1],
            // and right partition is in [I2+1..Idx2]
            //
            if( ap::fp_less(xy(bufs.idxbuf(i1),idxbest),tbest) )
            {
                i1 = i1+1;
                continue;
            }
            if( ap::fp_greater_eq(xy(bufs.idxbuf(i2),idxbest),tbest) )
            {
                i2 = i2-1;
                continue;
            }
            j = bufs.idxbuf(i1);
            bufs.idxbuf(i1) = bufs.idxbuf(i2);
            bufs.idxbuf(i2) = j;
            i1 = i1+1;
            i2 = i2-1;
        }
        oldnp = numprocessed;
        numprocessed = numprocessed+innernodewidth;
        dfbuildtreerec(xy, npoints, nvars, nclasses, nfeatures, nvarsinpool, flags, numprocessed, idx1, i1-1, bufs);
        bufs.treebuf(oldnp+2) = numprocessed;
        dfbuildtreerec(xy, npoints, nvars, nclasses, nfeatures, nvarsinpool, flags, numprocessed, i2+1, idx2, bufs);
    }
}


/*************************************************************************
Makes weak split on attribute
*************************************************************************/
static void dfweakspliti(ap::real_1d_array& x,
     ap::integer_1d_array& y,
     int n,
     int nclasses,
     int& info,
     double& threshold,
     double& e)
{
    int i;
    int neq;
    int nless;
    int ngreater;

    tagsortfasti(x, y, n);
    if( n%2==1 )
    {
        
        //
        // odd number of elements
        //
        threshold = x(n/2);
    }
    else
    {
        
        //
        // even number of elements.
        //
        // if two closest to the middle of the array are equal,
        // we will select one of them (to avoid possible problems with
        // floating point errors).
        // we will select halfsum otherwise.
        //
        if( ap::fp_eq(x(n/2-1),x(n/2)) )
        {
            threshold = x(n/2-1);
        }
        else
        {
            threshold = 0.5*(x(n/2-1)+x(n/2));
        }
    }
    neq = 0;
    nless = 0;
    ngreater = 0;
    for(i = 0; i <= n-1; i++)
    {
        if( ap::fp_less(x(i),threshold) )
        {
            nless = nless+1;
        }
        if( ap::fp_eq(x(i),threshold) )
        {
            neq = neq+1;
        }
        if( ap::fp_greater(x(i),threshold) )
        {
            ngreater = ngreater+1;
        }
    }
    if( nless==0&&ngreater==0 )
    {
        info = -3;
    }
    else
    {
        if( neq!=0 )
        {
            if( nless<ngreater )
            {
                threshold = 0.5*(x(nless+neq-1)+x(nless+neq));
            }
            else
            {
                threshold = 0.5*(x(nless-1)+x(nless));
            }
        }
        info = 1;
        e = 0;
    }
}


/*************************************************************************
Makes split on attribute
*************************************************************************/
static void dfsplitc(ap::real_1d_array& x,
     ap::integer_1d_array& c,
     ap::integer_1d_array& cntbuf,
     int n,
     int nc,
     int flags,
     int& info,
     double& threshold,
     double& e)
{
    int i;
    int neq;
    int nless;
    int ngreater;
    int q;
    int qmin;
    int qmax;
    int qcnt;
    double cursplit;
    int nleft;
    double v;
    double cure;
    double w;
    double sl;
    double sr;

    tagsortfasti(x, c, n);
    e = ap::maxrealnumber;
    threshold = 0.5*(x(0)+x(n-1));
    info = -3;
    if( flags/dfusestrongsplits%2==0 )
    {
        
        //
        // weak splits, split at half
        //
        qcnt = 2;
        qmin = 1;
        qmax = 1;
    }
    else
    {
        
        //
        // strong splits: choose best quartile
        //
        qcnt = 4;
        qmin = 1;
        qmax = 3;
    }
    for(q = qmin; q <= qmax; q++)
    {
        cursplit = x(n*q/qcnt);
        neq = 0;
        nless = 0;
        ngreater = 0;
        for(i = 0; i <= n-1; i++)
        {
            if( ap::fp_less(x(i),cursplit) )
            {
                nless = nless+1;
            }
            if( ap::fp_eq(x(i),cursplit) )
            {
                neq = neq+1;
            }
            if( ap::fp_greater(x(i),cursplit) )
            {
                ngreater = ngreater+1;
            }
        }
        ap::ap_error::make_assertion(neq!=0, "DFSplitR: NEq=0, something strange!!!");
        if( nless!=0||ngreater!=0 )
        {
            
            //
            // set threshold between two partitions, with
            // some tweaking to avoid problems with floating point
            // arithmetics.
            //
            // The problem is that when you calculates C = 0.5*(A+B) there
            // can be no C which lies strictly between A and B (for example,
            // there is no floating point number which is
            // greater than 1 and less than 1+eps). In such situations
            // we choose right side as theshold (remember that
            // points which lie on threshold falls to the right side).
            //
            if( nless<ngreater )
            {
                cursplit = 0.5*(x(nless+neq-1)+x(nless+neq));
                nleft = nless+neq;
                if( ap::fp_less_eq(cursplit,x(nless+neq-1)) )
                {
                    cursplit = x(nless+neq);
                }
            }
            else
            {
                cursplit = 0.5*(x(nless-1)+x(nless));
                nleft = nless;
                if( ap::fp_less_eq(cursplit,x(nless-1)) )
                {
                    cursplit = x(nless);
                }
            }
            info = 1;
            cure = 0;
            for(i = 0; i <= 2*nc-1; i++)
            {
                cntbuf(i) = 0;
            }
            for(i = 0; i <= nleft-1; i++)
            {
                cntbuf(c(i)) = cntbuf(c(i))+1;
            }
            for(i = nleft; i <= n-1; i++)
            {
                cntbuf(nc+c(i)) = cntbuf(nc+c(i))+1;
            }
            sl = nleft;
            sr = n-nleft;
            v = 0;
            for(i = 0; i <= nc-1; i++)
            {
                w = cntbuf(i);
                v = v+w*ap::sqr(w/sl-1);
                v = v+(sl-w)*ap::sqr(w/sl);
                w = cntbuf(nc+i);
                v = v+w*ap::sqr(w/sr-1);
                v = v+(sr-w)*ap::sqr(w/sr);
            }
            cure = sqrt(v/(nc*n));
            if( ap::fp_less(cure,e) )
            {
                threshold = cursplit;
                e = cure;
            }
        }
    }
}


/*************************************************************************
Makes split on attribute
*************************************************************************/
static void dfsplitr(ap::real_1d_array& x,
     ap::real_1d_array& y,
     int n,
     int flags,
     int& info,
     double& threshold,
     double& e)
{
    int i;
    int neq;
    int nless;
    int ngreater;
    int q;
    int qmin;
    int qmax;
    int qcnt;
    double cursplit;
    int nleft;
    double v;
    double cure;

    tagsortfastr(x, y, n);
    e = ap::maxrealnumber;
    threshold = 0.5*(x(0)+x(n-1));
    info = -3;
    if( flags/dfusestrongsplits%2==0 )
    {
        
        //
        // weak splits, split at half
        //
        qcnt = 2;
        qmin = 1;
        qmax = 1;
    }
    else
    {
        
        //
        // strong splits: choose best quartile
        //
        qcnt = 4;
        qmin = 1;
        qmax = 3;
    }
    for(q = qmin; q <= qmax; q++)
    {
        cursplit = x(n*q/qcnt);
        neq = 0;
        nless = 0;
        ngreater = 0;
        for(i = 0; i <= n-1; i++)
        {
            if( ap::fp_less(x(i),cursplit) )
            {
                nless = nless+1;
            }
            if( ap::fp_eq(x(i),cursplit) )
            {
                neq = neq+1;
            }
            if( ap::fp_greater(x(i),cursplit) )
            {
                ngreater = ngreater+1;
            }
        }
        ap::ap_error::make_assertion(neq!=0, "DFSplitR: NEq=0, something strange!!!");
        if( nless!=0||ngreater!=0 )
        {
            
            //
            // set threshold between two partitions, with
            // some tweaking to avoid problems with floating point
            // arithmetics.
            //
            // The problem is that when you calculates C = 0.5*(A+B) there
            // can be no C which lies strictly between A and B (for example,
            // there is no floating point number which is
            // greater than 1 and less than 1+eps). In such situations
            // we choose right side as theshold (remember that
            // points which lie on threshold falls to the right side).
            //
            if( nless<ngreater )
            {
                cursplit = 0.5*(x(nless+neq-1)+x(nless+neq));
                nleft = nless+neq;
                if( ap::fp_less_eq(cursplit,x(nless+neq-1)) )
                {
                    cursplit = x(nless+neq);
                }
            }
            else
            {
                cursplit = 0.5*(x(nless-1)+x(nless));
                nleft = nless;
                if( ap::fp_less_eq(cursplit,x(nless-1)) )
                {
                    cursplit = x(nless);
                }
            }
            info = 1;
            cure = 0;
            v = 0;
            for(i = 0; i <= nleft-1; i++)
            {
                v = v+y(i);
            }
            v = v/nleft;
            for(i = 0; i <= nleft-1; i++)
            {
                cure = cure+ap::sqr(y(i)-v);
            }
            v = 0;
            for(i = nleft; i <= n-1; i++)
            {
                v = v+y(i);
            }
            v = v/(n-nleft);
            for(i = nleft; i <= n-1; i++)
            {
                cure = cure+ap::sqr(y(i)-v);
            }
            cure = sqrt(cure/n);
            if( ap::fp_less(cure,e) )
            {
                threshold = cursplit;
                e = cure;
            }
        }
    }
}




