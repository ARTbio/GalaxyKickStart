/*************************************************************************
Copyright 2008 by Sergey Bochkanov (ALGLIB project).

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
#include "bdss.h"

static void dskfoldsplit(const ap::real_2d_array& xy,
     int npoints,
     int nclasses,
     int foldscount,
     bool stratifiedsplits,
     ap::integer_1d_array& folds);
static double xlny(double x, double y);
static double getcv(const ap::integer_1d_array& cnt, int nc);
static void tieaddc(const ap::integer_1d_array& c,
     const ap::integer_1d_array& ties,
     int ntie,
     int nc,
     ap::integer_1d_array& cnt);
static void tiesubc(const ap::integer_1d_array& c,
     const ap::integer_1d_array& ties,
     int ntie,
     int nc,
     ap::integer_1d_array& cnt);
static void tiegetc(const ap::integer_1d_array& c,
     const ap::integer_1d_array& ties,
     int ntie,
     int nc,
     ap::integer_1d_array& cnt);

/*************************************************************************
This set of routines (DSErrAllocate, DSErrAccumulate, DSErrFinish)
calculates different error functions (classification error, cross-entropy,
rms, avg, avg.rel errors).

1. DSErrAllocate prepares buffer.
2. DSErrAccumulate accumulates individual errors:
    * Y contains predicted output (posterior probabilities for classification)
    * DesiredY contains desired output (class number for classification)
3. DSErrFinish outputs results:
   * Buf[0] contains relative classification error (zero for regression tasks)
   * Buf[1] contains avg. cross-entropy (zero for regression tasks)
   * Buf[2] contains rms error (regression, classification)
   * Buf[3] contains average error (regression, classification)
   * Buf[4] contains average relative error (regression, classification)
   
NOTES(1):
    "NClasses>0" means that we have classification task.
    "NClasses<0" means regression task with -NClasses real outputs.

NOTES(2):
    rms. avg, avg.rel errors for classification tasks are interpreted as
    errors in posterior probabilities with respect to probabilities given
    by training/test set.

  -- ALGLIB --
     Copyright 11.01.2009 by Bochkanov Sergey
*************************************************************************/
void dserrallocate(int nclasses, ap::real_1d_array& buf)
{

    buf.setbounds(0, 7);
    buf(0) = 0;
    buf(1) = 0;
    buf(2) = 0;
    buf(3) = 0;
    buf(4) = 0;
    buf(5) = nclasses;
    buf(6) = 0;
    buf(7) = 0;
}


/*************************************************************************
See DSErrAllocate for comments on this routine.

  -- ALGLIB --
     Copyright 11.01.2009 by Bochkanov Sergey
*************************************************************************/
void dserraccumulate(ap::real_1d_array& buf,
     const ap::real_1d_array& y,
     const ap::real_1d_array& desiredy)
{
    int nclasses;
    int nout;
    int offs;
    int mmax;
    int rmax;
    int j;
    double v;
    double ev;

    offs = 5;
    nclasses = ap::round(buf(offs));
    if( nclasses>0 )
    {
        
        //
        // Classification
        //
        rmax = ap::round(desiredy(0));
        mmax = 0;
        for(j = 1; j <= nclasses-1; j++)
        {
            if( ap::fp_greater(y(j),y(mmax)) )
            {
                mmax = j;
            }
        }
        if( mmax!=rmax )
        {
            buf(0) = buf(0)+1;
        }
        if( ap::fp_greater(y(rmax),0) )
        {
            buf(1) = buf(1)-log(y(rmax));
        }
        else
        {
            buf(1) = buf(1)+log(ap::maxrealnumber);
        }
        for(j = 0; j <= nclasses-1; j++)
        {
            v = y(j);
            if( j==rmax )
            {
                ev = 1;
            }
            else
            {
                ev = 0;
            }
            buf(2) = buf(2)+ap::sqr(v-ev);
            buf(3) = buf(3)+fabs(v-ev);
            if( ap::fp_neq(ev,0) )
            {
                buf(4) = buf(4)+fabs((v-ev)/ev);
                buf(offs+2) = buf(offs+2)+1;
            }
        }
        buf(offs+1) = buf(offs+1)+1;
    }
    else
    {
        
        //
        // Regression
        //
        nout = -nclasses;
        rmax = 0;
        for(j = 1; j <= nout-1; j++)
        {
            if( ap::fp_greater(desiredy(j),desiredy(rmax)) )
            {
                rmax = j;
            }
        }
        mmax = 0;
        for(j = 1; j <= nout-1; j++)
        {
            if( ap::fp_greater(y(j),y(mmax)) )
            {
                mmax = j;
            }
        }
        if( mmax!=rmax )
        {
            buf(0) = buf(0)+1;
        }
        for(j = 0; j <= nout-1; j++)
        {
            v = y(j);
            ev = desiredy(j);
            buf(2) = buf(2)+ap::sqr(v-ev);
            buf(3) = buf(3)+fabs(v-ev);
            if( ap::fp_neq(ev,0) )
            {
                buf(4) = buf(4)+fabs((v-ev)/ev);
                buf(offs+2) = buf(offs+2)+1;
            }
        }
        buf(offs+1) = buf(offs+1)+1;
    }
}


/*************************************************************************
See DSErrAllocate for comments on this routine.

  -- ALGLIB --
     Copyright 11.01.2009 by Bochkanov Sergey
*************************************************************************/
void dserrfinish(ap::real_1d_array& buf)
{
    int nout;
    int offs;

    offs = 5;
    nout = abs(ap::round(buf(offs)));
    if( ap::fp_neq(buf(offs+1),0) )
    {
        buf(0) = buf(0)/buf(offs+1);
        buf(1) = buf(1)/buf(offs+1);
        buf(2) = sqrt(buf(2)/(nout*buf(offs+1)));
        buf(3) = buf(3)/(nout*buf(offs+1));
    }
    if( ap::fp_neq(buf(offs+2),0) )
    {
        buf(4) = buf(4)/buf(offs+2);
    }
}


/*************************************************************************

  -- ALGLIB --
     Copyright 19.05.2008 by Bochkanov Sergey
*************************************************************************/
void dsnormalize(ap::real_2d_array& xy,
     int npoints,
     int nvars,
     int& info,
     ap::real_1d_array& means,
     ap::real_1d_array& sigmas)
{
    int i;
    int j;
    ap::real_1d_array tmp;
    double mean;
    double variance;
    double skewness;
    double kurtosis;

    
    //
    // Test parameters
    //
    if( npoints<=0||nvars<1 )
    {
        info = -1;
        return;
    }
    info = 1;
    
    //
    // Standartization
    //
    means.setbounds(0, nvars-1);
    sigmas.setbounds(0, nvars-1);
    tmp.setbounds(0, npoints-1);
    for(j = 0; j <= nvars-1; j++)
    {
        ap::vmove(&tmp(0), 1, &xy(0, j), xy.getstride(), ap::vlen(0,npoints-1));
        calculatemoments(tmp, npoints, mean, variance, skewness, kurtosis);
        means(j) = mean;
        sigmas(j) = sqrt(variance);
        if( ap::fp_eq(sigmas(j),0) )
        {
            sigmas(j) = 1;
        }
        for(i = 0; i <= npoints-1; i++)
        {
            xy(i,j) = (xy(i,j)-means(j))/sigmas(j);
        }
    }
}


/*************************************************************************

  -- ALGLIB --
     Copyright 19.05.2008 by Bochkanov Sergey
*************************************************************************/
void dsnormalizec(const ap::real_2d_array& xy,
     int npoints,
     int nvars,
     int& info,
     ap::real_1d_array& means,
     ap::real_1d_array& sigmas)
{
    int j;
    ap::real_1d_array tmp;
    double mean;
    double variance;
    double skewness;
    double kurtosis;

    
    //
    // Test parameters
    //
    if( npoints<=0||nvars<1 )
    {
        info = -1;
        return;
    }
    info = 1;
    
    //
    // Standartization
    //
    means.setbounds(0, nvars-1);
    sigmas.setbounds(0, nvars-1);
    tmp.setbounds(0, npoints-1);
    for(j = 0; j <= nvars-1; j++)
    {
        ap::vmove(&tmp(0), 1, &xy(0, j), xy.getstride(), ap::vlen(0,npoints-1));
        calculatemoments(tmp, npoints, mean, variance, skewness, kurtosis);
        means(j) = mean;
        sigmas(j) = sqrt(variance);
        if( ap::fp_eq(sigmas(j),0) )
        {
            sigmas(j) = 1;
        }
    }
}


/*************************************************************************

  -- ALGLIB --
     Copyright 19.05.2008 by Bochkanov Sergey
*************************************************************************/
double dsgetmeanmindistance(const ap::real_2d_array& xy,
     int npoints,
     int nvars)
{
    double result;
    int i;
    int j;
    ap::real_1d_array tmp;
    ap::real_1d_array tmp2;
    double v;

    
    //
    // Test parameters
    //
    if( npoints<=0||nvars<1 )
    {
        result = 0;
        return result;
    }
    
    //
    // Process
    //
    tmp.setbounds(0, npoints-1);
    for(i = 0; i <= npoints-1; i++)
    {
        tmp(i) = ap::maxrealnumber;
    }
    tmp2.setbounds(0, nvars-1);
    for(i = 0; i <= npoints-1; i++)
    {
        for(j = i+1; j <= npoints-1; j++)
        {
            ap::vmove(&tmp2(0), 1, &xy(i, 0), 1, ap::vlen(0,nvars-1));
            ap::vsub(&tmp2(0), 1, &xy(j, 0), 1, ap::vlen(0,nvars-1));
            v = ap::vdotproduct(&tmp2(0), 1, &tmp2(0), 1, ap::vlen(0,nvars-1));
            v = sqrt(v);
            tmp(i) = ap::minreal(tmp(i), v);
            tmp(j) = ap::minreal(tmp(j), v);
        }
    }
    result = 0;
    for(i = 0; i <= npoints-1; i++)
    {
        result = result+tmp(i)/npoints;
    }
    return result;
}


/*************************************************************************

  -- ALGLIB --
     Copyright 19.05.2008 by Bochkanov Sergey
*************************************************************************/
void dstie(ap::real_1d_array& a,
     int n,
     ap::integer_1d_array& ties,
     int& tiecount,
     ap::integer_1d_array& p1,
     ap::integer_1d_array& p2)
{
    int i;
    int k;
    ap::integer_1d_array tmp;

    
    //
    // Special case
    //
    if( n<=0 )
    {
        tiecount = 0;
        return;
    }
    
    //
    // Sort A
    //
    tagsort(a, n, p1, p2);
    
    //
    // Process ties
    //
    tiecount = 1;
    for(i = 1; i <= n-1; i++)
    {
        if( ap::fp_neq(a(i),a(i-1)) )
        {
            tiecount = tiecount+1;
        }
    }
    ties.setbounds(0, tiecount);
    ties(0) = 0;
    k = 1;
    for(i = 1; i <= n-1; i++)
    {
        if( ap::fp_neq(a(i),a(i-1)) )
        {
            ties(k) = i;
            k = k+1;
        }
    }
    ties(tiecount) = n;
}


/*************************************************************************

  -- ALGLIB --
     Copyright 11.12.2008 by Bochkanov Sergey
*************************************************************************/
void dstiefasti(ap::real_1d_array& a,
     ap::integer_1d_array& b,
     int n,
     ap::integer_1d_array& ties,
     int& tiecount)
{
    int i;
    int k;
    ap::integer_1d_array tmp;

    
    //
    // Special case
    //
    if( n<=0 )
    {
        tiecount = 0;
        return;
    }
    
    //
    // Sort A
    //
    tagsortfasti(a, b, n);
    
    //
    // Process ties
    //
    ties(0) = 0;
    k = 1;
    for(i = 1; i <= n-1; i++)
    {
        if( ap::fp_neq(a(i),a(i-1)) )
        {
            ties(k) = i;
            k = k+1;
        }
    }
    ties(k) = n;
    tiecount = k;
}


/*************************************************************************
Optimal partition, internal subroutine.

  -- ALGLIB --
     Copyright 22.05.2008 by Bochkanov Sergey
*************************************************************************/
void dsoptimalsplit2(ap::real_1d_array a,
     ap::integer_1d_array c,
     int n,
     int& info,
     double& threshold,
     double& pal,
     double& pbl,
     double& par,
     double& pbr,
     double& cve)
{
    int i;
    int t;
    double s;
    ap::integer_1d_array ties;
    int tiecount;
    ap::integer_1d_array p1;
    ap::integer_1d_array p2;
    int k;
    int koptimal;
    double pak;
    double pbk;
    double cvoptimal;
    double cv;

    
    //
    // Test for errors in inputs
    //
    if( n<=0 )
    {
        info = -1;
        return;
    }
    for(i = 0; i <= n-1; i++)
    {
        if( c(i)!=0&&c(i)!=1 )
        {
            info = -2;
            return;
        }
    }
    info = 1;
    
    //
    // Tie
    //
    dstie(a, n, ties, tiecount, p1, p2);
    for(i = 0; i <= n-1; i++)
    {
        if( p2(i)!=i )
        {
            t = c(i);
            c(i) = c(p2(i));
            c(p2(i)) = t;
        }
    }
    
    //
    // Special case: number of ties is 1.
    //
    // NOTE: we assume that P[i,j] equals to 0 or 1,
    //       intermediate values are not allowed.
    //
    if( tiecount==1 )
    {
        info = -3;
        return;
    }
    
    //
    // General case, number of ties > 1
    //
    // NOTE: we assume that P[i,j] equals to 0 or 1,
    //       intermediate values are not allowed.
    //
    pal = 0;
    pbl = 0;
    par = 0;
    pbr = 0;
    for(i = 0; i <= n-1; i++)
    {
        if( c(i)==0 )
        {
            par = par+1;
        }
        if( c(i)==1 )
        {
            pbr = pbr+1;
        }
    }
    koptimal = -1;
    cvoptimal = ap::maxrealnumber;
    for(k = 0; k <= tiecount-2; k++)
    {
        
        //
        // first, obtain information about K-th tie which is
        // moved from R-part to L-part
        //
        pak = 0;
        pbk = 0;
        for(i = ties(k); i <= ties(k+1)-1; i++)
        {
            if( c(i)==0 )
            {
                pak = pak+1;
            }
            if( c(i)==1 )
            {
                pbk = pbk+1;
            }
        }
        
        //
        // Calculate cross-validation CE
        //
        cv = 0;
        cv = cv-xlny(pal+pak, (pal+pak)/(pal+pak+pbl+pbk+1));
        cv = cv-xlny(pbl+pbk, (pbl+pbk)/(pal+pak+1+pbl+pbk));
        cv = cv-xlny(par-pak, (par-pak)/(par-pak+pbr-pbk+1));
        cv = cv-xlny(pbr-pbk, (pbr-pbk)/(par-pak+1+pbr-pbk));
        
        //
        // Compare with best
        //
        if( ap::fp_less(cv,cvoptimal) )
        {
            cvoptimal = cv;
            koptimal = k;
        }
        
        //
        // update
        //
        pal = pal+pak;
        pbl = pbl+pbk;
        par = par-pak;
        pbr = pbr-pbk;
    }
    cve = cvoptimal;
    threshold = 0.5*(a(ties(koptimal))+a(ties(koptimal+1)));
    pal = 0;
    pbl = 0;
    par = 0;
    pbr = 0;
    for(i = 0; i <= n-1; i++)
    {
        if( ap::fp_less(a(i),threshold) )
        {
            if( c(i)==0 )
            {
                pal = pal+1;
            }
            else
            {
                pbl = pbl+1;
            }
        }
        else
        {
            if( c(i)==0 )
            {
                par = par+1;
            }
            else
            {
                pbr = pbr+1;
            }
        }
    }
    s = pal+pbl;
    pal = pal/s;
    pbl = pbl/s;
    s = par+pbr;
    par = par/s;
    pbr = pbr/s;
}


/*************************************************************************
Optimal partition, internal subroutine. Fast version.

Accepts:
    A       array[0..N-1]       array of attributes     array[0..N-1]
    C       array[0..N-1]       array of class labels
    TiesBuf array[0..N]         temporaries (ties)
    CntBuf  array[0..2*NC-1]    temporaries (counts)
    Alpha                       centering factor (0<=alpha<=1, recommended value - 0.05)
    
Output:
    Info    error code (">0"=OK, "<0"=bad)
    RMS     training set RMS error
    CVRMS   leave-one-out RMS error
    
Note:
    content of all arrays is changed by subroutine

  -- ALGLIB --
     Copyright 11.12.2008 by Bochkanov Sergey
*************************************************************************/
void dsoptimalsplit2fast(ap::real_1d_array& a,
     ap::integer_1d_array& c,
     ap::integer_1d_array& tiesbuf,
     ap::integer_1d_array& cntbuf,
     int n,
     int nc,
     double alpha,
     int& info,
     double& threshold,
     double& rms,
     double& cvrms)
{
    int i;
    int k;
    int cl;
    int tiecount;
    double cbest;
    double cc;
    int koptimal;
    int sl;
    int sr;
    double v;
    double w;
    double x;

    
    //
    // Test for errors in inputs
    //
    if( n<=0||nc<2 )
    {
        info = -1;
        return;
    }
    for(i = 0; i <= n-1; i++)
    {
        if( c(i)<0||c(i)>=nc )
        {
            info = -2;
            return;
        }
    }
    info = 1;
    
    //
    // Tie
    //
    dstiefasti(a, c, n, tiesbuf, tiecount);
    
    //
    // Special case: number of ties is 1.
    //
    if( tiecount==1 )
    {
        info = -3;
        return;
    }
    
    //
    // General case, number of ties > 1
    //
    for(i = 0; i <= 2*nc-1; i++)
    {
        cntbuf(i) = 0;
    }
    for(i = 0; i <= n-1; i++)
    {
        cntbuf(nc+c(i)) = cntbuf(nc+c(i))+1;
    }
    koptimal = -1;
    threshold = a(n-1);
    cbest = ap::maxrealnumber;
    sl = 0;
    sr = n;
    for(k = 0; k <= tiecount-2; k++)
    {
        
        //
        // first, move Kth tie from right to left
        //
        for(i = tiesbuf(k); i <= tiesbuf(k+1)-1; i++)
        {
            cl = c(i);
            cntbuf(cl) = cntbuf(cl)+1;
            cntbuf(nc+cl) = cntbuf(nc+cl)-1;
        }
        sl = sl+(tiesbuf(k+1)-tiesbuf(k));
        sr = sr-(tiesbuf(k+1)-tiesbuf(k));
        
        //
        // Calculate RMS error
        //
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
        v = sqrt(v/(nc*n));
        
        //
        // Compare with best
        //
        x = double(2*sl)/double(sl+sr)-1;
        cc = v*(1-alpha+alpha*ap::sqr(x));
        if( ap::fp_less(cc,cbest) )
        {
            
            //
            // store split
            //
            rms = v;
            koptimal = k;
            cbest = cc;
            
            //
            // calculate CVRMS error
            //
            cvrms = 0;
            for(i = 0; i <= nc-1; i++)
            {
                if( sl>1 )
                {
                    w = cntbuf(i);
                    cvrms = cvrms+w*ap::sqr((w-1)/(sl-1)-1);
                    cvrms = cvrms+(sl-w)*ap::sqr(w/(sl-1));
                }
                else
                {
                    w = cntbuf(i);
                    cvrms = cvrms+w*ap::sqr(double(1)/double(nc)-1);
                    cvrms = cvrms+(sl-w)*ap::sqr(double(1)/double(nc));
                }
                if( sr>1 )
                {
                    w = cntbuf(nc+i);
                    cvrms = cvrms+w*ap::sqr((w-1)/(sr-1)-1);
                    cvrms = cvrms+(sr-w)*ap::sqr(w/(sr-1));
                }
                else
                {
                    w = cntbuf(nc+i);
                    cvrms = cvrms+w*ap::sqr(double(1)/double(nc)-1);
                    cvrms = cvrms+(sr-w)*ap::sqr(double(1)/double(nc));
                }
            }
            cvrms = sqrt(cvrms/(nc*n));
        }
    }
    
    //
    // Calculate threshold.
    // Code is a bit complicated because there can be such
    // numbers that 0.5(A+B) equals to A or B (if A-B=epsilon)
    //
    threshold = 0.5*(a(tiesbuf(koptimal))+a(tiesbuf(koptimal+1)));
    if( ap::fp_less_eq(threshold,a(tiesbuf(koptimal))) )
    {
        threshold = a(tiesbuf(koptimal+1));
    }
}


/*************************************************************************
Automatic non-optimal discretization, internal subroutine.

  -- ALGLIB --
     Copyright 22.05.2008 by Bochkanov Sergey
*************************************************************************/
void dssplitk(ap::real_1d_array a,
     ap::integer_1d_array c,
     int n,
     int nc,
     int kmax,
     int& info,
     ap::real_1d_array& thresholds,
     int& ni,
     double& cve)
{
    int i;
    int j;
    int j1;
    int k;
    ap::integer_1d_array ties;
    int tiecount;
    ap::integer_1d_array p1;
    ap::integer_1d_array p2;
    ap::integer_1d_array cnt;
    double v2;
    int bestk;
    double bestcve;
    ap::integer_1d_array bestsizes;
    double curcve;
    ap::integer_1d_array cursizes;

    
    //
    // Test for errors in inputs
    //
    if( n<=0||nc<2||kmax<2 )
    {
        info = -1;
        return;
    }
    for(i = 0; i <= n-1; i++)
    {
        if( c(i)<0||c(i)>=nc )
        {
            info = -2;
            return;
        }
    }
    info = 1;
    
    //
    // Tie
    //
    dstie(a, n, ties, tiecount, p1, p2);
    for(i = 0; i <= n-1; i++)
    {
        if( p2(i)!=i )
        {
            k = c(i);
            c(i) = c(p2(i));
            c(p2(i)) = k;
        }
    }
    
    //
    // Special cases
    //
    if( tiecount==1 )
    {
        info = -3;
        return;
    }
    
    //
    // General case:
    // 0. allocate arrays
    //
    kmax = ap::minint(kmax, tiecount);
    bestsizes.setbounds(0, kmax-1);
    cursizes.setbounds(0, kmax-1);
    cnt.setbounds(0, nc-1);
    
    //
    // General case:
    // 1. prepare "weak" solution (two subintervals, divided at median)
    //
    v2 = ap::maxrealnumber;
    j = -1;
    for(i = 1; i <= tiecount-1; i++)
    {
        if( ap::fp_less(fabs(ties(i)-0.5*(n-1)),v2) )
        {
            v2 = fabs(ties(i)-0.5*n);
            j = i;
        }
    }
    ap::ap_error::make_assertion(j>0, "DSSplitK: internal error #1!");
    bestk = 2;
    bestsizes(0) = ties(j);
    bestsizes(1) = n-j;
    bestcve = 0;
    for(i = 0; i <= nc-1; i++)
    {
        cnt(i) = 0;
    }
    for(i = 0; i <= j-1; i++)
    {
        tieaddc(c, ties, i, nc, cnt);
    }
    bestcve = bestcve+getcv(cnt, nc);
    for(i = 0; i <= nc-1; i++)
    {
        cnt(i) = 0;
    }
    for(i = j; i <= tiecount-1; i++)
    {
        tieaddc(c, ties, i, nc, cnt);
    }
    bestcve = bestcve+getcv(cnt, nc);
    
    //
    // General case:
    // 2. Use greedy algorithm to find sub-optimal split in O(KMax*N) time
    //
    for(k = 2; k <= kmax; k++)
    {
        
        //
        // Prepare greedy K-interval split
        //
        for(i = 0; i <= k-1; i++)
        {
            cursizes(i) = 0;
        }
        i = 0;
        j = 0;
        while(j<=tiecount-1&&i<=k-1)
        {
            
            //
            // Rule: I-th bin is empty, fill it
            //
            if( cursizes(i)==0 )
            {
                cursizes(i) = ties(j+1)-ties(j);
                j = j+1;
                continue;
            }
            
            //
            // Rule: (K-1-I) bins left, (K-1-I) ties left (1 tie per bin); next bin
            //
            if( tiecount-j==k-1-i )
            {
                i = i+1;
                continue;
            }
            
            //
            // Rule: last bin, always place in current
            //
            if( i==k-1 )
            {
                cursizes(i) = cursizes(i)+ties(j+1)-ties(j);
                j = j+1;
                continue;
            }
            
            //
            // Place J-th tie in I-th bin, or leave for I+1-th bin.
            //
            if( ap::fp_less(fabs(cursizes(i)+ties(j+1)-ties(j)-double(n)/double(k)),fabs(cursizes(i)-double(n)/double(k))) )
            {
                cursizes(i) = cursizes(i)+ties(j+1)-ties(j);
                j = j+1;
            }
            else
            {
                i = i+1;
            }
        }
        ap::ap_error::make_assertion(cursizes(k-1)!=0&&j==tiecount, "DSSplitK: internal error #1");
        
        //
        // Calculate CVE
        //
        curcve = 0;
        j = 0;
        for(i = 0; i <= k-1; i++)
        {
            for(j1 = 0; j1 <= nc-1; j1++)
            {
                cnt(j1) = 0;
            }
            for(j1 = j; j1 <= j+cursizes(i)-1; j1++)
            {
                cnt(c(j1)) = cnt(c(j1))+1;
            }
            curcve = curcve+getcv(cnt, nc);
            j = j+cursizes(i);
        }
        
        //
        // Choose best variant
        //
        if( ap::fp_less(curcve,bestcve) )
        {
            for(i = 0; i <= k-1; i++)
            {
                bestsizes(i) = cursizes(i);
            }
            bestcve = curcve;
            bestk = k;
        }
    }
    
    //
    // Transform from sizes to thresholds
    //
    cve = bestcve;
    ni = bestk;
    thresholds.setbounds(0, ni-2);
    j = bestsizes(0);
    for(i = 1; i <= bestk-1; i++)
    {
        thresholds(i-1) = 0.5*(a(j-1)+a(j));
        j = j+bestsizes(i);
    }
}


/*************************************************************************
Automatic optimal discretization, internal subroutine.

  -- ALGLIB --
     Copyright 22.05.2008 by Bochkanov Sergey
*************************************************************************/
void dsoptimalsplitk(ap::real_1d_array a,
     ap::integer_1d_array c,
     int n,
     int nc,
     int kmax,
     int& info,
     ap::real_1d_array& thresholds,
     int& ni,
     double& cve)
{
    int i;
    int j;
    int s;
    int jl;
    int jr;
    double v2;
    ap::integer_1d_array ties;
    int tiecount;
    ap::integer_1d_array p1;
    ap::integer_1d_array p2;
    double cvtemp;
    ap::integer_1d_array cnt;
    ap::integer_1d_array cnt2;
    ap::real_2d_array cv;
    ap::integer_2d_array splits;
    int k;
    int koptimal;
    double cvoptimal;

    
    //
    // Test for errors in inputs
    //
    if( n<=0||nc<2||kmax<2 )
    {
        info = -1;
        return;
    }
    for(i = 0; i <= n-1; i++)
    {
        if( c(i)<0||c(i)>=nc )
        {
            info = -2;
            return;
        }
    }
    info = 1;
    
    //
    // Tie
    //
    dstie(a, n, ties, tiecount, p1, p2);
    for(i = 0; i <= n-1; i++)
    {
        if( p2(i)!=i )
        {
            k = c(i);
            c(i) = c(p2(i));
            c(p2(i)) = k;
        }
    }
    
    //
    // Special cases
    //
    if( tiecount==1 )
    {
        info = -3;
        return;
    }
    
    //
    // General case
    // Use dynamic programming to find best split in O(KMax*NC*TieCount^2) time
    //
    kmax = ap::minint(kmax, tiecount);
    cv.setbounds(0, kmax-1, 0, tiecount-1);
    splits.setbounds(0, kmax-1, 0, tiecount-1);
    cnt.setbounds(0, nc-1);
    cnt2.setbounds(0, nc-1);
    for(j = 0; j <= nc-1; j++)
    {
        cnt(j) = 0;
    }
    for(j = 0; j <= tiecount-1; j++)
    {
        tieaddc(c, ties, j, nc, cnt);
        splits(0,j) = 0;
        cv(0,j) = getcv(cnt, nc);
    }
    for(k = 1; k <= kmax-1; k++)
    {
        for(j = 0; j <= nc-1; j++)
        {
            cnt(j) = 0;
        }
        
        //
        // Subtask size J in [K..TieCount-1]:
        // optimal K-splitting on ties from 0-th to J-th.
        //
        for(j = k; j <= tiecount-1; j++)
        {
            
            //
            // Update Cnt - let it contain classes of ties from K-th to J-th
            //
            tieaddc(c, ties, j, nc, cnt);
            
            //
            // Search for optimal split point S in [K..J]
            //
            for(i = 0; i <= nc-1; i++)
            {
                cnt2(i) = cnt(i);
            }
            cv(k,j) = cv(k-1,j-1)+getcv(cnt2, nc);
            splits(k,j) = j;
            for(s = k+1; s <= j; s++)
            {
                
                //
                // Update Cnt2 - let it contain classes of ties from S-th to J-th
                //
                tiesubc(c, ties, s-1, nc, cnt2);
                
                //
                // Calculate CVE
                //
                cvtemp = cv(k-1,s-1)+getcv(cnt2, nc);
                if( ap::fp_less(cvtemp,cv(k,j)) )
                {
                    cv(k,j) = cvtemp;
                    splits(k,j) = s;
                }
            }
        }
    }
    
    //
    // Choose best partition, output result
    //
    koptimal = -1;
    cvoptimal = ap::maxrealnumber;
    for(k = 0; k <= kmax-1; k++)
    {
        if( ap::fp_less(cv(k,tiecount-1),cvoptimal) )
        {
            cvoptimal = cv(k,tiecount-1);
            koptimal = k;
        }
    }
    ap::ap_error::make_assertion(koptimal>=0, "DSOptimalSplitK: internal error #1!");
    if( koptimal==0 )
    {
        
        //
        // Special case: best partition is one big interval.
        // Even 2-partition is not better.
        // This is possible when dealing with "weak" predictor variables.
        //
        // Make binary split as close to the median as possible.
        //
        v2 = ap::maxrealnumber;
        j = -1;
        for(i = 1; i <= tiecount-1; i++)
        {
            if( ap::fp_less(fabs(ties(i)-0.5*(n-1)),v2) )
            {
                v2 = fabs(ties(i)-0.5*(n-1));
                j = i;
            }
        }
        ap::ap_error::make_assertion(j>0, "DSOptimalSplitK: internal error #2!");
        thresholds.setbounds(0, 0);
        thresholds(0) = 0.5*(a(ties(j-1))+a(ties(j)));
        ni = 2;
        cve = 0;
        for(i = 0; i <= nc-1; i++)
        {
            cnt(i) = 0;
        }
        for(i = 0; i <= j-1; i++)
        {
            tieaddc(c, ties, i, nc, cnt);
        }
        cve = cve+getcv(cnt, nc);
        for(i = 0; i <= nc-1; i++)
        {
            cnt(i) = 0;
        }
        for(i = j; i <= tiecount-1; i++)
        {
            tieaddc(c, ties, i, nc, cnt);
        }
        cve = cve+getcv(cnt, nc);
    }
    else
    {
        
        //
        // General case: 2 or more intervals
        //
        thresholds.setbounds(0, koptimal-1);
        ni = koptimal+1;
        cve = cv(koptimal,tiecount-1);
        jl = splits(koptimal,tiecount-1);
        jr = tiecount-1;
        for(k = koptimal; k >= 1; k--)
        {
            thresholds(k-1) = 0.5*(a(ties(jl-1))+a(ties(jl)));
            jr = jl-1;
            jl = splits(k-1,jl-1);
        }
    }
}


/*************************************************************************
Subroutine prepares K-fold split of the training set.

NOTES:
    "NClasses>0" means that we have classification task.
    "NClasses<0" means regression task with -NClasses real outputs.

  -- ALGLIB --
     Copyright 11.01.2009 by Bochkanov Sergey
*************************************************************************/
static void dskfoldsplit(const ap::real_2d_array& xy,
     int npoints,
     int nclasses,
     int foldscount,
     bool stratifiedsplits,
     ap::integer_1d_array& folds)
{
    int i;
    int j;
    int k;

    
    //
    // test parameters
    //
    ap::ap_error::make_assertion(npoints>0, "DSKFoldSplit: wrong NPoints!");
    ap::ap_error::make_assertion(nclasses>1||nclasses<0, "DSKFoldSplit: wrong NClasses!");
    ap::ap_error::make_assertion(foldscount>=2&&foldscount<=npoints, "DSKFoldSplit: wrong FoldsCount!");
    ap::ap_error::make_assertion(!stratifiedsplits, "DSKFoldSplit: stratified splits are not supported!");
    
    //
    // Folds
    //
    folds.setbounds(0, npoints-1);
    for(i = 0; i <= npoints-1; i++)
    {
        folds(i) = i*foldscount/npoints;
    }
    for(i = 0; i <= npoints-2; i++)
    {
        j = i+ap::randominteger(npoints-i);
        if( j!=i )
        {
            k = folds(i);
            folds(i) = folds(j);
            folds(j) = k;
        }
    }
}


/*************************************************************************
Internal function
*************************************************************************/
static double xlny(double x, double y)
{
    double result;

    if( ap::fp_eq(x,0) )
    {
        result = 0;
    }
    else
    {
        result = x*log(y);
    }
    return result;
}


/*************************************************************************
Internal function,
returns number of samples of class I in Cnt[I]
*************************************************************************/
static double getcv(const ap::integer_1d_array& cnt, int nc)
{
    double result;
    int i;
    double s;

    s = 0;
    for(i = 0; i <= nc-1; i++)
    {
        s = s+cnt(i);
    }
    result = 0;
    for(i = 0; i <= nc-1; i++)
    {
        result = result-xlny(double(cnt(i)), cnt(i)/(s+nc-1));
    }
    return result;
}


/*************************************************************************
Internal function, adds number of samples of class I in tie NTie to Cnt[I]
*************************************************************************/
static void tieaddc(const ap::integer_1d_array& c,
     const ap::integer_1d_array& ties,
     int ntie,
     int nc,
     ap::integer_1d_array& cnt)
{
    int i;

    for(i = ties(ntie); i <= ties(ntie+1)-1; i++)
    {
        cnt(c(i)) = cnt(c(i))+1;
    }
}


/*************************************************************************
Internal function, subtracts number of samples of class I in tie NTie to Cnt[I]
*************************************************************************/
static void tiesubc(const ap::integer_1d_array& c,
     const ap::integer_1d_array& ties,
     int ntie,
     int nc,
     ap::integer_1d_array& cnt)
{
    int i;

    for(i = ties(ntie); i <= ties(ntie+1)-1; i++)
    {
        cnt(c(i)) = cnt(c(i))-1;
    }
}


/*************************************************************************
Internal function,
returns number of samples of class I in Cnt[I]
*************************************************************************/
static void tiegetc(const ap::integer_1d_array& c,
     const ap::integer_1d_array& ties,
     int ntie,
     int nc,
     ap::integer_1d_array& cnt)
{
    int i;

    for(i = 0; i <= nc-1; i++)
    {
        cnt(i) = 0;
    }
    for(i = ties(ntie); i <= ties(ntie+1)-1; i++)
    {
        cnt(c(i)) = cnt(c(i))+1;
    }
}




