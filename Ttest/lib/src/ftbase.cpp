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
#include "ftbase.h"

static const int ftbaseplanentrysize = 8;
static const int ftbasecffttask = 0;
static const int ftbaserfhttask = 1;
static const int ftbaserffttask = 2;
static const int fftcooleytukeyplan = 0;
static const int fftbluesteinplan = 1;
static const int fftcodeletplan = 2;
static const int fhtcooleytukeyplan = 3;
static const int fhtcodeletplan = 4;
static const int fftrealcooleytukeyplan = 5;
static const int fftemptyplan = 6;
static const int fhtn2plan = 999;
static const int ftbaseupdatetw = 4;
static const int ftbasecodeletmax = 5;
static const int ftbasecodeletrecommended = 5;
static const double ftbaseinefficiencyfactor = 1.3;
static const int ftbasemaxsmoothfactor = 5;

static void ftbasegenerateplanrec(int n,
     int tasktype,
     ftplan& plan,
     int& plansize,
     int& precomputedsize,
     int& planarraysize,
     int& tmpmemsize,
     int& stackmemsize,
     int stackptr);
static void ftbaseprecomputeplanrec(ftplan& plan,
     int entryoffset,
     int stackptr);
static void ffttwcalc(ap::real_1d_array& a, int aoffset, int n1, int n2);
static void internalcomplexlintranspose(ap::real_1d_array& a,
     int m,
     int n,
     int astart,
     ap::real_1d_array& buf);
static void internalreallintranspose(ap::real_1d_array& a,
     int m,
     int n,
     int astart,
     ap::real_1d_array& buf);
static void ffticltrec(ap::real_1d_array& a,
     int astart,
     int astride,
     ap::real_1d_array& b,
     int bstart,
     int bstride,
     int m,
     int n);
static void fftirltrec(ap::real_1d_array& a,
     int astart,
     int astride,
     ap::real_1d_array& b,
     int bstart,
     int bstride,
     int m,
     int n);
static void ftbasefindsmoothrec(int n, int seed, int leastfactor, int& best);
static void fftarrayresize(ap::integer_1d_array& a, int& asize, int newasize);
static void reffht(ap::real_1d_array& a, int n, int offs);

/*************************************************************************
This subroutine generates FFT plan - a decomposition of a N-length FFT to
the more simpler operations. Plan consists of the root entry and the child
entries.

Subroutine parameters:
    N               task size
    
Output parameters:
    Plan            plan

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
void ftbasegeneratecomplexfftplan(int n, ftplan& plan)
{
    int planarraysize;
    int plansize;
    int precomputedsize;
    int tmpmemsize;
    int stackmemsize;
    int stackptr;

    planarraysize = 1;
    plansize = 0;
    precomputedsize = 0;
    stackmemsize = 0;
    stackptr = 0;
    tmpmemsize = 2*n;
    plan.plan.setlength(planarraysize);
    ftbasegenerateplanrec(n, ftbasecffttask, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr);
    ap::ap_error::make_assertion(stackptr==0, "Internal error in FTBaseGenerateComplexFFTPlan: stack ptr!");
    plan.stackbuf.setlength(ap::maxint(stackmemsize, 1));
    plan.tmpbuf.setlength(ap::maxint(tmpmemsize, 1));
    plan.precomputed.setlength(ap::maxint(precomputedsize, 1));
    stackptr = 0;
    ftbaseprecomputeplanrec(plan, 0, stackptr);
    ap::ap_error::make_assertion(stackptr==0, "Internal error in FTBaseGenerateComplexFFTPlan: stack ptr!");
}


/*************************************************************************
Generates real FFT plan
*************************************************************************/
void ftbasegeneraterealfftplan(int n, ftplan& plan)
{
    int planarraysize;
    int plansize;
    int precomputedsize;
    int tmpmemsize;
    int stackmemsize;
    int stackptr;

    planarraysize = 1;
    plansize = 0;
    precomputedsize = 0;
    stackmemsize = 0;
    stackptr = 0;
    tmpmemsize = 2*n;
    plan.plan.setlength(planarraysize);
    ftbasegenerateplanrec(n, ftbaserffttask, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr);
    ap::ap_error::make_assertion(stackptr==0, "Internal error in FTBaseGenerateRealFFTPlan: stack ptr!");
    plan.stackbuf.setlength(ap::maxint(stackmemsize, 1));
    plan.tmpbuf.setlength(ap::maxint(tmpmemsize, 1));
    plan.precomputed.setlength(ap::maxint(precomputedsize, 1));
    stackptr = 0;
    ftbaseprecomputeplanrec(plan, 0, stackptr);
    ap::ap_error::make_assertion(stackptr==0, "Internal error in FTBaseGenerateRealFFTPlan: stack ptr!");
}


/*************************************************************************
Generates real FHT plan
*************************************************************************/
void ftbasegeneraterealfhtplan(int n, ftplan& plan)
{
    int planarraysize;
    int plansize;
    int precomputedsize;
    int tmpmemsize;
    int stackmemsize;
    int stackptr;

    planarraysize = 1;
    plansize = 0;
    precomputedsize = 0;
    stackmemsize = 0;
    stackptr = 0;
    tmpmemsize = n;
    plan.plan.setlength(planarraysize);
    ftbasegenerateplanrec(n, ftbaserfhttask, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr);
    ap::ap_error::make_assertion(stackptr==0, "Internal error in FTBaseGenerateRealFHTPlan: stack ptr!");
    plan.stackbuf.setlength(ap::maxint(stackmemsize, 1));
    plan.tmpbuf.setlength(ap::maxint(tmpmemsize, 1));
    plan.precomputed.setlength(ap::maxint(precomputedsize, 1));
    stackptr = 0;
    ftbaseprecomputeplanrec(plan, 0, stackptr);
    ap::ap_error::make_assertion(stackptr==0, "Internal error in FTBaseGenerateRealFHTPlan: stack ptr!");
}


/*************************************************************************
This subroutine executes FFT/FHT plan.

If Plan is a:
* complex FFT plan  -   sizeof(A)=2*N,
                        A contains interleaved real/imaginary values
* real FFT plan     -   sizeof(A)=2*N,
                        A contains real values interleaved with zeros
* real FHT plan     -   sizeof(A)=2*N,
                        A contains real values interleaved with zeros

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
void ftbaseexecuteplan(ap::real_1d_array& a,
     int aoffset,
     int n,
     ftplan& plan)
{
    int stackptr;

    stackptr = 0;
    ftbaseexecuteplanrec(a, aoffset, plan, 0, stackptr);
}


/*************************************************************************
Recurrent subroutine for the FTBaseExecutePlan

Parameters:
    A           FFT'ed array
    AOffset     offset of the FFT'ed part (distance is measured in doubles)

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
void ftbaseexecuteplanrec(ap::real_1d_array& a,
     int aoffset,
     ftplan& plan,
     int entryoffset,
     int stackptr)
{
    int i;
    int j;
    int k;
    int n1;
    int n2;
    int n;
    int m;
    int offs;
    int offs1;
    int offs2;
    int offsa;
    int offsb;
    int offsp;
    double hk;
    double hnk;
    double x;
    double y;
    double bx;
    double by;
    ap::real_1d_array emptyarray;
    double a0x;
    double a0y;
    double a1x;
    double a1y;
    double a2x;
    double a2y;
    double a3x;
    double a3y;
    double v0;
    double v1;
    double v2;
    double v3;
    double t1x;
    double t1y;
    double t2x;
    double t2y;
    double t3x;
    double t3y;
    double t4x;
    double t4y;
    double t5x;
    double t5y;
    double m1x;
    double m1y;
    double m2x;
    double m2y;
    double m3x;
    double m3y;
    double m4x;
    double m4y;
    double m5x;
    double m5y;
    double s1x;
    double s1y;
    double s2x;
    double s2y;
    double s3x;
    double s3y;
    double s4x;
    double s4y;
    double s5x;
    double s5y;
    double c1;
    double c2;
    double c3;
    double c4;
    double c5;
    ap::real_1d_array tmp;

    if( plan.plan(entryoffset+3)==fftemptyplan )
    {
        return;
    }
    if( plan.plan(entryoffset+3)==fftcooleytukeyplan )
    {
        
        //
        // Cooley-Tukey plan
        // * transposition
        // * row-wise FFT
        // * twiddle factors:
        //   - TwBase is a basis twiddle factor for I=1, J=1
        //   - TwRow is a twiddle factor for a second element in a row (J=1)
        //   - Tw is a twiddle factor for a current element
        // * transposition again
        // * row-wise FFT again
        //
        n1 = plan.plan(entryoffset+1);
        n2 = plan.plan(entryoffset+2);
        internalcomplexlintranspose(a, n1, n2, aoffset, plan.tmpbuf);
        for(i = 0; i <= n2-1; i++)
        {
            ftbaseexecuteplanrec(a, aoffset+i*n1*2, plan, plan.plan(entryoffset+5), stackptr);
        }
        ffttwcalc(a, aoffset, n1, n2);
        internalcomplexlintranspose(a, n2, n1, aoffset, plan.tmpbuf);
        for(i = 0; i <= n1-1; i++)
        {
            ftbaseexecuteplanrec(a, aoffset+i*n2*2, plan, plan.plan(entryoffset+6), stackptr);
        }
        internalcomplexlintranspose(a, n1, n2, aoffset, plan.tmpbuf);
        return;
    }
    if( plan.plan(entryoffset+3)==fftrealcooleytukeyplan )
    {
        
        //
        // Cooley-Tukey plan
        // * transposition
        // * row-wise FFT
        // * twiddle factors:
        //   - TwBase is a basis twiddle factor for I=1, J=1
        //   - TwRow is a twiddle factor for a second element in a row (J=1)
        //   - Tw is a twiddle factor for a current element
        // * transposition again
        // * row-wise FFT again
        //
        n1 = plan.plan(entryoffset+1);
        n2 = plan.plan(entryoffset+2);
        internalcomplexlintranspose(a, n2, n1, aoffset, plan.tmpbuf);
        for(i = 0; i <= n1/2-1; i++)
        {
            
            //
            // pack two adjacent smaller real FFT's together,
            // make one complex FFT,
            // unpack result
            //
            offs = aoffset+2*i*n2*2;
            for(k = 0; k <= n2-1; k++)
            {
                a(offs+2*k+1) = a(offs+2*n2+2*k+0);
            }
            ftbaseexecuteplanrec(a, offs, plan, plan.plan(entryoffset+6), stackptr);
            plan.tmpbuf(0) = a(offs+0);
            plan.tmpbuf(1) = 0;
            plan.tmpbuf(2*n2+0) = a(offs+1);
            plan.tmpbuf(2*n2+1) = 0;
            for(k = 1; k <= n2-1; k++)
            {
                offs1 = 2*k;
                offs2 = 2*n2+2*k;
                hk = a(offs+2*k+0);
                hnk = a(offs+2*(n2-k)+0);
                plan.tmpbuf(offs1+0) = +0.5*(hk+hnk);
                plan.tmpbuf(offs2+1) = -0.5*(hk-hnk);
                hk = a(offs+2*k+1);
                hnk = a(offs+2*(n2-k)+1);
                plan.tmpbuf(offs2+0) = +0.5*(hk+hnk);
                plan.tmpbuf(offs1+1) = +0.5*(hk-hnk);
            }
            ap::vmove(&a(offs), 1, &plan.tmpbuf(0), 1, ap::vlen(offs,offs+2*n2*2-1));
        }
        if( n1%2!=0 )
        {
            ftbaseexecuteplanrec(a, aoffset+(n1-1)*n2*2, plan, plan.plan(entryoffset+6), stackptr);
        }
        ffttwcalc(a, aoffset, n2, n1);
        internalcomplexlintranspose(a, n1, n2, aoffset, plan.tmpbuf);
        for(i = 0; i <= n2-1; i++)
        {
            ftbaseexecuteplanrec(a, aoffset+i*n1*2, plan, plan.plan(entryoffset+5), stackptr);
        }
        internalcomplexlintranspose(a, n2, n1, aoffset, plan.tmpbuf);
        return;
    }
    if( plan.plan(entryoffset+3)==fhtcooleytukeyplan )
    {
        
        //
        // Cooley-Tukey FHT plan:
        // * transpose                    \
        // * smaller FHT's                |
        // * pre-process                  |
        // * multiply by twiddle factors  | corresponds to multiplication by H1
        // * post-process                 |
        // * transpose again              /
        // * multiply by H2 (smaller FHT's)
        // * final transposition
        //
        // For more details see Vitezslav Vesely, "Fast algorithms
        // of Fourier and Hartley transform and their implementation in MATLAB",
        // page 31.
        //
        n1 = plan.plan(entryoffset+1);
        n2 = plan.plan(entryoffset+2);
        n = n1*n2;
        internalreallintranspose(a, n1, n2, aoffset, plan.tmpbuf);
        for(i = 0; i <= n2-1; i++)
        {
            ftbaseexecuteplanrec(a, aoffset+i*n1, plan, plan.plan(entryoffset+5), stackptr);
        }
        for(i = 0; i <= n2-1; i++)
        {
            for(j = 0; j <= n1-1; j++)
            {
                offsa = aoffset+i*n1;
                hk = a(offsa+j);
                hnk = a(offsa+(n1-j)%n1);
                offs = 2*(i*n1+j);
                plan.tmpbuf(offs+0) = -0.5*(hnk-hk);
                plan.tmpbuf(offs+1) = +0.5*(hk+hnk);
            }
        }
        ffttwcalc(plan.tmpbuf, 0, n1, n2);
        for(j = 0; j <= n1-1; j++)
        {
            a(aoffset+j) = plan.tmpbuf(2*j+0)+plan.tmpbuf(2*j+1);
        }
        if( n2%2==0 )
        {
            offs = 2*(n2/2)*n1;
            offsa = aoffset+n2/2*n1;
            for(j = 0; j <= n1-1; j++)
            {
                a(offsa+j) = plan.tmpbuf(offs+2*j+0)+plan.tmpbuf(offs+2*j+1);
            }
        }
        for(i = 1; i <= (n2+1)/2-1; i++)
        {
            offs = 2*i*n1;
            offs2 = 2*(n2-i)*n1;
            offsa = aoffset+i*n1;
            for(j = 0; j <= n1-1; j++)
            {
                a(offsa+j) = plan.tmpbuf(offs+2*j+1)+plan.tmpbuf(offs2+2*j+0);
            }
            offsa = aoffset+(n2-i)*n1;
            for(j = 0; j <= n1-1; j++)
            {
                a(offsa+j) = plan.tmpbuf(offs+2*j+0)+plan.tmpbuf(offs2+2*j+1);
            }
        }
        internalreallintranspose(a, n2, n1, aoffset, plan.tmpbuf);
        for(i = 0; i <= n1-1; i++)
        {
            ftbaseexecuteplanrec(a, aoffset+i*n2, plan, plan.plan(entryoffset+6), stackptr);
        }
        internalreallintranspose(a, n1, n2, aoffset, plan.tmpbuf);
        return;
    }
    if( plan.plan(entryoffset+3)==fhtn2plan )
    {
        
        //
        // Cooley-Tukey FHT plan
        //
        n1 = plan.plan(entryoffset+1);
        n2 = plan.plan(entryoffset+2);
        n = n1*n2;
        reffht(a, n, aoffset);
        return;
    }
    if( plan.plan(entryoffset+3)==fftcodeletplan )
    {
        n1 = plan.plan(entryoffset+1);
        n2 = plan.plan(entryoffset+2);
        n = n1*n2;
        if( n==2 )
        {
            a0x = a(aoffset+0);
            a0y = a(aoffset+1);
            a1x = a(aoffset+2);
            a1y = a(aoffset+3);
            v0 = a0x+a1x;
            v1 = a0y+a1y;
            v2 = a0x-a1x;
            v3 = a0y-a1y;
            a(aoffset+0) = v0;
            a(aoffset+1) = v1;
            a(aoffset+2) = v2;
            a(aoffset+3) = v3;
            return;
        }
        if( n==3 )
        {
            offs = plan.plan(entryoffset+7);
            c1 = plan.precomputed(offs+0);
            c2 = plan.precomputed(offs+1);
            a0x = a(aoffset+0);
            a0y = a(aoffset+1);
            a1x = a(aoffset+2);
            a1y = a(aoffset+3);
            a2x = a(aoffset+4);
            a2y = a(aoffset+5);
            t1x = a1x+a2x;
            t1y = a1y+a2y;
            a0x = a0x+t1x;
            a0y = a0y+t1y;
            m1x = c1*t1x;
            m1y = c1*t1y;
            m2x = c2*(a1y-a2y);
            m2y = c2*(a2x-a1x);
            s1x = a0x+m1x;
            s1y = a0y+m1y;
            a1x = s1x+m2x;
            a1y = s1y+m2y;
            a2x = s1x-m2x;
            a2y = s1y-m2y;
            a(aoffset+0) = a0x;
            a(aoffset+1) = a0y;
            a(aoffset+2) = a1x;
            a(aoffset+3) = a1y;
            a(aoffset+4) = a2x;
            a(aoffset+5) = a2y;
            return;
        }
        if( n==4 )
        {
            a0x = a(aoffset+0);
            a0y = a(aoffset+1);
            a1x = a(aoffset+2);
            a1y = a(aoffset+3);
            a2x = a(aoffset+4);
            a2y = a(aoffset+5);
            a3x = a(aoffset+6);
            a3y = a(aoffset+7);
            t1x = a0x+a2x;
            t1y = a0y+a2y;
            t2x = a1x+a3x;
            t2y = a1y+a3y;
            m2x = a0x-a2x;
            m2y = a0y-a2y;
            m3x = a1y-a3y;
            m3y = a3x-a1x;
            a(aoffset+0) = t1x+t2x;
            a(aoffset+1) = t1y+t2y;
            a(aoffset+4) = t1x-t2x;
            a(aoffset+5) = t1y-t2y;
            a(aoffset+2) = m2x+m3x;
            a(aoffset+3) = m2y+m3y;
            a(aoffset+6) = m2x-m3x;
            a(aoffset+7) = m2y-m3y;
            return;
        }
        if( n==5 )
        {
            offs = plan.plan(entryoffset+7);
            c1 = plan.precomputed(offs+0);
            c2 = plan.precomputed(offs+1);
            c3 = plan.precomputed(offs+2);
            c4 = plan.precomputed(offs+3);
            c5 = plan.precomputed(offs+4);
            t1x = a(aoffset+2)+a(aoffset+8);
            t1y = a(aoffset+3)+a(aoffset+9);
            t2x = a(aoffset+4)+a(aoffset+6);
            t2y = a(aoffset+5)+a(aoffset+7);
            t3x = a(aoffset+2)-a(aoffset+8);
            t3y = a(aoffset+3)-a(aoffset+9);
            t4x = a(aoffset+6)-a(aoffset+4);
            t4y = a(aoffset+7)-a(aoffset+5);
            t5x = t1x+t2x;
            t5y = t1y+t2y;
            a(aoffset+0) = a(aoffset+0)+t5x;
            a(aoffset+1) = a(aoffset+1)+t5y;
            m1x = c1*t5x;
            m1y = c1*t5y;
            m2x = c2*(t1x-t2x);
            m2y = c2*(t1y-t2y);
            m3x = -c3*(t3y+t4y);
            m3y = c3*(t3x+t4x);
            m4x = -c4*t4y;
            m4y = c4*t4x;
            m5x = -c5*t3y;
            m5y = c5*t3x;
            s3x = m3x-m4x;
            s3y = m3y-m4y;
            s5x = m3x+m5x;
            s5y = m3y+m5y;
            s1x = a(aoffset+0)+m1x;
            s1y = a(aoffset+1)+m1y;
            s2x = s1x+m2x;
            s2y = s1y+m2y;
            s4x = s1x-m2x;
            s4y = s1y-m2y;
            a(aoffset+2) = s2x+s3x;
            a(aoffset+3) = s2y+s3y;
            a(aoffset+4) = s4x+s5x;
            a(aoffset+5) = s4y+s5y;
            a(aoffset+6) = s4x-s5x;
            a(aoffset+7) = s4y-s5y;
            a(aoffset+8) = s2x-s3x;
            a(aoffset+9) = s2y-s3y;
            return;
        }
    }
    if( plan.plan(entryoffset+3)==fhtcodeletplan )
    {
        n1 = plan.plan(entryoffset+1);
        n2 = plan.plan(entryoffset+2);
        n = n1*n2;
        if( n==2 )
        {
            a0x = a(aoffset+0);
            a1x = a(aoffset+1);
            a(aoffset+0) = a0x+a1x;
            a(aoffset+1) = a0x-a1x;
            return;
        }
        if( n==3 )
        {
            offs = plan.plan(entryoffset+7);
            c1 = plan.precomputed(offs+0);
            c2 = plan.precomputed(offs+1);
            a0x = a(aoffset+0);
            a1x = a(aoffset+1);
            a2x = a(aoffset+2);
            t1x = a1x+a2x;
            a0x = a0x+t1x;
            m1x = c1*t1x;
            m2y = c2*(a2x-a1x);
            s1x = a0x+m1x;
            a(aoffset+0) = a0x;
            a(aoffset+1) = s1x-m2y;
            a(aoffset+2) = s1x+m2y;
            return;
        }
        if( n==4 )
        {
            a0x = a(aoffset+0);
            a1x = a(aoffset+1);
            a2x = a(aoffset+2);
            a3x = a(aoffset+3);
            t1x = a0x+a2x;
            t2x = a1x+a3x;
            m2x = a0x-a2x;
            m3y = a3x-a1x;
            a(aoffset+0) = t1x+t2x;
            a(aoffset+1) = m2x-m3y;
            a(aoffset+2) = t1x-t2x;
            a(aoffset+3) = m2x+m3y;
            return;
        }
        if( n==5 )
        {
            offs = plan.plan(entryoffset+7);
            c1 = plan.precomputed(offs+0);
            c2 = plan.precomputed(offs+1);
            c3 = plan.precomputed(offs+2);
            c4 = plan.precomputed(offs+3);
            c5 = plan.precomputed(offs+4);
            t1x = a(aoffset+1)+a(aoffset+4);
            t2x = a(aoffset+2)+a(aoffset+3);
            t3x = a(aoffset+1)-a(aoffset+4);
            t4x = a(aoffset+3)-a(aoffset+2);
            t5x = t1x+t2x;
            v0 = a(aoffset+0)+t5x;
            a(aoffset+0) = v0;
            m2x = c2*(t1x-t2x);
            m3y = c3*(t3x+t4x);
            s3y = m3y-c4*t4x;
            s5y = m3y+c5*t3x;
            s1x = v0+c1*t5x;
            s2x = s1x+m2x;
            s4x = s1x-m2x;
            a(aoffset+1) = s2x-s3y;
            a(aoffset+2) = s4x-s5y;
            a(aoffset+3) = s4x+s5y;
            a(aoffset+4) = s2x+s3y;
            return;
        }
    }
    if( plan.plan(entryoffset+3)==fftbluesteinplan )
    {
        
        //
        // Bluestein plan:
        // 1. multiply by precomputed coefficients
        // 2. make convolution: forward FFT, multiplication by precomputed FFT
        //    and backward FFT. backward FFT is represented as
        //
        //        invfft(x) = fft(x')'/M
        //
        //    for performance reasons reduction of inverse FFT to
        //    forward FFT is merged with multiplication of FFT components
        //    and last stage of Bluestein's transformation.
        // 3. post-multiplication by Bluestein factors
        //
        n = plan.plan(entryoffset+1);
        m = plan.plan(entryoffset+4);
        offs = plan.plan(entryoffset+7);
        for(i = stackptr+2*n; i <= stackptr+2*m-1; i++)
        {
            plan.stackbuf(i) = 0;
        }
        offsp = offs+2*m;
        offsa = aoffset;
        offsb = stackptr;
        for(i = 0; i <= n-1; i++)
        {
            bx = plan.precomputed(offsp+0);
            by = plan.precomputed(offsp+1);
            x = a(offsa+0);
            y = a(offsa+1);
            plan.stackbuf(offsb+0) = x*bx-y*(-by);
            plan.stackbuf(offsb+1) = x*(-by)+y*bx;
            offsp = offsp+2;
            offsa = offsa+2;
            offsb = offsb+2;
        }
        ftbaseexecuteplanrec(plan.stackbuf, stackptr, plan, plan.plan(entryoffset+5), stackptr+2*2*m);
        offsb = stackptr;
        offsp = offs;
        for(i = 0; i <= m-1; i++)
        {
            x = plan.stackbuf(offsb+0);
            y = plan.stackbuf(offsb+1);
            bx = plan.precomputed(offsp+0);
            by = plan.precomputed(offsp+1);
            plan.stackbuf(offsb+0) = x*bx-y*by;
            plan.stackbuf(offsb+1) = -(x*by+y*bx);
            offsb = offsb+2;
            offsp = offsp+2;
        }
        ftbaseexecuteplanrec(plan.stackbuf, stackptr, plan, plan.plan(entryoffset+5), stackptr+2*2*m);
        offsb = stackptr;
        offsp = offs+2*m;
        offsa = aoffset;
        for(i = 0; i <= n-1; i++)
        {
            x = +plan.stackbuf(offsb+0)/m;
            y = -plan.stackbuf(offsb+1)/m;
            bx = plan.precomputed(offsp+0);
            by = plan.precomputed(offsp+1);
            a(offsa+0) = x*bx-y*(-by);
            a(offsa+1) = x*(-by)+y*bx;
            offsp = offsp+2;
            offsa = offsa+2;
            offsb = offsb+2;
        }
        return;
    }
}


/*************************************************************************
Returns good factorization N=N1*N2.

Usually N1<=N2 (but not always - small N's may be exception).
if N1<>1 then N2<>1.

Factorization is chosen depending on task type and codelets we have.

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
void ftbasefactorize(int n, int tasktype, int& n1, int& n2)
{
    int j;

    n1 = 0;
    n2 = 0;
    
    //
    // try to find good codelet
    //
    if( n1*n2!=n )
    {
        for(j = ftbasecodeletrecommended; j >= 2; j--)
        {
            if( n%j==0 )
            {
                n1 = j;
                n2 = n/j;
                break;
            }
        }
    }
    
    //
    // try to factorize N
    //
    if( n1*n2!=n )
    {
        for(j = ftbasecodeletrecommended+1; j <= n-1; j++)
        {
            if( n%j==0 )
            {
                n1 = j;
                n2 = n/j;
                break;
            }
        }
    }
    
    //
    // looks like N is prime :(
    //
    if( n1*n2!=n )
    {
        n1 = 1;
        n2 = n;
    }
    
    //
    // normalize
    //
    if( n2==1&&n1!=1 )
    {
        n2 = n1;
        n1 = 1;
    }
}


/*************************************************************************
Is number smooth?

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
bool ftbaseissmooth(int n)
{
    bool result;
    int i;

    for(i = 2; i <= ftbasemaxsmoothfactor; i++)
    {
        while(n%i==0)
        {
            n = n/i;
        }
    }
    result = n==1;
    return result;
}


/*************************************************************************
Returns smallest smooth (divisible only by 2, 3, 5) number that is greater
than or equal to max(N,2)

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
int ftbasefindsmooth(int n)
{
    int result;
    int best;

    best = 2;
    while(best<n)
    {
        best = 2*best;
    }
    ftbasefindsmoothrec(n, 1, 2, best);
    result = best;
    return result;
}


/*************************************************************************
Returns  smallest  smooth  (divisible only by 2, 3, 5) even number that is
greater than or equal to max(N,2)

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
int ftbasefindsmootheven(int n)
{
    int result;
    int best;

    best = 2;
    while(best<n)
    {
        best = 2*best;
    }
    ftbasefindsmoothrec(n, 2, 2, best);
    result = best;
    return result;
}


/*************************************************************************
Returns estimate of FLOP count for the FFT.

It is only an estimate based on operations count for the PERFECT FFT
and relative inefficiency of the algorithm actually used.

N should be power of 2, estimates are badly wrong for non-power-of-2 N's.

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
double ftbasegetflopestimate(int n)
{
    double result;

    result = ftbaseinefficiencyfactor*(4*n*log(double(n))/log(double(2))-6*n+8);
    return result;
}


/*************************************************************************
Recurrent subroutine for the FFTGeneratePlan:

PARAMETERS:
    N                   plan size
    IsReal              whether input is real or not.
                        subroutine MUST NOT ignore this flag because real
                        inputs comes with non-initialized imaginary parts,
                        so ignoring this flag will result in corrupted output
    HalfOut             whether full output or only half of it from 0 to
                        floor(N/2) is needed. This flag may be ignored if
                        doing so will simplify calculations
    Plan                plan array
    PlanSize            size of used part (in integers)
    PrecomputedSize     size of precomputed array allocated yet
    PlanArraySize       plan array size (actual)
    TmpMemSize          temporary memory required size
    BluesteinMemSize    temporary memory required size

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbasegenerateplanrec(int n,
     int tasktype,
     ftplan& plan,
     int& plansize,
     int& precomputedsize,
     int& planarraysize,
     int& tmpmemsize,
     int& stackmemsize,
     int stackptr)
{
    int k;
    int m;
    int n1;
    int n2;
    int esize;
    int entryoffset;

    
    //
    // prepare
    //
    if( plansize+ftbaseplanentrysize>planarraysize )
    {
        fftarrayresize(plan.plan, planarraysize, 8*planarraysize);
    }
    entryoffset = plansize;
    esize = ftbaseplanentrysize;
    plansize = plansize+esize;
    
    //
    // if N=1, generate empty plan and exit
    //
    if( n==1 )
    {
        plan.plan(entryoffset+0) = esize;
        plan.plan(entryoffset+1) = -1;
        plan.plan(entryoffset+2) = -1;
        plan.plan(entryoffset+3) = fftemptyplan;
        plan.plan(entryoffset+4) = -1;
        plan.plan(entryoffset+5) = -1;
        plan.plan(entryoffset+6) = -1;
        plan.plan(entryoffset+7) = -1;
        return;
    }
    
    //
    // generate plans
    //
    ftbasefactorize(n, tasktype, n1, n2);
    if( tasktype==ftbasecffttask||tasktype==ftbaserffttask )
    {
        
        //
        // complex FFT plans
        //
        if( n1!=1 )
        {
            
            //
            // Cooley-Tukey plan (real or complex)
            //
            // Note that child plans are COMPLEX
            // (whether plan itself is complex or not).
            //
            tmpmemsize = ap::maxint(tmpmemsize, 2*n1*n2);
            plan.plan(entryoffset+0) = esize;
            plan.plan(entryoffset+1) = n1;
            plan.plan(entryoffset+2) = n2;
            if( tasktype==ftbasecffttask )
            {
                plan.plan(entryoffset+3) = fftcooleytukeyplan;
            }
            else
            {
                plan.plan(entryoffset+3) = fftrealcooleytukeyplan;
            }
            plan.plan(entryoffset+4) = 0;
            plan.plan(entryoffset+5) = plansize;
            ftbasegenerateplanrec(n1, ftbasecffttask, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr);
            plan.plan(entryoffset+6) = plansize;
            ftbasegenerateplanrec(n2, ftbasecffttask, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr);
            plan.plan(entryoffset+7) = -1;
            return;
        }
        else
        {
            if( n==2||n==3||n==4||n==5 )
            {
                
                //
                // hard-coded plan
                //
                plan.plan(entryoffset+0) = esize;
                plan.plan(entryoffset+1) = n1;
                plan.plan(entryoffset+2) = n2;
                plan.plan(entryoffset+3) = fftcodeletplan;
                plan.plan(entryoffset+4) = 0;
                plan.plan(entryoffset+5) = -1;
                plan.plan(entryoffset+6) = -1;
                plan.plan(entryoffset+7) = precomputedsize;
                if( n==3 )
                {
                    precomputedsize = precomputedsize+2;
                }
                if( n==5 )
                {
                    precomputedsize = precomputedsize+5;
                }
                return;
            }
            else
            {
                
                //
                // Bluestein's plan
                //
                // Select such M that M>=2*N-1, M is composite, and M's
                // factors are 2, 3, 5
                //
                k = 2*n2-1;
                m = ftbasefindsmooth(k);
                tmpmemsize = ap::maxint(tmpmemsize, 2*m);
                plan.plan(entryoffset+0) = esize;
                plan.plan(entryoffset+1) = n2;
                plan.plan(entryoffset+2) = -1;
                plan.plan(entryoffset+3) = fftbluesteinplan;
                plan.plan(entryoffset+4) = m;
                plan.plan(entryoffset+5) = plansize;
                stackptr = stackptr+2*2*m;
                stackmemsize = ap::maxint(stackmemsize, stackptr);
                ftbasegenerateplanrec(m, ftbasecffttask, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr);
                stackptr = stackptr-2*2*m;
                plan.plan(entryoffset+6) = -1;
                plan.plan(entryoffset+7) = precomputedsize;
                precomputedsize = precomputedsize+2*m+2*n;
                return;
            }
        }
    }
    if( tasktype==ftbaserfhttask )
    {
        
        //
        // real FHT plans
        //
        if( n1!=1 )
        {
            
            //
            // Cooley-Tukey plan
            //
            //
            tmpmemsize = ap::maxint(tmpmemsize, 2*n1*n2);
            plan.plan(entryoffset+0) = esize;
            plan.plan(entryoffset+1) = n1;
            plan.plan(entryoffset+2) = n2;
            plan.plan(entryoffset+3) = fhtcooleytukeyplan;
            plan.plan(entryoffset+4) = 0;
            plan.plan(entryoffset+5) = plansize;
            ftbasegenerateplanrec(n1, tasktype, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr);
            plan.plan(entryoffset+6) = plansize;
            ftbasegenerateplanrec(n2, tasktype, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr);
            plan.plan(entryoffset+7) = -1;
            return;
        }
        else
        {
            
            //
            // N2 plan
            //
            plan.plan(entryoffset+0) = esize;
            plan.plan(entryoffset+1) = n1;
            plan.plan(entryoffset+2) = n2;
            plan.plan(entryoffset+3) = fhtn2plan;
            plan.plan(entryoffset+4) = 0;
            plan.plan(entryoffset+5) = -1;
            plan.plan(entryoffset+6) = -1;
            plan.plan(entryoffset+7) = -1;
            if( n==2||n==3||n==4||n==5 )
            {
                
                //
                // hard-coded plan
                //
                plan.plan(entryoffset+0) = esize;
                plan.plan(entryoffset+1) = n1;
                plan.plan(entryoffset+2) = n2;
                plan.plan(entryoffset+3) = fhtcodeletplan;
                plan.plan(entryoffset+4) = 0;
                plan.plan(entryoffset+5) = -1;
                plan.plan(entryoffset+6) = -1;
                plan.plan(entryoffset+7) = precomputedsize;
                if( n==3 )
                {
                    precomputedsize = precomputedsize+2;
                }
                if( n==5 )
                {
                    precomputedsize = precomputedsize+5;
                }
                return;
            }
            return;
        }
    }
}


/*************************************************************************
Recurrent subroutine for precomputing FFT plans

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbaseprecomputeplanrec(ftplan& plan,
     int entryoffset,
     int stackptr)
{
    int i;
    int idx;
    int n1;
    int n2;
    int n;
    int m;
    int offs;
    double v;
    ap::real_1d_array emptyarray;
    double bx;
    double by;

    if( plan.plan(entryoffset+3)==fftcooleytukeyplan||plan.plan(entryoffset+3)==fftrealcooleytukeyplan||plan.plan(entryoffset+3)==fhtcooleytukeyplan )
    {
        ftbaseprecomputeplanrec(plan, plan.plan(entryoffset+5), stackptr);
        ftbaseprecomputeplanrec(plan, plan.plan(entryoffset+6), stackptr);
        return;
    }
    if( plan.plan(entryoffset+3)==fftcodeletplan||plan.plan(entryoffset+3)==fhtcodeletplan )
    {
        n1 = plan.plan(entryoffset+1);
        n2 = plan.plan(entryoffset+2);
        n = n1*n2;
        if( n==3 )
        {
            offs = plan.plan(entryoffset+7);
            plan.precomputed(offs+0) = cos(2*ap::pi()/3)-1;
            plan.precomputed(offs+1) = sin(2*ap::pi()/3);
            return;
        }
        if( n==5 )
        {
            offs = plan.plan(entryoffset+7);
            v = 2*ap::pi()/5;
            plan.precomputed(offs+0) = (cos(v)+cos(2*v))/2-1;
            plan.precomputed(offs+1) = (cos(v)-cos(2*v))/2;
            plan.precomputed(offs+2) = -sin(v);
            plan.precomputed(offs+3) = -(sin(v)+sin(2*v));
            plan.precomputed(offs+4) = sin(v)-sin(2*v);
            return;
        }
    }
    if( plan.plan(entryoffset+3)==fftbluesteinplan )
    {
        ftbaseprecomputeplanrec(plan, plan.plan(entryoffset+5), stackptr);
        n = plan.plan(entryoffset+1);
        m = plan.plan(entryoffset+4);
        offs = plan.plan(entryoffset+7);
        for(i = 0; i <= 2*m-1; i++)
        {
            plan.precomputed(offs+i) = 0;
        }
        for(i = 0; i <= n-1; i++)
        {
            bx = cos(ap::pi()*ap::sqr(double(i))/n);
            by = sin(ap::pi()*ap::sqr(double(i))/n);
            plan.precomputed(offs+2*i+0) = bx;
            plan.precomputed(offs+2*i+1) = by;
            plan.precomputed(offs+2*m+2*i+0) = bx;
            plan.precomputed(offs+2*m+2*i+1) = by;
            if( i>0 )
            {
                plan.precomputed(offs+2*(m-i)+0) = bx;
                plan.precomputed(offs+2*(m-i)+1) = by;
            }
        }
        ftbaseexecuteplanrec(plan.precomputed, offs, plan, plan.plan(entryoffset+5), stackptr);
        return;
    }
}


/*************************************************************************
Twiddle factors calculation

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ffttwcalc(ap::real_1d_array& a, int aoffset, int n1, int n2)
{
    int i;
    int j;
    int n;
    int idx;
    int offs;
    double x;
    double y;
    double twxm1;
    double twy;
    double twbasexm1;
    double twbasey;
    double twrowxm1;
    double twrowy;
    double tmpx;
    double tmpy;
    double v;

    n = n1*n2;
    v = -2*ap::pi()/n;
    twbasexm1 = -2*ap::sqr(sin(0.5*v));
    twbasey = sin(v);
    twrowxm1 = 0;
    twrowy = 0;
    for(i = 0; i <= n2-1; i++)
    {
        twxm1 = 0;
        twy = 0;
        for(j = 0; j <= n1-1; j++)
        {
            idx = i*n1+j;
            offs = aoffset+2*idx;
            x = a(offs+0);
            y = a(offs+1);
            tmpx = x*twxm1-y*twy;
            tmpy = x*twy+y*twxm1;
            a(offs+0) = x+tmpx;
            a(offs+1) = y+tmpy;
            
            //
            // update Tw: Tw(new) = Tw(old)*TwRow
            //
            if( j<n1-1 )
            {
                if( j%ftbaseupdatetw==0 )
                {
                    v = -2*ap::pi()*i*(j+1)/n;
                    twxm1 = -2*ap::sqr(sin(0.5*v));
                    twy = sin(v);
                }
                else
                {
                    tmpx = twrowxm1+twxm1*twrowxm1-twy*twrowy;
                    tmpy = twrowy+twxm1*twrowy+twy*twrowxm1;
                    twxm1 = twxm1+tmpx;
                    twy = twy+tmpy;
                }
            }
        }
        
        //
        // update TwRow: TwRow(new) = TwRow(old)*TwBase
        //
        if( i<n2-1 )
        {
            if( j%ftbaseupdatetw==0 )
            {
                v = -2*ap::pi()*(i+1)/n;
                twrowxm1 = -2*ap::sqr(sin(0.5*v));
                twrowy = sin(v);
            }
            else
            {
                tmpx = twbasexm1+twrowxm1*twbasexm1-twrowy*twbasey;
                tmpy = twbasey+twrowxm1*twbasey+twrowy*twbasexm1;
                twrowxm1 = twrowxm1+tmpx;
                twrowy = twrowy+tmpy;
            }
        }
    }
}


/*************************************************************************
Linear transpose: transpose complex matrix stored in 1-dimensional array

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void internalcomplexlintranspose(ap::real_1d_array& a,
     int m,
     int n,
     int astart,
     ap::real_1d_array& buf)
{

    ffticltrec(a, astart, n, buf, 0, m, m, n);
    ap::vmove(&a(astart), 1, &buf(0), 1, ap::vlen(astart,astart+2*m*n-1));
}


/*************************************************************************
Linear transpose: transpose real matrix stored in 1-dimensional array

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void internalreallintranspose(ap::real_1d_array& a,
     int m,
     int n,
     int astart,
     ap::real_1d_array& buf)
{

    fftirltrec(a, astart, n, buf, 0, m, m, n);
    ap::vmove(&a(astart), 1, &buf(0), 1, ap::vlen(astart,astart+m*n-1));
}


/*************************************************************************
Recurrent subroutine for a InternalComplexLinTranspose

Write A^T to B, where:
* A is m*n complex matrix stored in array A as pairs of real/image values,
  beginning from AStart position, with AStride stride
* B is n*m complex matrix stored in array B as pairs of real/image values,
  beginning from BStart position, with BStride stride
stride is measured in complex numbers, i.e. in real/image pairs.

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ffticltrec(ap::real_1d_array& a,
     int astart,
     int astride,
     ap::real_1d_array& b,
     int bstart,
     int bstride,
     int m,
     int n)
{
    int i;
    int j;
    int idx1;
    int idx2;
    int m2;
    int m1;
    int n1;

    if( m==0||n==0 )
    {
        return;
    }
    if( ap::maxint(m, n)<=8 )
    {
        m2 = 2*bstride;
        for(i = 0; i <= m-1; i++)
        {
            idx1 = bstart+2*i;
            idx2 = astart+2*i*astride;
            for(j = 0; j <= n-1; j++)
            {
                b(idx1+0) = a(idx2+0);
                b(idx1+1) = a(idx2+1);
                idx1 = idx1+m2;
                idx2 = idx2+2;
            }
        }
        return;
    }
    if( n>m )
    {
        
        //
        // New partition:
        //
        // "A^T -> B" becomes "(A1 A2)^T -> ( B1 )
        //                                  ( B2 )
        //
        n1 = n/2;
        if( n-n1>=8&&n1%8!=0 )
        {
            n1 = n1+(8-n1%8);
        }
        ap::ap_error::make_assertion(n-n1>0, "");
        ffticltrec(a, astart, astride, b, bstart, bstride, m, n1);
        ffticltrec(a, astart+2*n1, astride, b, bstart+2*n1*bstride, bstride, m, n-n1);
    }
    else
    {
        
        //
        // New partition:
        //
        // "A^T -> B" becomes "( A1 )^T -> ( B1 B2 )
        //                     ( A2 )
        //
        m1 = m/2;
        if( m-m1>=8&&m1%8!=0 )
        {
            m1 = m1+(8-m1%8);
        }
        ap::ap_error::make_assertion(m-m1>0, "");
        ffticltrec(a, astart, astride, b, bstart, bstride, m1, n);
        ffticltrec(a, astart+2*m1*astride, astride, b, bstart+2*m1, bstride, m-m1, n);
    }
}


/*************************************************************************
Recurrent subroutine for a InternalRealLinTranspose


  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void fftirltrec(ap::real_1d_array& a,
     int astart,
     int astride,
     ap::real_1d_array& b,
     int bstart,
     int bstride,
     int m,
     int n)
{
    int i;
    int j;
    int idx1;
    int idx2;
    int m1;
    int n1;

    if( m==0||n==0 )
    {
        return;
    }
    if( ap::maxint(m, n)<=8 )
    {
        for(i = 0; i <= m-1; i++)
        {
            idx1 = bstart+i;
            idx2 = astart+i*astride;
            for(j = 0; j <= n-1; j++)
            {
                b(idx1) = a(idx2);
                idx1 = idx1+bstride;
                idx2 = idx2+1;
            }
        }
        return;
    }
    if( n>m )
    {
        
        //
        // New partition:
        //
        // "A^T -> B" becomes "(A1 A2)^T -> ( B1 )
        //                                  ( B2 )
        //
        n1 = n/2;
        if( n-n1>=8&&n1%8!=0 )
        {
            n1 = n1+(8-n1%8);
        }
        ap::ap_error::make_assertion(n-n1>0, "");
        fftirltrec(a, astart, astride, b, bstart, bstride, m, n1);
        fftirltrec(a, astart+n1, astride, b, bstart+n1*bstride, bstride, m, n-n1);
    }
    else
    {
        
        //
        // New partition:
        //
        // "A^T -> B" becomes "( A1 )^T -> ( B1 B2 )
        //                     ( A2 )
        //
        m1 = m/2;
        if( m-m1>=8&&m1%8!=0 )
        {
            m1 = m1+(8-m1%8);
        }
        ap::ap_error::make_assertion(m-m1>0, "");
        fftirltrec(a, astart, astride, b, bstart, bstride, m1, n);
        fftirltrec(a, astart+m1*astride, astride, b, bstart+m1, bstride, m-m1, n);
    }
}


/*************************************************************************
recurrent subroutine for FFTFindSmoothRec

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbasefindsmoothrec(int n, int seed, int leastfactor, int& best)
{

    ap::ap_error::make_assertion(ftbasemaxsmoothfactor<=5, "FTBaseFindSmoothRec: internal error!");
    if( seed>=n )
    {
        best = ap::minint(best, seed);
        return;
    }
    if( leastfactor<=2 )
    {
        ftbasefindsmoothrec(n, seed*2, 2, best);
    }
    if( leastfactor<=3 )
    {
        ftbasefindsmoothrec(n, seed*3, 3, best);
    }
    if( leastfactor<=5 )
    {
        ftbasefindsmoothrec(n, seed*5, 5, best);
    }
}


/*************************************************************************
Internal subroutine: array resize

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void fftarrayresize(ap::integer_1d_array& a, int& asize, int newasize)
{
    ap::integer_1d_array tmp;
    int i;

    tmp.setlength(asize);
    for(i = 0; i <= asize-1; i++)
    {
        tmp(i) = a(i);
    }
    a.setlength(newasize);
    for(i = 0; i <= asize-1; i++)
    {
        a(i) = tmp(i);
    }
    asize = newasize;
}


/*************************************************************************
Reference FHT stub
*************************************************************************/
static void reffht(ap::real_1d_array& a, int n, int offs)
{
    ap::real_1d_array buf;
    int i;
    int j;
    double v;

    ap::ap_error::make_assertion(n>0, "RefFHTR1D: incorrect N!");
    buf.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        v = 0;
        for(j = 0; j <= n-1; j++)
        {
            v = v+a(offs+j)*(cos(2*ap::pi()*i*j/n)+sin(2*ap::pi()*i*j/n));
        }
        buf(i) = v;
    }
    for(i = 0; i <= n-1; i++)
    {
        a(offs+i) = buf(i);
    }
}




