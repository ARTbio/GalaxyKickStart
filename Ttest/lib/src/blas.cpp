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
#include "blas.h"

double vectornorm2(const ap::real_1d_array& x, int i1, int i2)
{
    double result;
    int n;
    int ix;
    double absxi;
    double scl;
    double ssq;

    n = i2-i1+1;
    if( n<1 )
    {
        result = 0;
        return result;
    }
    if( n==1 )
    {
        result = fabs(x(i1));
        return result;
    }
    scl = 0;
    ssq = 1;
    for(ix = i1; ix <= i2; ix++)
    {
        if( ap::fp_neq(x(ix),0) )
        {
            absxi = fabs(x(ix));
            if( ap::fp_less(scl,absxi) )
            {
                ssq = 1+ssq*ap::sqr(scl/absxi);
                scl = absxi;
            }
            else
            {
                ssq = ssq+ap::sqr(absxi/scl);
            }
        }
    }
    result = scl*sqrt(ssq);
    return result;
}


int vectoridxabsmax(const ap::real_1d_array& x, int i1, int i2)
{
    int result;
    int i;
    double a;

    result = i1;
    a = fabs(x(result));
    for(i = i1+1; i <= i2; i++)
    {
        if( ap::fp_greater(fabs(x(i)),fabs(x(result))) )
        {
            result = i;
        }
    }
    return result;
}


int columnidxabsmax(const ap::real_2d_array& x, int i1, int i2, int j)
{
    int result;
    int i;
    double a;

    result = i1;
    a = fabs(x(result,j));
    for(i = i1+1; i <= i2; i++)
    {
        if( ap::fp_greater(fabs(x(i,j)),fabs(x(result,j))) )
        {
            result = i;
        }
    }
    return result;
}


int rowidxabsmax(const ap::real_2d_array& x, int j1, int j2, int i)
{
    int result;
    int j;
    double a;

    result = j1;
    a = fabs(x(i,result));
    for(j = j1+1; j <= j2; j++)
    {
        if( ap::fp_greater(fabs(x(i,j)),fabs(x(i,result))) )
        {
            result = j;
        }
    }
    return result;
}


double upperhessenberg1norm(const ap::real_2d_array& a,
     int i1,
     int i2,
     int j1,
     int j2,
     ap::real_1d_array& work)
{
    double result;
    int i;
    int j;

    ap::ap_error::make_assertion(i2-i1==j2-j1, "UpperHessenberg1Norm: I2-I1<>J2-J1!");
    for(j = j1; j <= j2; j++)
    {
        work(j) = 0;
    }
    for(i = i1; i <= i2; i++)
    {
        for(j = ap::maxint(j1, j1+i-i1-1); j <= j2; j++)
        {
            work(j) = work(j)+fabs(a(i,j));
        }
    }
    result = 0;
    for(j = j1; j <= j2; j++)
    {
        result = ap::maxreal(result, work(j));
    }
    return result;
}


void copymatrix(const ap::real_2d_array& a,
     int is1,
     int is2,
     int js1,
     int js2,
     ap::real_2d_array& b,
     int id1,
     int id2,
     int jd1,
     int jd2)
{
    int isrc;
    int idst;

    if( is1>is2||js1>js2 )
    {
        return;
    }
    ap::ap_error::make_assertion(is2-is1==id2-id1, "CopyMatrix: different sizes!");
    ap::ap_error::make_assertion(js2-js1==jd2-jd1, "CopyMatrix: different sizes!");
    for(isrc = is1; isrc <= is2; isrc++)
    {
        idst = isrc-is1+id1;
        ap::vmove(&b(idst, jd1), 1, &a(isrc, js1), 1, ap::vlen(jd1,jd2));
    }
}


void inplacetranspose(ap::real_2d_array& a,
     int i1,
     int i2,
     int j1,
     int j2,
     ap::real_1d_array& work)
{
    int i;
    int j;
    int ips;
    int jps;
    int l;

    if( i1>i2||j1>j2 )
    {
        return;
    }
    ap::ap_error::make_assertion(i1-i2==j1-j2, "InplaceTranspose error: incorrect array size!");
    for(i = i1; i <= i2-1; i++)
    {
        j = j1+i-i1;
        ips = i+1;
        jps = j1+ips-i1;
        l = i2-i;
        ap::vmove(&work(1), 1, &a(ips, j), a.getstride(), ap::vlen(1,l));
        ap::vmove(&a(ips, j), a.getstride(), &a(i, jps), 1, ap::vlen(ips,i2));
        ap::vmove(&a(i, jps), 1, &work(1), 1, ap::vlen(jps,j2));
    }
}


void copyandtranspose(const ap::real_2d_array& a,
     int is1,
     int is2,
     int js1,
     int js2,
     ap::real_2d_array& b,
     int id1,
     int id2,
     int jd1,
     int jd2)
{
    int isrc;
    int jdst;

    if( is1>is2||js1>js2 )
    {
        return;
    }
    ap::ap_error::make_assertion(is2-is1==jd2-jd1, "CopyAndTranspose: different sizes!");
    ap::ap_error::make_assertion(js2-js1==id2-id1, "CopyAndTranspose: different sizes!");
    for(isrc = is1; isrc <= is2; isrc++)
    {
        jdst = isrc-is1+jd1;
        ap::vmove(&b(id1, jdst), b.getstride(), &a(isrc, js1), 1, ap::vlen(id1,id2));
    }
}


void matrixvectormultiply(const ap::real_2d_array& a,
     int i1,
     int i2,
     int j1,
     int j2,
     bool trans,
     const ap::real_1d_array& x,
     int ix1,
     int ix2,
     double alpha,
     ap::real_1d_array& y,
     int iy1,
     int iy2,
     double beta)
{
    int i;
    double v;

    if( !trans )
    {
        
        //
        // y := alpha*A*x + beta*y;
        //
        if( i1>i2||j1>j2 )
        {
            return;
        }
        ap::ap_error::make_assertion(j2-j1==ix2-ix1, "MatrixVectorMultiply: A and X dont match!");
        ap::ap_error::make_assertion(i2-i1==iy2-iy1, "MatrixVectorMultiply: A and Y dont match!");
        
        //
        // beta*y
        //
        if( ap::fp_eq(beta,0) )
        {
            for(i = iy1; i <= iy2; i++)
            {
                y(i) = 0;
            }
        }
        else
        {
            ap::vmul(&y(iy1), 1, ap::vlen(iy1,iy2), beta);
        }
        
        //
        // alpha*A*x
        //
        for(i = i1; i <= i2; i++)
        {
            v = ap::vdotproduct(&a(i, j1), 1, &x(ix1), 1, ap::vlen(j1,j2));
            y(iy1+i-i1) = y(iy1+i-i1)+alpha*v;
        }
    }
    else
    {
        
        //
        // y := alpha*A'*x + beta*y;
        //
        if( i1>i2||j1>j2 )
        {
            return;
        }
        ap::ap_error::make_assertion(i2-i1==ix2-ix1, "MatrixVectorMultiply: A and X dont match!");
        ap::ap_error::make_assertion(j2-j1==iy2-iy1, "MatrixVectorMultiply: A and Y dont match!");
        
        //
        // beta*y
        //
        if( ap::fp_eq(beta,0) )
        {
            for(i = iy1; i <= iy2; i++)
            {
                y(i) = 0;
            }
        }
        else
        {
            ap::vmul(&y(iy1), 1, ap::vlen(iy1,iy2), beta);
        }
        
        //
        // alpha*A'*x
        //
        for(i = i1; i <= i2; i++)
        {
            v = alpha*x(ix1+i-i1);
            ap::vadd(&y(iy1), 1, &a(i, j1), 1, ap::vlen(iy1,iy2), v);
        }
    }
}


double pythag2(double x, double y)
{
    double result;
    double w;
    double xabs;
    double yabs;
    double z;

    xabs = fabs(x);
    yabs = fabs(y);
    w = ap::maxreal(xabs, yabs);
    z = ap::minreal(xabs, yabs);
    if( ap::fp_eq(z,0) )
    {
        result = w;
    }
    else
    {
        result = w*sqrt(1+ap::sqr(z/w));
    }
    return result;
}


void matrixmatrixmultiply(const ap::real_2d_array& a,
     int ai1,
     int ai2,
     int aj1,
     int aj2,
     bool transa,
     const ap::real_2d_array& b,
     int bi1,
     int bi2,
     int bj1,
     int bj2,
     bool transb,
     double alpha,
     ap::real_2d_array& c,
     int ci1,
     int ci2,
     int cj1,
     int cj2,
     double beta,
     ap::real_1d_array& work)
{
    int arows;
    int acols;
    int brows;
    int bcols;
    int crows;
    int ccols;
    int i;
    int j;
    int k;
    int l;
    int r;
    double v;

    
    //
    // Setup
    //
    if( !transa )
    {
        arows = ai2-ai1+1;
        acols = aj2-aj1+1;
    }
    else
    {
        arows = aj2-aj1+1;
        acols = ai2-ai1+1;
    }
    if( !transb )
    {
        brows = bi2-bi1+1;
        bcols = bj2-bj1+1;
    }
    else
    {
        brows = bj2-bj1+1;
        bcols = bi2-bi1+1;
    }
    ap::ap_error::make_assertion(acols==brows, "MatrixMatrixMultiply: incorrect matrix sizes!");
    if( arows<=0||acols<=0||brows<=0||bcols<=0 )
    {
        return;
    }
    crows = arows;
    ccols = bcols;
    
    //
    // Test WORK
    //
    i = ap::maxint(arows, acols);
    i = ap::maxint(brows, i);
    i = ap::maxint(i, bcols);
    work(1) = 0;
    work(i) = 0;
    
    //
    // Prepare C
    //
    if( ap::fp_eq(beta,0) )
    {
        for(i = ci1; i <= ci2; i++)
        {
            for(j = cj1; j <= cj2; j++)
            {
                c(i,j) = 0;
            }
        }
    }
    else
    {
        for(i = ci1; i <= ci2; i++)
        {
            ap::vmul(&c(i, cj1), 1, ap::vlen(cj1,cj2), beta);
        }
    }
    
    //
    // A*B
    //
    if( !transa&&!transb )
    {
        for(l = ai1; l <= ai2; l++)
        {
            for(r = bi1; r <= bi2; r++)
            {
                v = alpha*a(l,aj1+r-bi1);
                k = ci1+l-ai1;
                ap::vadd(&c(k, cj1), 1, &b(r, bj1), 1, ap::vlen(cj1,cj2), v);
            }
        }
        return;
    }
    
    //
    // A*B'
    //
    if( !transa&&transb )
    {
        if( arows*acols<brows*bcols )
        {
            for(r = bi1; r <= bi2; r++)
            {
                for(l = ai1; l <= ai2; l++)
                {
                    v = ap::vdotproduct(&a(l, aj1), 1, &b(r, bj1), 1, ap::vlen(aj1,aj2));
                    c(ci1+l-ai1,cj1+r-bi1) = c(ci1+l-ai1,cj1+r-bi1)+alpha*v;
                }
            }
            return;
        }
        else
        {
            for(l = ai1; l <= ai2; l++)
            {
                for(r = bi1; r <= bi2; r++)
                {
                    v = ap::vdotproduct(&a(l, aj1), 1, &b(r, bj1), 1, ap::vlen(aj1,aj2));
                    c(ci1+l-ai1,cj1+r-bi1) = c(ci1+l-ai1,cj1+r-bi1)+alpha*v;
                }
            }
            return;
        }
    }
    
    //
    // A'*B
    //
    if( transa&&!transb )
    {
        for(l = aj1; l <= aj2; l++)
        {
            for(r = bi1; r <= bi2; r++)
            {
                v = alpha*a(ai1+r-bi1,l);
                k = ci1+l-aj1;
                ap::vadd(&c(k, cj1), 1, &b(r, bj1), 1, ap::vlen(cj1,cj2), v);
            }
        }
        return;
    }
    
    //
    // A'*B'
    //
    if( transa&&transb )
    {
        if( arows*acols<brows*bcols )
        {
            for(r = bi1; r <= bi2; r++)
            {
                for(i = 1; i <= crows; i++)
                {
                    work(i) = 0.0;
                }
                for(l = ai1; l <= ai2; l++)
                {
                    v = alpha*b(r,bj1+l-ai1);
                    k = cj1+r-bi1;
                    ap::vadd(&work(1), 1, &a(l, aj1), 1, ap::vlen(1,crows), v);
                }
                ap::vadd(&c(ci1, k), c.getstride(), &work(1), 1, ap::vlen(ci1,ci2));
            }
            return;
        }
        else
        {
            for(l = aj1; l <= aj2; l++)
            {
                k = ai2-ai1+1;
                ap::vmove(&work(1), 1, &a(ai1, l), a.getstride(), ap::vlen(1,k));
                for(r = bi1; r <= bi2; r++)
                {
                    v = ap::vdotproduct(&work(1), 1, &b(r, bj1), 1, ap::vlen(1,k));
                    c(ci1+l-aj1,cj1+r-bi1) = c(ci1+l-aj1,cj1+r-bi1)+alpha*v;
                }
            }
            return;
        }
    }
}




