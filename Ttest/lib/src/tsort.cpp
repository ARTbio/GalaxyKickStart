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
#include "tsort.h"

void tagsort(ap::real_1d_array& a,
     int n,
     ap::integer_1d_array& p1,
     ap::integer_1d_array& p2)
{
    int i;
    ap::integer_1d_array pv;
    ap::integer_1d_array vp;
    int lv;
    int lp;
    int rv;
    int rp;

    
    //
    // Special cases
    //
    if( n<=0 )
    {
        return;
    }
    if( n==1 )
    {
        p1.setbounds(0, 0);
        p2.setbounds(0, 0);
        p1(0) = 0;
        p2(0) = 0;
        return;
    }
    
    //
    // General case, N>1: prepare permutations table P1
    //
    p1.setbounds(0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        p1(i) = i;
    }
    
    //
    // General case, N>1: sort, update P1
    //
    tagsortfasti(a, p1, n);
    
    //
    // General case, N>1: fill permutations table P2
    //
    // To fill P2 we maintain two arrays:
    // * PV, Position(Value). PV[i] contains position of I-th key at the moment
    // * VP, Value(Position). VP[i] contains key which has position I at the moment
    //
    // At each step we making permutation of two items:
    //   Left, which is given by position/value pair LP/LV
    //   and Right, which is given by RP/RV
    // and updating PV[] and VP[] correspondingly.
    //
    pv.setbounds(0, n-1);
    vp.setbounds(0, n-1);
    p2.setbounds(0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        pv(i) = i;
        vp(i) = i;
    }
    for(i = 0; i <= n-1; i++)
    {
        
        //
        // calculate LP, LV, RP, RV
        //
        lp = i;
        lv = vp(lp);
        rv = p1(i);
        rp = pv(rv);
        
        //
        // Fill P2
        //
        p2(i) = rp;
        
        //
        // update PV and VP
        //
        vp(lp) = rv;
        vp(rp) = lv;
        pv(lv) = rp;
        pv(rv) = lp;
    }
}


void tagsortfasti(ap::real_1d_array& a, ap::integer_1d_array& b, int n)
{
    int i;
    int k;
    int t;
    double tmp;
    int tmpi;

    
    //
    // Special cases
    //
    if( n<=1 )
    {
        return;
    }
    
    //
    // General case, N>1: sort, update B
    //
    i = 2;
    do
    {
        t = i;
        while(t!=1)
        {
            k = t/2;
            if( ap::fp_greater_eq(a(k-1),a(t-1)) )
            {
                t = 1;
            }
            else
            {
                tmp = a(k-1);
                a(k-1) = a(t-1);
                a(t-1) = tmp;
                tmpi = b(k-1);
                b(k-1) = b(t-1);
                b(t-1) = tmpi;
                t = k;
            }
        }
        i = i+1;
    }
    while(i<=n);
    i = n-1;
    do
    {
        tmp = a(i);
        a(i) = a(0);
        a(0) = tmp;
        tmpi = b(i);
        b(i) = b(0);
        b(0) = tmpi;
        t = 1;
        while(t!=0)
        {
            k = 2*t;
            if( k>i )
            {
                t = 0;
            }
            else
            {
                if( k<i )
                {
                    if( ap::fp_greater(a(k),a(k-1)) )
                    {
                        k = k+1;
                    }
                }
                if( ap::fp_greater_eq(a(t-1),a(k-1)) )
                {
                    t = 0;
                }
                else
                {
                    tmp = a(k-1);
                    a(k-1) = a(t-1);
                    a(t-1) = tmp;
                    tmpi = b(k-1);
                    b(k-1) = b(t-1);
                    b(t-1) = tmpi;
                    t = k;
                }
            }
        }
        i = i-1;
    }
    while(i>=1);
}


void tagsortfastr(ap::real_1d_array& a, ap::real_1d_array& b, int n)
{
    int i;
    int k;
    int t;
    double tmp;
    double tmpr;

    
    //
    // Special cases
    //
    if( n<=1 )
    {
        return;
    }
    
    //
    // General case, N>1: sort, update B
    //
    i = 2;
    do
    {
        t = i;
        while(t!=1)
        {
            k = t/2;
            if( ap::fp_greater_eq(a(k-1),a(t-1)) )
            {
                t = 1;
            }
            else
            {
                tmp = a(k-1);
                a(k-1) = a(t-1);
                a(t-1) = tmp;
                tmpr = b(k-1);
                b(k-1) = b(t-1);
                b(t-1) = tmpr;
                t = k;
            }
        }
        i = i+1;
    }
    while(i<=n);
    i = n-1;
    do
    {
        tmp = a(i);
        a(i) = a(0);
        a(0) = tmp;
        tmpr = b(i);
        b(i) = b(0);
        b(0) = tmpr;
        t = 1;
        while(t!=0)
        {
            k = 2*t;
            if( k>i )
            {
                t = 0;
            }
            else
            {
                if( k<i )
                {
                    if( ap::fp_greater(a(k),a(k-1)) )
                    {
                        k = k+1;
                    }
                }
                if( ap::fp_greater_eq(a(t-1),a(k-1)) )
                {
                    t = 0;
                }
                else
                {
                    tmp = a(k-1);
                    a(k-1) = a(t-1);
                    a(t-1) = tmp;
                    tmpr = b(k-1);
                    b(k-1) = b(t-1);
                    b(t-1) = tmpr;
                    t = k;
                }
            }
        }
        i = i-1;
    }
    while(i>=1);
}


void tagsortfast(ap::real_1d_array& a, int n)
{
    int i;
    int k;
    int t;
    double tmp;

    
    //
    // Special cases
    //
    if( n<=1 )
    {
        return;
    }
    
    //
    // General case, N>1: sort, update B
    //
    i = 2;
    do
    {
        t = i;
        while(t!=1)
        {
            k = t/2;
            if( ap::fp_greater_eq(a(k-1),a(t-1)) )
            {
                t = 1;
            }
            else
            {
                tmp = a(k-1);
                a(k-1) = a(t-1);
                a(t-1) = tmp;
                t = k;
            }
        }
        i = i+1;
    }
    while(i<=n);
    i = n-1;
    do
    {
        tmp = a(i);
        a(i) = a(0);
        a(0) = tmp;
        t = 1;
        while(t!=0)
        {
            k = 2*t;
            if( k>i )
            {
                t = 0;
            }
            else
            {
                if( k<i )
                {
                    if( ap::fp_greater(a(k),a(k-1)) )
                    {
                        k = k+1;
                    }
                }
                if( ap::fp_greater_eq(a(t-1),a(k-1)) )
                {
                    t = 0;
                }
                else
                {
                    tmp = a(k-1);
                    a(k-1) = a(t-1);
                    a(t-1) = tmp;
                    t = k;
                }
            }
        }
        i = i-1;
    }
    while(i>=1);
}


/*************************************************************************
Heap operations: adds element to the heap

PARAMETERS:
    A       -   heap itself, must be at least array[0..N]
    B       -   array of integer tags, which are updated according to
                permutations in the heap
    N       -   size of the heap (without new element).
                updated on output
    VA      -   value of the element being added
    VB      -   value of the tag

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void tagheappushi(ap::real_1d_array& a,
     ap::integer_1d_array& b,
     int& n,
     double va,
     int vb)
{
    int j;
    int k;
    double v;

    if( n<0 )
    {
        return;
    }
    
    //
    // N=0 is a special case
    //
    if( n==0 )
    {
        a(0) = va;
        b(0) = vb;
        n = n+1;
        return;
    }
    
    //
    // add current point to the heap
    // (add to the bottom, then move up)
    //
    // we don't write point to the heap
    // until its final position is determined
    // (it allow us to reduce number of array access operations)
    //
    j = n;
    n = n+1;
    while(j>0)
    {
        k = (j-1)/2;
        v = a(k);
        if( ap::fp_less(v,va) )
        {
            
            //
            // swap with higher element
            //
            a(j) = v;
            b(j) = b(k);
            j = k;
        }
        else
        {
            
            //
            // element in its place. terminate.
            //
            break;
        }
    }
    a(j) = va;
    b(j) = vb;
}


/*************************************************************************
Heap operations: replaces top element with new element
(which is moved down)

PARAMETERS:
    A       -   heap itself, must be at least array[0..N-1]
    B       -   array of integer tags, which are updated according to
                permutations in the heap
    N       -   size of the heap
    VA      -   value of the element which replaces top element
    VB      -   value of the tag

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void tagheapreplacetopi(ap::real_1d_array& a,
     ap::integer_1d_array& b,
     int n,
     double va,
     int vb)
{
    int j;
    int k1;
    int k2;
    double v;
    double v1;
    double v2;

    if( n<1 )
    {
        return;
    }
    
    //
    // N=1 is a special case
    //
    if( n==1 )
    {
        a(0) = va;
        b(0) = vb;
        return;
    }
    
    //
    // move down through heap:
    // * J  -   current element
    // * K1 -   first child (always exists)
    // * K2 -   second child (may not exists)
    //
    // we don't write point to the heap
    // until its final position is determined
    // (it allow us to reduce number of array access operations)
    //
    j = 0;
    k1 = 1;
    k2 = 2;
    while(k1<n)
    {
        if( k2>=n )
        {
            
            //
            // only one child.
            //
            // swap and terminate (because this child
            // have no siblings due to heap structure)
            //
            v = a(k1);
            if( ap::fp_greater(v,va) )
            {
                a(j) = v;
                b(j) = b(k1);
                j = k1;
            }
            break;
        }
        else
        {
            
            //
            // two childs
            //
            v1 = a(k1);
            v2 = a(k2);
            if( ap::fp_greater(v1,v2) )
            {
                if( ap::fp_less(va,v1) )
                {
                    a(j) = v1;
                    b(j) = b(k1);
                    j = k1;
                }
                else
                {
                    break;
                }
            }
            else
            {
                if( ap::fp_less(va,v2) )
                {
                    a(j) = v2;
                    b(j) = b(k2);
                    j = k2;
                }
                else
                {
                    break;
                }
            }
            k1 = 2*j+1;
            k2 = 2*j+2;
        }
    }
    a(j) = va;
    b(j) = vb;
}


/*************************************************************************
Heap operations: pops top element from the heap

PARAMETERS:
    A       -   heap itself, must be at least array[0..N-1]
    B       -   array of integer tags, which are updated according to
                permutations in the heap
    N       -   size of the heap, N>=1

On output top element is moved to A[N-1], B[N-1], heap is reordered, N is
decreased by 1.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void tagheappopi(ap::real_1d_array& a, ap::integer_1d_array& b, int& n)
{
    double va;
    int vb;

    if( n<1 )
    {
        return;
    }
    
    //
    // N=1 is a special case
    //
    if( n==1 )
    {
        n = 0;
        return;
    }
    
    //
    // swap top element and last element,
    // then reorder heap
    //
    va = a(n-1);
    vb = b(n-1);
    a(n-1) = a(0);
    b(n-1) = b(0);
    n = n-1;
    tagheapreplacetopi(a, b, n, va, vb);
}




