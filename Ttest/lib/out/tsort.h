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

#ifndef _tsort_h
#define _tsort_h

#include "ap.h"
#include "ialglib.h"

void tagsort(ap::real_1d_array& a,
     int n,
     ap::integer_1d_array& p1,
     ap::integer_1d_array& p2);


void tagsortfasti(ap::real_1d_array& a, ap::integer_1d_array& b, int n);


void tagsortfastr(ap::real_1d_array& a, ap::real_1d_array& b, int n);


void tagsortfast(ap::real_1d_array& a, int n);


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
     int vb);


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
     int vb);


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
void tagheappopi(ap::real_1d_array& a, ap::integer_1d_array& b, int& n);


#endif

