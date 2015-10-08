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

#ifndef _ftbase_h
#define _ftbase_h

#include "ap.h"
#include "ialglib.h"

struct ftplan
{
    ap::integer_1d_array plan;
    ap::real_1d_array precomputed;
    ap::real_1d_array tmpbuf;
    ap::real_1d_array stackbuf;
};




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
void ftbasegeneratecomplexfftplan(int n, ftplan& plan);


/*************************************************************************
Generates real FFT plan
*************************************************************************/
void ftbasegeneraterealfftplan(int n, ftplan& plan);


/*************************************************************************
Generates real FHT plan
*************************************************************************/
void ftbasegeneraterealfhtplan(int n, ftplan& plan);


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
void ftbaseexecuteplan(ap::real_1d_array& a, int aoffset, int n, ftplan& plan);


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
     int stackptr);


/*************************************************************************
Returns good factorization N=N1*N2.

Usually N1<=N2 (but not always - small N's may be exception).
if N1<>1 then N2<>1.

Factorization is chosen depending on task type and codelets we have.

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
void ftbasefactorize(int n, int tasktype, int& n1, int& n2);


/*************************************************************************
Is number smooth?

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
bool ftbaseissmooth(int n);


/*************************************************************************
Returns smallest smooth (divisible only by 2, 3, 5) number that is greater
than or equal to max(N,2)

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
int ftbasefindsmooth(int n);


/*************************************************************************
Returns  smallest  smooth  (divisible only by 2, 3, 5) even number that is
greater than or equal to max(N,2)

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
int ftbasefindsmootheven(int n);


/*************************************************************************
Returns estimate of FLOP count for the FFT.

It is only an estimate based on operations count for the PERFECT FFT
and relative inefficiency of the algorithm actually used.

N should be power of 2, estimates are badly wrong for non-power-of-2 N's.

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
double ftbasegetflopestimate(int n);


#endif

