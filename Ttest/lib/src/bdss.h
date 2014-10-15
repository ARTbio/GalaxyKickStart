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

#ifndef _bdss_h
#define _bdss_h

#include "ap.h"
#include "ialglib.h"

#include "tsort.h"
#include "descriptivestatistics.h"


struct cvreport
{
    double relclserror;
    double avgce;
    double rmserror;
    double avgerror;
    double avgrelerror;
};




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
void dserrallocate(int nclasses, ap::real_1d_array& buf);


/*************************************************************************
See DSErrAllocate for comments on this routine.

  -- ALGLIB --
     Copyright 11.01.2009 by Bochkanov Sergey
*************************************************************************/
void dserraccumulate(ap::real_1d_array& buf,
     const ap::real_1d_array& y,
     const ap::real_1d_array& desiredy);


/*************************************************************************
See DSErrAllocate for comments on this routine.

  -- ALGLIB --
     Copyright 11.01.2009 by Bochkanov Sergey
*************************************************************************/
void dserrfinish(ap::real_1d_array& buf);


/*************************************************************************

  -- ALGLIB --
     Copyright 19.05.2008 by Bochkanov Sergey
*************************************************************************/
void dsnormalize(ap::real_2d_array& xy,
     int npoints,
     int nvars,
     int& info,
     ap::real_1d_array& means,
     ap::real_1d_array& sigmas);


/*************************************************************************

  -- ALGLIB --
     Copyright 19.05.2008 by Bochkanov Sergey
*************************************************************************/
void dsnormalizec(const ap::real_2d_array& xy,
     int npoints,
     int nvars,
     int& info,
     ap::real_1d_array& means,
     ap::real_1d_array& sigmas);


/*************************************************************************

  -- ALGLIB --
     Copyright 19.05.2008 by Bochkanov Sergey
*************************************************************************/
double dsgetmeanmindistance(const ap::real_2d_array& xy,
     int npoints,
     int nvars);


/*************************************************************************

  -- ALGLIB --
     Copyright 19.05.2008 by Bochkanov Sergey
*************************************************************************/
void dstie(ap::real_1d_array& a,
     int n,
     ap::integer_1d_array& ties,
     int& tiecount,
     ap::integer_1d_array& p1,
     ap::integer_1d_array& p2);


/*************************************************************************

  -- ALGLIB --
     Copyright 11.12.2008 by Bochkanov Sergey
*************************************************************************/
void dstiefasti(ap::real_1d_array& a,
     ap::integer_1d_array& b,
     int n,
     ap::integer_1d_array& ties,
     int& tiecount);


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
     double& cve);


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
     double& cvrms);


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
     double& cve);


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
     double& cve);


#endif

