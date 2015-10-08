/*************************************************************************
This file is a part of ALGLIB project.

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

#ifndef _safesolve_h
#define _safesolve_h

#include "ap.h"
#include "ialglib.h"

/*************************************************************************
Real implementation of CMatrixScaledTRSafeSolve

  -- ALGLIB routine --
     21.01.2010
     Bochkanov Sergey
*************************************************************************/
bool rmatrixscaledtrsafesolve(const ap::real_2d_array& a,
     double sa,
     int n,
     ap::real_1d_array& x,
     bool isupper,
     int trans,
     bool isunit,
     double maxgrowth);


/*************************************************************************
Internal subroutine for safe solution of

    SA*op(A)=b
    
where  A  is  NxN  upper/lower  triangular/unitriangular  matrix, op(A) is
either identity transform, transposition or Hermitian transposition, SA is
a scaling factor such that max(|SA*A[i,j]|) is close to 1.0 in magnutude.

This subroutine  limits  relative  growth  of  solution  (in inf-norm)  by
MaxGrowth,  returning  False  if  growth  exceeds MaxGrowth. Degenerate or
near-degenerate matrices are handled correctly (False is returned) as long
as MaxGrowth is significantly less than MaxRealNumber/norm(b).

  -- ALGLIB routine --
     21.01.2010
     Bochkanov Sergey
*************************************************************************/
bool cmatrixscaledtrsafesolve(const ap::complex_2d_array& a,
     double sa,
     int n,
     ap::complex_1d_array& x,
     bool isupper,
     int trans,
     bool isunit,
     double maxgrowth);


#endif

