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

#ifndef _gq_h
#define _gq_h

#include "ap.h"
#include "ialglib.h"

#include "hblas.h"
#include "reflections.h"
#include "creflections.h"
#include "sblas.h"
#include "ablasf.h"
#include "ablas.h"
#include "ortfac.h"
#include "blas.h"
#include "rotations.h"
#include "hsschur.h"
#include "evd.h"
#include "gammafunc.h"


/*************************************************************************
Computation of nodes and weights for a Gauss quadrature formula

The algorithm generates the N-point Gauss quadrature formula  with  weight
function given by coefficients alpha and beta  of  a  recurrence  relation
which generates a system of orthogonal polynomials:

P-1(x)   =  0
P0(x)    =  1
Pn+1(x)  =  (x-alpha(n))*Pn(x)  -  beta(n)*Pn-1(x)

and zeroth moment Mu0

Mu0 = integral(W(x)dx,a,b)

INPUT PARAMETERS:
    Alpha   –   array[0..N-1], alpha coefficients
    Beta    –   array[0..N-1], beta coefficients
                Zero-indexed element is not used and may be arbitrary.
                Beta[I]>0.
    Mu0     –   zeroth moment of the weight function.
    N       –   number of nodes of the quadrature formula, N>=1

OUTPUT PARAMETERS:
    Info    -   error code:
                * -3    internal eigenproblem solver hasn't converged
                * -2    Beta[i]<=0
                * -1    incorrect N was passed
                *  1    OK
    X       -   array[0..N-1] - array of quadrature nodes,
                in ascending order.
    W       -   array[0..N-1] - array of quadrature weights.

  -- ALGLIB --
     Copyright 2005-2009 by Bochkanov Sergey
*************************************************************************/
void gqgeneraterec(const ap::real_1d_array& alpha,
     const ap::real_1d_array& beta,
     double mu0,
     int n,
     int& info,
     ap::real_1d_array& x,
     ap::real_1d_array& w);


/*************************************************************************
Computation of nodes and weights for a Gauss-Lobatto quadrature formula

The algorithm generates the N-point Gauss-Lobatto quadrature formula  with
weight function given by coefficients alpha and beta of a recurrence which
generates a system of orthogonal polynomials.

P-1(x)   =  0
P0(x)    =  1
Pn+1(x)  =  (x-alpha(n))*Pn(x)  -  beta(n)*Pn-1(x)

and zeroth moment Mu0

Mu0 = integral(W(x)dx,a,b)

INPUT PARAMETERS:
    Alpha   –   array[0..N-2], alpha coefficients
    Beta    –   array[0..N-2], beta coefficients.
                Zero-indexed element is not used, may be arbitrary.
                Beta[I]>0
    Mu0     –   zeroth moment of the weighting function.
    A       –   left boundary of the integration interval.
    B       –   right boundary of the integration interval.
    N       –   number of nodes of the quadrature formula, N>=3
                (including the left and right boundary nodes).

OUTPUT PARAMETERS:
    Info    -   error code:
                * -3    internal eigenproblem solver hasn't converged
                * -2    Beta[i]<=0
                * -1    incorrect N was passed
                *  1    OK
    X       -   array[0..N-1] - array of quadrature nodes,
                in ascending order.
    W       -   array[0..N-1] - array of quadrature weights.

  -- ALGLIB --
     Copyright 2005-2009 by Bochkanov Sergey
*************************************************************************/
void gqgenerategausslobattorec(ap::real_1d_array alpha,
     ap::real_1d_array beta,
     double mu0,
     double a,
     double b,
     int n,
     int& info,
     ap::real_1d_array& x,
     ap::real_1d_array& w);


/*************************************************************************
Computation of nodes and weights for a Gauss-Radau quadrature formula

The algorithm generates the N-point Gauss-Radau  quadrature  formula  with
weight function given by the coefficients alpha and  beta  of a recurrence
which generates a system of orthogonal polynomials.

P-1(x)   =  0
P0(x)    =  1
Pn+1(x)  =  (x-alpha(n))*Pn(x)  -  beta(n)*Pn-1(x)

and zeroth moment Mu0

Mu0 = integral(W(x)dx,a,b)

INPUT PARAMETERS:
    Alpha   –   array[0..N-2], alpha coefficients.
    Beta    –   array[0..N-1], beta coefficients
                Zero-indexed element is not used.
                Beta[I]>0
    Mu0     –   zeroth moment of the weighting function.
    A       –   left boundary of the integration interval.
    N       –   number of nodes of the quadrature formula, N>=2
                (including the left boundary node).

OUTPUT PARAMETERS:
    Info    -   error code:
                * -3    internal eigenproblem solver hasn't converged
                * -2    Beta[i]<=0
                * -1    incorrect N was passed
                *  1    OK
    X       -   array[0..N-1] - array of quadrature nodes,
                in ascending order.
    W       -   array[0..N-1] - array of quadrature weights.


  -- ALGLIB --
     Copyright 2005-2009 by Bochkanov Sergey
*************************************************************************/
void gqgenerategaussradaurec(ap::real_1d_array alpha,
     ap::real_1d_array beta,
     double mu0,
     double a,
     int n,
     int& info,
     ap::real_1d_array& x,
     ap::real_1d_array& w);


/*************************************************************************
Returns nodes/weights for Gauss-Legendre quadrature on [-1,1] with N
nodes.

INPUT PARAMETERS:
    N           -   number of nodes, >=1

OUTPUT PARAMETERS:
    Info        -   error code:
                    * -4    an  error   was   detected   when  calculating
                            weights/nodes.  N  is  too  large   to  obtain
                            weights/nodes  with  high   enough   accuracy.
                            Try  to   use   multiple   precision  version.
                    * -3    internal eigenproblem solver hasn't  converged
                    * -1    incorrect N was passed
                    * +1    OK
    X           -   array[0..N-1] - array of quadrature nodes,
                    in ascending order.
    W           -   array[0..N-1] - array of quadrature weights.


  -- ALGLIB --
     Copyright 12.05.2009 by Bochkanov Sergey
*************************************************************************/
void gqgenerategausslegendre(int n,
     int& info,
     ap::real_1d_array& x,
     ap::real_1d_array& w);


/*************************************************************************
Returns  nodes/weights  for  Gauss-Jacobi quadrature on [-1,1] with weight
function W(x)=Power(1-x,Alpha)*Power(1+x,Beta).

INPUT PARAMETERS:
    N           -   number of nodes, >=1
    Alpha       -   power-law coefficient, Alpha>-1
    Beta        -   power-law coefficient, Beta>-1

OUTPUT PARAMETERS:
    Info        -   error code:
                    * -4    an  error  was   detected   when   calculating
                            weights/nodes. Alpha or  Beta  are  too  close
                            to -1 to obtain weights/nodes with high enough
                            accuracy, or, may be, N is too large.  Try  to
                            use multiple precision version.
                    * -3    internal eigenproblem solver hasn't converged
                    * -1    incorrect N/Alpha/Beta was passed
                    * +1    OK
    X           -   array[0..N-1] - array of quadrature nodes,
                    in ascending order.
    W           -   array[0..N-1] - array of quadrature weights.


  -- ALGLIB --
     Copyright 12.05.2009 by Bochkanov Sergey
*************************************************************************/
void gqgenerategaussjacobi(int n,
     double alpha,
     double beta,
     int& info,
     ap::real_1d_array& x,
     ap::real_1d_array& w);


/*************************************************************************
Returns  nodes/weights  for  Gauss-Laguerre  quadrature  on  [0,+inf) with
weight function W(x)=Power(x,Alpha)*Exp(-x)

INPUT PARAMETERS:
    N           -   number of nodes, >=1
    Alpha       -   power-law coefficient, Alpha>-1

OUTPUT PARAMETERS:
    Info        -   error code:
                    * -4    an  error  was   detected   when   calculating
                            weights/nodes. Alpha is too  close  to  -1  to
                            obtain weights/nodes with high enough accuracy
                            or, may  be,  N  is  too  large.  Try  to  use
                            multiple precision version.
                    * -3    internal eigenproblem solver hasn't converged
                    * -1    incorrect N/Alpha was passed
                    * +1    OK
    X           -   array[0..N-1] - array of quadrature nodes,
                    in ascending order.
    W           -   array[0..N-1] - array of quadrature weights.


  -- ALGLIB --
     Copyright 12.05.2009 by Bochkanov Sergey
*************************************************************************/
void gqgenerategausslaguerre(int n,
     double alpha,
     int& info,
     ap::real_1d_array& x,
     ap::real_1d_array& w);


/*************************************************************************
Returns  nodes/weights  for  Gauss-Hermite  quadrature on (-inf,+inf) with
weight function W(x)=Exp(-x*x)

INPUT PARAMETERS:
    N           -   number of nodes, >=1

OUTPUT PARAMETERS:
    Info        -   error code:
                    * -4    an  error  was   detected   when   calculating
                            weights/nodes.  May be, N is too large. Try to
                            use multiple precision version.
                    * -3    internal eigenproblem solver hasn't converged
                    * -1    incorrect N/Alpha was passed
                    * +1    OK
    X           -   array[0..N-1] - array of quadrature nodes,
                    in ascending order.
    W           -   array[0..N-1] - array of quadrature weights.


  -- ALGLIB --
     Copyright 12.05.2009 by Bochkanov Sergey
*************************************************************************/
void gqgenerategausshermite(int n,
     int& info,
     ap::real_1d_array& x,
     ap::real_1d_array& w);


#endif

