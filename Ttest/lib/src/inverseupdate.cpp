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
#include "inverseupdate.h"

/*************************************************************************
Inverse matrix update by the Sherman-Morrison formula

The algorithm updates matrix A^-1 when adding a number to an element
of matrix A.

Input parameters:
    InvA    -   inverse of matrix A.
                Array whose indexes range within [0..N-1, 0..N-1].
    N       -   size of matrix A.
    UpdRow  -   row where the element to be updated is stored.
    UpdColumn - column where the element to be updated is stored.
    UpdVal  -   a number to be added to the element.


Output parameters:
    InvA    -   inverse of modified matrix A.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************/
void rmatrixinvupdatesimple(ap::real_2d_array& inva,
     int n,
     int updrow,
     int updcolumn,
     double updval)
{
    ap::real_1d_array t1;
    ap::real_1d_array t2;
    int i;
    double lambda;
    double vt;

    ap::ap_error::make_assertion(updrow>=0&&updrow<n, "RMatrixInvUpdateSimple: incorrect UpdRow!");
    ap::ap_error::make_assertion(updcolumn>=0&&updcolumn<n, "RMatrixInvUpdateSimple: incorrect UpdColumn!");
    t1.setbounds(0, n-1);
    t2.setbounds(0, n-1);
    
    //
    // T1 = InvA * U
    //
    ap::vmove(&t1(0), 1, &inva(0, updrow), inva.getstride(), ap::vlen(0,n-1));
    
    //
    // T2 = v*InvA
    //
    ap::vmove(&t2(0), 1, &inva(updcolumn, 0), 1, ap::vlen(0,n-1));
    
    //
    // Lambda = v * InvA * U
    //
    lambda = updval*inva(updcolumn,updrow);
    
    //
    // InvA = InvA - correction
    //
    for(i = 0; i <= n-1; i++)
    {
        vt = updval*t1(i);
        vt = vt/(1+lambda);
        ap::vsub(&inva(i, 0), 1, &t2(0), 1, ap::vlen(0,n-1), vt);
    }
}


/*************************************************************************
Inverse matrix update by the Sherman-Morrison formula

The algorithm updates matrix A^-1 when adding a vector to a row
of matrix A.

Input parameters:
    InvA    -   inverse of matrix A.
                Array whose indexes range within [0..N-1, 0..N-1].
    N       -   size of matrix A.
    UpdRow  -   the row of A whose vector V was added.
                0 <= Row <= N-1
    V       -   the vector to be added to a row.
                Array whose index ranges within [0..N-1].

Output parameters:
    InvA    -   inverse of modified matrix A.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************/
void rmatrixinvupdaterow(ap::real_2d_array& inva,
     int n,
     int updrow,
     const ap::real_1d_array& v)
{
    ap::real_1d_array t1;
    ap::real_1d_array t2;
    int i;
    int j;
    double lambda;
    double vt;

    t1.setbounds(0, n-1);
    t2.setbounds(0, n-1);
    
    //
    // T1 = InvA * U
    //
    ap::vmove(&t1(0), 1, &inva(0, updrow), inva.getstride(), ap::vlen(0,n-1));
    
    //
    // T2 = v*InvA
    // Lambda = v * InvA * U
    //
    for(j = 0; j <= n-1; j++)
    {
        vt = ap::vdotproduct(&v(0), 1, &inva(0, j), inva.getstride(), ap::vlen(0,n-1));
        t2(j) = vt;
    }
    lambda = t2(updrow);
    
    //
    // InvA = InvA - correction
    //
    for(i = 0; i <= n-1; i++)
    {
        vt = t1(i)/(1+lambda);
        ap::vsub(&inva(i, 0), 1, &t2(0), 1, ap::vlen(0,n-1), vt);
    }
}


/*************************************************************************
Inverse matrix update by the Sherman-Morrison formula

The algorithm updates matrix A^-1 when adding a vector to a column
of matrix A.

Input parameters:
    InvA        -   inverse of matrix A.
                    Array whose indexes range within [0..N-1, 0..N-1].
    N           -   size of matrix A.
    UpdColumn   -   the column of A whose vector U was added.
                    0 <= UpdColumn <= N-1
    U           -   the vector to be added to a column.
                    Array whose index ranges within [0..N-1].

Output parameters:
    InvA        -   inverse of modified matrix A.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************/
void rmatrixinvupdatecolumn(ap::real_2d_array& inva,
     int n,
     int updcolumn,
     const ap::real_1d_array& u)
{
    ap::real_1d_array t1;
    ap::real_1d_array t2;
    int i;
    double lambda;
    double vt;

    t1.setbounds(0, n-1);
    t2.setbounds(0, n-1);
    
    //
    // T1 = InvA * U
    // Lambda = v * InvA * U
    //
    for(i = 0; i <= n-1; i++)
    {
        vt = ap::vdotproduct(&inva(i, 0), 1, &u(0), 1, ap::vlen(0,n-1));
        t1(i) = vt;
    }
    lambda = t1(updcolumn);
    
    //
    // T2 = v*InvA
    //
    ap::vmove(&t2(0), 1, &inva(updcolumn, 0), 1, ap::vlen(0,n-1));
    
    //
    // InvA = InvA - correction
    //
    for(i = 0; i <= n-1; i++)
    {
        vt = t1(i)/(1+lambda);
        ap::vsub(&inva(i, 0), 1, &t2(0), 1, ap::vlen(0,n-1), vt);
    }
}


/*************************************************************************
Inverse matrix update by the Sherman-Morrison formula

The algorithm computes the inverse of matrix A+u*v’ by using the given matrix
A^-1 and the vectors u and v.

Input parameters:
    InvA    -   inverse of matrix A.
                Array whose indexes range within [0..N-1, 0..N-1].
    N       -   size of matrix A.
    U       -   the vector modifying the matrix.
                Array whose index ranges within [0..N-1].
    V       -   the vector modifying the matrix.
                Array whose index ranges within [0..N-1].

Output parameters:
    InvA - inverse of matrix A + u*v'.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************/
void rmatrixinvupdateuv(ap::real_2d_array& inva,
     int n,
     const ap::real_1d_array& u,
     const ap::real_1d_array& v)
{
    ap::real_1d_array t1;
    ap::real_1d_array t2;
    int i;
    int j;
    double lambda;
    double vt;

    t1.setbounds(0, n-1);
    t2.setbounds(0, n-1);
    
    //
    // T1 = InvA * U
    // Lambda = v * T1
    //
    for(i = 0; i <= n-1; i++)
    {
        vt = ap::vdotproduct(&inva(i, 0), 1, &u(0), 1, ap::vlen(0,n-1));
        t1(i) = vt;
    }
    lambda = ap::vdotproduct(&v(0), 1, &t1(0), 1, ap::vlen(0,n-1));
    
    //
    // T2 = v*InvA
    //
    for(j = 0; j <= n-1; j++)
    {
        vt = ap::vdotproduct(&v(0), 1, &inva(0, j), inva.getstride(), ap::vlen(0,n-1));
        t2(j) = vt;
    }
    
    //
    // InvA = InvA - correction
    //
    for(i = 0; i <= n-1; i++)
    {
        vt = t1(i)/(1+lambda);
        ap::vsub(&inva(i, 0), 1, &t2(0), 1, ap::vlen(0,n-1), vt);
    }
}


void shermanmorrisonsimpleupdate(ap::real_2d_array& inva,
     int n,
     int updrow,
     int updcolumn,
     double updval)
{
    ap::real_1d_array t1;
    ap::real_1d_array t2;
    int i;
    double lambda;
    double vt;

    t1.setbounds(1, n);
    t2.setbounds(1, n);
    
    //
    // T1 = InvA * U
    //
    ap::vmove(&t1(1), 1, &inva(1, updrow), inva.getstride(), ap::vlen(1,n));
    
    //
    // T2 = v*InvA
    //
    ap::vmove(&t2(1), 1, &inva(updcolumn, 1), 1, ap::vlen(1,n));
    
    //
    // Lambda = v * InvA * U
    //
    lambda = updval*inva(updcolumn,updrow);
    
    //
    // InvA = InvA - correction
    //
    for(i = 1; i <= n; i++)
    {
        vt = updval*t1(i);
        vt = vt/(1+lambda);
        ap::vsub(&inva(i, 1), 1, &t2(1), 1, ap::vlen(1,n), vt);
    }
}


void shermanmorrisonupdaterow(ap::real_2d_array& inva,
     int n,
     int updrow,
     const ap::real_1d_array& v)
{
    ap::real_1d_array t1;
    ap::real_1d_array t2;
    int i;
    int j;
    double lambda;
    double vt;

    t1.setbounds(1, n);
    t2.setbounds(1, n);
    
    //
    // T1 = InvA * U
    //
    ap::vmove(&t1(1), 1, &inva(1, updrow), inva.getstride(), ap::vlen(1,n));
    
    //
    // T2 = v*InvA
    // Lambda = v * InvA * U
    //
    for(j = 1; j <= n; j++)
    {
        vt = ap::vdotproduct(&v(1), 1, &inva(1, j), inva.getstride(), ap::vlen(1,n));
        t2(j) = vt;
    }
    lambda = t2(updrow);
    
    //
    // InvA = InvA - correction
    //
    for(i = 1; i <= n; i++)
    {
        vt = t1(i)/(1+lambda);
        ap::vsub(&inva(i, 1), 1, &t2(1), 1, ap::vlen(1,n), vt);
    }
}


void shermanmorrisonupdatecolumn(ap::real_2d_array& inva,
     int n,
     int updcolumn,
     const ap::real_1d_array& u)
{
    ap::real_1d_array t1;
    ap::real_1d_array t2;
    int i;
    double lambda;
    double vt;

    t1.setbounds(1, n);
    t2.setbounds(1, n);
    
    //
    // T1 = InvA * U
    // Lambda = v * InvA * U
    //
    for(i = 1; i <= n; i++)
    {
        vt = ap::vdotproduct(&inva(i, 1), 1, &u(1), 1, ap::vlen(1,n));
        t1(i) = vt;
    }
    lambda = t1(updcolumn);
    
    //
    // T2 = v*InvA
    //
    ap::vmove(&t2(1), 1, &inva(updcolumn, 1), 1, ap::vlen(1,n));
    
    //
    // InvA = InvA - correction
    //
    for(i = 1; i <= n; i++)
    {
        vt = t1(i)/(1+lambda);
        ap::vsub(&inva(i, 1), 1, &t2(1), 1, ap::vlen(1,n), vt);
    }
}


void shermanmorrisonupdateuv(ap::real_2d_array& inva,
     int n,
     const ap::real_1d_array& u,
     const ap::real_1d_array& v)
{
    ap::real_1d_array t1;
    ap::real_1d_array t2;
    int i;
    int j;
    double lambda;
    double vt;

    t1.setbounds(1, n);
    t2.setbounds(1, n);
    
    //
    // T1 = InvA * U
    // Lambda = v * T1
    //
    for(i = 1; i <= n; i++)
    {
        vt = ap::vdotproduct(&inva(i, 1), 1, &u(1), 1, ap::vlen(1,n));
        t1(i) = vt;
    }
    lambda = ap::vdotproduct(&v(1), 1, &t1(1), 1, ap::vlen(1,n));
    
    //
    // T2 = v*InvA
    //
    for(j = 1; j <= n; j++)
    {
        vt = ap::vdotproduct(&v(1), 1, &inva(1, j), inva.getstride(), ap::vlen(1,n));
        t2(j) = vt;
    }
    
    //
    // InvA = InvA - correction
    //
    for(i = 1; i <= n; i++)
    {
        vt = t1(i)/(1+lambda);
        ap::vsub(&inva(i, 1), 1, &t2(1), 1, ap::vlen(1,n), vt);
    }
}




