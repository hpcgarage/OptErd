/*
 * Copyright (c) 2003-2010 University of Florida
 * Copyright (c) 2013-2015 Georgia Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * The GNU Lesser General Public License is included in this distribution
 * in the file COPYING.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <immintrin.h>

#include "erd.h"
#include "erdutil.h"

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__2D_PQ_INTEGRALS */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation calculates a full table of 2D PQ X,Y,Z */
/*                integrals using the Rys vertical recurrence scheme */
/*                VRR explained below. */
/*                The Rys weight is multiplied to the 2DX PQ integral */
/*                to reduce overall FLOP count. Note, that the Rys weight */
/*                factor needs to be introduced only three times for the */
/*                starting 2DX PQ integrals for the recurrence scheme, */
/*                namely to the (0,0), (1,0) and (0,1) elements. The */
/*                weight factor is then automatically propagated */
/*                through the vertical transfer equations (see below). */
/*                The recurrence scheme VRR is due to Rys, Dupuis and */
/*                King, J. Comp. Chem. 4, p.154-157 (1983). */
/*                   INT2D (0,0) = 1.D0    (* WEIGHT for the 2DX case) */
/*                   INT2D (1,0) = C00     (* WEIGHT for the 2DX case) */
/*                   INT2D (0,1) = D00     (* WEIGHT for the 2DX case) */
/*                   For I = 1,...,SHELLP-1 */
/*                       INT2D (I+1,0) = I * B10 * INT2D (I-1,0) */
/*                                         + C00 * INT2D (I,0) */
/*                   For K = 1,...,SHELLQ-1 */
/*                       INT2D (0,K+1) = K * B01 * INT2D (0,K-1) */
/*                                         + D00 * INT2D (0,K) */
/*                   For I = 1,...,SHELLP */
/*                       INT2D (I,1)   = I * B00 * INT2D (I-1,0) */
/*                                         + D00 * INT2D (I,0) */
/*                   For K = 2,...,SHELLQ */
/*                       INT2D (1,K)   = K * B00 * INT2D (0,K-1) */
/*                                         + C00 * INT2D (0,K) */
/*                   For K = 2,...,SHELLQ */
/*                   For I = 2,...,SHELLP */
/*                       INT2D (I,K)   = (I-1) * B10 * INT2D (I-2,K) */
/*                                         + K * B00 * INT2D (I-1,K-1) */
/*                                             + C00 * INT2D (I-1,K) */
/*                The 2D PQ integrals are calculated for all roots (info */
/*                already present in transmitted VRR coefficients!) and */
/*                for all exponent quadruples simultaneously and placed */
/*                into a 3-dimensional array. */
/*                  Input: */
/*                    SHELLx      =  maximum shell type for electrons */
/*                                   1 and 2 (x = P,Q) */
/*                    NGQEXQ      =  product of # of gaussian quadrature */
/*                                   points times exponent quadruplets */
/*                    WTS         =  all quadrature weights */
/*                    B00,B01,B10 =  VRR expansion coefficients */
/*                                   (cartesian coordinate independent) */
/*                    C00x,D00x   =  cartesian coordinate dependent */
/*                                   VRR expansion coefficients */
/*                                   (x = X,Y,Z) */
/*                    CASE2D      =  logical flag for simplifications */
/*                                   in 2D integral evaluation for */
/*                                   low quantum numbers */
/*                  Output: */
/*                    INT2Dx      =  all 2D PQ integrals for each */
/*                                   cartesian component (x = X,Y,Z) */
/* ------------------------------------------------------------------------ */
int erd__2d_pq_integrals (int shellp, int shellq,
                          int ngqexq, double *b00,
                          double *b01, double *b10, double *c00x,
                          double *c00y, double *c00z,
                          double *d00x, double *d00y,
                          double *d00z, int case2d,
                          double *int2dx, double *int2dy,
                          double *int2dz)
{
    int i, k, n, n1;
    double b0, b1;
    double weight;
/*             ...jump according to the 4 different cases that can arise: */
/*                  P-shell = s- or higher angular momentum */
/*                  Q-shell = s- or higher angular momentum */
/*                each leading to simplifications in the VRR formulas. */
/*                The case present has been evaluated outside this */
/*                routine and is transmitted via argument. */

    switch (case2d)
    {
        case 1:
            goto L1;
        case 2:
            goto L3;
        case 3:
            goto L3;
        case 4:
            goto L2;
        case 5:
            goto L4;
        case 6:
            goto L4;
        case 7:
            goto L2;
        case 8:
            goto L4;
        case 9:
            goto L4;
    }

/*             ...the case P = s-shell and Q = s-shell. */
  L1:
    for (n = 0; n < ngqexq; n+=SIMDW)
    {
        #pragma vector aligned
        #pragma simd
        for (n1 = 0; n1 < SIMDW; n1++)
        {
            int2dy[n + n1] = 1.0;
            int2dz[n + n1] = 1.0;
        }
    }

    return 0;

/*             ...the cases P = s-shell and Q >= p-shell. */
/*                Evaluate I=0 and K=0,1. */
  L2:
    for (n = 0; n < ngqexq; n+=SIMDW)
    {
        ERD_SIMD_ALIGN double int2dx_0[SIMDW], int2dx_1[SIMDW], int2dx_2[SIMDW];
        ERD_SIMD_ALIGN double int2dy_0[SIMDW], int2dy_1[SIMDW], int2dy_2[SIMDW];
        ERD_SIMD_ALIGN double int2dz_0[SIMDW], int2dz_1[SIMDW], int2dz_2[SIMDW];
        #pragma vector aligned
        #pragma simd
        for(n1 = 0; n1 < SIMDW; n1++)
        {
            int2dx_2[n1] = int2dx[n + n1];          
            int2dy[n + n1] = int2dy_2[n1] = 1.;
            int2dz[n + n1] = int2dz_2[n1] = 1.;
            int2dx[n + n1 + ngqexq] = int2dx_1[n1] = d00x[n + n1] * int2dx[n + n1];
            int2dy[n + n1 + ngqexq] = int2dy_1[n1] = d00y[n + n1];
            int2dz[n + n1 + ngqexq] = int2dz_1[n1] = d00z[n + n1];
        }

/*             ...evaluate I=0 and K=2,SHELLQ (if any). */
        for (k = 2; k <= shellq; ++k)
        {
            double k1 = k - 1;

            #pragma vector aligned
            #pragma simd
            for(n1 = 0; n1 < SIMDW; n1++)
            {
                b1 = k1 * b01[n + n1];
                int2dx_0[n1] = b1 * int2dx_2[n1] + d00x[n + n1] * int2dx_1[n1];
                int2dx_2[n1] = int2dx_1[n1];
                int2dx_1[n1] = int2dx_0[n1];
                int2dx[n + n1 + k * ngqexq] = int2dx_0[n1];

                int2dy_0[n1] = b1 * int2dy_2[n1] + d00y[n + n1] * int2dy_1[n1];
                int2dy_2[n1] = int2dy_1[n1];
                int2dy_1[n1] = int2dy_0[n1];
                int2dy[n + n1 + k * ngqexq] = int2dy_0[n1];

                int2dz_0[n1] = b1 * int2dz_2[n1] + d00z[n + n1] * int2dz_1[n1];
                int2dz_2[n1] = int2dz_1[n1];
                int2dz_1[n1] = int2dz_0[n1];
                int2dz[n + n1 + k * ngqexq] = int2dz_0[n1];
            }
        }
    }
    return 0;


/*             ...the cases P >= p-shell and Q = s-shell. */
/*                Evaluate I=0,1 and K=0. */
  L3:
    for (n = 0; n < ngqexq; n+=SIMDW)
    {
        ERD_SIMD_ALIGN double int2dx_0[SIMDW], int2dx_1[SIMDW], int2dx_2[SIMDW];
        ERD_SIMD_ALIGN double int2dy_0[SIMDW], int2dy_1[SIMDW], int2dy_2[SIMDW];
        ERD_SIMD_ALIGN double int2dz_0[SIMDW], int2dz_1[SIMDW], int2dz_2[SIMDW];

        #pragma vector aligned
        #pragma simd
        for(n1 = 0; n1 < SIMDW; n1++)
        {
            int2dx_2[n1] = int2dx[n + n1];
            int2dy[n + n1] = int2dy_2[n1] = 1.;
            int2dz[n + n1] = int2dz_2[n1] = 1.;
            int2dx[n + n1 + ngqexq] = int2dx_1[n1] = c00x[n + n1] * int2dx[n + n1];
            int2dy[n + n1 + ngqexq] = int2dy_1[n1] = c00y[n + n1];
            int2dz[n + n1 + ngqexq] = int2dz_1[n1] = c00z[n + n1];
        }
/*             ...evaluate I=2,SHELLP (if any) and K=0. */

        for (i = 2; i <= shellp; ++i)
        {
            double i1 = i - 1;

            #pragma vector aligned
            #pragma simd
            for(n1 = 0; n1 < SIMDW; n1++)
            {
                b1 = i1 * b10[n + n1];
                int2dx_0[n1] = b1 * int2dx_2[n1] + c00x[n + n1] * int2dx_1[n1];
                int2dx_2[n1] = int2dx_1[n1];
                int2dx_1[n1] = int2dx_0[n1];
                int2dx[n + n1 + i * ngqexq] = int2dx_0[n1];

                int2dy_0[n1] = b1 * int2dy_2[n1] + c00y[n + n1] * int2dy_1[n1];
                int2dy_2[n1] = int2dy_1[n1];
                int2dy_1[n1] = int2dy_0[n1];
                int2dy[n + n1 + i * ngqexq] = int2dy_0[n1];

                int2dz_0[n1] = b1 * int2dz_2[n1] + c00z[n + n1] * int2dz_1[n1];
                int2dz_2[n1] = int2dz_1[n1];
                int2dz_1[n1] = int2dz_0[n1];
                int2dz[n + n1 + i * ngqexq] = int2dz_0[n1];
            }
        }
    }
    return 0;


/*             ...the cases P >= p-shell and Q >= p-shell. */
/*                Evaluate I=0,SHELLP       I=0 */
/*                         K=0        and   K=0,SHELLQ */
  L4:
    for (n = 0; n < ngqexq; n+=SIMDW)
    {
        ERD_SIMD_ALIGN double int2dx_0[SIMDW], int2dx_i1[SIMDW], int2dx_k1[SIMDW], int2dx_2[SIMDW];
        ERD_SIMD_ALIGN double int2dy_0[SIMDW], int2dy_i1[SIMDW], int2dy_k1[SIMDW], int2dy_2[SIMDW];
        ERD_SIMD_ALIGN double int2dz_0[SIMDW], int2dz_i1[SIMDW], int2dz_k1[SIMDW], int2dz_2[SIMDW];

        #pragma vector aligned
        #pragma simd
        for(n1 = 0; n1 < SIMDW; n1++)
        {
            int2dx_2[n1] = int2dx[n + n1];
            int2dy[n + n1] = int2dy_2[n1] = 1.;
            int2dz[n + n1] = int2dz_2[n1] = 1.;
            int2dx[n + n1 + ngqexq] = int2dx_i1[n1] = c00x[n + n1] * int2dx[n + n1];
            int2dy[n + n1 + ngqexq] = int2dy_i1[n1] = c00y[n + n1];
            int2dz[n + n1 + ngqexq] = int2dz_i1[n1] = c00z[n + n1];
        }

        for (i = 2; i <= shellp; ++i)
        {
            double i1 = i - 1;

            #pragma vector aligned
            #pragma simd
            for(n1 = 0; n1 < SIMDW; n1++)
            {
                b1 = i1 * b10[n + n1];
                int2dx_0[n1] = b1 * int2dx_2[n1] + c00x[n + n1] * int2dx_i1[n1];
                int2dx_2[n1] = int2dx_i1[n1];
                int2dx_i1[n1] = int2dx_0[n1];
                int2dx[n + n1 + i * ngqexq] = int2dx_0[n1];

                int2dy_0[n1] = b1 * int2dy_2[n1] + c00y[n + n1] * int2dy_i1[n1];
                int2dy_2[n1] = int2dy_i1[n1];
                int2dy_i1[n1] = int2dy_0[n1];
                int2dy[n + n1 + i * ngqexq] = int2dy_0[n1];

                int2dz_0[n1] = b1 * int2dz_2[n1] + c00z[n + n1] * int2dz_i1[n1];
                int2dz_2[n1] = int2dz_i1[n1];
                int2dz_i1[n1] = int2dz_0[n1];
                int2dz[n + n1 + i * ngqexq] = int2dz_0[n1];
            }
        }

        #pragma vector aligned
        #pragma simd
        for(n1 = 0; n1 < SIMDW; n1++)
        {
            weight = int2dx[n + n1];
            int2dx_2[n1] = weight;
            int2dy_2[n1] = 1.;
            int2dz_2[n1] = 1.;
            int2dx[n + n1 + (shellp + 1) * ngqexq] = int2dx_k1[n1] = d00x[n + n1] * weight;
            int2dy[n + n1 + (shellp + 1) * ngqexq] = int2dy_k1[n1] = d00y[n + n1];
            int2dz[n + n1 + (shellp + 1) * ngqexq] = int2dz_k1[n1] = d00z[n + n1];
        }

        for (k = 2; k <= shellq; ++k)
        {
            double k1 = k - 1;

            #pragma vector aligned
            #pragma simd
            for(n1 = 0; n1 < SIMDW; n1++)
            {
                b1 = k1 * b01[n + n1];
                int2dx_0[n1] = b1 * int2dx_2[n1] + d00x[n + n1] * int2dx_k1[n1];
                int2dx_2[n1] = int2dx_k1[n1];
                int2dx_k1[n1] = int2dx_0[n1];
                int2dx[n + n1 + k * (shellp + 1) * ngqexq] = int2dx_0[n1];

                int2dy_0[n1] = b1 * int2dy_2[n1] + d00y[n + n1] * int2dy_k1[n1];
                int2dy_2[n1] = int2dy_k1[n1];
                int2dy_k1[n1] = int2dy_0[n1];
                int2dy[n + n1 + k * (shellp + 1) * ngqexq] = int2dy_0[n1];

                int2dz_0[n1] = b1 * int2dz_2[n1] + d00z[n + n1] * int2dz_k1[n1];
                int2dz_2[n1] = int2dz_k1[n1];
                int2dz_k1[n1] = int2dz_0[n1];
                int2dz[n + n1 + k * (shellp + 1) * ngqexq] = int2dz_0[n1];
            }
        }
    }


/*             ...evaluate I=1,SHELLP and K=1,SHELLQ (if any) */
/*                in most economical way. */


    if (shellq <= shellp)
    {
        for (n = 0; n < ngqexq; n+=SIMDW)
        {
            ERD_SIMD_ALIGN double int2dx_00[SIMDW], int2dx_10[SIMDW], int2dx_20[SIMDW], int2dx_11[SIMDW];
            ERD_SIMD_ALIGN double int2dy_00[SIMDW], int2dy_10[SIMDW], int2dy_20[SIMDW], int2dy_11[SIMDW];
            ERD_SIMD_ALIGN double int2dz_00[SIMDW], int2dz_10[SIMDW], int2dz_20[SIMDW], int2dz_11[SIMDW];

            for (k = 1; k <= shellq; ++k)
            {
                int k1 = k - 1;

                #pragma vector aligned
                #pragma simd
                for(n1 = 0; n1 < SIMDW; n1++)
                {
                    b0 = k * b00[n + n1];

                    int2dx_10[n1] = b0 * int2dx[n + n1 + k1 * (shellp + 1) * ngqexq] +
                        c00x[n + n1] * int2dx[n + n1 + k * (shellp + 1) * ngqexq];

                    int2dy_10[n1] = b0 * int2dy[n + n1 + k1 * (shellp + 1) * ngqexq] +
                        c00y[n + n1] * int2dy[n + n1 + k * (shellp + 1) * ngqexq];

                    int2dz_10[n1] = b0 * int2dz[n + n1 + k1 * (shellp + 1) * ngqexq] +
                        c00z[n + n1] * int2dz[n + n1 + k * (shellp + 1) * ngqexq];

                    int2dx_20[n1] = int2dx[n + n1 + (k * (shellp + 1)) * ngqexq];
                    int2dy_20[n1] = int2dy[n + n1 + (k * (shellp + 1)) * ngqexq];
                    int2dz_20[n1] = int2dz[n + n1 + (k * (shellp + 1)) * ngqexq];
                }

                #pragma vector aligned
                #pragma simd
                for(n1 = 0; n1 < SIMDW; n1++)
                {
                    int2dx[n + n1 + (k * (shellp + 1) + 1) * ngqexq] = int2dx_10[n1];
                    int2dy[n + n1 + (k * (shellp + 1) + 1) * ngqexq] = int2dy_10[n1];
                    int2dz[n + n1 + (k * (shellp + 1) + 1) * ngqexq] = int2dz_10[n1];
                }
                for (i = 2; i <= shellp; ++i)
                {
                    int i1 = i - 1;
                    #pragma vector aligned
                    #pragma simd
                    for(n1 = 0; n1 < SIMDW; n1++)
                    {
                        b0 = k * b00[n + n1];
                        b1 = i1 * b10[n + n1];
                        int2dx_11[n1] = int2dx[n + n1 + (i1 + k1 * (shellp + 1)) * ngqexq];
                        int2dx_00[n1] = b0 * int2dx_11[n1] + b1 * int2dx_20[n1] + c00x[n + n1] * int2dx_10[n1];
                        int2dx_20[n1] = int2dx_10[n1];
                        int2dx_10[n1] = int2dx_00[n1];

                        int2dy_11[n1] = int2dy[n + n1 + (i1 + k1 * (shellp + 1)) * ngqexq];
                        int2dy_00[n1] = b0 * int2dy_11[n1] + b1 * int2dy_20[n1] + c00y[n + n1] * int2dy_10[n1];
                        int2dy_20[n1] = int2dy_10[n1];
                        int2dy_10[n1] = int2dy_00[n1];

                        int2dz_11[n1] = int2dz[n + n1 + (i1 + k1 * (shellp + 1)) * ngqexq];
                        int2dz_00[n1] = b0 * int2dz_11[n1] + b1 * int2dz_20[n1] + c00z[n + n1] * int2dz_10[n1];
                        int2dz_20[n1] = int2dz_10[n1];
                        int2dz_10[n1] = int2dz_00[n1];
                    }
                    #pragma vector aligned
                    #pragma simd
                    for(n1 = 0; n1 < SIMDW; n1++)
                    {
                        int2dx[n + n1 + (i + k * (shellp + 1)) * ngqexq] = int2dx_00[n1];
                        int2dy[n + n1 + (i + k * (shellp + 1)) * ngqexq] = int2dy_00[n1];
                        int2dz[n + n1 + (i + k * (shellp + 1)) * ngqexq] = int2dz_00[n1];
                    }
                }
            }
        }
    }
    else
    {
        for (n = 0; n < ngqexq; n+=SIMDW)
        {
            ERD_SIMD_ALIGN double int2dx_00[SIMDW], int2dx_01[SIMDW], int2dx_02[SIMDW], int2dx_11[SIMDW];
            ERD_SIMD_ALIGN double int2dy_00[SIMDW], int2dy_01[SIMDW], int2dy_02[SIMDW], int2dy_11[SIMDW];
            ERD_SIMD_ALIGN double int2dz_00[SIMDW], int2dz_01[SIMDW], int2dz_02[SIMDW], int2dz_11[SIMDW];

            for (i = 1; i <= shellp; ++i)
            {
                int i1 = i - 1;

                #pragma vector aligned
                #pragma simd
                for(n1 = 0; n1 < SIMDW; n1++)
                {
                    b0 = i * b00[n + n1];
                    int2dx_01[n1] = b0 * int2dx[n + n1 + i1 * ngqexq]
                        + d00x[n + n1] * int2dx[n + n1 + i * ngqexq];

                    int2dy_01[n1] = b0 * int2dy[n + n1 + i1 * ngqexq]
                        + d00y[n + n1] * int2dy[n + n1 + i * ngqexq];

                    int2dz_01[n1] = b0 * int2dz[n + n1 + i1 * ngqexq]
                        + d00z[n + n1] * int2dz[n + n1 + i * ngqexq];

                    int2dx_02[n1] = int2dx[n + n1 + (i) * ngqexq];
                    int2dy_02[n1] = int2dy[n + n1 + (i) * ngqexq];
                    int2dz_02[n1] = int2dz[n + n1 + (i) * ngqexq];
                }
                #pragma vector aligned
                #pragma simd
                for(n1 = 0; n1 < SIMDW; n1++)
                {
                    int2dx[n + n1 + (i + (shellp + 1)) * ngqexq] = int2dx_01[n1];
                    int2dy[n + n1 + (i + (shellp + 1)) * ngqexq] = int2dy_01[n1];
                    int2dz[n + n1 + (i + (shellp + 1)) * ngqexq] = int2dz_01[n1];
                }

                for (k = 2; k <= shellq; ++k)
                {
                    int k1 = k - 1;
                    #pragma vector aligned
                    #pragma simd
                    for(n1 = 0; n1 < SIMDW; n1++)
                    {
                        b0 = i * b00[n + n1];
                        b1 = k1 * b01[n + n1];

                        int2dx_11[n1] = int2dx[n + n1 + (i1 + k1 * (shellp + 1)) * ngqexq];
                        int2dx_00[n1] = b0 * int2dx_11[n1] + b1 * int2dx_02[n1] + d00x[n + n1] * int2dx_01[n1];
                        int2dx_02[n1] = int2dx_01[n1];
                        int2dx_01[n1] = int2dx_00[n1];

                        int2dy_11[n1] = int2dy[n + n1 + (i1 + k1 * (shellp + 1)) * ngqexq];
                        int2dy_00[n1] = b0 * int2dy_11[n1] + b1 * int2dy_02[n1] + d00y[n + n1] * int2dy_01[n1];
                        int2dy_02[n1] = int2dy_01[n1];
                        int2dy_01[n1] = int2dy_00[n1];

                        int2dz_11[n1] = int2dz[n + n1 + (i1 + k1 * (shellp + 1)) * ngqexq];
                        int2dz_00[n1] = b0 * int2dz_11[n1] + b1 * int2dz_02[n1] + d00z[n + n1] * int2dz_01[n1];
                        int2dz_02[n1] = int2dz_01[n1];
                        int2dz_01[n1] = int2dz_00[n1];
                    }
                    #pragma vector aligned
                    #pragma simd
                    for(n1 = 0; n1 < SIMDW; n1++)
                    {
                        int2dx[n + n1 + (i + k * (shellp + 1)) * ngqexq] = int2dx_00[n1];
                        int2dy[n + n1 + (i + k * (shellp + 1)) * ngqexq] = int2dy_00[n1];
                        int2dz[n + n1 + (i + k * (shellp + 1)) * ngqexq] = int2dz_00[n1];
                    }
                }
            }
        }
    }

    return 0;
}


#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
