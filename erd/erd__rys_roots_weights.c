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
#include <string.h>
#include <assert.h>
#include <math.h>

#include "erd.h"
#include "erdutil.h"
#include "boys.h"

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__RYS_ROOTS_WEIGHTS */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : ERD__RYS_1_ROOTS_WEIGHTS */
/*                ERD__RYS_2_ROOTS_WEIGHTS */
/*                ERD__RYS_3_ROOTS_WEIGHTS */
/*                ERD__RYS_4_ROOTS_WEIGHTS */
/*                ERD__RYS_5_ROOTS_WEIGHTS */
/*                ERD__RYS_X_ROOTS_WEIGHTS */
/*  DESCRIPTION : This routine calculates NGQP-point Gaussian quadrature */
/*                rules on [0,1] over the Rys weight functions: */
/*                                       exp(-T*x) */
/*                             W   (x) = --------- */
/*                              Rys      2*sqrt(x) */
/*                for a set of NT T-exponents. Special interpolation */
/*                routines are provided for low number of roots and */
/*                weigths (NGQP < 6). On exit, NT x NGQP = NTGQP roots */
/*                and weights have been produced. */
/*                  Input: */
/*                    NT           =  # of T-exponents */
/*                    NTGQP        =  # of roots times # of T-exponents */
/*                    NGQP         =  # of gaussian quadrature points */
/*                                    (roots) */
/*                    NMOM         =  # of necessary moment integrals */
/*                                    to calculate the quadrature roots */
/*                    TVAL         =  the T-exponents */
/*                    RYSZERO      =  will hold the zeroth Rys moments */
/*                                    for all T-exponents */
/*                    FTABLE       =  Fm (T) table for interpolation */
/*                                    in low T region */
/*                    MGRID        =  maximum m in Fm (T) table */
/*                    NGRID        =  # of T's for which Fm (T) table */
/*                                    was set up */
/*                    TMAX         =  maximum T in Fm (T) table */
/*                    TSTEP        =  difference between two consecutive */
/*                                    T's in Fm (T) table */
/*                    TVSTEP       =  Inverse of TSTEP */
/*                    A,B          =  will contain the recurrence */
/*                                    coefficients for the auxillary */
/*                                    polynomials */
/*                    MOM          =  will contain the normed auxillary */
/*                                    polynomial modified moments */
/*                    DIA,OFF      =  will contain the diagonal and */
/*                                    offdiagonal elements of the */
/*                                    tridiagonal symmetric terminal */
/*                                    matrix */
/*                    ROW1,ROW2    =  first,second row intermediates. */
/*                                    Will be used to evaluate the */
/*                                    tridiagonal elements of the */
/*                                    symmetric terminal matrix in an */
/*                                    efficient way using Sack and */
/*                                    Donovan's method */
/*                  Output: */
/*                    RTS          =  the roots array */
/*                    WTS          =  the weights array */
/* ------------------------------------------------------------------------ */
ERD_OFFLOAD void erd__rys_roots_weights(uint32_t nt, uint32_t ngqp, uint32_t nmom,
                            const double tval[restrict],
                            double rts[restrict], double wts[restrict])
{
    switch (ngqp) {
#if 1
        case 1:
            erd__rys_1_roots_weights(nt, tval, rts, wts);
            return;
        case 2:
            erd__rys_2_roots_weights(nt, tval, rts, wts);
            return;
        case 3:
            erd__rys_3_roots_weights(nt, tval, rts, wts);
            return;
        case 4:
            erd__rys_4_roots_weights(nt, tval, rts, wts);
            return;
        case 5:
            erd__rys_5_roots_weights(nt, tval, rts, wts);
            return;
#endif
        default:
        {
            ERD_SIMD_ALIGN double ryszero[nt];
            /* ...# of roots and weights >= 6. Accumulate all zeroth Rys moments and call the general routine. */
            for (int n = 0; n < nt; n++) {
                const double t = tval[n];
                if (t == 0.0) {
                    ryszero[n] = 1.0;
                } else if (t <= tmax) {
                    const int tgrid = lround(t * tvstep);
                    const double delta = tgrid * tstep - t;
                    ryszero[n] = (((((boys_table[tgrid][6] * delta * 0.166666666666667 +
                                      boys_table[tgrid][5]) * delta * 0.2 +
                                      boys_table[tgrid][4]) * delta * 0.25 +
                                      boys_table[tgrid][3]) * delta * 0.333333333333333 +
                                      boys_table[tgrid][2]) * delta * 0.5 +
                                      boys_table[tgrid][1]) * delta +
                                      boys_table[tgrid][0];
                } else {
                    ryszero[n] = sqrt (3.141592653589793 / t) * .5;
                }
            }
            int ntgqp = nt * ngqp;
            erd__rys_x_roots_weights(nt, ntgqp, ngqp, nmom, tval, ryszero, rts, wts);
            return;
        }
    }

}
