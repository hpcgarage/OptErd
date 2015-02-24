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

#include "erd.h"


#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__INT2D_TO_E000 */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This routine assembles the set of batches of cartesian */
/*                eris: */
/*                        [E0|00] or [00|E0] , E = A to P */
/*                adding up all the contributions from all the */
/*                respective 2D integrals: */
/*                                   P0 or 0P */
/*                Simplified version of the general E0F0 routine to */
/*                reduce loop overheads for those cases where there is */
/*                at least one s-shell on the bra or ket side. For */
/*                comments and details see the general E0F0 routine. */
/*                  Input: */
/*                    SHELLx      =  shell types for individual csh */
/*                                   x=A and csh sum P=A+B */
/*                    NGQP        =  # of gaussian quadrature points */
/*                                   (roots) */
/*                    NEXQ        =  current # of exponent quadruplets */
/*                    NGQEXQ      =  product of # of gaussian quadrature */
/*                                   points times exponent quadruplets */
/*                    NXYZET      =  sum of # of cartesian monomials */
/*                                   for all shells in the range */
/*                                   E = A,...,P=A+B */
/*                    NXYZP       =  # of cartesian monomials for the */
/*                                   P=A+B shell */
/*                    INT2Dx      =  all current 2D P0/0P integrals for */
/*                                   each cartesian component */
/*                                   (x = X,Y,Z) */
/*                    TEMP1(2)    =  scratch arrays holding intermediate */
/*                                   2D P0/0P integral products */
/*                    SCALE       =  the NGQEXQ scaling factors */
/*                  Output: */
/*                    BATCH       =  batch of primitive cartesian */
/*                                   [E0|00] integrals corresponding */
/*                                   to all current exponent quadruplets */
/* ------------------------------------------------------------------------ */
int erd__int2d_to_e000 (int shella, int shellp, int ngqp, int nexq, int ngqexq,
                        int nxyzet, int nxyzp,
                        double *int2dx, double *int2dy, double *int2dz,
                        double *temp1, double *temp2,
                        double *scale, double *batch)
{
    int i, k, m, n, se, xe, ye, ze, xep;
    double sum;
    int xye, xyep, seend, yeend, xemax, nxyze;

    xep = nxyzet + 2;
    for (xe = 0; xe <= shellp; ++xe)
    {
        xep = xep + xe - 2;
        xemax = xe * shellp;
        yeend = shellp - xe;
        for (m = 0; m < ngqexq; ++m)
        {
            temp1[m] = scale[m] * int2dx[m + xe * ngqexq];
        }
/*             ...middle loops over y-contributions. Skip multiplication */
/*                of y-contributions, if we have a 0-type, as then the */
/*                2DY integral is equal to 1. */
        xyep = xep - xemax;
        for (ye = 0; ye <= yeend; ++ye)
        {
            xye = xe + ye;
            --xyep;
            seend = MAX (shella, xye);
            if (ye == 0)
            {
                for (n = 0; n < ngqexq; ++n)
                {
                    temp2[n] = temp1[n];
                }
            }
            else
            {
                for (n = 0; n < ngqexq; ++n)
                {
                    temp2[n] = temp1[n] * int2dy[n + ye * ngqexq];
                }
            }
/*             ...inner loops over E-pairs. Skip multiplication */
/*                of z-contributions, if we have a 0-type, as */
/*                then the 2DZ integral is equal to 1. */
            i = xyep;
            nxyze = nxyzp;
            for (se = shellp; se >= seend; --se)
            {
                ze = se - xye;
/*             ...all info concerning all x-, y- and z-contributions */
/*                have been collected for all exponent quadruplets at */
/*                once. Sum up the 2D X,Y,Z integral products to the */
/*                appropriate place of the batch. */
                batch[i] = 0.0;
                if (ze == 0)
                {
                    k = 0;
                    for (m = 0; m < nexq; ++m)
                    {
                        sum = 0.;
                        for (n = 0; n < ngqp; ++n)
                        {
                            sum += temp2[k + n];
                        }
                        k += ngqp;
                        batch[i] += sum;
                    }
                }
                else
                {
                    k = 0;
                    for (m = 0; m < nexq; ++m)
                    {
                        sum = 0.;
                        for (n = 0; n < ngqp; ++n)
                        {
                            sum += temp2[k + n] *
                                int2dz[k + n + ze * ngqexq];
                        }
                        k += ngqp;
                        batch[i] += sum;
                    }
                }
/*             ...next z-contribution. */
                i = i - nxyze + xe;
                nxyze = nxyze - se - 1;
            }
/*             ...next y- and x-contribution. */
        }
    }

    return 0;
} 

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
