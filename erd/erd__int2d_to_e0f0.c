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

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif


/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__INT2D_TO_E0F0 */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This routine assembles the set of batches of cartesian */
/*                eris [E0|F0] , E = A to P, F = C to Q, adding up all */
/*                the contributions from all the 2D PQ integrals. */
/*                The routine uses the reduced Rys multiplication scheme */
/*                as suggested in R.Lindh, U.Ryu and B.Liu, JCP 95, 5889. */
/*                This scheme reuses intermediate products between */
/*                2DX and 2DY integrals, which can be achieved by having */
/*                the outer loops run over all possible x and y monomial */
/*                parts and the inner loops over all allowed E and F */
/*                shell combinations. The price to pay for such loop */
/*                ordering is the scattered addressing of locations */
/*                within the batch array, which has its rows and columns */
/*                ordered such that the E and F shells are increasing */
/*                and within each E and F shell the monomials are */
/*                ordered such that x>y>z in exponents. */
/*                An example follows: */
/*                ------------------- */
/*                Let E = 0,2 and F = 0,1. Then we have the left and */
/*                right hand of the batch array ordered as follows: */
/*                          left xyz        right xyz */
/*                             000             000 */
/*                             ---             --- */
/*                             100             100 */
/*                             010             010 */
/*                             001             001 */
/*                             --- */
/*                             200 -> -5 */
/*                             110 */
/*                             101 -> -3 */
/*                             020 */
/*                             011 */
/*                             002 -> 0 */
/*                The batch would thus have dimensions 10 x 4. For */
/*                the left side (and analogous for the right side) the */
/*                reduced multiplication scheme would have its outer */
/*                most loop run over x=0,2, followed by the next loop */
/*                y=0,2-x. The innermost loop would then run over the */
/*                allowed shells E=E(max),max(E(min),x+y). In this */
/*                case all x,y-pairs can be reused for all appropriate */
/*                shell combinations. */
/*                To find the address of a specific x,y,z,E combination */
/*                inside the batch array, we first note that the z-part */
/*                is dependent on the x,y-parts and is hence not needed. */
/*                Lets look at the E-part first. The E-part is evaluated */
/*                from its dimension formula (E+1)*(E+2)/2. Organizing */
/*                the inner E-loop to run from E(max) always, the */
/*                dimension for E(max) is passed as the argument NXYZP */
/*                and all lower E dimensions are calculated by the */
/*                formula relating dimensions between E and E+1: */
/*                           dim(E+1) = dim(E) + E + 2 */
/*                In this way multiplications in the E-part can be */
/*                entirely avoided. The x,y-part is defined as the */
/*                part which has to be subtracted from dim(E) to */
/*                reach the xyz monomial position inside the E shell. */
/*                It can be divided into an x-part and a y-part. The */
/*                x-part is given by the fomula: */
/*                          x-part = - x*E + x(x-3)/2 */
/*                and for the example above has been given for E=2 */
/*                and marked with arrows ->. The last term of the x-part */
/*                involves 1 multiplication and division, however it */
/*                can be changed to: */
/*                                               x-1 */
/*                          x-part = - x*E - x + sum i */
/*                                               i=0 */
/*                and clever additions inside the x-loop avoid the use */
/*                of multiplications and divisions. The y-part is trivial */
/*                and is simply equal to -y. The overall conclusion is */
/*                thus that the location of a specific x,y,z,E quadruple */
/*                inside the batch comes at the cost of one x*E(max) */
/*                multiplication in the outermost x-loops, since the */
/*                other x*E ones can again be reached via stepwise */
/*                subtraction of x from x*E(max). */
/*                Due to the very computational intensive steps inside */
/*                the x,y,z,E loops, special sections of identical */
/*                x,y,z,E loop structers have been given for each */
/*                # of roots =< 9, thus saving considerable computing */
/*                time over the general case. */
/*                For comments on how the x,y,z,E loop structures are */
/*                coded please refer to the general root case. */
/*                  Input: */
/*                    SHELLx      =  shell types for individual csh */
/*                                   x=A,C and csh sums P=A+B,Q=C+D */
/*                    NGQP        =  # of gaussian quadrature points */
/*                                   (roots) */
/*                    NEXQ        =  current # of exponent quadruplets */
/*                    NGQEXQ      =  product of # of gaussian quadrature */
/*                                   points times exponent quadruplets */
/*                    NXYZE(F)T   =  sum of # of cartesian monomials */
/*                                   for all shells in the range */
/*                                   E = A,...,P=A+B and in the range */
/*                                   F = C,...,Q=C+D */
/*                    NXYZy       =  # of cartesian monomials for */
/*                                   y = P,Q shells */
/*                    INT2Dx      =  all current 2D PQ integrals for */
/*                                   each cartesian component */
/*                                   (x = X,Y,Z) */
/*                    TEMP1(2)    =  scratch arrays holding intermediate */
/*                                   2D PQ integral products */
/*                    SCALE       =  the NGQEXQ scaling factors */
/*                  Output: */
/*                    BATCH       =  batch of primitive cartesian */
/*                                   [E0|F0] integrals corresponding */
/*                                   to all current exponent quadruplets */
/* ------------------------------------------------------------------------ */

int erd__int2d_to_e0f0 (int shella, int shellp, int shellc, int shellq,
                        int ngqexq, int nxyzet, int nxyzft,
                        double *int2dx, double *int2dy, double *int2dz,
                        int **vrrtab, double *batch)
{
    uint64_t m;
    int *tabe;
    int *tabf;
    int ke;
    int kf;
    uint64_t indx;
    uint64_t indy;
    uint64_t indz;
    uint64_t indb;
               
    tabe = vrrtab[shella];
    tabf = vrrtab[shellc];
    indb = 0;
#if defined (__MIC__)
    __m512i zero512 = _mm512_setzero_epi32();
    int indxyz[16] __attribute__((aligned(64)));
    __m512i ngqexq512 = _mm512_set1_epi32(ngqexq);
    __m512i shellp_plus_one512 = _mm512_set1_epi32(shellp + 1);
    __m512i indxyz512;

    for (kf = 0; kf < nxyzft; kf++)
    {
        __m512i xyzf512 = _mm512_extloadunpacklo_epi32(zero512, &tabf[kf * 4], _MM_UPCONV_EPI32_NONE, _MM_HINT_NONE);
        for(ke = 0; ke < nxyzet; ke++)
        {
            __m512i xyze512 = _mm512_extloadunpacklo_epi32(zero512, &tabe[ke * 4], _MM_UPCONV_EPI32_NONE, _MM_HINT_NONE);
            indxyz512 = _mm512_fmadd_epi32(xyzf512, shellp_plus_one512, xyze512);
            indxyz512 = _mm512_mullo_epi32(indxyz512, ngqexq512);
            _mm512_store_epi32(indxyz, indxyz512);
            indx = indxyz[0];
            indy = indxyz[1];
            indz = indxyz[2];
            //indb = kf * nxyzet + ke;
            __m512d sum512 = _mm512_setzero_pd();
            for(m = 0; m < ngqexq; m+=SIMDW)
            {
                __m512d int2dx512 = _mm512_load_pd(&int2dx[indx + m]);
                __m512d int2dy512 = _mm512_load_pd(&int2dy[indy + m]);
                __m512d int2dz512 = _mm512_load_pd(&int2dz[indz + m]);

                __m512d temp512 = _mm512_mul_pd(int2dx512, int2dy512);
                sum512 = _mm512_fmadd_pd(temp512, int2dz512, sum512);
            }
            batch[indb++] = _mm512_reduce_add_pd(sum512);
        }
    }
#else
    int xe, ye, ze, xf, yf, zf;   
    int m1;
    double sum = 0;

    for (kf = 0; kf < nxyzft; kf++)
    {
        xf = tabf[kf * 4 + 0];
        yf = tabf[kf * 4 + 1];
        zf = tabf[kf * 4 + 2];
        for(ke = 0; ke < nxyzet; ke++)
        {
            xe = tabe[ke * 4 + 0];
            ye = tabe[ke * 4 + 1];
            ze = tabe[ke * 4 + 2];
            indx = (xe + xf * (shellp + 1)) * ngqexq;
            indy = (ye + yf * (shellp + 1)) * ngqexq;
            indz = (ze + zf * (shellp + 1)) * ngqexq;
            //indb = kf * nxyzet + ke;
            sum = 0.0;
            for(m = 0; m < ngqexq; m+=SIMDW)
            {
#pragma vector aligned
                for(m1 = 0; m1 < SIMDW; m1++)
                {
                    sum += int2dx[m + m1 + indx]
                        * int2dy[m + m1 + indy]
                        * int2dz[m + m1 + indz];
                }
            }
            batch[indb++] = sum;
        }
    }
#endif
    return 0;
}


#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
