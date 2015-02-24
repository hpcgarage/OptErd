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

#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "erd.h"
#include "erdutil.h"

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__SET_ABCD */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This routine handles the logistics on how to evaluate */
/*                the (12|34) integral batch in the most efficient way. */
/*                It performs the label map: */
/*                             (12|34) --> (AB|CD) */
/*                according to certain criteria that have to be met for */
/*                efficiency. The freedom we have in making the internal */
/*                association 1,2,3,4 -> A,B,C,D follows from the 8-fold */
/*                permutational symmetry of the integrals in (12|34): */
/*                     (12|34) = (21|34) = (12|43) = (21|43) = */
/*                     (34|12) = (43|12) = (34|21) = (43|21) */
/*                where the first line has only 12 -> 21 and 34 -> 43 */
/*                switches and the second line is obtained from the */
/*                first by bra <-> ket transpositions. */
/*                The type of switch to be applied is simply governed */
/*                by the demand that the final A,B,C,D shell labels */
/*                obey the relations A>=B and C>=D, since this means */
/*                the least amount of work (# of steps) for the HRR */
/*                procedure. The decision to perform a bra <-> ket */
/*                transposition comes from handling memory allocations */
/*                during both HRR on the bra and ket sides. */
/*                  Input (x = 1,2,3 and 4): */
/*                   NCGTOx        =  # of contractions for csh x */
/*                   NPGTOx        =  # of primitives per contraction */
/*                                    for csh x */
/*                   SHELLx        =  the shell type for csh x */
/*                   Xy,Yy,Zy      =  the x,y,z-coordinates for centers */
/*                                    y = 1,2,3 and 4 */
/*                   EXPx          =  primitive exponents for csh x */
/*                   CCx           =  contraction coeffs for csh x */
/*                   SPHERIC       =  is true, if spherical integrals */
/*                                    are wanted, false if cartesian */
/*                                    ones are wanted */
/*                  Output (x = A,B,C and D): */
/*                   NCGTOx        =  # of contractions for csh x */
/*                   NPGTOx        =  # of primitives per contraction */
/*                                    for csh x */
/*                   SHELLx        =  the shell type for csh x */
/*                   SHELLy        =  the shell sums: y = P,Q,T = */
/*                                    A+B,C+D,P+Q */
/*                   MXSHELL       =  the largest (maximum) shell type */
/*                   Xy,Yy,Zy      =  the x,y,z-coordinates for centers */
/*                                    y = A,B,C and D */
/*                   ATOMIC        =  indicates, if purely atomic */
/*                                    integrals will be evaluated */
/*                   ATOMxy        =  indicates, if atoms x and y are */
/*                                    considered to be equal for the */
/*                                    pairs xy = AB and CD */
/*                   EQUALxy       =  indicates, if csh x and csh y are */
/*                                    considered to be equal for the */
/*                                    pairs xy = AB and CD */
/*                   xxX,xxY,xxZ   =  the x,y,z-coordinate differences */
/*                                    for centers xx = AB and CD */
/*                   NxxCOOR       =  # of non-zero x,y,z-coordinate */
/*                                    differences for centers xx = AB */
/*                                    and CD */
/*                   RNxxSQ        =  square of the magnitude of the */
/*                                    distance between centers xx = AB */
/*                                    and CD */
/*                   SPNORM        =  normalization factor due to */
/*                                    presence of s- and p-type shells. */
/*                                    For each s-type shell there is a */
/*                                    factor of 1 and for each p-type */
/*                                    shell a factor of 2 */
/*                   NXYZx         =  # of cartesian monomials for csh x */
/*                   NXYZE(F)T     =  sum of # of cartesian monomials */
/*                                    for all shells in the range */
/*                                    E = A,...,A+B and in the range */
/*                                    F = C,...,C+D */
/*                   NXYZy         =  # of cartesian monomials for */
/*                                    y = P,Q shells */
/*                   NRYx          =  # of spherical functions for csh x */
/*                   INDEXx        =  index A,B,C,D -> 1,2,3,4 map */
/*                   SWAPxy        =  is .true. for xy = 12 and 34, if */
/*                                    a swap 1 <-> 2 and 3 <-> 4 has */
/*                                    been performed */
/*                   SWAPRS(TU)    =  is set .true. if the contraction */
/*                                    order of the primitives pair AB(CD) */
/*                                    will be performed in reverse order */
/*                                    BA(DC) for efficiency reasons */
/*                   TR1234        =  is .true., if a bra <-> ket */
/*                                    transposition has been applied */
/*                   LEXPx         =  pointers to locate appropriate */
/*                                    section of the exponent array */
/*                                    corresponding to csh x */
/*                   LCCx          =  pointers to locate appropriate */
/*                                    section of the contraction coeff */
/*                                    array corresponding to csh x */
/*                   LCCSEGx       =  pointers to locate appropriate */
/*                                    section of the lowest and highest */
/*                                    primitive index array defining */
/*                                    segmented contraction boundaries */
/*                                    for csh x */
/*                   NXYZHRR       =  maximum dimension of cartesian and */
/*                                    spherical components during the */
/*                                    entire HRR contraction procedure */
/*                   NCOLHRR       =  maximum # of HRR rotation matrix */
/*                                    columns needed to generate the */
/*                                    final HRR rotation matrix */
/*                   NROTHRR       =  maximum # of HRR rotation matrix */
/*                                    elements needed to generate the */
/*                                    final HRR rotation matrix */
/*                   EMPTY         =  int flag, indicating if an */
/*                                    empty batch of integrals is */
/*                                    expected. */
/* ------------------------------------------------------------------------ */
ERD_OFFLOAD void erd__set_abcd(
    uint32_t A_ptr[restrict static 1], uint32_t B_ptr[restrict static 1], uint32_t C_ptr[restrict static 1], uint32_t D_ptr[restrict static 1],
    const uint32_t shell[restrict static 1], const double xyz0[restrict static 1],
    bool spheric,
    uint32_t indexa_ptr[restrict static 1], uint32_t indexb_ptr[restrict static 1], uint32_t indexc_ptr[restrict static 1], uint32_t indexd_ptr[restrict static 1],
    uint32_t *restrict nxyza_ptr, uint32_t *restrict nxyzb_ptr, uint32_t *restrict nxyzc_ptr, uint32_t *restrict nxyzd_ptr,
    uint32_t *restrict nxyzet_ptr, uint32_t *restrict nxyzft_ptr,
    uint32_t nrya_ptr[restrict static 1], uint32_t nryb_ptr[restrict static 1], uint32_t nryc_ptr[restrict static 1], uint32_t nryd_ptr[restrict static 1],
    uint32_t nabcoor_ptr[restrict static 1], uint32_t ncdcoor_ptr[restrict static 1],
    uint32_t ncolhrr[restrict static 1], uint32_t nrothrr[restrict static 1],
    uint32_t nxyzhrr[restrict static 1], bool empty[restrict static 1])
{
/*             ...generate all 1,2,3,4 data. Decide as early as */
/*                possible, if a zero batch of integrals is expected. */
    *empty = false;
    uint32_t A = *A_ptr, B = *B_ptr, C = *C_ptr, D = *D_ptr;
    uint32_t shell1 = shell[A], shell2 = shell[B], shell3 = shell[C], shell4 = shell[D];

    const uint32_t preshellp = shell1 + shell2;
    const uint32_t preshellq = shell3 + shell4;
    const uint32_t shellt = preshellp + preshellq;

    const uint32_t mxshell = max4x32u(shell1, shell2, shell3, shell4);
    const bool atomic = ((A ^ B) | (B ^ C) | (C ^ D)) == 0;
    const bool case1 = (shellt % 2 == 1);
    const bool case2 = spheric && (2 * mxshell > shellt);
    if (atomic && (case1 || case2)) {
        *empty = true;
        return;
    }

/*             ...determine csh equality between center pairs 1,2 */
/*                and 3,4 in increasing order of complexity: */
/*                 centers -> shells -> exponents -> ctr coefficients */
/*             ...set the cartesian and spherical dimensions. In case */
/*                no spherical transformations are wanted, set the */
/*                corresponding dimensions equal to the cartesian ones. */
    uint32_t nxyz[4] = {
        ((shell1 + 1) * (shell1 + 2)) / 2,
        ((shell2 + 1) * (shell2 + 2)) / 2,
        ((shell3 + 1) * (shell3 + 2)) / 2,
        ((shell4 + 1) * (shell4 + 2)) / 2
    };
    uint32_t nry[4] = {
        2 * shell1 + 1,
        2 * shell2 + 1,
        2 * shell3 + 1,
        2 * shell4 + 1
    };
    if (!spheric) {
        nry[0] = nxyz[0];
        nry[1] = nxyz[1];
        nry[2] = nxyz[2];
        nry[3] = nxyz[3];
    }

/*             ...decide on the 1 <-> 2 and/or 3 <-> 4 swapping. */

/*             ...calculate NXYZHRR for the two possible HRR */
/*                and (if any) cartesian -> spherical transformation */
/*                sequences: */
/*                    i) initial dimension:  NXYZE * NXYZF */
/*                   ii) perform HRR on 34:  NXYZE * NXYZ3 * NXYZ4 */
/*                  iii) cart -> sph on 34:  NXYZE * NRY3 * NRY4 */
/*                   iv) perform HRR on 12:  NXYZ1 * NXYZ2 * NRY3 * NRY4 */
/*                    v) cart -> sph on 12:  NRY1 * NRY2 * NRY3 * NRY4 */
/*                    i) initial dimension:  NXYZE * NXYZF */
/*                   ii) perform HRR on 12:  NXYZ1 * NXYZ2 * NXYZF */
/*                  iii) cart -> sph on 12:  NRY1 * NRY2 * NXYZF */
/*                   iv) perform HRR on 34:  NRY1 * NRY2 * NXYZ3 * NXYZ4 */
/*                    v) cart -> sph on 34:  NRY1 * NRY2 * NRY3 * NRY4 */
/*                The only dimension increasing steps are ii) and iv) */
/*                in both cases. Hence we first find the maximum */
/*                between ii) and iv) for both sequences and then we */
/*                take the overall minimum of these two maxima. */
/*                Since the order of sequence of the HRR on the A,B,C,D */
/*                labels is CD followed by AB, the overall minimum */
/*                will decide if to perform a bra <-> ket transposition */
/*                on the 12/34 labels. */
    const uint32_t preshella = max32u(shell1, shell2);
    const uint32_t preshellb = min32u(shell1, shell2);
    const uint32_t preshellc = max32u(shell3, shell4);
    const uint32_t preshelld = min32u(shell3, shell4);
    const uint32_t nxyze = (preshellp + 1) * (preshellp + 2) * (preshellp + 3) / 6 - preshella * (preshella + 1) * (preshella + 2) / 6;
    const uint32_t nxyzf = (preshellq + 1) * (preshellq + 2) * (preshellq + 3) / 6 - preshellc * (preshellc + 1) * (preshellc + 2) / 6;

    bool transpose = false;
    if ((preshellb | preshelld) == 0) {
        *nxyzhrr = nxyze * nxyzf;
    } else {
        const uint32_t nhrr1st = max32u(nxyze * nxyz[2] * nxyz[3], nxyz[0] * nxyz[1] * nry[2] * nry[3]);
        const uint32_t nhrr2nd = max32u(nxyzf * nxyz[0] * nxyz[1], nxyz[2] * nxyz[3] * nry[0] * nry[1]);
        *nxyzhrr = min32u(nhrr1st, nhrr2nd);
        transpose = nhrr1st > nhrr2nd;
    }

/*             ...according to the previously gathered info, set the */
/*                new A,B,C,D shells, # of primitives + contraction */
/*                coeffs as well as pointers to the alpha exponents */
/*                and contraction coefficients. Also set the info for */
/*                evaluation of the [e0|f0] batches and for the HRR */
/*                steps later on. */
    uint32_t indexa = 0, indexb = 1, indexc = 2, indexd = 3;
    uint32_t nxyzet = nxyze, nxyzft = nxyzf;
    if (transpose) {
        ERD_SWAP(nxyzet, nxyzft);

        ERD_SWAP(A, C);
        ERD_SWAP(B, D);
        
        indexa = 2;
        indexb = 3;
        indexc = 0;
        indexd = 1;
    }
    *nxyzet_ptr = nxyzet;
    *nxyzft_ptr = nxyzft;
    
    if (shell[A] < shell[B]) {
        ERD_SWAP(A, B);
        ERD_SWAP(indexa, indexb);
    }
    if (shell[C] < shell[D]) {
        ERD_SWAP(C, D);
        ERD_SWAP(indexc, indexd);
    }
    const double xa = xyz0[A*4], xb = xyz0[B*4], xc = xyz0[C*4], xd = xyz0[D*4];
    const double ya = xyz0[A*4+1], yb = xyz0[B*4+1], yc = xyz0[C*4+1], yd = xyz0[D*4+1];
    const double za = xyz0[A*4+2], zb = xyz0[B*4+2], zc = xyz0[C*4+2], zd = xyz0[D*4+2];
    *A_ptr = A;
    *B_ptr = B;
    *C_ptr = C;
    *D_ptr = D;
    *indexa_ptr = indexa;
    *indexb_ptr = indexb;
    *indexc_ptr = indexc;
    *indexd_ptr = indexd;
    *nxyza_ptr = nxyz[indexa];
    *nxyzb_ptr = nxyz[indexb];
    *nxyzc_ptr = nxyz[indexc];
    *nxyzd_ptr = nxyz[indexd];
    *nrya_ptr = nry[indexa];
    *nryb_ptr = nry[indexb];
    *nryc_ptr = nry[indexc];
    *nryd_ptr = nry[indexd];
    const uint8_t nabcoor = (uint8_t)(xa != xb) + (uint8_t)(ya != yb) + (uint8_t)(za != zb);
    const uint8_t ncdcoor = (uint8_t)(xc != xd) + (uint8_t)(yc != yd) + (uint8_t)(zc != zd);
    const uint32_t shella = shell[A], shellb = shell[B], shellc = shell[C], shelld = shell[D];
    const uint32_t shellp = shella + shellb;
    const uint32_t shellq = shellc + shelld;
    const uint32_t nxyzp = (shellp + 1) * (shellp + 2) / 2;
    const uint32_t nxyzq = (shellq + 1) * (shellq + 2) / 2;
    
    *ncolhrr = 0;
    *nrothrr = 0;
    if (shellb != 0) {
        const uint32_t ngh = nxyzet;
        uint32_t nxyzgo = nxyzet;
        uint32_t nxyzho = 1;
        uint32_t nxyzi = nxyzp;
        uint32_t shellg = shellp;
        uint32_t nrow = 1;
        uint32_t ncol = ngh;
        uint32_t nrot = ngh;
        for (uint32_t shellh = 1; shellh <= shellb; shellh++) {
            nxyzgo -= nxyzi;
            nxyzho += shellh + 1;
            const uint32_t ngho = nxyzgo * nxyzho;
            switch (nabcoor) {
                case 3:
                {
                    const uint32_t m = shellh / 3 + 1;
                    nrow += m * m;
                    if (shellh % 3 == 2)
                        nrow += m;
                    break;
                }
                case 2:
                    nrow += shellh / 2 + 1;
                    break;
                case 1:
                    ++nrow;
                    break;
            }
            ncol = max32u(ngho, ncol);
            nrot = max32u(nrow * ngho, nrot);
            nxyzi = nxyzi - shellg - 1;
            --shellg;
        }
        *ncolhrr = ncol;
        *nrothrr = nrot;
    }
/*             ...next find maximum values for the HRR on the CD-part */
/*                and set overall maximum values. */
    if (shelld != 0) {
        const uint32_t ngh = nxyzft;
        uint32_t nxyzgo = nxyzft;
        uint32_t nxyzho = 1;
        uint32_t nxyzi = nxyzq;
        uint32_t shellg = shellq;
        uint32_t nrow = 1;
        uint32_t ncol = ngh;
        uint32_t nrot = ngh;
        for (uint32_t shellh = 1; shellh <= shelld; shellh++) {
            nxyzgo -= nxyzi;
            nxyzho += shellh + 1;
            const uint32_t ngho = nxyzgo * nxyzho;
            switch (ncdcoor) {
                case 3:
                {
                    const uint32_t m = shellh / 3 + 1;
                    nrow += m * m;
                    if (shellh % 3 == 2)
                        nrow += m;
                    break;
                }
                case 2:
                    nrow += shellh / 2 + 1;
                    break;
                case 1:
                    nrow += 1;
                    break;
            }
            ncol = max32u(ngho, ncol);
            nrot = max32u(nrow * ngho, nrot);
            nxyzi -= shellg + 1;
            --shellg;
        }
        *ncolhrr = max32u(ncol, *ncolhrr);
        *nrothrr = max32u(nrot, *nrothrr);
    }
    *nabcoor_ptr = nabcoor;
    *ncdcoor_ptr = ncdcoor;
}
