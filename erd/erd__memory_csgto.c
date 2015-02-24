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

#include "erd.h"
#include "erdutil.h"

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__MEMORY_CSGTO */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : ERD__SET_ABCD */
/*                ERD__E0F0_DEF_BLOCKS */
/*  DESCRIPTION : This operation calculates the minimum and optimum */
/*                int/flp memory needed for evaluating a batch */
/*                of contracted electron repulsion integrals on up to */
/*                four different centers between cartesian or spherical */
/*                gaussian type shells. */
/*                  Input (x = 1,2,3 and 4): */
/*                    NALPHA       =  total # of exponents */
/*                    NCOEFF       =  total # of contraction coeffs */
/*                    NCGTOx       =  # of contractions for csh x */
/*                    NPGTOx       =  # of primitives per contraction */
/*                                    for csh x */
/*                    SHELLx       =  the shell type for csh x */
/*                    Xy,Yy,Zy     =  the x,y,z-coordinates for centers */
/*                                    y = 1,2,3 and 4 */
/*                    ALPHA        =  primitive exponents for csh */
/*                                    1,2,3,4 in that order */
/*                    CC           =  contraction coefficient for csh */
/*                                    1,2,3,4 in that order, for each */
/*                                    csh individually such that an */
/*                                    (I,J) element corresponds to the */
/*                                    I-th primitive and J-th contraction. */
/*                    L1CACHE      =  Size of level 1 cache in units of */
/*                                    8 Byte */
/*                    NCTROW       =  minimum # of rows that are */
/*                                    accepted for blocked contractions */
/*                    SPHERIC      =  is true, if spherical integrals */
/*                                    are wanted, false if cartesian */
/*                                    ones are wanted */
/*                  Output: */
/*                    IMIN,IOPT    =  minimum/optimum int memory */
/*                    ZMIN,ZOPT    =  minimum/optimum flp memory */
/* ------------------------------------------------------------------------ */
ERD_OFFLOAD size_t erd__memory_csgto(uint32_t npgto1, uint32_t npgto2, uint32_t npgto3, uint32_t npgto4,
    uint32_t shell1, uint32_t shell2, uint32_t shell3, uint32_t shell4,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double x3, double y3, double z3,
    double x4, double y4, double z4,
    bool spheric)
{
    const uint32_t nxyz1 = ((shell1 + 1) * (shell1 + 2)) / 2;
    const uint32_t nxyz2 = ((shell2 + 1) * (shell2 + 2)) / 2;
    const uint32_t nxyz3 = ((shell3 + 1) * (shell3 + 2)) / 2;
    const uint32_t nxyz4 = ((shell4 + 1) * (shell4 + 2)) / 2;
    uint32_t nry1 = 2 * shell1 + 1;
    uint32_t nry2 = 2 * shell2 + 1;
    uint32_t nry3 = 2 * shell3 + 1;
    uint32_t nry4 = 2 * shell4 + 1;
    if (!spheric) {
        nry1 = nxyz1;
        nry2 = nxyz2;
        nry3 = nxyz3;
        nry4 = nxyz4;
    }

    const uint32_t shellp = shell1 + shell2;
    const uint32_t shellq = shell3 + shell4;
    const uint32_t shella = max32u(shell1, shell2);
    const uint32_t shellb = min32u(shell1, shell2);
    const uint32_t shellc = max32u(shell3, shell4);
    const uint32_t shelld = min32u(shell3, shell4);
    const uint32_t nxyze = (shellp + 1) * (shellp + 2) * (shellp + 3) / 6 - shella * (shella + 1) * (shella + 2) / 6;
    const uint32_t nxyzf = (shellq + 1) * (shellq + 2) * (shellq + 3) / 6 - shellc * (shellc + 1) * (shellc + 2) / 6;
    const uint32_t nxyzt = nxyze * nxyzf;

    uint32_t nxyzhrr;
    if ((shellb | shelld) == 0) {
        nxyzhrr = nxyze * nxyzf;
    } else {
        const uint32_t nhrr1st = max32u(nxyze * nxyz3 * nxyz4, nxyz1 * nxyz2 * nry3 * nry4);
        const uint32_t nhrr2nd = max32u(nxyzf * nxyz1 * nxyz2, nxyz3 * nxyz4 * nry1 * nry2);
        nxyzhrr = min32u(nhrr1st, nhrr2nd);
    }

    return PAD_LEN(nxyzt) + 2 * nxyzhrr;
}
