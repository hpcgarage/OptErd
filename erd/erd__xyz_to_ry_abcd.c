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
#include "erdutil.h"

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__XYZ_TO_RY_ABCD */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : ERD__XYZ_TO_RY_MATRIX */
/*  DESCRIPTION : This operation generates all the info needed for */
/*                cartesian -> spherical transformation of a batch */
/*                of cartesian integrals (AB|CD). No duplicate info */
/*                is generated, i.e. the above transmitted pointers */
/*                of the type Z00x, I0x1 and I0x2 with x=A,B,C,D may */
/*                coincide if the shell types among the A,B,C,D are */
/*                equal (see below). The transformation data is placed */
/*                in two arrays (int and flp) at the appropriate */
/*                locations. */
/*                Input: */
/*                    NXYZx = dimension of xyz-monomial basis */
/*                            corresponding to the shell SHELLx */
/*                            for x=A,B,C,D. */
/*                     NRYx = dimension of ry-spherical basis */
/*                            corresponding to the shell SHELLx */
/*                            for x=A,B,C,D. */
/*                   SHELLx = the shell types for A,B,C,D. */
/*                I(Z)START = Starting location for the int (flp) */
/*                            data. */
/*                Output: */
/*                    NROWx = maximum # of xyz-monomials contributing */
/*                            to the ry-components for x=A,B,C,D. */
/*                    NROTx = maximum # of elements in transformation */
/*                            matrix for x=A,B,C,D. This is equal to */
/*                            NROWx times NRYx. */
/*                     Z00x = pointer for the transformation matrix */
/*                            elements for x=A,B,C,D. */
/*                     I0x1 = pointer for the # of row indices leading */
/*                            to non-zero elements in the transformation */
/*                            matrix for x=A,B,C,D. */
/*                     I0x2 = pointer for the row index labels of the */
/*                            non-zero elements in the transformation */
/*                            matrix for x=A,B,C,D. */
/*                 I(Z)USED = # of int (flp) words used. */
/*                 I(Z)CORE = The int (flp) arrays holding the */
/*                            transformation data. */
/*                If any of the A,B,C,D shell labels are equal, their */
/*                offsets will be set equal in the following sequence, */
/*                governed by their order of usage: */
/*                    1) If shell C = D, then: */
/*                              Z00C = Z00D */
/*                              I0C1 = I0D1 */
/*                              I0C2 = I0D2 */
/*                    2) If shell B = C or D, then: */
/*                              Z00B = Z00C or Z00D */
/*                              I0B1 = I0C1 or I0D1 */
/*                              I0B2 = I0C2 or I0D2 */
/*                    3) If shell A = B or C or D, then: */
/*                              Z00A = Z00B or Z00C or Z00D */
/*                              I0A1 = I0B1 or I0C1 or I0D1 */
/*                              I0A2 = I0B2 or I0C2 or I0D2 */
/*                Only mutually different transformation matrices + */
/*                associated data are generated. */
/* ------------------------------------------------------------------------ */
ERD_OFFLOAD void erd__xyz_to_ry_abcd(uint32_t nxyza, uint32_t nxyzb, uint32_t nxyzc, uint32_t nxyzd,
    uint32_t nrya, uint32_t nryb, uint32_t nryc, uint32_t nryd,
    uint32_t shella, uint32_t shellb, uint32_t shellc, uint32_t shelld,
    uint32_t nrowa[restrict static 1], uint32_t nrowb[restrict static 1], uint32_t nrowc[restrict static 1], uint32_t nrowd[restrict static 1],
    uint32_t nrota[restrict static 1], uint32_t nrotb[restrict static 1], uint32_t nrotc[restrict static 1], uint32_t nrotd[restrict static 1],
    uint32_t z00a[restrict static 1], uint32_t z00b[restrict static 1], uint32_t z00c[restrict static 1], uint32_t z00d[restrict static 1],
    uint32_t i0a1[restrict static 1], uint32_t i0b1[restrict static 1], uint32_t i0c1[restrict static 1], uint32_t i0d1[restrict static 1],
    uint32_t i0a2[restrict static 1], uint32_t i0b2[restrict static 1], uint32_t i0c2[restrict static 1], uint32_t i0d2[restrict static 1],
    uint32_t icore[restrict static 1], double zcore[restrict static 1])
{
    /* ...shell D data. */
    if (shelld > 1) {
        *nrowd = ((shelld / 2 + 1) * (shelld / 2 + 2)) / 2;
        *nrotd = *nrowd * nryd;
        *z00d = 0;
        const uint32_t z0dt = *z00d + *nrotd;
        *i0d1 = 0;
        *i0d2 = *i0d1 + nryd;
        erd__xyz_to_ry_matrix(nxyzd, *nrowd, shelld, &icore[*i0d1], &icore[*i0d2], &zcore[*z00d]);
    } else {
        *nrowd = 0;
        *nrotd = 0;
        *z00d = 0;
        *i0d2 = 0;
    }

    /* ...shell C data. */
    bool cdata = false;
    if (shellc > 1) {
        if (shellc == shelld) {
            *z00c = *z00d;
            *i0c1 = *i0d1;
            *i0c2 = *i0d2;
            *nrowc = *nrowd;
            *nrotc = *nrotd;
        } else {
            *nrowc = ((shellc / 2 + 1) * (shellc / 2 + 2)) / 2;
            *nrotc = *nrowc * nryc;
            *z00c = *z00d + *nrotd;
            const uint32_t z0ct = *z00c + *nrotc;
            *i0c1 = *i0d2 + *nrotd;
            *i0c2 = *i0c1 + nryc;
            erd__xyz_to_ry_matrix(nxyzc, *nrowc, shellc, &icore[*i0c1], &icore[*i0c2], &zcore[*z00c]);
            cdata = true;
        }
    } else {
        *nrowc = 0;
        *nrotc = 0;
        *z00c = *z00d;
        *i0c2 = *i0d2;
    }

    /* ...shell B data (being careful, using SHELLD data if SHELLC data is not present!). */
    bool bdata = false;
    if (shellb > 1) {
        if (shellb == shellc) {
            *z00b = *z00c;
            *i0b1 = *i0c1;
            *i0b2 = *i0c2;
            *nrowb = *nrowc;
            *nrotb = *nrotc;
        } else if (shellb == shelld) {
            *z00b = *z00d;
            *i0b1 = *i0d1;
            *i0b2 = *i0d2;
            *nrowb = *nrowd;
            *nrotb = *nrotd;
        } else {
            *nrowb = ((shellb / 2 + 1) * (shellb / 2 + 2)) / 2;
            *nrotb = *nrowb * nryb;
            if (cdata) {
                *z00b = *z00c + *nrotc;
                *i0b1 = *i0c2 + *nrotc;
                *i0b2 = *i0b1 + nryb;
            } else {
                *z00b = *z00d + *nrotd;
                *i0b1 = *i0d2 + *nrotd;
                *i0b2 = *i0b1 + nryb;
            }
            erd__xyz_to_ry_matrix(nxyzb, *nrowb, shellb, &icore[*i0b1], &icore[*i0b2], &zcore[*z00b]);
            bdata = true;
        }
    } else {
        *nrowb = 0;
        *nrotb = 0;
        *z00b = *z00c;
        *i0b2 = *i0c2;
    }

    /* ...shell A data (being careful, using SHELLC data if SHELLB data is not present or using SHELLD data if also SHELLC data is not present!). */
    if (shella > 1) {
        if (shella == shellb) {
            *z00a = *z00b;
            *i0a1 = *i0b1;
            *i0a2 = *i0b2;
            *nrowa = *nrowb;
            *nrota = *nrotb;
        } else if (shella == shellc) {
            *z00a = *z00c;
            *i0a1 = *i0c1;
            *i0a2 = *i0c2;
            *nrowa = *nrowc;
            *nrota = *nrotc;
        } else if (shella == shelld) {
            *z00a = *z00d;
            *i0a1 = *i0d1;
            *i0a2 = *i0d2;
            *nrowa = *nrowd;
            *nrota = *nrotd;
        } else {
            *nrowa = ((shella / 2 + 1) * (shella / 2 + 2)) / 2;
            *nrota = *nrowa * nrya;
            if (bdata) {
                *z00a = *z00b + *nrotb;
                *i0a1 = *i0b2 + *nrotb;
                *i0a2 = *i0a1 + nrya;
            } else if (cdata) {
                *z00a = *z00c + *nrotc;
                *i0a1 = *i0c2 + *nrotc;
                *i0a2 = *i0a1 + nrya;
            } else {
                *z00a = *z00d + *nrotd;
                *i0a1 = *i0d2 + *nrotd;
                *i0a2 = *i0a1 + nrya;
            }
            erd__xyz_to_ry_matrix(nxyza, *nrowa, shella, &icore[*i0a1], &icore[*i0a2], &zcore[*z00a]);
        }
    }
}
