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
#include <math.h>

#include "erd.h"
#include "erdutil.h"

/**
 * @brief Constructs the whole HRR transformation matrix
 * @detailes This operation constructs the whole HRR transformation matrix T, which will contain the information for the complete transformation in xyz-basis:
 *                             (e0| --> (ab|  ;  e = a+b ; a>=b 
 *   The matrix is constructed stepwise, starting from a unit matrix and operating on its columns with elementary HRR steps.
 * @param[in]  nrothrr            Maximum # of elements of @a t and @a row matrix expected during construction of final @a t and @a row.
 * @param[in]  ncolhrr            Maximum # of columns of @a t and @a row matrix expected during construction of final @a t and @a row.
 * @param[in]  nxyzet             Monomial dimension of the (e0| part.
 * @param[in]  nxyza              Monomial dimension for shell a.
 * @param[in]  nxyzp              Monomial dimension for shell p=a+b.
 * @param[in]  shella             Sheel type for shell A.
 * @param[in]  shellb             Sheel type for shell B.
 * @param[in]  shellp             Sheel type for shell P.
 * @param[in]  nabcoor            # of nonzero coordinate differences between nuclear centers A and B.
 * @param[in]  abx                The x-coordinate differences between nuclear centers A and B.
 * @param[in]  aby                The y-coordinate differences between nuclear centers A and B.
 * @param[in]  abz                The z-coordinate differences between nuclear centers A and B.
 * @param[out] in1                Starting index position of final @a nrow vector in the big @a nrow, @a t and @a row arrays (which are 2x the maximum size).
 * @param[out] in2                Starting index position final @a t and @a row matrices (@a in2) in the big @a nrow, @a t and @a row arrays (which are 2x the maximum size).
 * @param[out] nrowout            Maximum # of nonzero row labels of matrix @a t and @a row.
 * @param[out] nrow               Vector containing # of nonzero entries in columns of @a t and @a row matrix
 * @param[out] row                The nonzero row labels of matrix @a t.
 * @param[out] t                  The HRR transformation matrix. The xyz-basis for the a- and b-parts in columns of matrix T will be ordered such that a preceeds b.
 */
ERD_OFFLOAD void erd__hrr_matrix(uint32_t nrothrr, uint32_t ncolhrr,
    uint32_t nxyzet, uint32_t nxyza, uint32_t nxyzp,
    uint32_t shella, uint32_t shellb, uint32_t shellp,
    uint32_t nabcoor, double abx, double aby, double abz,
    uint32_t in1_ptr[restrict static 1], uint32_t in2_ptr[restrict static 1],
    uint32_t nrowout_ptr[restrict static 1], uint32_t nrow[restrict],
    uint32_t row[restrict], double t[restrict])
{
    /* ...accumulate T. */
    uint32_t in1 = 0;
    uint32_t out1 = in1 + ncolhrr;
    uint32_t in2 = 0;
    uint32_t out2 = in2 + nrothrr;

    /* ...form initial 'unit' T. */
    for (uint32_t i = 0; i < nxyzet; i++) {
        nrow[i] = 1;
        row[i] = i + 1;
        t[i] = 1.0;
    }
    /* ...build up the HRR transformation matrix + data. */
    uint32_t nxyzg = nxyzet;
    uint32_t nxyzh = 1;
    uint32_t nxyzgo = nxyzet;
    uint32_t nxyzho = 1;
    uint32_t nxyzi = nxyzp;
    uint32_t shellg = shellp;
    uint32_t nrowin = 1;
    uint32_t nrowout = 1;
    for (uint32_t shellh = 1; shellh <= shellb; shellh++) {
        nxyzgo -= nxyzi;
        nxyzho += shellh + 1;
        const uint32_t ngho = nxyzgo * nxyzho;
        switch (nabcoor) {
            case 3:
            {
                const uint32_t m = shellh / 3 + 1;
                nrowout += m * (m + ((shellh % 3) == 2));
                break;
            }
            case 2:
                nrowout += shellh / 2 + 1;
                break;
            case 1:
                nrowout += 1;
                break;            
        }
        erd__hrr_step(ngho, nrowin, nrowout,
                       nxyza, nxyzg, nxyzh, nxyzgo,
                       shella, shellg, shellh - 1,
                       abx, aby, abz,
                       &nrow[in1], &row[in2], &t[in2],
                       &nrow[out1], &row[out2], &t[out2]);
        nxyzh = nxyzho;
        if (shellh != shellb) {
            nxyzg = nxyzgo;
            nxyzi -= shellg + 1;
            nrowin = nrowout;
            shellg--;
        }
        ERD_SWAP(out1, in1);
        ERD_SWAP(out2, in2);
    }

/*             ...the resulting T matrix of dimension NROWOUT x (NXYZA */
/*                x NXYZH) has the following structure: the columns */
/*                are such that they come in NXYZA identical copies */
/*                NXYZH times. Hence we can condense the T matrix into */
/*                only NROWOUT x NXYZH distinct elements. The same */
/*                applies to the array NROW of size NXYZA x NXYZH, which */
/*                also allows for condensation into only NXYZH elements. */
    uint32_t l = 1;
    uint32_t m = 0;
    uint32_t n = 0;
    nrow = &nrow[in1];
    t = &t[in2];
    const uint32_t tleap = nrowout * nxyza;
    for (uint32_t j = 0; j < nxyzh; ++j) {
        nrow[j] = nrow[l];
        for (uint32_t i = 0; i < nrowout; ++i) {
            t[m + i] = t[n + i];
        }
        l += nxyza;
        m += nrowout;
        n += tleap;
    }

    *in1_ptr = in1;
    *in2_ptr = in2;
    *nrowout_ptr = nrowout;
}
