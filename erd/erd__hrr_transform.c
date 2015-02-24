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
/*  OPERATION   : ERD__HRR_TRANSFORM */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation performs a HRR transformation on a */
/*                batch of contracted cartesian gaussian integrals: */
/*                                      nxyzet */
/*                     y (m,nxyzab)  =   sum   x (m,i) * rot (i,nxyzab) */
/*                                       i=1 */
/*                where rot is the HRR transformation matrix. Due to */
/*                the very sparse nature of this matrix, only those */
/*                i indices in the summation are addressed which */
/*                correspond to nonzero HRR transformation matrix */
/*                elements. */
/*                  Input: */
/*                     M          =  # of elements not involved in the */
/*                                   transformation (invariant indices) */
/*                     NROW       =  maximum # of nonzero row elements */
/*                                   per column in transformation matrix */
/*                     NXYZET     =  dimension of the cartesian e0-part */
/*                     NXYZAB     =  dimension of the resulting cartesian */
/*                                   ab-part */
/*                     NXYZA      =  dimension of the cartesian part due */
/*                                   to the a-shell */
/*                     NXYZB      =  dimension of the cartesian part due */
/*                                   to the b-shell */
/*                     LROW (N)   =  # of nonzero entries in column N of */
/*                                   ROT and ROW matrix. N ranges from */
/*                                   1 to NXYZB */
/*                     ROW (I,N)  =  I-th nonzero row label of column N */
/*                                   in ROT matrix. N ranges from 1 to */
/*                                   NXYZAB */
/*                     ROT (I,N)  =  I-th nonzero HRR transformation */
/*                                   matrix element of column N. N ranges */
/*                                   from 1 to NXYZB */
/*                     X          =  batch of untransformed integrals */
/*                                   (m,e0) */
/*                  Output: */
/*                     Y          =  batch of HRR transformed integrals */
/*                                   (m,ab) */
/* ------------------------------------------------------------------------ */
ERD_OFFLOAD void erd__hrr_transform(uint32_t m, uint32_t nrow,
    uint32_t nxyza, uint32_t nxyzb,
    const uint32_t lrow[restrict static nxyzb], const uint32_t row[restrict],
    const double rot[restrict], const double x[restrict], double y[restrict])
{
/*             ...perform the HRR transformation. One of the main */
/*                properties of this transformation is that the */
/*                last nonzero element of the HRR transformation */
/*                matrix is always equal to 1. Hence we can skip */
/*                the multiplication with that element. */
/*                Use basic row grouping of the transformation */
/*                to improve cache line reusing. */
    uint32_t n = 0;
    for (uint32_t b = 0; b < nxyzb; b++) {
        const uint32_t mrow = lrow[b];
        for (uint32_t a = 0; a < nxyza; a++) {
            for (uint32_t j = 0; j < m; j++) {
                y[j + n * m] = 0.0;
            }           
            for (uint32_t i = 0; i < mrow; i++) {
                const uint32_t xcol1 = row[i + n * nrow] - 1;
                const double rot1 = rot[i + b * nrow];
                for (uint32_t j = 0; j < m; j++) {
                    y[j + n * m] += rot1 * x[j + xcol1 * m];
                }
            }
            n++;
        }      
    }
}
