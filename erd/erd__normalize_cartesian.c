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
/*  OPERATION   : ERD__NORMALIZE_CARTESIAN */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation normalizes a cartesian monomial part */
/*                of a batch of contracted cartesian gaussian integrals. */
/*                The normalization factors are xyz-exponent dependent */
/*                and are given for a specific monomial as: */
/*                                    ______________________ */
/*                     l m n         /       2^(2L) */
/*                    x y z -->     / ----------------------- */
/*                                \/ (2l-1)!!(2m-1)!!(2n-1)!! */
/*                where L = l+m+n. The best way to deal with these */
/*                factors for a complete set of monomials for fixed L */
/*                is to split up each factor into its l-,m- and n- */
/*                component: */
/*                       _______        _______        _______ */
/*                      / 2^(2l)       / 2^(2m)       / 2^(2n) */
/*                     / -------  *   / -------  *   / ------- */
/*                   \/ (2l-1)!!    \/ (2m-1)!!    \/ (2n-1)!! */
/*                These factors are passed in argument as NORM. */
/*                  Input: */
/*                    M           =  # of elements not involved in the */
/*                                   normalization (invariant indices) */
/*                    NXYZ        =  # of monomials for shell L */
/*                    L           =  the shell type for which the */
/*                                   normalization will be done */
/*                    NORM        =  individual normalization factors */
/*                                   for each monomial exponent */
/*                    BATCH       =  batch of unnormalized integrals */
/*                  Output: */
/*                    BATCH       =  batch of normalized integrals */
/* ------------------------------------------------------------------------ */
ERD_OFFLOAD void erd__normalize_cartesian(uint32_t m, uint32_t l, const double norm[restrict static l+1], double batch[restrict]) {
    for (ssize_t x = l; x >= 0; x--) {
        const double xnorm = norm[x];
        const ssize_t ybeg = l - x;
        for (ssize_t y = ybeg; y >= 0; y--) {
            const double scalar = xnorm * norm[y] * norm[ybeg - y];
            for (uint32_t i = 0; i < m; i++) {
                *batch++ *= scalar;
            }
        }
    }
}
