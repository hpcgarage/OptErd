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

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__CARTESIAN_NORMS */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation generates all partial cartesian */
/*                normalization factors. The cartesian normalization */
/*                factors are xyz-exponent dependent and are given for */
/*                a specific monomial as: */
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
/*                and to precalculate each possible partial factor */
/*                only once for the range l = 0,1,...,L, where L denotes */
/*                the maximum shell quantum number that can possibly */
/*                arise during integral evaluation. */
/*                Note, that the first two factors for l = 0 and 1 */
/*                are equal to 1 and 2 and in fact these are the only */
/*                ones that are needed for s- and p-functions. */
/*                  Input: */
/*                       L         =  maximum shell quantum number */
/*                  Output: */
/*                       NORM (I)  =  cartesian norms from I=0,1,...,L */
/* ------------------------------------------------------------------------ */
ERD_OFFLOAD void erd__cartesian_norms(uint32_t length, double norm[restrict static length+1]) {
    double odd = 1.0;
    norm[0] = 1.0;
    norm[1] = 2.0;
    for (uint32_t i = 2; i <= length; i++) {
        odd += 2.0;
        norm[i] = norm[i - 1] * 2.0 / sqrt(odd);
    }
}
