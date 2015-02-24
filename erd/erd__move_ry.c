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
/*  OPERATION   : ERD__MOVE_RY */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : ERD__TRANSPOSE_BATCH */
/*                ERD__MAP_IJKL_TO_IKJL */
/*  DESCRIPTION : This routine moves all ry-components located on the */
/*                far right in array X to a specific position to the left */
/*                in array Y: */
/*                     X (1,2,3,...,RY) ---> Y ( 1, 2, 3,...) */
/*                                   |          |  |  |   | */
/*                                    --->------^--^--^...^ */
/*                The part of X that is not moved (i.e. the # of */
/*                invariant indices to the left in X) has been calculated */
/*                beforehand and is transmitted through argument. */
/*                  Input: */
/*                       NINTGRL    =  total # of integrals */
/*                       NINDEX     =  total # of possible target places */
/*                       NOTMOVE    =  inactive # of indices */
/*                       MOVE       =  # of indices that will be moved */
/*                       NRY        =  # of ry-components to be moved */
/*                       INDEX      =  place 1,2,3,... to which the */
/*                                     ry-components will be placed in */
/*                                     array Y (must be within the range */
/*                                     1 =< INDEX =< NINDEX) */
/*                       TILE       =  the level 1 cache tile */
/*                       X          =  initial set of integrals */
/*                       IXOFF (I)  =  NINDEX-fold array indicating total */
/*                                     # of elements preceeding place */
/*                                     I = 1,2,3,... before the move */
/*                                     (with IXOFF(1)=1) */
/*                  Output: */
/*                       IXOFF (I)  =  updated NINDEX-fold array */
/*                                     indicating total # of elements */
/*                                     preceeding place I = 1,2,3,... */
/*                                     after the move */
/*                       Y          =  final set of integrals */
/* ------------------------------------------------------------------------ */
void erd__move_ry(uint32_t nindex, uint32_t notmove, uint32_t move, uint32_t nry,
    uint32_t index, const double x[restrict static notmove * move * nry],
    uint32_t ixoff[restrict static nindex], double y[restrict static notmove * move * nry])
{
/*             ...check, if the move is simply a transposition */
/*                or a more general ijkl -> ikjl move. */
    if (notmove == 1) {
        for (uint32_t j = 0; j < move; j++) {
            for (uint32_t i = 0; i < nry; i++) {
                y[j * nry + i] = x[i * move + j];
            }
        }
    } else {
        for (uint32_t k = 0; k < nry; k++) {
            for (uint32_t j = 0; j < move; j++) {
                for (uint32_t i = 0; i < notmove; i++) {
                    y[(j * nry + k) * notmove + i] = x[(k * move + j) * notmove + i];
                }
            }
        }
    }

    for (uint32_t i = index; i < nindex; i++) {
        ixoff[i] *= nry;
    }
}

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
