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
/*  OPERATION   : ERD__SPHERICAL_TRANSFORM */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation performs a simple cartesian to spherical */
/*                transformation step on a batch of contracted cartesian */
/*                gaussian integrals: */
/*                         y (m,r) = sum  mat (i,r) * x (m,i) */
/*                                    i */
/*                where i is an element of the cartesian basis, r an */
/*                element of the spherical basis, mat is the cartesian */
/*                -> spherical transformation matrix and m are the bra */
/*                elements not involved. Ordering of the cartesian */
/*                monomial basis is assumed to be: */
/*                         e f g     a b c */
/*                        X Y Z  >  X Y Z   , if e > a */
/*                                            for e=a if f > b */
/*                                            for e=a and f=b if g > c */
/*                  Input: */
/*                    M           =  # of elements not involved in the */
/*                                   transformation (invariant indices) */
/*                    NROW        =  maximum # of nonzero rows in the */
/*                                   transformation matrix */
/*                    NXYZ        =  # of cartesian monomials xyz for */
/*                                   the shell to be transformed */
/*                    NRY         =  # of spherical functions ry for */
/*                                   the shell to be transformed */
/*                    LROW (R)    =  # of xyz-monomials contributing */
/*                                   to the R-th ry-component */
/*                    ROW (I,R)   =  I-th xyz-monomial row index */
/*                                   containing nonzero contribution to */
/*                                   the R-th ry-component */
/*                    ROT (I,R)   =  I-th nonzero xyz-monomial to R-th */
/*                                   ry-component transformation matrix */
/*                                   element */
/*                    X           =  input batch of cartesian integrals */
/*                  Output: */
/*                    Y           =  output batch of spherical integrals */
/* ------------------------------------------------------------------------ */
ERD_OFFLOAD void erd__spherical_transform(uint32_t m, uint32_t nrow, uint32_t nry, uint32_t lrow[restrict static nry], uint32_t row[restrict], const double rot[restrict], const double x[restrict], double y[restrict]) {   
/*             ...perform the cartesian -> spherical transformation. */
/*                Use basic row grouping of the transformation */
/*                to improve cache line reusing. */
    for (uint32_t r = 0; r < nry; r++) {
        const uint32_t mrow = lrow[r];

        for (uint32_t n = 0; n < m; n++) {
            y[r * m + n] = 0.0;
        }
        
        for (uint32_t i = 0; i < mrow; i++) {
            const uint32_t xcol1 = row[r * nrow + i] - 1;
            const double rot1 = rot[r * nrow + i];
            for (uint32_t n = 0; n < m; n++) {
                y[r * m + n] += rot1 * x[xcol1 * m + n];
            }
        }
    }
}
