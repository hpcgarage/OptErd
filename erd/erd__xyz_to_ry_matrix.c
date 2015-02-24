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
/*  OPERATION   : ERD__XYZ_TO_RY_MATRIX */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation generates the transformation matrix for */
/*                performing a primitive gaussian integral batch */
/*                transformation from cartesian to spherical gaussian */
/*                type orbitals. */
/*                Input: */
/*                     NXYZ = dimension of xyz-monomial basis */
/*                            corresponding to the specified shell */
/*                            quantum number L. */
/*                            Must at least be equal to (L+1)*(L+2)/2. */
/*                      NRY = dimension of ry-spherical basis */
/*                            corresponding to the specified shell */
/*                            quantum number L. */
/*                            Must at least be equal to 2L+1. */
/*                   NROWMX = maximum number of xyz-monomials in one */
/*                            ry-solid harmonic function corresponding */
/*                            to the specified shell quantum number L. */
/*                            Must at least be equal to (L/2+1)*(L/2+2)/2. */
/*                        L = Shell quantum number. */
/*                     TEMP = temporary scratch array. */
/*                Output: */
/*                     NROW (I) = Number of xyz-monomials contributing */
/*                                to the I-th ry-component. */
/*                    ROW (K,I) = K-th xyz-monomial position number */
/*                                containing nonzero contribution to */
/*                                the I-th ry-component. */
/*                   TMAT (K,I) = K-th nonzero xyz-monomial coefficient */
/*                                to the I-th ry-component. */
/*                Each ry-component for a specific shell type corresponds */
/*                to a sherical harmonic function Y (L,M) with M ranging */
/*                from -L to +L. The TMAT matrix ry-component indices are */
/*                ordered such that +L corresponds to the 1st index and */
/*                -L to the last one. */
/*                For the ordering of the xyz-monomials used in the */
/*                expansion of the ry-components we use the usual */
/*                preference ordering p>q>r in the x-,y- and z-exponents */
/*                of x^p y^q z^r. For L=3 we would thus have the 10 */
/*                monomials arranged as follows: */
/*                                #  |  x  y  z */
/*                              ----------------- */
/*                                1  |  3  0  0 */
/*                                2  |  2  1  0 */
/*                                3  |  2  0  1 */
/*                                4  |  1  2  0 */
/*                                5  |  1  1  1 */
/*                                6  |  1  0  2 */
/*                                7  |  0  3  0 */
/*                                8  |  0  2  1 */
/*                                9  |  0  1  2 */
/*                               10  |  0  0  3 */
/*                This ordering is also assumed when evaluating batches */
/*                of cartesian gaussian integrals, hence the present */
/*                routine and the cartesian gaussian integral generation */
/*                routines go hand in hand and must be viewed together */
/*                when performing xyz-monomial basis ordering changes! */
/*                The ultimate goal of using the present routine is */
/*                to achieve achieve a transformation matrix T relating */
/*                cartesian gaussian functions to spherical ones. */
/*                At his stage it is convenient to introduce all */
/*                factors independent! of the gaussian exponents that */
/*                enter the final normalization constant for each */
/*                spherical gaussian function. Such a function is given */
/*                by the following expression: */
/*                    LM                     L    LM                2 */
/*                 GTO  (r,t,p) = N (L,a) * r  * Y  (t,p) * exp (-ar ) */
/*                    a */
/*                where t = theta and p = phi and N (L,a) denotes a */
/*                normalization factor such that */
/*                              LM       *    LM */
/*                 integral {GTO  (r,t,p)  GTO  (r,t,p) dr dt dp} = 1 */
/*                              a             a */
/*                Since the spherical harmonics are normalized, we can */
/*                integrate out over the angles t and p and arrive at */
/*                the following defining expression for N (L,a): */
/*                        2  r=oo  2L+2          2 */
/*                 N (L,a)   int  r     exp (-2ar ) dr  =  1 */
/*                           r=0 */
/*                which leads to: */
/*                               ____________________________ */
/*                              / 2^(2L+3+1/2) * a^((2L+3)/2) */
/*                 N (L,a) =   / ----------------------------- */
/*                           \/     (2L+1)!! * sqrt (pi) */
/*                The following part of N (L,a) has already been */
/*                taken care of during the cartesian integrals */
/*                generation: */
/*                               ____________________________ */
/*                              / 2^(1+1/2) * a^((2L+3)/2) */
/*                 N (L,a) =   / ----------------------------- */
/*                           \/          sqrt (pi) */
/*                hence we are left only with the following part */
/*                of N (L,a) to deal with here: */
/*                                   __________ */
/*                                  / 2^(2L+2) */
/*                                 / ----------- */
/*                               \/   (2L+1)!! */
/*                which is designated RNORM (to reflect its origin) */
/*                and is being calculated in logarithmic form to avoid */
/*                int overflow for the double factorial. The */
/*                normalization constant for the spherical harmonic */
/*                part, depending on both L and M quantum numbers, */
/*                is designated YNORM and is evaluated stepwise */
/*                from its base value at M = -L or L downwards to */
/*                -1 or 0. Note, that YNORM should also include a factor */
/*                in terms of pi, namely 1/sqrt(pi), but this has also */
/*                already been incorporated into the cartesian integrals */
/*                evaluation routine. */

/* ------------------------------------------------------------------------ */
ERD_OFFLOAD void erd__xyz_to_ry_matrix(
    uint32_t nxyz,
    uint32_t nrowmx,
    uint32_t l,
    uint32_t nrow[restrict static 2*l+1],
    uint32_t row[restrict static nrowmx*(2*l+1)],
    double tmat[restrict static nrowmx*(2*l+1)])
{
    double a, b, c, d;
    ERD_SIMD_ALIGN double temp[PAD_LEN(nxyz)];

/*             ...calculate the RNORM and the initial value of */
/*                YNORM at M = -L. Form initial combined norm RYNORM. */
    const double two_log10 = 0x1.34413509F79FFp-2; /* log10(2.0) */
    double rnorm = (double)(l * 2 + 2) * two_log10;
    for (uint32_t i = 0; i < l; i++) {
        rnorm -= log10((double)(int32_t)(2*i + 3));
    }
    rnorm = pow(10, rnorm * 0.5);
    double ynorm = 0.5;
    for (uint32_t i = 1; i <= l; i++) {
        ynorm = ynorm * (double)(int32_t)(2*i + 1) / (double)(int32_t)(2*i);
    }
    ynorm = sqrt(ynorm);
    double rynorm = rnorm * ynorm;

/*             ...evaluate monomial expansion for spherical harmonics */
/*                Y(L,M) for the cases M < 0. Note that the loop below */
/*                uses M, which goes from L down to 1. This M represents */
/*                in fact the absolute value of M, i.e. |M|. */
/*                There are two ways of coding this: */
/*                    1) symbolically, in which the nominators and */
/*                       denominators of all fractions are handled */
/*                       separately, performing common divisor cuts */
/*                       at every single step. This procedure is */
/*                       slower than the direct evaluation but more */
/*                       accurate. Also the symbollic procedure suffers */
/*                       from int overflow at about L = 16, so */
/*                       there is a limit on the shell size one can */
/*                       safely handle. */
/*                    2) the direct mode, in which all fractions are */
/*                       multiplied together as flp's. */
/*                This routine uses the direct mode. */
    const double xl = (double)l;
    const double xll = (double)(2 * l);
    for (uint32_t m = l; m >= 1; m--) {
        for (uint32_t xyz = 0; xyz < nxyz; xyz++) {
            temp[xyz] = 0.0;
        }
        const uint32_t lmm = l - m;
        const uint32_t lpm = l + m;
        const double xm = (double)(int32_t)m;
        const double xlpm = (double)(int32_t)lpm;
        const double xlmm = (double)(int32_t)lmm;
        const uint32_t kend = lmm / 2;
        const uint32_t pend = (m - 1) / 2;
        if (m == l) {
            a = xl;
        } else {
            if (lmm % 2 == 0) {
                const double x = xlmm / (xlpm + 1.0);
                rynorm *= sqrt(x);
                a = a * (xm / (xm + 1.0)) / x;
            } else {
                rynorm *= sqrt(xlmm * xlpm + xlmm);
                a = a * (xm / (xm + 1.0)) / xlmm;
            }
        }
        b = a;
        for (uint32_t k = 0; k <= kend; k++) {
            const double xk = (double)(int32_t) k;
            if (k != 0) {
                const double xkk = (double)(int32_t) (2 * k);
                const double x = xlmm - xkk + 1.0;
                b = -b * (x / (xll - xkk + 1.0)) * ((x + 1.0) / 2.0);
            }
            c = b;
            for (uint32_t p = 0; p <= pend; p++) {
                const uint32_t pp = (double)(int32_t) 2 * p;
                if (p != 0) {
                    const double xpp = (double)(int32_t) pp;
                    const double x = xm - xpp;
                    c = -c * (x / xpp) * ((x + 1.0) / (xpp + 1.0));
                }
                d = c;
                for (int32_t t = k; t >= 0; --t) {
                    const int tt = 2 * t;
                    const double xt = (double) t;
                    if (t != k) {
                        d /= xk - xt;
                    }
                    double xcol = d;
                    for (uint32_t v = 0; v <= t; ++v) {
                        const uint32_t vv = 2 * v;
                        if (v == 0) {
                            for (uint32_t i = 2; i <= t; ++i) {
                                xcol /= (double)(int32_t) i;
                            }
                        } else {
                            const double xv = (double)(int32_t) v;
                            xcol *= (xt - xv + 1.0) / xv;
                        }
                        const uint32_t xyz = (lmm - vv + pp + 2u) * (lmm - vv + pp + 1u) / 2u + lmm - tt;
                        temp[xyz] += xcol;
                    }
                }
            }
        }
        uint32_t i = 0;
        for (uint32_t xyz = 0; xyz < nxyz; xyz++) {
            if (temp[xyz] != 0.0) {
                row[i + lpm * nrowmx] = xyz+1;
                tmat[i + lpm * nrowmx] = rynorm * temp[xyz];
                i++;
            }
        }
        nrow[lpm] = i;
    }


/*             ...evaluate monomial expansion for spherical harmonics */
/*                Y(L,M) for the cases M >= 0. */
    rynorm = rnorm * ynorm;
    for (int32_t m = l; m >= 0; --m) {
        for (uint32_t xyz = 0; xyz < nxyz; xyz++) {
            temp[xyz] = 0.;
        }
        const uint32_t lmm = l - m;
        const uint32_t lpm = l + m;
        const double xm = (double)(int32_t)m;
        const double xlpm = (double)(int32_t)lpm;
        const double xlmm = (double)(int32_t)lmm;
        const uint32_t kend = lmm / 2;
        const uint32_t pend = m / 2;
        if (m == 0) {
            const double sqrt_half = 0x1.6A09E667F3BCDp-1;
            rynorm *= sqrt_half;
        }
        if (m == l) {
            a = 1.0;
        } else {
            if (lmm % 2 == 0) {
                const double x = xlmm / (xlpm + 1.0);
                rynorm *= sqrt(x);
                a /= x;
            } else {
                rynorm *= sqrt(xlmm * xlpm + xlmm);
                a /= xlmm;
            }
        }
        b = a;
        for (uint32_t k = 0; k <= kend; ++k) {
            const double xk = (double)(int32_t) k;
            if (k != 0) {
                const double xkk = (double)(int32_t) (2 * k);
                const double x = xlmm - xkk + 1.0;
                b = -b * (x / (xll - xkk + 1.0)) * ((x + 1.0) / 2.0);
            }
            c = b;
            for (uint32_t p = 0; p <= pend; ++p) {
                const uint32_t pp = 2 * p;
                if (p != 0) {
                    const double xpp = (double)(int32_t) pp;
                    const double x = xm - xpp + 1.0;
                    c = -c * (x / xpp) * ((x + 1.0) / (xpp - 1.0));
                }
                d = c;
                for (int32_t t = k; t >= 0; --t) {
                    const int32_t tt = 2 * t;
                    const double xt = (double) t;
                    if (t != k) {
                        d /= xk - xt;
                    }
                    double xcol = d;
                    for (uint32_t v = 0; v <= t; ++v) {
                        const uint32_t vv = 2 * v;
                        if (v == 0) {
                            for (uint32_t i = 2; i <= t; ++i) {
                                xcol /= (double)(int32_t) i;
                            }
                        } else {
                            const double xv = (double)(int32_t) v;
                            xcol *= (xt - xv + 1.0) / xv;
                        }
                        const uint32_t xyz = (lmm - vv + pp + 1) * (lmm - vv + pp) / 2 + lmm - tt;
                        temp[xyz] += xcol;
                    }
                }
            }
        }
        uint32_t i = 0;
        for (uint32_t xyz = 0; xyz < nxyz; xyz++) {
            if (temp[xyz] != 0.0) {
                row[i + lmm * nrowmx] = xyz+1;
                tmat[i + lmm * nrowmx] = rynorm * temp[xyz];
                i++;
            }
        }
        nrow[lmm] = i;
    }
}
