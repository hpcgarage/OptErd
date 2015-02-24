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
 * @brief Performs a single horizontal recurrence relation (HRR) step
 * @details This operation performs a single horizontal recurrence relation (HRR) step on the input transformation matrix WIN of formal
 *   dimension NXYZET x NAB, where NXYZET is the dimension of the starting shell combination (e0|.
 *   The columns of the matrix WIN are combined such that: 
 *     (a(b+1i)| = ((a+1i)b| + (Ai-Bi)(ab| ; i=x,y,z
 *   Note, that the shell symbol a stands for a range of shells, starting from shell x and ending at shell p-1 for (a(b+1i)| and (ab| and
 *   ending at shell p for ((a+1i)b|.
 *   Since the input transformation matrix WIN contains already all the HRR info from previous steps, the output transformation matrix will
 *   have the complete information for performing a HRR of the following kind:
 *                           (e0| --> (a(b+1i)|
 *                Strategy used to perform the HRR step:
 *                -------------------------------------
 *                The HRR is split into two parts: part I) deals with the (Ai-Bi)(ab| part and part II) with the ((a+1i)b| part.
 *                Part I):  In this case we have the following scheme, with a-shell range a = x to p-1:
 *                           W = WIN (xa,ya,za,xb,yb,zb) ->
 *                               + ABX*W to WOUT (xa,ya,za,xb+1,yb,zb)
 *                               + ABY*W to WOUT (xa,ya,za,xb,yb+1,zb)
 *                               + ABZ*W to WOUT (xa,ya,za,xb,yb,zb+1)
 *                Part II): Here we have the following scheme, with a-shell range a = x+1 to p:
 *                           W = WIN (xa,ya,za,xb,yb,zb) ->
 *                               + W to WOUT (xa-1,ya,za,xb+1,yb,zb)
 *                               + W to WOUT (xa,ya-1,za,xb,yb+1,zb)
 *                               + W to WOUT (xa,ya,za-1,xb,yb,zb+1)
 *                How to get from the b- to the (b+1)-shell monomials:
 *                ---------------------------------------------------
 *                To perform parts I) and II) of the HRR we need an algorithm which creates a unique set of (b+1)-shell monomials from
 *                those of the b-shell. The following strategy is adopted:
 *                     1) if x-exponent in b-shell monomial is > 0,
 *                        add +1 to x-exponent.
 *                     2) if x-exponent in b-shell monomial is = 0 and y-exponent is > 0, add +1 to x-exponent and y-exponent.
 *                     3) if x-exponent and y-exponent in b-shell monomial is = 0, add +1 to all exponents.
 *                Example for b=2 --> b+1=3 case:
 *                            200 --> 300
 *                            110 --> 210
 *                            101 --> 201
 *                            020 --> 120,030
 *                            011 --> 111,021
 *                            002 --> 102,012,003
 *
 * @param[in]  nabo         Total # of monomials of the output (a(b+1i)| part, i.e. # of columns of output transformation matrix
 * @param[in]  nrowin       Maximum # of nonzero row elements per column in input transformation matrix
 * @param[in]  mrowout      Maximum # of nonzero row elements per column in output transformation matrix
 * @param[in]  nxyzx        Total # of monomials for shell x
 * @param[in]  nxyzp        Total # of monomials for shell p
 * @param[in]  nxyza        Sum of total # of monomials for shell range a = x to p
 * @param[in]  nxyzb        Total # of monomials for shell b
 * @param[in]  nxyzao       Sum of total # of monomials for shell range a = x to p-1
 * @param[in]  shellx       Shell type for shell X
 * @param[in]  shellp       Shell type for shell P
 * @param[in]  shellb       Shell type for shell B
 * @param[in]  abx          The x-coordinate differences between centers A and B
 * @param[in]  aby          The y-coordinate differences between centers A and B
 * @param[in]  abz          The z-coordinate differences between centers A and B
 * @param[in]  nrowin(j)    # of nonzero row elements of J-th input transformation matrix column
 * @param[in]  rowin(i,j)   Nonzero row indices of J-th input transformation matrix column
 * @param[in]  win(i,j)     Nonzero elements of the input transformation matrix corresponding to nonzero row index I and column J
 * @param[out] nrowout(j)   # of nonzero row elements of J-th output transformation matrix column
 * @param[out] rowout(i,j)  Nonzero row indices of J-th output transformation matrix column
 * @param[out] wout(i,j)    nonzero elements of the output transformation matrix corresponding to nonzero row index I and column J
 */
ERD_OFFLOAD void erd__hrr_step(uint32_t nabo, uint32_t mrowin,
    uint32_t mrowout, uint32_t nxyzx,
    uint32_t nxyza, uint32_t nxyzb, uint32_t nxyzao,
    uint32_t shellx, uint32_t shellp, uint32_t shellb,
    double abx, double aby, double abz,
    const uint32_t nrowin[restrict], const uint32_t rowin[restrict], const double win[restrict],
    uint32_t nrowout[restrict static nabo], uint32_t rowout[restrict], double wout[restrict])
{
    /*
     * Hold the pair of column indices of the input transformation matrix that
     * have to be combined to form each output transformation matrix column
     */
    uint32_t cpair[2*nabo];
    const int32_t nxbgt0 = nxyzb - shellb - 1;
    const uint32_t offybo = nxyzao * (shellb + 1);
    /*
     * ...determine the column pairs to be 'added' to define the output tranformation matrix.
     * ------- Part I) section ---------
     * Outer loop over b- to (b+1)-shell monomials with xb-exponent > 0.
     * Inner loop over a-part contractions for part I) of HRR for a-shell range a = x to p-1.
     */
    uint32_t x = -nxyza;
    uint32_t xo = -nxyzao;
    for (int32_t b = 0; b < nxbgt0; b++) {
        x += nxyza;
        xo += nxyzao;
        if (abx != 0.0) {
            for (uint32_t j = 0; j < nxyzao; j++) {
                cpair[xo + j] = x + j;
            }
        }
    }
    /*
     * ...outer loop over b- to (b+1)-shell monomials with xb-exponent = 0 and yb-exponent > 0.
     * Inner loop over a-part contractions for part I) of HRR for a-shell range a = x to p-1.
     */
    uint32_t yo = xo + offybo;
    for (uint32_t b = 0; b < shellb; b++) {
        x += nxyza;
        xo += nxyzao;
        yo += nxyzao;
        if (abx != 0.0) {
            for (uint32_t j = 0; j < nxyzao; j++) {
                cpair[xo + j] = x + j;
            }
        }
        if (aby != 0.0) {
            for (uint32_t j = 0; j < nxyzao; j++) {
                cpair[yo + j] = x + j;
            }
        }
    }
    /*
     * ...last b- to (b+1)-shell monomial with xb-exponent = 0 and yb-exponent = 0.
     * Inner loop over a-part contractions for part I) of HRR for a-shell range a = x to p-1.
     */
    x += nxyza;
    xo += nxyzao;
    yo += nxyzao;
    uint32_t zo = yo + nxyzao;
    if (abx != 0.0) {
        for (uint32_t j = 0; j < nxyzao; j++) {
            cpair[xo + j] = x + j;
        }
    }
    if (aby != 0.0) {
        for (uint32_t j = 0; j < nxyzao; j++) {
            cpair[yo + j] = x + j;
        }
    }
    if (abz != 0.0) {
        for (uint32_t j = 0; j < nxyzao; j++) {
            cpair[zo + j] = x + j;
        }
    }

    /*
     * ------- Part II) section --------- 
     * ...outer loop over b- to (b+1)-shell monomials with xb-exponent > 0.
     * Inner loop over a-part contractions for part II) of HRR for a-shell range a = x+1 to p.
     */
    x = 0;
    xo = 0;
    for (int32_t b = 0; b < nxbgt0; b++) {
        x += nxyzx;
        for (uint32_t shella = shellx; shella < shellp; shella++) {
            for (int32_t xa = shella; xa >= 0; xa--) {
                for (int32_t ya = shella - xa; ya >= 0; ya--) {
                    cpair[xo + nabo] = x;
                    x++;
                    xo++;
                }
            }
            x += shella + 2;
        }
    }

    /*
     * ...outer loop over b- to (b+1)-shell monomials with xb-exponent = 0 and yb-exponent > 0.
     * Inner loop over a-part contractions for part II) of HRR for a-shell range a = x+1 to p.
     */
    yo = xo + offybo;
    for (uint32_t b = 0; b < shellb; b++) {
        x += nxyzx;
        for (uint32_t shella = shellx; shella < shellp; ++shella) {
            uint32_t offya = 0;
            uint32_t y;
            for (int32_t xa = shella; xa >= 0; xa--) {
                ++offya;
                for (uint32_t ya = offya; ya != 0; ya--) {
                    ++x;
                    cpair[xo + nabo] = x - 1;
                    cpair[yo + nabo] = x + offya - 1;
                    y = x + offya;
                    ++xo;
                    ++yo;
                }
            }
            x = y + 1;
        }
    }

    /*
     * ...last b- to (b+1)-shell monomial with xb-exponent = 0 and yb-exponent = 0.
     * Inner loop over a-part contractions for part II) of HRR for a-shell range a = x+1 to p.
     */
    x += nxyzx;
    zo = yo + nxyzao;
    for (uint32_t shella = shellx; shella < shellp; shella++) {
        uint32_t offya = 0;
        uint32_t z;
        for (int32_t xa = shella; xa >= 0; xa--) {
            ++offya;
            for (uint32_t ya = offya; ya != 0; ya--) {
                ++x;
                const uint32_t y = x + offya;
                z = y + 1;
                cpair[xo + nabo] = x - 1;
                cpair[yo + nabo] = y - 1;
                cpair[zo + nabo] = z - 1;
                ++xo;
                ++yo;
                ++zo;
            }
        }
        x = z;
    }

    /* ...the column pairs are ready. Construct the new transformation matrix. */
    for (uint32_t n = 0; n < nabo; n++) {
        double coeff;
        if (n < xo) {
            coeff = abx;
        } else if (n < yo) {
            coeff = aby;
        } else {
            coeff = abz;
        }
        if (coeff != 0.0) {
            const uint32_t c1 = cpair[n];
            const uint32_t c2 = cpair[n + nabo];
            const uint32_t nrow1 = nrowin[c1];
            const uint32_t nrow2 = nrowin[c2];
            const uint32_t m = min32u(nrow1, nrow2);
            uint32_t i1 = 0;
            uint32_t i2 = 0;
            uint32_t nout = 0;
            for (uint32_t i = 1; i < m; i++) {
                const uint32_t row1 = rowin[i1 + c1 * mrowin];
                const uint32_t row2 = rowin[i2 + c2 * mrowin];
                if (row1 == row2) {
                    wout[nout + n * mrowout] = coeff * win[i1 + c1 * mrowin] + win[i2 + c2 * mrowin];
                    rowout[nout + n * mrowout] = row1;
                    nout++;
                    i1++;
                    i2++;
                } else if (row1 < row2) {
                    wout[nout + n * mrowout] = coeff * win[i1 + c1 * mrowin];
                    rowout[nout + n * mrowout] = row1;
                    nout++;
                    i1++;
                } else {
                    wout[nout + n * mrowout] = win[i2 + c2 * mrowin];
                    rowout[nout + n * mrowout] = row2;
                    nout++;
                    i2++;
                }
            }
            for (;;) {
                const uint32_t row1 = rowin[i1 + c1 * mrowin];
                const uint32_t row2 = rowin[i2 + c2 * mrowin];
                if (row1 == row2) {
                    wout[nout + n * mrowout] = coeff * win[i1 + c1 * mrowin] + win[i2 + c2 * mrowin];
                    rowout[nout + n * mrowout] = row1;
                    nout++;
                    i1++;
                    i2++;
                } else if (row1 < row2) {
                    wout[nout + n * mrowout] = coeff * win[i1 + c1 * mrowin];
                    rowout[nout + n * mrowout] = row1;
                    nout++;
                    i1++;
                } else {
                    wout[nout + n * mrowout] = win[i2 + c2 * mrowin];
                    rowout[nout + n * mrowout] = row2;
                    nout++;
                    i2++;
                }
                if (i1 >= nrow1) {
                    for (uint32_t i = i2; i < nrow2; i++) {
                        const uint32_t row2 = rowin[i + c2 * mrowin];
                        wout[nout + n * mrowout] = win[i + c2 * mrowin];
                        rowout[nout + n * mrowout] = row2;
                        nout++;
                    }
                    break;
                } else if (i2 >= nrow2) {
                    for (uint32_t i = i1; i < nrow1; i++) {
                        const uint32_t row1 = rowin[i + c1 * mrowin];
                        wout[nout + n * mrowout] = coeff * win[i + c1 * mrowin];
                        rowout[nout + n * mrowout] = row1;
                        nout++;
                    }
                    break;
                }
            }
            nrowout[n] = nout;
        } else {
            const uint32_t c2 = cpair[n + nabo];
            const uint32_t nrow2 = nrowin[c2];
            for (uint32_t i2 = 0; i2 < nrow2; i2++) {
                wout[i2 + n * mrowout] = win[i2 + c2 * mrowin];
                rowout[i2 + n * mrowout] = rowin[i2 + c2 * mrowin];
            }
            nrowout[n] = nrow2;
        }
    }
}
