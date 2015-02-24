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

/* erd__rys_x_roots_weights.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
    on Microsoft Windows system, link with libf2c.lib;
    on Linux or Unix systems, link with .../path/to/libf2c.a -lm
    or, if you install libf2c.a in a standard place, with -lf2c -lm
    -- in that order, at the end of the command line, as in
        cc *.o -lf2c -lm
    Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

        http://www.netlib.org/f2c/libf2c.zip
*/

#include "jacobi.h"
#include <stdint.h>
#include <stddef.h>
#include <math.h>
#include <assert.h>
#include "erd.h"
#include "erdutil.h"

ERD_OFFLOAD void erd__rys_x_roots_weights(int nt, int ntgqp, int ngqp,
    int nmom, const double tval[restrict],
    const double ryszero[restrict],
    double rts[restrict], double wts[restrict])
{
    ERD_SIMD_ALIGN double a[nmom];
    ERD_SIMD_ALIGN double b[nmom-1];
    ERD_SIMD_ALIGN double mom[nmom];
    ERD_SIMD_ALIGN double dia[ngqp];
    ERD_SIMD_ALIGN double off[ngqp];
    ERD_SIMD_ALIGN double row1[nmom];
    ERD_SIMD_ALIGN double row2[nmom];

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__RYS_X_ROOTS_WEIGHTS */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : Master routine to evaluate roots and weights in the */
/*                interval [0,1] over the Rys weight function: */


/*                                       exp(-T*x) */
/*                             W   (x) = --------- */
/*                              Rys      2*sqrt(x) */


/*                using the general Gaussian Quadrature technique */
/*                consisting of the following basic steps: */

/*                   1) Calculate the auxillary polynomial moments */
/*                   2) Calculate the auxillary polynomial coefficients */
/*                   3) Set up the tridiagonal symmetric terminal */
/*                      matrix */
/*                   4) Solve the tridiagonal symmetric terminal */
/*                      matrix for roots and weights. */

/*                In this routine all T-values are treated all at once. */
/*                The value of the T-parameter dictates which type of */
/*                auxillary set of polynomials is to be used for the */
/*                modified Chebyshev algorithm. The auxillary polynomials */
/*                are often chosen (as is the case here) to be orthogonal */
/*                relative to some classical weight function. */

/*                The classical weights used to establish the auxillary */
/*                polynomials are the following (notation according to */
/*                M.Abramowitz and I.A.Stegun, Handbook of Mathematical */
/*                Functions, 1964): */


/*                i) Range of validity: 0 =< T =< 30 */

/*                   Here we use shifted Jacobi weights and polynomials: */


/*                        W      (p,q,x)  =  (1-x)^(p-q) * x^(q-1) */
/*                         Jacobi */

/*                        i-th Jacobi polynomial  =  G (p,q,x) */
/*                                                    i */

/*                   with conditions: x in interval [0,1] */
/*                                    p-q > -1 , q > 0 */
/*                                    i = 0,1,2,... */


/*                ii) Range of validity: 1 < T =< oo (infinity) */

/*                   Here we use generalized Laguerre weights and */
/*                   polynomials: */


/*                        W        (a,x)  =  exp^(-x) * x^a */
/*                         Laguerre */

/*                        i-th Laguerre polynomial  =  L (a,x) */
/*                                                      i */

/*                   with conditions: x in interval [0,inf] */
/*                                    a > -1 */
/*                                    i = 0,1,2,... */


/*                Range of validity means that for all the specified */
/*                T's within the range the resulting moment integrals */
/*                are accurate to within 1.D-16. */


/*                  Input: */

/*                    NT           =  # of T-values */
/*                    NTGQP        =  # of roots times # of T-values */
/*                    NGQP         =  # of gaussian quadrature points */
/*                                    (roots) */
/*                    NMOM         =  # of necessary moment integrals */
/*                                    to calculate the quadrature roots */
/*                    TVAL         =  the T-values */
/*                    RYSZERO      =  the zeroth Rys moments for all */
/*                                    T-values */
/*                    A,B          =  will contain the recurrence */
/*                                    coefficients for the auxillary */
/*                                    polynomials */
/*                    MOM          =  will contain the normed auxillary */
/*                                    polynomial modified moments */
/*                    DIA,OFF      =  will contain the diagonal and */
/*                                    offdiagonal elements of the */
/*                                    tridiagonal symmetric terminal */
/*                                    matrix */
/*                    ROW1,ROW2    =  will be used to evaluate the */
/*                                    tridiagonal elements of the */
/*                                    symmetric terminal matrix in an */
/*                                    efficient way using Sack and */
/*                                    Donovan's method */

/*                  Output: */

/*                    RTS          =  the roots array */
/*                    WTS          =  the weights array */


/*  AUTHOR      : Norbert Flocke */
/* ------------------------------------------------------------------------ */


/*             ...include files and declare variables. */




/* ------------------------------------------------------------------------ */


/*             ...main loop over all T values. Check which T case */
/*                applies. */


/* ------------------------------------------------------------------------ */
/*  INCLUDE FILE: ERD__JACOBI */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  DESCRIPTION : Include file containing: */


/*                a) The Jacobi moment recurence coefficients R2 and */
/*                   SINV to evaluate the Jacobi moments via: */

/*                    MOM (i-1) = SINV (i) * (MOM (i+1) - R(i) * MOM (i)) */

/*                   with: */

/*                    R (i)    = (2i+1) / 2T  +  R2 (i) */
/*                    R2 (i)   = (2i+1)/(4i-1)(4i+3) */
/*                    SINV (i) = (4i-3)(4i+1)(4i-1)**2/2i(2i+1)(2i-1)**2 */

/*                   All R2 and SINV values are tabulated from i=1,100 */
/*                   with all significant figures up to quadruple */
/*                   precision. */


/*                b) The Jacobi proportionality coefficients CSMALL for */
/*                   evaluating Jacobi moments corresponding to very */
/*                   small T values. */


/*                c) Recurrence coefficients AJAC and BJAC for the */
/*                   shifted Jacobi polynomials: */


/*                        AJAC (i+1) = 4i*(2i+1)-1 / (4i+3)(4i-1) */

/*                        BJAC (i+1) = (4i**2)((4i**2)-4i+1) / */
/*                                     (4i-3)(4i+1)((4i-1)**2) */

/*                   in the range i=0,99 for A and i=1,99 for B. */


/*  AUTHOR      : Norbert Flocke */
/* ------------------------------------------------------------------------ */


/*             ...declare and set the data. */


    /* Function Body */
    int nrts = 0;
    for (int n = 0; n < nt; n += 1) {
        const double t = tval[n];
        const double momzero = ryszero[n];
        if (t <= 15.0) {
            /* 
             * ...The Jacobi section. Check first, if the number of moments wanted is within limits for the Jacobi case.
             *    If ok, we proceed by calculating two steps at the same time: 
             *
             *    1) normed shifted Jacobi modified moments defined as integrals on [0,1] over the Rys weight function: 
             *
             *                           exp(-Tx) 
             *                W   (x) = --------- 
             *                 Rys      2*sqrt(x) 
             *
             *       and shifted Jacobi polynomials: 
             *
             *                    G (0.5,0.5,x) 
             *                     i 
             *
             *       A 3-term recursion relation can be given for these moments as follows: 
             *
             *           MOM (i+1) = R * MOM (i)  +  S * MOM (i-1) 
             *
             *       where 
             *
             *           R  =  (2i+1)/2T  +  (2i+1)/(4i-1)(4i+3) 
             *           S  =  2i(2i+1)(2i-1)**2 / (4i-3)(4i+1)(4i-1)**2 
             *
             *       However, in order to evaluate the moments to the required accuracy,
             *       a downward recursion must be applied starting with some tiny seed values for MOM (i+1) and MOM (i): 
             *
             *           MOM (i-1) = 1/S * ( MOM (i+1) - R * MOM (i)) 
             *
             *       The sequence of moments thus obtained is not normalized. Hence the downward recursion is performed down to
             *       the zeroth moment ZMOM and the relevant moments are normalized by subsequent division with ZMOM. 
             *
             *       The above downward recursion technique is extremely unstable if T is very small,
             *       in which case even the small seed values are not small enough to prevent overflow of the lowest not normalized moments. 
             *       In such a case it is better to build the normed sequence starting with the first moment after the following considerations.
             *       An expansion of MOM (i) in terms of powers of T shows that the leading term is in the i-th power of T, i.e. 
             *
             *                        inf           k 
             *           MOM (i)  =   sum  c (i,k) T 
             *                       k = i 
             *
             *       with leading coefficient c (i,i) getting smaller with increasing i. Hence for very small T it is sufficient to set: 
             *
             *                                 i 
             *           MOM (i)  =   c (i,i) T     (for T very small) 
             *
             *
             *       an equation which is valid in terms of computer accuray if T is less or equal
             *       to the minimum possible nonzero number within the precision of the mantissa allowed
             *       (here in the present code this is double precision, hence T less or equal 1.0e-16).
             *       Rather than evaluating each leading coefficient c (i,i) individualy
             *       (they are given by very! complicated expressions in terms of sums of double factorials),
             *       the procedure adopted here was to predetermine each c (i,i) numerically by running the routine
             *       with T = 1.D-16 and presetting the resulting coefficients in a data array CSMALL
             *       to be used whenever T =< 1.D-16. 

             *       If T is in the range 1.D-16 to TMAX = 30.D0, then the routine was calibrated such that the mantissa 
             *       of the resulting moments are accurate to the 16th decimal place. Calibration in this context means the 
             *       predetermination of the maximum moment (controled in the code by its index IMAX) that has to be generated
             *       to perform the downward recurrence relation starting with the tiny seeds.
             *       Hence the calibration depends on four things: 
             *
             *            i) The total maximum number MOMMAX of normed moments wanted. 
             *
             *           ii) The mantissa accuray wanted for the normed moments. 
             *
             *          iii) The tiny values of the two seeds (the nonzero seed should correpond to the smallest possible nonzero number
             *               representable on the present computer). 
             *
             *           iv) T range. Each T range requires a different maximum moment (i.e. a different IMAX value) to start with. 
             *
             *       In order to perform a calibration for a different machine, one has to write a separate small program using this
             *       routine with different seeds, T- and IMAX values. Once the seeds have been set, this little program should recalculate
             *       the MOMMAX normed moments for increasing IMAX values starting with IMAX= MOMMAX. The IMAX value finally taken for a 
             *       particular T-value should then obey the following inequality: 
             *
             *
             *        | MOM (i,IMAX) - MOM (i,IMAX-1)| 
             *        | -----------------------------| < mantissa accuracy 
             *        |         MOM (i,IMAX-1)       | 
             *
             *
             *            COMMENTS FOR FAST EVALUATION OF THE MOMENTS 
             *           --------------------------------------------- 
             *       As can be seen from the downward recursion formula for the moments, all we need are the values of the R and 1/S recursion parameters,
             *       which are conveniently precomputed (with R being split as R = R1/2T + R2) and supplied via an include table 'erd__jacobi_table.c'. 
             *       Hence 'erd__jacobi_table.c' will have the values of the complicated R2 and SINV = 1/S expressions in the range from i = 1 to 100.
             *       The R1 values are simply calculated by using the appropriate decrement value of 2. The include file 'jacobi.h' will also contain
             *       the CSMALL array for the moments of very small T's. 
             *
             *
             *    2) recurrence coefficients for the shifted Jacobi polynomials G (0.5,0.5,x), denoted simply by G (x): 
             *
             *
             *             G (x)  =  1 
             *              0 
             *
             *             G (x)  =  (x - A ) 
             *              1              1 
             *
             *             G   (x)  =  (x - A   ) * G (x) - B   G   (x) 
             *              i+1              i+1     i       i+1 i-1 
             *
             *
             *       The result consists of the recurrence coefficients: 
             *
             *                   A (I) , I = 1,NMOM 
             *                   B (I) , I = 2,NMOM 
             *
             *       whose values are given by the following expressions: 
             *
             *
             *            A (i+1) = 4i*(2i+1)-1 / (4i+3)(4i-1) 
             *
             *            B (i+1) = (4i**2)((4i**2)-4i+1) / 
             *                      (4i-3)(4i+1)((4i-1)**2) 
             *
             *       Since these are complicated and time consuming expressions,
             *       they are precalculated and included into the include file 'jacobi.h'.
             */

            assert(nmom <= 30);


            /* ...the very small T case. */


            if (t <= 1.0e-16) {
                const int imax = (nmom < 16) ? nmom : 16;
                a[0] = ajac[0];
                mom[0] = csmall[0] * t;
                double tpower = t;
                for (int i = 2; i <= imax; ++i) {
                    tpower *= t;
                    a[i-1] = ajac[i - 1];
                    b[i-2] = bjac[i - 2];
                    mom[i-1] = csmall[i - 1] * tpower;
                }
                for (int i = imax + 1; i <= nmom; ++i) {
                    a[i-1] = ajac[i - 1];
                    b[i-2] = bjac[i - 2];
                    mom[i-1] = 0.;
                }
            } else {
                /*
                 * ...the general Jacobi case.
                 *    Set maximum number of moments necessary to get required moments to an accuracy of at least 1.D-16.
                 *    See above in the routine description for details and calibration of the setting.
                 */

                int imax;
                if (nmom <= 5) {
                    if (t < 1.0e-6) {
                        imax = nmom + 1;
                    } else if (t < .1) {
                        imax = nmom + 3;
                    } else if (t < 2.) {
                        imax = nmom + 7;
                    } else if (t < 10.) {
                        imax = nmom + 13;
                    } else {
                        imax = nmom + 22;
                    }
                } else {
                    if (t < 1.0e-6) {
                        imax = nmom;
                    } else if (t < .1) {
                        imax = nmom + 2;
                    } else if (t < 2.) {
                        imax = nmom + 4;
                    } else if (t < 10.) {
                        imax = nmom + 8;
                    } else {
                        imax = nmom + 16;
                    }
                }


                /*
                 * ...proceed by setting seed values for downward recursion to obtain minimal solution
                 *    and start recursive evaluation for all Jacobi moments down to first moment
                 *    (two loops are necessary here due to calculation of higher moments than actually needed later on).
                 */


                double momi = 1.0e-300;
                double momip1 = 0.0;
                const double tinvhf = .5 / t;
                double r1 = (double) ((imax << 1) + 5);
                for (int i = imax + 1; i >= nmom + 2; --i) {
                    r1 -= 2.0;
                    const double r = r1 * tinvhf + r2[i - 1];
                    const double momim1 = sinv[i - 1] * (momip1 - r * momi);
                    momip1 = momi;
                    momi = momim1;
                }
                for (int i = nmom + 1; i >= 2; --i) {
                    r1 -= 2.0;
                    const double r = r1 * tinvhf + r2[i - 1];
                    const double momim1 = sinv[i - 1] * (momip1 - r * momi);
                    mom[i - 2] = momim1;
                    momip1 = momi;
                    momi = momim1;
                }

                /*
                 * ...evaluate zeroth moment and normalize sequence.
                 *    If the absolute zeroth moment is less than the approximate absolute nonzero minimum
                 *    (set here equal to 1.D-300), the normalization looses its meaning and the calculation must be stopped.
                 *    Set also here the recurrence relation coefficients A and B for the shifted Jacobi polynomials.
                 */

                const double r = tinvhf * 3.0 + r2[0];
                const double zmom = sinv[0] * (momip1 - r * momi);
                assert(fabs(zmom) >= 1.0e-300);
                a[0] = ajac[0];
                const double zinv = 1. / zmom;
                mom[0] *= zinv;
                for (int i = 2; i <= nmom; ++i) {
                    a[i-1] = ajac[i - 1];
                    b[i-2] = bjac[i - 2];
                    mom[i-1] *= zinv;
                }
            }
        } else {

            /*
             * ...The Laguerre section. As for the Jacobi case,
             *    we calculate two steps at the same time:
             *
             *    1) normed Laguerre modified moments defined as
             *       integrals on [0,1] over the Rys weight function:
             *
             *                            exp(-Tx)
             *                 W   (x) = ---------
             *                  Rys      2*sqrt(x)
             *
             *       and monic generalized 'scaled' Laguerre polynomials:
             *
             *                    L (-0.5,Tx)
             *                     i
             *
             *       A closed formula can be given for these moments in
             *       terms of monic generalized 'scaled' Laguerre
             *       polynomials with generalization +0.5:
             *
             *           MOM (i) = SCALE * L   (+0.5,T)
             *                              i-1
             *       where
             *
             *                   - exp(-T)
             *           SCALE = ---------   ;  F0(T) = Rys zero moment
             *                    2T*F0(T)
             *
             *       The recursion relation for the +0.5 polynomials is:
             *
             *
             *         L   (+0.5,T) = R * L   (+0.5,T) - S * L   (+0.5,T)
             *          i-1                i-2                i-3
             *
             *       where
             *
             *              R  =  (T - 2i + 5/2) / T
             *              S  =  (i - 2)(i - 3/2) / T*T
             *
             *
             *       All moments MOM (i);i=1,2,...NMOM using the above
             *       outlined algorithm are evaluated.
             *
             *
             *    2) recurrence coefficients for T-scaled generalized
             *       monic Laguerre polynomials L (-0.5,Tx), denoted
             *       simply by L (Tx):
             *
             *
             *           L (Tx)    =   1
             *            0
             *
             *           L (Tx)    =   (x - A )
             *            1                  1
             *
             *           L   (Tx)  =   (x - A   ) * L (Tx) - B   L   (Tx)
             *            i+1                i+1     i        i+1 i-1
             *
             *
             *       The result consists of the recurrence coefficients:
             *
             *                 A (I) , I = 1,NMOM
             *                 B (I) , I = 2,NMOM
             *
             *       whose values are given by the following expressions:
             *
             *               A (1) = 1 / 2T
             *
             *             A (i+1) = (2i+1/2) / T
             *
             *             B (i+1) = i*(i-1/2) / (T*T)
             */

            const double texp = exp(-t);
            const double tinv = 1.0 / t;
            const double tinv2 = tinv * 2.;
            const double tinvhf = tinv * .5;
            const double tinvsq = tinv * tinv;
            const double scale = -tinvhf * texp / momzero;
            if (nmom == 1) {
                a[0] = tinvhf;
                mom[0] = scale;
            } else {
                a[0] = tinvhf;
                a[1] = tinvhf + tinv2;
                b[0] = tinvsq * .5;
                mom[0] = scale;
                double r = 1. - tinv * 1.5;
                mom[1] = scale * r;
                double s = 0.0;
                double binc = 0.5;
                double sinc = -0.5;
                double lim2 = r;
                double lim3 = 1.0;
                for (int i = 3; i <= nmom; ++i) {
                    binc += 2.;
                    a[i-1] = a[i-2] + tinv2;
                    b[i-2] = b[i - 3] + binc * tinvsq;
                    sinc += 2.;
                    r -= tinv2;
                    s += sinc * tinvsq;
                    const double lim1 = r * lim2 - s * lim3;
                    mom[i-1] = scale * lim1;
                    lim3 = lim2;
                    lim2 = lim1;
                }
            }
        }

        /*
         * ...This section calculates a symmetric terminal matrix
         *    using the normed modified moments MOM and related
         *        recurrence coefficients A and B of monic polynomials
         *        established in the previous section. The algorithm is
         *        a transcription of the LQMD algorithm described by Sack
         *        and Donovan in Numer. Mathematik 18, 465-478 (1972).
         *
         *        The needed data is as follows (NMOM = 2*NGQP-1):
         *
         *           1) Normed modified moments: MOM (I),I=1,NMOM
         *
         *           2) Recurrence coefficients of the corresponding
         *              monic polynomials: A (I),I=1,NMOM and B (I),
         *              I=2,NMOM. The recurrence relation for the monic
         *              polynomials is defined as:
         *
         *                 P (x)  =  1
         *                  0
         *
         *                 P (x)  =  (x - A )
         *                  1              1
         *
         *                 P   (x)  =  (x - A   ) * P (x) - B   P   (x)
         *                  i+1              i+1     i       i+1 i-1
         *
         *        The result will consist in the diagonal and offdiagonal
         *        parts of the symmetric terminal matrix:
         *
         *                    DIA (I) , I = 1,NGQP
         *                    OFF (I) , I = 1,NGQP-1
         *
         *        Proceed now with Sack and Donovan's algorithm.
         *        Handle low number (1 or 2) of quadrature points first.
         */

        if (ngqp == 1) {
            dia[0] = mom[0] + a[0];
        } else if (ngqp == 2) {
            const double sigma = mom[0] + a[0];
            dia[0] = sigma;
            const double theta = (a[1] - sigma) * mom[0] + mom[1] + b[0];
            off[0] = sqrt(theta);
            dia[1] = ((a[2] - sigma) * mom[1] + mom[2] + b[1] * mom[0]) / theta - mom[0] + a[1];
        } else {
            /*
             * ...Handle case for number of quadrature points > 2.
             *    Set maximum values for I and J and evaluate first diagonal element.
             */

            const int imax = ngqp - 1;
            /*
             * IF WE REMOVE STATIC, THE PROGRAM CRASHES IN RUN-TIME.
             * THIS LIKELY INDICATES A BUG SOMEWHERE
             */
            static int jmax = 0;
            jmax = ngqp + imax;
            for (int j = 1; j <= jmax; ++j) {
                row1[j-1] = mom[j-1];
            }
            double sigma = row1[0] + a[0];
            dia[0] = sigma;


            /* ...evaluate 2nd row of terminal matrix. */


            row2[0] = (a[1] - sigma) * row1[0] + row1[1] + b[0];
            double theta = row2[0];
            off[0] = sqrt(theta);
            --jmax;
            for (int j = 2; j <= jmax; ++j) {
                row2[j-1] = (a[j] - sigma) * row1[j-1] + row1[j] + b[j-1] * row1[j - 2];
            }
            sigma = row2[1] / theta - row1[0] + a[1];
            dia[1] = sigma;


            /* ...proceed with higher rows. */


            for (int i = 2; i <= imax; ++i) {
                --jmax;
                if (i % 2 == 0) {
                    for (int j = i; j <= jmax; ++j) {
                        row1[j-1] = (a[j] - sigma) * row2[j-1] + row2[j] + b[j-1] * row2[j - 2] - theta * row1[j-1];
                    }
                    sigma = a[i] - row2[i-1] / row2[i - 2] + row1[i] / row1[i-1];
                    theta = row1[i-1] / row2[i - 2];
                } else {
                    for (int j = i; j <= jmax; ++j) {
                        row2[j-1] = (a[j] - sigma) * row1[j-1] + row1[j] + b[j-1] * row1[j - 2] - theta * row2[j-1];
                    }
                    sigma = a[i] - row1[i-1] / row1[i - 2] + row2[i] / row2[i-1];
                    theta = row2[i-1] / row1[i - 2];
                }
                dia[i] = sigma;
                off[i-1] = sqrt(theta);
            }
        }

        /*
         *     ...The last section computes the gaussian quadrature roots
         *        and weights from a previously established symmetric
         *        terminal matrix by the Golub-Welsch algorithm (see
         *        G.H. Golub and J.H. Welsch, Math. of Computation 23,
         *        p. 221-230 and A1-A10, 1969), which is based on
         *        a result of Wilf (see H.S. Wilf, Mathematics for the
         *        Physical Sciences, New York: Wiley, Problem 9, p. 80).
         *
         *        Wilf has shown that if Z (K,I) is the K-th element of
         *        the I-th normalized eigenvector of the terminal matrix
         *        corresponding to the I-th eigenvalue D (I), then the
         *        roots RTS (i.e. the zeros of the N-th orthogonal
         *        monic polynomial) and weights WTS to be used for the
         *        gaussian quadrature are given by:
         *
         *                   RTS (I) = D (I)
         *                   TSi (I) = MOMZERO * (Z (1,I)**2)
         *
         *        where MOMZERO is the value of the definite integral
         *        over the weight function W (x) alone. In our case
         *        it is equal to the value of the zeroth Rys moment.
         *
         *        The present section performs hence a diagonalization
         *        of the tridiagonal symmetric terminal matrix keeping
         *        only the first component of the eigenvectors and sets
         *        the roots and weights equal to the above relations.
         *        The diagonalization code was derived from the routine
         *        IMTQL2 in the EISPACK collection and uses the implicit
         *        QL method.
         *
         *        Note, that the original diagonals DIA and offdiagonals
         *        OFF of the terminal matrix are destroyed during the
         *        diagonalization process !!!
         *
         *        Handle special case, if order of terminal matrix is 1.
         */


        if (ngqp == 1) {
            ++nrts;
            rts[nrts-1] = dia[0];
            wts[nrts-1] *= momzero;
        } else {
            /*
             * ...initialize vector for collecting first component of eigenvectors.
             *    To save space, array A is used, which can be done safely,
             *    since its dimension NMOM (# of moments) is always >= # of roots NGQP.
             */


            a[0] = 1.0;
            for (int j = 2; j <= ngqp; ++j) {
                a[j-1] = 0.0;
            }


            /* ...QL iterations. */


            off[ngqp-1] = 0.0;
            int m, iter = 0;
            for (int j = 1; j <= ngqp; ++j) {
next_iteration:
                for (m = j; m < ngqp; ++m) {
                    const double test1 = fabs(dia[m-1]) + fabs(dia[m]);
                    const double test2 = test1 + fabs(off[m-1]);
                    if (test2 == test1) {
                        break;
                    }
                }
                double p = dia[j-1];
                if (m != j) {
                    /* Root not converged */
                    assert(iter != 30);
                    ++iter;
                    double g = (dia[j] - p) / (off[j-1] * 2.);
                    double r = sqrt(g * g + 1.);
                    g = dia[m-1] - p + off[j-1] / (g + copysign(r, g));
                    double s = 1.0;
                    double c = 1.0;
                    p = 0.0;
                    for (int i = m - 1; i >= j; --i) {
                        double f = s * off[i-1];
                        const double d = c * off[i-1];
                        r = sqrt(f * f + g * g);
                        off[i] = r;
                        if (r == 0.0) {
                            dia[i] -= p;
                            off[m-1] = 0.0;
                            goto next_iteration;
                        }
                        s = f / r;
                        c = g / r;
                        g = dia[i] - p;
                        r = (dia[i-1] - g) * s + c * 2. * d;
                        p = s * r;
                        dia[i] = g + p;
                        g = c * r - d;
                        f = a[i];
                        a[i] = s * a[i-1] + c * f;
                        a[i-1] = c * a[i-1] - s * f;
                    }
                    dia[j-1] -= p;
                    off[j-1] = g;
                    off[m-1] = 0.0;
                    goto next_iteration;
                }
            }


            /*
             * ...calculate roots and weights.
             *    Since it is known that the roots must lay between 0 and 1,
             *    a check is made on them to see if they are actually within this range.
             */


            for (int i = 1; i <= ngqp; ++i) {
                const double root = dia[i-1];
                /* Quadrature root not in range 0-1 */
                assert((root >= 0.0) && (root <= 1.0));
                rts[nrts+i-1] = root;
                /* Computing 2nd power */
                const double ai = a[i-1];
                wts[nrts+i-1] *= momzero * (ai * ai);
            }
            nrts += ngqp;
        }

        /* ...next T-value */

    }
}
