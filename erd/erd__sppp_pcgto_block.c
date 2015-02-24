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
#include "boys.h"
#include "erd.h"
#include "erdutil.h"

/* ------------------------------------------------------------------------ */
/*  OPERATION   : ERD__SPPP_PCGTO_BLOCK */
/*  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT */
/*  MODULE-ID   : ERD */
/*  SUBROUTINES : none */
/*  DESCRIPTION : This operation is designed to provide ultrafast block */
/*                evaluation of a batch of normalized electron repulsion */
/*                integrals between s-shell and p-shell primitive */
/*                spherical gaussian type orbitals. */
/*                A batch is defined here as containing all possible */
/*                integrals, that is its dimension is determined by */
/*                the total number of primitive functions (here = 27) */
/*                times the total number of ij and kl exponent pair */
/*                combinations. */
/*                The integrals are ordered in the batch the following */
/*                way (first index varying fastest): */
/*                    batch (nxyz1,nxyz2,nxyz3,nxyz4,kl,ij) */
/*                where ij and kl indicates alpha exponent pairs */
/*                defining the present block. */
/*                The present routine evaluates batches of the type: */
/*                            sppp , pspp , ppsp , ppps */
/*                The cartesian primitive integrals are evaluated using */
/*                the auxillary functions technique, described by */
/*                V.R. Saunders, "An Introduction to Molecular Integral */
/*                Evaluation" in "Computational Techniques in Quantum */
/*                Chemistry and Molecular Physics" edited by */
/*                GHF Diercksen, BT Sutcliff and A Veillard, D. Reidel */
/*                Publ. Comp. Dordrecht (1975), p. 347. */
/*                The cartesian primitive integrals are each evaluated */
/*                explicitely using common multipliers if possible and */
/*                avoiding array addresses. */
/*                  Input: */
/*                    NBATCH       =  size of the primitive cartesian */
/*                                    sppp/pspp/ppsp/ppps integral batch */
/*                    ATOMIC       =  indicates, if purely atomic */
/*                                    integrals will be evaluated */
/*                    ATOM12(34)   =  indicates, if centers 1 and 2 */
/*                                    (3 and 4) coincide */
/*                    MIJ(KL)      =  current # of ij (kl) primitive */
/*                                    index pairs corresponding to */
/*                                    the contracted shell pairs 1,2 */
/*                                    (3,4) */
/*                    NIJ          =  total # of ij primitive index */
/*                                    pairs for the contracted shell */
/*                                    pair 1,2 */
/*                    NIJBEG(END)  =  first(last) ij primitive index */
/*                                    defining the ij block */
/*                    NKL          =  total # of kl primitive index */
/*                                    pairs for the contracted shell */
/*                                    pair 3,4 */
/*                    NKLBEG(END)  =  first(last) kl primitive index */
/*                                    defining the kl block */
/*                    NPGTOx       =  # of primitives per contraction */
/*                                    for contraction shells x = 1,2,3,4 */
/*                    SHELLx       =  the shell type for contraction */
/*                                    shells x = 1,3,P=1+2 */
/*                    Xx,Yx,Zx     =  the x,y,z-coordinates for centers */
/*                                    x = 1,2,3,4 */
/*                    Xxx,Yxx,Zxx  =  the x,y,z-coordinate differences */
/*                                    between centers xx = 12 and 34 */
/*                    ALPHAx       =  the primitive exponents for */
/*                                    contraction shells x = 1,2,3,4 */
/*                    FTABLE       =  Fm (T) table for interpolation */
/*                                    in low T region */
/*                    MGRID        =  maximum m in Fm (T) table */
/*                    NGRID        =  # of T's for which Fm (T) table */
/*                                    was set up */
/*                    TMAX         =  maximum T in Fm (T) table */
/*                    TSTEP        =  difference between two consecutive */
/*                                    T's in Fm (T) table */
/*                    TVSTEP       =  Inverse of TSTEP */
/*                    PRIMx        =  i,j,k,l labels of primitives for */
/*                                    the respective contraction shells */
/*                                    x = 1,2,3,4 */
/*                    NORMx        =  the normalization factors due to */
/*                                    the primitive exponents for the */
/*                                    contraction shells x = 1,2,3,4 */
/*                    RHO12(34)    =  the complete set of NIJ (NKL) */
/*                                    exponential prefactors between */
/*                                    contraction shells 1 and 2 */
/*                                    (3 and 4) */
/*                    P            =  will hold current MIJ exponent */
/*                                    sums for contraction shells 1 */
/*                                    and 2 */
/*                    Px           =  will hold current MIJ coordinates */
/*                                    x=X,Y,Z for the gaussian product */
/*                                    centers P=1+2 */
/*                    SCALEP       =  will hold current MIJ values of */
/*                                    scaling factors related to point P */
/*                    Q            =  will hold current MKL exponent */
/*                                    sums for contraction shells 3 */
/*                                    and 4 */
/*                    Qx           =  will hold current MKL coordinates */
/*                                    x=X,Y,Z for the gaussian product */
/*                                    centers Q=3+4 */
/*                    SCALEQ       =  will hold current MKL values of */
/*                                    scaling factors related to point Q */
/*                  Output: */
/*                    BATCH        =  current batch of primitive */
/*                                    cartesian sppp/pspp/ppsp/ppps */
/*                                    integrals */
/* ------------------------------------------------------------------------ */
ERD_OFFLOAD void erd__sppp_pcgto_block(uint32_t nij, uint32_t nkl,
    uint32_t shell1, uint32_t shell3, uint32_t shellp,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double x3, double y3, double z3,
    double x4, double y4, double z4,
    const double *restrict alpha1, const double *restrict alpha2, const double *restrict alpha3, const double *restrict alpha4,
    const double *restrict cc1, const double *restrict cc2, const double *restrict cc3, const double *restrict cc4,
    const uint32_t *restrict prim1, const uint32_t *restrict prim2, const uint32_t *restrict prim3, const uint32_t *restrict prim4,
    const double *restrict norm1, const double *restrict norm2, const double *restrict norm3, const double *restrict norm4,
    const double *restrict rho12, const double *restrict rho34,
    double *restrict cbatch)
{
    const double x12 = x1 - x2;
    const double y12 = y1 - y2;
    const double z12 = z1 - z2;
    const double x34 = x3 - x4;
    const double y34 = y3 - y4;
    const double z34 = z3 - z4;

    ERD_SIMD_ALIGN double p[PAD_LEN(nij)], px[PAD_LEN(nij)], py[PAD_LEN(nij)], pz[PAD_LEN(nij)], scalep[PAD_LEN(nij)];
    for (uint32_t ij = 0; ij < nij; ij += 1) {
        const uint32_t i = prim1[ij];
        const uint32_t j = prim2[ij];
        const double exp1 = alpha1[i];
        const double exp2 = alpha2[j];
        double pval = exp1 + exp2;
        p[ij] = pval;
        pval = exp1 / pval;
        px[ij] = pval * x12 + x2;
        py[ij] = pval * y12 + y2;
        pz[ij] = pval * z12 + z2;
        scalep[ij] = cc1[i] * cc2[j] * norm1[i] * norm2[j] * rho12[ij];
    }

    ERD_SIMD_ALIGN double q[PAD_LEN(nkl)], qx[PAD_LEN(nkl)], qy[PAD_LEN(nkl)], qz[PAD_LEN(nkl)], scaleq[PAD_LEN(nkl)];
    for (uint32_t kl = 0; kl < nkl; kl += 1) {
        const uint32_t k = prim3[kl];
        const uint32_t l = prim4[kl];
        const double exp3 = alpha3[k];
        const double exp4 = alpha4[l];
        double qval = exp3 + exp4;
        q[kl] = qval;
        qval = exp3 / qval;
        qx[kl] = qval * x34 + x4;
        qy[kl] = qval * y34 + y4;
        qz[kl] = qval * z34 + z4;
        scaleq[kl] = cc3[k] * cc4[l] * norm3[k] * norm4[l] * rho34[kl];
    }

    if (shellp == 1) {
        // 1     5   |  (AB|CD)  4-center   sppp and pspp
        double pxsub, pysub, pzsub;
        if (shell1 == 1) {
            pxsub = x1;
            pysub = y1;
            pzsub = z1;
        } else {
            pxsub = x2;
            pysub = y2;
            pzsub = z2;
        }
        for (uint32_t ij = 0; ij < nij; ij += 1) {
            const double pval = p[ij];
            const double pxval = px[ij];
            const double pyval = py[ij];
            const double pzval = pz[ij];
            const double pscale = scalep[ij];
            const double xpss1 = pxval - pxsub;
            const double ypss1 = pyval - pysub;
            const double zpss1 = pzval - pzsub;
            for (uint32_t kl = 0; kl < nkl; kl += 1) {
                const double qval = q[kl];
                const double qinv = 1. / qval;
                const double qxval = qx[kl];
                const double qyval = qy[kl];
                const double qzval = qz[kl];
                const double pqmult = pval * qval;
                const double pqplus = pval + qval;
                const double pqpinv = 1. / pqplus;
                const double pqx = pxval - qxval;
                const double pqy = pyval - qyval;
                const double pqz = pzval - qzval;
                const double t = (pqx * pqx + pqy * pqy + pqz * pqz) * pqmult * pqpinv;
                const double scale = pscale * scaleq[kl] / (pqmult * __builtin_sqrt(pqplus));
                const struct Boys0123 boys = boys0123(t, scale);
                const double u0 = pval * pqpinv;
                const double u1 = -qval * pqpinv;
                const double u2 = pqpinv * .5;
                const double u3 = u2 + pqpinv;
                const double u4 = qinv * .5;
                const double u5 = u0 * u4;

                /* ...the X-terms. */
                const double xssp1 = qxval - x4;
                const double xsps1 = qxval - x3;
                const double xssp2 = pqx * u0;
                const double xpss2 = pqx * u1;
                double a = xsps1 + xssp1;
                double b = xpss1 * xssp2 + u2;
                const double xspp1 = xsps1 * xssp1 + u4;
                const double xspp2 = a * xssp2 - u5;
                const double xspp3 = xssp2 * xssp2;
                const double xpsp1 = xpss1 * xssp1;
                const double xpsp2 = xssp1 * xpss2 + b;
                const double xpsp3 = xssp2 * xpss2;
                const double xpps1 = xpss1 * xsps1;
                const double xpps2 = xsps1 * xpss2 + b;
                const double xppp1 = xpss1 * xspp1;
                const double xppp2 = xpss1 * xspp2 + xspp1 * xpss2 + a * u2;
                const double xppp3 = xpss1 * xspp3 + a * xpsp3 + u3 * xssp2;
                const double xppp4 = xpss2 * xspp3;

                /* ...the Y-terms. */
                const double yssp1 = qyval - y4;
                const double ysps1 = qyval - y3;
                const double yssp2 = pqy * u0;
                const double ypss2 = pqy * u1;
                a = ysps1 + yssp1;
                b = ypss1 * yssp2 + u2;
                const double yspp1 = ysps1 * yssp1 + u4;
                const double yspp2 = a * yssp2 - u5;
                const double yspp3 = yssp2 * yssp2;
                const double ypsp1 = ypss1 * yssp1;
                const double ypsp2 = yssp1 * ypss2 + b;
                const double ypsp3 = yssp2 * ypss2;
                const double ypps1 = ypss1 * ysps1;
                const double ypps2 = ysps1 * ypss2 + b;
                const double yppp1 = ypss1 * yspp1;
                const double yppp2 = ypss1 * yspp2 + yspp1 * ypss2 + a * u2;
                const double yppp3 = ypss1 * yspp3 + a * ypsp3 + u3 * yssp2;
                const double yppp4 = ypss2 * yspp3;

                /* ...the Z-terms. */
                const double zssp1 = qzval - z4;
                const double zsps1 = qzval - z3;
                const double zssp2 = pqz * u0;
                const double zpss2 = pqz * u1;
                a = zsps1 + zssp1;
                b = zpss1 * zssp2 + u2;
                const double zspp1 = zsps1 * zssp1 + u4;
                const double zspp2 = a * zssp2 - u5;
                const double zspp3 = zssp2 * zssp2;
                const double zpsp1 = zpss1 * zssp1;
                const double zpsp2 = zssp1 * zpss2 + b;
                const double zpsp3 = zssp2 * zpss2;
                const double zpps1 = zpss1 * zsps1;
                const double zpps2 = zsps1 * zpss2 + b;
                const double zppp1 = zpss1 * zspp1;
                const double zppp2 = zpss1 * zspp2 + zspp1 * zpss2 + a * u2;
                const double zppp3 = zpss1 * zspp3 + a * zpsp3 + u3 * zssp2;
                const double zppp4 = zpss2 * zspp3;

                /* ...assemble the 4-center (AB|CD) type integrals. */
                const double gxxx = xppp1 * boys.f0 + xppp2 * boys.f1 + xppp3 * boys.f2 + xppp4 * boys.f3;
                const double gyyy = yppp1 * boys.f0 + yppp2 * boys.f1 + yppp3 * boys.f2 + yppp4 * boys.f3;
                const double gzzz = zppp1 * boys.f0 + zppp2 * boys.f1 + zppp3 * boys.f2 + zppp4 * boys.f3;
                double aa = xpsp3 * boys.f2;
                double bb = xpsp3 * boys.f3;
                a = xpps1 * boys.f0 + xpps2 * boys.f1 + aa;
                b = xpps1 * boys.f1 + xpps2 * boys.f2 + bb;
                double c = xpsp1 * boys.f0 + xpsp2 * boys.f1 + aa;
                double d = xpsp1 * boys.f1 + xpsp2 * boys.f2 + bb;
                double e = xspp1 * boys.f0 + xspp2 * boys.f1 + xspp3 * boys.f2;
                double f = xspp1 * boys.f1 + xspp2 * boys.f2 + xspp3 * boys.f3;
                const double gxxy = yssp1 * a + yssp2 * b;
                const double gxxz = zssp1 * a + zssp2 * b;
                const double gxyx = ysps1 * c + yssp2 * d;
                const double gxzx = zsps1 * c + zssp2 * d;
                const double gyxx = ypss1 * e + ypss2 * f;
                const double gzxx = zpss1 * e + zpss2 * f;
                aa = ypsp3 * boys.f2;
                bb = ypsp3 * boys.f3;
                a = ypps1 * boys.f0 + ypps2 * boys.f1 + aa;
                b = ypps1 * boys.f1 + ypps2 * boys.f2 + bb;
                c = ypsp1 * boys.f0 + ypsp2 * boys.f1 + aa;
                d = ypsp1 * boys.f1 + ypsp2 * boys.f2 + bb;
                e = yspp1 * boys.f0 + yspp2 * boys.f1 + yspp3 * boys.f2;
                f = yspp1 * boys.f1 + yspp2 * boys.f2 + yspp3 * boys.f3;
                const double gyyx = xssp1 * a + xssp2 * b;
                const double gyyz = zssp1 * a + zssp2 * b;
                const double gyxy = xsps1 * c + xssp2 * d;
                const double gyzy = zsps1 * c + zssp2 * d;
                const double gxyy = xpss1 * e + xpss2 * f;
                const double gzyy = zpss1 * e + zpss2 * f;
                aa = zpsp3 * boys.f2;
                bb = zpsp3 * boys.f3;
                a = zpps1 * boys.f0 + zpps2 * boys.f1 + aa;
                b = zpps1 * boys.f1 + zpps2 * boys.f2 + bb;
                c = zpsp1 * boys.f0 + zpsp2 * boys.f1 + aa;
                d = zpsp1 * boys.f1 + zpsp2 * boys.f2 + bb;
                e = zspp1 * boys.f0 + zspp2 * boys.f1 + zspp3 * boys.f2;
                f = zspp1 * boys.f1 + zspp2 * boys.f2 + zspp3 * boys.f3;
                const double gzzx = xssp1 * a + xssp2 * b;
                const double gzzy = yssp1 * a + yssp2 * b;
                const double gzxz = xsps1 * c + xssp2 * d;
                const double gzyz = ysps1 * c + yssp2 * d;
                const double gxzz = xpss1 * e + xpss2 * f;
                const double gyzz = ypss1 * e + ypss2 * f;
                a = xpss1 * boys.f0 + xpss2 * boys.f1;
                b = xpss1 * boys.f1 + xpss2 * boys.f2;
                c = xpss1 * boys.f2 + xpss2 * boys.f3;
                d = ypss1 * boys.f0 + ypss2 * boys.f1;
                e = ypss1 * boys.f1 + ypss2 * boys.f2;
                f = ypss1 * boys.f2 + ypss2 * boys.f3;
                double g = zpss1 * boys.f0 + zpss2 * boys.f1;
                double h = zpss1 * boys.f1 + zpss2 * boys.f2;
                double r = zpss1 * boys.f2 + zpss2 * boys.f3;
                const double gxyz = ysps1 * (zssp1 * a + zssp2 * b) + yssp2 * (zssp1 * b + zssp2 * c);
                const double gxzy = zsps1 * (yssp1 * a + yssp2 * b) + zssp2 * (yssp1 * b + yssp2 * c);
                const double gyxz = xsps1 * (zssp1 * d + zssp2 * e) + xssp2 * (zssp1 * e + zssp2 * f);
                const double gyzx = zsps1 * (xssp1 * d + xssp2 * e) + zssp2 * (xssp1 * e + xssp2 * f);
                const double gzxy = xsps1 * (yssp1 * g + yssp2 * h) + xssp2 * (yssp1 * h + yssp2 * r);
                const double gzyx = ysps1 * (xssp1 * g + xssp2 * h) + yssp2 * (xssp1 * h + xssp2 * r);
                cbatch[0] += gxxx;
                cbatch[1] += gyxx;
                cbatch[2] += gzxx;
                cbatch[3] += gxyx;
                cbatch[4] += gyyx;
                cbatch[5] += gzyx;
                cbatch[6] += gxzx;
                cbatch[7] += gyzx;
                cbatch[8] += gzzx;
                cbatch[9] += gxxy;
                cbatch[10] += gyxy;
                cbatch[11] += gzxy;
                cbatch[12] += gxyy;
                cbatch[13] += gyyy;
                cbatch[14] += gzyy;
                cbatch[15] += gxzy;
                cbatch[16] += gyzy;
                cbatch[17] += gzzy;
                cbatch[18] += gxxz;
                cbatch[19] += gyxz;
                cbatch[20] += gzxz;
                cbatch[21] += gxyz;
                cbatch[22] += gyyz;
                cbatch[23] += gzyz;
                cbatch[24] += gxzz;
                cbatch[25] += gyzz;
                cbatch[26] += gzzz;
            }
        }
    } else {
        // 2     5   |  (AB|CD)  4-center   ppsp and ppps
        double qxsub, qysub, qzsub;
        if (shell3 == 1) {
            qxsub = x3;
            qysub = y3;
            qzsub = z3;
        } else {
            qxsub = x4;
            qysub = y4;
            qzsub = z4;
        }
        for (uint32_t ij = 0; ij < nij; ij += 1) {
            const double pval = p[ij];
            const double pinv = 1. / pval;
            const double pxval = px[ij];
            const double pyval = py[ij];
            const double pzval = pz[ij];
            const double pscale = scalep[ij];
            const double u4 = pinv * .5;
            const double xsps1 = pxval - x2;
            const double xpss1 = pxval - x1;
            const double ysps1 = pyval - y2;
            const double ypss1 = pyval - y1;
            const double zsps1 = pzval - z2;
            const double zpss1 = pzval - z1;
            for (uint32_t kl = 0; kl < nkl; kl += 1) {
                const double qval = q[kl];
                const double qxval = qx[kl];
                const double qyval = qy[kl];
                const double qzval = qz[kl];
                const double pqmult = pval * qval;
                const double pqplus = pval + qval;
                const double pqpinv = 1. / pqplus;
                const double pqx = pxval - qxval;
                const double pqy = pyval - qyval;
                const double pqz = pzval - qzval;
                const double t = (pqx * pqx + pqy * pqy + pqz * pqz) * pqmult * pqpinv;
                const double scale = pscale * scaleq[kl] / (pqmult * __builtin_sqrt(pqplus));
                const struct Boys0123 boys = boys0123(t, scale);
                const double u0 = pval * pqpinv;
                const double u1 = -qval * pqpinv;
                const double u2 = pqpinv * .5;
                const double u3 = u2 + pqpinv;
                const double u5 = u1 * u4;

                /* ...the X-terms. */
                const double xssp1 = qxval - qxsub;
                const double xssp2 = pqx * u0;
                const double xsps2 = pqx * u1;
                double a = xpss1 + xsps1;
                double b = xssp1 * xsps2 + u2;
                const double xspp1 = xsps1 * xssp1;
                const double xspp2 = xsps1 * xssp2 + b;
                const double xspp3 = xssp2 * xsps2;
                const double xpsp1 = xpss1 * xssp1;
                const double xpsp2 = xpss1 * xssp2 + b;
                const double xpps1 = xpss1 * xsps1 + u4;
                const double xpps2 = a * xsps2 + u5;
                const double xpps3 = xsps2 * xsps2;
                const double xppp1 = xssp1 * xpps1;
                const double xppp2 = xssp1 * xpps2 + xpps1 * xssp2 + a * u2;
                const double xppp3 = xssp1 * xpps3 + a * xspp3 + u3 * xsps2;
                const double xppp4 = xssp2 * xpps3;

                /* ...the Y-terms. */
                const double yssp1 = qyval - qysub;
                const double yssp2 = pqy * u0;
                const double ysps2 = pqy * u1;
                a = ypss1 + ysps1;
                b = yssp1 * ysps2 + u2;
                const double yspp1 = ysps1 * yssp1;
                const double yspp2 = ysps1 * yssp2 + b;
                const double yspp3 = yssp2 * ysps2;
                const double ypsp1 = ypss1 * yssp1;
                const double ypsp2 = ypss1 * yssp2 + b;
                const double ypps1 = ypss1 * ysps1 + u4;
                const double ypps2 = a * ysps2 + u5;
                const double ypps3 = ysps2 * ysps2;
                const double yppp1 = yssp1 * ypps1;
                const double yppp2 = yssp1 * ypps2 + ypps1 * yssp2 + a * u2;
                const double yppp3 = yssp1 * ypps3 + a * yspp3 + u3 * ysps2;
                const double yppp4 = yssp2 * ypps3;

                /* ...the Z-terms. */
                const double zssp1 = qzval - qzsub;
                const double zssp2 = pqz * u0;
                const double zsps2 = pqz * u1;
                a = zpss1 + zsps1;
                b = zssp1 * zsps2 + u2;
                const double zspp1 = zsps1 * zssp1;
                const double zspp2 = zsps1 * zssp2 + b;
                const double zspp3 = zssp2 * zsps2;
                const double zpsp1 = zpss1 * zssp1;
                const double zpsp2 = zpss1 * zssp2 + b;
                const double zpps1 = zpss1 * zsps1 + u4;
                const double zpps2 = a * zsps2 + u5;
                const double zpps3 = zsps2 * zsps2;
                const double zppp1 = zssp1 * zpps1;
                const double zppp2 = zssp1 * zpps2 + zpps1 * zssp2 + a * u2;
                const double zppp3 = zssp1 * zpps3 + a * zspp3 + u3 * zsps2;
                const double zppp4 = zssp2 * zpps3;

                /* ...assemble the 4-center (AB|CD) type integrals. */
                const double gxxx = xppp1 * boys.f0 + xppp2 * boys.f1 + xppp3 * boys.f2 + xppp4 * boys.f3;
                const double gyyy = yppp1 * boys.f0 + yppp2 * boys.f1 + yppp3 * boys.f2 + yppp4 * boys.f3;
                const double gzzz = zppp1 * boys.f0 + zppp2 * boys.f1 + zppp3 * boys.f2 + zppp4 * boys.f3;
                double aa = xspp3 * boys.f2;
                double bb = xspp3 * boys.f3;
                a = xpps1 * boys.f0 + xpps2 * boys.f1 + xpps3 * boys.f2;
                b = xpps1 * boys.f1 + xpps2 * boys.f2 + xpps3 * boys.f3;
                double c = xpsp1 * boys.f0 + xpsp2 * boys.f1 + aa;
                double d = xpsp1 * boys.f1 + xpsp2 * boys.f2 + bb;
                double e = xspp1 * boys.f0 + xspp2 * boys.f1 + aa;
                double f = xspp1 * boys.f1 + xspp2 * boys.f2 + bb;
                const double gxxy = yssp1 * a + yssp2 * b;
                const double gxxz = zssp1 * a + zssp2 * b;
                const double gxyx = ysps1 * c + ysps2 * d;
                const double gxzx = zsps1 * c + zsps2 * d;
                const double gyxx = ypss1 * e + ysps2 * f;
                const double gzxx = zpss1 * e + zsps2 * f;
                aa = yspp3 * boys.f2;
                bb = yspp3 * boys.f3;
                a = ypps1 * boys.f0 + ypps2 * boys.f1 + ypps3 * boys.f2;
                b = ypps1 * boys.f1 + ypps2 * boys.f2 + ypps3 * boys.f3;
                c = ypsp1 * boys.f0 + ypsp2 * boys.f1 + aa;
                d = ypsp1 * boys.f1 + ypsp2 * boys.f2 + bb;
                e = yspp1 * boys.f0 + yspp2 * boys.f1 + aa;
                f = yspp1 * boys.f1 + yspp2 * boys.f2 + bb;
                const double gyyx = xssp1 * a + xssp2 * b;
                const double gyyz = zssp1 * a + zssp2 * b;
                const double gyxy = xsps1 * c + xsps2 * d;
                const double gyzy = zsps1 * c + zsps2 * d;
                const double gxyy = xpss1 * e + xsps2 * f;
                const double gzyy = zpss1 * e + zsps2 * f;
                aa = zspp3 * boys.f2;
                bb = zspp3 * boys.f3;
                a = zpps1 * boys.f0 + zpps2 * boys.f1 + zpps3 * boys.f2;
                b = zpps1 * boys.f1 + zpps2 * boys.f2 + zpps3 * boys.f3;
                c = zpsp1 * boys.f0 + zpsp2 * boys.f1 + aa;
                d = zpsp1 * boys.f1 + zpsp2 * boys.f2 + bb;
                e = zspp1 * boys.f0 + zspp2 * boys.f1 + aa;
                f = zspp1 * boys.f1 + zspp2 * boys.f2 + bb;
                const double gzzx = xssp1 * a + xssp2 * b;
                const double gzzy = yssp1 * a + yssp2 * b;
                const double gzxz = xsps1 * c + xsps2 * d;
                const double gzyz = ysps1 * c + ysps2 * d;
                const double gxzz = xpss1 * e + xsps2 * f;
                const double gyzz = ypss1 * e + ysps2 * f;
                a = xpss1 * boys.f0 + xsps2 * boys.f1;
                b = xpss1 * boys.f1 + xsps2 * boys.f2;
                c = xpss1 * boys.f2 + xsps2 * boys.f3;
                d = ypss1 * boys.f0 + ysps2 * boys.f1;
                e = ypss1 * boys.f1 + ysps2 * boys.f2;
                f = ypss1 * boys.f2 + ysps2 * boys.f3;
                double g = zpss1 * boys.f0 + zsps2 * boys.f1;
                double h = zpss1 * boys.f1 + zsps2 * boys.f2;
                double r = zpss1 * boys.f2 + zsps2 * boys.f3;
                const double gxyz = ysps1 * (zssp1 * a + zssp2 * b) + ysps2 * (zssp1 * b + zssp2 * c);
                const double gxzy = zsps1 * (yssp1 * a + yssp2 * b) + zsps2 * (yssp1 * b + yssp2 * c);
                const double gyxz = xsps1 * (zssp1 * d + zssp2 * e) + xsps2 * (zssp1 * e + zssp2 * f);
                const double gyzx = zsps1 * (xssp1 * d + xssp2 * e) + zsps2 * (xssp1 * e + xssp2 * f);
                const double gzxy = xsps1 * (yssp1 * g + yssp2 * h) + xsps2 * (yssp1 * h + yssp2 * r);
                const double gzyx = ysps1 * (xssp1 * g + xssp2 * h) + ysps2 * (xssp1 * h + xssp2 * r);
                cbatch[0] += gxxx;
                cbatch[1] += gyxx;
                cbatch[2] += gzxx;
                cbatch[3] += gxyx;
                cbatch[4] += gyyx;
                cbatch[5] += gzyx;
                cbatch[6] += gxzx;
                cbatch[7] += gyzx;
                cbatch[8] += gzzx;
                cbatch[9] += gxxy;
                cbatch[10] += gyxy;
                cbatch[11] += gzxy;
                cbatch[12] += gxyy;
                cbatch[13] += gyyy;
                cbatch[14] += gzyy;
                cbatch[15] += gxzy;
                cbatch[16] += gyzy;
                cbatch[17] += gzzy;
                cbatch[18] += gxxz;
                cbatch[19] += gyxz;
                cbatch[20] += gzxz;
                cbatch[21] += gxyz;
                cbatch[22] += gyyz;
                cbatch[23] += gzyz;
                cbatch[24] += gxzz;
                cbatch[25] += gyzz;
                cbatch[26] += gzzz;
            }
        }
    }
}
